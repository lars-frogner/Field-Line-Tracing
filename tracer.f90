module tracer_mod
implicit none

private

integer, parameter :: SP = kind(0.0)

public :: trace

contains

subroutine trace(x_initial, y_initial, z_initial, &
                 n_output_points,                 &
                 x_out, y_out, z_out,             &
                 aux_out,                         &
                 fieldline_length,                &
                 decoupling_index)

    use tracer_params_mod, only : direction,         &
                                  nonuniform_output, &
                                  decoupling_rate,   &
                                  slope_correction,  &
                                  max_output_points, &
                                  ds_out,            &
                                  n_aux,             &
                                  z_max

    use tracer_base_mod, only : initialize_tracer, &
                                adjust_stepsize,   &
                                compute_aux_output

    use stepper_mod, only : initialize_stepper,        &
                            step,                      &
                            compute_error,             &
                            prepare_outputting,        &
                            get_interpolated_position, &
                            update_directions

    use mesh_mod, only : zm, u_l

    real(SP), intent(in) :: x_initial, y_initial, z_initial

    integer,               intent(out) :: n_output_points
    real(SP), allocatable, intent(out) :: x_out(:), y_out(:), z_out(:)
    real(SP), allocatable, intent(out) :: aux_out(:, :)
    real(SP),              intent(out) :: fieldline_length
    integer,               intent(out) :: decoupling_index

    real(SP) :: s_current
    real(SP) :: s_old

    real(SP) :: error_current
    real(SP) :: error_old

    real(SP) :: ds_current
    real(SP) :: ds_next

    real(SP) :: x_current, y_current, z_current
    real(SP) :: x_old, y_old, z_old

    integer  :: x_idx_current, y_idx_current, z_idx_current
    integer  :: x_idx_old, y_idx_old, z_idx_old

    real(SP) :: dx_ds_initial, dy_ds_initial, dz_ds_initial

    logical  :: decoupled
    real(SP) :: s_first_decoupling
    real(SP) :: B_weight, z_weight

    integer  :: n_extension_points

    real(SP), allocatable :: out_temp(:)
    real(SP), allocatable :: out_aux_temp(:, :)

    integer  :: n
    real(SP) :: s_out
    integer  :: x_idx_out, y_idx_out, z_idx_out
    logical  :: terminated
    logical  :: accepted

    decoupling_index = -1

    ! ************************ Allocate output arrays ************************

    if (allocated(x_out))   deallocate(x_out)
    if (allocated(y_out))   deallocate(y_out)
    if (allocated(z_out))   deallocate(z_out)
    if (allocated(aux_out)) deallocate(aux_out)

    allocate(x_out(max_output_points), &
             y_out(max_output_points), &
             z_out(max_output_points), &
             aux_out(max_output_points, n_aux))

    ! ************************** Initialize stepper **************************

    call initialize_tracer(x_initial, y_initial, z_initial, &
                           ds_current,                      &
                           error_old,                       &
                           s_old,                           &
                           x_old, y_old, z_old,             &
                           x_idx_old, y_idx_old, z_idx_old, &
                           dx_ds_initial, dy_ds_initial, dz_ds_initial)

    call initialize_stepper(dx_ds_initial, dy_ds_initial, dz_ds_initial)

    ! ************************** Set initial values **************************

    x_out(1) = x_initial
    y_out(1) = y_initial
    z_out(1) = z_initial

    call compute_aux_output(x_out(1), y_out(1), z_out(1),    &
                            x_idx_old, y_idx_old, z_idx_old, &
                            0.0,                             &
                            aux_out(1, :),                   &
                            terminated)

    ! *************************** Perform stepping ***************************

    n = 2
    s_out = ds_out
    accepted = .true.
    terminated = .false.

    s_first_decoupling = -1.0
    B_weight = direction
    z_weight = 0.0
    decoupled = .false.

    tracing: do

        step_attempts: do

            ! Perform a trial step with the current step size
            call step(x_old, y_old, z_old,                         &
                      x_idx_old, y_idx_old, z_idx_old,             &
                      ds_current,                                  &
                      B_weight, z_weight,                          &
                      x_current, y_current, z_current,             &
                      x_idx_current, y_idx_current, z_idx_current, &
                      terminated)

            ! Stop if we hit a non-periodic boundary
            if (terminated) then
                n_output_points = n - 1
                fieldline_length = s_old
                exit tracing
            end if

            ! Find the estimated error
            call compute_error(x_current, y_current, z_current, &
                               x_old, y_old, z_old,             &
                               error_current)

            ! Evaluate error and adjust step size accordingly
            call adjust_stepsize(ds_current,    &
                                 error_current, &
                                 error_old,     &
                                 accepted,      &
                                 ds_next)

            ! If accepted, the step is valid and we continue
            if (accepted) exit step_attempts

            ! Otherwise, apply the new step size and try again
            ds_current = ds_next

        end do step_attempts

        ! Find new position along fieldline
        s_current = s_old + ds_current

        if (nonuniform_output) then

            x_out(n) = x_current
            y_out(n) = y_current
            z_out(n) = z_current

            x_idx_out = x_idx_current
            y_idx_out = y_idx_current
            z_idx_out = z_idx_current

            ! Get auxiliary output
            call compute_aux_output(x_out(n), y_out(n), z_out(n),    &
                                    x_idx_out, y_idx_out, z_idx_out, &
                                    z_weight,                        &
                                    aux_out(n, :),                   &
                                    decoupled)

            ! Adjust weights for smooth decoupling transition
            if (decoupled) then

                if (s_first_decoupling < 0.0) then
                    s_first_decoupling = s_current
                    decoupling_index = n
                end if

                z_weight = (s_current - s_first_decoupling)*decoupling_rate
                B_weight = direction*(1 - z_weight)

                ! Terminate when the direction is straight down
                if (z_weight >= 1.0) then

                    n_output_points = n
                    fieldline_length = s_current
                    exit tracing

                end if

            end if

            ! Increment number of output values
            n = n + 1

            ! Stop if the output array is full
            if (n > max_output_points) then
                n_output_points = max_output_points
                fieldline_length = s_current
                terminated = .true.
                exit tracing
            end if

        else if (s_current >= s_out) then

            ! We have passed the next output position

            ! Compute interpolation coefficients for output
            call prepare_outputting(ds_current, &
                                    x_old, y_old, z_old)

            output_steps: do

                ! Write output until we reach the computed position

                ! Get interpolated output position
                call get_interpolated_position(s_out,                           &
                                               s_old,                           &
                                               ds_current,                      &
                                               x_idx_old, y_idx_old, z_idx_old, &
                                               x_out(n), y_out(n), z_out(n),    &
                                               x_idx_out, y_idx_out, z_idx_out)

                ! Get auxiliary output
                call compute_aux_output(x_out(n), y_out(n), z_out(n),    &
                                        x_idx_out, y_idx_out, z_idx_out, &
                                        z_weight,                        &
                                        aux_out(n, :),                   &
                                        decoupled)

                ! Adjust weights for smooth decoupling transition
                if (decoupled) then

                    if (s_first_decoupling < 0.0) then
                        s_first_decoupling = s_out
                        decoupling_index = n
                    end if

                    z_weight = (s_out - s_first_decoupling)*decoupling_rate
                    B_weight = direction*(1 - z_weight)

                    ! Terminate when the direction is straight down
                    if (z_weight >= 1.0) then

                        n_output_points = n
                        fieldline_length = s_out
                        exit tracing

                    end if

                end if

                ! Increment number of output values
                n = n + 1

                ! Stop if the output array is full
                if (n > max_output_points) then
                    n_output_points = max_output_points
                    fieldline_length = s_out
                    terminated = .true.
                    exit tracing
                end if

                ! Set next output position
                s_out = s_out + ds_out

                ! When the next output position is past the newest position, we can continue
                if (s_out > s_current) exit output_steps

            end do output_steps

        end if

        ! Move to the new position

        s_old = s_current

        x_old = x_current
        y_old = y_current
        z_old = z_current

        x_idx_old = x_idx_current
        y_idx_old = y_idx_current
        z_idx_old = z_idx_current

        call update_directions()

        ! Use the adjusted step size for the next step
        ds_current = ds_next

    end do tracing

    ! ****************** Extend straight down to the bottom ******************

    if (.not. terminated) then

        n_extension_points = floor((z_max - z_out(n_output_points))/ds_out)

        if (n_output_points + n_extension_points > max_output_points) then
            n_extension_points = max_output_points - n_output_points
        end if

        do n = n_output_points + 1, n_output_points + n_extension_points

            x_out(n) = x_out(n_output_points)
            y_out(n) = y_out(n_output_points)
            z_out(n) = z_out(n-1) + ds_out

            if (zm(z_idx_out+1) <= z_out(n)) then
                do
                    z_idx_out = z_idx_out + 1
                    if (zm(z_idx_out+1) > z_out(n)) exit
                end do
            end if

            call compute_aux_output(x_out(n), y_out(n), z_out(n),    &
                                    x_idx_out, y_idx_out, z_idx_out, &
                                    1.0,                             &
                                    aux_out(n, :),                   &
                                    decoupled)

        end do

        n_output_points = n_output_points + n_extension_points
        fieldline_length = fieldline_length + n_extension_points*ds_out

    end if

    ! ************************* Resize output arrays *************************

    if (n_output_points < max_output_points) then

        allocate(out_temp(n_output_points))

        out_temp = x_out(1:n_output_points)
        deallocate(x_out)
        allocate(x_out(n_output_points))
        x_out = u_l*out_temp

        out_temp = y_out(1:n_output_points)
        deallocate(y_out)
        allocate(y_out(n_output_points))
        y_out = u_l*out_temp

        out_temp = z_out(1:n_output_points)
        deallocate(z_out)
        allocate(z_out(n_output_points))
        z_out = u_l*out_temp

        deallocate(out_temp)

        allocate(out_aux_temp(n_output_points, n_aux))

        out_aux_temp = aux_out(1:n_output_points, :)
        deallocate(aux_out)
        allocate(aux_out(n_output_points, n_aux))
        aux_out = out_aux_temp

    end if

end subroutine trace

end module tracer_mod
