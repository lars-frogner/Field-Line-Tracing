module stepper_mod
implicit none

! The Bogacki–Shampine scheme, a third-order Runge-Kutta method with error
! estimation through an embedded second-order step.

private

integer, parameter :: SP = kind(0.0)

! ***************************** Coefficients *****************************

real(SP) :: a2(1) = [ 1.0/2.0   ]
real(SP) :: a3(2) = [ 0.0     , &
                      3.0/4.0   ]
real(SP) :: a4(3) = [ 2.0/9.0 , &
                      1.0/3.0 , &
                      4.0/9.0   ]

real(SP) ::  e(4) = [-5.0/72.0, &
                      1.0/12.0, &
                      1.0/9.0 , &
                     -1.0/8.0   ]

real(SP) :: cx(4)
real(SP) :: cy(4)
real(SP) :: cz(4)

! ****************************** Variables *******************************

real(SP) :: dx_ds_current(4)
real(SP) :: dy_ds_current(4)
real(SP) :: dz_ds_current(4)

real(SP) :: dx_ds_old(4)
real(SP) :: dy_ds_old(4)
real(SP) :: dz_ds_old(4)

real(SP) :: dx, dy, dz
real(SP) :: delta_x, delta_y, delta_z


public :: initialize_stepper, &
          step,               &
          update_directions,  &
          compute_error,      &
          prepare_outputting, &
          get_interpolated_position

contains

subroutine initialize_stepper(dx_ds_initial, dy_ds_initial, dz_ds_initial)

    use tracer_params_mod, only : beta, alpha

    real(SP), intent(in) :: dx_ds_initial, dy_ds_initial, dz_ds_initial

    alpha = 1.0/3.0 - 0.75*beta

    dx_ds_old(4) = dx_ds_initial
    dy_ds_old(4) = dy_ds_initial
    dz_ds_old(4) = dz_ds_initial

end subroutine initialize_stepper


subroutine step(x_old, y_old, z_old,                         &
                x_idx_old, y_idx_old, z_idx_old,             &
                ds_current,                                  &
                B_weight, z_weight,                          &
                x_current, y_current, z_current,             &
                x_idx_current, y_idx_current, z_idx_current, &
                terminate)

    use tracer_base_mod, only : apply_boundary_conditions, &
                                update_position_indices,   &
                                interpolate_field

    real(SP), intent(in) :: x_old, y_old, z_old
    integer,  intent(in) :: x_idx_old, y_idx_old, z_idx_old
    real(SP), intent(in) :: ds_current
    real(SP), intent(in) :: B_weight, z_weight

    real(SP), intent(out) :: x_current, y_current, z_current
    integer,  intent(out) :: x_idx_current, y_idx_current, z_idx_current
    logical,  intent(out) :: terminate

    real(SP) :: norm
    real(SP) :: Bx, By, Bz

    x_idx_current = x_idx_old
    y_idx_current = y_idx_old
    z_idx_current = z_idx_old

    ! ******************************** Step 1 ********************************

    dx_ds_current(1) = dx_ds_old(4)
    dy_ds_current(1) = dy_ds_old(4)
    dz_ds_current(1) = dz_ds_old(4)

    x_current = x_old + ds_current*a2(1)*dx_ds_current(1)
    y_current = y_old + ds_current*a2(1)*dy_ds_current(1)
    z_current = z_old + ds_current*a2(1)*dz_ds_current(1)

    call apply_boundary_conditions(x_current, y_current, z_current,             &
                                   x_idx_current, y_idx_current, z_idx_current, &
                                   terminate)

    call update_position_indices(x_current, y_current, z_current, &
                                 x_idx_current, y_idx_current, z_idx_current)

    call interpolate_field(x_idx_current, y_idx_current, z_idx_current, &
                           x_current, y_current, z_current,             &
                           Bx, By, Bz)

    ! ******************************** Step 2 ********************************

    Bx = B_weight*Bx
    By = B_weight*By
    Bz = B_weight*Bz + z_weight

    norm = 1.0/sqrt(Bx*Bx + By*By + Bz*Bz)

    dx_ds_current(2) = Bx*norm
    dy_ds_current(2) = By*norm
    dz_ds_current(2) = Bz*norm

    x_current = x_old + ds_current*(a3(2)*dx_ds_current(2))
    y_current = y_old + ds_current*(a3(2)*dy_ds_current(2))
    z_current = z_old + ds_current*(a3(2)*dz_ds_current(2))

    call apply_boundary_conditions(x_current, y_current, z_current,             &
                                   x_idx_current, y_idx_current, z_idx_current, &
                                   terminate)

    call update_position_indices(x_current, y_current, z_current, &
                                 x_idx_current, y_idx_current, z_idx_current)

    call interpolate_field(x_idx_current, y_idx_current, z_idx_current, &
                           x_current, y_current, z_current,             &
                           Bx, By, Bz)

    ! ******************************** Step 3 ********************************

    Bx = B_weight*Bx
    By = B_weight*By
    Bz = B_weight*Bz + z_weight

    norm = 1.0/sqrt(Bx*Bx + By*By + Bz*Bz)

    dx_ds_current(3) = Bx*norm
    dy_ds_current(3) = By*norm
    dz_ds_current(3) = Bz*norm

    dx = ds_current*(a4(1)*dx_ds_current(1) + &
                     a4(2)*dx_ds_current(2) + &
                     a4(3)*dx_ds_current(3))

    dy = ds_current*(a4(1)*dy_ds_current(1) + &
                     a4(2)*dy_ds_current(2) + &
                     a4(3)*dy_ds_current(3))

    dz = ds_current*(a4(1)*dz_ds_current(1) + &
                     a4(2)*dz_ds_current(2) + &
                     a4(3)*dz_ds_current(3))

    x_current = x_old + dx
    y_current = y_old + dy
    z_current = z_old + dz

    ! Set actual ds_current that was used in the final step (might deviate from original ds_current)
    !ds_current = sqrt(dx*dx + dy*dy)

    call apply_boundary_conditions(x_current, y_current, z_current,             &
                                   x_idx_current, y_idx_current, z_idx_current, &
                                   terminate)

    call update_position_indices(x_current, y_current, z_current, &
                                 x_idx_current, y_idx_current, z_idx_current)

    call interpolate_field(x_idx_current, y_idx_current, z_idx_current, &
                           x_current, y_current, z_current,             &
                           Bx, By, Bz)

    ! ************************** Prepare next step ***************************

    Bx = B_weight*Bx
    By = B_weight*By
    Bz = B_weight*Bz + z_weight

    norm = 1.0/sqrt(Bx*Bx + By*By + Bz*Bz)

    dx_ds_current(4) = Bx*norm
    dy_ds_current(4) = By*norm
    dz_ds_current(4) = Bz*norm

    ! **************************** Compute error *****************************

    delta_x = ds_current*(e(1)*dx_ds_current(1) + &
                          e(2)*dx_ds_current(2) + &
                          e(3)*dx_ds_current(3) + &
                          e(4)*dx_ds_current(4))

    delta_y = ds_current*(e(1)*dy_ds_current(1) + &
                          e(2)*dy_ds_current(2) + &
                          e(3)*dy_ds_current(3) + &
                          e(4)*dy_ds_current(4))

    delta_z = ds_current*(e(1)*dz_ds_current(1) + &
                          e(2)*dz_ds_current(2) + &
                          e(3)*dz_ds_current(3) + &
                          e(4)*dz_ds_current(4))

end subroutine step


subroutine update_directions()

    dx_ds_old = dx_ds_current
    dy_ds_old = dy_ds_current
    dz_ds_old = dz_ds_current

end subroutine update_directions


subroutine compute_error(x_current, y_current, z_current, &
                         x_old, y_old, z_old,             &
                         error_current)

    use tracer_params_mod, only : abs_tolerance, &
                                  rel_tolerance

    real(SP), intent(in) :: x_current, y_current, z_current
    real(SP), intent(in) :: x_old, y_old, z_old

    real(SP), intent(out) :: error_current

    real(SP) :: error_x, error_y, error_z

    error_x = delta_x/(abs_tolerance + &
                       rel_tolerance*max(abs(x_old), abs(x_current)))

    error_y = delta_y/(abs_tolerance + &
                       rel_tolerance*max(abs(y_old), abs(y_current)))

    error_z = delta_z/(abs_tolerance + &
                       rel_tolerance*max(abs(z_old), abs(z_current)))

    error_current = sqrt(0.5*(error_x*error_x + error_y*error_y + error_z*error_z))

end subroutine compute_error


subroutine prepare_outputting(ds_current, &
                              x_old, y_old, z_old)

    ! Computes coefficients for the third-order interpolating polynomial
    ! between s_old and s.

    real(SP), intent(in) :: ds_current
    real(SP), intent(in) :: x_old, y_old, z_old

    cx(1) = x_old
    cy(1) = y_old
    cz(1) = z_old

    cx(2) = dx
    cy(2) = dy
    cz(2) = dz

    cx(3) = ds_current*dx_ds_current(1)
    cy(3) = ds_current*dy_ds_current(1)
    cz(3) = ds_current*dz_ds_current(1)

    cx(4) = ds_current*dx_ds_current(4)
    cy(4) = ds_current*dy_ds_current(4)
    cz(4) = ds_current*dz_ds_current(4)

end subroutine prepare_outputting


subroutine get_interpolated_position(s_out,                           &
                                     s_old,                           &
                                     ds_current,                      &
                                     x_idx_old, y_idx_old, z_idx_old, &
                                     x_out, y_out, z_out,             &
                                     x_idx_out, y_idx_out, z_idx_out)

    use tracer_base_mod, only : apply_boundary_conditions, &
                                update_position_indices

    real(SP), intent(in) :: s_out
    real(SP), intent(in) :: s_old
    real(SP), intent(in) :: ds_current
    integer,  intent(in) :: x_idx_old, y_idx_old, z_idx_old

    real(SP), intent(out) :: x_out, y_out, z_out
    integer,  intent(out) :: x_idx_out, y_idx_out, z_idx_out

    real(SP) :: t, t_minus_one, t_times_t_minus_one, one_minus_two_t
    logical  :: terminate

    t = (s_out - s_old)/ds_current
    t_minus_one = t - 1.0
    t_times_t_minus_one = t*t_minus_one
    one_minus_two_t = -(t + t_minus_one)

    x_out = cx(1) + t*cx(2) +                                    &
                    t_times_t_minus_one*(one_minus_two_t*cx(2) + &
                                         t_minus_one*cx(3) +     &
                                         t*cx(4))

    y_out = cy(1) + t*cy(2) +                                    &
                    t_times_t_minus_one*(one_minus_two_t*cy(2) + &
                                         t_minus_one*cy(3) +     &
                                         t*cy(4))

    z_out = cz(1) + t*cz(2) +                                    &
                    t_times_t_minus_one*(one_minus_two_t*cz(2) + &
                                         t_minus_one*cz(3) +     &
                                         t*cz(4))

    x_idx_out = x_idx_old
    y_idx_out = y_idx_old
    z_idx_out = z_idx_old

    call apply_boundary_conditions(x_out, y_out, z_out,             &
                                   x_idx_out, y_idx_out, z_idx_out, &
                                   terminate)

    call update_position_indices(x_out, y_out, z_out, &
                                 x_idx_out, y_idx_out, z_idx_out)

end subroutine get_interpolated_position

end module stepper_mod
