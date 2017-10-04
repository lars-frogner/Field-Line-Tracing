module tracer_base_mod
implicit none

private

integer, parameter :: SP = kind(0.0)

public :: initialize_tracer,    &
          adjust_stepsize,      &
          update_position_data, &
          interpolate_field,    &
          compute_aux_output

contains

subroutine initialize_tracer(x_initial, y_initial, z_initial, &
                             ds_current,                      &
                             error_old,                       &
                             s_old,                           &
                             x_old, y_old, z_old,             &
                             x_idx_old, y_idx_old, z_idx_old, &
                             dx_ds_initial, dy_ds_initial, dz_ds_initial)

    use tracer_params_mod, only : direction,             &
                                  interp_order,          &
                                  interp_bias,           &
                                  s_initial,             &
                                  ds_initial,            &
                                  error_initial,         &
                                  n_interp,              &
                                  start_offset,          &
                                  x_min, x_max, x_range, &
                                  y_min, y_max, y_range, &
                                  z_min, z_max, z_range

    use mesh_mod, only : xs, xe, &
                         ys, ye, &
                         zs, ze, &
                         xm, ym, zm

    real(SP), intent(in)  :: x_initial, y_initial, z_initial

    real(SP), intent(out) :: ds_current
    real(SP), intent(out) :: error_old
    real(SP), intent(out) :: s_old
    real(SP), intent(out) :: x_old, y_old, z_old
    integer,  intent(out) :: x_idx_old, y_idx_old, z_idx_old
    real(SP), intent(out) :: dx_ds_initial, dy_ds_initial, dz_ds_initial

    integer  :: minloc_arr(1)
    real(SP) :: norm
    real(SP) :: Bx, By, Bz

    ! ************************ Set derived parameters ************************

    n_interp = interp_order + 1
    start_offset = interp_bias + 1 - n_interp/2

    x_min = xm(xs)
    y_min = ym(ys)
    z_min = zm(zs)

    x_max = xm(xe)
    y_max = ym(ye)
    z_max = zm(ze)

    x_range = x_max - x_min
    y_range = y_max - y_min
    z_range = z_max - z_min

    ! ************************* Initialize step size *************************

    ds_current = ds_initial
    error_old = error_initial

    ! ************************* Initialize position **************************

    s_old = s_initial

    x_old = x_initial
    y_old = y_initial
    z_old = z_initial

    ! ************************** Find start indices **************************

    minloc_arr = minloc(xm(xs:xe) - x_old, mask=(xm(xs:xe) > x_old))
    x_idx_old = xs + minloc_arr(1) - 2

    minloc_arr = minloc(ym(ys:ye) - y_old, mask=(ym(ys:ye) > y_old))
    y_idx_old = ys + minloc_arr(1) - 2

    minloc_arr = minloc(zm(zs:ze) - z_old, mask=(zm(zs:ze) > z_old))
    z_idx_old = zs + minloc_arr(1) - 2

    ! ************************ Find initial direction ************************

    call interpolate_field(x_idx_old, y_idx_old, z_idx_old, &
                           x_old, y_old, z_old,             &
                           Bx, By, Bz)

    norm = direction/sqrt(Bx*Bx + By*By + Bz*Bz)

    dx_ds_initial = Bx*norm
    dy_ds_initial = By*norm
    dz_ds_initial = Bz*norm

end subroutine initialize_tracer


subroutine adjust_stepsize(ds_current,    &
                           error_current, &
                           error_old,     &
                           accepted,      &
                           ds_next)

    use tracer_params_mod, only : constant_stepsize, &
                                  scale_min,         &
                                  scale_max,         &
                                  safety_factor,     &
                                  alpha, beta

    real(SP), intent(in)    :: ds_current
    real(SP), intent(in)    :: error_current
    real(SP), intent(inout) :: error_old
    logical,  intent(inout) :: accepted
    real(SP), intent(out)   :: ds_next

    real(SP)               :: scale

    if (constant_stepsize) then

        accepted = .true.

    else if (error_current <= 1.0) then

        if (error_current < 1.0e-9) then ! Avoid division by zero
            scale = scale_max
        else

            scale = safety_factor*(error_old**beta)/(error_current**alpha)

            if (scale < scale_min) then
                scale = scale_min
            else if (scale > scale_max) then
                scale = scale_max
            end if

        end if

        if (.not. accepted) then
            ! Don't increase ds if the previous attempt was rejected
            ds_next = ds_current*min(scale, 1.0)
        else
            ds_next = ds_current*scale
        end if

        accepted = .true.

    else

        ds_next = max(safety_factor/(error_current**alpha), &
                      scale_min)*ds_current
        accepted = .false.

    end if

    error_old = error_current

end subroutine adjust_stepsize


subroutine update_position_data(x, y, z,             &
                                x_idx, y_idx, z_idx, &
                                terminated,          &
                                terminated_here)

    use tracer_params_mod, only : x_min, x_max, x_range,  &
                                  y_min, y_max, y_range,  &
                                  z_min, z_max, z_range

    use mesh_mod, only : xs, xe,     &
                         ys, ye,     &
                         zs, ze,     &
                         xm, ym, zm, &
                         periodic_x, periodic_y, periodic_z

    real(SP), intent(inout) :: x, y, z
    integer,  intent(inout) :: x_idx, y_idx, z_idx
    logical,  intent(inout) :: terminated
    logical,  intent(out)   :: terminated_here

    if ((.not. periodic_x .and. (x < x_min .or. x >= x_max)) .or. &
        (.not. periodic_y .and. (y < y_min .or. y >= y_max)) .or. &
        (.not. periodic_z .and. (z < z_min .or. z >= z_max))) then

        terminated = .true.
        terminated_here = .true.

    else

        terminated_here = .false.

    end if

    if (.not. terminated) then

        if (periodic_x) then

            if (x < x_min) then
                x = x_min + modulo(x - x_min,  x_range)
                x_idx = xe - 1
            else if (x >= x_max) then
                x = x_min + modulo(x - x_min,  x_range)
                x_idx = xs
            end if

        end if

        if (periodic_y) then

            if (y < y_min) then
                y = y_min + modulo(y - y_min,  y_range)
                y_idx = ye - 1
            else if (y >= y_max) then
                y = y_min + modulo(y - y_min,  y_range)
                y_idx = ys
            end if

        end if

        if (periodic_z) then

            if (z < z_min) then
                z = z_min + modulo(z - z_min,  z_range)
                z_idx = ze - 1
            else if (z >= z_max) then
                z = z_min + modulo(z - z_min,  z_range)
                z_idx = zs
            end if

        end if

        if (xm(x_idx+1) <= x) then
            do
                x_idx = x_idx + 1
                if (xm(x_idx+1) > x) exit
            end do
        else if (xm(x_idx) > x) then
            do
                x_idx = x_idx - 1
                if (xm(x_idx) <= x) exit
            end do
        end if

        if (ym(y_idx+1) <= y) then
            do
                y_idx = y_idx + 1
                if (ym(y_idx+1) > y) exit
            end do
        else if (ym(y_idx) > y) then
            do
                y_idx = y_idx - 1
                if (ym(y_idx) <= y) exit
            end do
        end if

        if (zm(z_idx+1) <= z) then
            do
                z_idx = z_idx + 1
                if (zm(z_idx+1) > z) exit
            end do
        else if (zm(z_idx) > z) then
            do
                z_idx = z_idx - 1
                if (zm(z_idx) <= z) exit
            end do
        end if

    end if

end subroutine update_position_data


subroutine interpolate_field(x_idx, y_idx, z_idx, &
                             x, y, z,             &
                             Bx, By, Bz)

    use tracer_params_mod, only : n_interp, &
                                  start_offset

    use mesh_mod, only : xsb, xeb,         &
                         ysb, yeb,         &
                         zsb, zeb,         &
                         xm, ym, zm,       &
                         xmdn, ymdn, zmdn, &
                         bxm, bym, bzm

    use poly_interpolation_mod, only : poly_interpolate_3d

    integer,  intent(in)  :: x_idx, y_idx, z_idx
    real(SP), intent(in)  :: x, y, z
    real(SP), intent(out) :: Bx, By, Bz

    integer :: x_start_idx, y_start_idx, z_start_idx

    x_start_idx = x_idx + start_offset
    y_start_idx = y_idx + start_offset
    z_start_idx = z_idx + start_offset

    Bx = poly_interpolate_3d(xsb, xeb,                              &
                             ysb, yeb,                              &
                             zsb, zeb,                              &
                             xmdn, ym, zm,                          &
                             bxm,                                   &
                             x, y, z,                               &
                             x_start_idx, y_start_idx, z_start_idx, &
                             n_interp)

    By = poly_interpolate_3d(xsb, xeb,                              &
                             ysb, yeb,                              &
                             zsb, zeb,                              &
                             xm, ymdn, zm,                          &
                             bym,                                   &
                             x, y, z,                               &
                             x_start_idx, y_start_idx, z_start_idx, &
                             n_interp)

    Bz = poly_interpolate_3d(xsb, xeb,                              &
                             ysb, yeb,                              &
                             zsb, zeb,                              &
                             xm, ym, zmdn,                          &
                             bzm,                                   &
                             x, y, z,                               &
                             x_start_idx, y_start_idx, z_start_idx, &
                             n_interp)

end subroutine interpolate_field


subroutine compute_aux_output(x_out, y_out, z_out,             &
                              x_idx_out, y_idx_out, z_idx_out, &
                              z_weight,                        &
                              aux_out,                         &
                              decoupled)

    use tracer_params_mod, only : slope_correction, &
                                  n_aux, &
                                  aux_names

    real(SP), intent(in)  :: x_out, y_out, z_out
    integer,  intent(in)  :: x_idx_out, y_idx_out, z_idx_out
    real(SP), intent(in)  :: z_weight

    real(SP), intent(out)   :: aux_out(n_aux)
    logical,  intent(inout) :: decoupled

    integer :: i

    do i = 1, n_aux

        select case (trim(aux_names(i)))

            case ('r')
                aux_out(i) = get_density(x_out, y_out, z_out, &
                                         x_idx_out, y_idx_out, z_idx_out)
            case ('r_eff')
                aux_out(i) = get_density(x_out, y_out, z_out, &
                                         x_idx_out, y_idx_out, z_idx_out)*(1 + slope_correction*z_weight)
            case ('e')
                aux_out(i) = get_internal_energy(x_out, y_out, z_out, &
                                                 x_idx_out, y_idx_out, z_idx_out)
            case ('Pg')
                aux_out(i) = get_gas_pressure(x_out, y_out, z_out, &
                                              x_idx_out, y_idx_out, z_idx_out)
            case ('PB')
                aux_out(i) = get_magnetic_pressure(x_out, y_out, z_out, &
                                                   x_idx_out, y_idx_out, z_idx_out)
            case ('beta')
                aux_out(i) = get_plasma_beta(x_out, y_out, z_out, &
                                             x_idx_out, y_idx_out, z_idx_out)
                if (.not. decoupled) then
                    decoupled = evaluate_decoupling_condition(z_out, aux_out(i))
                end if

        end select

    end do

end subroutine compute_aux_output


function get_density(x_out, y_out, z_out, &
                     x_idx_out, y_idx_out, z_idx_out) result(r_cgs)

    use mesh_mod, only : xsb, xeb,         &
                         ysb, yeb,         &
                         zsb, zeb,         &
                         u_r,              &
                         xmdn, ymdn, zmdn, &
                         rm

    use tracer_params_mod, only : n_interp, &
                                  start_offset

    use poly_interpolation_mod, only : poly_interpolate_3d

    real(SP), intent(in)  :: x_out, y_out, z_out
    integer,  intent(in)  :: x_idx_out, y_idx_out, z_idx_out
    real(SP)              :: r_cgs

    r_cgs = u_r*poly_interpolate_3d(xsb, xeb,                 &
                                    ysb, yeb,                 &
                                    zsb, zeb,                 &
                                    xmdn, ymdn, zmdn,         &
                                    rm,                       &
                                    x_out, y_out, z_out,      &
                                    x_idx_out + start_offset, &
                                    y_idx_out + start_offset, &
                                    z_idx_out + start_offset, &
                                    n_interp)

end function get_density


function get_internal_energy(x_out, y_out, z_out, &
                             x_idx_out, y_idx_out, z_idx_out) result(e_cgs)

    use mesh_mod, only : xsb, xeb,   &
                         ysb, yeb,   &
                         zsb, zeb,   &
                         u_e,        &
                         xm, ym, zm, &
                         em

    use tracer_params_mod, only : n_interp, &
                                  start_offset

    use poly_interpolation_mod, only : poly_interpolate_3d

    real(SP), intent(in)  :: x_out, y_out, z_out
    integer,  intent(in)  :: x_idx_out, y_idx_out, z_idx_out
    real(SP)              :: e_cgs

    e_cgs = u_e*poly_interpolate_3d(xsb, xeb,                 &
                                    ysb, yeb,                 &
                                    zsb, zeb,                 &
                                    xm, ym, zm,               &
                                    em,                       &
                                    x_out, y_out, z_out,      &
                                    x_idx_out + start_offset, &
                                    y_idx_out + start_offset, &
                                    z_idx_out + start_offset, &
                                    n_interp)

end function get_internal_energy


function get_gas_pressure(x_out, y_out, z_out, &
                          x_idx_out, y_idx_out, z_idx_out) result(Pg_cgs)

    use mesh_mod, only : n_eos_bins_ee, n_eos_bins_r, &
                         ln_eem_cgs, ln_rm_cgs,       &
                         ln_Pgm_cgs

    use poly_interpolation_mod, only : poly_interpolate_2d

    real(SP), intent(in)  :: x_out, y_out, z_out
    integer,  intent(in)  :: x_idx_out, y_idx_out, z_idx_out
    real(SP)              :: Pg_cgs

    integer,  parameter :: n_interp_eos = 3
    integer,  parameter :: start_offset_eos = 1 - n_interp_eos/2

    real(SP) :: r_cgs, e_cgs
    real(SP) :: ln_ee_cgs, ln_r_cgs
    real(SP) :: ln_Pg_cgs
    integer  :: minloc_arr(1)
    integer  :: ee_idx, r_idx

    r_cgs = get_density(x_out, y_out, z_out, &
                        x_idx_out, y_idx_out, z_idx_out)

    e_cgs = get_internal_energy(x_out, y_out, z_out, &
                                x_idx_out, y_idx_out, z_idx_out)

    ln_ee_cgs = log(e_cgs/r_cgs)
    ln_r_cgs = log(r_cgs)

    minloc_arr = minloc(ln_eem_cgs - ln_ee_cgs, mask=(ln_eem_cgs > ln_ee_cgs))
    ee_idx = minloc_arr(1) - 1

    minloc_arr = minloc(ln_rm_cgs - ln_r_cgs, mask=(ln_rm_cgs > ln_r_cgs))
    r_idx = minloc_arr(1) - 1

    if (ee_idx + start_offset_eos < 1) ee_idx = 1 - start_offset_eos

    ln_Pg_cgs = poly_interpolate_2d(1, n_eos_bins_ee,          &
                                    1, n_eos_bins_r,           &
                                    ln_eem_cgs, ln_rm_cgs,     &
                                    ln_Pgm_cgs,                &
                                    ln_ee_cgs, ln_r_cgs,       &
                                    ee_idx + start_offset_eos, &
                                    r_idx  + start_offset_eos, &
                                    n_interp_eos)

    Pg_cgs = exp(ln_Pg_cgs)

end function get_gas_pressure


function get_magnetic_pressure(x_out, y_out, z_out, &
                               x_idx_out, y_idx_out, z_idx_out) result(PB_cgs)

    use mesh_mod, only : u_B

    real(SP), intent(in)  :: x_out, y_out, z_out
    integer,  intent(in)  :: x_idx_out, y_idx_out, z_idx_out
    real(SP)              :: PB_cgs

    real(SP), parameter :: ONE_OVER_8_PI = 3.978873577e-2

    real(SP) :: Bx, By, Bz
    real(SP) :: B_cgs

    call interpolate_field(x_idx_out, y_idx_out, z_idx_out, &
                           x_out, y_out, z_out,             &
                           Bx, By, Bz)

    B_cgs = u_B*sqrt(Bx*Bx + By*By + Bz*Bz)

    PB_cgs = B_cgs*B_cgs*ONE_OVER_8_PI

end function get_magnetic_pressure


function get_plasma_beta(x_out, y_out, z_out, &
                         x_idx_out, y_idx_out, z_idx_out) result(beta)

    real(SP), intent(in)  :: x_out, y_out, z_out
    integer,  intent(in)  :: x_idx_out, y_idx_out, z_idx_out
    real(SP)              :: beta

    real(SP) :: Pg_cgs, PB_cgs

    Pg_cgs = get_gas_pressure(x_out, y_out, z_out, &
                              x_idx_out, y_idx_out, z_idx_out)

    PB_cgs = get_magnetic_pressure(x_out, y_out, z_out, &
                                   x_idx_out, y_idx_out, z_idx_out)

    beta = Pg_cgs/PB_cgs

end function get_plasma_beta


function evaluate_decoupling_condition(z, beta) result(decouple)

    use tracer_params_mod, only : decoupling_z, &
                                  decoupling_beta

    real(SP), intent(in) :: z, beta
    logical              :: decouple

    decouple = (z > decoupling_z .and. beta > decoupling_beta)

end function evaluate_decoupling_condition

end module tracer_base_mod
