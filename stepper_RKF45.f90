module stepper_mod
implicit none

! The Dormand-Prince scheme, a fifth-order Runge-Kutta method with error
! estimation through an embedded fourth-order step.

private

integer, parameter :: SP = kind(0.0)

! ***************************** Coefficients *****************************

real(SP) :: a2(1) = [     1.0/5.0       ]
real(SP) :: a3(2) = [     3.0/40.0    , &
                          9.0/40.0      ]
real(SP) :: a4(3) = [    44.0/45.0    , &
                        -56.0/15.0    , &
                         32.0/9.0       ]
real(SP) :: a5(4) = [ 19372.0/6561.0  , &
                     -25360.0/2187.0  , &
                      64448.0/6561.0  , &
                       -212.0/729.0     ]
real(SP) :: a6(5) = [  9017.0/3168.0  , &
                       -355.0/33.0    , &
                      46732.0/5247.0  , &
                         49.0/176.0   , &
                      -5103.0/18656.0   ]
real(SP) :: a7(6) = [    35.0/384.0   , &
                          0.0         , &
                        500.0/1113.0  , &
                        125.0/192.0   , &
                      -2187.0/6784.0  , &
                         11.0/84.0      ]

real(SP) ::  e(7) = [    71.0/57600.0 , &
                          0.0         , &
                        -71.0/16695.0 , &
                         71.0/1920.0  , &
                     -17253.0/339200.0, &
                         22.0/525.0   , &
                         -1.0/40.0      ]

real(SP) ::  d(7) = [-12715105075.0/11282082432.0 , &
                                0.0               , &
                      87487479700.0/32700410799.0 , &
                     -10690763975.0/1880347072.0  , &
                     701980252875.0/199316789632.0, &
                      -1453857185.0/822651844.0   , &
                         69997945.0/29380423.0      ]

real(SP) :: cx(5)
real(SP) :: cy(5)
real(SP) :: cz(5)

! ****************************** Variables *******************************

real(SP) :: dx_ds_current(7)
real(SP) :: dy_ds_current(7)
real(SP) :: dz_ds_current(7)

real(SP) :: dx_ds_old(7)
real(SP) :: dy_ds_old(7)
real(SP) :: dz_ds_old(7)

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

    alpha = 1.0/5.0 - 0.75*beta

    dx_ds_old(7) = dx_ds_initial
    dy_ds_old(7) = dy_ds_initial
    dz_ds_old(7) = dz_ds_initial

end subroutine initialize_stepper


subroutine step(x_old, y_old, z_old,                         &
                x_idx_old, y_idx_old, z_idx_old,             &
                ds_current,                                  &
                x_current, y_current, z_current,             &
                x_idx_current, y_idx_current, z_idx_current, &
                terminate)

    use tracer_params_mod, only : direction

    use tracer_base_mod, only : apply_boundary_conditions, &
                                update_position_indices,   &
                                interpolate_field

    real(SP), intent(in) :: x_old, y_old, z_old
    integer,  intent(in) :: x_idx_old, y_idx_old, z_idx_old
    real(SP), intent(in) :: ds_current

    real(SP), intent(out) :: x_current, y_current, z_current
    integer,  intent(out) :: x_idx_current, y_idx_current, z_idx_current
    logical,  intent(out) :: terminate

    real(SP) :: norm
    real(SP) :: Bx, By, Bz

    x_idx_current = x_idx_old
    y_idx_current = y_idx_old
    z_idx_current = z_idx_old

    ! ******************************** Step 1 ********************************

    dx_ds_current(1) = dx_ds_old(7)
    dy_ds_current(1) = dy_ds_old(7)
    dz_ds_current(1) = dz_ds_old(7)

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

    norm = direction/sqrt(Bx*Bx + By*By + Bz*Bz)
    dx_ds_current(2) = Bx*norm
    dy_ds_current(2) = By*norm
    dz_ds_current(2) = Bz*norm

    x_current = x_old + ds_current*(a3(1)*dx_ds_current(1) + &
                                    a3(2)*dx_ds_current(2))

    y_current = y_old + ds_current*(a3(1)*dy_ds_current(1) + &
                                    a3(2)*dy_ds_current(2))

    z_current = z_old + ds_current*(a3(1)*dz_ds_current(1) + &
                                    a3(2)*dz_ds_current(2))

    call apply_boundary_conditions(x_current, y_current, z_current,             &
                                   x_idx_current, y_idx_current, z_idx_current, &
                                   terminate)

    call update_position_indices(x_current, y_current, z_current, &
                                 x_idx_current, y_idx_current, z_idx_current)

    call interpolate_field(x_idx_current, y_idx_current, z_idx_current, &
                           x_current, y_current, z_current,             &
                           Bx, By, Bz)

    ! ******************************** Step 3 ********************************

    norm = direction/sqrt(Bx*Bx + By*By + Bz*Bz)
    dx_ds_current(3) = Bx*norm
    dy_ds_current(3) = By*norm
    dz_ds_current(3) = Bz*norm

    x_current = x_old + ds_current*(a4(1)*dx_ds_current(1) + &
                                    a4(2)*dx_ds_current(2) + &
                                    a4(3)*dx_ds_current(3))

    y_current = y_old + ds_current*(a4(1)*dy_ds_current(1) + &
                                    a4(2)*dy_ds_current(2) + &
                                    a4(3)*dy_ds_current(3))

    z_current = z_old + ds_current*(a4(1)*dz_ds_current(1) + &
                                    a4(2)*dz_ds_current(2) + &
                                    a4(3)*dz_ds_current(3))

    call apply_boundary_conditions(x_current, y_current, z_current,             &
                                   x_idx_current, y_idx_current, z_idx_current, &
                                   terminate)

    call update_position_indices(x_current, y_current, z_current, &
                                 x_idx_current, y_idx_current, z_idx_current)

    call interpolate_field(x_idx_current, y_idx_current, z_idx_current, &
                           x_current, y_current, z_current,             &
                           Bx, By, Bz)

    ! ******************************** Step 4 ********************************

    norm = direction/sqrt(Bx*Bx + By*By + Bz*Bz)
    dx_ds_current(4) = Bx*norm
    dy_ds_current(4) = By*norm
    dz_ds_current(4) = Bz*norm

    x_current = x_old + ds_current*(a5(1)*dx_ds_current(1) + &
                                    a5(2)*dx_ds_current(2) + &
                                    a5(3)*dx_ds_current(3) + &
                                    a5(4)*dx_ds_current(4))

    y_current = y_old + ds_current*(a5(1)*dy_ds_current(1) + &
                                    a5(2)*dy_ds_current(2) + &
                                    a5(3)*dy_ds_current(3) + &
                                    a5(4)*dy_ds_current(4))

    z_current = z_old + ds_current*(a5(1)*dz_ds_current(1) + &
                                    a5(2)*dz_ds_current(2) + &
                                    a5(3)*dz_ds_current(3) + &
                                    a5(4)*dz_ds_current(4))

    call apply_boundary_conditions(x_current, y_current, z_current,             &
                                   x_idx_current, y_idx_current, z_idx_current, &
                                   terminate)

    call update_position_indices(x_current, y_current, z_current, &
                                 x_idx_current, y_idx_current, z_idx_current)

    call interpolate_field(x_idx_current, y_idx_current, z_idx_current, &
                           x_current, y_current, z_current,             &
                           Bx, By, Bz)

    ! ******************************** Step 5 ********************************

    norm = direction/sqrt(Bx*Bx + By*By + Bz*Bz)
    dx_ds_current(5) = Bx*norm
    dy_ds_current(5) = By*norm
    dz_ds_current(5) = Bz*norm

    x_current = x_old + ds_current*(a6(1)*dx_ds_current(1) + &
                                    a6(2)*dx_ds_current(2) + &
                                    a6(3)*dx_ds_current(3) + &
                                    a6(4)*dx_ds_current(4) + &
                                    a6(5)*dx_ds_current(5))

    y_current = y_old + ds_current*(a6(1)*dy_ds_current(1) + &
                                    a6(2)*dy_ds_current(2) + &
                                    a6(3)*dy_ds_current(3) + &
                                    a6(4)*dy_ds_current(4) + &
                                    a6(5)*dy_ds_current(5))

    z_current = z_old + ds_current*(a6(1)*dz_ds_current(1) + &
                                    a6(2)*dz_ds_current(2) + &
                                    a6(3)*dz_ds_current(3) + &
                                    a6(4)*dz_ds_current(4) + &
                                    a6(5)*dz_ds_current(5))

    call apply_boundary_conditions(x_current, y_current, z_current,             &
                                   x_idx_current, y_idx_current, z_idx_current, &
                                   terminate)

    call update_position_indices(x_current, y_current, z_current, &
                                 x_idx_current, y_idx_current, z_idx_current)

    call interpolate_field(x_idx_current, y_idx_current, z_idx_current, &
                           x_current, y_current, z_current,             &
                           Bx, By, Bz)

    ! ******************************** Step 6 ********************************

    norm = direction/sqrt(Bx*Bx + By*By + Bz*Bz)
    dx_ds_current(6) = Bx*norm
    dy_ds_current(6) = By*norm
    dz_ds_current(6) = Bz*norm

    dx = ds_current*(a7(1)*dx_ds_current(1) + &
                     a7(3)*dx_ds_current(3) + &
                     a7(4)*dx_ds_current(4) + &
                     a7(5)*dx_ds_current(5) + &
                     a7(6)*dx_ds_current(6))

    dy = ds_current*(a7(1)*dy_ds_current(1) + &
                     a7(3)*dy_ds_current(3) + &
                     a7(4)*dy_ds_current(4) + &
                     a7(5)*dy_ds_current(5) + &
                     a7(6)*dy_ds_current(6))

    dz = ds_current*(a7(1)*dz_ds_current(1) + &
                     a7(3)*dz_ds_current(3) + &
                     a7(4)*dz_ds_current(4) + &
                     a7(5)*dz_ds_current(5) + &
                     a7(6)*dz_ds_current(6))

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

    norm = direction/sqrt(Bx*Bx + By*By + Bz*Bz)
    dx_ds_current(7) = Bx*norm
    dy_ds_current(7) = By*norm
    dz_ds_current(7) = Bz*norm

    ! **************************** Compute error *****************************

    delta_x = ds_current*(e(1)*dx_ds_current(1) + &
                          e(3)*dx_ds_current(3) + &
                          e(4)*dx_ds_current(4) + &
                          e(5)*dx_ds_current(5) + &
                          e(6)*dx_ds_current(6) + &
                          e(7)*dx_ds_current(7))

    delta_y = ds_current*(e(1)*dy_ds_current(1) + &
                          e(3)*dy_ds_current(3) + &
                          e(4)*dy_ds_current(4) + &
                          e(5)*dy_ds_current(5) + &
                          e(6)*dy_ds_current(6) + &
                          e(7)*dy_ds_current(7))

    delta_z = ds_current*(e(1)*dz_ds_current(1) + &
                          e(3)*dz_ds_current(3) + &
                          e(4)*dz_ds_current(4) + &
                          e(5)*dz_ds_current(5) + &
                          e(6)*dz_ds_current(6) + &
                          e(7)*dz_ds_current(7))

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

    ! Computes coefficients for the fourth-order interpolating polynomial
    ! between s_old and s.

    cx(1) = x_old
    cy(1) = y_old
    cz(1) = z_old

    cx(2) = dx
    cy(2) = dy
    cz(2) = dz

    cx(3) = ds_current*dx_ds_current(1) - cx(2)
    cy(3) = ds_current*dy_ds_current(1) - cy(2)
    cz(3) = ds_current*dz_ds_current(1) - cz(2)

    cx(4) = cx(2) - ds_current*dx_ds_current(7) - cx(3)
    cy(4) = cy(2) - ds_current*dy_ds_current(7) - cy(3)
    cz(4) = cz(2) - ds_current*dz_ds_current(7) - cz(3)

    cx(5) = ds_current*(d(1)*dx_ds_current(1) + &
                        d(3)*dx_ds_current(3) + &
                        d(4)*dx_ds_current(4) + &
                        d(5)*dx_ds_current(5) + &
                        d(6)*dx_ds_current(6) + &
                        d(7)*dx_ds_current(7))

    cy(5) = ds_current*(d(1)*dy_ds_current(1) + &
                        d(3)*dy_ds_current(3) + &
                        d(4)*dy_ds_current(4) + &
                        d(5)*dy_ds_current(5) + &
                        d(6)*dy_ds_current(6) + &
                        d(7)*dy_ds_current(7))

    cz(5) = ds_current*(d(1)*dz_ds_current(1) + &
                        d(3)*dz_ds_current(3) + &
                        d(4)*dz_ds_current(4) + &
                        d(5)*dz_ds_current(5) + &
                        d(6)*dz_ds_current(6) + &
                        d(7)*dz_ds_current(7))

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

    real(SP) :: t, one_minus_t
    logical  :: terminate

    t = (s_out - s_old)/ds_current
    one_minus_t = 1.0 - t

    x_out = cx(1) + t*(cx(2) +                 &
                       one_minus_t*(cx(3) +    &
                                    t*(cx(4) + &
                                       one_minus_t*cx(5))))

    y_out = cy(1) + t*(cy(2) +                 &
                       one_minus_t*(cy(3) +    &
                                    t*(cy(4) + &
                                       one_minus_t*cy(5))))

    z_out = cz(1) + t*(cz(2) +                 &
                       one_minus_t*(cz(3) +    &
                                    t*(cz(4) + &
                                       one_minus_t*cz(5))))

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
