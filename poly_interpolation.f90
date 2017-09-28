module poly_interpolation_mod
implicit none

private

integer, parameter :: SP = kind(0.0)

public :: poly_interpolate_2d, poly_interpolate_3d

contains

function poly_interpolate_2d(xsb, xeb,                 &
                             ysb, yeb,                 &
                             xm, ym,                   &
                             fm,                       &
                             x, y,                     &
                             x_start_idx, y_start_idx, &
                             n_interp) result(P_xy)

    integer,  intent(in) :: xsb, xeb
    integer,  intent(in) :: ysb, yeb
    real(SP), intent(in) :: xm(xsb:xeb), ym(ysb:yeb)
    real(SP), intent(in) :: fm(xsb:xeb, ysb:yeb)
    real(SP), intent(in) :: x, y
    integer,  intent(in) :: x_start_idx, y_start_idx
    integer,  intent(in) :: n_interp
    real(SP)             :: P_xy

    real(SP) :: x_interp(n_interp)
    real(SP) :: y_interp(n_interp)
    real(SP) :: f_interp(n_interp, n_interp)

    real(SP) :: P_x(n_interp)

    real(SP) :: C(n_interp)
    real(SP) :: D(n_interp)

    integer :: x_end_idx
    integer :: y_end_idx
    integer :: order

    real(SP) :: correction_factor
    integer :: n, i, j

    order = n_interp - 1

    x_end_idx = x_start_idx + order
    y_end_idx = y_start_idx + order

    ! ************************ Validate index ranges *************************

    if (x_start_idx < xsb .or. x_end_idx > xeb) then
        write(*, '(A, A, I0, A, I0, A, I0, A, I0, A)') &
            'Warning in poly_interpolate: x-index range outside limits ', &
            '(x_start_idx = ', x_start_idx, ', x_end_idx = ', x_end_idx, &
            ', limits = [', xsb, ', ', xeb, '])'
    end if

    if (y_start_idx < ysb .or. y_end_idx > yeb) then
        write(*, '(A, A, I0, A, I0, A, I0, A, I0, A)') &
            'Warning in poly_interpolate: y-index range outside limits ', &
            '(y_start_idx = ', y_start_idx, ', y_end_idx = ', y_end_idx, &
            ', limits = [', ysb, ', ', yeb, '])'
    end if

    ! **************** Extract relevant coordinates and values ***************

    x_interp = xm(x_start_idx:x_end_idx)
    y_interp = ym(y_start_idx:y_end_idx)
    f_interp = fm(x_start_idx:x_end_idx, y_start_idx:y_end_idx)

    ! ********************* Validate interpolation point *********************

    if (x < x_interp(1) .or. x > x_interp(n_interp)) then
        write(*, '(A, A, G0.3, A, G0.3, A, G0.3, A)') &
            'Warning in poly_interpolate: x outside interpolation range ', &
            '(x = ', x, ', range = [', x_interp(1), ', ', x_interp(n_interp), '])'
    end if

    if (y < y_interp(1) .or. y > y_interp(n_interp)) then
        write(*, '(A, A, G0.3, A, G0.3, A, G0.3, A)') &
            'Warning in poly_interpolate: y outside interpolation range ', &
            '(y = ', y, ', range = [', y_interp(1), ', ', y_interp(n_interp), '])'
    end if

    ! ********************** Interpolate in x-direction **********************

    do j = 1, n_interp

        C = f_interp(:, j)
        D = C

        P_x(j) = C(1)

        do n = 1, order

            do i = 1, n_interp-n

                ! Update corrections
                correction_factor = (C(i+1) - D(i))/(x_interp(i+n) - x_interp(i))
                C(i) = (x - x_interp(i))*correction_factor
                D(i) = (x - x_interp(i+n))*correction_factor

            end do

            ! Correct the current approximation
            P_x(j) = P_x(j) + C(1)

        end do

    end do

    ! ********************** Interpolate in y-direction **********************

    C = P_x
    D = C

    P_xy = C(1)

    do n = 1, order

        do j = 1, n_interp-n

            ! Update corrections
            correction_factor = (C(j+1) - D(j))/(y_interp(j+n) - y_interp(j))
            C(j) = (y - y_interp(j))*correction_factor
            D(j) = (y - y_interp(j+n))*correction_factor

        end do

        ! Correct the current approximation
        P_xy = P_xy + C(1)

    end do

end function poly_interpolate_2d

function poly_interpolate_3d(xsb, xeb,                              &
                             ysb, yeb,                              &
                             zsb, zeb,                              &
                             xm, ym, zm,                            &
                             fm,                                    &
                             x, y, z,                               &
                             x_start_idx, y_start_idx, z_start_idx, &
                             n_interp) result(P_xyz)

    integer,  intent(in) :: xsb, xeb
    integer,  intent(in) :: ysb, yeb
    integer,  intent(in) :: zsb, zeb

    real(SP), intent(in) :: xm(xsb:xeb)
    real(SP), intent(in) :: ym(ysb:yeb)
    real(SP), intent(in) :: zm(zsb:zeb)
    real(SP), intent(in) :: fm(xsb:xeb, ysb:yeb, zsb:zeb)

    real(SP), intent(in) :: x, y, z
    integer,  intent(in) :: x_start_idx, y_start_idx, z_start_idx

    integer,  intent(in) :: n_interp

    real(SP)             :: P_xyz

    real(SP) :: x_interp(n_interp)
    real(SP) :: y_interp(n_interp)
    real(SP) :: z_interp(n_interp)
    real(SP) :: f_interp(n_interp, n_interp, n_interp)

    real(SP) :: P_x(n_interp, n_interp)
    real(SP) :: P_xy(n_interp)

    real(SP) :: C(n_interp)
    real(SP) :: D(n_interp)

    integer  :: order

    integer  :: x_end_idx, y_end_idx, z_end_idx

    real(SP) :: correction_factor
    integer  :: n, i, j, k

    order = n_interp - 1

    x_end_idx = x_start_idx + order
    y_end_idx = y_start_idx + order
    z_end_idx = z_start_idx + order

    ! ************************ Validate index ranges *************************

    if (x_start_idx < xsb .or. x_end_idx > xeb) then
        write(*, '(A, A, I0, A, I0, A, I0, A, I0, A)') &
            'Warning in poly_interpolate: x-index range outside limits ', &
            '(x_start_idx = ', x_start_idx, ', x_end_idx = ', x_end_idx, &
            ', limits = [', xsb, ', ', xeb, '])'
    end if

    if (y_start_idx < ysb .or. y_end_idx > yeb) then
        write(*, '(A, A, I0, A, I0, A, I0, A, I0, A)') &
            'Warning in poly_interpolate: y-index range outside limits ', &
            '(y_start_idx = ', y_start_idx, ', y_end_idx = ', y_end_idx, &
            ', limits = [', ysb, ', ', yeb, '])'
    end if

    if (z_start_idx < zsb .or. z_end_idx > zeb) then
        write(*, '(A, A, I0, A, I0, A, I0, A, I0, A)') &
            'Warning in poly_interpolate: z-index range outside limits ', &
            '(z_start_idx = ', z_start_idx, ', z_end_idx = ', z_end_idx, &
            ', limits = [', zsb, ', ', zeb, '])'
    end if

    ! **************** Extract relevant coordinates and values ***************

    x_interp = xm(x_start_idx:x_end_idx)
    y_interp = ym(y_start_idx:y_end_idx)
    z_interp = zm(z_start_idx:z_end_idx)
    f_interp = fm(x_start_idx:x_end_idx, y_start_idx:y_end_idx, z_start_idx:z_end_idx)

    ! ********************* Validate interpolation point *********************

    if (x < x_interp(1) .or. x > x_interp(n_interp)) then
        write(*, '(A, A, G0.3, A, G0.3, A, G0.3, A)') &
            'Warning in poly_interpolate: x outside interpolation range ', &
            '(x = ', x, ', range = [', x_interp(1), ', ', x_interp(n_interp), '])'
    end if

    if (y < y_interp(1) .or. y > y_interp(n_interp)) then
        write(*, '(A, A, G0.3, A, G0.3, A, G0.3, A)') &
            'Warning in poly_interpolate: y outside interpolation range ', &
            '(y = ', y, ', range = [', y_interp(1), ', ', y_interp(n_interp), '])'
    end if

    if (z < z_interp(1) .or. z > z_interp(n_interp)) then
        write(*, '(A, A, G0.3, A, G0.3, A, G0.3, A)') &
            'Warning in poly_interpolate: z outside interpolation range ', &
            '(z = ', z, ', range = [', z_interp(1), ', ', z_interp(n_interp), '])'
    end if

    ! ********************** Interpolate in x-direction **********************

    do k = 1, n_interp
        do j = 1, n_interp

            C = f_interp(:, j, k)
            D = C

            P_x(j, k) = C(1)

            do n = 1, order

                do i = 1, n_interp-n

                    ! Update corrections
                    correction_factor = (C(i+1) - D(i))/(x_interp(i+n) - x_interp(i))
                    C(i) = (x - x_interp(i))*correction_factor
                    D(i) = (x - x_interp(i+n))*correction_factor

                end do

                ! Correct the current approximation
                P_x(j, k) = P_x(j, k) + C(1)

            end do

        end do
    end do

    ! ********************** Interpolate in y-direction **********************

    do k = 1, n_interp

        C = P_x(:, k)
        D = C

        P_xy(k) = C(1)

        do n = 1, order

            do j = 1, n_interp-n

                ! Update corrections
                correction_factor = (C(j+1) - D(j))/(y_interp(j+n) - y_interp(j))
                C(j) = (y - y_interp(j))*correction_factor
                D(j) = (y - y_interp(j+n))*correction_factor

            end do

            ! Correct the current approximation
            P_xy(k) = P_xy(k) + C(1)

        end do

    end do

    ! ********************** Interpolate in z-direction **********************

    C = P_xy
    D = C

    P_xyz = C(1)

    do n = 1, order

        do k = 1, n_interp-n

            ! Update corrections
            correction_factor = (C(k+1) - D(k))/(z_interp(k+n) - z_interp(k))
            C(k) = (z - z_interp(k))*correction_factor
            D(k) = (z - z_interp(k+n))*correction_factor

        end do

        ! Correct the current approximation
        P_xyz = P_xyz + C(1)

    end do

end function poly_interpolate_3d

end module poly_interpolation_mod
