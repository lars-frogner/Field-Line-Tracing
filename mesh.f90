module mesh_mod
implicit none

private

integer, parameter :: SP = kind(0.0)

integer,               public :: mb
integer,               public :: mx, my, mz
logical,               public :: periodic_x, periodic_y, periodic_z
integer,               public :: xs, ys, zs
integer,               public :: xe, ye, ze
integer,               public :: xsb, ysb, zsb
integer,               public :: xeb, yeb, zeb

real(SP), allocatable, public :: xm(:), ym(:), zm(:)
real(SP), allocatable, public :: xmdn(:), ymdn(:), zmdn(:)
real(SP), allocatable, public :: dxidxup(:), dyidyup(:), dzidzup(:)
real(SP), allocatable, public :: dxidxdn(:), dyidydn(:), dzidzdn(:)
real(SP), allocatable, public :: bxm(:, :, :), bym(:, :, :), bzm(:, :, :)
real(SP), allocatable, public :: rm(:, :, :), em(:, :, :)

real(SP),              public :: u_l, u_B, u_r, u_e, u_ee

integer,               public :: n_eos_bins_ee, n_eos_bins_r

real(SP), allocatable, public :: ln_eem_cgs(:), ln_rm_cgs(:)
real(SP), allocatable, public :: ln_Pgm_cgs(:, :)

logical,               public :: mesh_initialized = .false.

real(SP) :: eos_ee_cgs_min, eos_ee_cgs_max, eos_r_cgs_min, eos_r_cgs_max

public :: initialize_mesh, &
          free_mesh

contains

subroutine initialize_mesh(idl_file_path)

    character(len=*), intent(in) :: idl_file_path

    character(len=120) :: path = ' '
    character(len=60)  :: idl_file = ' ', eos_file = ' ', mesh_file = ' ', snap_file = ' '
    integer :: pos, len

    len = len_trim(idl_file_path)
    pos = index(idl_file_path, '/', .true.)

    if (pos /= 0) then
        path(1:pos) = idl_file_path(1:pos)
        idl_file(1:len-pos) = idl_file_path(pos+1:len)
    end if

    call read_params(path, idl_file, eos_file, mesh_file, snap_file)
    call read_meshgrid(trim(path)//trim(mesh_file))
    call read_snap(trim(path)//trim(snap_file), 'bx')
    call read_snap(trim(path)//trim(snap_file), 'by')
    call read_snap(trim(path)//trim(snap_file), 'bz')
    call set_boundaries()

    mesh_initialized = .true.

end subroutine initialize_mesh


subroutine read_params(path, idl_file, eos_file, mesh_file, snap_file)

    character(len=*),  intent(in)  :: path, idl_file
    character(len=60), intent(out) :: eos_file, mesh_file, snap_file

    integer :: unit, iostat
    character(len=70) :: buffer
    character(len=15) :: param
    character(len=1) :: eq
    character(len=3) :: snap_str
    integer :: pos, temp, isnap

    open(newunit=unit, file=trim(path)//trim(idl_file), status='old', action='read', iostat=iostat)
    if (iostat /= 0) then
        write(*, '(A, A)') 'Error: could not open ', trim(path)//trim(idl_file)
        stop
    end if

    do
        read(unit, '(A)', iostat=iostat) buffer
        if (iostat /= 0) exit

        buffer = adjustl(buffer)

        if (len_trim(buffer) > 0 .and. buffer(1:1) /= ';') then

            read(buffer, *) param

            select case (param)

                case ('mx')
                    read(buffer, *) param, eq, mx
                case ('my')
                    read(buffer, *) param, eq, my
                case ('mz')
                    read(buffer, *) param, eq, mz
                case ('mb')
                    read(buffer, *) param, eq, mb
                case ('periodic_x')
                    read(buffer, *) param, eq, temp
                    periodic_x = (temp == 1)
                case ('periodic_y')
                    read(buffer, *) param, eq, temp
                    periodic_y = (temp == 1)
                case ('periodic_z')
                    read(buffer, *) param, eq, temp
                    periodic_z = (temp == 1)
                case ('meshfile')
                    read(buffer, *) param, eq, mesh_file
                    mesh_file = adjustl(mesh_file)
                    pos = index(mesh_file, '"')
                    mesh_file(1:pos-2) = mesh_file(2:pos-1)
                case ('snapname')
                    read(buffer, *) param, eq, snap_file
                    snap_file = adjustl(snap_file)
                    pos = index(snap_file, '"')
                    snap_file(1:pos-2) = snap_file(2:pos-1)
                case ('isnap')
                    read(buffer, *) param, eq, isnap
                case ('u_l')
                    read(buffer, *) param, eq, u_l
                case ('u_B')
                    read(buffer, *) param, eq, u_B
                case ('u_r')
                    read(buffer, *) param, eq, u_r
                case ('u_e')
                    read(buffer, *) param, eq, u_e
                case ('u_ee')
                    read(buffer, *) param, eq, u_ee

            end select

        end if

    end do

    close(unit)

    write(snap_str, '(I3.3)') isnap
    pos = index(snap_file, ' ')
    snap_file(pos:pos+8) = '_'//snap_str//'.snap'

    open(newunit=unit, file=trim(path)//'tabparam.in', status='old', action='read', iostat=iostat)
    if (iostat /= 0) then
        write(*, '(A, A)') 'Error: could not open ', trim(path)//'tabparam.in'
        stop
    end if

    do
        read(unit, '(A)', iostat=iostat) buffer
        if (iostat /= 0) exit

        buffer = adjustl(buffer)

        read(buffer, *) param

        select case (param)

            case ('EOSTableFile')
                read(buffer, *) param, eq, eos_file
                eos_file = adjustl(eos_file)
                pos = index(eos_file, "'")
                eos_file(1:pos-2) = eos_file(2:pos-1)
            case ('nRhoBin')
                read(buffer, *) param, eq, n_eos_bins_r
            case ('nEiBin')
                read(buffer, *) param, eq, n_eos_bins_ee
            case ('RhoMin')
                read(buffer, *) param, eq, eos_r_cgs_min
            case ('RhoMax')
                read(buffer, *) param, eq, eos_r_cgs_max
            case ('EiMin')
                read(buffer, *) param, eq, eos_ee_cgs_min
            case ('EiMax')
                read(buffer, *) param, eq, eos_ee_cgs_max

        end select

    end do

    close(unit)

end subroutine read_params


subroutine read_eos_table(eos_file_path)

    character(len=*), intent(in) :: eos_file_path

    integer  :: unit, iostat, i
    real(SP) :: d_ln_ee_cgs, d_ln_r_cgs

    open(newunit=unit, file=eos_file_path, status='old', form='unformatted', access='stream', iostat=iostat)
    if (iostat /= 0) then
        write(*, '(A, A)') 'Error: could not open ', eos_file_path
        stop
    end if

    if (allocated(ln_Pgm_cgs)) deallocate(ln_Pgm_cgs)
    allocate(ln_Pgm_cgs(n_eos_bins_ee, n_eos_bins_r))

    read(unit) ln_Pgm_cgs

    close(unit)

    d_ln_ee_cgs = (log(eos_ee_cgs_max) - log(eos_ee_cgs_min))/real(n_eos_bins_ee - 1, SP)
    d_ln_r_cgs  = (log(eos_r_cgs_max)  - log(eos_r_cgs_min)) /real(n_eos_bins_r  - 1, SP)

    if (allocated(ln_eem_cgs)) deallocate(ln_eem_cgs)
    if (allocated(ln_rm_cgs)) deallocate(ln_rm_cgs)
    allocate(ln_eem_cgs(n_eos_bins_ee), ln_rm_cgs(n_eos_bins_r))

    ln_eem_cgs = [(log(eos_ee_cgs_min) + d_ln_ee_cgs*(i-1), i = 1, n_eos_bins_ee)]
    ln_rm_cgs  = [(log(eos_r_cgs_min)  + d_ln_r_cgs *(i-1), i = 1, n_eos_bins_r)]

end subroutine read_eos_table


subroutine read_meshgrid(mesh_file_path)

    character(len=*), intent(in) :: mesh_file_path
    integer  :: unit, iostat

    if (allocated(xm))      deallocate(xm)
    if (allocated(xmdn))    deallocate(xmdn)
    if (allocated(dxidxup)) deallocate(dxidxup)
    if (allocated(dxidxdn)) deallocate(dxidxdn)

    if (allocated(ym))      deallocate(ym)
    if (allocated(ymdn))    deallocate(ymdn)
    if (allocated(dyidyup)) deallocate(dyidyup)
    if (allocated(dyidydn)) deallocate(dyidydn)

    if (allocated(zm))      deallocate(zm)
    if (allocated(zmdn))    deallocate(zmdn)
    if (allocated(dzidzup)) deallocate(dzidzup)
    if (allocated(dzidzdn)) deallocate(dzidzdn)

    open(newunit=unit, file=mesh_file_path, status='old', action='read', iostat=iostat)
    if (iostat /= 0) then
        write(*, '(A, A)') 'Error: could not open ', mesh_file_path
        stop
    end if

    read(unit, *) mx

    xs  = 1
    xe  = mx
    xsb = xs - mb
    xeb = xe + mb

    allocate(xm(xsb:xeb),      &
             xmdn(xsb:xeb),    &
             dxidxup(xsb:xeb), &
             dxidxdn(xsb:xeb))

    read(unit, *) xm(xs:xe)
    read(unit, *) xmdn(xs:xe)
    read(unit, *) dxidxup(xs:xe)
    read(unit, *) dxidxdn(xs:xe)

    read(unit, *) my

    ys  = 1
    ye  = my
    ysb = ys - mb
    yeb = ye + mb

    allocate(ym(ysb:yeb),      &
             ymdn(ysb:yeb),    &
             dyidyup(ysb:yeb), &
             dyidydn(ysb:yeb))

    read(unit, *) ym(ys:ye)
    read(unit, *) ymdn(ys:ye)
    read(unit, *) dyidyup(ys:ye)
    read(unit, *) dyidydn(ys:ye)

    read(unit, *) mz

    zs  = 1
    ze  = mz
    zsb = zs - mb
    zeb = ze + mb

    allocate(zm(zsb:zeb),      &
             zmdn(zsb:zeb),    &
             dzidzup(zsb:zeb), &
             dzidzdn(zsb:zeb))

    read(unit, *) zm(zs:ze)
    read(unit, *) zmdn(zs:ze)
    read(unit, *) dzidzup(zs:ze)
    read(unit, *) dzidzdn(zs:ze)

    close(unit)

end subroutine read_meshgrid


subroutine read_snap(snap_file_path, var_name)

    character(len=*), intent(in) :: snap_file_path
    character(len=*), intent(in) :: var_name

    integer, parameter :: dsize = 4

    integer :: unit, iostat
    integer(8) :: var_size

    var_size = mx*my*mz*dsize

    open(newunit=unit, file=snap_file_path, status='old', form='unformatted', access='stream', iostat=iostat)
    if (iostat /= 0) then
        write(*, '(A, A)') 'Error: could not open ', snap_file_path
        stop
    end if

    select case (trim(var_name))

        case ('rm')

            write(*, *) 'rm'

            if (allocated(rm)) deallocate(rm)
            allocate(rm(xsb:xeb, ysb:yeb, zsb:zeb))
            read(unit, pos=1) rm(xs:xe, ys:ye, zs:ze)

        case ('em')

            write(*, *) 'em'

            if (allocated(em)) deallocate(em)
            allocate(em(xsb:xeb, ysb:yeb, zsb:zeb))
            read(unit, pos=(1 + var_size*4)) em(xs:xe, ys:ye, zs:ze)

        case ('bx')

            write(*, *) 'bx'

            if (allocated(bxm)) deallocate(bxm)
            allocate(bxm(xsb:xeb, ysb:yeb, zsb:zeb))
            read(unit, pos=(1 + var_size*5)) bxm(xs:xe, ys:ye, zs:ze)

        case ('by')

            write(*, *) 'by'

            if (allocated(bym)) deallocate(bym)
            allocate(bym(xsb:xeb, ysb:yeb, zsb:zeb))
            read(unit, pos=(1 + var_size*6)) bym(xs:xe, ys:ye, zs:ze)

        case ('bz')

            write(*, *) 'bz'

            if (allocated(bzm)) deallocate(bzm)
            allocate(bzm(xsb:xeb, ysb:yeb, zsb:zeb))
            read(unit, pos=(1 + var_size*7)) bzm(xs:xe, ys:ye, zs:ze)

    end select

    close(unit)

end subroutine read_snap


subroutine set_boundaries()

    real(SP) :: x_range, y_range, z_range
    real(SP) :: dx, dy, dz
    integer  :: i, j, k

    if (periodic_x) then

        x_range = xm(xe) - xm(xs)

        xm(xsb:xs-1) = xm(xe-mb:xe-1) - x_range
        xm(xe+1:xeb) = xm(xs+1:xs+mb) + x_range

        xmdn(xsb:xs-1) = xmdn(xe-mb:xe-1) - x_range
        xmdn(xe+1:xeb) = xmdn(xs+1:xs+mb) + x_range

        if (allocated(bxm)) then
            bxm(xsb:xs-1, :, :) = bxm(xe-mb:xe-1, :, :)
            bxm(xe+1:xeb, :, :) = bxm(xs+1:xs+mb, :, :)
        end if
        if (allocated(bym)) then
            bym(xsb:xs-1, :, :) = bym(xe-mb:xe-1, :, :)
            bym(xe+1:xeb, :, :) = bym(xs+1:xs+mb, :, :)
        end if
        if (allocated(bzm)) then
            bzm(xsb:xs-1, :, :) = bzm(xe-mb:xe-1, :, :)
            bzm(xe+1:xeb, :, :) = bzm(xs+1:xs+mb, :, :)
        end if
        if (allocated(rm)) then
            rm(xsb:xs-1, :, :) = rm(xe-mb:xe-1, :, :)
            rm(xe+1:xeb, :, :) = rm(xs+1:xs+mb, :, :)
        end if
        if (allocated(em)) then
            em(xsb:xs-1, :, :) = em(xe-mb:xe-1, :, :)
            em(xe+1:xeb, :, :) = em(xs+1:xs+mb, :, :)
        end if

    else

        dx = xm(xs+1) - xm(xs)
        xm(xsb:xs-1) = [(xm(xs) + (i-mb-1)*dx, i = 1, mb)]
        dx = xm(xe) - xm(xe-1)
        xm(xe+1:xeb) = [(xm(xe) +        i*dx, i = 1, mb)]

        dx = xmdn(xs+1) - xmdn(xs)
        xmdn(xsb:xs-1) = [(xmdn(xs) + (i-mb-1)*dx, i = 1, mb)]
        dx = xmdn(xe) - xmdn(xe-1)
        xmdn(xe+1:xeb) = [(xmdn(xe) +        i*dx, i = 1, mb)]

        if (allocated(bxm)) then
            bxm(xsb:xs-1, :, :) = 0.0
            bxm(xe+1:xeb, :, :) = 0.0
        end if
        if (allocated(bym)) then
            bym(xsb:xs-1, :, :) = 0.0
            bym(xe+1:xeb, :, :) = 0.0
        end if
        if (allocated(bzm)) then
            bzm(xsb:xs-1, :, :) = 0.0
            bzm(xe+1:xeb, :, :) = 0.0
        end if
        if (allocated(rm)) then
            rm(xsb:xs-1, :, :) = 0.0
            rm(xe+1:xeb, :, :) = 0.0
        end if
        if (allocated(em)) then
            em(xsb:xs-1, :, :) = 0.0
            em(xe+1:xeb, :, :) = 0.0
        end if

    end if

    dxidxup(xsb:xs-1) = dxidxup(xs)
    dxidxup(xe+1:xeb) = dxidxup(xe)

    dxidxdn(xsb:xs-1) = dxidxdn(xs)
    dxidxdn(xe+1:xeb) = dxidxdn(xe)

    if (periodic_y) then

        y_range = ym(ye) - ym(ys)

        ym(ysb:ys-1) = ym(ye-mb:ye-1) - y_range
        ym(ye+1:yeb) = ym(ys+1:ys+mb) + y_range

        ymdn(ysb:ys-1) = ymdn(ye-mb:ye-1) - y_range
        ymdn(ye+1:yeb) = ymdn(ys+1:ys+mb) + y_range

        if (allocated(bxm)) then
            bxm(:, ysb:ys-1, :) = bxm(:, ye-mb:ye-1, :)
            bxm(:, ye+1:yeb, :) = bxm(:, ys+1:ys+mb, :)
        end if
        if (allocated(bym)) then
            bym(:, ysb:ys-1, :) = bym(:, ye-mb:ye-1, :)
            bym(:, ye+1:yeb, :) = bym(:, ys+1:ys+mb, :)
        end if
        if (allocated(bzm)) then
            bzm(:, ysb:ys-1, :) = bzm(:, ye-mb:ye-1, :)
            bzm(:, ye+1:yeb, :) = bzm(:, ys+1:ys+mb, :)
        end if
        if (allocated(rm)) then
            rm(:, ysb:ys-1, :) = rm(:, ye-mb:ye-1, :)
            rm(:, ye+1:yeb, :) = rm(:, ys+1:ys+mb, :)
        end if
        if (allocated(em)) then
            em(:, ysb:ys-1, :) = em(:, ye-mb:ye-1, :)
            em(:, ye+1:yeb, :) = em(:, ys+1:ys+mb, :)
        end if

    else

        dy = ym(ys+1) - ym(ys)
        ym(ysb:ys-1) = [(ym(ys) + (j-mb-1)*dy, j = 1, mb)]
        dy = ym(ye) - ym(ye-1)
        ym(ye+1:yeb) = [(ym(ye) +        j*dy, j = 1, mb)]

        dy = ymdn(ys+1) - ymdn(ys)
        ymdn(ysb:ys-1) = [(ymdn(ys) + (j-mb-1)*dy, j = 1, mb)]
        dy = ymdn(ye) - ymdn(ye-1)
        ymdn(ye+1:yeb) = [(ymdn(ye) +        j*dy, j = 1, mb)]

        if (allocated(bxm)) then
            bxm(:, ysb:ys-1, :) = 0.0
            bxm(:, ye+1:yeb, :) = 0.0
        end if
        if (allocated(bym)) then
            bym(:, ysb:ys-1, :) = 0.0
            bym(:, ye+1:yeb, :) = 0.0
        end if
        if (allocated(bzm)) then
            bzm(:, ysb:ys-1, :) = 0.0
            bzm(:, ye+1:yeb, :) = 0.0
        end if
        if (allocated(rm)) then
            rm(:, ysb:ys-1, :) = 0.0
            rm(:, ye+1:yeb, :) = 0.0
        end if
        if (allocated(em)) then
            em(:, ysb:ys-1, :) = 0.0
            em(:, ye+1:yeb, :) = 0.0
        end if

    end if

    dyidyup(ysb:ys-1) = dyidyup(ys)
    dyidyup(ye+1:yeb) = dyidyup(ye)

    dyidydn(ysb:ys-1) = dyidydn(ys)
    dyidydn(ye+1:yeb) = dyidydn(ye)

    if (periodic_z) then

        z_range = zm(ze) - zm(zs)

        zm(zsb:zs-1) = zm(ze-mb:ze-1) - z_range
        zm(ze+1:zeb) = zm(zs+1:zs+mb) + z_range

        zmdn(zsb:zs-1) = zmdn(ze-mb:ze-1) - z_range
        zmdn(ze+1:zeb) = zmdn(zs+1:zs+mb) + z_range

        if (allocated(bxm)) then
            bxm(:, :, zsb:zs-1) = bxm(:, :, ze-mb:ze-1)
            bxm(:, :, ze+1:zeb) = bxm(:, :, zs+1:zs+mb)
        end if
        if (allocated(bym)) then
            bym(:, :, zsb:zs-1) = bym(:, :, ze-mb:ze-1)
            bym(:, :, ze+1:zeb) = bym(:, :, zs+1:zs+mb)
        end if
        if (allocated(bzm)) then
            bzm(:, :, zsb:zs-1) = bzm(:, :, ze-mb:ze-1)
            bzm(:, :, ze+1:zeb) = bzm(:, :, zs+1:zs+mb)
        end if
        if (allocated(rm)) then
            rm(:, :, zsb:zs-1) = rm(:, :, ze-mb:ze-1)
            rm(:, :, ze+1:zeb) = rm(:, :, zs+1:zs+mb)
        end if
        if (allocated(em)) then
            em(:, :, zsb:zs-1) = em(:, :, ze-mb:ze-1)
            em(:, :, ze+1:zeb) = em(:, :, zs+1:zs+mb)
        end if

    else

        dz = zm(zs+1) - zm(zs)
        zm(zsb:zs-1) = [(zm(zs) + (k-mb-1)*dz, k = 1, mb)]
        dz = zm(ze) - zm(ze-1)
        zm(ze+1:zeb) = [(zm(ze) +        k*dz, k = 1, mb)]

        dz = zmdn(zs+1) - zmdn(zs)
        zmdn(zsb:zs-1) = [(zmdn(zs) + (k-mb-1)*dz, k = 1, mb)]
        dz = zmdn(ze) - zmdn(ze-1)
        zmdn(ze+1:zeb) = [(zmdn(ze) +        k*dz, k = 1, mb)]

        if (allocated(bxm)) then
            bxm(:, :, zsb:zs-1) = 0.0
            bxm(:, :, ze+1:zeb) = 0.0
        end if
        if (allocated(bym)) then
            bym(:, :, zsb:zs-1) = 0.0
            bym(:, :, ze+1:zeb) = 0.0
        end if
        if (allocated(bzm)) then
            bzm(:, :, zsb:zs-1) = 0.0
            bzm(:, :, ze+1:zeb) = 0.0
        end if
        if (allocated(rm)) then
            rm(:, :, zsb:zs-1) = 0.0
            rm(:, :, ze+1:zeb) = 0.0
        end if
        if (allocated(em)) then
            em(:, :, zsb:zs-1) = 0.0
            em(:, :, ze+1:zeb) = 0.0
        end if

    end if

    dzidzup(zsb:zs-1) = dzidzup(zs)
    dzidzup(ze+1:zeb) = dzidzup(ze)

    dzidzdn(zsb:zs-1) = dzidzdn(zs)
    dzidzdn(ze+1:zeb) = dzidzdn(ze)

end subroutine set_boundaries


subroutine free_mesh()

    if (allocated(xm))         deallocate(xm)
    if (allocated(xmdn))       deallocate(xmdn)
    if (allocated(dxidxup))    deallocate(dxidxup)
    if (allocated(dxidxdn))    deallocate(dxidxdn)
    if (allocated(bxm))        deallocate(bxm)

    if (allocated(ym))         deallocate(ym)
    if (allocated(ymdn))       deallocate(ymdn)
    if (allocated(dyidyup))    deallocate(dyidyup)
    if (allocated(dyidydn))    deallocate(dyidydn)
    if (allocated(bym))        deallocate(bym)

    if (allocated(zm))         deallocate(zm)
    if (allocated(zmdn))       deallocate(zmdn)
    if (allocated(dzidzup))    deallocate(dzidzup)
    if (allocated(dzidzdn))    deallocate(dzidzdn)
    if (allocated(bzm))        deallocate(bzm)

    if (allocated(rm))         deallocate(rm)
    if (allocated(em))         deallocate(em)

    if (allocated(ln_eem_cgs)) deallocate(ln_eem_cgs)
    if (allocated(ln_rm_cgs))  deallocate(ln_rm_cgs)
    if (allocated(ln_Pgm_cgs)) deallocate(ln_Pgm_cgs)

    mesh_initialized = .false.

end subroutine free_mesh


end module mesh_mod
