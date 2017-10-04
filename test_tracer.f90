program test_tracer

use tracer_params_mod, only : direction,              &
                              nonuniform_output,      &
                              constant_stepsize,      &
                              interp_order,           &
                              interp_bias,            &
                              abs_tolerance,          &
                              rel_tolerance,          &
                              ds_out,                 &
                              ds_initial,             &
                              error_initial,          &
                              beta,                   &
                              safety_factor,          &
                              scale_min,              &
                              scale_max,              &
                              decoupling_z,           &
                              decoupling_beta,        &
                              decoupling_rate,        &
                              slope_correction,       &
                              max_output_points,      &
                              n_aux,                  &
                              aux_names

use mesh_mod, only : initialize_mesh, &
                     free_mesh,       &
                     u_l,             &
                     xs, xe,          &
                     ys, ye,          &
                     zs, ze,          &
                     xm, ym, zm

use tracer_mod, only : trace

implicit none

integer,  parameter :: SP = kind(0.0)
integer,  parameter :: DP = kind(0.0d0)

integer               :: n_output_points
real(SP), allocatable :: x_out(:), y_out(:), z_out(:)
real(SP), allocatable :: aux_out(:, :)
real(SP)              :: fieldline_length
integer               :: decoupling_index

real(SP), allocatable :: x_initials(:), y_initials(:), z_initials(:)

real(SP) :: directions

character(len=23) :: command_arg
integer :: unit
integer :: i, j
integer :: n_args, n_params, n_initials

integer(8) :: start_count
integer(8) :: end_count
integer(8) :: count_rate
real(DP)   :: elapsed_time

write(*, '(A/)') '*************************** test_tracer.f90 ****************************'
write(*, '(A)') 'Parameters:'

n_args = command_argument_count()

i = 1
call get_command_argument(i, command_arg)
read(command_arg, '(G23.0)')   directions
write(*, '(A, G0.3)') 'directions = ', directions
i = i + 1
call get_command_argument(i, command_arg)
read(command_arg, '(L23)')   nonuniform_output
write(*, '(A, L1)') 'nonuniform_output = ', nonuniform_output
i = i + 1
call get_command_argument(i, command_arg)
read(command_arg, '(L23)')   constant_stepsize
write(*, '(A, L1)') 'constant_stepsize = ', constant_stepsize
i = i + 1
call get_command_argument(i, command_arg)
read(command_arg, '(I23)')   interp_order
write(*, '(A, I0)') 'interp_order = ', interp_order
i = i + 1
call get_command_argument(i, command_arg)
read(command_arg, '(I23)')   interp_bias
write(*, '(A, I0)') 'interp_bias = ', interp_bias
i = i + 1
call get_command_argument(i, command_arg)
read(command_arg, '(G23.0)') abs_tolerance
write(*, '(A, G0.3)') 'abs_tolerance = ', abs_tolerance
i = i + 1
call get_command_argument(i, command_arg)
read(command_arg, '(G23.0)') rel_tolerance
write(*, '(A, G0.3)') 'rel_tolerance = ', rel_tolerance
i = i + 1
call get_command_argument(i, command_arg)
read(command_arg, '(G23.0)') ds_out
write(*, '(A, G0.3)') 'ds_out = ', ds_out
i = i + 1
call get_command_argument(i, command_arg)
read(command_arg, '(G23.0)') ds_initial
write(*, '(A, G0.3)') 'ds_initial = ', ds_initial
i = i + 1
call get_command_argument(i, command_arg)
read(command_arg, '(G23.0)') error_initial
write(*, '(A, G0.3)') 'error_initial = ', error_initial
i = i + 1
call get_command_argument(i, command_arg)
read(command_arg, '(G23.0)') beta
write(*, '(A, G0.3)') 'beta = ', beta
i = i + 1
call get_command_argument(i, command_arg)
read(command_arg, '(G23.0)') safety_factor
write(*, '(A, G0.3)') 'safety_factor = ', safety_factor
i = i + 1
call get_command_argument(i, command_arg)
read(command_arg, '(G23.0)') scale_min
write(*, '(A, G0.3)') 'scale_min = ', scale_min
i = i + 1
call get_command_argument(i, command_arg)
read(command_arg, '(G23.0)') scale_max
write(*, '(A, G0.3)') 'scale_max = ', scale_max
i = i + 1
call get_command_argument(i, command_arg)
read(command_arg, '(G23.0)') decoupling_z
write(*, '(A, G0.3)') 'decoupling_z = ', decoupling_z
i = i + 1
call get_command_argument(i, command_arg)
read(command_arg, '(G23.0)') decoupling_beta
write(*, '(A, G0.3)') 'decoupling_beta = ', decoupling_beta
i = i + 1
call get_command_argument(i, command_arg)
read(command_arg, '(G23.0)') decoupling_rate
write(*, '(A, G0.3)') 'decoupling_rate = ', decoupling_rate
i = i + 1
call get_command_argument(i, command_arg)
read(command_arg, '(G23.0)') slope_correction
write(*, '(A, G0.3)') 'slope_correction = ', slope_correction
i = i + 1
call get_command_argument(i, command_arg)
read(command_arg, '(I23)')   max_output_points
write(*, '(A, I0)') 'max_output_points = ', max_output_points
i = i + 1
call get_command_argument(i, command_arg)
read(command_arg, '(I23)')   n_aux
write(*, '(A, I0)') 'n_aux = ', n_aux

n_params = i

allocate(aux_names(n_aux))

do i = 1, n_aux
    call get_command_argument(n_params + i, aux_names(i))
    aux_names(i) = adjustl(aux_names(i))
end do

n_params = n_params + n_aux
n_initials = (n_args - n_params)/3

allocate(x_initials(n_initials), &
         y_initials(n_initials), &
         z_initials(n_initials))

do i = 1, n_initials

    call get_command_argument(n_params + 3*i - 2, command_arg)
    read(command_arg, '(G23.0)') x_initials(i)

    call get_command_argument(n_params + 3*i - 1, command_arg)
    read(command_arg, '(G23.0)') y_initials(i)

    call get_command_argument(n_params + 3*i, command_arg)
    read(command_arg, '(G23.0)') z_initials(i)

end do

call initialize_mesh('/home/lars/Data/en024031_emer3.0med/en024031_emer3.0med_438.idl')

open(newunit=unit, file='fieldlines.dat', status='replace', action='write')

write(unit, *) xs, xe, ys, ye, zs, ze, n_output_points
write(unit, *) u_l*xm(xs:xe)
write(unit, *) u_l*ym(ys:ye)
write(unit, *) u_l*zm(zs:ze)

elapsed_time = 0.0

if (directions == 0.0) then

    do i = 1, n_initials

        write(*, '(G0.3, A, G0.3, A, G0.3)') x_initials(i), ', ', y_initials(i), ', ', z_initials(i)

        direction = 1.0

        call system_clock(start_count, count_rate)

        call trace(x_initials(i), y_initials(i), z_initials(i), &
                   n_output_points,                             &
                   x_out, y_out, z_out,                         &
                   aux_out,                                     &
                   fieldline_length,                            &
                   decoupling_index)

        call system_clock(end_count)
        elapsed_time = elapsed_time + (end_count - start_count)/real(count_rate, DP)

        write(unit, *) fieldline_length
        write(unit, *) decoupling_index
        write(unit, *) x_out
        write(unit, *) y_out
        write(unit, *) z_out
        do j = 1, n_aux
            write(unit, *) aux_out(:, j)
        end do

        direction = -1.0

        call system_clock(start_count, count_rate)

        call trace(x_initials(i), y_initials(i), z_initials(i), &
                   n_output_points,                             &
                   x_out, y_out, z_out,                         &
                   aux_out,                                     &
                   fieldline_length,                            &
                   decoupling_index)

        call system_clock(end_count)
        elapsed_time = elapsed_time + (end_count - start_count)/real(count_rate, DP)

        write(unit, *) fieldline_length
        write(unit, *) decoupling_index
        write(unit, *) x_out
        write(unit, *) y_out
        write(unit, *) z_out
        do j = 1, n_aux
            write(unit, *) aux_out(:, j)
        end do

    end do

else

    direction = directions

    do i = 1, n_initials

        write(*, '(G0.3, A, G0.3, A, G0.3)') x_initials(i), ', ', y_initials(i), ', ', z_initials(i)

        call system_clock(start_count, count_rate)

        call trace(x_initials(i), y_initials(i), z_initials(i), &
                   n_output_points,                             &
                   x_out, y_out, z_out,                         &
                   aux_out,                                     &
                   fieldline_length,                            &
                   decoupling_index)

        call system_clock(end_count)
        elapsed_time = elapsed_time + (end_count - start_count)/real(count_rate, DP)

        write(unit, *) fieldline_length
        write(unit, *) decoupling_index
        write(unit, *) x_out
        write(unit, *) y_out
        write(unit, *) z_out
        do j = 1, n_aux
            write(unit, *) aux_out(:, j)
        end do

    end do

end if

write(unit, *) elapsed_time/real(n_initials)

close(unit)

call free_mesh()

deallocate(x_out)
deallocate(y_out)
deallocate(z_out)
deallocate(aux_out)

deallocate(x_initials)
deallocate(y_initials)
deallocate(z_initials)

deallocate(aux_names)

end program test_tracer
