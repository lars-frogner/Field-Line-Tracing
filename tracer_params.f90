module tracer_params_mod
implicit none

private

integer, parameter :: SP = kind(0.0)

! ********************** User-modifiable parameters **********************

real(SP), public :: direction         = 1.0

logical,  public :: nonuniform_output = .false.
logical,  public :: constant_stepsize = .false.

integer,  public :: interp_order      = 3
integer,  public :: interp_bias       = 0

real(SP), public :: abs_tolerance     = 1.0e-6
real(SP), public :: rel_tolerance     = 1.0e-6

real(SP), public :: ds_out            = 1.0e-3

real(SP), public :: s_initial         = 0.0
real(SP), public :: ds_initial        = 1.0e-6
real(SP), public :: error_initial     = 1.0e-4

! For PI control, set beta = ~0.4/k (k = scheme order, e.g. k=5 for RKF45)
real(SP), public :: beta              = 0.0
real(SP), public :: safety_factor     = 0.9
real(SP), public :: scale_min         = 0.2
real(SP), public :: scale_max         = 10.0

real(SP), public :: decoupling_z      = -1.3
real(SP), public :: decoupling_beta   = 5.0
real(SP), public :: decoupling_rate   = 1.0
real(SP), public :: slope_correction  = 1.0

integer,  public :: max_output_points = 5000

integer,  public :: n_aux             = 1

character(len=5), allocatable, public :: aux_names(:)

! ************************** Derived parameters **************************

integer,  public :: n_interp
integer,  public :: start_offset

real(SP), public :: x_min, x_max, x_range
real(SP), public :: y_min, y_max, y_range
real(SP), public :: z_min, z_max, z_range

real(SP), public :: alpha

end module tracer_params_mod
