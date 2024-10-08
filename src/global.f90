module global
    use stdlib_constants
    use minpack_module
    use stdlib_io, only: loadtxt
    use stdlib_math, only: linspace
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    integer, parameter :: t_slices = 10
    integer, parameter :: z_slices = 5
    integer, parameter :: z_pos_slices = 101
    real(dp), parameter :: tau = 8.0e-9_dp
    real(dp), parameter :: wid = tau/(2*(log(2.0_dp))**0.5)
    real(dp), parameter :: w0 = 14e-6_dp
    real(dp), parameter :: M2 = 1.0_dp
    real(dp), parameter :: ep = 154.e-9_dp
    real(dp), parameter :: phi_pk = 0.94_dp*ep/(tau)
    real(dp), parameter :: wavelength = 532.0e-9_dp
    real(dp), parameter :: zr = (PI_dp*w0**2)/(M2*wavelength)
    real(dp), parameter :: frq = c/wavelength
    real(dp), allocatable :: t(:)
    real(dp), allocatable :: intensities(:, :)
    real(dp), allocatable :: final_intensities(:, :)
    real(dp), allocatable :: tot(:)
    real(dp), allocatable :: int0(:)
    real(dp), allocatable :: normals(:)
    real(dp), allocatable :: z_pos(:)
    real(dp), allocatable :: z_sample(:)
    real(dp), allocatable :: exp(:, :)
    real(dp), allocatable :: z_exp(:)
    real(dp), allocatable :: transmission(:)
    real(dp), allocatable :: pop(:, :)
    real(dp), allocatable :: fvec(:)
    real(dp), dimension(15) :: current_pop
    real(dp), dimension(5) :: initial_pop
    real(dp), dimension(15) :: current_intensity
    integer :: i, j, pos_idx
    real(dp) :: z, z0, zout, t0, tout, start, finish

    public

contains

end module global

