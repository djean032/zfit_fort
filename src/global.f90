module global
    use stdlib_constants
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    integer, parameter :: t_slices = 10
    integer, parameter :: z_pos_slices = 101
    integer, parameter :: z_slices = 20
    real(dp), parameter :: tau = 8.0e-9_dp
    real(dp), parameter :: wid = tau/(2*(log(2.0_dp))**0.5)
    real(dp), parameter :: w0 = 30.45e-6_dp
    real(dp), parameter :: M2 = 1.0_dp
    real(dp), parameter :: ep = 400.0e-9_dp
    real(dp), parameter :: phi_pk = 0.94_dp * ep / (tau)
    real(dp), parameter :: wavelength = 532.0e-9_dp
    real(dp), parameter :: zr = (PI_dp*w0**2)/(M2*wavelength)
    real(dp), parameter :: frq = c/wavelength
    real(dp), allocatable :: pop(:, :)
    real(dp), dimension(15) :: current_pop
    real(dp), dimension(5) :: initial_pop
    real(dp), dimension(15) :: current_intensity

    public

contains

end module global

