module global
    use stdlib_constants
    use stdlib_math, only: linspace
    use odepack_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none

    integer, parameter :: neq_rates = 5
    integer, parameter :: neq_intensity = 1
    integer, parameter :: t_slices = 101
    integer, parameter :: z_pos_slices = 101
    integer, parameter :: z_slices = 50
    real(dp), parameter :: tau = 10.0e-9_dp
    real(dp), parameter :: wid = tau/(2*(log(2.0_dp))**0.5)
    real(dp), parameter :: w0 = 10.0e-5_dp
    real(dp), parameter :: M2 = 1.2_dp
    real(dp), parameter :: ep = 0.9e-3_dp
    real(dp), parameter :: phi_pk = 0.94_dp*ep/(tau)
    real(dp), parameter :: wavelength = 532e-9_dp
    real(dp), parameter :: frq = c/wavelength
    real(dp), dimension(10) :: par
    real(dp), allocatable :: pop(:, :)
    real(dp), dimension(5) :: current_pop

    public

contains

end module global

