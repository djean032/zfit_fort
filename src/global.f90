module global
   use stdlib_constants
   use odepack_mod
   use, intrinsic :: iso_fortran_env, only: dp => real64
   implicit none

   integer, parameter :: neq_rates = 5
   integer, parameter :: neq_intensity = 1
   real(dp), parameter :: tau = 10.0e-9_dp
   real(dp), parameter :: wid = tau/(2*(log(2.0_dp))**0.5)
   real(dp), parameter :: w0 = 10.0e-5_dp
   real(dp), parameter :: M2 = 1.2_dp
   real(dp), parameter :: ep = 0.9e-3_dp
   real(dp), parameter :: phi_pk = 0.94_dp * ep / (tau)
   real(dp), parameter :: wavelength = 532e-9_dp
   real(dp), parameter :: frq = c/wavelength
   real(dp), dimension(10) :: par
   real(dp), dimension(5, 100) :: pop
   real(dp), dimension(5) :: current_pop

   public

contains
   subroutine linspace(from, to, array)
    real(dp), intent(in) :: from, to
    real(dp), intent(out) :: array(:)
    real(dp) :: range
    integer :: n, i
    n = size(array)
    range = to - from

    if (n == 0) return

    if (n == 1) then
        array(1) = from
        return
    end if


    do i=1, n
        array(i) = from + range * (i - 1) / (n - 1)
    end do
   end subroutine

end module global

