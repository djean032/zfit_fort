program main
  use odepack_mod
  use minpack_module
  use stdlib_constants
  use laser
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  type(lsoda_class) :: ls
  integer :: neq_rates, neq_intensity, ierr_rates, ierr_intensity, &
             istate, itask, i, j
  real(dp) :: z, zout, t, tout, rtol, atol(1), y(5), par(10), &
              sum, wavelength, frq, phi, start, finish, t_fin, Wz

  ! Set the parameters
  neq_rates = 5
  wavelength = 532e-9_dp
  frq = c/wavelength
  par(:) = [1000._dp, frq, 3.2e-18_dp, 16.e-18_dp, 14.e-18_dp, &
            1.0e-15_dp, 1.25e-9_dp, 1.0e-15_dp, 250.e-6_dp, 1.0e-15_dp]

  !call ls%initialize(rhs_rates, neq_rates, istate=istate)

  ! Set the initial conditions
  y(:) = [1.75e18_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]

  ! Set solver parameters
  t_fin = 10.0e-9_dp
  t = 0.0_dp
  tout = 10.0e-9_dp
  rtol = 1.0e-8_dp
  atol(1) = 1.0e-8_dp
  itask = 1
  istate = 1
  call cpu_time(start)
  Wz = W(1.0e-3_dp, 0.0_dp, 1.0_dp, wavelength)
  print *, Wz
!  do
!    call ls%integrate(y, t, tout, rtol, atol, itask, istate)
!    if (t .ge. t_fin) exit
!    tout = tout + 10.0e-10_dp
!    j = j + 1
!  end do

contains
  subroutine rhs_rates(self, neq, t, y, ydot, ierr)
    class(lsoda_class), intent(inout) :: self
    integer, intent(in) :: neq
    real(8), intent(in) :: t
    real(8), dimension(neq), intent(in) :: y
    real(8), dimension(neq), intent(out) :: ydot
    integer, intent(out) :: ierr
    ydot(1) = -(par(3)*y(1)*par(1))/(h*par(2)) &
              + (y(2)/par(6)) &
              + (y(4)/par(9))
    ydot(2) = (par(3)*y(1)*par(1))/(h*par(2)) &
              - (y(2)/par(6)) &
              - (par(4)*y(2)*par(1))/(h*par(2)) &
              + (y(3)/par(8)) &
              - (y(4)/par(9))
    ydot(3) = (par(4)*y(2)*par(1))/(h*par(2)) &
              - (y(3)/par(8))
    ydot(4) = -(par(5)*y(4)*par(1))/(h*par(2)) &
              + (y(5)/par(10)) &
              + (y(2)/par(7)) &
              - (y(4)/par(9))
    ydot(5) = (par(5)*y(4)*par(1))/(h*par(2)) &
              - (y(5)/par(10))
    ierr = 0
  end subroutine rhs_rates
end program main
