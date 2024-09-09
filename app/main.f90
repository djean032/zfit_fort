program main
  use minpack_module
  use global
  use laser
  use chem
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

  call ls%initialize(rhs_rates, neq_rates, istate=istate)

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
    call cpu_time(finish)
  print *, Wz
    print *, 'CPU time = ', finish - start
!  do
!    call ls%integrate(y, t, tout, rtol, atol, itask, istate)
!    if (t .ge. t_fin) exit
!    tout = tout + 10.0e-10_dp
!    j = j + 1
!  end do

contains
end program main
