program main
   use minpack_module
   use global
   use laser
   use chem
   implicit none
   type(lsoda_class) :: ls_rates, ls_intensity
   integer :: istate, itask, i, j
   real(dp) :: z, zout, tout, rtol, atol(1), start, finish
   real(dp), dimension(t_slices) :: t, intensities
   real(dp), dimension(z_slices) :: z_pos

   call ls_rates%initialize(rhs_rates, neq_rates, istate=istate)
   call ls_intensity%initialize(rhs_intensity, neq_intensity, istate=istate)

   rtol = 1.0e-8_dp
   atol(1) = 1.0e-8_dp
   itask = 1
   istate = 1

   ! Set the parameters
   par(:) = [0._dp, frq, 3.2e-18_dp, 16.e-18_dp, 14.e-18_dp, &
             1.0e-15_dp, 1.25e-9_dp, 1.0e-15_dp, 250.e-6_dp, 1.0e-15_dp]

   ! Set the initial conditions
   call cpu_time(start)
   pop(1:5, 1) = [1.75e18_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
   t = linspace(-5.0e-8_dp, 5.0e-8_dp, t_slices)
   z_pos = linspace(-20.0_dp, 20.0_dp, z_slices)
   z = 0.0_dp
   do concurrent(i = 1:1000)
        intensities(i) = irradiance(0.0_dp, t(i), z)
   end do
   call cpu_time(finish)
   print *, z_pos
   print *, "Time taken: ", finish - start


!   do
!      call ls%integrate(y, t, tout, rtol, atol, itask, istate)
!      if (t .ge. t_fin) exit
!      tout = tout + 10.0e-10_dp
!   end do

contains
end program main
