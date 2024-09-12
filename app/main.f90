program main
    use minpack_module
    use global
    use laser
    use chem
    implicit none
    type(lsoda_class) :: ls_rates, ls_intensity
    integer :: istate, itask, i, j, pos_idx
    real(dp) :: z, zout, tout, rtol, atol(1), start, finish, y(5)
    real(dp), allocatable :: t(:)
    real(dp), allocatable :: intensities(:, :)
    real(dp), allocatable :: final_intensities(:, :)
    real(dp) :: current_intensity(1)
    real(dp), allocatable :: z_pos(:)
    real(dp), allocatable :: z_sample(:)
    allocate(t(t_slices))
    allocate(z_sample(z_slices))
    allocate(intensities(z_pos_slices, t_slices))
    allocate(final_intensities(z_pos_slices, t_slices))
    allocate(z_pos(z_pos_slices))
    allocate(pop(neq_rates, z_pos_slices))


    rtol = 1.0e-8_dp
    atol(1) = 1.0e-8_dp
    itask = 1
    istate = 1

    ! Set the parameters
    par(:) = [0._dp, frq, 3.2e-18_dp, 16.e-18_dp, 14.e-18_dp, &
              1.0e-15_dp, 1.25e-9_dp, 1.0e-15_dp, 250.e-6_dp, 1.0e-15_dp]

    ! Set the initial conditions
    call cpu_time(start)
    current_pop(:) = [1.75e18_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    t = linspace(-5.0e-8_dp, 5.0e-8_dp, t_slices)
    z_pos = linspace(-20.0_dp, 20.0_dp, z_pos_slices)
    z_sample = linspace(0.0_dp, 1.0e-3_dp, z_slices)
    do concurrent(j=1:z_slices)
        pop(:, j) = current_pop
    end do
    do concurrent (j=1:z_pos_slices)
        z = z_pos(j)
        do concurrent(i=1:t_slices)
            intensities(j, i) = irradiance(0.0_dp, t(i), z)
        end do
    end do
    do pos_idx=1, z_pos_slices
        do i=1, t_slices
            current_intensity = intensities(pos_idx, i)
            par(1) = intensities(pos_idx, i)
            !print *, par(1)
            do j=1, z_slices
                current_pop = pop(:, j)
                call ls_rates%initialize(rhs_rates, neq_rates, istate=istate)
                call ls_intensity%initialize(rhs_intensity, neq_intensity, istate=istate)
                call ls_rates%integrate(current_pop, t(1), t(2), rtol, atol, itask, istate)
                call ls_intensity%integrate(current_intensity, z_sample(1), z_sample(2), rtol, atol, itask, istate)
                itask = 1
                istate = 1
                pop(:, j) = current_pop
            end do
            final_intensities(pos_idx, i) = current_intensity(1)
        end do
    end do
    call cpu_time(finish)
    print *, "Time taken: ", finish - start
!   do
!      call ls%integrate(y, t, tout, rtol, atol, itask, istate)
!      if (t .ge. t_fin) exit
!      tout = tout + 10.0e-10_dp
!   end do

contains
end program main
