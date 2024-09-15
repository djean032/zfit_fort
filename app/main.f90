program main
    use minpack_module
    use global
    use laser
    use chem
    use stdlib_math, only: linspace
    implicit none
    integer :: i, j, pos_idx
    real(dp) :: z, z0, zout, t0, tout, start, finish
    real(dp), allocatable :: t(:)
    real(dp), allocatable :: intensities(:, :)
    real(dp), allocatable :: final_intensities(:, :)
    real(dp), allocatable :: tot(:)
    real(dp), allocatable :: int0(:)
    real(dp), allocatable :: z_pos(:)
    real(dp), allocatable :: z_sample(:)
    allocate(t(t_slices))
    allocate(z_sample(z_slices))
    allocate(intensities(z_pos_slices, t_slices))
    allocate(final_intensities(z_pos_slices, t_slices))
    allocate(tot(z_pos_slices))
    allocate(int0(z_pos_slices))
    allocate(z_pos(z_pos_slices))
    allocate(pop(5, z_pos_slices))



    ! Set the initial conditions
    call cpu_time(start)
    ! Initialize y vectors for DLSODA
    initial_pop(:) = [1.75e18_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    current_pop(:) = [1.75e18_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                      0._dp, frq, 3.2e-18_dp, 16.e-18_dp, 14.e-18_dp, &
                      1.0e-15_dp, 1.25e-9_dp, 1.0e-15_dp, 250.e-6_dp, &
                      1.0e-15_dp]
    current_intensity(:) = [0._dp, 1.75e18_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                            0.0_dp, frq, 3.2e-18_dp, 16.e-18_dp, 14.e-18_dp, &
                            1.0e-15_dp, 1.25e-9_dp, 1.0e-15_dp, 250.e-6_dp, &
                            1.0e-15_dp]
    ! Initialize t and z values for DLSODA
    t0 = 0.0_dp
    tout = 1.0e-8_dp
    z_pos = linspace(-30.0e-2_dp, 30.0e-2_dp, z_pos_slices)
    z_sample = linspace(0.0_dp, 1.0e-3_dp, z_slices)
    t = linspace(-5.0e-8_dp, 5.0e-8_dp, t_slices)
    z = 0.0_dp
    zout = z_sample(2) - z_sample(1)
    do concurrent (j=1:z_pos_slices)
        z = z_pos(j)
        do concurrent(i=1:t_slices)
            intensities(j, i) = irradiance(0.0_dp, t(i), z)
        end do
    end do
    do concurrent (pos_idx=1:z_pos_slices)
        do concurrent(j=1:z_slices)
            pop(:, j) = initial_pop
        end do
        do i=1, t_slices
            current_intensity(1) = intensities(pos_idx, i)
            do j=1, z_slices
                 current_pop(1:5) = pop(:, j)
                 current_pop(6) = current_intensity(1)
                 current_pop = solve_rates(current_pop, t0, tout)
                 current_intensity(2:6) = current_pop(1:5)
                 pop(:, j) = current_pop(1:5)
                 current_intensity = solve_intensity(current_intensity, z0, zout)
            end do
            final_intensities(pos_idx, i) = current_intensity(1)
        end do
    end do
    call cpu_time(finish)
    tot = sum(final_intensities, dim=2)
    int0 = sum(intensities, dim=2)
    open(1, file="intensities.dat", status="unknown")
    do i=1, z_pos_slices
        write(1, *) z_pos(i), tot(i)/int0(i)
    end do
    print *, "Time taken: ", finish - start

contains
end program main


