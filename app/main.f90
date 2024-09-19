program main
    use global
    use chem
    use laser
    implicit none

    real(dp) :: params(1)
    integer :: loc(1), idx, zdx, tdx
    real(dp) :: initial_params(5)
    real(dp) :: rnorm(5)
    real(dp) :: final_params(5)

    call cpu_time(start)
    ! Allocate arrays
    allocate (t(t_slices))
    allocate (z_sample(z_slices))
    allocate (intensities(z_pos_slices, t_slices))
    allocate (final_intensities(z_pos_slices, t_slices))
    allocate (tot(z_pos_slices))
    allocate (int0(z_pos_slices))
    allocate (normals(z_pos_slices))
    allocate (fvec(z_pos_slices))
    allocate (pop(5, z_pos_slices))

    ! Load experimental data
    call loadtxt('./data/6-Br/6-Br-pbt-532nm-154.2nJ.all', exp, skiprows=10, fmt="*")
    z_pos = exp(:, 1)*1.0e-3_dp
    transmission = exp(:, 3)/exp(:, 2)
    transmission = transmission + (1.0_dp - maxval(transmission))
    loc = minloc(transmission)
    z_pos = z_pos - z_pos(loc(1))

    ! Initialize t and z values for DLSODA
    z_sample = linspace(0.0_dp, 1.0e-1_dp, z_slices)
    t = linspace(-4.0e-9_dp, 4.0e-9_dp, t_slices)

    ! Initialize intensities
    do concurrent(zdx=1:z_pos_slices)
        z = z_pos(zdx)
        do concurrent(tdx=1:t_slices)
            intensities(zdx, tdx) = irradiance(0.0_dp, t(tdx), z)
        end do
    end do

    ! Initialize y vectors for DLSODA
    initial_pop(:) = [2.90869e17_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
    current_pop(:) = [2.90869e17_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                      0.0_dp, frq, 5.48e-18_dp, 16.0e-18_dp, 28.1e-18_dp, &
                      1.0e-12_dp, 1.0e-12_dp, 1.0e-12_dp, 1.05e-7_dp, &
                      1.0e-12_dp]
    current_intensity(:) = [0.0_dp, 2.90869e17_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                            frq, 5.48e-18_dp, 16.0e-18_dp, 28.1e-18_dp, &
                            1.0e-12_dp, 1.0e-12_dp, 1.0e-12_dp, 1.05e-7_dp, &
                            1.0e-12_dp]
    initial_params = [2.8e-14_dp, 2.8e-15_dp, 2.8e-16_dp, 2.8e-17_dp, 2.8e-18_dp]

    ! Solve system
    do idx = 1, size(initial_params)
        params(1) = initial_params(idx)
        call fit_scan(z_pos, transmission, solve_system, params, fvec)
        final_params(idx) = params(1)
        rnorm(idx) = enorm(z_pos_slices, fvec)
    end do
    print *, "Final parameters: ", final_params
    print *, "Final residuals: ", rnorm
    loc = minloc(rnorm)
    print *, "Best fit: ", final_params(loc(1))
    params(1) = final_params(loc(1))
    tot = solve_system(z_pos, params)

    call cpu_time(finish)
    int0 = sum(intensities, dim=2)
    open (1, file="intensities.dat", status="unknown")
    open (2, file="exp.dat", status="unknown")
    do i = 1, z_pos_slices
        write (1, *) z_pos(i), tot(i)
        write (2, *) z_pos(i), transmission(i)
    end do
    print *, "Time taken: ", finish - start
contains
end program main
