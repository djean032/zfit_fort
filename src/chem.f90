module chem
    use global
    implicit none

    private

    abstract interface
        function expr_f(x, pars) result(y)
            import :: dp
            implicit none
            real(dp), intent(in) :: x(:)
            real(dp), intent(in) :: pars(:)
            real(dp) :: y(size(x))
        end function expr_f
    end interface

    abstract interface
        pure subroutine odepack_f(neq, t, y, ydot)
            import :: dp
            implicit none
            integer, intent(in) :: neq
            real(dp), intent(in) :: t
            real(dp), intent(in) :: y(15)
            real(dp), intent(out) :: ydot(neq)
        end subroutine

        pure subroutine odepack_jac()
        end subroutine
    end interface
    interface
        pure subroutine DLSODA(f, neq, y, t, tout, itol, rtol, atol, itask, &
                               istate, iopt, rwork, lrw, iwork, liw, jac, jt)
            import :: odepack_f, odepack_jac
            import :: dp
            implicit none
            integer, intent(in) :: neq, itol, itask, iopt, lrw, liw, jt
            integer, intent(inout) :: istate
            real(dp), intent(in) :: rtol, atol(1)
            real(dp), intent(inout) :: t
            real(dp), intent(in) :: tout
            real(dp), intent(inout) :: y(neq)
            real(dp), intent(inout) :: rwork(lrw)
            integer, intent(inout) :: iwork(liw)
            procedure(odepack_f) :: f
            procedure(odepack_jac) :: jac
        end subroutine

    end interface

    public :: rhs_rates, rhs_intensity, jdum, solve_rates, solve_intensity, &
              solve_system, fit_scan

contains
    pure subroutine rhs_rates(neq, t, y, ydot)
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp), dimension(15), intent(in) :: y
        real(dp), dimension(neq), intent(out) :: ydot
        ydot(1) = -(y(8)*y(1)*y(6))/(h*y(7)) &
                  + (y(2)/y(11)) &
                  + (y(4)/y(14))
        ydot(2) = (y(8)*y(1)*y(6))/(h*y(7)) &
                  - (y(2)/y(11)) &
                  - (y(9)*y(2)*y(6))/(h*y(7)) &
                  + (y(3)/y(13)) &
                  - (y(2)/y(12))
        ydot(3) = (y(9)*y(2)*y(6))/(h*y(7)) &
                  - (y(3)/y(13))
        ydot(4) = -(y(10)*y(4)*y(6))/(h*y(7)) &
                  + (y(5)/y(15)) &
                  + (y(2)/y(12)) &
                  - (y(4)/y(14))
        ydot(5) = (y(10)*y(4)*y(6))/(h*y(7)) &
                  - (y(5)/y(15))
    end subroutine rhs_rates

    pure subroutine rhs_intensity(neq, t, y, ydot)
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp), dimension(15), intent(in) :: y
        real(dp), dimension(neq), intent(out) :: ydot
        ydot(1) = -y(8)*y(2)*y(1) &
                  - y(9)*y(3)*y(1) &
                  - y(10)*y(5)*y(1)
    end subroutine rhs_intensity

    pure subroutine jdum()
    end subroutine jdum

    pure function solve_rates(y, t, tout) result(y_ret)
        integer :: neq, itol, itask, iopt, lrw, &
                   liw, jt, istate
        real(dp), intent(in) :: t, tout
        real(dp) :: t_ret, tout_ret, rtol, atol(1), rwork(102)
        real(dp), intent(in) :: y(15)
        real(dp) :: y_ret(15)
        integer :: iwork(25)
        t_ret = t
        tout_ret = tout
        neq = 5
        lrw = 102
        liw = 25
        rtol = 1.0e-8_dp
        atol(1) = 1.0e-8_dp
        itol = 1
        itask = 1
        istate = 1
        iopt = 0
        lrw = 102
        liw = 25
        jt = 2
        y_ret = y
        call DLSODA(rhs_rates, neq, y_ret, t_ret, tout_ret, itol, &
                    rtol, atol(1), itask, istate, iopt, rwork, lrw, iwork, liw, jdum, jt)
    end function solve_rates

    pure function solve_intensity(y, t, tout) result(y_ret)
        integer :: neq, itol, itask, iopt, lrw, &
                   liw, jt, istate
        real(dp), intent(in) :: t, tout
        real(dp) :: t_ret, tout_ret, rtol, atol(1), rwork(102)
        real(dp), intent(in) :: y(15)
        real(dp) :: y_ret(15)
        integer :: iwork(25)
        t_ret = t
        tout_ret = tout
        neq = 1
        lrw = 102
        liw = 25
        rtol = 1.0e-6_dp
        atol(1) = 1.0e-6_dp
        itol = 1
        itask = 1
        istate = 1
        iopt = 0
        lrw = 102
        liw = 25
        jt = 2
        y_ret = y
        call DLSODA(rhs_intensity, neq, y_ret, t_ret, tout_ret, itol, &
                    rtol, atol(1), itask, istate, iopt, rwork, lrw, iwork, liw, jdum, jt)
    end function solve_intensity

    function solve_system(x, pars) result(y)
        real(dp), intent(in) :: x(:)
        real(dp), intent(in) :: pars(:)
        real(dp) :: y(size(x)), normals(size(x))
        initial_pop(:) = [2.90869e17_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
        current_pop(1:5) = initial_pop
        current_intensity(2:6) = initial_pop
        current_intensity(10) = pars(1)
        current_pop(10) = pars(1)
        t0 = 0.0_dp
        tout = t(2) - t(1)
        z0 = 0.0_dp
        zout = z_sample(2)
        do concurrent(pos_idx=1:z_pos_slices)
            do concurrent(j=1:z_slices)
                pop(:, j) = initial_pop
            end do
            do i = 1, t_slices
                current_intensity(1) = intensities(pos_idx, i)
                do j = 1, z_slices
                    current_pop(1:5) = pop(:, j)
                    current_pop(6) = current_intensity(1)
                    current_pop = solve_rates(current_pop, t0, tout)
                    pop(:, j) = current_pop(1:5)
                    current_intensity(2:6) = current_pop(1:5)
                    current_intensity = solve_intensity(current_intensity, z0, zout)
                end do
                final_intensities(pos_idx, i) = current_intensity(1)
            end do
        end do
        y = sum(final_intensities, dim=2)
        int0 = sum(intensities, dim=2)
        normals = y/int0
        y = y/int0 + (1.0_dp - maxval(normals))
    end function solve_system

    subroutine fit_scan(data_x, data_y, expr, pars, fvec)
        real(dp), intent(in) :: data_x(:)
        real(dp), intent(in) :: data_y(:)
        procedure(expr_f) :: expr
        real(dp), intent(out) :: pars(:)

        real(dp) :: tol
        real(dp), intent(out) :: fvec(size(data_x))
        integer :: iwa(size(pars)), info, m, n
        real(dp), allocatable :: wa(:)

        tol = sqrt(epsilon(1.0_dp))
        m = size(fvec)
        n = size(pars)
        allocate (wa(2*m*n + 5*n + m))
        call lmdif1(fcn, m, n, pars, fvec, tol, info, iwa, wa, size(wa))

    contains

        subroutine fcn(m, n, x, fvec, iflag)
            integer, intent(in) :: m, n
            integer, intent(inout) :: iflag
            real(dp), intent(in) :: x(n)
            real(dp), intent(out) :: fvec(m)
            real(dp) :: y(size(data_x))
            fvec(1) = iflag
            y = expr(data_x, x)
            fvec = (data_y - y)
        end subroutine fcn
    end subroutine fit_scan

end module chem

