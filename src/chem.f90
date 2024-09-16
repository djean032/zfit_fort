module chem
    use global
    implicit none
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

    private

    public :: rhs_rates, rhs_intensity, jdum, solve_rates, solve_intensity

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

end module chem

