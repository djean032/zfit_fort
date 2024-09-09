module chem
  use global
  use odepack_mod
  implicit none

  private

  public :: rhs_rates 

contains
  subroutine rhs_rates(self, neq, t, y, ydot, ierr)
    class(lsoda_class), intent(inout) :: self
    integer, intent(in) :: neq
    real(dp), intent(in) :: t
    real(dp), dimension(10) :: par
    real(dp), dimension(neq), intent(in) :: y
    real(dp), dimension(neq), intent(out) :: ydot
    integer, intent(out) :: ierr
    ydot(1) = -(par(3)*y(1)*par(1))/(h*par(2)) &
              + (y(2)/par(6)) &
              + (y(4)/par(9))
    ydot(2) = (par(3)*y(1)*par(1))/(h*par(2)) &
              - (y(2)/par(6)) &
              - (par(4)*y(2)*par(1))/(h*par(2)) &
              + (y(3)/par(dp)) &
              - (y(4)/par(9))
    ydot(3) = (par(4)*y(2)*par(1))/(h*par(2)) &
              - (y(3)/par(dp))
    ydot(4) = -(par(5)*y(4)*par(1))/(h*par(2)) &
              + (y(5)/par(10)) &
              + (y(2)/par(7)) &
              - (y(4)/par(9))
    ydot(5) = (par(5)*y(4)*par(1))/(h*par(2)) &
              - (y(5)/par(10))
    ierr = 0
  end subroutine rhs_rates
end module chem

