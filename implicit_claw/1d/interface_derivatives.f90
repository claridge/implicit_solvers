! These routines evaluate 1-dimensional derivatives at an interface by using
! 4-cell centered stencils.  Mostly intended to use arbitrary boundary
! conditions when a true solution is available.
!
! They will work in a pinch anywhere you need a derivative, but performance
! takes a big hit if you make pointwise function calls at every grid cell.

double precision function derivative0(q)
    implicit none
    double precision, dimension(4), intent(in) :: q
    derivative0 = -0.0625d0 * q(1) + 0.5625d0 * q(2)  &
        + 0.5625d0 * q(3) - 0.0625d0 * q(4)
end function derivative0


double precision function derivative1(q, dx)
    implicit none
    double precision, dimension(4), intent(in) :: q
    double precision, intent(in) :: dx
    derivative1 = (4.166666666666667d-2 * q(1) - 1.125d0 * q(2)  &
        + 1.125d0 * q(3) - 4.166666666666667d-2 * q(4)) / dx
end function derivative1


double precision function derivative2(q, dx)
    implicit none
    double precision, dimension(4), intent(in) :: q
    double precision, intent(in) :: dx
    derivative2 = (.5d0 * q(1) - .5d0 * q(2) - .5d0 * q(3) + .5d0 * q(4)) / dx**2
end function derivative2


double precision function derivative3(q, dx)
    implicit none
    double precision, dimension(4), intent(in) :: q
    double precision, intent(in) :: dx
    derivative3 = (-q(1) + 3.d0 * q(2) - 3.d0 * q(3) + q(4)) / dx**3
end function derivative3


double precision function gen_derivative(i, q, dx)
    implicit none
    integer, intent(in) :: i
    double precision, dimension(4), intent(in) :: q
    double precision, intent(in) :: dx
    double precision, external :: derivative0, derivative1, derivative2, derivative3
    
    if (i == 0) then
        gen_derivative = derivative0(q)
    else if (i == 1) then
        gen_derivative = derivative1(q, dx)
    else if (i == 2) then
        gen_derivative = derivative2(q, dx)
    else if (i == 3) then
        gen_derivative = derivative3(q, dx)
    end if

end function gen_derivative
