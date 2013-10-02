subroutine qinit(maxmx, maxmy, meqn, mbc, mx, my, x_lower, y_lower, dx, dy, q)

! Sets initial condition.

    !$ use omp_lib
    implicit none

    integer, intent(in) :: maxmx, maxmy, meqn, mbc, mx, my
    double precision, intent(in) :: x_lower, y_lower, dx, dy
    double precision, intent(out) :: q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)

    integer :: ix, iy
    double precision :: x, y
    double precision, external :: true_solution

    !$omp parallel do private(x, y, ix)
    do iy = 1, my
        y = y_lower + (iy-0.5d0)*dy
        do ix = 1, mx
            x = x_lower + (ix-0.5d0)*dx
            q(ix, iy, 1) = true_solution(x, y, 0.d0)
        end do
    end do

end subroutine qinit
