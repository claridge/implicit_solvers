subroutine qinit(maxmx, meqn, mbc, mx, x_lower, dx, q, maux, aux)

! Sets up the initial condition.

    implicit none

    integer, intent(in) :: maxmx, meqn, mbc, mx, maux
    double precision, intent(in) :: x_lower, dx, aux
    double precision, intent(out) :: q(1-mbc:mx+mbc,meqn)
    
    integer :: ix
    double precision :: x
    double precision, external :: true_solution


    do ix = 1, mx
        x = x_lower + (ix - 0.5d0) * dx
        q(ix, 1) = 1.d0
        if (abs(x) < 1.d0) then
            q(ix, 2) = 2 * (1.d0 - abs(x)**10)
        else
            q(ix, 2) = 0.d0
        end if
    end do

end subroutine qinit
