subroutine qinit(maxmx, meqn, mbc, mx, x_lower, dx, q, maux, aux)

    implicit none

    integer, intent(in) :: maxmx, meqn, mbc, mx, maux
    double precision, intent(in) :: x_lower, dx, aux
    double precision, intent(out) :: q(1-mbc:mx+mbc,meqn)
    
    integer :: i
    double precision :: x


    do i = 1, mx
        x = x_lower + (i-0.5d0)*dx
        call set_true_solution(x, 0.d0, q(i, :))
    end do

end subroutine qinit
