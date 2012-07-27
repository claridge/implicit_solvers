subroutine qinit(maxmx, meqn, mbc, mx, x_low, dx, q, maux, aux)

    implicit none

    double precision :: h_left, h_right
    common /boundary_config/ h_left, h_right

    integer, intent(in) :: maxmx, meqn, mbc, mx, maux
    double precision, intent(in) :: x_low, dx, aux
    double precision, intent(out) :: q(1-mbc:mx+mbc,meqn)
    
    integer :: i
    double precision :: x


    do i = 1, mx
        x = x_low + (i-0.5d0)*dx
        q(i,1) = h_left + (h_right - h_left) * (1.d0 + tanh(x/(.2)))/2
    end do

end subroutine qinit
