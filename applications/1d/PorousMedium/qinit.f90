subroutine qinit(maxmx, meqn, mbc, mx, x_lower, dx, q, maux, aux)

    implicit none

    integer, intent(in) :: maxmx, meqn, mbc, mx, maux
    double precision, intent(in) :: x_lower, dx, aux
    double precision, intent(out) :: q(1-mbc:mx+mbc,meqn)
    
    integer :: i
    double precision :: x, h_left, h_right
    double precision, external :: true_solution

    integer :: mx_common, mbc_common, meqn_common
    double precision :: x_lower_common, dx_common
    common /claw_config/ mx_common, mbc_common, x_lower_common, dx_common,  &
        meqn_common

    ! Copy configuration constants into common block.  This only actually needs
    ! to happen once, but I don't think there's a great place for it.  And it's
    ! cheap.
    mx_common = mx
    mbc_common = mbc
    x_lower_common = x_lower
    dx_common = dx
    meqn_common = meqn

    do i = 1,mx
        x = x_lower + (i-0.5d0)*dx
        q(i,1) = true_solution(x, 0.d0)
    end do

end subroutine qinit
