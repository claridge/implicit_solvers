subroutine qinit(maxmx, meqn, mbc, mx, x_lower, dx, q, maux, aux)

    implicit none

    integer, intent(in) :: maxmx, meqn, mbc, mx, maux
    double precision, intent(in) :: x_lower, dx, aux
    double precision, intent(out) :: q(1-mbc:mx+mbc,meqn)
    
    integer :: i
    double precision :: x
    double precision, external :: true_solution

    ! BEGIN: DON'T REMOVE
    ! This common block needs to get built before bc1 executes, as it
    ! potentially calls apply_bcs.  There's no appropriate place for user-
    ! defined code to do so, but this works.
    integer :: mx_common, mbc_common, meqn_common
    double precision :: x_lower_common, dx_common
    common /claw_config/ mx_common, mbc_common, x_lower_common, dx_common,  &
        meqn_common
    mx_common = mx
    mbc_common = mbc
    x_lower_common = x_lower
    dx_common = dx
    meqn_common = meqn
    ! END: DON'T REMOVE


    do i = 1, mx
        x = x_lower + (i-0.5d0)*dx
        q(i, 1) = true_solution(x, 0.d0)
        ! q(i, 1) = 2.d0 + sin(x)
    end do

end subroutine qinit
