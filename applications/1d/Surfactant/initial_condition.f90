subroutine qinit(maxmx, meqn, mbc, mx, x_lower, dx, q, maux, aux)

    implicit none

    integer, intent(in) :: maxmx, meqn, mbc, mx, maux
    double precision, intent(in) :: x_lower, dx, aux
    double precision, intent(out) :: q(1-mbc:mx+mbc,meqn)
    
    integer :: ix
    double precision :: x
    double precision, external :: true_solution

    ! BEGIN: DON'T REMOVE
    ! This common block needs to get built before bc1 executes, as it
    ! potentially calls apply_bcs.  There's no appropriate place for user-
    ! defined code to do so, but this works.
    ! TODO: Do this in driver.f90 by directly reading claw.data.
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


    do ix = 1, mx
        x = x_lower + (ix - 0.5d0) * dx
        q(ix, 1) = 1.d0
        if (abs(x) < 1.d0) then
            q(ix, 2) = 2 * (1.d0 - abs(x)**10)
        else
            q(ix, 2) = 0.d0
        end if
!         q(ix, 2) = tanh(3 * (x + 1)) - tanh(3 * (x - 1)) + .1d0
    end do

end subroutine qinit
