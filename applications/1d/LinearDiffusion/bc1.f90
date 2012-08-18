subroutine bc1(maxmx, meqn, mbc, mx, x_low, dx, q, maux, aux, t, dt, mthbc)

    implicit none
    
    !---> Variables --->    
    ! Inputs
    integer, intent(in) :: maxmx, meqn, mbc, mx, mthbc(2), maux
    double precision, intent(in) :: x_low, dx, aux, t, dt
    
    ! Input/output
    double precision, intent(inout) :: q(1-mbc:mx+mbc,meqn)
    !<--- Variables <---
    
    call apply_bcs(t, q(:,1))

end subroutine bc1
