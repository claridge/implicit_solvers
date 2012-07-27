subroutine bc1(maxmx, meqn, mbc, mx, x_low, dx, q, maux, aux, t, dt, mthbc)

    implicit none
    
    integer, intent(in) :: maxmx, meqn, mbc, mx, mthbc(2), maux
    double precision, intent(in) :: x_low, dx, aux, t, dt    
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(inout) :: q

    call apply_bcs(t, q)

end subroutine bc1
