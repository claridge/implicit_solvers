subroutine rp1(maxmx, meqn, mwaves, mbc, mx,  &
               q_left, q_right, aux_left, aux_right,  &
               waves, wave_speeds, a_minus_dq, a_plus_dq)

! Riemann solver for hyperbolic part of the driven thin film equation.

    implicit none

    integer, intent(in) :: maxmx, meqn, mwaves, mbc, mx
    double precision, intent(in), dimension(1-mbc:mx+mbc,meqn) ::  &
        q_left, q_right, aux_left, aux_right
    double precision, intent(out), dimension(1-mbc:maxmx+mbc, mwaves) :: wave_speeds
    double precision, intent(out), dimension(1-mbc:maxmx+mbc, meqn, mwaves) :: waves
    double precision, intent(out), dimension(1-mbc:maxmx+mbc, meqn) :: a_minus_dq, a_plus_dq
    
    integer :: i
    double precision :: ql, qr
    
    do i=2-mbc,mx+mbc
        qr = q_right(i,1)
        ql = q_left(i-1,1)
        waves(i,1,1) = (qr**2 - qr**3) - (ql**2 - ql**3)
        wave_speeds(i,1) = (qr + ql) - (qr**2 + ql**2)
        a_minus_dq(i,1) = 0.d0
        a_plus_dq(i,1) = waves(i,1,1)
    end do

end subroutine rp1
