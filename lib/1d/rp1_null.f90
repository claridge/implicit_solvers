subroutine rp1(maxmx, meqn, mwaves, mbc, mx,  &
               q_left, q_right, aux_left, aux_right,  &
               waves, wave_speeds, a_minus_dq, a_plus_dq)

!----------------------
! Null Riemann solver.
!----------------------

    implicit none

    !==== Variables ====
    
    !---- Inputs ----
    integer, intent(in) :: maxmx, meqn, mwaves, mbc, mx
    double precision, intent(in), dimension(1-mbc:mx+mbc,meqn) ::  &
        q_left, q_right, aux_left, aux_right
    
    !---- Input/output ----
    double precision, intent(out), dimension(1-mbc:maxmx+mbc, mwaves) :: wave_speeds
    double precision, intent(out), dimension(1-mbc:maxmx+mbc, meqn, mwaves) :: waves
    double precision, intent(out), dimension(1-mbc:maxmx+mbc, meqn) :: a_minus_dq, a_plus_dq
    
    !---- Local ----
    integer :: i
    double precision :: ql, qr
    
    do i=2-mbc,mx+mbc
        qr = 0.d0
        ql = 0.d0
        waves(i,1,1) = 0.d0
        wave_speeds(i,1) = 0.d0
        a_minus_dq(i,1) = 0.d0
        a_plus_dq(i,1) = 0.d0
    end do

end
