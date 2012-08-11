subroutine rp1(maxmx, meqn, mwaves, mbc, mx,  &
               q_left, q_right, aux_left, aux_right,  &
               waves, wave_speeds, a_minus_dq, a_plus_dq)

! Null Riemann solver.
!
! TODO: Possible to make sense for mwaves == 0?
! TODO: Change file to null_riemann_solver; include b4step1.

    implicit none

    integer, intent(in) :: maxmx, meqn, mwaves, mbc, mx
    double precision, intent(in), dimension(1-mbc:mx+mbc,meqn) ::  &
        q_left, q_right, aux_left, aux_right
    
    double precision, intent(out), dimension(1-mbc:maxmx+mbc, mwaves) :: wave_speeds
    double precision, intent(out), dimension(1-mbc:maxmx+mbc, meqn, mwaves) :: waves
    double precision, intent(out), dimension(1-mbc:maxmx+mbc, meqn) :: a_minus_dq, a_plus_dq
    
    integer :: i
    
    waves = 0.d0
    wave_speeds = 0.d0
    a_minus_dq = 0.d0
    a_plus_dq = 0.d0

end subroutine rp1
