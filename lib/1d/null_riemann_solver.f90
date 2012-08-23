subroutine bc1(maxmx, meqn, mbc, mx, x_low, dx, q, maux, aux, t, dt, mthbc)
! Fills ghost cell values; called prior to rp1.
!
! This is a dummy version that does nothing.
    implicit none
    integer, intent(in) :: maxmx, meqn, mbc, mx, mthbc(2), maux
    double precision, intent(in) :: x_low, dx, aux, t, dt
    double precision, intent(inout) :: q(1-mbc:mx+mbc,meqn)
end subroutine bc1


subroutine b4step1(maxmx, mbc, mx, meqn, q, x_lower, dx, t, dt, maux, aux)
! Called once before each call to step1 (1-dimensional time-stepping
! routine).  Use to set time-dependent aux arrays, or perform other
! tasks that must be done every time step.
!
! This is a dummy verison that does nothing.
    implicit none
    integer, intent(in) :: maxmx, mbc, mx, meqn, maux
    double precision, intent(in) :: x_lower, dx, t, dt
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(in) :: q
    double precision, dimension(1-mbc:mx+mbc, maux), intent(in) :: aux 
end subroutine b4step1




subroutine rp1(maxmx, meqn, mwaves, mbc, mx,  &
               q_left, q_right, aux_left, aux_right,  &
               waves, wave_speeds, a_minus_dq, a_plus_dq)
! Riemann solver.
!
! This is a dummy version that will not affect the solution.
!
! TODO: Possible to make sense for mwaves == 0?

    implicit none

    integer, intent(in) :: maxmx, meqn, mwaves, mbc, mx
    double precision, intent(in), dimension(1-mbc:mx+mbc,meqn) ::  &
        q_left, q_right, aux_left, aux_right    
    double precision, intent(out), dimension(1-mbc:maxmx+mbc, mwaves) :: wave_speeds
    double precision, intent(out), dimension(1-mbc:maxmx+mbc, meqn, mwaves) :: waves
    double precision, intent(out), dimension(1-mbc:maxmx+mbc, meqn) :: a_minus_dq, a_plus_dq

    waves = 0.d0
    wave_speeds = 0.d0
    a_minus_dq = 0.d0
    a_plus_dq = 0.d0

end subroutine rp1
