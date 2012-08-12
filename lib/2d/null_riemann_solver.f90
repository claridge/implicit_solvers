subroutine b4step2(maxmx, maxmy, mbc, mx, my, meqn, q,  &
                   x_low, y_low, dx, dy, t, dt, maux, aux)

! Called from claw2 before each call to step2.  Used to set the 
! aux arrays with velocity field data.
!
! "step2" means "2-dimensional step routine", not "step number 2".

    integer, intent(in) :: maxmx, maxmy, mbc, mx, my, meqn, maux
    double precision, intent(in) :: x_low, y_low, dx, dy, t, dt
    double precision, dimension(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn),  &
        intent(in) :: q

    ! TODO: This probably crashes if maux == 0.
    double precision, dimension(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux),  &
        intent(out) :: aux

    ! Null version; does nothing.

end subroutine b4step2


subroutine rpn2(slice_direction, maxm, meqn, mwaves, mbc, mx,  &
                q_left, q_right, aux_left, aux_right,  &
                waves, wave_speeds, a_minus_dq, a_plus_dq)

! TODO: Add description.

    implicit none

    integer, intent(in) :: slice_direction, maxm, meqn, mwaves, mbc, mx
    double precision, dimension(1-mbc:maxm+mbc,meqn), intent(in) :: q_left, q_right
    double precision, dimension(1-mbc:maxm+mbc, *), intent(in) :: aux_left, aux_right
    double precision, dimension(1-mbc:maxm+mbc, meqn, mwaves), intent(out) :: waves
    double precision, dimension(1-mbc:maxm+mbc, mwaves), intent(out) :: wave_speeds
    double precision, dimension(1-mbc:maxm+mbc, meqn), intent(out) ::  &
        a_plus_dq, a_minus_dq

    waves = 0.d0
    wave_speeds = 0.d0
    a_minus_dq = 0.d0
    a_plus_dq = 0.d0

end subroutine rpn2


subroutine rpt2(slice_direction, maxm, meqn, mwaves, mbc, mx, ql, qr, &
                aux_below, aux_aligned, aux_above,        &
                imp, asdq, bmasdq, bpasdq)

! TODO: Add description

    implicit none

    integer, intent(in) :: slice_direction, maxm, meqn, mwaves, mbc, mx, imp    
    double precision, intent(in) :: ql(1-mbc:maxm+mbc, meqn)
    double precision, intent(in) :: qr(1-mbc:maxm+mbc, meqn)
    double precision, intent(in), dimension(1-mbc:maxm+mbc,*) ::  &
        aux_below, aux_aligned, aux_above
    double precision, intent(in) :: asdq(1-mbc:maxm+mbc, meqn)
    double precision, intent(out) :: bmasdq(1-mbc:maxm+mbc, meqn)
    double precision, intent(out) :: bpasdq(1-mbc:maxm+mbc, meqn)

    bmasdq = 0.d0
    bpasdq = 0.d0

end subroutine rpt2
