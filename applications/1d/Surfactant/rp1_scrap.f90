subroutine rp1(maxmx, meqn, mwaves, mbc, mx,  &
               q, q_unused, aux, aux_unused,  &
               waves, wave_speeds, a_minus_dq, a_plus_dq)

! Riemann solver for hyperbolic part of the thin film/surfactant problem.
!
! Args:
!   maxmx: Max number of interior cells.  Remnant of F77, where array dimensions
!     had to be specified at compile time.  If using provided driver, it's safe
!     to assume that maxmx == mx.
!   meqn: Number of equations in the PDE.
!   mwaves: Number of waves.
!   mbc: Number of ghost cells at each end of the domain.
!   mx: Actual number of grid cells.
!   q: PDE solution.
!   q_unused: Same as q for our purposes.  (In some applications, q and q_unused
!     would specify values at the left and right ends of a cell, respectively.)
!   aux: Auxiliary array.  If used, set it in b4step2.
!   aux_unused: Same as aux for our purposes.  See comment on q_unused.
!   waves: Wave values.
!   wave_speeds: Wave speeds.
!   a_minus_dq: Left-moving fluctuation.
!   a_plus_dq: Right-moving fluctuation.
!
! For interface ix, the lower cell is ix-1, and the upper cell is ix.

    implicit none

    integer, intent(in) :: mx, meqn, mwaves, mbc, mx
    double precision, intent(in), dimension(1-mbc:mx+mbc,meqn) :: q, q_unused
        q, q_unused
    double precision, intent(in), dimension(1-mbc:mx+mbc,*) :: aux, aux_unused
    double precision, intent(out), dimension(1-mbc:mx+mbc, mwaves) :: wave_speeds
    double precision, intent(out), dimension(1-mbc:mx+mbc, meqn, mwaves) :: waves
    double precision, intent(out), dimension(1-mbc:mx+mbc, meqn) :: a_minus_dq, a_plus_dq
    
    integer :: ix
    double precision :: ql, qr

    ! aux(:, 1) := cell-centered surface tension.
    call compute_surface_tension(mx + 2*mbc, q(:, 2), aux(:, 1))
    
    ! Wave 1: Film height
    do ix = 2-mbc, mx+mbc
        
        velocity_lower = q(ix-1, 1)**2/2 *  &
            (aux(ix + 1, 1) - aux(ix - 2, 1)) / (2 * dx)
        velocity_upper = q(ix, 1)**2/2 *  &
            (aux(ix + 2, 1) - aux(ix - 1, 1)) / (2 * dx)
        
        q_lower = q(i-1, 1)
        q_upper = q(i, 1)
        
        wave_values(i, 2, 1) = 0.d0

        !---- Speed and fluctuations ----
        if ((velocity_lower>=0.d0) .and. (velocity_upper>=0.d0)) then
            !-------------------------------------------------------------------
            ! All motion is to the right, so only the right-moving fluctuation
            ! a_plus_dq(i) is nonzero.  
            !
            ! The wave_speeds is determined by the Rankine-Hugoniot condition applied
            ! to the flux function
            !       f(q,x) = 1/2 * u(x) * q**2
            !-------------------------------------------------------------------
            wave_values(ix,1,1) = velocity_upper * q(ix,1)**2 / 2 -  &
                velocity_lower * q(ix-1, 1)**2 / 2
            wave_speeds(i,1) = .5d0 * (velocity_upper * q(ix, 1) +  &
                               sqrt(velocity_upper * velocity_lower) * q(ix-1, 1))
            a_minus_dq(i,1) = 0.d0
            a_plus_dq(i,1) = wave_values(ix, 1, 1)

        else if ((velocity_lower<=0.d0) .and. (velocity_upper<=0.d0)) then
            waves(i,1,1) = velocity_upper*q_left(i,1)**2/2.d0 -  &
                           velocity_lower*q_right(i-1,1)**2/2.d0
            wave_speeds(i,1) = .5d0 * (-sqrt(velocity_upper*velocity_lower)*q_left(i,1)  &
                               + velocity_lower*q_right(i-1,1))
            a_minus_dq(i,1) = waves(i,1,1)
            a_plus_dq(i,1) = 0.d0

        else  ! Velocities are in different directions
            !---- Velocities are in different directions ----
            ! In this case, regardless of whether the velocities are moving
            ! towards or away from the interface, there is no flux across the
            ! interface itself.  However, to maintain conservation, a_minus_dq must
            ! take a value that compares the left state in cell i-1 to a right
            ! state of zero flux, and similarly for a_plus_dq regarding cell i.
            !-------------------------------------------------------------------
            waves(i,1,1) = 0.d0
            wave_speeds(i,1) = 0.d0
            a_minus_dq(i,1) = -velocity_lower * q_right(i-1,1)**2/2.d0
            a_plus_dq(i,1) = velocity_upper * q_left(i,1)**2/2.d0
        end if
    end do

end subroutine rp1
