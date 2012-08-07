subroutine b4step1(maxmx, mbc, mx, meqn, q, x_lower, dx, t, dt, maux, aux)

    implicit none

    double precision :: beta, kappa, delta, mu
    common /physics_config/ beta, kappa, delta, mu

    integer, intent(in) :: maxmx, mbc, mx, meqn, maux
    double precision, intent(in) :: x_lower, dx, t, dt
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(in) :: q
    double precision, dimension(1-mbc:mx+mbc, maux), intent(out) :: aux

    double precision, dimension(1-mbc:mx+mbc) :: surface_tension
    integer :: ix

    call compute_surface_tension(mx + 2*mbc, q(:, 2), surface_tension)
    do ix = 1, mx
        aux(ix, 1) = q(ix, 1)**2 / 2 *  &
            (surface_tension(ix+1) - surface_tension(ix-1)) / (2 * dx)
    end do
    call extend_to_ghost_cells('odd', 'odd', aux(:, 1))

    ! velocity := -1/2 * beta * h^2 * h_x(i,j)
    do ix = 1, mx
        aux(ix, 2) = -.5d0 * beta * q(ix, 1)**2 *  &
            (q(ix+1, 1) - q(ix-1, 1)) / (2 * dx)
    end do
    call extend_to_ghost_cells('odd', 'odd', aux(:, 2))

end subroutine b4step1



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

    integer, intent(in) :: maxmx, meqn, mwaves, mbc, mx
    double precision, dimension(1-mbc:maxmx+mbc,meqn), intent(in) :: q, q_unused
    double precision, dimension(1-mbc:maxmx+mbc, *), intent(inout) :: aux, aux_unused
    double precision, dimension(1-mbc:maxmx+mbc, meqn, mwaves), intent(out) :: waves
    double precision, dimension(1-mbc:maxmx+mbc, mwaves), intent(out) :: wave_speeds
    double precision, dimension(1-mbc:maxmx+mbc, meqn), intent(out) ::  &
        a_plus_dq, a_minus_dq

    integer :: i
    double precision :: velocity_left, velocity_right


    !---> Film height --->
    do i=2-mbc,mx+mbc
        velocity_left = aux(i-1, 1)
        velocity_right = aux(i, 1)
        
        ! This wave doesn't affect the surfactant
        waves(i, 2, 1) = 0.d0

        !---- Speed and fluctuations ----
        if ((velocity_left >= 0.d0) .and. (velocity_right >= 0.d0)) then
            !-------------------------------------------------------------------
            ! All motion is to the right, so only the right-moving fluctuation
            ! a_plus_dq(i) is nonzero.  
            !
            ! The wave_speeds is determined by the Rankine-Hugoniot condition applied
            ! to the flux function
            !       f(q,x) = 1/2 * u(x) * q**2
            !-------------------------------------------------------------------
            waves(i,1,1) = velocity_right*q(i,1)**2/2.d0 -  &
                           velocity_left*q(i-1,1)**2/2.d0
            wave_speeds(i,1) = .5d0 * (velocity_right*q(i,1) +  &
                               sqrt(velocity_right*velocity_left)*q(i-1,1))
            a_minus_dq(i,1) = 0.d0
            a_plus_dq(i,1) = waves(i,1,1)

        else if ((velocity_left<=0.d0) .and. (velocity_right<=0.d0)) then
            waves(i,1,1) = velocity_right*q(i,1)**2/2.d0 -  &
                           velocity_left*q(i-1,1)**2/2.d0
            wave_speeds(i,1) = .5d0 * (-sqrt(velocity_right*velocity_left)*q(i,1)  &
                               + velocity_left*q(i-1,1))
            a_minus_dq(i,1) = waves(i,1,1)
            a_plus_dq(i,1) = 0.d0

        else
            !---- Velocities are in different directions ----
            ! In this case, regardless of whether the velocities are moving
            ! towards or away from the interface, there is no flux across the
            ! interface itself.  However, to maintain conservation, a_minus_dq must
            ! take a value that compares the left state in cell i-1 to a right
            ! state of zero flux, and similarly for a_plus_dq regarding cell i.
            !-------------------------------------------------------------------
            waves(i,1,1) = 0.d0
            wave_speeds(i,1) = 0.d0
            a_minus_dq(i,1) = -velocity_left * q(i-1,1)**2/2.d0
            a_plus_dq(i,1) = velocity_right * q(i,1)**2/2.d0
        end if
    end do
    !<--- Film height <---


    !---> Surfactant concentration --->
    do i = 2-mbc, mx+mbc
        velocity_left = aux(i-1, 2)
        velocity_right = aux(i, 2)

        ! This wave doesn't affect the film height
        waves(i,1,2) = 0.d0
        
        !---- Speed and fluctuations ----
        if ((velocity_left>=0.d0) .and. (velocity_right>=0.d0)) then
            !------------------------------------------------------------------
            ! All motion is to the right, so only the right-moving fluctuation
            ! a_plus_dq(i) is nonzero.  The wave_speeds is determined by the lower cell,
            ! as its mass is moving into the interface.
            !------------------------------------------------------------------
            waves(i,2,2) = velocity_right*q(i,2) - velocity_left*q(i-1,2) 
            wave_speeds(i,2) = velocity_right
            a_minus_dq(i,2)  = 0.d0
            a_plus_dq(i,2)   = waves(i,2,2)
            
        else if ((velocity_left<=0.d0) .and. (velocity_right<=0.d0)) then
            !-----------------------------------------------------------------
            ! All motion is to the left, so only the left-moving fluctuation
            ! a_minus_dq(i) is nonzero.  The wave_speeds is determined by the upper cell,
            ! as its mass is moving into the interface.
            !-----------------------------------------------------------------
            waves(i,2,2) = velocity_right*q(i,2) - velocity_left*q(i-1,2)
            wave_speeds(i,2) = velocity_left
            a_minus_dq(i,2) = waves(i,2,2)
            a_plus_dq(i,2) = 0.d0

        else  ! Velocities are in different directions
            !-------------------------------------------------------------------
            ! In this case, regardless of whether the velocities are moving
            ! towards or away from the interface, there is no flux across the
            ! interface itself.  However, to maintain conservation, a_minus_dq must
            ! take a value that compares the left state in cell i-1 to a right
            ! state of zero flux, and similarly for a_plus_dq regarding cell i.
            !-------------------------------------------------------------------
            waves(i,2,2) = 0.d0
            wave_speeds(i,2) = 0.d0    !(velocity_left + velocity_right) / 2.d0
            a_minus_dq(i,2) = -velocity_left * q(i-1,2)
            a_plus_dq(i,2) = velocity_right * q(i,2)
        end if
    end do
    !<--- Surfactant concentration <---

end subroutine rp1
