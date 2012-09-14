subroutine rpn2(slice_direction, maxm, meqn, mwaves, mbc, mx,  &
                q_left, q_right, aux_left, aux_right,  &
                waves, wave_speeds, a_minus_dq, a_plus_dq)

!--------------------------------------------------------------------------
! This solves Riemann problems for the thin film/surfactant model.
! The film (h=q(:,:,1)) and surfactant (Gamma=q(:,:,2)) act independently.
! On input, q_left contains the state vector at the left edge of each cell,
! and q_right contains the state vector at the right edge of each cell.
!
! The film advects with a Burgers-type equation,
!     h_t + (1/2 * h**2 * sigma_x)_x + (1/2 * h**2 * sigma_y)_y = 0.
!
! In b4step2, sigma_x and sigma_y have been stored on each cell, with
!
!     aux(i,j,1) = sigma_x(i,j), and
!     aux(i,j,2) = sigma_y(i,j).
!
!
! The surfactant evolves by simple advection,
!     Gamma_t + (-1/2 * beta * h^2 * h_x * Gamma)_x 
!             + (-1/2 * beta * h^2 * h_y * Gamma)_y = 0.
!
! The velocities have been stored in the aux arrays, with
!
!     aux(i,j,3) = -1/2 * beta * h^2 * h_x(i,j)
!     aux(i,j,4) = -1/2 * beta * h^2 * h_y(i,j).
!--------------------------------------------------------------------------

    implicit none

    !---> Variables --->
    ! Input
    integer, intent(in) :: slice_direction, maxm, meqn, mwaves, mbc, mx
    double precision, dimension(1-mbc:maxm+mbc,meqn), intent(in) :: q_left, q_right
    double precision, dimension(1-mbc:maxm+mbc, *), intent(in) :: aux_left, aux_right

    ! Output
    double precision, dimension(1-mbc:maxm+mbc, meqn, mwaves), intent(out) :: waves
    double precision, dimension(1-mbc:maxm+mbc, mwaves), intent(out) :: wave_speeds
    double precision, dimension(1-mbc:maxm+mbc, meqn), intent(out) ::  &
        a_plus_dq, a_minus_dq

    ! Local
    integer :: i
    double precision :: velocity_left, velocity_right
    !<--- Variables <---


    !---> Film height --->
    do i=2-mbc,mx+mbc
        !---- Set velocity on either side of interface ----
        if (slice_direction==1) then
            velocity_left = aux_right(i-1,1)
            velocity_right = aux_left(i,1)
        else
            velocity_left = aux_right(i-1,2)
            velocity_right = aux_left(i,2)
        end if
        
        !---- This waves doesn't affect the surfactant ----
        waves(i,2,1) = 0.d0

        !---- Speed and fluctuations ----
        if ((velocity_left>=0.d0) .and. (velocity_right>=0.d0)) then
            !-------------------------------------------------------------------
            ! All motion is to the right, so only the right-moving fluctuation
            ! a_plus_dq(i) is nonzero.  
            !
            ! The wave_speeds is determined by the Rankine-Hugoniot condition applied
            ! to the flux function
            !       f(q,x) = 1/2 * u(x) * q**2
            !-------------------------------------------------------------------
            waves(i,1,1) = velocity_right*q_left(i,1)**2/2.d0 -  &
                           velocity_left*q_right(i-1,1)**2/2.d0
            wave_speeds(i,1) = .5d0 * (velocity_right*q_left(i,1) +  &
                               sqrt(velocity_right*velocity_left)*q_right(i-1,1))
            a_minus_dq(i,1) = 0.d0
            a_plus_dq(i,1) = waves(i,1,1)

        else if ((velocity_left<=0.d0) .and. (velocity_right<=0.d0)) then
            waves(i,1,1) = velocity_right*q_left(i,1)**2/2.d0 -  &
                           velocity_left*q_right(i-1,1)**2/2.d0
            wave_speeds(i,1) = .5d0 * (-sqrt(velocity_right*velocity_left)*q_left(i,1)  &
                               + velocity_left*q_right(i-1,1))
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
            wave_speeds(i,1) = 0.d0    !(velocity_left + velocity_right)/2.d0 * (q_left(i,1) + q_right(i-1,1))/2.d0
            a_minus_dq(i,1) = -velocity_left * q_right(i-1,1)**2/2.d0
            a_plus_dq(i,1) = velocity_right * q_left(i,1)**2/2.d0
        end if
    end do
    !<--- Film height <---


    !---> Surfactant concentration --->
    do i = 2-mbc, mx+mbc
        !---- Select velocity for this cell interface ----          
        if (slice_direction==1) then
            velocity_left = aux_right(i-1,3)
            velocity_right = aux_left(i,3)
        else
            velocity_left = aux_right(i-1,4)
            velocity_right = aux_left(i,4)
        end if

        !---- This waves doesn't affect the film height ----
        waves(i,1,2) = 0.d0
        
        !---- Speed and fluctuations ----
        if ((velocity_left>=0.d0) .and. (velocity_right>=0.d0)) then
            !------------------------------------------------------------------
            ! All motion is to the right, so only the right-moving fluctuation
            ! a_plus_dq(i) is nonzero.  The wave_speeds is determined by the lower cell,
            ! as its mass is moving into the interface.
            !------------------------------------------------------------------
            waves(i,2,2) = velocity_right*q_left(i,2) - velocity_left*q_right(i-1,2) 
            wave_speeds(i,2) = velocity_right
            a_minus_dq(i,2)  = 0.d0
            a_plus_dq(i,2)   = waves(i,2,2)
            
        else if ((velocity_left<=0.d0) .and. (velocity_right<=0.d0)) then
            !-----------------------------------------------------------------
            ! All motion is to the left, so only the left-moving fluctuation
            ! a_minus_dq(i) is nonzero.  The wave_speeds is determined by the upper cell,
            ! as its mass is moving into the interface.
            !-----------------------------------------------------------------
            waves(i,2,2) = velocity_right*q_left(i,2) - velocity_left*q_right(i-1,2)
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
            a_minus_dq(i,2) = -velocity_left * q_right(i-1,2)
            a_plus_dq(i,2) = velocity_right * q_left(i,2)
        end if
    end do
    !<--- Surfactant concentration <---



    return
end subroutine rpn2