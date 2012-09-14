subroutine rpt2(slice_direction, maxm, meqn, mwaves, mbc, mx, ql, qr, &
                aux_below, aux_aligned, aux_above,        &
                imp, asdq, bmasdq, bpasdq)

!---------------------------------------------------------------
! Transverse Riemann solver
!---------------------------------------------------------------


    implicit none


    !---> Variables --->
    
    !---- Inputs ----
    integer, intent(in) :: slice_direction, maxm, meqn, mwaves, mbc, mx, imp    
    double precision, intent(in) ::     ql(1-mbc:maxm+mbc, meqn)
    double precision, intent(in) ::     qr(1-mbc:maxm+mbc, meqn)
    double precision, intent(in), dimension(1-mbc:maxm+mbc,*) :: aux_below, aux_aligned, aux_above
    double precision, intent(in) ::   asdq(1-mbc:maxm+mbc, meqn)

    !--- Outputs ----
    double precision, intent(out) :: bmasdq(1-mbc:maxm+mbc, meqn)
    double precision, intent(out) :: bpasdq(1-mbc:maxm+mbc, meqn)

    !---- Locals ----
    integer :: i, i_velocity, iaux_transverse, iaux_normal
    double precision :: s_transverse, u_upper, u_lower, q_wave
    double precision :: small_value
    parameter( small_value = 1.d-8 )

    !<--- Variables <---



    !---> Film height --->

    !---- Determine which part of the aux array to use ----
    
    if (slice_direction==1) then
        iaux_normal     = 1
        iaux_transverse = 2
    else
        iaux_normal     = 2
        iaux_transverse = 1
    end if
    
    
    do i = 2-mbc, mx+mbc
      
        !---- Choose the cell whose velocity we use ----
        !---------------------------------------------------------------
        ! If imp==1, then we're looking at the left-moving fluctuation,
        ! amdq.  Thus, we need the velocity in the left cell.
        ! If imp==2, we're looking at the right-moving fluctuation, so
        ! we need the velocity in the right cell.
        !---------------------------------------------------------------

        if (imp .eq. 1) then
            i_velocity = i-1
        else
            i_velocity = i
        end if


        !---> Approximate the transverse velocity --->
        !-----------------------------------------------------------------
        ! s_transverse approximates g'(q), the transverse speed of
        ! the wave.  This is based on velocity field in the cell from
        ! which the wave originates.  Note that the transverse velocities
        ! in the origin cell and the destination cell must agree in order
        ! for there to be any transverse motion.
        !-----------------------------------------------------------------
        
        !---- Approximate q on the wave ---
        u_upper = aux_aligned(i-1,iaux_normal)
        u_lower = aux_aligned(i,iaux_normal)
        
        
        !---- Checking the downwind velocity against small_value to prevent division by zero (JLC, 7/21/2011) ----
        
        if ( (u_upper>small_value) .and. (u_lower>0.d0) ) then
            q_wave = ( sqrt(u_lower/u_upper)*qr(i-1,1) + ql(i,1) ) / 2.d0
        else if ((u_upper<0.d0) .and. (u_lower)<-small_value) then
            q_wave = ( qr(i-1,1) + sqrt(u_upper/u_lower)*ql(i,1) ) / 2.d0
        else
            q_wave = 0.d0
        end if        
        
        
        
        if (aux_aligned(i_velocity,iaux_transverse) >= 0.d0) then
            !---- Transverse motion is up ----
          
            if (aux_above(i_velocity,iaux_transverse) > 0.d0) then
                s_transverse = aux_aligned(i_velocity,iaux_transverse) * q_wave
            else
                s_transverse = 0.d0
            end if
            
            bmasdq(i,1) = 0.d0
            bpasdq(i,1) = s_transverse * asdq(i,1)
        
        else
            !---- Transverse motion is down ----
            
            if (aux_below(i_velocity,iaux_transverse) < 0.d0) then
                s_transverse = aux_aligned(i_velocity,iaux_transverse) * q_wave
            else
                s_transverse = 0.d0
            end if
            
            bmasdq(i,1) = s_transverse * asdq(i,1)
            bpasdq(i,1) = 0.d0
        end if
        
        !<--- Approximate the transverse velocity <---

    end do
    
    !<--- Film height <---




    !---> Surfactant concentration --->

    !---- Determine which part of the aux array to use ----
    
    if (slice_direction==1) then
        iaux_normal     = 3
        iaux_transverse = 4
    else
        iaux_normal     = 4
        iaux_transverse = 3
    end if


    do i = 2-mbc, mx+mbc
      
        !---- Choose the cell whose velocity we use ----
        !---------------------------------------------------------------
        ! If imp==1, then we're looking at the left-moving fluctuation,
        ! amdq.  Thus, we need the velocity in the left cell.
        ! If imp==2, we're looking at the right-moving fluctuation, so
        ! we need the velocity in the right cell.
        !---------------------------------------------------------------

        if (imp .eq. 1) then
            i_velocity = i-1
        else
            i_velocity = i
        end if


        !---- Approximate the transverse velocity ----
        !------------------------------------------------------------
        ! In this case, g'(q) is simply the velocity.  (See the film
        ! calculation above for more details.)
        !------------------------------------------------------------

        if (aux_aligned(i_velocity,iaux_transverse) >= 0.d0) then
          
            if (aux_above(i_velocity,iaux_transverse) > 0.d0) then
                s_transverse = aux_aligned(i_velocity,iaux_transverse)
            else
                s_transverse = 0.d0
            end if
    
            bmasdq(i,2) = 0.d0
            bpasdq(i,2) = s_transverse * asdq(i,2)

        else
          
            if (aux_below(i_velocity,iaux_transverse) < 0.d0) then
                s_transverse = aux_aligned(i_velocity,iaux_transverse)
            else
                s_transverse = 0.d0
            end if
    
            bmasdq(i,2) = s_transverse * asdq(i,2)
            bpasdq(i,2) = 0.d0
            
        end if
        
    end do
    
    !<--- Surfactant concentration <---
    
    

    return
end subroutine rpt2
