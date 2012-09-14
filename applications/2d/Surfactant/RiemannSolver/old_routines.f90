subroutine surface_tension_d1(tension_d1, mu, mx)
    implicit none
    integer, intent(in) :: mx
    double precision, intent(in) :: mu
    double precision, intent(inout) :: tension_d1(mx)
    integer :: i
    
    do i=1,mx
        tension_d1(i) = -3.d0*mu*(1.d0+mu*tension_d1(i))**(-4)
    end do
end subroutine surface_tension_d1


subroutine film_apply_bcs(mx, my, mbc,  &
                          film)
    implicit none
    integer, intent(in) :: mx, my, mbc
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc), intent(out) :: film
    call extend_to_ghost_cells(mx, my, mbc,  &
                               'even', 'periodic',  &
                               film)
end subroutine film_apply_bcs

subroutine film_x_velocity_apply_bcs(mx, my, mbc,  &
                                     velocity)
    implicit none
    integer, intent(in) :: mx, my, mbc
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc), intent(out) :: velocity
    call extend_to_ghost_cells(mx, my, mbc,  &
                               'odd', 'periodic',  &
                               velocity)
end subroutine film_x_velocity_apply_bcs

subroutine film_y_velocity_apply_bcs(mx, my, mbc,  &
                                     velocity)
    implicit none
    integer, intent(in) :: mx, my, mbc
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc), intent(out) :: velocity
    call extend_to_ghost_cells(mx, my, mbc,  &
                               'even', 'odd',  &
                               velocity)
end subroutine film_y_velocity_apply_bcs

subroutine film_laplacian_apply_bcs(mx, my, mbc,  &
                                    film_laplacian)
    implicit none
    integer, intent(in) :: mx, my, mbc
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc), intent(out) ::  &
        film_laplacian
    call extend_to_ghost_cells(mx, my, mbc,  &
                               'even', 'even',  &
                               film_laplacian)
end subroutine film_laplacian_apply_bcs

subroutine surfactant_apply_bcs(mx, my, mbc,  &
                                surfactant)
    implicit none    
    integer, intent(in) :: mx, my, mbc
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc), intent(out) :: surfactant
    call extend_to_ghost_cells(mx, my, mbc,  &
                               'even', 'even',  &
                               surfactant)
end subroutine surfactant_apply_bcs

subroutine surfactant_x_velocity_apply_bcs(mx, my, mbc,  &
                                           velocity)
    implicit none
    integer, intent(in) :: mx, my, mbc
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc), intent(out) :: velocity
    call extend_to_ghost_cells(mx, my, mbc,  &
                               'odd', 'even',  &
                               velocity)
end subroutine surfactant_x_velocity_apply_bcs

subroutine surfactant_y_velocity_apply_bcs(mx, my, mbc,  &
                                           velocity)
    implicit none
    integer, intent(in) :: mx, my, mbc
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc), intent(out) :: velocity
    call extend_to_ghost_cells(mx, my, mbc,  &
                               'even', 'odd',  &
                               velocity)
end subroutine surfactant_y_velocity_apply_bcs


subroutine extend_to_ghost_cells (mx, my, mbc,  &
                                  x_extension_type, y_extension_type,  &
                                  q)
                                  
!--------------------------------------------------------------------
! Extends the interior data of q into its ghost cells, based on
! either even, odd, or periodic or odd extension in both directions.
!--------------------------------------------------------------------

    !$ use omp_lib
    implicit none
    
    !---> Variables --->
    ! Input
    integer, intent(in) :: mx, my, mbc
    character*1 :: x_extension_type, y_extension_type
    
    ! In/out
    double precision, intent(inout) :: q(1-mbc:mx+mbc, 1-mbc:my+mbc)

    ! Local
    integer :: ix, iy, iy_target
    logical :: lsame
    !$ integer :: nthreads
    !$ common /omp_block/ nthreads
    !<--- Variables <---

    !$ call omp_set_num_threads(nthreads)

    !---> Perform x-extension --->
    if ( lsame(x_extension_type, 'e') ) then  ! Even extension
        !$omp parallel do
        do iy = 1, my
            q(1-mbc:0,iy)     = q(mbc:1:-1,      iy)
            q(mx+1:mx+mbc,iy) = q(mx:mx-mbc+1:-1,iy)
        end do

    else if ( lsame(x_extension_type, 'o') ) then  ! Odd extension
        !$omp parallel do
        do iy = 1, my       
            q(1-mbc:0,iy)     = -q(mbc:1:-1,      iy)
            q(mx+1:mx+mbc,iy) = -q(mx:mx-mbc+1:-1,iy)
        end do
        
    else if ( lsame(x_extension_type, 'p') ) then  ! Periodic extension
        !$omp parallel do
        do iy = 1, my
            q(1-mbc:0,     iy) = q(mx-mbc+1:mx, iy)
            q(mx+1:mx+mbc, iy) = q(1:mbc,       iy)
        end do
        
    else
        print *, "Error (extend_to_ghost_cells): Provided invalid x extension type."
        stop
    end if
    !<--- Perform x-extension <---
        

    !---> Perform y-extension --->
    if ( lsame(y_extension_type, 'e') ) then  ! Even extension
        do iy = 1-mbc, 0
            iy_target = 1-iy
            !$omp parallel do
            do ix = 1-mbc, mx+mbc
                q(ix, iy) = q(ix, iy_target)
            end do
        end do
        
        do iy = my+1, my+mbc
            iy_target = my - (iy-my) + 1            
            !$omp parallel do
            do ix = 1-mbc, mx+mbc
                q(ix, iy) = q(ix, iy_target)
            end do
        end do
    
    else if ( lsame(y_extension_type, 'o') ) then  ! Odd extension
        do iy = 1-mbc, 0
            iy_target = 1-iy
            !$omp parallel do
            do ix = 1-mbc, mx+mbc
                q(ix, iy) = -q(ix, iy_target)
            end do
        end do
        
        do iy = my+1, my+mbc
            iy_target = my - (iy-my) + 1
            !$omp parallel do
            do ix = 1-mbc, mx+mbc
                q(ix, iy) = -q(ix, iy_target)
            end do
        end do
    
    else if ( lsame(y_extension_type, 'p') ) then  ! Periodic extension
        do iy = 1-mbc, 0
            iy_target = iy + my
            !$omp parallel do
            do ix = 1-mbc, mx+mbc
                q(ix, iy) = q(ix, iy_target)
            end do
        end do
        
        do iy = my+1, my+mbc
            iy_target = iy - my
            !$omp parallel do
            do ix = 1-mbc, mx+mbc
                q(ix, iy) = q(ix, iy_target)
            end do
        end do
            
    else
        print *, "Error (extend_to_ghost_cells): Provided invalid y extension type."
        stop
    end if
    !<--- y-extension <---
    
end subroutine extend_to_ghost_cells