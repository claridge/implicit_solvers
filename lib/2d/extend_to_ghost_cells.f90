subroutine extend_to_ghost_cells(x_lower_extension, x_upper_extension,  &
                                 y_lower_extension, y_upper_extension,  &
                                 q)
                                  
! Extends the interior data of q into its ghost cells.  Possible extension types
! are "even", "odd", "periodic", or "none".  (Only the first letter of the
! argument actually matters.)
!
! Corner ghost cells are filled by first extending in the x-direction only for
! interior y-coordinates, and then extending in the y-direction for both
! interior and exterior x-coordinates.
!
! Unlike the other library routines, this acts only on a single component of the
! PDE solution.
!
! Args:
!   x_lower_extension: Extension type, as above.
!   x_upper_extension: Extension type, as above.
!   y_lower_extension: Extension type, as above.
!   y_upper_extension: Extension type, as above.
!   q: Array of values being extended.

    !$ use omp_lib
    implicit none
    
    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn

    character*1, intent(in) :: x_lower_extension, x_upper_extension,  &
        y_lower_extension, y_upper_extension
    double precision, intent(inout) :: q(1-mbc:mx+mbc, 1-mbc:my+mbc)

    integer :: ix, iy, iy_target
    logical :: lsame


    ! TODO: What about corners?

    if (x_lower_extension == "e") then
        !$omp parallel do
        do iy = 1, my
            q(1-mbc:0, iy) = q(mbc:1:-1, iy)
        end do
    else if (x_lower_extension == "o") then
        !$omp parallel do
        do iy = 1, my
            q(1-mbc:0, iy) = -q(mbc:1:-1, iy)
        end do
    else if (x_lower_extension == "p") then
        !$omp parallel do
        do iy = 1, my
            q(1-mbc:0, iy) = q(mx-mbc+1:mx, iy)
        end do
    else if (x_lower_extension == "n") then
        ! Nothing
    else
        print *, "extend_to_ghost_cells: Error: ",  &
            "Invalid value for x_lower_extension."
        stop
    end if


    if (x_upper_extension == "e") then
        !$omp parallel do
        do iy = 1, my
            q(mx+1:mx+mbc, iy) = q(mx:mx-mbc+1:-1, iy)
        end do
    else if (x_upper_extension == "o") then
        !$omp parallel do
        do iy = 1, my
            q(mx+1:mx+mbc, iy) = -q(mx:mx-mbc+1:-1, iy)
        end do
    else if (x_upper_extension == "p") then
        !$omp parallel do
        do iy = 1, my
            q(mx+1:mx+mbc, iy) = q(1:mbc, iy)
        end do
    else if (x_upper_extension == "n") then
        ! Nothing
    else
        print *, "extend_to_ghost_cells: Error: ",  &
            "Invalid value for x_upper_extension."
        stop
    end if
    

    ! y-extensions are implemented differently because a strip of ghost cells
    ! in the y-direction are not stored contiguously in memory.
    if (y_lower_extension == "e") then
        do iy = 1-mbc, 0
            iy_target = 1-iy
            !$omp parallel do
            do ix = 1-mbc, mx+mbc
                q(ix, iy) = q(ix, iy_target)
            end do
        end do    
    else if (y_lower_extension == "o") then
        do iy = 1-mbc, 0
            iy_target = 1-iy
            !$omp parallel do
            do ix = 1-mbc, mx+mbc
                q(ix, iy) = -q(ix, iy_target)
            end do
        end do
    else if (y_lower_extension == "p") then
        do iy = 1-mbc, 0
            iy_target = iy + my
            !$omp parallel do
            do ix = 1-mbc, mx+mbc
                q(ix, iy) = q(ix, iy_target)
            end do
        end do
    else if (y_lower_extension == "n") then
        ! Nothing
    else
        print *, "extend_to_ghost_cells: Error: ",  &
            "Invalid value for y_lower_extension."
        stop
    end if


    if (y_upper_extension == "e") then
        do iy = my+1, my+mbc
            iy_target = my - (iy-my) + 1            
            !$omp parallel do
            do ix = 1-mbc, mx+mbc
                q(ix, iy) = q(ix, iy_target)
            end do
        end do
    else if (y_upper_extension == "o") then
        do iy = my+1, my+mbc
            iy_target = my - (iy-my) + 1
            !$omp parallel do
            do ix = 1-mbc, mx+mbc
                q(ix, iy) = -q(ix, iy_target)
            end do
        end do
    else if (y_upper_extension == "p") then
        do iy = my+1, my+mbc
            iy_target = iy - my
            !$omp parallel do
            do ix = 1-mbc, mx+mbc
                q(ix, iy) = q(ix, iy_target)
            end do
        end do
    else if (y_upper_extension == "n") then
        ! Nothing
    else
        print *, "extend_to_ghost_cells: Error: ",  &
            "Invalid value for y_upper_extension."
        stop
    end if
    
end subroutine extend_to_ghost_cells
