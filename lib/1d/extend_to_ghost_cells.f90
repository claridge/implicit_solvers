subroutine extend_to_ghost_cells (lower_extension, upper_extension, q)

! Extends the interior data of q into its ghost cells.  Possible extension types
! are "even", "odd", "periodic", or "none".  (Only the first letter of the
! argument actually matters.)
!
! Args:
!   lower_extension: Extension type, as above.
!   upper_extension: Extension type, as above.
!   q: Array of values being extended.


    implicit none
    
    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    character*1, intent(in) :: lower_extension, upper_extension
    double precision, intent(inout) :: q(1-mbc:mx+mbc)

    integer :: ix
    logical :: lsame


    if (lower_extension == "e") then
        q(1-mbc:0) = q(mbc:1:-1)
    else if (lower_extension == "o") then
        q(1-mbc:0) = -q(mbc:1:-1)
    else if (lower_extension == "p") then
        q(1-mbc:0) = q(mx-mbc+1:mx)
    else if (lower_extension == "n") then
        ! Nothing
    else
        print *, "extend_to_ghost_cells: Error: ",  &
            "Invalid value for lower_extnsion_type."
    end if

    if (upper_extension == "e") then
        q(mx+1:mx+mbc) = q(mx:mx-mbc+1:-1)
    else if (upper_extension == "o") then
        q(mx+1:mx+mbc) = -q(mx:mx-mbc+1:-1)
    else if (upper_extension == "p") then
        q(mx+1:mx+mbc) = q(1:mbc)
    else if (upper_extension == "n") then
        ! Nothing
    else
        print *, "extend_to_ghost_cells: Error: ",  &
            "Invalid value for upper_extnsion_type."
    end if
    
end subroutine extend_to_ghost_cells
