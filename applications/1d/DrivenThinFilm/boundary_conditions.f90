subroutine apply_homogeneous_bcs(q)

! Apply the relevant homogeneous boundary condition operator.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision, dimension(1-mbc:mx+mbc, meqn), intent(inout) :: q

    call extend_to_ghost_cells(mx, mbc, 'even', q(:,1))

end subroutine apply_homogeneous_bcs


subroutine apply_bcs(t, q)

! Fill ghost cells with appropriate boundary values.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision :: h_left, h_right
    common /boundary_config/ h_left, h_right

    double precision, intent(in) :: t
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(inout) :: q

    double precision :: boundary_value, delta
    integer :: i

    call extend_to_ghost_cells(mx, mbc, 'even', q(:, 1))

    ! ! Left boundary
    ! delta = 2.d0 * (h_left - q(1, 1))
    ! do i = 1, mbc
    !     q(1 - i, 1) = q(i, 1) + i * delta
    ! end do

    ! ! Right boundary
    ! delta = 2.d0 * (h_right - q(mx, 1))
    ! do i = 1, mbc
    !     q(mx + i, 1) = q(mx, 1) + i * delta
    ! end do

end subroutine apply_bcs


subroutine extend_to_ghost_cells (mx, mbc,  &
                                  extension_type,  &
                                  q)
                                  
!---------------------------------------------------------------
! Extends the interior data of q into its ghost cells, based on
! either even, odd, or periodic or odd extension.
!---------------------------------------------------------------

    implicit none
    
    !---> Variables --->
    ! Boilerplate input
    integer, intent(in) :: mx, mbc

    ! Input
    character*1 :: extension_type
    
    ! In/out
    double precision, intent(inout) :: q(1-mbc:mx+mbc)

    ! Local
    integer :: ix
    logical :: lsame
    !<--- Variables <---


    if (lsame(extension_type, 'e')) then
        ! Even extension.
        q(1-mbc:0) = q(mbc:1:-1)
        q(mx+1:mx+mbc) = q(mx:mx-mbc+1:-1)

    else if (lsame(extension_type, 'o')) then
        ! Odd extension.
        q(1-mbc:0) = -q(mbc:1:-1)
        q(mx+1:mx+mbc) = -q(mx:mx-mbc+1:-1)
        
    else if (lsame(extension_type, 'p')) then
        ! Periodic extension.
        q(1-mbc:0) = q(mx-mbc+1:mx)
        q(mx+1:mx+mbc) = q(1:mbc)
        
    else
        print *, "Error (extend_to_ghost_cells): Provided invalid extension type."
        stop
    end if
    
end subroutine extend_to_ghost_cells
