subroutine fill_ghost_cells_homogeneous(bc_options, q_component)

! Fill ghost cells with homogeneous boundary conditions.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    character(len=2), dimension(2), intent(in) :: bc_options
    double precision, dimension(1-mbc:mx+mbc), intent(inout) :: q_component
    double precision, dimension(2) :: zeros
    integer :: i

    zeros = 0.d0

    call fill_ghost_cells(bc_options, zeros, zeros, q_component)

end subroutine fill_ghost_cells_homogeneous


subroutine fill_ghost_cells(bc_options, lower_values, upper_values, q)

! Fills ghost cells on either side of the domain for a single solution
! component.
!
! 'option' arguments indicate the condition used to fill the cell.
! Allowable values are:
!   'p': periodic
!   '0' or '1': Match derivative 0 or 1
!   '00', '01', '02', '03', '12', '13', or '23': Match derivatives of the
!     two specified orders.
!
! If either option is 'p', then both must be.
!
! 'value' arguments are double precision, indicating the value to be
! matched.  If the corresponding option is 'p', this value is unused.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    character(len=2), dimension(2) :: bc_options
    double precision, intent(in), dimension(2) :: lower_values, upper_values
    double precision, intent(inout) :: q(1-mbc:mx+mbc)

    double precision, dimension(2, 2) :: boundary_coeffs
    double precision, dimension(3, 2) :: cell_coeffs

    ! TODO: Make sure periodic options are consistent on initialization.

    ! if (bc_options(1) == 'p' .or. bc_options(2) == 'p') then
    !     if (bc_options(1) /= bc_options(2)) then
    !         print *, "Error (fill_ghost_cells): If one option is 'p', ",  &
    !             "then both must be."
    !         stop
    !     end if

    !     q(1-mbc:0) = q(mx+1-mbc:mx)
    !     q(mx+1:mx+mbc) = q(1:mbc)
    !     return
    ! end if

    call get_lower_bc_coefficients(bc_options(1), dx, boundary_coeffs, cell_coeffs)
    if (bc_options(1) == 'p') then
        q(1-mbc:0) = q(mx+1-mbc:mx)
    else if (len_trim(bc_options(1)) == 1) then
        q(-1:0) = lower_values(1) * boundary_coeffs(1, :) + matmul(q(1:3), cell_coeffs)
    else
        q(-1:0) = matmul(lower_values, boundary_coeffs) + matmul(q(1:3), cell_coeffs)
    end if

    call get_upper_bc_coefficients(bc_options(2), dx, boundary_coeffs, cell_coeffs)
    if (bc_options(2) == 'p') then
        q(mx+1:mx+mbc) = q(1:mbc)
    else if (len_trim(bc_options(2)) == 1) then
        q(mx+1:mx+2) = upper_values(1) * boundary_coeffs(1, :) + matmul(q(mx-2:mx), cell_coeffs)
    else
        q(mx+1:mx+2) = matmul(upper_values, boundary_coeffs) + matmul(q(mx-2:mx), cell_coeffs)
    end if
    
end subroutine fill_ghost_cells
