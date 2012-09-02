subroutine fill_2_ghost_cells(lower_option, lower_values,  &
                              upper_option, upper_values,  &
                              q)

! Fills 1 ghost cell on either side of the domain.
!
! 'option' arguments indicate the condition used to fill the cell.
! Allowable values are:
!   'p': periodic
!
! TODO: Specify allowable values.
!
! If either option is 'p', then both must be.
!
! 'value' arguments are double precision, indicating the value to be
! matched.  If the corresponding option is 'p', this value is unused.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    character*2, intent(in) :: lower_option, upper_option
    double precision, intent(in), dimension(2) :: lower_values, upper_values
    double precision, intent(inout) :: q(1-mbc:mx+mbc)

    integer, dimension(2) :: derivative_orders, ipiv
    integer :: info, i
    double precision, dimension(2,2) ::boundary_coeffs, cell_coeffs

    if (lower_option(1:1) == 'p' .or. upper_option(1:1) == 'p') then
        if (lower_option(1:1) /= upper_option(1:1)) then
            print *, "Error (fill_2_ghost_cells): If one option is 'p', ",  &
                "then both must be."
            stop
        end if

        q(-1:0) = q(mx-1:mx)
        q(mx+1:mx+2) = q(1:2)
        return
    end if


    call get_lower_bc_coefficients(lower_option, dx, boundary_coeffs, cell_coeffs)
    q(-1) = boundary_coeffs(1, 1) * lower_values(1) + boundary_coeffs(2, 1) * lower_values(2)  &
        + cell_coeffs(1, 1) * q(1) + cell_coeffs(2, 1) * q(2)
    q(0) = boundary_coeffs(1, 2) * lower_values(1) + boundary_coeffs(2, 2) * lower_values(2)  &
        + cell_coeffs(1, 2) * q(1) + cell_coeffs(2, 2) * q(2)
    

    call get_upper_bc_coefficients(upper_option, dx, boundary_coeffs, cell_coeffs)
    q(mx+1) = boundary_coeffs(1, 1) * upper_values(1) + boundary_coeffs(2, 1) * upper_values(2)  &
        + cell_coeffs(1, 1) * q(mx-1) + cell_coeffs(2, 1) * q(mx)
    q(mx+2) = boundary_coeffs(1, 2) * upper_values(1) + boundary_coeffs(2, 2) * upper_values(2)  &
        + cell_coeffs(1, 2) * q(mx-1) + cell_coeffs(2, 2) * q(mx)
    
end subroutine fill_2_ghost_cells
