subroutine apply_homogeneous_bcs(q)

! Apply the relevant homogeneous boundary condition operator.

    implicit none

    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn

    character(len=2), dimension(4, 10) :: bc_options
    common /bc_config/ bc_options

    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn), intent(inout) :: q
    ! Dimensions are correct; there are my horizontal strips and mx vertical strips.
    double precision, dimension(2, my) :: x_zeros
    double precision, dimension(2, mx) :: y_zeros
    integer :: i


    x_zeros = 0.d0
    y_zeros = 0.d0
    do i = 1, meqn
        call fill_ghost_cells(bc_options(:, i), x_zeros, x_zeros,  &
                              y_zeros, y_zeros, q(:, :, i))
    end do

end subroutine apply_homogeneous_bcs


subroutine apply_bcs(t, q)

! Fill ghost cells with appropriate boundary values.

    implicit none

    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn

    character(len=2), dimension(4, 10) :: bc_options
    common /bc_config/ bc_options

    double precision, intent(in) :: t
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn), intent(inout) :: q

    double precision, dimension(2, my, meqn) :: x_lower_values, x_upper_values
    double precision, dimension(2, mx, meqn) :: y_lower_values, y_upper_values
    integer :: i
    

    call set_implicit_boundary_data(t, x_lower_values, x_upper_values, y_lower_values, y_upper_values)

    do i = 1, meqn
        call fill_ghost_cells(bc_options(:, i), x_lower_values, x_upper_values,  &
                              y_lower_values, y_upper_values, q(:, :, i))
    end do

end subroutine apply_bcs


subroutine fill_ghost_cells(bc_options, x_lower_values, x_upper_values,  &
                            y_lower_values, y_upper_values, q)

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

    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn

    character(len=2), dimension(4) :: bc_options
    double precision, intent(in), dimension(2, my) :: x_lower_values, x_upper_values
    double precision, intent(in), dimension(2, mx) :: y_lower_values, y_upper_values
    double precision, intent(inout) :: q(1-mbc:mx+mbc, 1-mbc:my+mbc)

    double precision, dimension(2,2) ::boundary_coeffs, cell_coeffs
    integer :: num_cells, ix, iy


    if (bc_options(1) == 'p' .or. bc_options(2) == 'p') then
        if (bc_options(1) /= bc_options(2)) then
            print *, "Error (fill_ghost_cells): If one option is 'p', ",  &
                "then both must be."
            stop
        end if

        do iy = 1, my
            q(1-mbc:0, iy) = q(mx+1-mbc:mx, iy)
            q(mx+1:mx+mbc, iy) = q(1:mbc, iy)
        end do
    else
        call get_lower_bc_coefficients(bc_options(1), dx, boundary_coeffs, cell_coeffs)
        num_cells = len_trim(bc_options(1))
        if (num_cells == 1) then
            do iy = 1, my
                q(0, iy) = boundary_coeffs(1, 1) * x_lower_values(1, iy)  &
                    + cell_coeffs(1, 1) * q(1, iy) + cell_coeffs(2, 1) * q(2, iy)
            end do
        else if (num_cells == 2) then
            do iy = 1, my
                q(-1, iy) = boundary_coeffs(1, 1) * x_lower_values(1, iy) + boundary_coeffs(2, 1) * x_lower_values(2, iy)  &
                    + cell_coeffs(1, 1) * q(1, iy) + cell_coeffs(2, 1) * q(2, iy)
                    q(0, iy) = boundary_coeffs(1, 2) * x_lower_values(1, iy) + boundary_coeffs(2, 2) * x_lower_values(2, iy)  &
                    + cell_coeffs(1, 2) * q(1, iy) + cell_coeffs(2, 2) * q(2, iy)
            end do
        end if
        
        call get_upper_bc_coefficients(bc_options(2), dx, boundary_coeffs, cell_coeffs)
        num_cells = len_trim(bc_options(2))
        if (num_cells == 1) then
            do iy = 1, mx
                q(mx+1, iy) = boundary_coeffs(1, 1) * x_upper_values(1, iy)  &
                    + cell_coeffs(1, 1) * q(mx-1, iy) + cell_coeffs(2, 1) * q(mx, iy)
            end do
        else if (num_cells == 2) then
            do iy = 1, my
                q(mx+1, iy) = boundary_coeffs(1, 1) * x_upper_values(1, iy) + boundary_coeffs(2, 1) * x_upper_values(2, iy)  &
                    + cell_coeffs(1, 1) * q(mx-1, iy) + cell_coeffs(2, 1) * q(mx, iy)
                q(mx+2, iy) = boundary_coeffs(1, 2) * x_upper_values(1, iy) + boundary_coeffs(2, 2) * x_upper_values(2, iy)  &
                    + cell_coeffs(1, 2) * q(mx-1, iy) + cell_coeffs(2, 2) * q(mx, iy)
            end do
        end if
    end if
    

    if (bc_options(3) == 'p' .or. bc_options(4) == 'p') then
        if (bc_options(3) /= bc_options(4)) then
            print *, "Error (fill_ghost_cells): If one option is 'p', ",  &
                "then both must be."
            stop
        end if

        do ix = 1, mx
            q(ix, 1-mbc:0) = q(ix, my+1-mbc:my)
            q(my+1:my+mbc, ix) = q(1:mbc, ix)
        end do
    else
        call get_lower_bc_coefficients(bc_options(3), dy, boundary_coeffs, cell_coeffs)
        num_cells = len_trim(bc_options(3))
        if (num_cells == 1) then
            do ix = 1, mx
                q(ix, 0) = boundary_coeffs(1, 1) * y_lower_values(1, ix)  &
                    + cell_coeffs(1, 1) * q(ix, 1) + cell_coeffs(2, 1) * q(ix, 2)
            end do
        else if (num_cells == 2) then
            do ix = 1, mx
                q(ix, -1) = boundary_coeffs(1, 1) * y_lower_values(1, ix) + boundary_coeffs(2, 1) * y_lower_values(2, ix)  &
                    + cell_coeffs(1, 1) * q(ix, 1) + cell_coeffs(2, 1) * q(ix, 2)
                q(ix, 0) = boundary_coeffs(1, 2) * y_lower_values(1, ix) + boundary_coeffs(2, 2) * y_lower_values(2, ix)  &
                    + cell_coeffs(1, 2) * q(ix, 1) + cell_coeffs(2, 2) * q(ix, 2)
            end do
        end if
        
        call get_upper_bc_coefficients(bc_options(4), dy, boundary_coeffs, cell_coeffs)
        num_cells = len_trim(bc_options(4))
        if (num_cells == 1) then
            do ix = 1, mx
                q(ix, my+1) = boundary_coeffs(1, 1) * y_upper_values(1, ix)  &
                    + cell_coeffs(1, 1) * q(ix, my-1) + cell_coeffs(2, 1) * q(ix, my)
            end do
        else if (num_cells == 2) then
            do ix = 1, mx
                q(ix, my+1) = boundary_coeffs(1, 1) * y_upper_values(1, ix) + boundary_coeffs(2, 1) * y_upper_values(2, ix)  &
                    + cell_coeffs(1, 1) * q(ix, my-1) + cell_coeffs(2, 1) * q(ix, my)
                q(ix, my+2) = boundary_coeffs(1, 2) * y_upper_values(1, ix) + boundary_coeffs(2, 2) * y_upper_values(2, ix)  &
                    + cell_coeffs(1, 2) * q(ix, my-1) + cell_coeffs(2, 2) * q(ix, my)
            end do
        end if
    end if
    
    ! TODO: Handle corner ghost cells
    
end subroutine fill_ghost_cells