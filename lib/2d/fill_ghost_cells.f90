subroutine fill_ghost_cells_homogeneous(bc_options, q)

! Apply the relevant homogeneous boundary condition operator.

    implicit none

    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn

    character(len=2), dimension(4), intent(in) :: bc_options
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn), intent(inout) :: q
    double precision, dimension(:, :), allocatable :: zeros
    integer :: i

    allocate(zeros(2, max(mx, my)))
    zeros = 0.d0

    do i = 1, meqn
        ! Dimensions are correct; there are my horizontal strips and mx vertical strips.
        call fill_ghost_cells(bc_options, zeros(:, 1:my), zeros(:, 1:my),  &
                              zeros(:, 1:mx), zeros(:, 1:mx), q(:, :, i))
    end do

end subroutine fill_ghost_cells_homogeneous



! TODO: Move to separate file, and add dedicated homogeneous method for
! convenience.
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

    character(len=2), dimension(4), intent(in) :: bc_options
    double precision, intent(in), dimension(2, my) :: x_lower_values, x_upper_values
    double precision, intent(in), dimension(2, mx) :: y_lower_values, y_upper_values
    double precision, intent(inout) :: q(1-mbc:mx+mbc, 1-mbc:my+mbc)

    ! TODO: Clean up unused variables
    double precision, dimension(2, 2) :: boundary_coeffs, cell_coeffs
    double precision, dimension(2, 2) :: corner
    integer :: num_cells, ix, iy, i
    double precision, dimension(2, mbc) :: zeros
    double precision, dimension(2, 1) :: extrap_boundary_values
    double precision :: temp_corner


    zeros = 0.d0

    ! TODO: Make sure there's a safety check on the periodic option.
    call set_x_lower_ghost_cells(bc_options(1), x_lower_values, 1, my)
    call set_x_upper_ghost_cells(bc_options(2), x_upper_values, 1, my)
    call set_y_lower_ghost_cells(bc_options(3), y_lower_values, 1, mx)
    call set_y_upper_ghost_cells(bc_options(4), y_upper_values, 1, mx)

    ! Handle corners if one or both dimensions are periodic.
    if (bc_options(1) == 'p') then
        call set_x_lower_ghost_cells(bc_options(1), zeros, 1-mbc, 0)
        call set_x_lower_ghost_cells(bc_options(1), zeros, my+1, my+mbc)
        call set_x_upper_ghost_cells(bc_options(2), zeros, 1-mbc, 0)
        call set_x_upper_ghost_cells(bc_options(2), zeros, my+1, my+mbc)
    else if (bc_options(3) == 'p') then
        call set_y_lower_ghost_cells(bc_options(3), zeros, 1-mbc, 0)
        call set_y_lower_ghost_cells(bc_options(3), zeros, mx+1, mx+mbc)
        call set_y_upper_ghost_cells(bc_options(4), zeros, 1-mbc, 0)
        call set_y_upper_ghost_cells(bc_options(4), zeros, mx+1, mx+mbc)
    else
        q(0, 0) = (quartic_extrap(q(1:5, 0)) + quartic_extrap(q(0, 1:5))) / 2.d0
        q(mx+1, 0) = (quartic_extrap(q(mx:mx-4:-1, 0)) + quartic_extrap(q(mx+1, 1:5))) / 2.d0
        q(0, my+1) = (quartic_extrap(q(1:5, my+1)) + quartic_extrap(q(0, my:my-4:-1))) / 2.d0
        q(mx+1, my+1) = (quartic_extrap(q(mx:mx-4:-1, my+1)) + quartic_extrap(q(mx+1, my:my-4:-1))) / 2.d0        
    ! TODO: Clear this out when I'm sure I don't want it.
    !         do i = 1, 2
    !             extrap_boundary_values(i, 1) = cubic_extrap(x_lower_values(i, 1:4))
    !         end do
    !         call set_x_lower_ghost_cells(bc_options(1), extrap_boundary_values, 0, 0)
    !         do i = 1, 2
    !             extrap_boundary_values(i, 1) = cubic_extrap(y_lower_values(i, 1:4))
    !         end do
    !         call set_y_lower_ghost_cells(bc_options(3), extrap_boundary_values, 0, 0)
    !         q(0, 0) = (q(0, 0) + temp_corner) / 2.d0
    end if


    contains


    double precision function cubic_extrap(u)
        double precision, dimension(4) :: u
        cubic_extrap = 4.d0 * u(1) - 6.d0 * u(2) + 4.d0 * u(3) - u(4)
    end function cubic_extrap

    double precision function quartic_extrap(u)
        double precision, dimension(5) :: u
        quartic_extrap = 5.d0 * u(1) - 10.d0 * u(2) + 10.d0 * u(3) - 5.d0 * u(4) + u(5)
    end function quartic_extrap    
    

    subroutine set_x_lower_ghost_cells(options, boundary_values, iy_low, iy_high)
        implicit none
        character(len=2), intent(in) :: options
        integer, intent(in) :: iy_low, iy_high
        double precision, dimension(2, iy_high-iy_low+1), intent(in) :: boundary_values
        double precision, dimension(2, 2) :: boundary_coeffs, cell_ceffs

        if (options == 'p') then
            do iy = iy_low, iy_high
                q(1-mbc:0, iy) = q(mx+1-mbc:mx, iy)
            end do
        else
            call get_lower_bc_coefficients(options, dx, boundary_coeffs, cell_coeffs)
            num_cells = len_trim(options)
            if (num_cells == 1) then
                do iy = iy_low, iy_high
                    q(0, iy) = boundary_coeffs(1, 1) * boundary_values(1, iy) + dot_product(q(1:2, iy), cell_coeffs(1:2, 1))
                end do
            else if (num_cells == 2) then
                do iy = iy_low, iy_high
                    q(-1:0, iy) = matmul(boundary_values(:, iy), boundary_coeffs) + matmul(q(1:2, iy), cell_coeffs)
                end do
            end if
        end if
    end subroutine set_x_lower_ghost_cells


    subroutine set_x_upper_ghost_cells(options, boundary_values, iy_low, iy_high)
        implicit none
        character(len=2), intent(in) :: options
        integer, intent(in) :: iy_low, iy_high
        double precision, dimension(2, iy_high-iy_low+1), intent(in) :: boundary_values
        double precision, dimension(2, 2) :: boundary_coeffs, cell_ceffs

        if (options == 'p') then
            do iy = iy_low, iy_high
                q(mx+1:mx+mbc, iy) = q(1:mbc, iy)
            end do
        else
            call get_upper_bc_coefficients(options, dx, boundary_coeffs, cell_coeffs)
            num_cells = len_trim(options)
            if (num_cells == 1) then
                do iy = iy_low, iy_high
                    q(mx+1, iy) = boundary_coeffs(1, 1) * boundary_values(1, iy)  &
                        + dot_product(q(mx-1:mx, iy), cell_coeffs(:, 1))
                end do
            else if (num_cells == 2) then
                do iy = iy_low, iy_high
                    q(mx+1:mx+2, iy) = matmul(boundary_values(:, iy), boundary_coeffs)  &
                        + matmul(q(mx-1:mx, iy), cell_coeffs)
                end do
            end if
        end if
    end subroutine set_x_upper_ghost_cells


    subroutine set_y_lower_ghost_cells(options, boundary_values, ix_low, ix_high)
        implicit none
        character(len=2), intent(in) :: options
        integer, intent(in) :: ix_low, ix_high
        double precision, dimension(2, ix_high-ix_low+1), intent(in) :: boundary_values
        double precision, dimension(2, 2) :: boundary_coeffs, cell_ceffs

        if (options == 'p') then
            do ix = ix_low, ix_high
                q(ix, 1-mbc:0) = q(ix, my+1-mbc:my)
            end do
        else
            call get_lower_bc_coefficients(bc_options(3), dy, boundary_coeffs, cell_coeffs)
            num_cells = len_trim(bc_options(3))
            if (num_cells == 1) then
                do ix = ix_low, ix_high
                    q(ix, 0) = boundary_coeffs(1, 1) * boundary_values(1, ix)  &
                        + dot_product(q(ix, 1:2), cell_coeffs(:, 1))
                end do
            else if (num_cells == 2) then
                do ix = ix_low, ix_high
                    q(ix, -1:0) = matmul(boundary_values(:, ix), boundary_coeffs) + matmul(q(ix, 1:2), cell_coeffs)
                end do
            end if
        end if
    end subroutine set_y_lower_ghost_cells
    

    subroutine set_y_upper_ghost_cells(options, boundary_values, ix_low, ix_high)
        implicit none
        character(len=2), intent(in) :: options
        integer, intent(in) :: ix_low, ix_high
        double precision, dimension(2, ix_high-ix_low+1), intent(in) :: boundary_values
        double precision, dimension(2, 2) :: boundary_coeffs, cell_ceffs

        if (options == 'p') then
            do ix = ix_low, ix_high
                q(ix, my+1:my+mbc) = q(ix, 1:mbc)
            end do
        else
            call get_upper_bc_coefficients(bc_options(4), dy, boundary_coeffs, cell_coeffs)
            num_cells = len_trim(bc_options(4))
            if (num_cells == 1) then
                do ix = 1, mx
                    q(ix, my+1) = boundary_coeffs(1, 1) * boundary_values(1, ix)  &
                        + dot_product(q(ix, my-1:my), cell_coeffs(:, 1))
                end do
            else if (num_cells == 2) then
                do ix = 1, mx
                    q(ix, my+1:my+2) = matmul(boundary_values(:, ix), boundary_coeffs) + matmul(q(ix, my-1:my), cell_coeffs)
                end do
            end if
        end if
    end subroutine set_y_upper_ghost_cells

end subroutine fill_ghost_cells