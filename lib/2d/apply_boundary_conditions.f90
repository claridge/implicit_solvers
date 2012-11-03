subroutine apply_homogeneous_bcs(q)

! Apply the relevant homogeneous boundary condition operator.

    implicit none

    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn

    character(len=2), dimension(4, 10) :: bc_options
    common /bc_config/ bc_options

    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn), intent(inout) :: q
    integer :: i

    do i = 1, meqn
        call fill_ghost_cells_homogeneous(bc_options(:, i), q(:,:,i))
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

    call set_implicit_boundary_data(t, x_lower_values, x_upper_values,  &
                                    y_lower_values, y_upper_values)

    do i = 1, meqn
        call fill_ghost_cells(bc_options(:, i), x_lower_values(:, :, i), x_upper_values(:, :, i),  &
                              y_lower_values(:, :, i), y_upper_values(:, :, i), q(:, :, i))
    end do

end subroutine apply_bcs
