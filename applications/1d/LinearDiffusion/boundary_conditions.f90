subroutine set_implicit_boundary_data(t, lower_values, upper_values)

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    character(len=2), dimension(2, 10) :: bc_options
    common /bc_config/ bc_options

    double precision, dimension(4, 0:3) :: stencils
    common /stencil_config/ stencils

    double precision, intent(in) :: t
    double precision, dimension(2), intent(out) :: lower_values, upper_values

    double precision, external :: true_solution
    double precision, dimension(-1:2) :: cell_values
    integer :: i, order

    if (bc_options(1, 1) == 'p') return

    call get_true_solution_by_index_range(-1, 2, t, cell_values)
    if (bc_options(1, 1) == '0') then
        lower_values(1) = (cell_values(0) + cell_values(1)) / 2.d0
    else if (bc_options(1, 1) == '1') then
        lower_values(1) = (cell_values(1) - cell_values(0)) / dx
    end if

    call get_true_solution_by_index_range(mx-1, mx+2, t, cell_values)
    if (bc_options(1, 1) == '0') then
        upper_values(1) = (cell_values(0) + cell_values(1)) / 2.d0
    else if (bc_options(1, 1) == '1') then
        upper_values(1) = (cell_values(1) - cell_values(0)) / dx
    end if

end subroutine set_implicit_boundary_data
