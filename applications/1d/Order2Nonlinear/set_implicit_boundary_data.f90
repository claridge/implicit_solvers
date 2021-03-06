subroutine set_implicit_boundary_data(t, lower_values, upper_values)

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    character(len=2), dimension(2, 10) :: bc_options
    common /bc_config/ bc_options

    double precision, intent(in) :: t
    double precision, dimension(2), intent(out) :: lower_values, upper_values

    double precision, external :: true_solution, gen_derivative
    double precision, dimension(-1:2) :: cell_values
    integer :: i, order


    if (bc_options(1, 1) == 'p') return

    call get_true_solution_by_index_range(-1, 2, t, cell_values)
    read(bc_options(1, 1), '(I1)') order
    lower_values(1) = gen_derivative(order, cell_values, dx)

    call get_true_solution_by_index_range(mx-1, mx+2, t, cell_values)
    read(bc_options(2, 1), '(I1)') order
    upper_values(1) = gen_derivative(order, cell_values, dx)

end subroutine set_implicit_boundary_data
