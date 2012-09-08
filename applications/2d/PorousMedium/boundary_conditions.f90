subroutine set_implicit_boundary_data(t, x_lower_values, x_upper_values,  &
                                      y_lower_values, y_upper_values)

    implicit none
    
    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn
    
    double precision, intent(in) :: t
    double precision, dimension(2, my), intent(out) :: x_lower_values, x_upper_values
    double precision, dimension(2, mx), intent(out) :: y_lower_values, y_upper_values

    character(len=2), dimension(4, 10) :: bc_options
    common /bc_config/ bc_options

    integer :: order, ix, iy
    double precision, external :: true_solution, gen_derivative
    double precision, dimension(4, my) :: cell_values_my
    double precision, dimension(mx, 4) :: cell_values_mx
    double precision, dimension(4) :: cell_values


    x_lower_values = 0.d0
    x_upper_values = 0.d0
    y_lower_values = 0.d0
    y_upper_values = 0.d0

    if (bc_options(1, 1) /= 'p') then
        read(bc_options(1, 1), '(I1)') order
        do iy = 1, my
            call get_true_solution_by_index_range(-1, 2, iy, iy, t, cell_values)
            x_lower_values(1, iy) = gen_derivative(order, cell_values, dx)
        end do

        read(bc_options(2, 1), '(I1)') order
        do iy = 1, my
            call get_true_solution_by_index_range(mx-1, mx+2, iy, iy, t, cell_values)
            x_upper_values(1, iy) = gen_derivative(order, cell_values, dx)
        end do
    end if

    if (bc_options(3, 1) /= 'p') then
        read(bc_options(3, 1), '(I1)') order
        do ix = 1, my
            call get_true_solution_by_index_range(ix, ix, -1, 2, t, cell_values)
            y_lower_values(1, ix) = gen_derivative(order, cell_values, dy)
        end do

        read(bc_options(4, 1), '(I1)') order
        do ix = 1, my
            call get_true_solution_by_index_range(ix, ix, my-1, my+2, t, cell_values)
            y_upper_values(1, ix) = gen_derivative(order, cell_values, dy)
        end do
    end if
        
end subroutine set_implicit_boundary_data
