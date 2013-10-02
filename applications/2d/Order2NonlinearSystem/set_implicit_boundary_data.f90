subroutine set_implicit_boundary_data(t, x_lower_values, x_upper_values,  &
                                      y_lower_values, y_upper_values)

    implicit none
    
    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn
    
    double precision, intent(in) :: t
    double precision, dimension(2, my, meqn), intent(out) :: x_lower_values, x_upper_values
    double precision, dimension(2, mx, meqn), intent(out) :: y_lower_values, y_upper_values

    integer :: ix, iy
    double precision :: x_upper, y_upper, x, y


    x_lower_values = 0.d0
    x_upper_values = 0.d0
    y_lower_values = 0.d0
    y_upper_values = 0.d0

    x_upper = x_lower + mx * dx
    y_upper = y_lower + my * dy

    do iy = 1, my
        y = y_lower + (iy - .5d0) * dy
        call set_true_solution(x_lower, y, t, x_lower_values(1, iy, :))
        call set_true_solution(x_upper, y, t, x_upper_values(1, iy, :))
    end do
    
    do ix = 1, mx
        x = x_lower + (ix - .5d0) * dx
        call set_true_solution(x, y_lower, t, y_lower_values(1, ix, :))
        call set_true_solution(x, y_upper, t, y_upper_values(1, ix, :))
    end do
        
end subroutine set_implicit_boundary_data
