subroutine set_implicit_boundary_data(t, x_lower_values, x_upper_values,  &
                                      y_lower_values, y_upper_values)

    implicit none
    
    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn

    double precision :: beta, kappa, delta, mu, right_film_height
    common /surfactant_config/ beta, kappa, delta, mu, right_film_height
    
    double precision, intent(in) :: t
    double precision, dimension(2, my, meqn), intent(out) :: x_lower_values, x_upper_values
    double precision, dimension(2, mx, meqn), intent(out) :: y_lower_values, y_upper_values

    integer :: iy
    
    x_lower_values(:, :, 1) = 0.d0
    x_upper_values(1, :, 1) = right_film_height
    x_upper_values(2, :, 1) = 0.d0

    x_lower_values(:, :, 2) = 0.d0
    x_upper_values(:, :, 2) = 0.d0
    
    ! y is periodic, so we don't need to set the y boundary values.
    
end subroutine set_implicit_boundary_data
