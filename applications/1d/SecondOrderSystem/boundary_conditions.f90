subroutine set_implicit_boundary_data(t, lower_values, upper_values)

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    character(len=2), dimension(2, 10) :: bc_options
    common /bc_config/ bc_options

    double precision, intent(in) :: t
    double precision, dimension(2, meqn), intent(out) :: lower_values, upper_values

    call set_true_solution(x_lower, t, lower_values(1, :))
    call set_true_solution(x_lower + mx * dx, t, upper_values(1, :))

end subroutine set_implicit_boundary_data
