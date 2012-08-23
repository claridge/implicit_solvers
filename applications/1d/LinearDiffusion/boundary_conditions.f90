subroutine apply_homogeneous_bcs(q)

! Apply the relevant homogeneous boundary condition operator.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision, dimension(1-mbc:mx+mbc, meqn), intent(inout) :: q

    ! Dirichlet boundary conditions
    call fill_1_ghost_cell('0', 0.d0, '0', 0.d0, q)

    ! ! Neumann boundary conditions
    ! call fill_1_ghost_cell('1', 0.d0, '1', 0.d0, q)

end subroutine apply_homogeneous_bcs


subroutine apply_bcs(t, q)

! Fill ghost cells with appropriate boundary values.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision, intent(in) :: t
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(inout) :: q

    double precision :: x_upper, lower_value, upper_value
    integer :: i
    double precision, external :: true_solution

    ! Dirichlet boundary conditions
    call fill_1_ghost_cell('0', true_solution(x_lower, t),  &
                           '0', true_solution(x_lower + mx * dx, t),  &
                           q)

    ! ! Neumann boundary conditions
    ! lower_value = (true_solution(x_lower + dx/2, t) -  &
    !     true_solution(x_lower - dx/2, t)) / dx
    ! x_upper = x_lower + mx * dx
    ! upper_value = (true_solution(x_upper + dx/2, t) -  &
    !     true_solution(x_upper - dx/2, t)) / dx
    ! call fill_1_ghost_cell('1', lower_value, '1', upper_value, q)

end subroutine apply_bcs
