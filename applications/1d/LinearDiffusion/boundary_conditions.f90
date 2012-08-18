subroutine apply_homogeneous_bcs(q)

! Apply the relevant homogeneous boundary condition operator.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision, dimension(1-mbc:mx+mbc, meqn), intent(inout) :: q

    call extend_to_ghost_cells('odd', 'odd', q(:,1))

end subroutine apply_homogeneous_bcs


subroutine apply_bcs(t, q)

! Fill ghost cells with appropriate boundary values.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision, intent(in) :: t
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(inout) :: q

    double precision :: boundary_value, delta
    integer :: i
    double precision, external :: true_solution

    call extend_to_ghost_cells('odd', 'odd', q(:, 1))

    boundary_value = true_solution(x_lower, t)
    q(1-mbc:0, 1) = q(1-mbc:0, 1) + 2 * boundary_value

    boundary_value = true_solution(x_lower + mx * dx, t)
    q(mx+1:mx+mbc, 1) = q(mx+1:mx+mbc, 1) + 2 * boundary_value
    
    print *, q(0, 1), q(mx+1, 1), true_solution(x_lower + mx*dx, t), t

end subroutine apply_bcs
