subroutine apply_homogeneous_bcs(q)

! Apply the relevant homogeneous boundary condition operator.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision, dimension(1-mbc:mx+mbc, meqn), intent(inout) :: q

    call extend_to_ghost_cells('even', 'even', q(:, 1))
    call extend_to_ghost_cells('even', 'even', q(:, 2))

end subroutine apply_homogeneous_bcs


subroutine apply_bcs(t, q)

! Fill ghost cells with appropriate boundary values.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision, intent(in) :: t
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(inout) :: q

    call extend_to_ghost_cells('even', 'even', q(:, 1))
    call extend_to_ghost_cells('even', 'even', q(:, 2))

end subroutine apply_bcs
