subroutine apply_homogeneous_bcs(q)

! Apply the relevant homogeneous boundary condition operator.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision, dimension(1-mbc:mx+mbc, meqn), intent(inout) :: q

    call extend_to_ghost_cells('even', 'even', q(:,1))

end subroutine apply_homogeneous_bcs


subroutine apply_bcs(t, q)

! Fill ghost cells with appropriate boundary values.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision :: h_left, h_right
    common /boundary_config/ h_left, h_right

    double precision, intent(in) :: t
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(inout) :: q

    double precision :: boundary_value, delta
    integer :: i

    call extend_to_ghost_cells('even', 'even', q(:, 1))

    ! ! Left boundary
    ! delta = 2.d0 * (h_left - q(1, 1))
    ! do i = 1, mbc
    !     q(1 - i, 1) = q(i, 1) + i * delta
    ! end do

    ! ! Right boundary
    ! delta = 2.d0 * (h_right - q(mx, 1))
    ! do i = 1, mbc
    !     q(mx + i, 1) = q(mx, 1) + i * delta
    ! end do

end subroutine apply_bcs
