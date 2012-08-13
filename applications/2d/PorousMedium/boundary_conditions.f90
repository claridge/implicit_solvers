subroutine bc2(maxmx, maxmy, meqn, mbc, mx, my, x_low, y_low, dx, dy,  &
               q, maux, aux, t, dt, mthbc)

! Applies boundary conditions before each call to rpn2.
!
! This version doesn't do anything.

    implicit none

    integer, intent(in) :: maxmx, maxmy, meqn, mbc, mx, my, maux
    double precision, intent(in) :: x_low, y_low, dx, dy, t, dt
    double precision :: aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, *)
    integer, dimension(4), intent(in) :: mthbc(4)
    double precision, dimension(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn),  &
        intent(inout) :: q

end subroutine bc2


subroutine apply_homogeneous_bcs(q)

! Apply the relevant homogeneous boundary operator.

    implicit none

    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn

    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn), intent(inout) :: q

    call extend_to_ghost_cells('odd', 'odd', 'odd', 'odd', q(:, :, 1))

end subroutine apply_homogeneous_bcs


subroutine apply_bcs(t, q)

! Fill ghost cells with appropriate boundary values.
!
! Corner ghost cells are not filled appropriately.

    implicit none

    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn

    double precision, intent(in) :: t
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn), intent(inout) :: q
    double precision :: boundary_value, x, y
    double precision, external :: true_solution
    integer :: ix, iy


    call extend_to_ghost_cells('odd', 'odd', 'odd', 'odd', q(:, :, 1))

    do iy = 1, my
        y = y_lower + (iy - .5d0) * dy
        boundary_value = true_solution(x_lower, y, t)
        q(1-mbc:0, iy, 1) = q(1-mbc:0, iy, 1) + 2 * boundary_value
        boundary_value = true_solution(x_lower + mx * dx, y, t)
        q(mx+1:mx+mbc, iy, 1) = q(mx+1:mx+mbc, iy, 1) + 2 * boundary_value
    end do

    do ix = 1, my
        x = x_lower + (ix - .5d0) * dx
        boundary_value = true_solution(x, y_lower, t)
        q(ix, 1-mbc:0, 1) = q(ix, 1-mbc:0, 1) + 2 * boundary_value
        boundary_value = true_solution(x, y_lower + my * dy, t)
        q(ix, my+1:my+mbc, 1) = q(ix, my+1:my+mbc, 1) + 2 * boundary_value
    end do

end subroutine apply_bcs
