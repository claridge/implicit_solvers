subroutine apply_homogeneous_bcs(q)

! Apply the relevant homogeneous boundary condition operator.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision, dimension(1-mbc:mx+mbc, meqn), intent(inout) :: q

    double precision, dimension(4) :: zeros

    ! zeros = 0.d0
    ! call apply_bcs_internal(zeros, zeros, q)

    call extend_to_ghost_cells('odd', 'odd', q)

end subroutine apply_homogeneous_bcs


subroutine apply_bcs(t, q)

! Fill ghost cells with appropriate boundary values.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision, intent(in) :: t
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(inout) :: q

    double precision :: boundary_value, delta, x
    integer :: i

    double precision, dimension(-1:2) :: lower_values, upper_values
    double precision, external :: true_solution
    double precision, dimension(4) :: d0_stencil, d1_stencil
    integer :: info

    ! call extend_to_ghost_cells('odd', 'odd', q)
    ! boundary_value = true_solution(x_lower, t)
    ! q(1-mbc:0, 1) = q(1-mbc:0, 1) + 2 * boundary_value

    ! boundary_value = true_solution(x_lower + mx * dx, t)
    ! q(mx+1:mx+mbc, 1) = q(mx+1:mx+mbc, 1) + 2 * boundary_value


    do i = -1, 2
        x = x_lower + (i - .5d0) * dx
        lower_values(i) = true_solution(x, t)
    end do
    do i = -1, 2
        x = x_lower + (mx + i - .5d0) * dx
        upper_values(i) = true_solution(x, t)
    end do
    ! print *, ''
    ! print *, 'Before'
    ! print *, lower_values
    ! print *, q(-1:2, 1)
 
    call apply_bcs_internal(lower_values, upper_values, q)

    d0_stencil = (/ -0.0625d0,  0.5625d0,  0.5625d0, -0.0625d0 /)
    print *, ''
    print *, true_solution(x_lower, t), (q(0,1) + q(1,1)) / 2.d0, dot_product(d0_stencil, q(-1:2, 1))
    print *, lower_values
    print *, q(-1:2, 1)


    ! d1_stencil = (/ 1.d0/24, -1.125d0, 1.125d0, -1.d0/24 /)
    ! print *, ''
    ! print *, 'After'
    ! print *, lower_values
    ! print *, q(-1:2, 1)
    ! print *, dot_product(d0_stencil, q(-1:2, 1))
    ! print *, dot_product(d0_stencil, lower_values)
    ! print *, dot_product(d1_stencil, q(-1:2, 1))
    ! print *, dot_product(d1_stencil, lower_values)

end subroutine apply_bcs


subroutine apply_bcs_internal(lower_values, upper_values, q)

! Creates and solves the linear system needed to apply boundary values, matching
! q(0, t) and q_x(0, t) to the supplied data.
!
! Args:
!   lower_values: True solution in cells -1:2.
!   upper_values: True solution in cells mx-1:mx+2.
!   q: The function being supplied with boundary data.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision, dimension(4), intent(in) :: lower_values, upper_values
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(inout) :: q

    double precision, dimension(4) :: d0_stencil, d1_stencil
    double precision, dimension(2, 2) :: a, b
    integer, dimension(2) :: ipiv
    integer :: info


    ! ! -1/16, 9/16, 9/16, -1/16
    ! d0_stencil = (/ -0.0625d0,  0.5625d0,  0.5625d0, -0.0625d0 /)
    d0_stencil = (/ 0.d0,  .5d0,  .5d0, 0.d0 /)

    ! ! 1/24, -9/8, 9/8, -1/24
    d1_stencil = (/ 1.d0/24, -1.125d0, 1.125d0, -1.d0/24 /) / dx

    a(1, 1:2) = d0_stencil(1:2)
    b(1, 1:2) = d0_stencil(3:4)
    a(2, 1:2) = d1_stencil(1:2)
    b(2, 1:2) = d1_stencil(3:4)

    q(-1:0, 1) = matmul(a, lower_values(1:2)) +  &
        matmul(b, lower_values(3:4) - q(1:2, 1))
    ! print *, ''
    ! print *, 'q lower'
    ! print *, q(-1:0, 1)
    call dgesv(2, 1, a, 2, ipiv, q(-1:0, 1), 2, info)
    ! print *, q(-1:0, 1)

    a(1, 1:2) = d0_stencil(1:2)
    b(1, 1:2) = d0_stencil(3:4)
    a(2, 1:2) = d1_stencil(1:2)
    b(2, 1:2) = d1_stencil(3:4)

    q(mx+1:mx+2, 1) = matmul(a, upper_values(1:2) - q(mx-1:mx, 1)) +  &
        matmul(b, upper_values(3:4))
    ! print *, ''
    ! print *, 'q upper'
    ! print *, q(mx+1:mx+2, 1)
    call dgesv(2, 1, b, 2, ipiv, q(mx+1:mx+2, 1), 2, info)
    ! print *, q(mx+1:mx+2, 1)

end subroutine apply_bcs_internal
