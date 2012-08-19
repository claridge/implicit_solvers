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

    integer :: ix, iy
    double precision, dimension(4) :: zeros

    zeros = 0.d0

    do iy = 1, my
        call apply_bcs_internal(zeros, zeros, q(:, iy, 1), mx)
    end do

    do ix = 1, mx
        call apply_bcs_internal(zeros, zeros, q(ix, :, 1), my)
    end do

    ! Fill innermost corner cells.
    q(0, 0, 1) = - q(1, 0, 1) - q(0, 1, 1) - q(1, 1, 1)
    q(mx+1, 0, 1) = - q(mx, 0, 1) - q(mx+1, 1, 1) - q(mx, 1, 1)
    q(0, my+1, 1) = - q(1, my+1, 1) - q(0, my, 1) - q(1, my, 1)
    q(mx+1, my+1, 1) = - q(mx, my+1, 1) - q(mx+1, my, 1) - q(mx, my, 1)

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
    double precision :: x, y, x_upper, y_upper
    double precision, external :: true_solution
    integer :: ix, iy, i
    double precision, dimension(-1:2) :: lower_values, upper_values

    x_upper = x_lower + mx * dx
    y_upper = y_lower + my * dy

    do iy = 1, my
        y = y_lower + (iy - .5d0) * dy
        do i = -1, 2
            x = x_lower + (i - .5d0) * dx
            lower_values(i) = true_solution(x, y, t)
        end do        
        
        do i = -1, 2
            x = x_upper + (i - .5d0) * dx
            upper_values(i) = true_solution(x, y, t)
        end do

        call apply_bcs_internal(lower_values, upper_values, q(:, iy, 1), mx)
    end do

    do ix = 1, mx
        x = x_lower + (ix - .5d0) * dx
        do i = -1, 2
            y = y_lower + (i - .5d0) * dy
            lower_values(i) = true_solution(x, y, t)
        end do        
    
        do i = -1, 2
            y = y_upper + (i - .5d0) * dy
            upper_values(i) = true_solution(x, y, t)
        end do

        call apply_bcs_internal(lower_values, upper_values, q(ix, :, 1), my)
    end do
    
    ! Fill innermost corner cells.
    q(0, 0, 1) = 4 * true_solution(x_lower, y_lower, t)  &
        - q(1, 0, 1) - q(0, 1, 1) - q(1, 1, 1)
    q(mx+1, 0, 1) = 4 * true_solution(x_upper, y_lower, t)  &
        - q(mx, 0, 1) - q(mx+1, 1, 1) - q(mx, 1, 1)
    q(0, my+1, 1) = 4 * true_solution(x_lower, y_upper, t)  &
        - q(1, my+1, 1) - q(0, my, 1) - q(1, my, 1)
    q(mx+1, my+1, 1) = 4 * true_solution(x_upper, y_upper, t)  &
        - q(mx, my+1, 1) - q(mx+1, my, 1) - q(mx, my, 1)

end subroutine apply_bcs


subroutine apply_bcs_internal(lower_values, upper_values, q, m)

! TODO: This only works for x-parallel slices right now.

! Creates and solves the linear system needed to apply boundary values, matching
! q(0, t) and q_x(0, t) to the supplied data.
!
! Args:
!   lower_values: True solution in cells -1:2.
!   upper_values: True solution in cells mx-1:mx+2.
!   q: The function being supplied with boundary data.

    implicit none

    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn

    double precision, dimension(4) :: d0_stencil, d1_stencil, d2_stencil,  &
        d3_stencil
    common /stencil_config/ d0_stencil, d1_stencil, d2_stencil, d3_stencil

    double precision, dimension(4), intent(in) :: lower_values, upper_values
    integer, intent(in) :: m
    double precision, dimension(1-mbc:m+mbc), intent(inout) :: q

    double precision, dimension(2, 2) :: a, b
    integer, dimension(2) :: ipiv
    integer :: info


    a(1, 1:2) = d0_stencil(1:2)
    b(1, 1:2) = d0_stencil(3:4)
    a(2, 1:2) = d1_stencil(1:2) / dx**2
    b(2, 1:2) = d1_stencil(3:4) / dx**2

    q(-1:0) = matmul(a, lower_values(1:2)) +  &
        matmul(b, lower_values(3:4) - q(1:2))
    call dgesv(2, 1, a, 2, ipiv, q(-1:0), 2, info)

    a(1, 1:2) = d0_stencil(1:2)
    b(1, 1:2) = d0_stencil(3:4)
    a(2, 1:2) = d1_stencil(1:2) / dx**2
    b(2, 1:2) = d1_stencil(3:4) / dx**2

    q(m+1:m+2) = matmul(a, upper_values(1:2) - q(m-1:m)) +  &
        matmul(b, upper_values(3:4))
    call dgesv(2, 1, b, 2, ipiv, q(m+1:m+2), 2, info)

end subroutine apply_bcs_internal
