subroutine apply_homogeneous_bcs(q)

! Apply the relevant homogeneous boundary operator.

    implicit none

    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn

    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn), intent(inout) :: q

    integer :: ix, iy
    double precision, dimension(4) :: zeros, x_temp, y_temp

    zeros = 0.d0

    do iy = 1, my
        call fill_lower_values_12(zeros, q(-1:2, iy, 1))
        call fill_upper_values_12(zeros, q(mx-1:mx+2, iy, 1))
    end do

    do ix = 1, mx
        call fill_lower_values_12(zeros, q(ix, -1:2, 1))
        call fill_upper_values_12(zeros, q(ix, my-1:my+2, 1))
    end do
    
    ! Corner (0,0)
    x_temp = q(-1:2, 0, 1)
    call fill_lower_values_12(zeros, x_temp)
    y_temp = q(0, -1:2, 1)
    call fill_lower_values_12(zeros, y_temp)
    q(0, 0, 1) = (x_temp(2) + y_temp(2)) / 2.d0

    ! Corner (mx+1,0)
    x_temp = q(mx-1:mx+2, 0, 1)
    call fill_upper_values_12(zeros, x_temp)
    y_temp = q(mx+1, -1:2, 1)
    call fill_lower_values_12(zeros, y_temp)
    q(mx+1, 0, 1) = (x_temp(3) + y_temp(2)) / 2.d0

    ! Corner (0,my+1)
    x_temp = q(-1:2, my+1, 1)
    call fill_lower_values_12(zeros, x_temp)
    y_temp = q(0, my-1:my+2, 1)
    call fill_upper_values_12(zeros, y_temp)
    q(0, my+1, 1) = (x_temp(2) + y_temp(3)) / 2.d0

    ! Corner (mx+1,my+1)
    x_temp = q(mx-1:mx+2, my+1, 1)
    call fill_upper_values_12(zeros, x_temp)
    y_temp = q(mx+1, my-1:my+2, 1)
    call fill_upper_values_12(zeros, y_temp)
    q(mx+1, my+1, 1) = (x_temp(3) + y_temp(3)) / 2.d0

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
    double precision, dimension(4) :: x_temp, y_temp

    x_upper = x_lower + mx * dx
    y_upper = y_lower + my * dy

    do iy = 1, my
        y = y_lower + (iy - .5d0) * dy
        do i = -1, 2
            x = x_lower + (i - .5d0) * dx
            lower_values(i) = true_solution(x, y, t)
        end do        
        call fill_lower_values_12(lower_values, q(-1:2, iy, 1))
        
        do i = -1, 2
            x = x_upper + (i - .5d0) * dx
            upper_values(i) = true_solution(x, y, t)
        end do
        call fill_upper_values_12(upper_values, q(mx-1:mx+2, iy, 1))
    end do

    do ix = 1, mx
        x = x_lower + (ix - .5d0) * dx
        do i = -1, 2
            y = y_lower + (i - .5d0) * dy
            lower_values(i) = true_solution(x, y, t)
        end do     
        call fill_lower_values_12(lower_values, q(ix, -1:2, 1))   
    
        do i = -1, 2
            y = y_upper + (i - .5d0) * dy
            upper_values(i) = true_solution(x, y, t)
        end do
        call fill_upper_values_12(upper_values, q(ix, my-1:my+2, 1))
    end do
    
    ! Corner (0,0)
    y = y_lower - dy / 2.d0
    do i = -1, 2
        x = x_lower + (i - .5d0) * dx
        lower_values(i) = true_solution(x, y, t)
    end do
    x_temp = q(-1:2, 0, 1)
    call fill_lower_values_12(lower_values, x_temp)

    x = x_lower - dx / 2.d0
    do i = -1, 2
        y = y_lower + (i - .5d0) * dy
        lower_values(i) = true_solution(x, y, t)
    end do     
    y_temp = q(0, -1:2, 1)
    call fill_lower_values_12(lower_values, y_temp)

    q(0, 0, 1) = (x_temp(2) + y_temp(2)) / 2.d0
    
    
    ! Corner (mx+1,0)
    y = y_lower - dy / 2.d0
    do i = -1, 2
        x = x_upper + (i - .5d0) * dx
        upper_values(i) = true_solution(x, y, t)
    end do
    x_temp = q(mx-1:mx+2, 0, 1)
    call fill_upper_values_12(upper_values, x_temp)

    x = x_upper + dx / 2.d0
    do i = -1, 2
        y = y_lower + (i - .5d0) * dy
        lower_values(i) = true_solution(x, y, t)
    end do     
    y_temp = q(mx+1, -1:2, 1)
    call fill_lower_values_12(lower_values, y_temp)

    q(mx+1, 0, 1) = (x_temp(3) + y_temp(2)) / 2.d0
    
    
    ! Corner (0,my+1)
    y = y_upper + dy / 2.d0
    do i = -1, 2
        x = x_lower + (i - .5d0) * dx
        lower_values(i) = true_solution(x, y, t)
    end do
    x_temp = q(-1:2, my+1, 1)
    call fill_lower_values_12(lower_values, x_temp)
    
    x = x_lower - dx / 2.d0
    do i = -1, 2
        y = y_upper + (i - .5d0) * dy
        upper_values(i) = true_solution(x, y, t)
    end do
    y_temp = q(0, my-1:my+2, 1)
    call fill_upper_values_12(upper_values, y_temp)
    
    q(0, my+1, 1) = (x_temp(2) + y_temp(3)) / 2.d0
     

    ! Corner (mx+1,my+1)
    y = y_upper + dy / 2.d0
    do i = -1, 2
        x = x_upper + (i - .5d0) * dx
        upper_values(i) = true_solution(x, y, t)
    end do
    x_temp = q(mx-1:mx+2, my+1, 1)
    call fill_upper_values_12(upper_values, x_temp)
    
    x = x_upper + dx / 2.d0
    do i = -1, 2
        y = y_upper + (i - .5d0) * dy
        upper_values(i) = true_solution(x, y, t)
    end do
    y_temp = q(mx+1, my-1:my+2, 1)
    call fill_upper_values_12(upper_values, y_temp)
    
    q(mx+1, my+1, 1) = (x_temp(3) + y_temp(3)) / 2.d0
    
!     print *, q(0, 0, 1), true_solution(x_lower-dx/2, y_lower-dy/2, t)
    
end subroutine apply_bcs


subroutine fill_lower_values_01(true_values, q_values)

    implicit none
    
    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn
    
    double precision, dimension(4) :: d0_stencil, d1_stencil, d2_stencil,  &
        d3_stencil
    common /stencil_config/ d0_stencil, d1_stencil, d2_stencil, d3_stencil

    double precision, dimension(4), intent(in) :: true_values
    double precision, dimension(4), intent(inout) :: q_values
    
    double precision, dimension(2, 2) :: a, b
    integer, dimension(2) :: ipiv
    integer :: info
    
    a(1, 1:2) = d0_stencil(1:2)
    b(1, 1:2) = d0_stencil(3:4)
    a(2, 1:2) = d1_stencil(1:2) / dx**2
    b(2, 1:2) = d1_stencil(3:4) / dx**2

    q_values(1:2) = matmul(a, true_values(1:2)) + matmul(b, true_values(3:4) - q_values(3:4))
    call dgesv(2, 1, a, 2, ipiv, q_values(1:2), 2, info)
        
end subroutine fill_lower_values_01


subroutine fill_lower_values_12(true_values, q_values)

    implicit none
    
    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn
    
    double precision, dimension(4) :: d0_stencil, d1_stencil, d2_stencil,  &
        d3_stencil
    common /stencil_config/ d0_stencil, d1_stencil, d2_stencil, d3_stencil

    double precision, dimension(4), intent(in) :: true_values
    double precision, dimension(4), intent(inout) :: q_values
    
    double precision, dimension(2, 2) :: a, b
    integer, dimension(2) :: ipiv
    integer :: info
    
    a(1, 1:2) = d0_stencil(1:2)
    b(1, 1:2) = d0_stencil(3:4)
    a(2, 1:2) = d3_stencil(1:2)
    b(2, 1:2) = d3_stencil(3:4)

    q_values(1:2) = matmul(a, true_values(1:2)) + matmul(b, true_values(3:4) - q_values(3:4))
    call dgesv(2, 1, a, 2, ipiv, q_values(1:2), 2, info)
        
end subroutine fill_lower_values_12


subroutine fill_upper_values_01(true_values, q_values)

    implicit none

    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn

    double precision, dimension(4) :: d0_stencil, d1_stencil, d2_stencil,  &
        d3_stencil
    common /stencil_config/ d0_stencil, d1_stencil, d2_stencil, d3_stencil

    double precision, dimension(4), intent(in) :: true_values
    double precision, dimension(4), intent(inout) :: q_values

    double precision, dimension(2, 2) :: a, b
    integer, dimension(2) :: ipiv
    integer :: info
    
    a(1, 1:2) = d0_stencil(1:2)
    b(1, 1:2) = d0_stencil(3:4)
    a(2, 1:2) = d1_stencil(1:2) / dx**2
    b(2, 1:2) = d1_stencil(3:4) / dx**2

    q_values(3:4) = matmul(a, true_values(1:2) - q_values(1:2)) + matmul(b, true_values(3:4))
    call dgesv(2, 1, b, 2, ipiv, q_values(3:4), 2, info)
    
end subroutine fill_upper_values_01


subroutine fill_upper_values_12(true_values, q_values)

    implicit none

    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn

    double precision, dimension(4) :: d0_stencil, d1_stencil, d2_stencil,  &
        d3_stencil
    common /stencil_config/ d0_stencil, d1_stencil, d2_stencil, d3_stencil

    double precision, dimension(4), intent(in) :: true_values
    double precision, dimension(4), intent(inout) :: q_values

    double precision, dimension(2, 2) :: a, b
    integer, dimension(2) :: ipiv
    integer :: info
    
    a(1, 1:2) = d0_stencil(1:2)
    b(1, 1:2) = d0_stencil(3:4)
    a(2, 1:2) = d3_stencil(1:2)
    b(2, 1:2) = d3_stencil(3:4)

    q_values(3:4) = matmul(a, true_values(1:2) - q_values(1:2)) + matmul(b, true_values(3:4))
    call dgesv(2, 1, b, 2, ipiv, q_values(3:4), 2, info)
    
end subroutine fill_upper_values_12
