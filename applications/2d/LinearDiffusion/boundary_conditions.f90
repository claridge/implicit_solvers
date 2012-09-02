subroutine set_implicit_boundary_data(t, x_lower_values, x_upper_values,  &
                                      y_lower_values, y_upper_values)

    implicit none
    
    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn
    
    double precision, intent(in) :: t
    double precision, dimension(2, my), intent(out) :: x_lower_values, x_upper_values
    double precision, dimension(2, mx), intent(out) :: y_lower_values, y_upper_values

    character(len=2), dimension(4, 10) :: bc_options
    common /bc_config/ bc_options

    double precision :: x_upper, y_upper, x, y
    integer :: ix, iy
    double precision, external :: true_solution
    double precision, dimension(4, my) :: cell_values_my
    double precision, dimension(mx, 4) :: cell_values_mx
    double precision, dimension(4) :: d0_stencil

    x_upper = x_lower + mx * dx
    y_upper = y_lower + my * dy
    
    
    d0_stencil = (/ -0.0625d0,  0.5625d0,  0.5625d0, -0.0625d0 /)

    call get_true_solution_by_index_range(-1, 2, 1, my, t, cell_values_my)
    do iy = 1, my
        x_lower_values(1, iy) = dot_product(cell_values_my(:, iy), d0_stencil)
    end do

    call get_true_solution_by_index_range(mx-1, mx+2, 1, my, t, cell_values_my)
    do iy = 1, my
        x_upper_values(1, iy) = dot_product(cell_values_my(:, iy), d0_stencil)
    end do

    
    call get_true_solution_by_index_range(1, mx, -1, 2, t, cell_values_mx)
    do ix = 1, my
        y_lower_values(1, ix) = dot_product(cell_values_mx(ix, :), d0_stencil)
    end do

    call get_true_solution_by_index_range(1, mx, my-1, my+2, t, cell_values_mx)
    do ix = 1, my
        y_upper_values(1, ix) = dot_product(cell_values_mx(ix, :), d0_stencil)
    end do
        
end subroutine set_implicit_boundary_data



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
