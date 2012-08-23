subroutine apply_homogeneous_bcs(q)

! Apply the relevant homogeneous boundary condition operator.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(inout) :: q
    double precision, dimension(4) :: zeros

    zeros = 0.d0
    call fill_2_ghost_cells('12', zeros, '12', zeros, q(:, 1))

end subroutine apply_homogeneous_bcs


subroutine apply_bcs(t, q)

! Fill ghost cells with appropriate boundary values.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision, dimension(4) :: d0_stencil, d1_stencil, d2_stencil,  &
        d3_stencil, stencils(4, 0:3)
    common /stencil_config/ d0_stencil, d1_stencil, d2_stencil, d3_stencil, stencils

    double precision, intent(in) :: t
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(inout) :: q
    double precision :: boundary_value, delta, x
    integer :: i

    double precision, dimension(-1:2) :: cellwise_lower_values, cellwise_upper_values
    double precision, dimension(2) :: lower_values, upper_values
    double precision, external :: true_solution
    integer :: info, order
    character*2 :: lower_option, upper_option

    call get_true_solution_by_index_range(-1, 2, t, cellwise_lower_values)
    call get_true_solution_by_index_range(mx-1, mx+2, t, cellwise_upper_values)

    lower_option = '12'
    upper_option = '12'

    do i = 1, 2
        read(lower_option(i:i), '(I1)') order
        lower_values(order) = dot_product(stencils(:, order), cellwise_lower_values) / dx**order
    end do
    do i = 1, 2
        read(upper_option(i:i), '(I1)') order
        upper_values(order) = dot_product(stencils(:, order), cellwise_upper_values) / dx**order
    end do

    call fill_2_ghost_cells('12', lower_values, '12', upper_values, q(:, 1))

end subroutine apply_bcs
