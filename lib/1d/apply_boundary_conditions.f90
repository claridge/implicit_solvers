subroutine apply_linearized_bcs(r, p)

! Apply linearized boundary conditions to the Newton perturbation.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    character(len=2), dimension(2, 10) :: bc_options
    common /bc_config/ bc_options

    double precision, dimension(1-mbc:mx+mbc, meqn), intent(in) :: r
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(inout) :: p
    double precision, dimension(2) :: zeros
    integer :: i


    do i = 1, meqn
        call fill_ghost_cells_homogeneous(bc_options(:, i), p(:, i))
    end do

end subroutine apply_linearized_bcs


subroutine apply_bcs(t, q)

! Fill ghost cells with appropriate boundary values.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    character(len=2), dimension(2, 10) :: bc_options
    common /bc_config/ bc_options

    double precision, intent(in) :: t
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(inout) :: q

    double precision, dimension(2, meqn) :: lower_values, upper_values
    integer :: i
    

    call set_implicit_boundary_data(t, lower_values, upper_values)

    do i = 1, meqn
        call fill_ghost_cells(bc_options(:, i), lower_values(:, i), upper_values(:, i), q(:, i))
    end do

end subroutine apply_bcs
