double precision function d0x(q)
    implicit none
    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn
    double precision, dimension(4, 0:3) :: stencils
    common /stencil_config/ stencils
    double precision, dimension(4), intent(in) :: q
    d0x = dot_product(stencils(:, 0), q)
end function d0x


double precision function d1x(q)
    implicit none
    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn
    double precision, dimension(4, 0:3) :: stencils
    common /stencil_config/ stencils
    double precision, dimension(4), intent(in) :: q
    d1x = dot_product(stencils(:, 1), q) / dx
end function d1x


double precision function d2x(q)
    implicit none
    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn
    double precision, dimension(4, 0:3) :: stencils
    common /stencil_config/ stencils
    double precision, dimension(4), intent(in) :: q
    d2x = dot_product(stencils(:, 2), q) / dx**2
end function d2x


double precision function d3x(q)
    implicit none
    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn
    double precision, dimension(4, 0:3) :: stencils
    common /stencil_config/ stencils
    double precision, dimension(4), intent(in) :: q
    d3x = dot_product(stencils(:, 3), q) / dx**3
end function d3x
