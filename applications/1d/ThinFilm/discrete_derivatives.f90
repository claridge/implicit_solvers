double precision function d0x(q)
    implicit none
    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn
    double precision, dimension(4) :: d0_stencil, d1_stencil, d2_stencil,  &
        d3_stencil
    common /stencil_config/ d0_stencil, d1_stencil, d2_stencil, d3_stencil
    double precision, dimension(4), intent(in) :: q
    d0x = dot_product(d0_stencil, q)
end function d0x


double precision function d1x(q)
    implicit none
    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn
    double precision, dimension(4) :: d0_stencil, d1_stencil, d2_stencil,  &
        d3_stencil
    common /stencil_config/ d0_stencil, d1_stencil, d2_stencil, d3_stencil
    double precision, dimension(4), intent(in) :: q
    d1x = dot_product(d1_stencil, q) / dx
end function d1x


double precision function d2x(q)
    implicit none
    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn
    double precision, dimension(4) :: d0_stencil, d1_stencil, d2_stencil,  &
        d3_stencil
    common /stencil_config/ d0_stencil, d1_stencil, d2_stencil, d3_stencil
    double precision, dimension(4), intent(in) :: q
    d2x = dot_product(d2_stencil, q) / dx**2
end function d2x


double precision function d3x(q)
    implicit none
    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn
    double precision, dimension(4) :: d0_stencil, d1_stencil, d2_stencil,  &
        d3_stencil
    common /stencil_config/ d0_stencil, d1_stencil, d2_stencil, d3_stencil
    double precision, dimension(4), intent(in) :: q
    d3x = dot_product(d3_stencil, q) / dx**3
end function d3x
