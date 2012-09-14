subroutine apply_pde_operator(t, q, output)

! For the PDE q_t = g(q), calculates g(q).
!
! Args:
!   t: Time at which the operator is evaluated.
!   q: PDE solution at time t.
!   output: g(q), calculated here.

    implicit none
    
    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision :: beta, kappa, delta, mu
    common /physics_config/ beta, kappa, delta, mu
    
    double precision, intent(in) :: t
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(in), target :: q
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(out) :: output

    double precision, dimension(2) :: lower_flux, upper_flux
    integer :: ix
    double precision, external :: d0x, d1x, d2x, d3x
    double precision, dimension(1-mbc:mx+mbc) :: surface_tension

    call compute_surface_tension(mx + 2*mbc, q(:, 2), surface_tension)

    do ix = 1, mx
        call compute_flux(ix, lower_flux)
        call compute_flux(ix+1, upper_flux)
        output(ix, 1:2) = -(upper_flux - lower_flux) / dx
    end do


    contains
    
    subroutine compute_flux(ix, flux)
        implicit none
        integer, intent(in) :: ix
        double precision, dimension(2), intent(out) :: flux
        double precision :: h_face, g_face, temp

        h_face = d0x(q(ix-2:ix+1, 1))
        g_face = d0x(q(ix-2:ix+1, 2))
        flux(1) = h_face**3 / 3 * (-beta * d1x(q(ix-2:ix+1, 1)) + kappa * d3x(q(ix-2:ix+1, 1)))
        flux(2) = h_face * g_face * d1x(surface_tension(ix-2:ix+1)) -  &
            delta * d1x(q(ix-2:ix+1, 2))
    end subroutine compute_flux

end subroutine apply_pde_operator
