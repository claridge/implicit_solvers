subroutine apply_linearized_pde_operator(t, r, p, output)

! For the PDE q_t = g(r), calculates g'[r](p), i.e. the linearization
! of g, about r, applied to perturbation p
!
! Args:
!   t: Time at which the operator is evaluated.
!   r: Base function of the linearization.
!   p: Perturbation to which the linearized operator is applied.
!   output: g'[r](p), calculated here.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: dx, x_lower
    common /claw_config/ mx, mbc, x_lower, dx, meqn
    
    double precision, intent(in) :: t
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(in) :: r, p
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(out) :: output

    double precision, dimension(2) :: lower_flux, upper_flux
    integer :: ix


    do ix = 1, mx
        call compute_flux(ix, lower_flux)
        call compute_flux(ix+1, upper_flux)
        output(ix, :) = -(upper_flux - lower_flux) / dx
    end do


    contains


    subroutine compute_flux(ix, flux)
        implicit none
        integer, intent(in) :: ix
        double precision, dimension(2), intent(out) :: flux
        double precision :: r1_face, r2_face, r2_x, p1_face, p2_face, p2_x
    
        r1_face = (r(ix-1, 1) + r(ix, 1)) / 2.d0
        r2_face = (r(ix-1, 2) + r(ix, 2)) / 2.d0
        r2_x = (r(ix, 2) - r(ix-1, 2)) / dx
        p1_face = (p(ix-1, 1) + p(ix, 1)) / 2.d0
        p2_face = (p(ix-1, 2) + p(ix, 2)) / 2.d0
        p2_x = (p(ix, 2) - p(ix-1, 2)) / dx

        flux(1) = -(.5d0 * r1_face**2 * p2_x + r1_face * p1_face * r2_x)
        flux(2) = -((r1_face * p2_face + p1_face * r2_face) * r2_x  &
            + r1_face * r2_face * p2_x)
    end subroutine compute_flux

    
end subroutine apply_linearized_pde_operator
