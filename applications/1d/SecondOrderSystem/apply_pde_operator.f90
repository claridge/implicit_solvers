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
    
    double precision, intent(in) :: t
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(in) :: q
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
        double precision :: q1_face, q2_face, q2_x
        
        q1_face = (q(ix-1, 1) + q(ix, 1)) / 2.d0
        q2_face = (q(ix-1, 2) + q(ix, 2)) / 2.d0
        q2_x = (q(ix, 2) - q(ix-1, 2)) / dx

        flux(1) = -.5d0 * q1_face**2 * q2_x
!         flux(1) = -.5d0 * q(ix-1, 1)**2 * q2_x
        flux(2) = -q1_face * q2_face * q2_x
    end subroutine compute_flux

end subroutine apply_pde_operator
