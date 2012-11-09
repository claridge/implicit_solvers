subroutine apply_pde_operator(t, q, output)

! For the PDE q_t = g(q), calculates g(q).
!
! Args:
!   t: Time at which the operator is evaluated.
!   q: PDE solution at time t.
!   output: g(q), calculated here.

    !$ use omp_lib
    implicit none
    
    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn

    double precision :: beta, kappa, delta, mu, right_film_height
    common /surfactant_config/ beta, kappa, delta, mu, right_film_height
    
    double precision, intent(in) :: t
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn), intent(in) :: q
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn), intent(out) :: output

    integer :: ix, iy
    double precision, dimension(2) :: lower_flux, upper_flux



    !$omp parallel do private(ix, lower_flux, upper_flux)
    do iy = 1, my
        do ix = 1, mx
            call compute_x_flux(ix, iy, lower_flux)
            call compute_x_flux(ix+1, iy, upper_flux)
            output(ix, iy, :) = -(upper_flux - lower_flux) / dx
            
            call compute_y_flux(ix, iy, lower_flux)
            call compute_y_flux(ix, iy+1, upper_flux)
            output(ix, iy, :) = output(ix, iy, :) - (upper_flux - lower_flux) / dy
        end do
    end do
    
    contains


    subroutine compute_x_flux(ix, iy, flux)
        implicit none
        integer, intent(in) :: ix, iy
        double precision, dimension(2), intent(out) :: flux
        double precision :: q1_face, q2_face, q2_x
        
        q1_face = (q(ix-1, iy, 1) + q(ix, iy, 1)) / 2.d0
        q2_face = (q(ix-1, iy, 2) + q(ix, iy, 2)) / 2.d0
        q2_x = (q(ix, iy, 2) - q(ix-1, iy, 2)) / dx
        
        flux(1) = -.5d0 * q1_face**2 * q2_x
        flux(2) = -q1_face * q2_face * q2_x
    end subroutine compute_x_flux
    
    
    subroutine compute_y_flux(ix, iy, flux)
        implicit none
        integer, intent(in) :: ix, iy
        double precision, dimension(2), intent(out) :: flux
        double precision :: q1_face, q2_face, q2_y

        q1_face = (q(ix, iy-1, 1) + q(ix, iy, 1)) / 2.d0
        q2_face = (q(ix, iy-1, 2) + q(ix, iy, 2)) / 2.d0
        q2_y = (q(ix, iy, 2) - q(ix, iy-1, 2)) / dy
        
        flux(1) = -.5d0 * q1_face**2 * q2_y
        flux(2) = -q1_face * q2_face * q2_y
    end subroutine compute_y_flux
    
end subroutine apply_pde_operator