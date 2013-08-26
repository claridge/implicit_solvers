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

    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc) :: h_laplacian, surface_tension
    integer :: ix, iy
    double precision, external :: derivative0, derivative1
    double precision, dimension(2) :: lower_flux, upper_flux


    call get_laplacian(q(:,:,1), h_laplacian)
    call compute_surface_tension(mx+2*mbc, my+2*mbc, q(:, :, 2), surface_tension)

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
        ! Use 2-cell stencils for the Laplacian and for g, as they don't have
        ! enough ghost cells for 4-cell stencils.
        implicit none
        integer, intent(in) :: ix, iy
        double precision, dimension(2), intent(out) :: flux
        double precision :: h_face, h_x, h_laplacian_x, g_face, g_x, surface_tension_x

        h_face = derivative0(q(ix-2:ix+1, iy, 1))
        h_x = derivative1(q(ix-2:ix+1, iy, 1), dx)
        h_laplacian_x = (h_laplacian(ix, iy) - h_laplacian(ix-1, iy)) / dx
        flux(1) = h_face**3 / 3.d0 * (-beta * h_x + kappa * h_laplacian_x)

        g_face = (q(ix, iy, 2) + q(ix-1, iy, 2)) / 2.d0
        g_x = (q(ix, iy, 2) - q(ix-1, iy, 2)) / dx
        surface_tension_x = (surface_tension(ix, iy) - surface_tension(ix-1, iy)) / dx
        flux(2) = h_face * g_face * surface_tension_x - delta * g_x
    end subroutine compute_x_flux
    
    
    subroutine compute_y_flux(ix, iy, flux)
        implicit none
        integer, intent(in) :: ix, iy
        double precision, dimension(2), intent(out) :: flux
        double precision :: h_face, h_y, h_laplacian_y, g_face, g_y, surface_tension_y
    
        h_face = derivative0(q(ix, iy-2:iy+1, 1))
        h_y = derivative1(q(ix, iy-2:iy+1, 1), dy)
        h_laplacian_y = (h_laplacian(ix, iy) - h_laplacian(ix, iy-1)) / dy
        flux(1) = h_face**3 / 3.d0 * (-beta * h_y + kappa * h_laplacian_y)

        g_face = (q(ix, iy, 2) + q(ix, iy-1, 2)) / 2.d0
        g_y = (q(ix, iy, 2) - q(ix, iy-1, 2)) / dy
        surface_tension_y = (surface_tension(ix, iy) - surface_tension(ix, iy-1)) / dy        
        flux(2) = h_face * g_face * surface_tension_y - delta * g_y
    end subroutine compute_y_flux
    
end subroutine apply_pde_operator


subroutine apply_linearized_pde_operator(t, q, p, output)

! For the PDE q_t = g(q), calculates g'[q](p), i.e. the linearization
! of g, about q, applied to perturbation p
!
! Args:
!   t: Time at which the operator is evaluated.
!   q: Base function of the linearization.
!   p: Perturbation to which the linearized operator is applied.
!   output: g'[q](p), calculated here.

    !#$ use omp_lib
    implicit none

    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn

    double precision :: beta, kappa, delta, mu, right_film_height
    common /surfactant_config/ beta, kappa, delta, mu, right_film_height
    
    double precision, intent(in) :: t
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn), intent(in) :: q, p
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn), intent(out) :: output
    
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc) :: h_laplacian, ph_laplacian,  &
        surface_tension, surface_tension_d1
    integer :: ix, iy
    double precision, external :: derivative0, derivative1
    double precision, dimension(2) :: lower_flux, upper_flux


    call get_laplacian(q(:,:,1), h_laplacian)
    call get_laplacian(p(:,:,1), ph_laplacian)
    call compute_surface_tension(mx+2*mbc, my+2*mbc, q(:,:,2), surface_tension)
    call compute_surface_tension_d1(mx+2*mbc, my+2*mbc, q(:,:,2), surface_tension_d1)
    
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
        ! Use 2-cell stencils for the Laplacian and for g, as they don't have
        ! enough ghost cells for 4-cell stencils.
        implicit none
        integer, intent(in) :: ix, iy
        double precision, dimension(2), intent(out) :: flux
        double precision :: h_face, h_x, h_laplacian_x, ph_face, ph_x, ph_laplacian_x
        double precision :: g_face, g_x, pg_face, pg_x, surface_tension_x

        h_face = derivative0(q(ix-2:ix+1, iy, 1))
        h_x = derivative1(q(ix-2:ix+1, iy, 1), dx)
        h_laplacian_x = (h_laplacian(ix, iy) - h_laplacian(ix-1, iy)) / dx
        ph_face = derivative0(p(ix-2:ix+1, iy, 1))
        ph_x = derivative1(p(ix-2:ix+1, iy, 1), dx)
        ph_laplacian_x = (ph_laplacian(ix, iy) - ph_laplacian(ix-1, iy)) / dx
        flux(1) = h_face**3 / 3.d0 * (-beta * ph_x + kappa * ph_laplacian_x)  &
            + h_face**2 * (-beta * h_x + kappa * h_laplacian_x) * ph_face

        g_face = (q(ix, iy, 2) + q(ix-1, iy, 2)) / 2.d0
        g_x = (q(ix, iy, 2) - q(ix-1, iy, 2)) / dx
        pg_face = (p(ix, iy, 2) + p(ix-1, iy, 2)) / 2.d0
        pg_x = (p(ix, iy, 2) - p(ix-1, iy, 2)) / dx
        surface_tension_x = (surface_tension(ix, iy) - surface_tension(ix-1, iy)) / dx
        flux(2) = (ph_face * g_face + h_face * pg_face) * surface_tension_x  &
            + h_face * g_face * (p(ix, iy, 2) * surface_tension_d1(ix, iy)  &
                - p(ix-1, iy, 2) * surface_tension_d1(ix-1, iy)) / dx  &
            - delta * pg_x
    end subroutine compute_x_flux
    
    
    subroutine compute_y_flux(ix, iy, flux)
        implicit none
        integer, intent(in) :: ix, iy
        double precision, dimension(2), intent(inout) :: flux
        double precision :: h_face, h_y, h_laplacian_y, ph_face, ph_y, ph_laplacian_y
        double precision :: g_face, g_y, pg_face, pg_y, surface_tension_y
        
        h_face = derivative0(q(ix, iy-2:iy+1, 1))
        h_y = derivative1(q(ix, iy-2:iy+1, 1), dy)
        h_laplacian_y = (h_laplacian(ix, iy) - h_laplacian(ix-1, iy))
        ph_y = derivative1(p(ix, iy-2:iy+1, 1), dy)
        ph_face = derivative0(p(ix, iy-2:ix+1, 1))
        ph_laplacian_y = (ph_laplacian(ix, iy) - ph_laplacian(ix, iy-1)) / dy
        flux(1) = h_face**3 / 3.d0 * (-beta * ph_y + kappa * ph_laplacian_y)  &
            + h_face**2 * (-beta * h_y + kappa * h_laplacian_y) * ph_face
            
        g_face = (q(ix, iy, 2) + q(ix, iy-1, 2)) / 2.d0
        g_y = (q(ix, iy, 2) - q(ix, iy-1, 2)) / dy
        pg_face = (p(ix, iy, 2) + p(ix, iy-1, 2)) / 2.d0
        pg_y = (p(ix, iy, 2) - p(ix, iy-1, 2)) / dy
        surface_tension_y = (surface_tension(ix, iy) - surface_tension(ix, iy-1)) / dy
        flux(2) = (ph_face * g_face + h_face * pg_face) * surface_tension_y  &
            + h_face * g_face * (p(ix, iy, 2) * surface_tension_d1(ix, iy)  &
                - p(ix, iy-1, 2) * surface_tension_d1(ix, iy-1)) / dy  &
            - delta * pg_y
    end subroutine compute_y_flux
    
end subroutine apply_linearized_pde_operator
