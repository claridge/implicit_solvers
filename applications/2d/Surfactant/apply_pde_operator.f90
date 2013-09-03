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
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn) :: x_flux, y_flux
    integer :: ix, iy, ieqn
    double precision, dimension(2) :: lower_flux, upper_flux
    double precision :: h_face, h_x, h_laplacian_x, g_face, g_x, surface_tension_x
    double precision :: h_y, h_laplacian_y, g_y, surface_tension_y


    call get_laplacian(q(:,:,1), h_laplacian)
    call compute_surface_tension(mx+2*mbc, my+2*mbc, q(:, :, 2), surface_tension)


    !$omp parallel do private(ix, h_face, h_x, h_laplacian_x, g_face, g_x, surface_tension_x)
    do iy = 1, my
        do ix = 1, mx+1
            h_face = (q(ix-1, iy, 1) + q(ix, iy, 1)) / 2.d0
            h_x = (q(ix, iy, 1) - q(ix-1, iy, 1)) / dx
            h_laplacian_x = (h_laplacian(ix, iy) - h_laplacian(ix-1, iy)) / dx
    
            g_face = (q(ix, iy, 2) + q(ix-1, iy, 2)) / 2.d0
            g_x = (q(ix, iy, 2) - q(ix-1, iy, 2)) / dx
            surface_tension_x = (surface_tension(ix, iy) - surface_tension(ix-1, iy)) / dx

            x_flux(ix, iy, 1) = h_face**3 / 3.d0 * (-beta * h_x + kappa * h_laplacian_x)  &
              + h_face**2 / 2.d0 * surface_tension_x

            x_flux(ix, iy, 2) = h_face * g_face * surface_tension_x - delta * g_x  &
              + kappa * h_face**2 / 2.d0 * h_laplacian_x * g_face
        end do
    end do

    !$omp parallel do private(ix, h_face, h_y, h_laplacian_y, g_face, g_y, surface_tension_y)
    do iy = 1, my+1
        do ix = 1, mx
            h_face = (q(ix, iy-1, 1) + q(ix, iy, 1)) / 2.d0
            h_y = (q(ix, iy, 1) - q(ix, iy-1, 1)) / dy
            h_laplacian_y = (h_laplacian(ix, iy) - h_laplacian(ix, iy-1)) / dy
    
            g_face = (q(ix, iy, 2) + q(ix, iy-1, 2)) / 2.d0
            g_y = (q(ix, iy, 2) - q(ix, iy-1, 2)) / dy
            surface_tension_y = (surface_tension(ix, iy) - surface_tension(ix, iy-1)) / dy        

            y_flux(ix, iy, 1) = h_face**3 / 3.d0 * (-beta * h_y + kappa * h_laplacian_y)  &
              + h_face**2 / 2.d0 * surface_tension_y

            y_flux(ix, iy, 2) = h_face * g_face * surface_tension_y - delta * g_y  &
              + kappa * h_face**2 / 2.d0 * h_laplacian_y * g_face
       end do
    end do

        
    do ieqn = 1, meqn
        !$omp parallel do private(ix)
        do iy = 1, my
            do ix = 1, mx
                output(ix, iy, ieqn) = - (x_flux(ix+1, iy, ieqn) - x_flux(ix, iy, ieqn)) / dx  &
                    - (y_flux(ix, iy+1, ieqn) - y_flux(ix, iy, ieqn)) / dy
            end do
        end do
    end do
    
end subroutine apply_pde_operator
