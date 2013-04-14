! subroutine apply_pde_operator(t, q, output)
! 
! ! For the PDE q_t = g(q), calculates g(q).
! !
! ! Args:
! !   t: Time at which the operator is evaluated.
! !   q: PDE solution at time t.
! !   output: g(q), calculated here.
! 
!     !$ use omp_lib
!     implicit none
!     
!     integer :: mx, my, mbc, meqn
!     double precision :: x_lower, y_lower, dx, dy
!     common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn
! 
!     double precision :: beta, kappa, delta, mu, right_film_height
!     common /surfactant_config/ beta, kappa, delta, mu, right_film_height
!     
!     double precision, intent(in) :: t
!     double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn), intent(in) :: q
!     double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn), intent(out) :: output
!     double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn) :: test_output
! 
!     double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc) :: h_laplacian, surface_tension
!     integer :: ix, iy, ieqn
!     double precision, external :: derivative0, derivative1
!     double precision, dimension(2) :: lower_flux, upper_flux
!     double precision :: max_diff
! 
! 
!     call calculate_laplacian(q(:,:,1), h_laplacian)
!     call compute_surface_tension(mx+2*mbc, my+2*mbc, q(:, :, 2), surface_tension)
! 
!     !$omp parallel do private(ix, lower_flux, upper_flux)
!     do iy = 1, my
!         do ix = 1, mx
!             call compute_x_flux(ix, iy, lower_flux)
!             call compute_x_flux(ix+1, iy, upper_flux)
!             output(ix, iy, :) = -(upper_flux - lower_flux) / dx
!             
!             call compute_y_flux(ix, iy, lower_flux)
!             call compute_y_flux(ix, iy+1, upper_flux)
!             output(ix, iy, :) = output(ix, iy, :) - (upper_flux - lower_flux) / dy
!         end do
!     end do
!     
! !     call test_apply_pde_operator(t, q, test_output)
! !     max_diff = 0.d0
! !     do ieqn = 1, 2
! !         do iy = 1, my
! !             do ix = 1, mx
! !                 max_diff = max(max_diff, abs(output(ix, iy, ieqn) - test_output(ix, iy, ieqn)))
! !             end do
! !         end do
! !     end do
! !     print *, 'max diff: ', max_diff
!     
!     contains
! 
! 
!     subroutine compute_x_flux(ix, iy, flux)
!         ! Use 2-cell stencils for the Laplacian and for g, as they don't have
!         ! enough ghost cells for 4-cell stencils.
!         implicit none
!         integer, intent(in) :: ix, iy
!         double precision, dimension(2), intent(out) :: flux
!         double precision :: h_face, h_x, h_laplacian_x, g_face, g_x, surface_tension_x
! 
!         h_face = derivative0(q(ix-2:ix+1, iy, 1))
!         h_x = derivative1(q(ix-2:ix+1, iy, 1), dx)
!         h_laplacian_x = (h_laplacian(ix, iy) - h_laplacian(ix-1, iy)) / dx
!         flux(1) = h_face**3 / 3.d0 * (-beta * h_x + kappa * h_laplacian_x)
! 
!         g_face = (q(ix, iy, 2) + q(ix-1, iy, 2)) / 2.d0
!         g_x = (q(ix, iy, 2) - q(ix-1, iy, 2)) / dx
!         surface_tension_x = (surface_tension(ix, iy) - surface_tension(ix-1, iy)) / dx
!         flux(2) = h_face * g_face * surface_tension_x - delta * g_x
!     end subroutine compute_x_flux
!     
!     
!     subroutine compute_y_flux(ix, iy, flux)
!         implicit none
!         integer, intent(in) :: ix, iy
!         double precision, dimension(2), intent(out) :: flux
!         double precision :: h_face, h_y, h_laplacian_y, g_face, g_y, surface_tension_y
!     
!         h_face = derivative0(q(ix, iy-2:iy+1, 1))
!         h_y = derivative1(q(ix, iy-2:iy+1, 1), dy)
!         h_laplacian_y = (h_laplacian(ix, iy) - h_laplacian(ix, iy-1)) / dy
!         flux(1) = h_face**3 / 3.d0 * (-beta * h_y + kappa * h_laplacian_y)
! 
!         g_face = (q(ix, iy, 2) + q(ix, iy-1, 2)) / 2.d0
!         g_y = (q(ix, iy, 2) - q(ix, iy-1, 2)) / dy
!         surface_tension_y = (surface_tension(ix, iy) - surface_tension(ix, iy-1)) / dy        
!         flux(2) = h_face * g_face * surface_tension_y - delta * g_y
!     end subroutine compute_y_flux
!     
! end subroutine apply_pde_operator


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
    double precision, external :: derivative0, derivative1
    double precision, dimension(2) :: lower_flux, upper_flux
    double precision :: h_face, h_x, h_laplacian_x, g_face, g_x, surface_tension_x
    double precision :: h_y, h_laplacian_y, g_y, surface_tension_y


    call calculate_laplacian(q(:,:,1), h_laplacian)
    call compute_surface_tension(mx+2*mbc, my+2*mbc, q(:, :, 2), surface_tension)


    !$omp parallel do private(ix, h_face, h_x, h_laplacian_x, g_face, g_x, surface_tension_x)
    do iy = 1, my
        do ix = 1, mx+1
            h_face = (q(ix-1, iy, 1) + q(ix, iy, 1)) / 2.d0  ! derivative0(q(ix-2:ix+1, iy, 1))
            h_x = (q(ix, iy, 1) - q(ix-1, iy, 1)) / dx  ! derivative1(q(ix-2:ix+1, iy, 1), dx)
            h_laplacian_x = (h_laplacian(ix, iy) - h_laplacian(ix-1, iy)) / dx
    
            g_face = (q(ix, iy, 2) + q(ix-1, iy, 2)) / 2.d0
            g_x = (q(ix, iy, 2) - q(ix-1, iy, 2)) / dx
            surface_tension_x = (surface_tension(ix, iy) - surface_tension(ix-1, iy)) / dx

            x_flux(ix, iy, 1) = h_face**3 / 3.d0 * (-beta * h_x + kappa * h_laplacian_x)
            x_flux(ix, iy, 2) = h_face * g_face * surface_tension_x - delta * g_x
        end do
    end do

    !$omp parallel do private(ix, h_face, h_y, h_laplacian_y, g_face, g_y, surface_tension_y)
    do iy = 1, my+1
        do ix = 1, mx
            h_face = (q(ix, iy-1, 1) + q(ix, iy, 1)) / 2.d0  ! derivative0(q(ix, iy-2:iy+1, 1))
            h_y = (q(ix, iy, 1) - q(ix, iy-1, 1)) / dy  ! derivative1(q(ix, iy-2:iy+1, 1), dy)
            h_laplacian_y = (h_laplacian(ix, iy) - h_laplacian(ix, iy-1)) / dy
    
            g_face = (q(ix, iy, 2) + q(ix, iy-1, 2)) / 2.d0
            g_y = (q(ix, iy, 2) - q(ix, iy-1, 2)) / dy
            surface_tension_y = (surface_tension(ix, iy) - surface_tension(ix, iy-1)) / dy        

            y_flux(ix, iy, 1) = h_face**3 / 3.d0 * (-beta * h_y + kappa * h_laplacian_y)
            y_flux(ix, iy, 2) = h_face * g_face * surface_tension_y - delta * g_y
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


subroutine calculate_laplacian(q, q_laplacian)

! Calculate Laplacian of q at cell centers.  Fills ghost cells within the
! outermost layer.

    !# use omp_lib
    implicit none

    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn

    double precision, dimension(1-mbc:mx+mbc,1-mbc:my+mbc), intent(in) :: q
    double precision, dimension(1-mbc:mx+mbc,1-mbc:my+mbc), intent(out) :: q_laplacian

    integer :: ix, iy

    !#omp parallel do private(ix)
    do iy = 2-mbc,my+mbc-1
        do ix = 2-mbc,mx+mbc-1
            q_laplacian(ix,iy) = (q(ix-1,iy) - 2.d0*q(ix,iy) + q(ix+1,iy)) / dx**2
            q_laplacian(ix,iy) = q_laplacian(ix,iy) + (q(ix,iy-1) - 2.d0*q(ix,iy) +  &
                                 q(ix,iy+1)) / dy**2
        end do
    end do
    
end subroutine calculate_laplacian
