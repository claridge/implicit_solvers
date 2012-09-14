subroutine b4step2(maxmx, maxmy, mbc, mx, my, meqn, q,  &
                   x_low, y_low, dx, dy, t, dt, maux, aux)

!-----------------------------------------------------------------
! Called from claw2 before each call to step2.  Used to set the 
! aux arrays with velocity field data.
!
!    aux(:,:,1) = film velocity field, x-component
!    aux(:,:,2) = film velocity field, y-component
!    aux(:,:,3) = surfactant velocity field, x-component
!    aux(:,:,4) = surfactant velocity field, y-component
!-----------------------------------------------------------------

    !$ use omp_lib
    implicit none

    double precision :: beta, kappa, delta, mu, right_film_height
    common /surfactant_config/ beta, kappa, delta, mu, right_film_height

    integer, intent(in) :: maxmx, maxmy, mbc, mx, my, meqn, maux
    double precision, intent(in) :: x_low, y_low, dx, dy, t, dt
    double precision, dimension(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, meqn),  &
        intent(in) :: q

    double precision, dimension(1-mbc:maxmx+mbc,1-mbc:maxmy+mbc, maux),  &
        intent(out) :: aux

    integer :: ix, iy

    ! TODO: Could potentially use the same array here, at the expense of meaningful names.
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc) :: surface_tension, film_laplacian
    double precision, dimension(2, my) :: x_lower_values, x_upper_values
    double precision, dimension(2, mx) :: junk  ! For y's periodic boundary conditions
    character(len=2), dimension(4) :: bc_options
    double precision, external :: derivative0, derivative1

    double precision :: film_face, film_x, film_y ,film_laplacian_x, film_laplacian_y

    ! We'll fill ghost cells for the velocity fields by using direct
    ! function values in the x-direction and periodic conditions in y.
    bc_options = (/ '0', '0', 'p', 'p' /)

    
    call compute_surface_tension(mx+2*mbc, my+2*mbc, q(:, :, 2), surface_tension)

    do iy = 1, my
        x_lower_values(1, iy) = 0.d0 !derivative1(surface_tension(-1:2, iy), dx)
        do ix = 1, mx
            aux(ix, iy, 1) = centered_x_derivative(surface_tension, ix, iy)
        end do
        x_upper_values(1, iy) = 0.d0 !derivative1(surface_tension(mx-1:mx+2, iy), dx)
    end do
    call fill_ghost_cells(bc_options, x_lower_values, x_upper_values, junk, junk, aux(:, :, 1))

    do iy = 1, my
        x_lower_values(1, iy) = 0.d0 !y_derivative_at_x_interface(surface_tension, 0, iy)
        do ix = 1, mx
            aux(ix, iy, 2) = centered_y_derivative(surface_tension, ix, iy)
        end do
        x_upper_values(1, iy) = 0.d0 !y_derivative_at_x_interface(surface_tension, mx+1, iy)
    end do
    call fill_ghost_cells(bc_options, x_lower_values, x_upper_values, junk, junk, aux(:, :, 2))
    

    ! print *, x_lower_values(1, 10)
    ! print *, aux(-1:3, 10, 1)


    ! TODO: Move calculate_laplacian to library.
    call calculate_laplacian(q(:, :, 1), film_laplacian)

    ! TODO: Perhaps name derivative routines more descriptively?  interface_d1_2cell,
    ! interface_d1_4cell, etc.
    do iy = 1, my
        film_face = derivative0(q(-1:2, iy, 1))
        film_x = derivative1(q(-1:2, iy, 1), dx)
        film_laplacian_x = (film_laplacian(1, iy) - film_laplacian(0, iy)) / dx
        x_lower_values(1, iy) = 0.d0 !surfactant_velocity(film_face, film_x, film_laplacian_x)

        do ix = 1, mx
            film_x = centered_x_derivative(q(:, :, 1), ix, iy)
            film_laplacian_x = centered_x_derivative(film_laplacian, ix, iy)
            aux(ix, iy, 3) = surfactant_velocity(q(ix, iy, 1), film_x, film_laplacian_x)
        end do
        
        film_face = derivative0(q(mx-1:mx+2, iy, 1))
        film_x = derivative1(q(mx-1:mx+2, iy, 1), dx)
        film_laplacian_x = (film_laplacian(mx+1, iy) - film_laplacian(mx, iy)) / dx
        x_upper_values(1, iy) = 0.d0 !surfactant_velocity(film_face, film_x, film_laplacian_x)
    end do
    call fill_ghost_cells(bc_options, x_lower_values, x_upper_values, junk, junk, aux(:, :, 3))

    do iy = 1, my
        film_face = derivative0(q(-1:2, iy, 1))
        film_y = y_derivative_at_x_interface(q(:, :, 1), 0, iy)
        film_laplacian_y = y_derivative_at_x_interface(film_laplacian, 0, iy)
        x_lower_values(1, iy) = 0.d0 !surfactant_velocity(film_face, film_y, film_laplacian_y)

        do ix = 1, mx
            film_y = centered_y_derivative(q(:, :, 1), ix, iy)
            film_laplacian_y = centered_y_derivative(film_laplacian, ix, iy)
            aux(ix, iy, 4) = surfactant_velocity(q(ix, iy, 1), film_y, film_laplacian_y)
        end do

        film_face = derivative0(q(mx-1:mx+2, iy, 1))
        film_y = y_derivative_at_x_interface(q(:, : ,1), mx+1, iy)
        film_laplacian_y = y_derivative_at_x_interface(film_laplacian, mx+1, iy)
        x_upper_values(1, iy) = 0.d0 !surfactant_velocity(film_face, film_y, film_laplacian_y)
    end do
    call fill_ghost_cells(bc_options, x_lower_values, x_upper_values, junk, junk, aux(:, :, 4))


    contains


    ! TODO: Will these work as statement functions?
    double precision function y_derivative_at_x_interface(u, ix, iy)
        implicit none
        double precision, dimension(1-mbc:mx+mbc, 1-mbc:mx+mbc) :: u
        integer :: ix, iy
        y_derivative_at_x_interface =  &
            (u(ix-1, iy+1) + u(ix, iy+1) - u(ix-1, iy-1) - u(ix, iy-1)) / (4.d0 * dy)
    end function y_derivative_at_x_interface

    double precision function centered_x_derivative(u, ix, iy)
        implicit none
        double precision, dimension(1-mbc:mx+mbc, 1-mbc:mx+mbc) :: u
        integer :: ix, iy
        centered_x_derivative = (u(ix+1, iy) - u(ix-1, iy)) / (2.d0 * dx)
    end function centered_x_derivative

    double precision function centered_y_derivative(u, ix, iy)
        implicit none
        double precision, dimension(1-mbc:mx+mbc, 1-mbc:mx+mbc) :: u
        integer :: ix, iy
        centered_y_derivative = (u(ix, iy+1) - u(ix, iy-1)) / (2.d0 * dy)
    end function centered_y_derivative

    double precision function surfactant_velocity(film_face, d_film, d_film_laplacian)
      double precision :: film_face, d_film, d_film_laplacian
      surfactant_velocity = .5d0 * film_face**2 * (-beta * d_film + kappa * d_film_laplacian)
    end function surfactant_velocity

end subroutine b4step2
