subroutine compute_surface_tension(nx, ny, surfactant, surface_tension)
    implicit none
    double precision :: beta, kappa, delta, mu, right_film_height
    common /surfactant_config/ beta, kappa, delta, mu, right_film_height

    integer, intent(in) :: nx, ny
    double precision, dimension(nx, ny), intent(in) :: surfactant
    double precision, dimension(nx, ny), intent(out) :: surface_tension

    integer :: ix, iy

    do iy = 1, ny
        do ix = 1, nx
            ! surface_tension(ix, iy) = (1.d0 + mu * surfactant(ix, iy))**(-3)
            surface_tension(ix, iy) = 1.d0 - surfactant(ix, iy)
        end do
    end do
end subroutine compute_surface_tension


subroutine compute_surface_tension_d1(nx, ny, surfactant, surface_tension_d1)
    implicit none
    double precision :: beta, kappa, delta, mu, right_film_height
    common /surfactant_config/ beta, kappa, delta, mu, right_film_height

    integer, intent(in) :: nx, ny
    double precision, dimension(nx, ny), intent(in) :: surfactant
    double precision, dimension(nx, ny), intent(out) :: surface_tension_d1

    integer :: ix, iy

    do iy = 1, ny
        do ix = 1, nx
            ! surface_tension_d1(ix, iy) = -3.d0 * mu * (1.d0 + mu * surfactant(ix, iy))**(-4)
            surface_tension_d1(ix, iy) = -1.d0
        end do
    end do
end subroutine compute_surface_tension_d1
