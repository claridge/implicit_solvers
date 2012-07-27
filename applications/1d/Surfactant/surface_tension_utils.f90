subroutine compute_surface_tension(n, surfactant, surface_tension)
    double precision :: mu
    common /physics_config/ mu

    integer, intent(in) :: n
    double precision, dimension(n), intent(in) :: surfactant
    double precision, dimension(n), intent(out) :: surface_tension

    integer :: i

    do i = 1, n
        surface_tension(i) = (1.d0 + mu * surfactant(i))**(-3)
        ! surface_tension(i) = 1.d0 - surfactant(i)
    end do
end subroutine compute_surface_tension


subroutine compute_surface_tension_d1(n, surfactant, surface_tension_d1)
    double precision :: mu
    common /physics_config/ mu

    integer, intent(in) :: n
    double precision, dimension(n), intent(in) :: surfactant
    double precision, dimension(n), intent(out) :: surface_tension_d1

    integer :: i

    do i = 1, n
        surface_tension_d1(i) = -3.d0 * mu * (1.d0 + mu * surfactant(i))**(-4)
        ! surface_tension_d1(i) = -1.d0
    end do
end subroutine compute_surface_tension_d1
