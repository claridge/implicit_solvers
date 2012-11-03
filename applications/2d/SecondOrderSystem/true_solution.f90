subroutine set_true_solution(x, y, t, components)

    implicit none

    double precision, intent(in) :: x, y, t
    double precision, dimension(2), intent(out) :: components
    
    double precision :: t0, tau, rho
    parameter(t0 = 1.d0)
    
    tau = t + t0
    rho = sqrt(x**2 + y**2) / tau**.25d0

    components(1) = 2.d0 * rho**2
    components(2) = -1.d0 / (8.d0 * tau**.5d0) * log(rho)

end subroutine set_true_solution