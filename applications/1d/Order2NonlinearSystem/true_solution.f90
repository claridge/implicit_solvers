subroutine set_true_solution(x, t, components)
    double precision, intent(in) :: x, t
    double precision, dimension(2) :: components
    double precision :: rho, tau
    parameter(t0 = .3d0)
    
    tau = t + t0
    rho = x / tau**(1.d0/3.d0)
    components(1) = 2.d0 * rho
    components(2) = 1 / (6.d0 * tau**(1.d0/3.d0)) * (1.d0 - rho)    
end subroutine set_true_solution