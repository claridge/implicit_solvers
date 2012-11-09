double precision function true_solution(x, t)

! True solution for the porous medium equation,
!     q_t = (q**2)_xx.

    double precision :: x, t
    double precision :: mass, tau, t0
    parameter(mass = 1.d0)
    parameter(t0 = .3d0)

    tau = t + t0
    true_solution = (mass - 1.d0/12.d0 * x**2 / tau**(2.d0/3.d0)) / tau**(1.d0/3.d0)
    if (true_solution < 0.d0) true_solution = 0.d0
end function true_solution
