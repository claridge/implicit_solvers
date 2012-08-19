double precision function true_solution(x, y, t)

! True solution for the porous medium equation,
!     q_t = (q**2)_xx.

    double precision :: x, y, t
    double precision :: mass, tau, t0
    parameter(mass = 1.d0)
    parameter(t0 = .3d0)

    tau = t + t0
    true_solution = (mass - 1.d0/16.d0 * (x**2 + y**2) / sqrt(tau)) / sqrt(tau)
    if (true_solution < 0.d0) true_solution = 0.d0
end function true_solution
