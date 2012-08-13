double precision function true_solution(x, y, t)
    double precision :: x, y, t
    double precision :: mass, tmp, tau, t0
    parameter(mass = 1.d0)
    parameter(t0 = .3d0)

    tau = t + t0

    tmp = mass - 0.1875d0 * (x**2 + y**2) / sqrt(tau)
    if (tmp > 0.d0) then
        true_solution = tmp / sqrt(tau)
    else
        true_solution = 0.d0
    end if
end function true_solution
