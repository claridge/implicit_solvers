double precision function true_solution(x, y, t)
    implicit none
    double precision, intent(in) :: x, y, t

    double precision :: t0 = .3d0, tau, mass, eta
    parameter(mass = 2.d0)
    integer :: ndim
    parameter(ndim = 2)

    tau = ((ndim + 4) * (t + t0))**(1.d0 / (ndim + 4))
    eta = sqrt(x**2 + y**2) / tau

    if (eta < mass) then
        true_solution = 1.d0/(8 * (ndim + 2) * tau**ndim) * (mass**2 - eta**2)**2
    else
        true_solution = 0.d0
    end if
end function true_solution
