function true_solution(x, t)

! True solution for the thin film equation,
!     q_t + (qq_xxx)_x = 0.

    implicit none
    double precision, intent(in) :: x, t
    double precision :: true_solution

    double precision :: t0 = .3d0, tau, ell, eta
    parameter(ell = 2.d0)

    tau = (5.d0 * (t + t0))**.2d0
    eta = x / tau

    if (eta < ell) then
        true_solution = 1.d0/(24.d0 * tau) * (ell**2 - eta**2)**2
    else
        true_solution = 0.d0
    end if
    
end function true_solution
