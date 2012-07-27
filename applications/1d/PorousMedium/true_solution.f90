function true_solution(x, t)

! True solution for the porous medium equation,
!     q_t = (qq_x)_x.

    double precision, intent(in) :: x, t
    double precision :: true_solution
    double precision :: t0 = .3d0, tau, one_third
    parameter(one_third = 1.d0/3.d0)
    double precision :: a = (4.5d0)**one_third


    tau = t + t0

    if (abs(x) <= a * tau**one_third) then
        true_solution = 1.d0 / (6 * tau) * ((a * tau**one_third)**2 - x**2)
    else
        true_solution = 0.d0
    end if

end function true_solution
