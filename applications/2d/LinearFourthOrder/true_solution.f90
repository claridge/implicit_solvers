double precision function true_solution(x, y, t)
    double precision :: x, y, t
    true_solution = exp(-4*t) * sin(x) * sin(y)
!     true_solution = exp(-t) * sin(x)
end function true_solution
