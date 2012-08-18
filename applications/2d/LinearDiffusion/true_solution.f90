double precision function true_solution(x, y, t)
    double precision :: x, y, t

    ! true_solution = exp(-t) * sin(y)
    true_solution = exp(-2*t) * sin(x) * sin(y)
end function true_solution
