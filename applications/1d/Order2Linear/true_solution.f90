function true_solution(x, t)

    double precision, intent(in) :: x, t
    double precision :: true_solution

    true_solution = exp(-t) * sin(x)
    
end function true_solution
