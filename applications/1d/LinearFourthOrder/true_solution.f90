function true_solution(x, t)

! True solution for the thin film equation,
!     q_t + (qq_xxx)_x = 0.

    implicit none
    double precision, intent(in) :: x, t
    double precision :: true_solution

    true_solution = exp(-t) * sin(x)
    
end function true_solution
