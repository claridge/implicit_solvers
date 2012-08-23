double precision function true_solution(x, t)

! True solution for the thin film equation,
!     q_t + (qq_xxx)_x = 0.

    implicit none
    double precision, intent(in) :: x, t

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


subroutine get_true_solution_by_index_range(i_low, i_high, t, values)

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    integer, intent(in) :: i_low, i_high
    double precision, intent(in) :: t
    double precision, dimension(i_low:i_high) :: values

    integer :: i
    double precision :: x
    double precision, external :: true_solution

    do i = i_low, i_high
        x = x_lower + (i - .5d0) * dx
        values(i) = true_solution(x, t)
    end do

end subroutine get_true_solution_by_index_range
