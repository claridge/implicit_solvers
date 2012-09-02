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
