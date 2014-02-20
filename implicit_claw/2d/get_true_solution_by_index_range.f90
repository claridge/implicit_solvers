subroutine get_true_solution_by_index_range(ix_low, ix_high, iy_low, iy_high, t, values)

    implicit none

    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn

    integer, intent(in) :: ix_low, ix_high, iy_low, iy_high
    double precision, intent(in) :: t
    double precision, dimension(ix_low:ix_high, iy_low:iy_high) :: values

    integer :: ix, iy
    double precision :: x, y
    double precision, external :: true_solution

    do iy = iy_low, iy_high
        y = y_lower + (iy - .5d0) * dy
        do ix = ix_low, ix_high
            x = x_lower + (ix - .5d0) * dx
            values(ix, iy) = true_solution(x, y, t)
        end do
    end do

end subroutine get_true_solution_by_index_range
