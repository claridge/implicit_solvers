program bc_test

    implicit none
    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision :: q(-1:102, 1), x
    double precision, external :: true_solution
    integer :: ix

    mx = 100
    mbc = 2
    meqn = 1
    x_lower = .1
    dx = .01

    do ix = 1, mx
        x = x_lower + (ix - .5d0) * dx
        q(ix, 1) = true_solution(x, 0.d0)
    end do

    call apply_bcs(0.d0, q)

    print *, q(-1:2, 1)
    print *, ''
    print *, q(mx-1:mx+2, 1)

end program bc_test


double precision function true_solution(x, t)
    implicit none
    double precision, intent(in) :: x, t
    true_solution = sin(x)
end function true_solution
