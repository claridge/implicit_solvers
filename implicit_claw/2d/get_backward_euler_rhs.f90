subroutine get_backward_euler_rhs(t, dt, q, r, rhs)

! Calculates the RHS for use with backward Euler
!
! Args:
!   t: Time at the beginning of the current step.
!   dt: Length of the current time step.
!   q: Solution values at beginning of current step.
!   r: Current Newton iterate
!   rhs: Returns the RHS.

    !$ use omp_lib
    implicit none

    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn

    double precision, intent(in) :: t, dt
    double precision, intent(in), dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn) :: q, r
    double precision, intent(out), dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn) :: rhs
    integer :: ix, iy, ieqn

    call apply_bcs(t + dt, r)
    call apply_pde_operator(t, r, rhs)

    do ieqn = 1, meqn
        !$omp parallel do private(ix)
        do iy = 1, my
            do ix = 1, mx
                rhs(ix, iy, ieqn) = q(ix, iy, ieqn) - r(ix, iy, ieqn)  &
                    + dt * rhs(ix, iy, ieqn)
            end do
        end do
    end do

end subroutine get_backward_euler_rhs
