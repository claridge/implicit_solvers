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

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision, intent(in) :: t, dt
    double precision, intent(in), dimension(1-mbc:mx+mbc, meqn) :: q, r
    double precision, intent(out), dimension(1-mbc:mx+mbc, meqn) :: rhs
    integer :: ix, ieqn

    call apply_bcs(t + dt, r)
    call apply_pde_operator(t, r, rhs)

    do ieqn = 1, meqn
        !$omp parallel do
        do ix = 1, mx
            rhs(ix, ieqn) = q(ix, ieqn) - r(ix, ieqn) + dt * rhs(ix, ieqn)
        end do
    end do
end subroutine get_backward_euler_rhs
