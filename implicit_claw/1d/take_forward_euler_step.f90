subroutine take_forward_euler_step(t, dt, q)

! Take a time step using Forward Euler.
!
! Args:
!   t: Time before the step.
!   dt: Length of time step.
!   q: PDE solution at time t.

    !$ use omp_lib
    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision, intent(in) :: t, dt
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(inout) :: q

    double precision, dimension(1-mbc:mx+mbc, meqn) :: op_output
    integer :: ix, ieqn

    call apply_bcs(t, q)
    call apply_pde_operator(t, q, op_output)

    do ieqn = 1, meqn
        !$omp parallel do
        do ix = 1, mx
            q(ix, ieqn) = q(ix, ieqn) + dt * op_output(ix, ieqn)
        end do
    end do

end subroutine take_forward_euler_step
