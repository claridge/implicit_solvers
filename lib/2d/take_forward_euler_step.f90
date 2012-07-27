subroutine take_forward_euler_step(t, dt, q)

! Take a time step using Forward Euler.
!
! Args:
!   t: Time before the step.
!   dt: Length of time step.
!   q: PDE solution at time t.

    !$ use omp_lib
    implicit none

    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn

    double precision, intent(in) :: t, dt
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn),  &
        intent(inout) :: q

    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn) :: op_output
    integer :: ix, ieqn

    call apply_bcs(t, q)
    call apply_pde_operator(t, q, op_output)

    do ieqn = 1, meqn
        !$omp parallel do private(ix)
        do iy = 1, my
            do ix = 1, mx
                q(ix, iy, ieqn) = q(ix, iy, ieqn) + dt * op_output(ix, iy, ieqn)
            end do
        end do
    end do

end subroutine take_forward_euler_step
