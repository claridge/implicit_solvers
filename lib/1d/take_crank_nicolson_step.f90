subroutine take_crank_nicolson_step(t, dt, q, success)

! Take a time step using Crank-Nicolson.
!
! Args:
!   t: Time before the step.
!   dt: Length of time step.
!   q: PDE solution being advanced.
!   success: Indicates whether Newton's method in the backward Euler stage
!     converged.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision, intent(in) :: t, dt
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(inout) :: q
    logical, intent(out) :: success

    call take_forward_euler_step(t, dt/2, q)
    call take_backward_euler_step(t + dt/2, dt/2, q, success)

end subroutine take_crank_nicolson_step
