subroutine get_backward_euler_lhs(t, dt, r, p, output)

!----------------------------------------------------------------
! Calculates
!     M'[r](p) = p + dt * G'[r](p),
! as appears on the left-hand side of Newton's method.
!
! Args:
!   t: Most recent time value.
!   dt: Time step size.
!   r: Newton iterate.
!   p: Newton perturbation.
!   output: Result of the operation.
!------------------------------------------------------------------------------

    !$ use omp_lib
    implicit none
    
    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision, intent(in) :: t, dt
    double precision, intent(in), dimension(1-mbc:mx+mbc, meqn) :: r
    double precision, intent(inout), dimension(1-mbc:mx+mbc, meqn) :: p
    double precision, intent(out), dimension(1-mbc:mx+mbc, meqn) :: output

    integer :: ix, ieqn

    
    call apply_linearized_bcs(r, p)                                                         
    call apply_linearized_pde_operator(t, r, p, output)
    
    do ieqn = 1, meqn
        !$omp parallel do
        do ix = 1,mx
            output(ix, ieqn) = p(ix, ieqn) - dt * output(ix, ieqn)
        end do
    end do
    
end subroutine get_backward_euler_lhs
