subroutine get_backward_euler_lhs(t, dt, iterate, perturbation, output)

!----------------------------------------------------------------
! Calculates
!     M'[iterate](perturbation) = perturbation + dt * G'[iterate](perturbation),
! as appears on the left-hand side of Newton's method.
!
! Args:
!   t: Most recent time value.
!   dt: Time step size.
!   q: Cell-centered PDE solution.
!   perturbation: What we're applying the Newton operator to.  Homogeneous boundary
!       conditions are applied here.
!   output: Result of the operation.
!------------------------------------------------------------------------------

    !$ use omp_lib
    implicit none
    
    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision, intent(in) :: t, dt
    double precision, intent(in), dimension(1-mbc:mx+mbc, meqn) :: iterate
    double precision, intent(inout), dimension(1-mbc:mx+mbc, meqn) :: perturbation
    double precision, intent(out), dimension(1-mbc:mx+mbc, meqn) :: output

    integer :: ix, ieqn

    
    call apply_homogeneous_bcs(perturbation)                                                         
    call apply_linearized_pde_operator(t, iterate, perturbation, output)
    
    do ieqn = 1, meqn
        !$omp parallel do
        do ix = 1,mx
            output(ix, ieqn) = perturbation(ix, ieqn) - dt * output(ix, ieqn)
        end do
    end do
    
end subroutine get_backward_euler_lhs
