subroutine apply_newton_operator(t, dt, iterate, d_iterate, output)

!----------------------------------------------------------------
! Calculates
!     m'[iterate](d_iterate) = d_iterate + dt * div(f'[iterate](d_iterate)),
! as appears on the left-hand side of Newton's method.
!
! Args:
!   t: Most recent time value.
!   dt: Time step size.
!   q: Cell-centered PDE solution.
!   d_iterate: What we're applying the Newton operator to.  Homogeneous boundary
!       conditions are applied here.
!   output: Result of the operation.
!----------------------------------------------------------------

    implicit none
    
    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision, intent(in) :: t, dt
    double precision, intent(in), dimension(1-mbc:mx+mbc, meqn) :: iterate
    double precision, intent(inout), dimension(1-mbc:mx+mbc, meqn) :: d_iterate
    double precision, intent(out), dimension(1-mbc:mx+mbc, meqn) :: output

    integer :: ix, ieqn

    
    call apply_homogeneous_bcs(d_iterate)                                                         
    call apply_linearized_pde_operator(t, iterate, d_iterate, output)
    
    do ieqn = 1, meqn
        do ix = 1,mx
            output(ix, ieqn) = d_iterate(ix, ieqn) - dt * output(ix, ieqn)
        end do
    end do
    
end subroutine apply_newton_operator
