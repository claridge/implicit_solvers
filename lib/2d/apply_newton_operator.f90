subroutine apply_newton_operator(t, dt, iterate, d_iterate, output)

!----------------------------------------------------------------
! Calculates
!     m'[iterate](d_iterate) = d_iterate + dt * g'[iterate](d_iterate),
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

    !$ use omp_lib
    implicit none
    
    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn

    double precision, intent(in) :: t, dt
    double precision, intent(in), dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn) :: iterate
    double precision, intent(inout), dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn) :: d_iterate
    double precision, intent(out), dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn) :: output

    integer :: ix, iy, ieqn

    
    call apply_homogeneous_bcs(d_iterate)                                                         
    call apply_linearized_pde_operator(t, iterate, d_iterate, output)
    
    do ieqn = 1, meqn
        !$omp parallel do private(ix)
        do iy = 1, my
            do ix = 1,mx
                output(ix, iy, ieqn) = d_iterate(ix, iy, ieqn) - dt * output(ix, iy, ieqn)
            end do
        end do
    end do
    
end subroutine apply_newton_operator
