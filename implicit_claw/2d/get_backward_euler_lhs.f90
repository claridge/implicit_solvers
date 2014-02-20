subroutine get_backward_euler_lhs(t, dt, r, p, output)

!----------------------------------------------------------------
! Calculates
!     m'[r](p) = p + dt * g'[r](p),
! as appears on the left-hand side of Newton's method.
!
! Args:
!   t: Most recent time value.
!   dt: Time step size.
!   r: Current Newton iterate.
!   p: Current Newton perturbation.
!   output: Result of the operation.
!----------------------------------------------------------------

    !$ use omp_lib
    implicit none
    
    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn

    double precision, intent(in) :: t, dt
    double precision, intent(in), dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn) :: r, p    
    double precision, intent(out), dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn) :: output

    integer :: ix, iy, ieqn

    
    call apply_linearized_bcs(r, p)     
    call apply_linearized_pde_operator(t, r, p, output)
    
    do ieqn = 1, meqn
        !$omp parallel do private(ix)
        do iy = 1, my
            do ix = 1,mx
                output(ix, iy, ieqn) = p(ix, iy, ieqn) - dt * output(ix, iy, ieqn)
            end do
        end do
    end do
    
end subroutine get_backward_euler_lhs
