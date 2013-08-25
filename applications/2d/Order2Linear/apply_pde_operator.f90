subroutine apply_pde_operator(t, q, output)

! For the PDE q_t = g(q), calculates g(q).
!
! Args:
!   t: Time at which the operator is evaluated.
!   q: PDE solution at time t.
!   output: g(q), calculated here.

    !$ use omp_lib
    implicit none
    
    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn
    
    double precision, intent(in) :: t
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn), intent(in) :: q
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn), intent(out) :: output

    integer :: ix, iy

    !$omp parallel do private(ix)
    do iy = 1, my
        do ix = 1, mx
            output(ix, iy, 1) =  &
                (q(ix+1, iy, 1) - 2*q(ix, iy, 1) + q(ix-1, iy, 1)) / dx**2 +  &
                (q(ix, iy+1, 1) - 2*q(ix, iy, 1) + q(ix, iy-1, 1)) / dy**2
        end do
    end do

end subroutine apply_pde_operator
