subroutine apply_pde_operator(t, q, output)

! For the PDE q_t = g(q), calculates g(q).
!
! Args:
!   t: Time at which the operator is evaluated.
!   q: PDE solution at time t.
!   output: g(q), calculated here.

    implicit none
    
    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision, intent(in) :: t
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(in), target :: q
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(out) :: output

    double precision, external :: derivative3
    integer :: ix


    do ix = 1, mx
        output(ix, 1) = -(q(ix-2, 1) - 4 * q(ix-1, 1) + 6 * q(ix, 1) - 4 * q(ix+1, 1) + q(ix+2, 1)) / dx**4        
    end do

end subroutine apply_pde_operator