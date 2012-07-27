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

    double precision :: gamma
    common /physics_config/ gamma
    
    double precision, intent(in) :: t
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(in), target :: q
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(out) :: output

    double precision, dimension(2-mbc:mx+mbc-1) :: q1_xx

    integer :: ix


    do ix = 0, mx + 1
        q1_xx(ix) = (q(ix-1, 1) - 2 * q(ix, 1) + q(ix+1, 1)) / dx**2
    end do

    do ix = 1, mx
        output(ix, 1) = -(flux(ix+1) - flux(ix)) / dx
    end do


    contains
    
    double precision function flux(ix)
        implicit none
        integer :: ix
        double precision :: q1_face
        q1_face = (q(ix-1, 1) + q(ix, 1)) / 2
        flux = gamma * q1_face**3 * (q1_xx(ix) - q1_xx(ix-1)) / dx
    end function flux

end subroutine apply_pde_operator
