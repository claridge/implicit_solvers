subroutine apply_pde_operator(t, q, output)

! For the PDE q_t = g(q), calculates g(q).
!
! Args:
!   t: Time at which the operator is evaluated.
!   q: PDE solution at time t.
!   output: g(q), calculated here.

    !$ use omp_lib
    implicit none
    
    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn
    
    double precision, intent(in) :: t
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(in), target :: q
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(out) :: output

    integer :: ix
    double precision, dimension(1:mx+1) :: flux
    double precision :: q1_face, q1_xxx

    !$omp parallel do
    do ix = 1, mx+1
        q1_face = (q(ix-1, 1) + q(ix, 1)) / 2.d0
        q1_xxx =  &
            (-q(ix-2, 1) + 3 * q(ix-1, 1) - 3 * q(ix, 1) + q(ix+1, 1))  &
            / dx**3
        flux(ix) = q1_face * q1_xxx
    end do

    !$omp parallel do
    do ix = 1, mx
        output(ix, 1) = -(flux(ix+1) - flux(ix)) / dx
    end do

end subroutine apply_pde_operator
