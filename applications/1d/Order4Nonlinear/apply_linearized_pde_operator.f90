subroutine apply_linearized_pde_operator(t, q, p, output)

! For the PDE q_t = g(q), calculates g'[q](p), i.e. the linearization
! of g, about q, applied to perturbation p
!
! Args:
!   t: Time at which the operator is evaluated.
!   q: Base function of the linearization.
!   p: Perturbation to which the linearized operator is applied.
!   output: g'[q](p), calculated here.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: dx, x_lower
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision, intent(in) :: t
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(in) :: q, p
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(out) :: output

    integer :: ix
    double precision, dimension(1:mx+1) :: fprime  ! Linearized flux
    double precision :: q1_face, p1_face, q1_xxx, p1_xxx

    ! Compute linearized flux at interfaces.
    !$omp parallel do
    do ix = 1, mx+1
        q1_face = (q(ix-1, 1) + q(ix, 1)) / 2.d0
        p1_face = (p(ix-1, 1) + p(ix, 1)) / 2.d0
        q1_xxx =  &
            (-q(ix-2, 1) + 3 * q(ix-1, 1) - 3 * q(ix, 1) + q(ix+1, 1))  &
            / dx**3
        p1_xxx =  &
            (-p(ix-2, 1) + 3 * p(ix-1, 1) - 3 * p(ix, 1) + p(ix+1, 1))  &
            / dx**3
        fprime(ix) = q1_xxx * p1_face + q1_face * p1_xxx                
    end do

    !$omp parallel do
    do ix = 1, mx
        output(ix, 1) = -(fprime(ix+1) - fprime(ix)) / dx
    end do

end subroutine apply_linearized_pde_operator
