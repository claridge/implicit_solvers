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

    double precision, dimension(1-mbc:mx+mbc, meqn) :: temp1, temp2

    double precision :: epsilon, norm_p
    double precision, external :: l2_norm
    integer :: ix

    norm_p = l2_norm(p)
    epsilon = 1d-3 / norm_p

    temp1 = q + epsilon * p
    call apply_pde_operator(t, temp1, output)
    temp1 = q - epsilon * p
    call apply_pde_operator(t, temp1, temp2)
    do ix = 1, mx
        output(ix, 1) = (output(ix, 1) - temp2(ix, 1)) / (2 * epsilon)
    end do

end subroutine apply_linearized_pde_operator
