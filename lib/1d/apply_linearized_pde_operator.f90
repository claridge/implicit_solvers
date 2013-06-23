subroutine apply_linearized_pde_operator(t, r, p, output)

! For the PDE r_t = g(r), calculates g'[r](p), i.e. the linearization
! of g, about r, applied to perturbation p
!
! Args:
!   t: Time at which the operator is evaluated.
!   r: Base function of the linearization.
!   p: Perturbation to which the linearized operator is applied.
!   output: g'[r](p), calculated here.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: dx, x_lower
    common /claw_config/ mx, mbc, x_lower, dx, meqn
    
    double precision, intent(in) :: t
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(in) :: r, p
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(out) :: output

    double precision, dimension(1-mbc:mx+mbc, meqn) :: temp1, temp2

    double precision :: epsilon, norm_r, norm_p
    double precision, external :: l2_norm
    integer :: ix, ieqn

    norm_r = l2_norm(r)
    norm_p = l2_norm(p)
    
    epsilon = 1d-4 * norm_r / norm_p

    temp1 = r + epsilon * p
    call apply_pde_operator(t, temp1, output)
    temp1 = r - epsilon * p
    call apply_pde_operator(t, temp1, temp2)

    do ieqn = 1, meqn
        do ix = 1, mx
            output(ix, ieqn) = (output(ix, ieqn) - temp2(ix, ieqn)) / (2 * epsilon)
        end do
    end do

end subroutine apply_linearized_pde_operator
