subroutine apply_linearized_pde_operator(t, r, p, output)

! For the PDE q_t = g(r), calculates g'[r](p), i.e. the linearization
! of g, about r, applied to perturbation p
!
! Args:
!   t: Time at which the operator is evaluated.
!   r: Base function of the linearization.
!   p: Perturbation to which the linearized operator is applied.
!   output: g'[r](p), calculated here.

    !$ use omp_lib
    implicit none

    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn
    
    double precision, intent(in) :: t
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn), intent(in) :: r, p
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn), intent(out) :: output

    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn) :: temp1, temp2

    double precision :: epsilon, norm_r, norm_p
    double precision, external :: l2_norm
    
    integer :: ix, iy, ieqn

    norm_r = l2_norm(r)
    norm_p = l2_norm(p)
    epsilon = 1d-4 * norm_r / norm_p

    temp1 = r + epsilon * p
    call apply_pde_operator(t, temp1, output)
    temp1 = r - epsilon * p
    call apply_pde_operator(t, temp1, temp2)
    
    do ieqn = 1, meqn
        do iy = 1, my
            do ix = 1, mx
                output(ix, iy, ieqn) = (output(ix, iy, ieqn) - temp2(ix, iy, ieqn)) / (2 * epsilon)
            end do
        end do
    end do
        
end subroutine apply_linearized_pde_operator


double precision function l2_norm(r)
    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn) :: r
    double precision, external :: inner_product
    l2_norm = sqrt(inner_product(r, r))
end function l2_norm