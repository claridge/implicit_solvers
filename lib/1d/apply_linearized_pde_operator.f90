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

    double precision :: epsilon
    double precision, external :: l2_norm
    integer :: ix, ieqn

    double precision :: err, local_err
    integer :: i_max

!     epsilon = 1d-4 * l2_norm(r) / l2_norm(p)
    epsilon = sqrt((1 + l2_norm(r)) * 1d-16) / l2_norm(p)

!     epsilon = 0.d0
!     do ieqn = 1, meqn
!         do ix = 1, mx
!             epsilon = epsilon + abs(r(ix, ieqn))
!         end do
!     end do
!     epsilon = 1.d-4 * (epsilon / (mx * l2_norm(p)) + 1.d0)

    temp1 = r + epsilon * p
    call apply_pde_operator(t, temp1, output)
    temp1 = r - epsilon * p
    call apply_pde_operator(t, temp1, temp2)

    do ieqn = 1, meqn
        do ix = 1, mx
            output(ix, ieqn) = (output(ix, ieqn) - temp2(ix, ieqn)) / (2 * epsilon)
        end do
    end do
    
    
!     temp1 = r + epsilon * p
!     call apply_pde_operator(t, temp1, output)
!     temp1 = r
!     call apply_pde_operator(t, temp1, temp2)
!     
!     do ieqn = 1, meqn
!         do ix = 1, mx
!             output(ix, ieqn) = (output(ix, ieqn) - temp2(ix, ieqn)) / epsilon
!         end do
!     end do

!     ! Hack check for linear problem
!     call apply_pde_operator(t, p, temp1)
!     err = 0.d0
!     do ix = 1, mx
!         local_err = abs(temp1(ix, 1) - output(ix, 1))
!         if (local_err > err) then
!             err = local_err
!             i_max = ix
!         end if
!     end do
!     print *, mx, epsilon, i_max, err
!     print *, r(i_max-2:i_max+2, 1)
!     print *, p(i_max-2:i_max+2, 1)
!     print *, ''


end subroutine apply_linearized_pde_operator
