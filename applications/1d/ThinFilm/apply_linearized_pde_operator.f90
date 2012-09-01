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

    double precision, dimension(4) :: d0_stencil, d1_stencil, d2_stencil,  &
        d3_stencil
    common /stencil_config/ d0_stencil, d1_stencil, d2_stencil, d3_stencil

    double precision, intent(in) :: t
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(in) :: q, p
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(out) :: output

    integer :: ix
    double precision, external :: d0x, d3x


    do ix = 1, mx
        output(ix, 1) = -(fprime(ix+1) - fprime(ix)) / dx
    end do


    contains 

    double precision function fprime(ix)
        implicit none
        integer :: ix
        double precision :: q1_face, p1_face, q1_xxx, p1_xxx
        q1_face = d0x(q(ix-2:ix+1, 1))
        p1_face = d0x(p(ix-2:ix+1, 1))
        q1_xxx = d3x(q(ix-2:ix+1, 1))
        p1_xxx = d3x(p(ix-2:ix+1, 1))
        fprime = q1_xxx * p1_face + q1_face * p1_xxx
    end function fprime

end subroutine apply_linearized_pde_operator
