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

    double precision, external :: d3x
    integer :: ix


    do ix = 1, mx
        output(ix, 1) = -(flux(ix+1) - flux(ix)) / dx
    end do


    contains
    
    double precision function flux(ix)
        implicit none
        integer :: ix
        double precision :: q1_face, q1_xxx
        q1_xxx = d3x(q(ix-2:ix+1, 1))
        flux = q1_xxx
    end function flux

end subroutine apply_pde_operator


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

    double precision, external :: d3x    
    integer :: ix


    do ix = 1, mx
        output(ix, 1) = -(fprime(ix+1) - fprime(ix)) / dx
    end do


    contains 

    double precision function fprime(ix)
        implicit none
        integer :: ix
        double precision :: q1_face, p1_face, q1_xxx, p1_xxx
        p1_xxx = d3x(p(ix-2:ix+1, 1))
        fprime = p1_xxx
    end function fprime

end subroutine apply_linearized_pde_operator
