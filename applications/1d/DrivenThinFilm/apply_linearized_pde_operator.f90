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

    double precision :: gamma
    common /physics_config/ gamma
    
    double precision, intent(in) :: t
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(in) :: q, p
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(out) :: output
    
    integer :: ix
    double precision, dimension(2-mbc:mx+mbc-1) :: q1_xx, p1_xx

    do ix = 0, mx + 1
        q1_xx(ix) = (q(ix-1, 1) - 2 * q(ix, 1) + q(ix+1, 1)) / dx**2
    end do

    do ix = 0, mx + 1
        p1_xx(ix) = (p(ix-1, 1) - 2 * p(ix, 1) + p(ix+1, 1)) / dx**2
    end do

    do ix = 1, mx
        output(ix, 1) = -(fprime(ix+1) - fprime(ix)) / dx
    end do


    contains 

    double precision function fprime(ix)
        implicit none
        integer :: ix
        double precision :: q1_face, p1_face, q1_xxx, p1_xxx
        q1_face = (q(ix-1, 1) + q(ix, 1)) / 2
        p1_face = (p(ix-1, 1) + p(ix, 1)) / 2
        q1_xxx = (q1_xx(ix) - q1_xx(ix-1)) / dx
        p1_xxx = (p1_xx(ix) - p1_xx(ix-1)) / dx
        fprime = gamma * (3 * q1_face**2 * q1_xxx * p1_face + q1_face**3 * p1_xxx)
    end function fprime
    
    
end subroutine apply_linearized_pde_operator
