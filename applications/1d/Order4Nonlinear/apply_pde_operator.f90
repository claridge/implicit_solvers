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
    double precision, external :: derivative0, derivative3

    integer :: ix

    !$omp parallel do
    do ix = 1, mx
        output(ix, 1) = -(flux(ix+1) - flux(ix)) / dx
    end do


    contains
    
    double precision function flux(ix)
        implicit none
        integer :: ix
        double precision :: q1_face, q1_xxx
        q1_face = derivative0(q(ix-2:ix+1, 1))
        q1_xxx = derivative3(q(ix-2:ix+1, 1), dx)
        flux = q1_face * q1_xxx
    end function flux

end subroutine apply_pde_operator
