subroutine calculate_newton_rhs(dt, density, iterate, rhs)

!------------------------------------------------------------------------------
! Evaluates the RHS for Newton's method applied to the porous medium equation,
!     rhs = density - iterate - dt * div(porous_flux(iterate)).
!
! Args:
!   dt: Time step size
!   iterate: The most recent Newton iterate
!   density: Cell-centered density from most recent time step
!   rhs: Output
!
! TODO: I'd like to be more descrptive than 'iterate' and 'density', but the
! most obvious descriptors to me are 'current' and 'previous', which are
! ambiguous.
!------------------------------------------------------------------------------

    implicit none

    integer :: mx, mbc
    double precision :: x_lower, dx
    common /grid_data/ mx, mbc, x_lower, dx

    double precision, intent(in) :: dt
    double precision, intent(in), dimension(1-mbc:mx+mbc) :: iterate, density
    double precision, intent(out), dimension(1-mbc:mx+mbc) :: rhs
    
    integer :: ix


    call div_flux(iterate, rhs)

    do ix = 1, mx
        rhs(ix) = density(ix) - iterate(ix) - dt*rhs(ix)
    end do

end subroutine calculate_newton_rhs
