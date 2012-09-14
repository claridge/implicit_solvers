subroutine bc2(maxmx, maxmy, meqn, mbc, mx, my, x_low, y_low, dx, dy,  &
               q, maux, aux, t, dt, mthbc)

!-----------------------------------------------------------
! Applies boundary conditions before each call to rpn2.
!-----------------------------------------------------------

    implicit none

    integer, intent(in) :: maxmx, maxmy, meqn, mbc, mx, my, maux
    double precision, intent(in) :: x_low, y_low, dx, dy, t, dt
    double precision :: aux(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, *)
    integer, dimension(4), intent(in) :: mthbc(4)
    double precision, dimension(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn),  &
        intent(inout) :: q
    
    ! This fills ghost cells using the same conditions as the implicit solver.
    call apply_bcs(t, q)

end subroutine bc2
