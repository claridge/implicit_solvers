subroutine fill_2_ghost_cells(lower_option, lower_values,  &
                              upper_option, upper_values,  &
                              q)

! Fills 1 ghost cell on either side of the domain.
!
! 'option' arguments indicate the condition used to fill the cell.
! Allowable values are:
!   'p': periodic
!
! TODO: Specify allowable values.
!
! If either option is 'p', then both must be.
!
! 'value' arguments are double precision, indicating the value to be
! matched.  If the corresponding option is 'p', this value is unused.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    ! TODO: Get rid of old stencils.
    double precision, dimension(4) :: d0_stencil, d1_stencil, d2_stencil,  &
        d3_stencil
    double precision, dimension(4, 0:3) :: stencils
    common /stencil_config/ d0_stencil, d1_stencil, d2_stencil, d3_stencil, stencils

    ! TODO: Safety check for derivative orders?

    character*2, intent(in) :: lower_option, upper_option
    double precision, intent(in), dimension(2) :: lower_values, upper_values
    double precision, intent(inout) :: q(1-mbc:mx+mbc)

    integer, dimension(2) :: derivative_orders, ipiv
    integer :: info, i
    double precision, dimension(2,2) :: a, b

    if (lower_option(1:1) == 'p' .or. upper_option(1:1) == 'p') then
        if (lower_option(1:1) /= upper_option(1:1)) then
            print *, "Error (fill_2_ghost_cells): If one option is 'p', ",  &
                "then both must be."
        end if

        q(-1:0) = q(mx-1:mx)
        q(mx+1:mx+2) = q(1:2)
        return
    end if

    ! If we're still here, then we're matching derivatives.
    do i = 1, 2
        read(lower_option(i:i), '(I1)') derivative_orders(i)
        read(lower_option(i:i), '(I1)') derivative_orders(i)
    end do

    a(1, 1:2) = stencils(1:2, derivative_orders(1)) / dx**derivative_orders(1)
    b(1, 1:2) = stencils(3:4, derivative_orders(1)) / dx**derivative_orders(1)
    a(2, 1:2) = stencils(1:2, derivative_orders(2)) / dx**derivative_orders(2)
    b(2, 1:2) = stencils(3:4, derivative_orders(2)) / dx**derivative_orders(2)


    q(-1:0) = lower_values - matmul(b, q(1:2))
    call dgesv(2, 1, a, 2, ipiv, q(-1:0), 2, info)

    do i = 1, 2
        read(upper_option(i:i), '(I1)') derivative_orders(i)
        read(upper_option(i:i), '(I1)') derivative_orders(i)
    end do

    a(1, 1:2) = stencils(1:2, derivative_orders(1)) / dx**derivative_orders(1)
    b(1, 1:2) = stencils(3:4, derivative_orders(1)) / dx**derivative_orders(1)
    a(2, 1:2) = stencils(1:2, derivative_orders(2)) / dx**derivative_orders(2)
    b(2, 1:2) = stencils(3:4, derivative_orders(2)) / dx**derivative_orders(2)

    q(mx+1:mx+2) = upper_values - matmul(a, q(mx-1:mx))

    call dgesv(2, 1, b, 2, ipiv, q(mx+1:mx+2), 2, info)

end subroutine fill_2_ghost_cells
