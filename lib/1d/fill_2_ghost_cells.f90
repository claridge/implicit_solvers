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
            stop
        end if

        q(-1:0) = q(mx-1:mx)
        q(mx+1:mx+2) = q(1:2)
        return
    end if

    if (lower_option == '01') then
        q(0) = 8.888888888888888d-01*lower_values(1)  &
            - 1.333333333333333d+00*dx*lower_values(2)  &
            + 2.d+00*q(1)  &
            - 1.111111111111111d-01*q(2)
        q(-1) = 2.4d+01*lower_values(1)  &
            - 1.2d+01*dx*lower_values(2)  &
            + 2.7d+01*q(1)  &
            - 2.d+00*q(2)
    else if (lower_option == '02') then
        q(0) = 2.d+00*lower_values(1)  &
            + 2.5d-01*dx**2*lower_values(2)  &
            - 1.d+00*q(1)
        q(-1) = 2.d+00*lower_values(1)  &
            + 2.25d+00*dx**2*lower_values(2)  &
            - 1.d+00*q(2)
    else if (lower_option == '03') then
        q(0) = 2.666666666666667d+00*lower_values(1)  &
            - 1.666666666666667d-01*dx**3*lower_values(2)  &
            - 2.d+00*q(1)  &
            + 3.333333333333333d-01*q(2)
        q(-1) = 8.d+00*lower_values(1)  &
            - 1.5d+00*dx**3*lower_values(2)  &
            - 9.d+00*q(1)  &
            + 2.d+00*q(2)
    else if (lower_option == '12') then
        q(0) = 9.230769230769229d-01*dx*lower_values(1)  &
            + 7.692307692307693d-02*dx**2*lower_values(2)  &
            + 1.076923076923077d+00*q(1)  &
            - 7.692307692307693d-02*q(2)
        q(-1) = 9.230769230769227d-01*dx*lower_values(1)  &
            + 2.076923076923077d+00*dx**2*lower_values(2)  &
            + 2.076923076923077d+00*q(1)  &
            - 1.076923076923077d+00*q(2)
    else if (lower_option == '13') then
        q(0) = 1.d+00*dx*lower_values(1)  &
            - 4.166666666666666d-02*dx**3*lower_values(2)  &
            + 1.d+00*q(1)
        q(-1) = 3.d+00*dx*lower_values(1)  &
            - 1.125d+00*dx**3*lower_values(2)  &
            + 1.d+00*q(2)
    else if (lower_option == '23') then
        q(0) = 1.d+00*dx**2*lower_values(1)  &
            + 5.d-01*dx**3*lower_values(2)  &
            + 2.d+00*q(1)  &
            - 1.d+00*q(2)
        q(-1) = 3.d+00*dx**2*lower_values(1)  &
            + 5.d-01*dx**3*lower_values(2)  &
            + 3.d+00*q(1)  &
            - 2.d+00*q(2)
    end if

    if (upper_option == '01') then
        q(mx+1) = 8.888888888888884d-01*upper_values(1)  &
            + 1.333333333333333d+00*dx*upper_values(2)  &
            - 1.111111111111111d-01*q(mx-1)  &
            + 2.d+00*q(mx)
        q(mx+2) = 2.4d+01*upper_values(1)  &
            + 1.2d+01*dx*upper_values(2)  &
            - 2.d+00*q(mx-1)  &
            + 2.7d+01*q(mx)
    else if (upper_option == '02') then
        q(mx+1) = 2.d+00*upper_values(1)  &
            + 2.5d-01*dx**2*upper_values(2)  &
            - 1.d+00*q(mx)
        q(mx+2) = 2.d+00*upper_values(1)  &
            + 2.25d+00*dx**2*upper_values(2)  &
            - 1.d+00*q(mx-1)
    else if (upper_option == '03') then
        q(mx+1) = 2.666666666666666d+00*upper_values(1)  &
            + 1.666666666666666d-01*dx**3*upper_values(2)  &
            + 3.333333333333333d-01*q(mx-1)  &
            - 2.d+00*q(mx)
        q(mx+2) = 7.999999999999998d+00*upper_values(1)  &
            + 1.5d+00*dx**3*upper_values(2)  &
            + 2.d+00*q(mx-1)  &
            - 8.999999999999998d+00*q(mx)
    else if (upper_option == '12') then
        q(mx+1) = 9.230769230769229d-01*dx*upper_values(1)  &
            + 7.692307692307693d-02*dx**2*upper_values(2)  &
            - 7.692307692307691d-02*q(mx-1)  &
            + 1.076923076923077d+00*q(mx)
        q(mx+2) = 9.230769230769229d-01*dx*upper_values(1)  &
            + 2.076923076923077d+00*dx**2*upper_values(2)  &
            - 1.076923076923077d+00*q(mx-1)  &
            + 2.076923076923077d+00*q(mx)
    else if (upper_option == '13') then
        q(mx+1) = 1.d+00*dx*upper_values(1)  &
            + 4.166666666666666d-02*dx**3*upper_values(2)  &
            + 1.d+00*q(mx)
        q(mx+2) = 3.d+00*dx*upper_values(1)  &
            + 1.125d+00*dx**3*upper_values(2)  &
            + 1.d+00*q(mx-1)
    else if (upper_option == '23') then
        q(mx+1) = 9.999999999999998d-01*dx**2*upper_values(1)  &
            - 4.999999999999998d-01*dx**3*upper_values(2)  &
            - 9.999999999999998d-01*q(mx-1)  &
            + 2.d+00*q(mx)
        q(mx+2) = 3.d+00*dx**2*upper_values(1)  &
            - 5.d-01*dx**3*upper_values(2)  &
            - 2.d+00*q(mx-1)  &
            + 3.d+00*q(mx)
    end if

end subroutine fill_2_ghost_cells
