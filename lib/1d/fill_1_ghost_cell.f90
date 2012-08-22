subroutine fill_1_ghost_cell(lower_option, lower_value,  &
                             upper_option, upper_value,  &
                             q)

! Fills 1 ghost cell on either side of the domain.
!
! 'option' arguments indicate the condition used to fill the cell.
! Allowable values are:
!   'p': periodic
!   '0': 0th derivative
!   '1': 1st derivative
! If either option is 'p', then both must be.
!
! 'value' arguments are double precision, indicating the value to be
! matched.  If the corresponding option is 'p', this value is unused.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    character*1, intent(in) :: lower_option, upper_option
    double precision, intent(in) :: lower_value, upper_value
    double precision, intent(inout) :: q(1-mbc:mx+mbc)


    if (lower_option == 'p') then
        if (upper_option .ne. 'p') then
            print *, "Error (fill_1_ghost_cell): If one option is 'p', ",  &
                "then both must be."
            stop
        end if
        q(0) = q(mx)
    else if (lower_option == '0') then
        q(0) = 2 * lower_value - q(1)
    else if (lower_option == '1') then
        q(0) = q(1) - lower_value * dx
    else
        print *, "Error (fill_1_ghost_cell): Invalid value for 'lower_option'."
    end if

    if (upper_option == 'p') then
        if (lower_option .ne. 'p') then
            print *, "Error (fill_1_ghost_cell): If one option is 'p', ",  &
                "then both must be."
            stop
        end if
        q(mx+1) = q(1)
    else if (upper_option == '0') then
        q(mx+1) = 2 * upper_value - q(mx)
    else if (upper_option == '1') then
        q(mx+1) = q(mx) + dx * upper_value
    else
        print *, "Error (fill_1_ghost_cell): Invalid value for 'upper_option'."
    end if

end subroutine fill_1_ghost_cell
