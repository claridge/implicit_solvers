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

    character*2, intent(in) :: lower_option, upper_option
    double precision, intent(in), dimension(2) :: lower_values, upper_values
    double precision, intent(inout) :: q(1-mbc:mx+mbc)

    integer, dimension(2) :: derivative_orders, ipiv
    integer :: info, i
    double precision, dimension(2,2) ::boundary_coeffs, cell_coeffs

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


    call get_lower_bc_coefficients(lower_option, dx, boundary_coeffs, cell_coeffs)
    q(-1) = boundary_coeffs(1, 1) * lower_values(1) + boundary_coeffs(2, 1) * lower_values(2)  &
        + cell_coeffs(1, 1) * q(1) + cell_coeffs(2, 1) * q(2)
    q(0) = boundary_coeffs(1, 2) * lower_values(1) + boundary_coeffs(2, 2) * lower_values(2)  &
        + cell_coeffs(1, 2) * q(1) + cell_coeffs(2, 2) * q(2)
    

    call get_upper_bc_coefficients(upper_option, dx, boundary_coeffs, cell_coeffs)
    q(mx+1) = boundary_coeffs(1, 1) * upper_values(1) + boundary_coeffs(2, 1) * upper_values(2)  &
        + cell_coeffs(1, 1) * q(mx-1) + cell_coeffs(2, 1) * q(mx)
    q(mx+2) = boundary_coeffs(1, 2) * upper_values(1) + boundary_coeffs(2, 2) * upper_values(2)  &
        + cell_coeffs(1, 2) * q(mx-1) + cell_coeffs(2, 2) * q(mx)
    
end subroutine fill_2_ghost_cells


subroutine get_lower_bc_coefficients(orders, dx, boundary_coefficients, cell_coefficients)
    character*2, intent(in) :: orders
    double precision, intent(in) :: dx
    double precision, dimension(2, 2) :: boundary_coefficients
    double precision, dimension(2, 2) :: cell_coefficients
    
    if (orders == '01') then
        boundary_coefficients(1, 1) = -2.400000000000000d+01
        boundary_coefficients(2, 1) = -1.200000000000000d+01 * dx
        cell_coefficients(1, 1) = 2.700000000000000d+01
        cell_coefficients(2, 1) = -2.000000000000000d+00
        boundary_coefficients(1, 2) = -8.888888888888888d-01
        boundary_coefficients(2, 2) = -1.333333333333333d+00 * dx
        cell_coefficients(1, 2) = 2.000000000000000d+00
        cell_coefficients(2, 2) = -1.111111111111111d-01
    else if (orders == '02') then
        boundary_coefficients(1, 1) = 2.000000000000000d+00
        boundary_coefficients(2, 1) = 2.250000000000000d+00 * dx**2
        cell_coefficients(1, 1) = 0.000000000000000d+00
        cell_coefficients(2, 1) = -1.000000000000000d+00
        boundary_coefficients(1, 2) = 2.000000000000000d+00
        boundary_coefficients(2, 2) = 2.500000000000000d-01 * dx**2
        cell_coefficients(1, 2) = -1.000000000000000d+00
        cell_coefficients(2, 2) = 0.000000000000000d+00
    else if (orders == '03') then
        boundary_coefficients(1, 1) = 8.000000000000000d+00
        boundary_coefficients(2, 1) = -1.500000000000000d+00 * dx**3
        cell_coefficients(1, 1) = -9.000000000000000d+00
        cell_coefficients(2, 1) = 2.000000000000000d+00
        boundary_coefficients(1, 2) = 2.666666666666667d+00
        boundary_coefficients(2, 2) = -1.666666666666667d-01 * dx**3
        cell_coefficients(1, 2) = -2.000000000000000d+00
        cell_coefficients(2, 2) = 3.333333333333333d-01
    else if (orders == '12') then
        boundary_coefficients(1, 1) = -9.230769230769227d-01 * dx
        boundary_coefficients(2, 1) = 2.076923076923077d+00 * dx**2
        cell_coefficients(1, 1) = 2.076923076923077d+00
        cell_coefficients(2, 1) = -1.076923076923077d+00
        boundary_coefficients(1, 2) = -9.230769230769229d-01 * dx
        boundary_coefficients(2, 2) = 7.692307692307693d-02 * dx**2
        cell_coefficients(1, 2) = 1.076923076923077d+00
        cell_coefficients(2, 2) = -7.692307692307693d-02
    else if (orders == '13') then
        boundary_coefficients(1, 1) = -3.000000000000000d+00 * dx
        boundary_coefficients(2, 1) = -1.125000000000000d+00 * dx**3
        cell_coefficients(1, 1) = 0.000000000000000d+00
        cell_coefficients(2, 1) = 1.000000000000000d+00
        boundary_coefficients(1, 2) = -1.000000000000000d+00 * dx
        boundary_coefficients(2, 2) = -4.166666666666666d-02 * dx**3
        cell_coefficients(1, 2) = 1.000000000000000d+00
        cell_coefficients(2, 2) = 0.000000000000000d+00
    else if (orders == '23') then
        boundary_coefficients(1, 1) = 3.000000000000000d+00 * dx**2
        boundary_coefficients(2, 1) = 5.000000000000000d-01 * dx**3
        cell_coefficients(1, 1) = 3.000000000000000d+00
        cell_coefficients(2, 1) = -2.000000000000000d+00
        boundary_coefficients(1, 2) = 1.000000000000000d+00 * dx**2
        boundary_coefficients(2, 2) = 5.000000000000000d-01 * dx**3
        cell_coefficients(1, 2) = 2.000000000000000d+00
        cell_coefficients(2, 2) = -1.000000000000000d+00
    end if
    
end subroutine get_lower_bc_coefficients    


subroutine get_upper_bc_coefficients(orders, dx, boundary_coefficients, cell_coefficients)
    character*2, intent(in) :: orders
    double precision, intent(in) :: dx
    double precision, dimension(2, 2) :: boundary_coefficients
    double precision, dimension(2, 2) :: cell_coefficients
    
    if (orders == '01') then
        boundary_coefficients(1, 1) = -8.888888888888884d-01
        boundary_coefficients(2, 1) = 1.333333333333333d+00 * dx
        cell_coefficients(1, 1) = -1.111111111111111d-01
        cell_coefficients(2, 1) = 2.000000000000000d+00
        boundary_coefficients(1, 2) = -2.400000000000000d+01
        boundary_coefficients(2, 2) = 1.200000000000000d+01 * dx
        cell_coefficients(1, 2) = -2.000000000000000d+00
        cell_coefficients(2, 2) = 2.700000000000000d+01
    else if (orders == '02') then
        boundary_coefficients(1, 1) = 2.000000000000000d+00
        boundary_coefficients(2, 1) = 2.500000000000000d-01 * dx**2
        cell_coefficients(1, 1) = -0.000000000000000d+00
        cell_coefficients(2, 1) = -1.000000000000000d+00
        boundary_coefficients(1, 2) = 2.000000000000000d+00
        boundary_coefficients(2, 2) = 2.250000000000000d+00 * dx**2
        cell_coefficients(1, 2) = -1.000000000000000d+00
        cell_coefficients(2, 2) = -0.000000000000000d+00
    else if (orders == '03') then
        boundary_coefficients(1, 1) = 2.666666666666666d+00
        boundary_coefficients(2, 1) = 1.666666666666666d-01 * dx**3
        cell_coefficients(1, 1) = 3.333333333333333d-01
        cell_coefficients(2, 1) = -2.000000000000000d+00
        boundary_coefficients(1, 2) = 7.999999999999998d+00
        boundary_coefficients(2, 2) = 1.500000000000000d+00 * dx**3
        cell_coefficients(1, 2) = 2.000000000000000d+00
        cell_coefficients(2, 2) = -8.999999999999998d+00
    else if (orders == '12') then
        boundary_coefficients(1, 1) = 9.230769230769229d-01 * dx
        boundary_coefficients(2, 1) = 7.692307692307693d-02 * dx**2
        cell_coefficients(1, 1) = -7.692307692307691d-02
        cell_coefficients(2, 1) = 1.076923076923077d+00
        boundary_coefficients(1, 2) = 9.230769230769229d-01 * dx
        boundary_coefficients(2, 2) = 2.076923076923077d+00 * dx**2
        cell_coefficients(1, 2) = -1.076923076923077d+00
        cell_coefficients(2, 2) = 2.076923076923077d+00
    else if (orders == '13') then
        boundary_coefficients(1, 1) = 1.000000000000000d+00 * dx
        boundary_coefficients(2, 1) = 4.166666666666666d-02 * dx**3
        cell_coefficients(1, 1) = 0.000000000000000d+00
        cell_coefficients(2, 1) = 1.000000000000000d+00
        boundary_coefficients(1, 2) = 3.000000000000000d+00 * dx
        boundary_coefficients(2, 2) = 1.125000000000000d+00 * dx**3
        cell_coefficients(1, 2) = 1.000000000000000d+00
        cell_coefficients(2, 2) = 0.000000000000000d+00
    else if (orders == '23') then
        boundary_coefficients(1, 1) = 9.999999999999998d-01 * dx**2
        boundary_coefficients(2, 1) = -4.999999999999998d-01 * dx**3
        cell_coefficients(1, 1) = -9.999999999999998d-01
        cell_coefficients(2, 1) = 2.000000000000000d+00
        boundary_coefficients(1, 2) = 3.000000000000000d+00 * dx**2
        boundary_coefficients(2, 2) = -5.000000000000000d-01 * dx**3
        cell_coefficients(1, 2) = -2.000000000000000d+00
        cell_coefficients(2, 2) = 3.000000000000000d+00
    end if

end subroutine get_upper_bc_coefficients    