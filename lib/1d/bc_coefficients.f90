subroutine get_lower_bc_coefficients(orders, dx, boundary_coefficients, cell_coefficients)
    character*2, intent(in) :: orders
    double precision, intent(in) :: dx
    double precision, dimension(2, 2) :: boundary_coefficients
    double precision, dimension(2, 2) :: cell_coefficients
    
    if (orders == '0') then
        boundary_coefficients(1, 1) = 2.666666666666667d+00
        cell_coefficients(1, 1) = -2.000000000000000d+00
        cell_coefficients(2, 1) = 3.333333333333333d-01
    else if (orders == '1') then
        boundary_coefficients(1, 1) = -1.000000000000000d+00 * dx
        cell_coefficients(1, 1) = 1.000000000000000d+00
        cell_coefficients(2, 1) = 0.000000000000000d+00
    else if (orders == '01') then
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
    
    if (orders == '0') then
        boundary_coefficients(1, 1) = 2.666666666666666d+00
        cell_coefficients(1, 1) = 3.333333333333333d-01
        cell_coefficients(2, 1) = -2.000000000000000d+00
    else if (orders == '1') then
        boundary_coefficients(1, 1) = 1.000000000000000d+00 * dx
        cell_coefficients(1, 1) = 0.000000000000000d+00
        cell_coefficients(2, 1) = 1.000000000000000d+00
    else if (orders == '01') then
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