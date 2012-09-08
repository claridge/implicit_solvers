subroutine get_lower_bc_coefficients(orders, dx, boundary_coefficients, cell_coefficients)
    character*2, intent(in) :: orders
    double precision, intent(in) :: dx
    double precision, dimension(2, 2) :: boundary_coefficients
    double precision, dimension(3, 2) :: cell_coefficients
    

    if (orders == '0') then
        boundary_coefficients(1, 1) = 1.280000000000000d+01
        boundary_coefficients(2, 1) = 0.d0
        cell_coefficients(1, 1) = -1.800000000000000d+01
        cell_coefficients(2, 1) = 8.000000000000000d+00
        cell_coefficients(3, 1) = -1.800000000000000d+00
        boundary_coefficients(1, 2) = 3.200000000000000d+00
        boundary_coefficients(2, 2) = 0.d0
        cell_coefficients(1, 2) = -3.000000000000000d+00
        cell_coefficients(2, 2) = 1.000000000000000d+00
        cell_coefficients(3, 2) = -2.000000000000000d-01
    else if (orders == '1') then
        boundary_coefficients(1, 1) = -4.173913043478263d+00 * dx
        boundary_coefficients(2, 1) = 0.d0
        cell_coefficients(1, 1) = -2.347826086956524d+00
        cell_coefficients(2, 1) = 4.521739130434785d+00
        cell_coefficients(3, 1) = -1.173913043478261d+00
        boundary_coefficients(1, 2) = -1.043478260869565d+00 * dx
        boundary_coefficients(2, 2) = 0.d0
        cell_coefficients(1, 2) = 9.130434782608695d-01
        cell_coefficients(2, 2) = 1.304347826086957d-01
        cell_coefficients(3, 2) = -4.347826086956522d-02
    else if (orders == '01') then
        boundary_coefficients(1, 1) = -4.608000000000000d+01
        boundary_coefficients(2, 1) = -1.920000000000000d+01 * dx
        cell_coefficients(1, 1) = 5.400000000000000d+01
        cell_coefficients(2, 1) = -8.000000000000000d+00
        cell_coefficients(3, 1) = 1.080000000000000d+00
        boundary_coefficients(1, 2) = -1.706666666666667d+00
        boundary_coefficients(2, 2) = -1.600000000000000d+00 * dx
        cell_coefficients(1, 2) = 3.000000000000000d+00
        cell_coefficients(2, 2) = -3.333333333333333d-01
        cell_coefficients(3, 2) = 4.000000000000000d-02
    else if (orders == '02') then
        boundary_coefficients(1, 1) = -2.226086956521739d+00
        boundary_coefficients(2, 1) = 3.130434782608695d+00 * dx**2
        cell_coefficients(1, 1) = 7.043478260869565d+00
        cell_coefficients(2, 1) = -4.521739130434782d+00
        cell_coefficients(3, 1) = 7.043478260869566d-01
        boundary_coefficients(1, 2) = 1.947826086956522d+00
        boundary_coefficients(2, 2) = 2.608695652173913d-01 * dx**2
        cell_coefficients(1, 2) = -9.130434782608695d-01
        cell_coefficients(2, 2) = -4.347826086956522d-02
        cell_coefficients(3, 2) = 8.695652173913044d-03
    else if (orders == '03') then
        boundary_coefficients(1, 1) = 8.533333333333333d+00
        boundary_coefficients(2, 1) = -1.333333333333333d+00 * dx**3
        cell_coefficients(1, 1) = -1.000000000000000d+01
        cell_coefficients(2, 1) = 2.666666666666667d+00
        cell_coefficients(3, 1) = -2.000000000000000d-01
        boundary_coefficients(1, 2) = 2.844444444444445d+00
        boundary_coefficients(2, 2) = -1.111111111111111d-01 * dx**3
        cell_coefficients(1, 2) = -2.333333333333333d+00
        cell_coefficients(2, 2) = 5.555555555555556d-01
        cell_coefficients(3, 2) = -6.666666666666667d-02
    else if (orders == '12') then
        boundary_coefficients(1, 1) = 9.746192893401013d-01 * dx
        boundary_coefficients(2, 1) = 3.289340101522842d+00 * dx**2
        cell_coefficients(1, 1) = 4.659898477157360d+00
        cell_coefficients(2, 1) = -4.345177664974619d+00
        cell_coefficients(3, 1) = 6.852791878172588d-01
        boundary_coefficients(1, 2) = -8.527918781725888d-01 * dx
        boundary_coefficients(2, 2) = 1.218274111675127d-01 * dx**2
        cell_coefficients(1, 2) = 1.172588832487310d+00
        cell_coefficients(2, 2) = -1.979695431472081d-01
        cell_coefficients(3, 2) = 2.538071065989847d-02
    else if (orders == '13') then
        boundary_coefficients(1, 1) = -3.000000000000000d+00 * dx
        boundary_coefficients(2, 1) = -1.125000000000000d+00 * dx**3
        cell_coefficients(1, 1) = 0.000000000000000d+00
        cell_coefficients(2, 1) = 1.000000000000000d+00
        cell_coefficients(3, 1) = -0.000000000000000d+00
        boundary_coefficients(1, 2) = -1.000000000000000d+00 * dx
        boundary_coefficients(2, 2) = -4.166666666666666d-02 * dx**3
        cell_coefficients(1, 2) = 1.000000000000000d+00
        cell_coefficients(2, 2) = 0.000000000000000d+00
        cell_coefficients(3, 2) = -0.000000000000000d+00
    else if (orders == '23') then
        boundary_coefficients(1, 1) = 2.482758620689655d+00 * dx**2
        boundary_coefficients(2, 1) = -2.758620689655175d-01 * dx**3
        cell_coefficients(1, 1) = 3.517241379310344d+00
        cell_coefficients(2, 1) = -3.034482758620689d+00
        cell_coefficients(3, 1) = 5.172413793103448d-01
        boundary_coefficients(1, 2) = 8.275862068965518d-01 * dx**2
        boundary_coefficients(2, 2) = 2.413793103448277d-01 * dx**3
        cell_coefficients(1, 2) = 2.172413793103448d+00
        cell_coefficients(2, 2) = -1.344827586206897d+00
        cell_coefficients(3, 2) = 1.724137931034483d-01
    end if
    
end subroutine get_lower_bc_coefficients


subroutine get_upper_bc_coefficients(orders, dx, boundary_coefficients, cell_coefficients)
    character*2, intent(in) :: orders
    double precision, intent(in) :: dx
    double precision, dimension(2, 2) :: boundary_coefficients
    double precision, dimension(3, 2) :: cell_coefficients
    

    if (orders == '0') then
        boundary_coefficients(1, 1) = 3.200000000000001d+00
        boundary_coefficients(2, 1) = 0.d0
        cell_coefficients(1, 1) = -2.000000000000001d-01
        cell_coefficients(2, 1) = 1.000000000000000d+00
        cell_coefficients(3, 1) = -3.000000000000001d+00
        boundary_coefficients(1, 2) = 1.280000000000000d+01
        boundary_coefficients(2, 2) = 0.d0
        cell_coefficients(1, 2) = -1.800000000000000d+00
        cell_coefficients(2, 2) = 8.000000000000002d+00
        cell_coefficients(3, 2) = -1.800000000000000d+01
    else if (orders == '1') then
        boundary_coefficients(1, 1) = 1.043478260869565d+00 * dx
        boundary_coefficients(2, 1) = 0.d0
        cell_coefficients(1, 1) = -4.347826086956522d-02
        cell_coefficients(2, 1) = 1.304347826086956d-01
        cell_coefficients(3, 1) = 9.130434782608695d-01
        boundary_coefficients(1, 2) = 4.173913043478263d+00 * dx
        boundary_coefficients(2, 2) = 0.d0
        cell_coefficients(1, 2) = -1.173913043478261d+00
        cell_coefficients(2, 2) = 4.521739130434783d+00
        cell_coefficients(3, 2) = -2.347826086956522d+00
    else if (orders == '01') then
        boundary_coefficients(1, 1) = -1.706666666666667d+00
        boundary_coefficients(2, 1) = 1.600000000000000d+00 * dx
        cell_coefficients(1, 1) = 4.000000000000002d-02
        cell_coefficients(2, 1) = -3.333333333333335d-01
        cell_coefficients(3, 1) = 3.000000000000001d+00
        boundary_coefficients(1, 2) = -4.608000000000001d+01
        boundary_coefficients(2, 2) = 1.920000000000000d+01 * dx
        cell_coefficients(1, 2) = 1.080000000000000d+00
        cell_coefficients(2, 2) = -8.000000000000002d+00
        cell_coefficients(3, 2) = 5.400000000000001d+01
    else if (orders == '02') then
        boundary_coefficients(1, 1) = 1.947826086956522d+00
        boundary_coefficients(2, 1) = 2.608695652173913d-01 * dx**2
        cell_coefficients(1, 1) = 8.695652173913042d-03
        cell_coefficients(2, 1) = -4.347826086956520d-02
        cell_coefficients(3, 1) = -9.130434782608696d-01
        boundary_coefficients(1, 2) = -2.226086956521739d+00
        boundary_coefficients(2, 2) = 3.130434782608695d+00 * dx**2
        cell_coefficients(1, 2) = 7.043478260869566d-01
        cell_coefficients(2, 2) = -4.521739130434782d+00
        cell_coefficients(3, 2) = 7.043478260869565d+00
    else if (orders == '03') then
        boundary_coefficients(1, 1) = 2.844444444444444d+00
        boundary_coefficients(2, 1) = 1.111111111111111d-01 * dx**3
        cell_coefficients(1, 1) = -6.666666666666665d-02
        cell_coefficients(2, 1) = 5.555555555555555d-01
        cell_coefficients(3, 1) = -2.333333333333333d+00
        boundary_coefficients(1, 2) = 8.533333333333333d+00
        boundary_coefficients(2, 2) = 1.333333333333333d+00 * dx**3
        cell_coefficients(1, 2) = -2.000000000000000d-01
        cell_coefficients(2, 2) = 2.666666666666667d+00
        cell_coefficients(3, 2) = -1.000000000000000d+01
    else if (orders == '12') then
        boundary_coefficients(1, 1) = 8.527918781725889d-01 * dx
        boundary_coefficients(2, 1) = 1.218274111675127d-01 * dx**2
        cell_coefficients(1, 1) = 2.538071065989848d-02
        cell_coefficients(2, 1) = -1.979695431472081d-01
        cell_coefficients(3, 1) = 1.172588832487310d+00
        boundary_coefficients(1, 2) = -9.746192893401009d-01 * dx
        boundary_coefficients(2, 2) = 3.289340101522843d+00 * dx**2
        cell_coefficients(1, 2) = 6.852791878172589d-01
        cell_coefficients(2, 2) = -4.345177664974620d+00
        cell_coefficients(3, 2) = 4.659898477157362d+00
    else if (orders == '13') then
        boundary_coefficients(1, 1) = 1.000000000000000d+00 * dx
        boundary_coefficients(2, 1) = 4.166666666666666d-02 * dx**3
        cell_coefficients(1, 1) = 0.000000000000000d+00
        cell_coefficients(2, 1) = 0.000000000000000d+00
        cell_coefficients(3, 1) = 1.000000000000000d+00
        boundary_coefficients(1, 2) = 3.000000000000000d+00 * dx
        boundary_coefficients(2, 2) = 1.125000000000000d+00 * dx**3
        cell_coefficients(1, 2) = 2.273736754432321d-17
        cell_coefficients(2, 2) = 9.999999999999999d-01
        cell_coefficients(3, 2) = 0.000000000000000d+00
    else if (orders == '23') then
        boundary_coefficients(1, 1) = 8.275862068965518d-01 * dx**2
        boundary_coefficients(2, 1) = -2.413793103448276d-01 * dx**3
        cell_coefficients(1, 1) = 1.724137931034483d-01
        cell_coefficients(2, 1) = -1.344827586206897d+00
        cell_coefficients(3, 1) = 2.172413793103449d+00
        boundary_coefficients(1, 2) = 2.482758620689656d+00 * dx**2
        boundary_coefficients(2, 2) = 2.758620689655171d-01 * dx**3
        cell_coefficients(1, 2) = 5.172413793103449d-01
        cell_coefficients(2, 2) = -3.034482758620690d+00
        cell_coefficients(3, 2) = 3.517241379310345d+00
    end if

end subroutine get_upper_bc_coefficients
