subroutine qinit(maxmx, maxmy, meqn, mbc, mx, my, x_low, y_low, dx, dy, q)

! Sets initial conditions for:
!     film height = q(:,:,1);
!     surfactant concentration = q(:,:,2).

    implicit none

    integer, intent(in) :: maxmx, maxmy, meqn, mbc, mx, my
    double precision, intent(in) :: x_low, y_low, dx, dy
    double precision, intent(out) :: q(1-mbc:maxmx+mbc, 1-mbc:maxmy+mbc, meqn)

    integer :: i, j, imode, k(4)
    double precision :: x, y, a, b, big_b, c(4), pi=dacos(-1.d0)
    double precision :: heaviside_1mx

    a = 1.d-2
    b = 5.d-2
    big_b = 5.d0

    c(1) = 1.d0
    c(2) = 1.d0
    c(3) = 1.d0
    c(4) = 0.5d0

    k(1) = 2
    k(2) = 5
    k(3) = 7
    k(4) = 20

    do j = 1,my
        y = y_low + (j-0.5d0)*dy
        do i = 1,mx
            x = x_low + (i-0.5d0)*dx

            heaviside_1mx = 0.5d0*(1+tanh(20*(1-x)))

            q(i,j,1) = (1-x**2+b)*heaviside_1mx + b*(1.d0-heaviside_1mx)
           
            do imode = 1,4
                q(i,j,1) = q(i,j,1) + a * exp(-big_b*(x-1)**2) * c(imode) * cos(k(imode)*y)
            end do

           q(i,j,2) = heaviside_1mx

        end do
    end do
  
end subroutine qinit
