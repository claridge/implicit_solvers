double precision function inner_product(u, v)

! Calculates the inner product of cell-centered vectors u and v, using
! interior cells only.

    !$ use omp_lib
    implicit none
    
    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn

    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn) :: u, v

    integer :: ix, iy, ieqn


    inner_product = 0.d0
    do ieqn = 1, meqn
        !$omp parallel do reduction(+ : inner_product) private(ix)
        do iy = 1, my
            do ix = 1, mx
                inner_product = inner_product + u(ix, iy, ieqn) * v(ix, iy, ieqn)
            end do
        end do
    end do

end function inner_product
