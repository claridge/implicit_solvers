double precision function inner_product(u, v)

! Calculates the inner product of cell-centered vectors u and v, using
! interior cells only.

    !$ use omp_lib
    implicit none
    
    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision, dimension(1-mbc:mx+mbc, meqn) :: u, v

    integer :: ix, ieqn

    inner_product = 0.d0
    do ieqn = 1, meqn
        !$omp parallel do reduction(+ : inner_product)
        do ix = 1, mx
            inner_product = inner_product + u(ix, ieqn) * v(ix, ieqn)
        end do
    end do

end function inner_product
