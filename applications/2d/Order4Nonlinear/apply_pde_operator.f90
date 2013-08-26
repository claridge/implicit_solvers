subroutine apply_pde_operator(t, q, output)

! For the PDE q_t = g(q), calculates g(q).
!
! Args:
!   t: Time at which the operator is evaluated.
!   q: PDE solution at time t.
!   output: g(q), calculated here.

    !# use omp_lib
    implicit none
    
    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn
    
    double precision, intent(in) :: t
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn), intent(in) :: q
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn), intent(out) :: output

    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc) :: q_laplacian, x_flux, y_flux
    integer :: ix, iy
    double precision, external :: derivative0
    double precision :: q_face


    call get_laplacian(q(:,:,1), q_laplacian)

    do iy = 1, my
        do ix = 1, mx+1
            q_face = derivative0(q(ix-2:ix+1, iy, 1))
            x_flux(ix, iy) = q_face * (q_laplacian(ix, iy) - q_laplacian(ix-1, iy)) / dx
        end do
    end do

    do iy = 1, my+1
        do ix = 1, mx
            q_face = derivative0(q(ix, iy-2:iy+1, 1))
            y_flux(ix, iy) = q_face * (q_laplacian(ix, iy) - q_laplacian(ix, iy-1)) / dy
        end do
    end do    

    do iy = 1, my
        do ix = 1, mx
            output(ix, iy, 1) = - (x_flux(ix+1, iy) - x_flux(ix, iy)) / dx  &
                - (y_flux(ix, iy+1) - y_flux(ix, iy)) / dy
        end do
    end do
    
end subroutine apply_pde_operator
