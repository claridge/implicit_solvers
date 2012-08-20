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

    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc) :: q_laplacian
    integer :: ix, iy

    q_laplacian = 0.d0
    call calculate_laplacian(q(:,:,1), q_laplacian)

! print *, 'q(-1:3, -1:3)'
! print *, q(-1:3, 3, 1)
! print *, q(-1:3, 2, 1)
! print *, q(-1:3, 1, 1)
! print *, q(-1:3, 0, 1)
! print *, q(-1:3, -1, 1)
! print *, ''
! 
! print *, 'q_laplacian(-1:3, -1:3)'
! print *, q_laplacian(-1:3, 3)
! print *, q_laplacian(-1:3, 2)
! print *, q_laplacian(-1:3, 1)
! print *, q_laplacian(-1:3, 0)
! print *, q_laplacian(-1:3, -1)
! print *, ''

! print *, x_flux(0, 1), x_flux(1, 1)
! print *, y_flux(0, 1), y_flux(1, 1)
! print *, ''

    !#omp parallel do private(ix)
    do iy = 1, my
        do ix = 1, mx
            output(ix, iy, 1) = - (x_flux(ix+1, iy) - x_flux(ix, iy)) / dx  &
                - (y_flux(ix, iy+1) - y_flux(ix, iy)) / dy
        end do
    end do
    
!     print *, q(-1:3, 3, 1)
!     print *, q(-1:3, 2, 1)
!     print *, q(-1:3, 1, 1)
!     print *, q(-1:3, 0, 1)
!     print *, q(-1:3, -1, 1)
!     print *, ''
!     print *, output(-1:3, 1, 1)
!     print *, ''
!     print *, ''
    
    contains
    
    double precision function x_flux(ix, iy)
        integer :: ix, iy
        double precision :: q_face
        q_face = (q(ix-1, iy, 1) + q(ix, iy, 1)) / 2.d0
        x_flux = q_face * (q_laplacian(ix, iy) - q_laplacian(ix-1, iy)) / dx
    end function x_flux

    double precision function y_flux(ix, iy)
        integer :: ix, iy
        double precision :: q_face
        q_face = (q(ix, iy-1, 1) + q(ix, iy, 1)) / 2.d0
        y_flux = q_face * (q_laplacian(ix, iy) - q_laplacian(ix, iy-1)) / dy
!         if (abs(y_flux) > 1d-12) print *, ix, iy, q(ix, iy-2:iy+1, 1)
    end function y_flux
    
end subroutine apply_pde_operator


subroutine apply_linearized_pde_operator(t, q, p, output)

! For the PDE q_t = g(q), calculates g'[q](p), i.e. the linearization
! of g, about q, applied to perturbation p
!
! Args:
!   t: Time at which the operator is evaluated.
!   q: Base function of the linearization.
!   p: Perturbation to which the linearized operator is applied.
!   output: g'[q](p), calculated here.

    !# use omp_lib
    implicit none

    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn
    
    double precision, intent(in) :: t
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn), intent(in) :: q, p
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn), intent(out) :: output
    
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc) :: q_laplacian, p_laplacian
    integer :: ix, iy

    
    call calculate_laplacian(q(:,:,1), q_laplacian)
    call calculate_laplacian(p(:,:,1), p_laplacian)
    
    !#omp parallel do private(ix)
    do iy = 1, my
        do ix = 1, mx
            output(ix, iy, 1) = - (x_flux(ix+1, iy) - x_flux(ix, iy)) / dx  &
                - (y_flux(ix, iy+1) - y_flux(ix, iy)) / dy
        end do
    end do
    
    contains
    
    double precision function x_flux(ix, iy)
        integer :: ix, iy
        double precision :: q_face, p_face
        q_face = (q(ix-1, iy, 1) + q(ix, iy, 1)) / 2.d0
        p_face = (p(ix-1, iy, 1) + p(ix, iy, 1)) / 2.d0
        x_flux = (q_laplacian(ix, iy) - q_laplacian(ix-1, iy)) / dx * p_face  &
            + q_face * (p_laplacian(ix, iy) - p_laplacian(ix-1, iy)) / dx
    end function x_flux

    double precision function y_flux(ix, iy)
        integer :: ix, iy
        double precision :: q_face, p_face
        q_face = (q(ix, iy-1, 1) + q(ix, iy, 1)) / 2.d0
        p_face = (p(ix, iy-1, 1) + p(ix, iy, 1)) / 2.d0
        y_flux = (q_laplacian(ix, iy) - q_laplacian(ix, iy-1)) / dy * p_face  &
            + q_face * (p_laplacian(ix, iy) - p_laplacian(ix, iy-1)) / dy
    end function y_flux
    
    
end subroutine apply_linearized_pde_operator


subroutine calculate_laplacian(q, q_laplacian)

! Calculate Laplacian of q at cell centers.  Fills ghost cells within the
! outermost layer.

    !# use omp_lib
    implicit none

    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn

    double precision, dimension(1-mbc:mx+mbc,1-mbc:my+mbc), intent(in) :: q
    double precision, dimension(1-mbc:mx+mbc,1-mbc:my+mbc), intent(out) :: q_laplacian

    integer :: ix, iy

    !#omp parallel do private(ix)
    do iy = 2-mbc,my+mbc-1
        do ix = 2-mbc,mx+mbc-1
            q_laplacian(ix,iy) = (q(ix-1,iy) - 2.d0*q(ix,iy) + q(ix+1,iy)) / dx**2
            q_laplacian(ix,iy) = q_laplacian(ix,iy) + (q(ix,iy-1) - 2.d0*q(ix,iy) +  &
                                 q(ix,iy+1)) / dy**2
        end do
    end do
    
end subroutine calculate_laplacian