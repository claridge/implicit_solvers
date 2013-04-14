subroutine get_backward_euler_lhs(t, dt, r, p, output)

!----------------------------------------------------------------
! Calculates
!     m'[r](p) = p + dt * g'[r](p),
! as appears on the left-hand side of Newton's method.
!
! Args:
!   t: Most recent time value.
!   dt: Time step size.
!   r: Current Newton iterate.
!   p: Current Newton perturbation.
!   output: Result of the operation.
!----------------------------------------------------------------

    !$ use omp_lib
    implicit none
    
    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn
    
    double precision :: beta, kappa, delta, mu, right_film_height
    common /surfactant_config/ beta, kappa, delta, mu, right_film_height    

    double precision, intent(in) :: t, dt
    double precision, intent(in), dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn) :: r, p    
    double precision, intent(out), dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn) :: output

    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc) :: rh_laplacian, output_tmp
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, 2) :: p_tmp

    integer :: ix, iy, ieqn, ix_max, iy_max
    double precision :: max_diff

    call apply_linearized_bcs(r, p)
    call apply_linearized_pde_operator(t, r, p, output)    

    do ieqn = 2, 2
        !$omp parallel do private(ix)
        do iy = 1, my
            do ix = 1,mx
                output(ix, iy, ieqn) = p(ix, iy, ieqn) - dt * output(ix, iy, ieqn)
            end do
        end do
    end do
    
    call calculate_laplacian(r(:, :, 1), rh_laplacian) 
    call apply_adi_y_operator(r(:, :, 1), p(:, :, 1), p_tmp(:, :, 1))

    call apply_linearized_bcs(r, p_tmp)
    call apply_adi_x_operator(r(:, :, 1), p_tmp(:, :, 1), output(:, :, 1))


!     max_diff = 0.d0
!     do iy = 1, my
!         do ix = 1, mx
!             if (abs(output_tmp(ix, iy) - output(ix, iy, 1)) > max_diff) then
!                 max_diff = abs(output_tmp(ix, iy) - output(ix, iy, 1))
!                 ix_max = ix
!                 iy_max = iy
!             end if
!         end do
!     end do
!     print *, ix_max, iy_max, output(ix_max, iy_max, 1), output_tmp(ix_max, iy_max)
! 
!     output(:, :, 1) = output_tmp

    contains
        

    subroutine apply_adi_x_operator(rh, ph, op_ph)
        implicit none
        double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc), intent(in) :: rh, ph
        double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc), intent(out) :: op_ph
        double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc) :: flux
        double precision :: rh_face, rh_x, rh_laplacian_x, ph_face, ph_x, ph_xxx
    
        !#omp parallel do private(ix, rh_face, rh_x, rh_laplacian_x, ph_face, ph_x, ph_xxx)
        do iy = 1, my
            do ix = 1, mx+1
                rh_face = (rh(ix-1, iy) + rh(ix, iy)) / 2.d0
                rh_x = (rh(ix, iy) - rh(ix-1, iy)) / dx
                rh_laplacian_x = (rh_laplacian(ix, iy) - rh_laplacian(ix-1, iy)) / dx

                ph_face = (ph(ix-1, iy) + ph(ix, iy)) / 2.d0
                ph_x = (ph(ix, iy) - ph(ix-1, iy)) / dx
                ph_xxx = (-ph(ix-2, iy) + 3.d0 * ph(ix-1, iy)  &
                           - 3.d0 * ph(ix, iy) + ph(ix+1, iy)) / dx**3
                flux(ix, iy) = rh_face**3 / 3.d0 * (-beta * ph_x + kappa * ph_xxx)  &
                    + rh_face**2 * (-beta * rh_x + kappa * rh_laplacian_x) * ph_face
            end do
        end do
        
        !$omp parallel do private(ix)
        do iy = 1, my
            do ix = 1, mx
                op_ph(ix, iy) = ph(ix, iy) + dt * (flux(ix+1, iy) - flux(ix, iy)) / dx
            end do
        end do
    end subroutine apply_adi_x_operator
    
    
    subroutine apply_adi_y_operator(rh, ph, op_ph)
        implicit none
        double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc), intent(in) :: rh, ph
        double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc), intent(out) :: op_ph
        double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc) :: flux
        double precision :: rh_face, rh_y, rh_laplacian_y, ph_face, ph_y, ph_yyy
    
        !$omp parallel do private(ix, rh_face, rh_y, rh_laplacian_y, ph_face, ph_y, ph_yyy)
        do iy = 1, my+1
            do ix = 1, mx
                rh_face = (rh(ix, iy-1) + rh(ix, iy)) / 2.d0
                rh_y = (rh(ix, iy) - rh(ix, iy-1)) / dy
                rh_laplacian_y = (rh_laplacian(ix, iy) - rh_laplacian(ix, iy-1)) / dy

                ph_face = (ph(ix, iy-1) + ph(ix, iy)) / 2.d0
                ph_y = (ph(ix, iy) - ph(ix, iy-1)) / dy
                ph_yyy = (-ph(ix, iy-2) + 3.d0 * ph(ix, iy-1)  &
                           - 3.d0 * ph(ix, iy) + ph(ix, iy+1)) / dy**3
                flux(ix, iy) = rh_face**3 / 3.d0 * (-beta * ph_y + kappa * ph_yyy)  &
                    + rh_face**2 * (-beta * rh_y + kappa * rh_laplacian_y) * ph_face
            end do
        end do
        
        !$omp parallel do private(ix)
        do iy = 1, my
            do ix = 1, mx
                op_ph(ix, iy) = ph(ix, iy) + dt * (flux(ix, iy+1) - flux(ix, iy)) / dy
            end do
        end do
    end subroutine apply_adi_y_operator



    
end subroutine get_backward_euler_lhs
