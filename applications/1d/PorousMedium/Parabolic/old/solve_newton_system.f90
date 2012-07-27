subroutine solve_newton_system(t, dt, iterate, d_iterate)

!----------------------------------------------------------
! Solves the system
!     m'[iterate](d_iterate) = -m(iterate)
! that arises in each iteration of Newton's method for
! the porous medium equation.
!
! Upon input, d_iterate stores the right-hand side.
!----------------------------------------------------------

    implicit none
    
    integer :: mx, mbc
    double precision :: x_lower, dx
    common /grid_data/ mx, mbc, x_lower, dx
    
    double precision, intent(in) :: t, dt
    double precision, intent(in), dimension(1-mbc:mx+mbc) :: iterate
    double precision, intent(inout), dimension(1-mbc:mx+mbc) :: d_iterate    

    double precision, dimension(1-mbc:mx+mbc) ::  &
        residual, search_direction, search_direction_image
    double precision :: descent_length, residual_ratio, residual_norm,  &
        old_residual_norm
    integer :: ix, iter, max_iter
    logical :: verbose=.true.

    double precision :: cg_tolerance
    common /cg_config/ cg_tolerance

    double precision, external :: inner_product

    
    print '(A)', "Beginning density_solve_newton_system"
    
    ! Initialization
    do ix = 1,mx
        residual(ix) = d_iterate(ix)
        search_direction(ix) = residual(ix)
        d_iterate(ix) = 0.d0
    end do

    residual_norm = sqrt(inner_product(residual, residual))
    print '(A,E16.10)', "  Initial residual_norm = ", residual_norm
    
    if (residual_norm <= cg_tolerance) then
        print '(A)', "  Finished without iterating"
        return
    end if

    old_residual_norm = residual_norm
    
    iter = 0
    max_iter = mx
    
    do iter = 1, max_iter
        
        if (verbose) print '(A,I4)', "  Iteration ", iter
        
        call apply_newton_operator(t, dt, iterate, search_direction,  &
                                   search_direction_image)
        
        descent_length = residual_norm**2  &
                         / inner_product(search_direction, search_direction_image)
        
        do ix = 1,mx
            d_iterate(ix) = d_iterate(ix) + descent_length * search_direction(ix)
            residual(ix) = residual(ix) - descent_length * search_direction_image(ix)
        end do
        
        old_residual_norm = residual_norm
        residual_norm = sqrt(inner_product(residual, residual))
                
        if (residual_norm <= cg_tolerance) then
            print '(A,E16.10)', "  Finished with residual_norm = ", residual_norm 
            return
        end if
            
        residual_ratio = residual_norm**2 / old_residual_norm**2
        
        do ix = 1,mx
            search_direction(ix) = residual(ix) + residual_ratio * search_direction(ix)
        end do
        
        if (verbose) print '(A,E16.10)', "  residual_norm = ", residual_norm
        
    end do
    
    if (residual_norm > cg_tolerance) then
        print *, "Error: density_solve_newton_system"
        print *, "Conjugate gradient method failed to converge."
        stop
    end if
    
end subroutine solve_newton_system
