subroutine solve_newton_system(t, dt, iterate, d_iterate)

! Solve the Netwon system using BiCGStab.  d_iterate stores RHS upon entry.
! TODO: Return success/failure.
! TODO: Make max_iter a global parameter.
! TODO: Add reference to Saad's book.

    !$ use omp_lib
    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision :: cg_tolerance
    integer :: cg_verbosity
    common /cg_config/ cg_tolerance, cg_verbosity

    double precision, intent(in) :: t, dt
    double precision, intent(in), dimension(1-mbc:mx+mbc, meqn) :: iterate
    double precision, intent(inout), dimension(1-mbc:mx+mbc, meqn) :: d_iterate

    integer :: max_iter, iter, ix, ieqn
    double precision, dimension(1-mbc:mx+mbc, meqn) :: r, r_star, p, s, Ap, As
    double precision :: alpha, omega, beta, r_dot_r_star, r_dot_r_star_old
    double precision, external :: inner_product
    double precision :: residual_norm


    max_iter = 10 * mx

    do ieqn = 1, meqn
        !$omp parallel do
        do ix = 1, mx
            r(ix, ieqn) = d_iterate(ix, ieqn)
            r_star(ix, ieqn) = r(ix, ieqn)
            p(ix, ieqn) = r(ix, ieqn)
            d_iterate(ix, ieqn) = 0.d0
        end do
    end do
    r_dot_r_star = inner_product(r, r_star)
    residual_norm = sqrt(r_dot_r_star)  ! Using r_star = r

    if (residual_norm <= cg_tolerance) then
        if (cg_verbosity > 0) then
            print '(A)', 'BiCGStab finished without iterating'
        end if
        return
    end if

    do iter = 1, max_iter
        if (cg_verbosity > 1) print '(A,I4)', '  Iteration ', iter

        call apply_newton_operator(t, dt, iterate, p, Ap)
        alpha = inner_product(r, r_star) / inner_product(Ap, r_star)

        do ieqn = 1, meqn
            !$omp parallel do
            do ix = 1, mx
                s(ix, ieqn) = r(ix, ieqn) - alpha * Ap(ix, ieqn)
            end do
        end do

        call apply_newton_operator(t, dt, iterate, s, As)
        omega = inner_product(As, s) / inner_product(As, As)
        do ieqn = 1, meqn
            !$omp parallel do
            do ix = 1, mx
                d_iterate(ix, ieqn) = d_iterate(ix, ieqn) +  &
                    alpha * p(ix, ieqn) + omega * s(ix, ieqn)
                r(ix, ieqn) = s(ix, ieqn) - omega * As(ix, ieqn)
            end do
        end do

        residual_norm = sqrt(inner_product(r, r))
        if (residual_norm <= cg_tolerance) then
            if (cg_verbosity > 0) then
                print '(A,I5,A,E16.10)', 'BiCGStab completed after ', iter,  &
                    ' iterations with norm(residual) = ', residual_norm
            end if
            return
        end if

        r_dot_r_star_old = r_dot_r_star
        r_dot_r_star = inner_product(r, r_star)
        
        beta = r_dot_r_star / r_dot_r_star_old * alpha / omega

        do ieqn = 1, meqn
            !$omp parallel do
            do ix = 1, mx
                p(ix, ieqn) = r(ix, ieqn) + beta * (p(ix, ieqn) -  &
                    omega * Ap(ix, ieqn))
            end do
        end do

        ! TODO: Merge with logging of iteration number
        if (cg_verbosity > 1) print '(A,E16.10)', "  residual_norm = ",  &
            residual_norm
    end do

    print *, 'Error: BiCGStab failed to converge'
    stop

end subroutine solve_newton_system
