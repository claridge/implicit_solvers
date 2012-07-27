subroutine solve_newton_system(t, dt, iterate, d_iterate)

! Solve the Netwon system using BiCGStab.  d_iterate stores RHS upon entry.

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision, intent(in) :: t, dt
    double precision, intent(in), dimension(1-mbc:mx+mbc, meqn) :: iterate
    double precision, intent(inout), dimension(1-mbc:mx+mbc, meqn) :: d_iterate

    integer :: max_iter, iter
    double precision, dimension(1-mbc:mx+mbc, meqn) :: r, r_star, p, s, Ap, As
    double precision :: alpha, omega, beta, r_dot_r_star, r_dot_r_star_old
    double precision, external :: inner_product
    double precision :: residual_norm

    double precision :: cg_tolerance
    integer :: cg_verbosity
    common /cg_config/ cg_tolerance, cg_verbosity


    max_iter = 10 * mx

    r = d_iterate
    r_star = r
    p = r
    d_iterate = 0.d0
    r_dot_r_star = inner_product(r, r_star)
    residual_norm = sqrt(inner_product(r, r))

    if (residual_norm <= cg_tolerance) then
        if (cg_verbosity > 0) then
            print '(A)', "BiCGStab finished without iterating"
        end if
        return
    end if

    do iter = 1, max_iter
        if (cg_verbosity > 1) print '(A,I4)', "  Iteration ", iter

        call apply_newton_operator(t, dt, iterate, p, Ap)
        alpha = inner_product(r, r_star) / inner_product(Ap, r_star)
        s = r - alpha * Ap
        call apply_newton_operator(t, dt, iterate, s, As)
        omega = inner_product(As, s) / inner_product(As, As)
        d_iterate = d_iterate + alpha * p + omega * s

        r_dot_r_star_old = r_dot_r_star
        r = s - omega * As
        residual_norm = sqrt(inner_product(r, r))
        if (residual_norm <= cg_tolerance) then
            if (cg_verbosity > 0) then
                print '(A,I5,A,E16.10)', 'BiCGStab completed after ', iter,  &
                    ' iterations with norm(residual) = ', residual_norm
            end if
            return
        end if

        r_dot_r_star = inner_product(r, r_star)
        
        beta = r_dot_r_star / r_dot_r_star_old * alpha / omega
        p = r + beta * (p - omega * Ap)

        ! TODO: Merge with logging of iteration number
        if (cg_verbosity > 1) print '(A,E16.10)', "  residual_norm = ",  &
            residual_norm
    end do

    print *, 'Error: BiCGStab failed to converge'
    stop

end subroutine solve_newton_system
