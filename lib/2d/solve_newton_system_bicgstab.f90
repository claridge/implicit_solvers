subroutine solve_newton_system(t, dt, iterate, d_iterate, success)

! Solve the Netwon system,
!     newton_operator(d_iterate) = newton_rhs,
! using BiCGStab.
!
! The implementation is taken from Algorithm 7.7 of Y. Saad, Iterative Methods
! for Sparse Linear Systems, 2nd Ed.
!
! Args:
!   t: Time at beginning of the current step.
!   dt: Length of the current time step.
!   iterate: The current Newton iterate.
!   d_iterate: iterate + d_iterate will yield the next Newton iterate.
!     Stores newton_rhs on input.
!   success: Returns .true. if BiCGStab converges, .false. if not.

! TODO: Make max_iter a global parameter?

    !$ use omp_lib
    implicit none

    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn

    double precision :: cg_tolerance
    integer :: cg_verbosity
    common /cg_config/ cg_tolerance, cg_verbosity

    double precision, intent(in) :: t, dt
    double precision, intent(in), dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn) :: iterate
    double precision, intent(inout), dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn) :: d_iterate
    logical, intent(out) :: success

    integer :: max_iter, iter, ix, iy, ieqn
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn) :: r, r_star, p, s, Ap, As
    double precision :: alpha, omega, beta, r_dot_r_star, r_dot_r_star_old
    double precision, external :: inner_product
    double precision :: residual_norm

    success = .false.
    max_iter = 10 * mx

    do ieqn = 1, meqn
        !$omp parallel do private(ix)
        do iy = 1, my
            do ix = 1, mx
                r(ix, iy, ieqn) = d_iterate(ix, iy, ieqn)
                r_star(ix, iy, ieqn) = r(ix, iy, ieqn)
                p(ix, iy, ieqn) = r(ix, iy, ieqn)
                d_iterate(ix, iy, ieqn) = 0.d0
            end do
        end do
    end do
    r_dot_r_star = inner_product(r, r_star)
    residual_norm = sqrt(r_dot_r_star)  ! Using r_star = r

    if (residual_norm <= cg_tolerance) then
        success = .true.
        if (cg_verbosity > 0) then
            print '(AE16.10)', 'BiCGStab finished without iterating; residual_norm = ', residual_norm
        end if
        return
    end if

    do iter = 1, max_iter
        call apply_newton_operator(t, dt, iterate, p, Ap)
        alpha = inner_product(r, r_star) / inner_product(Ap, r_star)

        do ieqn = 1, meqn
            !$omp parallel do private(ix)
            do iy = 1, my
                do ix = 1, mx
                    s(ix, iy, ieqn) = r(ix, iy, ieqn) - alpha * Ap(ix, iy, ieqn)
                end do
            end do
        end do

        call apply_newton_operator(t, dt, iterate, s, As)
        omega = inner_product(As, s) / inner_product(As, As)
        do ieqn = 1, meqn
            !$omp parallel do private(ix)
            do iy = 1, my
                do ix = 1, mx
                    d_iterate(ix, iy, ieqn) = d_iterate(ix, iy, ieqn) +  &
                        alpha * p(ix, iy, ieqn) + omega * s(ix, iy, ieqn)
                    r(ix, iy, ieqn) = s(ix, iy, ieqn) - omega * As(ix, iy, ieqn)
                end do
            end do
        end do

        residual_norm = sqrt(inner_product(r, r))
        if (residual_norm <= cg_tolerance) then
            success = .true.
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
            !$omp parallel do private(ix)
            do iy = 1, my
                do ix = 1, mx
                    p(ix, iy, ieqn) = r(ix, iy, ieqn) + beta * (p(ix, iy, ieqn) -  &
                        omega * Ap(ix, iy, ieqn))
                end do
            end do
        end do

        if (cg_verbosity > 1) then
            print '(A,I4,A,E16.10)', 'Iteration ', iter, ': residual_norm = ',  &
                residual_norm
        end if
    end do

end subroutine solve_newton_system
