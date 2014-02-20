subroutine solve_backward_euler_system(t, dt, iterate, perturbation, residual_norm)

! Solve the Netwon system,
!     newton_operator(perturbation) = newton_rhs,
! using BiCGStab.
!
! The implementation is taken from Algorithm 7.7 of Y. Saad, Iterative Methods
! for Sparse Linear Systems, 2nd Ed.
!
! The algorithm will continue until either norm(residual) < linear_solver_tolerance,
! or max_iter iterations have been performed.
!
! Args:
!   t: Time at beginning of the current step.
!   dt: Length of the current time step.
!   iterate: The current Newton iterate.
!   perturbation: iterate + perturbation will yield the next Newton iterate.
!     Stores newton_rhs on input.
!   residual_norm: Returns final norm of the residual

! TODO: Make max_iter a global parameter?

    !$ use omp_lib
    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision :: linear_solver_tolerance
    integer :: linear_solver_verbosity
    common /cg_config/ linear_solver_tolerance, linear_solver_verbosity

    double precision, intent(in) :: t, dt
    double precision, intent(in), dimension(1-mbc:mx+mbc, meqn) :: iterate
    double precision, intent(inout), dimension(1-mbc:mx+mbc, meqn) :: perturbation
    double precision, intent(out) :: residual_norm

    integer :: max_iter, iter, ix, ieqn
    double precision, dimension(1-mbc:mx+mbc, meqn) :: r, r_star, p, s, Ap, As
    double precision :: alpha, omega, beta, r_dot_r_star, r_dot_r_star_old
    double precision, external :: inner_product
    double precision :: denominator

    max_iter = 10 * mx

    do ieqn = 1, meqn
        !$omp parallel do
        do ix = 1, mx
            r(ix, ieqn) = perturbation(ix, ieqn)
            r_star(ix, ieqn) = r(ix, ieqn)
            p(ix, ieqn) = r(ix, ieqn)
            perturbation(ix, ieqn) = 0.d0
        end do
    end do
    r_dot_r_star = inner_product(r, r_star)
    residual_norm = sqrt(r_dot_r_star)  ! Using r_star = r

    if (residual_norm <= linear_solver_tolerance) then
        if (linear_solver_verbosity > 0) then
            print '(A)', 'BiCGStab finished without iterating'
        end if
        return
    end if

    do iter = 1, max_iter
        call get_backward_euler_lhs(t, dt, iterate, p, Ap)
        denominator = inner_product(Ap, r_star)
        alpha = inner_product(r, r_star) / denominator
        
        ! TODO: Not sure if this is reasonable.
        if (denominator == 0.d0 .or. alpha == 0.d0) then
            if (linear_solver_verbosity > 0) then
                print '(A,A,E16.10,A)', 'BiCGStab reached degenerate condition, ',  &
                    'but reporting success anyway with residual_norm = ', residual_norm, '.'
            end if
            return
        end if

        do ieqn = 1, meqn
            !$omp parallel do
            do ix = 1, mx
                s(ix, ieqn) = r(ix, ieqn) - alpha * Ap(ix, ieqn)
            end do
        end do

        call get_backward_euler_lhs(t, dt, iterate, s, As)

        denominator = inner_product(As, As)        
        omega = inner_product(As, s) / denominator
        if (denominator == 0.d0 .or. omega == 0.d0) then
            if (linear_solver_verbosity > 0) then
                print '(A,A,E16.10,A)', 'BiCGStab reached degenerate condition, ',  &
                    'but reporting success anyway with residual_norm = ', residual_norm, '.'
            end if
            return
        end if        

        do ieqn = 1, meqn
            !$omp parallel do
            do ix = 1, mx
                perturbation(ix, ieqn) = perturbation(ix, ieqn) +  &
                    alpha * p(ix, ieqn) + omega * s(ix, ieqn)
                r(ix, ieqn) = s(ix, ieqn) - omega * As(ix, ieqn)
            end do
        end do

        residual_norm = sqrt(inner_product(r, r))
        if (residual_norm <= linear_solver_tolerance) then
            if (linear_solver_verbosity > 0) then
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

        if (linear_solver_verbosity > 1) then
            print '(A,I4,A,E16.10)', 'Iteration ', iter, ': residual_norm = ',  &
                residual_norm
        end if
    end do
    
    if (linear_solver_verbosity > 0) then
        print '(A,I5,A,E16.10)', 'BiCGStab gave up after ', max_iter,  &
            ' iterations with norm(residual) = ', residual_norm
    end if

end subroutine solve_backward_euler_system
