subroutine take_backward_euler_step(t, dt, q, success)

! Takes a time step using backward Euler.
!
! Args:
!   t: Time at the beginning of the current step.
!   dt: Length of the current time step.
!   q: Discrete PDE solution being updated.  Holds value at beginning of step
!     on entry.  Contains value at end of step on return.
!   success: Indicates whether time step was successful.  (Newton's method
!     or the linear solver could have failed to converge.)

    !$ use omp_lib
    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    double precision, intent(in) :: t, dt
    double precision, dimension(1-mbc:mx+mbc, meqn), intent(inout) :: q
    logical, intent(out) :: success

    integer :: newton_max_iter, newton_verbosity
    double precision :: newton_reduction_factor, newton_tolerance
    common /newton_config/ newton_max_iter, newton_reduction_factor,  &
         newton_tolerance, newton_verbosity

    double precision, dimension(1-mbc:mx+mbc, meqn) :: d_iterate, iterate
    double precision :: norm_d_iterate, old_norm_d_iterate, residual_norm
    integer :: iter, ix, ieqn
    double precision, external :: inner_product


    success = .false.
    norm_d_iterate = newton_tolerance + 1.d0
    old_norm_d_iterate = norm_d_iterate
    iter = 0

    do ieqn = 1, meqn
        !$omp parallel do
        do ix = 1, mx
            iterate(ix, ieqn) = q(ix, ieqn)
        end do
    end do

    do while (norm_d_iterate > newton_tolerance)

        iter = iter + 1

        call get_backward_euler_rhs(t, dt, q, iterate, d_iterate)
        ! residual_norm is currently unused.
        call solve_backward_euler_system(t, dt, iterate, d_iterate, residual_norm)

        old_norm_d_iterate = norm_d_iterate
        norm_d_iterate = 0.d0

        do ieqn = 1, meqn
            !$omp parallel do
            do ix = 1, mx
                iterate(ix, ieqn) = iterate(ix, ieqn) + d_iterate(ix, ieqn)
            end do
        end do
        norm_d_iterate = sqrt(inner_product(d_iterate, d_iterate))

        if (newton_verbosity > 1) then
           print '(A,I2,A,E16.10)', 'Newton iteration ', iter,  &
                ', norm(d_iterate) = ', norm_d_iterate
        end if

        ! Give up if norm_d_iterate gets too big, as this indicates nonconvergence
        if (norm_d_iterate > 1e6) exit

        ! If we've converged, then finish.
        if (norm_d_iterate < newton_tolerance) then
            do ieqn = 1, meqn
                !$omp parallel do
                do ix = 1,mx
                    q(ix, ieqn) = iterate(ix, ieqn)
                end do
            end do

            success = .true.
            if (newton_verbosity > 0) then
               print '(A,I2,A,E16.10)', 'Newton''s method finished after ',  &
                    iter, ' iterations with norm(d_iterate) = ', norm_d_iterate
            end if
            exit
        end if

        ! If we've hit the max # of iterations, give up unless we're making
        ! sufficient progress.
        if ((iter >= newton_max_iter) .and.  &
            (norm_d_iterate > newton_reduction_factor *  &
             old_norm_d_iterate)) then
            exit
        end if

    end do

end subroutine take_backward_euler_step
