subroutine src1(maxmx, meqn, mbc, mx, x_lower, dx, q, maux, aux, t, dt)

    implicit none

    integer, intent(in) :: maxmx, meqn, mbc, mx, maux
    double precision, intent(in) :: x_lower, dx, t, dt
    double precision, intent(in) :: aux(1-mbc:maxmx+mbc, maux)    
    double precision, intent(inout), dimension(1-mbc:maxmx+mbc, meqn) :: q

    integer :: n_splits, max_splits, i
    double precision :: t_new, t_local, dt_local
    logical :: success
        
    max_splits = 0
    
    ! Clawpack currently gives src1 the time at the *end* of the step.
    t_new = t

    t_local  = t_new-dt
    dt_local = dt
    n_splits = 0
    success  = .false.


    do while (t_new - t_local > 1d-8)
        call take_backward_euler_step(t_local, dt_local, q, success)

        if (success) then
            t_local  = t_local + dt_local
        else
            dt_local = dt_local / 2.d0
            n_splits = n_splits + 1
            print *, 'Splitting time step'

            if (n_splits > max_splits) then
                print '(A)', 'Parabolic step failed with maximum number of dt splits; aborting.'
                stop
            end if
        end if
    
    end do
    
end subroutine src1
