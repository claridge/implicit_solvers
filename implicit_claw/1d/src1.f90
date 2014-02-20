subroutine src1(maxmx, meqn, mbc, mx, x_lower, dx, q, maux, aux, t, dt)

    implicit none

    integer, intent(in) :: maxmx, meqn, mbc, mx, maux
    double precision, intent(in) :: x_lower, dx, t, dt
    double precision, intent(in) :: aux(1-mbc:maxmx+mbc, maux)
    double precision, intent(inout), dimension(1-mbc:maxmx+mbc, meqn) :: q

    integer :: n_splits, i
    double precision :: t_new, t_local, dt_local
    logical :: success

    character :: implicit_integration_scheme
    common /implicit_config/ implicit_integration_scheme

    ! Not currently exposed.  If >0. we'll taking two halved time steps should
    ! Newton's method fail to converge.
    integer :: max_time_step_splits=0
        
    
    ! Clawpack currently gives src1 the time at the *end* of the step.
    t_new = t
    t_local  = t_new-dt
    dt_local = dt
    n_splits = 0
    success  = .false.

    do while (t_new - t_local > 1d-8)
        if (implicit_integration_scheme .eq. 'b' .or.  &
            implicit_integration_scheme .eq. 'B') then
            call take_backward_euler_step(t_local, dt_local, q, success)
        else if (implicit_integration_scheme .eq. 'c' .or.  &
                 implicit_integration_scheme .eq. 'C') then
            call take_crank_nicolson_step(t_local, dt_local, q, success)
        else if (implicit_integration_scheme .eq. 'f' .or.  &
                 implicit_integration_scheme .eq. 'F') then
            call take_forward_euler_step(t_local, dt_local, q, success)
        else
            print *, "src1: Error: Invalid implicit integration scheme specified."
            stop
        end if
 
        if (success) then
            t_local  = t_local + dt_local
        else if (max_time_step_splits > 0) then
            n_splits = n_splits + 1
            if (n_splits > max_time_step_splits) then
                print *, "src1: Implicit step failed with maximum number of ",  &
                    "dt splits; aborting."
                stop
            end if
            dt_local = dt_local / 2.d0
            print *, "src1: Implicit step failed; retrying with dt = ", dt_local
        else
            print *, "src1: Implicit step failed; aborting."
            stop
        end if
    end do
    
end subroutine src1
