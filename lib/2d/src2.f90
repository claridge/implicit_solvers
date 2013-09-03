subroutine src2(maxmx, maxmy, meqn, mbc, mx, my, x_lower, y_lower, dx, dy, q,  &
                maux, aux, t, dt)

    implicit none

    integer, intent(in) :: maxmx, maxmy, meqn, mbc, mx, my, maux
    double precision, intent(in) :: x_lower, y_lower, dx, dy, t, dt
    double precision, intent(in) :: aux(1-mbc:maxmx+mbc, maux)
    double precision, dimension(1-mbc:mx+mbc, 1-mbc:my+mbc, meqn),  &
        intent(inout) :: q

    integer :: n_splits, i
    double precision :: t_new, t_local, dt_local
    logical :: success

    character :: implicit_integration_scheme
    common /implicit_config/ implicit_integration_scheme

    ! Not currently exposed.  If >0. we"ll taking two halved time steps should
    ! Newton"s method fail to converge.
    integer :: max_time_step_splits=0
        
    ! src2 receives the time value at the *end* of the step, unlike src1.
    t_new = t
    t_local = t_new - dt
    dt_local = dt
    n_splits = 0
    success  = .false.

    do while (t_new - t_local > 1d-8)
        if (implicit_integration_scheme .eq. "b" .or.  &
            implicit_integration_scheme .eq. "B") then
            call take_backward_euler_step(t_local, dt_local, q, success)
        else if (implicit_integration_scheme .eq. "c" .or.  &
                 implicit_integration_scheme .eq. "C") then
            call take_crank_nicolson_step(t_local, dt_local, q, success)
        else if (implicit_integration_scheme .eq. "f" .or.  &
                 implicit_integration_scheme .eq. "F") then
            call take_forward_euler_step(t_local, dt_local, q, success)
        else
            print *, "src2: Error: Invalid implicit integration scheme specified."
            stop
        end if
 
        if (success) then
            t_local  = t_local + dt_local
        else if (max_time_step_splits > 0) then
            n_splits = n_splits + 1
            if (n_splits > max_time_step_splits) then
                print *, "src2: Implicit step failed with maximum number of ",  &
                    "dt splits; aborting."
                stop
            end if
            dt_local = dt_local / 2.d0
            print *, "src2: Implicit step failed; retrying with dt = ", dt_local
        else
            print *, "src2: Implicit step failed; aborting."
            stop
        end if    
    end do
    
end subroutine src2
