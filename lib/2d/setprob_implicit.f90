subroutine setprob_implicit(file_id)

    !$ use omp_lib
    implicit none
    
    integer :: mx, my, mbc, meqn
    double precision :: x_lower, y_lower, dx, dy
    common /claw_config/ mx, my, mbc, x_lower, y_lower, dx, dy, meqn

    character :: implicit_integration_scheme
    integer :: max_time_step_splits
    common /implicit_config/ implicit_integration_scheme, max_time_step_splits
    integer, intent(in) :: file_id
    
    integer :: newton_max_iter, newton_verbosity
    double precision :: newton_reduction_factor, newton_tolerance
    common /newton_config/ newton_max_iter, newton_reduction_factor,  &
        newton_tolerance, newton_verbosity

    double precision :: cg_tolerance
    integer :: cg_verbosity
    common /cg_config/ cg_tolerance, cg_verbosity

    ! The second dimension needs to be >= meqn.  It has to be specified
    ! as a compile-time constant, though, to allow use in a common block.
    character(len=2), dimension(4, 10) :: bc_options
    common /bc_config/ bc_options
    
    integer :: i, num_threads
    
    read(file_id, *) implicit_integration_scheme
    read(file_id, *) max_time_step_splits

    read(file_id, *) newton_max_iter
    read(file_id, *) newton_reduction_factor
    read(file_id, *) newton_tolerance
    read(file_id, *) newton_verbosity

    read(file_id, *) cg_tolerance
    read(file_id, *) cg_verbosity

    read(file_id, *) num_threads
    !$ call omp_set_num_threads(num_threads)
    
    do i = 1, meqn
        read(file_id, *) bc_options(:, i)
    end do
    
end subroutine setprob_implicit