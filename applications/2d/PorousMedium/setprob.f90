subroutine setprob

    !$ use omp_lib
    implicit none

    character :: implicit_integration_scheme
    integer :: max_time_step_splits
    common /implicit_config/ implicit_integration_scheme, max_time_step_splits

    integer :: newton_max_iter, newton_verbosity
    double precision :: newton_reduction_factor, newton_tolerance
    common /newton_config/ newton_max_iter, newton_reduction_factor,  &
        newton_tolerance, newton_verbosity

    double precision :: cg_tolerance
    integer :: cg_verbosity
    common /cg_config/ cg_tolerance, cg_verbosity

    character*12 fname
    integer :: iunit, num_threads


    iunit = 7
    fname = 'setprob.data'

    call opendatafile(iunit, fname)

    read(7, *) implicit_integration_scheme
    read(7, *) max_time_step_splits

    read(7, *) newton_max_iter
    read(7, *) newton_reduction_factor
    read(7, *) newton_tolerance
    read(7, *) newton_verbosity

    read(7, *) cg_tolerance
    read(7, *) cg_verbosity

    read(7, *) num_threads

    !$ call omp_set_num_threads(num_threads)
    
end subroutine setprob
