subroutine setprob

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

    double precision, dimension(4) :: d0_stencil, d1_stencil, d2_stencil,  &
        d3_stencil
    common /stencil_config/ d0_stencil, d1_stencil, d2_stencil, d3_stencil

    character*12 fname
    integer :: iunit


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

    ! -1/16, 9/16, 9/16, -1/16
    d0_stencil = (/ -0.0625d0,  0.5625d0,  0.5625d0, -0.0625d0 /)

    ! 1/24, -9/8, 9/8, -1/24
    d1_stencil = (/ 1.d0/24, -1.125d0, 1.125d0, -1.d0/24 /)

    ! 1/2, -1/2, -1/2, 1/2
    d2_stencil = (/ .5d0, -.5d0, -.5d0, .5d0 /)

    ! -1, 3, -3, 1
    d3_stencil = (/ -1.d0, 3.d0, -3.d0, 1.d0 /)

end subroutine setprob
