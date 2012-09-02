subroutine setprob

    implicit none

    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

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

    ! The second dimension needs to be >= meqn.  It has to be specified
    ! as a compile-time constant, though, to allow use in a common block.
    character(len=2), dimension(2, 10) :: bc_options
    common /bc_config/ bc_options

    double precision, dimension(4, 0:3) :: stencils
    common /stencil_config/ stencils

    character*12 fname
    integer :: iunit, i


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

    do i = 1, meqn
        read(7, *) bc_options(:, i)
    end do


    ! -1/16, 9/16, 9/16, -1/16
    ! 1/24, -9/8, 9/8, -1/24
    ! 1/2, -1/2, -1/2, 1/2
    ! -1, 3, -3, 1
    stencils(:, 0) = (/ -0.0625d0,  0.5625d0,  0.5625d0, -0.0625d0 /)
    stencils(:, 1) = (/ 1.d0/24, -1.125d0, 1.125d0, -1.d0/24 /)
    stencils(:, 2) = (/ .5d0, -.5d0, -.5d0, .5d0 /)
    stencils(:, 3) = (/ -1.d0, 3.d0, -3.d0, 1.d0 /)

end subroutine setprob
