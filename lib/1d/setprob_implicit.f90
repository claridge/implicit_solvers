subroutine setprob_implicit(file_id)

    implicit none
    
    character :: implicit_integration_scheme
    common /implicit_config/ implicit_integration_scheme
    integer, intent(in) :: file_id
    
    integer :: mx, mbc, meqn
    double precision :: x_lower, dx
    common /claw_config/ mx, mbc, x_lower, dx, meqn

    integer :: newton_max_iter, newton_verbosity
    double precision :: newton_tolerance
    common /newton_config/ newton_max_iter, newton_tolerance, newton_verbosity

    double precision :: linear_solver_tolerance
    integer :: linear_solver_verbosity
    common /cg_config/ linear_solver_tolerance, linear_solver_verbosity

    ! The second dimension needs to be >= meqn.  It has to be specified
    ! as a compile-time constant, though, to allow use in a common block.
    character(len=2), dimension(2, 10) :: bc_options
    common /bc_config/ bc_options
    
    integer :: i
    
    read(file_id, *) implicit_integration_scheme

    read(file_id, *) newton_max_iter
    read(file_id, *) newton_tolerance
    read(file_id, *) newton_verbosity

    read(file_id, *) linear_solver_tolerance
    read(file_id, *) linear_solver_verbosity

    do i = 1, meqn
        read(file_id, *) bc_options(:, i)
    end do
    
end subroutine setprob_implicit