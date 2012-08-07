program driver

! Driver program for a 1-dimensional Clawpack application.  We perform a few
! supplemental tasks in addition to the most basic driver functionality:
!   - Arrays are allocated dynamically, ensuring that maxmx==mx.
!   - Basic parameters of the spatial discretization are read into the common
!     block claw_config.
     
    implicit none

    common /claw_config/ mx, mbc, x_lower, dx, meqn
        
    character :: skip_me
    integer :: mx, mbc, meqn, mwaves, maux, i, mwork
    integer, dimension(:), allocatable :: mthlim
    double precision :: x_lower, x_upper, dx
    double precision, dimension(:), allocatable :: work
    double precision, dimension(:,:), allocatable :: q, aux


    call opendatafile(55,'claw.data')
    read(55, *) skip_me
    read(55, *) mx    
    do i = 1,14
        read(55, *) skip_me
    end do
    read(55, *) maux
    read(55, *) meqn
    read(55, *) mwaves
    do i = 1,2
        read(55, *) skip_me
    end do
    read(55, *) x_lower
    read(55, *) x_upper
    read(55, *) mbc
    close(55)

    dx = (x_upper - x_lower) / mx

    allocate( q(1-mbc:mx+mbc, meqn) )
    
    if (maux > 0) then
        allocate(aux(1-mbc:mx+mbc, maux))
    else
        allocate(aux(1,1))
    end if
    
    allocate( mthlim(mwaves) )

    mwork = (mx + 2 * mbc) * (2 + 4 * meqn + mwaves + meqn * mwaves)
    allocate(work(mwork))

    call claw1ez(mx, meqn, mwaves, mbc, maux, mwork, mthlim, q, work, aux)

    stop 

end program driver
