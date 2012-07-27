program driver
          
     
    implicit none
    
    
    character :: skip_me
    integer :: mx, mbc, meqn, mwaves, maux

    integer, dimension(:), allocatable            :: mthlim
    double precision, dimension(:,:), allocatable :: q, aux
    double precision, dimension(:), allocatable   :: work

    integer :: i, mwork


    !==== Read size info from claw.data ====
    call opendatafile(55,'claw.data')
    read(55,*) skip_me
    read(55,*) mx
    
    do i = 1,14
        read(55,*) skip_me
    end do
    
    read(55,*) maux
    read(55,*) meqn
    read(55,*) mwaves

    do i = 1,4
        read(55,*) skip_me
    end do
    
    read(55,*) mbc
    close(55)



    !==== Allocate memory ====
    allocate( q(1-mbc:mx+mbc, meqn) )
    
    if (maux > 0) then
        allocate( aux(1-mbc:mx+mbc, maux) )
    else
        allocate( aux(1,1) )
    end if
    
    allocate( mthlim(mwaves) )

    mwork = (mx + 2*mbc) * (2 + 4*meqn + mwaves + meqn*mwaves)
    allocate( work(mwork) )



    call claw1ez(mx, meqn, mwaves, mbc, maux, mwork, mthlim, q, work, aux)

    stop 


end program driver
