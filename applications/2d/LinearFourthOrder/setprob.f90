subroutine setprob

    implicit none
    integer :: file_id

    file_id = file_id
    call opendatafile(file_id, 'setprob.data')
    call setprob_implicit(file_id)
    
end subroutine setprob
