subroutine setprob

    implicit none
    
    double precision :: beta, kappa, delta, mu, right_film_height
    common /surfactant_config/ beta, kappa, delta, mu, right_film_height
    
    integer :: file_id

    file_id = 7
    call opendatafile(file_id, 'setprob.data')
    call setprob_implicit(file_id)
    
    read(file_id, *), beta
    read(file_id, *), kappa
    read(file_id, *), delta
    read(file_id, *), mu
    read(file_id, *), right_film_height
    
end subroutine setprob
