module nchiFunction

    use const_def, only: dp
    use setup
    real(dp), allocatable :: n_chi(:)
    
    contains

        subroutine nchiFunc(TchiIn)
                 
            real(dp), intent(in) :: TchiIn
            integer :: k
            
            allocate(n_chi(1:numzones))

            do k = 1, numzones
                n_chi(k) = EXP(-U(k)/(kB * TchiIn))
            end do

        end subroutine nchiFunc

end module nchiFunction