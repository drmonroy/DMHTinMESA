module nchiFunction

    use const_def, only: dp
    use setup
    real(dp) :: n_chi(1:50000)
    
    contains

        subroutine nchiFunc(TchiIn)
                 
            real(dp), intent(in) :: TchiIn
            integer :: k

            do k = 1, numzones
                n_chi(k) = EXP(-U(k)/(kB * TchiIn))
            end do

        end subroutine nchiFunc

end module nchiFunction