module nchiCoefficient

    use const_def, only: dp
    use setup
    use nchiFunction
    
    contains

        real(dp) function nchiCoeff(TchiIn)

            real(dp), intent(in) :: TchiIn
            real(dp) :: flowerbound,fupperbound
            real(dp) :: integralIncrement
            integer :: k
            real(dp) :: x

            !sets values in array nchi
            call nchiFunc(TchiIn)

            x = 0.D0
            
            do k = 1, numzones-1

                flowerbound = n_chi(k+1)/density(k+1)
                fupperbound = n_chi(k)/density(k)
                !print*, "function1", n_chi(k+1)/density(k+1)
                !print*, "function2", n_chi(k),density(k)

                integralIncrement = (mass(k) - mass(k+1))/2.D0 * (flowerbound + fupperbound)
                if (flowerbound + fupperbound > HUGE(x)) then
                    print*, "inf error, k value:", k
                    exit
                end if
                x = x + integralIncrement

            end do
            nchiCoeff = x
            
        end function nchiCoeff

end module nchiCoefficient