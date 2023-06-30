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

            !sets values in array nchi
            call nchiFunc(TchiIn)

            nchiCoeff = 0.D0
            do k = 1, numzones

                flowerbound = n_chi(k+1)/density(k+1)
                fupperbound = n_chi(k)/density(k)
            
                integralIncrement = (mass(k) - mass(k+1))/2.D0 * (flowerbound + fupperbound)
                nchiCoeff = nchiCoeff + integralIncrement
            end do
            
        end function nchiCoeff

end module nchiCoefficient