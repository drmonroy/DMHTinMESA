module luminosity
    
    use star_lib
    use const_def, only: dp
    use setup

    contains

        real(dp) function lumin(TchiIn)

            real(dp), intent(in) :: TchiIn
            real(dp) :: flowerbound, fupperbound
            real(dp) :: lumin_spec(1:numspecies)
            real(dp) :: integralIncrement
            integer :: j, k

            lumin_spec = 0.D0

            do j = 1, numspecies-1
            
                do k = 1, numzones

                !upper bound
                tempFactor = kB * (Tchi - Temp(k)) * SQRT(kB * (Temp(k)/m_spec_GeV(j) - Tchi/mchi))
                interaction = EXP(-U(k)/(kB * Tchi)) * n_spec(j,k) * sigma(j)
                fupperbound = tempFactor * interaction /  density(k)

                !lower bound
                tempFactor = kB * (Tchi - Temp(k+1)) * SQRT(kB * (Temp(k+1)/m_spec_GeV(j) - Tchi/mchi))
                interaction = EXP(-U(k+1)/(kB * Tchi)) * n_spec(j,k) * sigma(j)
                flowerbound = tempFactor * interaction /  density(k+1)

                !calculate integral
                integralIncrement = (mass(k) - mass(k+1)) * (flowerbound + fupperbound)/2.D0
                lumin_spec(j) = lumin_spec(j) + integralIncrement

                end do
            end do
            lumin = SUM(lumin_spec)
        end function lumin

end module luminosity