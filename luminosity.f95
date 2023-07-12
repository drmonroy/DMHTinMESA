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
            real(dp) :: tempFactor, interaction, massFactor
            integer :: j, k

            lumin_spec = 0.D0


            !print*, "Tchi", TchiIn

            if (.not. spindep) then
 
                do k = 1, numzones - 1
                
                    do j = 1, numspecies

                        if (A_spec(j) == 1) then

                        else

                            !upper bound
                            tempFactor = kB * (TchiIn - Temp(k)) * SQRT(kB * (Temp(k)/m_spec_GeV(j) + TchiIn/mchi))
                            interaction = EXP(-U(k)/(kB * TchiIn)) * n_spec(j,k) * sigma(j)
                            massFactor = (mchi * m_spec_GeV(j))/(mchi + m_spec_GeV(j))**2._dp
                            fupperbound = massFactor * tempFactor * interaction /  density(k)

                            !lower bound
                            tempFactor = kB * (TchiIn - Temp(k+1)) * SQRT(kB * (Temp(k+1)/m_spec_GeV(j) + TchiIn/mchi))
                            interaction = EXP(-U(k+1)/(kB * TchiIn)) * n_spec(j,k+1) * sigma(j)
                            massFactor = (mchi * m_spec_GeV(j))/(mchi + m_spec_GeV(j))**2._dp
                            flowerbound = massFactor * tempFactor * interaction /  density(k+1)

                            !calculate integral
                            integralIncrement = (mass(k) - mass(k+1)) * (flowerbound + fupperbound)/2.D0
                            lumin_spec(j) = lumin_spec(j) + integralIncrement
                        
                        end if

                    end do
                    
                end do
                lumin = SUM(lumin_spec)

            else

                lumin = 0.D0

                do k = 1, numzones - 1

                    !upper bound
                    tempFactor = kB * (TchiIn - Temp(k)) * SQRT(kB * (Temp(k)/m_prot + TchiIn/mchi))
                    interaction = EXP(-U(k)/(kB * TchiIn)) * n_H(k) * sigma_p
                    massFactor = (mchi * m_prot)/(mchi + m_prot)**2._dp
                    fupperbound = massFactor * tempFactor * interaction /  density(k)
                    !print*, "Tfact", TchiIn

                    !lower bound
                    tempFactor = kB * (TchiIn - Temp(k+1)) * SQRT(kB * (Temp(k+1)/m_prot + TchiIn/mchi))
                    interaction = EXP(-U(k+1)/(kB * TchiIn)) * n_H(k+1) * sigma_p
                    massFactor = (mchi * m_prot)/(mchi + m_prot)**2._dp
                    flowerbound = massFactor * tempFactor * interaction /  density(k+1)

                    !calculate integral
                    integralIncrement = (mass(k) - mass(k+1)) * (flowerbound + fupperbound)/2.D0
                    lumin = lumin  + integralIncrement

                end do

            end if
        end function lumin

end module luminosity