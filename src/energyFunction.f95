module energyFunction

    use const_def, only: dp
    use setup
    use nchiCoefficient

    real(dp) :: heat_transfer(1:50000)
    real(dp) :: temporary_heat_transfer
    real(dp) :: nchi_Integral

    contains

        subroutine energyFunc(TchiIn)
            
            real(dp) :: energy_spec(1:numspecies)
            real(dp), intent(in) :: TchiIn

            nchi_Integral = nchiCoeff(TchiIn)
            !print*, "captureRate, starAge", captureRate, starAge

            !energy transferred per unit stellar mass
            !in units of GeV/(g*s)
            heat_transfer = 0.D0
            if (.not. spindep) then
                do k = 1, numzones

                    do j = 1, numspecies

                        if (A_spec(j) == 1) then
                        
                        else

                            energy_spec(j) = 8.D0/density(k) * SQRT(2.D0/pi) * & 
                            ((mchi * m_spec_GeV(j))/(mchi + m_spec_GeV(j))**2.D0) * &
                            n_chi(k) * n_spec(j,k) * sigma(j) * c_vac * kB * (TchiIn - Temp(k)) * &
                            SQRT(kb*Temp(k)/m_spec_GeV(j) + kb*TchiIn/mchi)

                            heat_transfer(k) = heat_transfer(k) + N_DM/nchi_Integral * GeV2erg * energy_spec(j)
                            
                        end if

                    end do
                    
                end do
            else

                do k = 1, numzones
                    temporary_heat_transfer = 8.D0/density(k) * SQRT(2.D0/pi) * & 
                    (mchi * m_prot)/((mchi + m_prot)**2) * &
                    n_chi(k) * n_H(k) * sigma_p * c_vac * kB * (TchiIn - Temp(k)) * &
                    SQRT(kb*Temp(k)/m_prot + kb*TchiIn/mchi)

                    heat_transfer(k) = temporary_heat_transfer * N_DM/nchi_Integral * GeV2erg
                end do

            end if

            !heat_transfer = 0.D0
        end subroutine energyFunc

end module energyFunction