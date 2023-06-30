module energyFunction

    use const_def, only: dp
    use setup
    use nchiCoefficient

    real(dp), allocatable :: heat_transfer(:)

    contains

        subroutine energyFunc(TchiIn)
            
            real(dp) :: energy_spec(1:numspecies)
            real(dp), intent(in) :: TchiIn
            real(dp) :: nchi_Integral, N_DM

            nchi_Integral = nchiCoeff(TchiIn)
            N_DM = captureRate * starAge

            !energy transferred per unit stellar mass
            !in units of GeV/(g*s)
            allocate(heat_transfer(1:numzones))
            do k = 1, numzones
                do j = 1, numspecies
                    energy_spec(j) = 8._dp/density(k) * SQRT(2._dp/pi) * & 
                    ((mchi * m_spec_GeV(j))/(mchi + m_spec_GeV(j))**2._dp) * &
                    n_chi(k) * n_spec(j,k) * sigma(j) * c_vac * kB * (TchiIn - Temp(k)) * &
                    SQRT(kb*Temp(k)/m_spec_GeV(j) + kb*TchiIn/mchi)
                end do
                
                heat_transfer(k) = N_DM/nchi_Integral * GeV2erg * SUM(energy_spec)
            end do

            !energyFunc = 0.D0
        end subroutine energyFunc

end module energyFunction