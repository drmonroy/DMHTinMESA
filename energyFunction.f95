module energyFunction

    use const_def, only: dp
    use setup

    contains

        real(dp) function energyFunc(mchiIn, mspeciesIn, &
                                        nchiIn, &
                                        nspeciesIn, &
                                        rhoIn, &
                                        tempIn, &
                                        TchiIn)

            real(dp), intent(in) :: mchiIn, mspeciesIn, nchiIn, nspeciesIn, rhoIn, tempIn, TchiIn
            real(dp) :: sigma, sigmaconversion, energyconv
            real(dp) :: kB = 8.67*(10._dp**(-14._dp))

            !interaction cross section in GeV^-2
            if (ABS(mspeciesIn - 0.93911) < 0.001) then
                sigma = 2.568D-10
            else
                sigma = 2.568D-10 * mspeciesIn**4._dp * (mchiIn + 0.938272)**2._dp/ &
                ((mchiIn+mspeciesIn)**2._dp * 0.938272**2._dp)
            end if


            ! Factor of hbar^2 * c^3 in GeV^2 * cm^3 / s
            ! Converts cross section from GeV^-2 to cm^3/s
            sigmaconversion = 1.16733D-17
            sigma = sigmaconversion * sigma

            ! energy transferred per unit stellar mass
            ! in units of GeV/(g*s)
            energyFunc = 8._dp/rhoIn * SQRT(2._dp/pi) * & 
            ((mchiIn * mspeciesIn)/(mchiIn + mspeciesIn)**2._dp) * &
            nchiIn * nspeciesIn * sigma * kB * (TchiIn - tempIn) * &
            SQRT(kb*tempIn/mspeciesIn + kb*TchiIn/mchiIn)

            ! Conversion factor from GeV/(g*s) to ergs/(g*s)
            energyconv = 1.602D-3

            energyFunc = energyFunc*energyconv

            !print*, "energytest", energyFunc

            !energyFunc = 0.D0
        end function energyFunc

end module energyFunction