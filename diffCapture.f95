module diffCapture
    
    use const_def, only: dp
    use setup
    use chiFunction

    contains

        real(dp) function diffCap(mchiInput, mspeciesInput, vescapeInput, nucDensInput)
    
            real(dp), intent(in) :: mchiInput, mspeciesInput, vescapeInput, nucDensInput
            real(dp) :: mu, muplus, muminus
            real(dp) :: vbar, sigma, sigmaconversion
            real(dp) :: Apara2, Apara, Aplus, Aminus
            real(dp) :: vtilde, eta2, eta
            real(dp) :: dmDens
            
            !!!!INPUT PARAMETERS!!!!
            !speed of sun relative to dark matter distribution in km/s
            vtilde = 220._dp
            !velocity dispersion of the dark matter in km/s
            vbar = 270._dp
            !interaction cross section in GeV^-2
            if (ABS(mspeciesInput - 0.93911) < 0.001) then
                sigma = 2.568D-16
            else
                sigma = 2.568D-16 * mchiInput**2._dp * mspeciesInput**4._dp * (mchiInput + 0.938272)**2._dp/ &
                 ((mchiInput+mspeciesInput)**2._dp * mchiInput**2._dp * 0.938272**2._dp)
            end if
            !Dark matter density in GeV/cc
            dmDens = 0.4_dp
            !Dark matter density in #/cc
            dmDens = dmDens/mchiInput
            
            !!!!!!!!CONSTANTS!!!!!!!!
            !Conversion of cross section to cm^6/km*Rsol^3
            sigmaconversion = 7.74242D10
            
            !unitless parameters for masses
            mu = mspeciesInput/mchiInput
            muplus = (mu + 1._dp)/2._dp
            muminus = (mu - 1._dp)/2._dp
            !print*, "mu, mu+, mu-"
            !print*, mu, muplus, muminus
            
            !nondimensionalized version speed relative to dark matter
            eta2 = 1.5_dp * vtilde**2._dp/vbar**2._dp
            eta = SQRT(eta2)
            !print*, "eta"
            !print*, eta
            
            !unitless parameter
            Apara2 = (1.5_dp*vescapeInput**2._dp/vbar**2._dp) * mu/(muminus**2._dp)
            Apara = SQRT(Apara2)
            Aplus = Apara + eta
            Aminus = Apara - eta
            !print*, "A, A+, A-"
            !print*, Apara, Aplus, Aminus
            
            !convert the cross section
            sigma = sigma * sigmaconversion
            
            !!!!!!!CALCULATE DIFFERENTIAL CAPTURE RATE!!!!!!!
            diffCap = SQRT(6._dp/pi) * sigma * nucDensInput * dmDens * vbar * &
                      vescapeInput**2._dp/vbar**2._dp / (2._dp*eta*Apara2) * &
                      ( &
                          (Aplus*Aminus - 0.5_dp) * &
                          (chiFunc(-1._dp*eta,eta) - chiFunc(Aminus,Aplus)) + &
                          0.5_dp * Aplus*EXP(-1._dp*Aminus**2._dp) - &
                          0.5_dp * Aminus*EXP(-1._dp*Aplus**2._dp) - &
                          eta * EXP(-1._dp*eta2) &
                      )
            
        end function diffCap

end module diffCapture