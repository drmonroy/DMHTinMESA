module setup

    use star_lib
    use const_def, only: dp
    use chem_def
    implicit none

    !value of pi
    real(dp), parameter :: pi = 4.D0*ATAN(1.D0)
    !Speed of light in a vacuum in cm/s
    real(dp), parameter :: c_vac = 2.99792D10
    !Newtonian gravitational constant in cm^3/g/s^2
    real(dp), parameter :: G_newt = 6.674D-8 
    !boltzmann constant in GeV/K
    real(dp), parameter :: kB = 8.617D-14
    !proton mass in GeV
    real(dp), parameter :: m_prot = 0.9382720813D0
    !conversion from amu to GeV
    real(dp), parameter :: amu2GeV = 0.931494102D0
    !conversion from GeV to g
    real(dp), parameter :: GeV2g = 1.783D-24
    !conversion from GeV to ergs
    real(dp), parameter :: GeV2erg = 1.602D-3
    !conversion from g to solar masses
    real(dp), parameter :: g2Msol = 1.988D33

    !proton-DM cross section in cm^2
    real(dp), parameter :: sigma_p = 1.D-43
    !Dark matter mass in GeV
    real(dp), parameter :: mchi = 10.D0
    !speed of sun relative to dark matter distribution in km/s
    real(dp) :: vtilde = 220.D5
    !velocity dispersion of the dark matter in km/s
    real(dp) :: vbar = 270.D5
    !Dark matter density in #/cc
    real(dp) :: dmDens = 0.4D0/mchi

    !whether the calculation is spin-dependent or not
    logical :: spindep = .true.

    !model number of star
    integer :: modelNum
    !number of cells in star model
    integer :: numzones
    !number of isotopes in current model
    integer :: numspecies
    !index of species in chem_isos
    integer :: chemj
    !age of star in seconds
    real(dp) :: starAge
    !characteristic radius of the dark matter in cm
    real(dp) :: rchi
    !cell number of the characteristic radius
    integer :: rchi_Index
    !max temperature? maximum of what?
    real(dp) :: maxT
    !total integrated capture rate
    real(dp) :: captureRate
    !potential energy, escape speed, gravitational g, differential capture rate, radius, mass, density, temperature
    !in: GeV, cm/s, cm/s^2, #/cc/s, cm, g, g/cc, K
    real(dp) :: U(1:50000), v_esc(1:50000), g_value(1:50000), diffCap(1:50000)
    real(dp) :: radius(1:50000), mass(1:50000), density(1:50000), Temp(1:50000), n_H(1:50000)
    !number density of each species, at each point in the star
    real(dp) :: n_spec(1:50,1:50000)
    !mass in GeV, mass in g, x-sec w/DM in cm^2, and atomic number of each species
    real(dp) :: m_spec_GeV(1:50), m_spec_g(1:50), sigma(1:50), A_spec(1:50)

    !variable to store values per species
    !when calculating differential capture rate
    real(dp) :: diffCap_spec(1:50)

    !variable to store data in star pointer
    real(dp) :: X_CTRL(1:10)

    contains

    subroutine set_vars(id,ierr)

        integer :: j, k
        real(dp) :: mu, muminus
        real(dp) :: Apara2, Apara, Aplus, Aminus
        real(dp) :: eta2, eta
        real(dp) :: flower, fupper, intInc

        type (star_info), pointer :: starptr
        integer, intent(in) :: id
        integer, intent(out) :: ierr

        call star_ptr(id, starptr, ierr)
        if (ierr /= 0) return

        numzones = starptr% nz
        numspecies = starptr% species

        !star age in seconds
        starAge = starptr% star_age * 3.154D7

        do k = 1, numzones
            ! what the hell is this? apparently in cm/s^2
            ! check that this matches GM_enc/r^2
            g_value(k) = starptr% grav(k)
            ! the radius values at each cell in cm
            radius(k) = starptr% r(k)
            ! the mass values at each cell in g
            mass(k) = starptr% m(k)
            ! the mass per volume at each cell in g/cc
            density(k) = starptr% rho(k)
            ! the temperature at each cell in K
            Temp(k) = starptr% T(k)
        end do


        !model number of star
        modelNum = starptr% model_number

        !characteristic radius of dark matter in cm
        rchi = SQRT(3*kB*Temp(numzones)*c_vac**2/(2*pi*G_newt*density(numzones)*mchi))

        rchi_Index = 1
        do while (radius(rchi_Index) > rchi)
            rchi_Index = rchi_Index + 1
        end do

        !maximum temperature in K (of what, though?)
        maxT = 10.D0**(starptr% log_max_temperature)

        !escape speed in cm/s
        v_esc(1) = 0.D0
        do k = 1, numzones-1
            flower = g_value(k+1)/(radius(k+1)**2 * density(k+1))
            fupper = g_value(k)/(radius(k)**2 * density(k))

            intInc = (mass(k) - mass(k+1))/2._dp * (flower + fupper)
            v_esc(k+1) = v_esc(k) + intInc
            v_esc(k+1) = SQRT(v_esc(k+1)/(2*pi) + 2*radius(1)*g_value(1))
        end do
        v_esc(1) = SQRT(2*radius(1)*g_value(1))

        !potential energy in GeV
        do k = 1, numzones
            U(k) = 0.5D0 * mchi * (v_esc(numzones)**2 - v_esc(k)**2)/c_vac**2
        end do


        do j = 1, numspecies
            !index of species j
            chemj = starptr% chem_id(j)
            !mass of species j in GeV
            m_spec_GeV(j) = chem_isos% W(chemj) * amu2GeV
            !mass of species j in g
            m_spec_g(j) = m_spec_GeV(j) * GeV2g
            !atomic # of species j
            A_spec(j) = chem_isos% Z_plus_N(chemj)
            !cross section with dark matter in cm^2
            
            if (chem_isos% Z_plus_N(chemj) == 1) then
                
            else
    
                sigma(j) = sigma_p * (A_spec(j) * m_spec_GeV(j)/m_prot * (mchi + m_prot)/(mchi + m_spec_GeV(j)))**2

            end if
            
            do k = 1, numzones
                !mass fraction of each species at each point in the star
                n_spec(j,k) = starptr% xa(j,k)
                !number density of each species at each point in star
                n_spec(j,k) = n_spec(j,k) * starptr% rho(k) / m_spec_g(j)
            end do
        end do

        do k = 1, numzones
            !mass fraction of hydrogen at each point in the star
            n_H(k) = starptr% X(k)
            !number density of hydrogen at each point in the star
            n_H(k) = n_H(k) * starptr% rho(k) / m_prot
        end do


        !!!!!!!CALCULATE CAPTURE RATE!!!!!!!

        if (.not. spindep) then

            !nondimensionalized version of speed relative to dark matter
            eta2 = 1.5D0 * (vtilde/vbar)**2.D0
            eta = SQRT(eta2)
            !print*, "eta"
            !print*, eta

            diffCap = 0.D0
            do k = 1, numzones
                    
                do j = 1, numspecies

                    if (A_spec(j) == 1) then

                    else
            
                        !unitless parameters for masses
                        mu = m_spec_GeV(j)/mchi
                        muminus = (mu - 1.D0)/2.D0
                        !print*, "mu, mu+, mu-"
                        !print*, mu, muplus, muminus

                        Apara2 = (1.5D0*(v_esc(k)/vbar)**2.D0) * mu/(muminus**2.D0)
                        Apara = SQRT(Apara2)
                        Aplus = Apara + eta
                        Aminus = Apara - eta
                        !print*, "A, A+, A-"
                        !print*, Apara, Aplus, Aminus
                        
                        diffCap_spec(j) = SQRT(6.D0/pi) * sigma(j) * n_spec(j,k) * dmDens * vbar * &
                                muminus**2.D0 / (3.D0*eta*mu) * &
                                ( &
                                    (Aplus*Aminus - 0.5D0) * &
                                    (chiFunc(-eta,eta) - chiFunc(Aminus,Aplus)) + &
                                    0.5D0 * Aplus*EXP(-Aminus**2.D0) - &
                                    0.5D0 * Aminus*EXP(-Aplus**2.D0) - &
                                    eta * EXP(-eta2) &
                                )
                        diffCap(k) = diffCap(k) + diffCap_spec(j)
                    
                    end if

                end do

            end do

            captureRate = 0.D0          

            do k = 1, numzones-1

                fupper = diffCap(k)/density(k)
                flower = diffCap(k+1)/density(k+1)

                intInc = (mass(k) - mass(k+1)) * (flower + fupper) / 2.D0
                captureRate = captureRate + intInc

            end do
        else
            captureRate = 25.D21 * (dmDens/0.4D0) * (sigma_p/1.D-43) * (v_esc(1)/618.D5) * (270.D5/vbar) * (mass(1)/g2Msol)
        end if

    end subroutine set_vars


    real(dp) function chiFunc(chiFuncInputa,chiFuncInputb)
        real(dp), intent(in) :: chiFuncInputa, chiFuncInputb

        chiFunc = (SQRT(pi)/2)* (ERF(chiFuncInputb) - ERF(chiFuncInputa))

    end function chiFunc

    
end module setup