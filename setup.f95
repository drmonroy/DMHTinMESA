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

    !proton-DM cross section in cm^2
    real(dp), parameter :: sigma_p = 1.D-43
    !Dark matter mass in GeV
    real(dp), parameter :: mchi = 10.D0
    !speed of sun relative to dark matter distribution in km/s
    real(dp) :: vtilde = 220.D0
    !velocity dispersion of the dark matter in km/s
    real(dp) :: vbar = 270.D0
    !Dark matter density in GeV/cc
    real(dp) :: dmDens = 0.4D0

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
    real(dp), allocatable :: U(:), v_esc(:), g_value(:), diffCap(:), radius(:), mass(:), density(:), Temp(:)
    !number density of each species, at each point in the star
    real(dp), allocatable :: n_spec(:,:)
    !mass in GeV, mass in g, x-sec w/DM in cm^2, and atomic number of each species
    real(dp), allocatable :: m_spec_GeV(:), m_spec_g(:), sigma(:), A_spec(:)
    

    !variable to store values per species
    !when calculating differential capture rate
    real(dp), allocatable :: diffCap_spec(:)
    


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
        allocate(U(1:numzones),v_esc(1:numzones))
        allocate(g_value(1:numzones),radius(1:numzones),mass(1:numzones),density(1:numzones),Temp(numzones))

        numspecies = starptr% species
        allocate(m_spec_GeV(1:numspecies),m_spec_g(1:numspecies),sigma(1:numspecies),A_spec(1:numspecies))
        allocate(n_spec(1:numspecies,1:numzones))
        allocate(diffCap_spec(1:numspecies),diffCap(1:numzones))

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
        v_esc(1) = 2*radius(1)*g_value(1)


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
            !check against Troy w/spin-dependent case
            sigma(j) = sigma_p * (A_spec(j) * m_spec_GeV(j)/m_prot * (mchi + m_prot)/(mchi + m_spec_GeV(j)))**2
            
            do k = 1, numzones
                !mass fraction of each species at each point in the star
                n_spec(j,k) = starptr% xa(j,k)
                !number density of each species at each point in star
                n_spec(j,k) = n_spec(j,k) * starptr% rho(k) / m_spec_g(j)
            end do
        end do


        !!!!!!!CALCULATE DIFFERENTIAL CAPTURE RATE!!!!!!!

        !!!!CONVERSIONS!!!!
        vtilde = vtilde * 1.D5 !in cm/s
        vbar = vbar * 1.D5 !in cm/s
        dmDens = dmDens/mchi !in #/cc

        !nondimensionalized version speed relative to dark matter
        eta2 = 1.5D0 * (vtilde/vbar)**2.D0
        eta = SQRT(eta2)
        !print*, "eta"
        !print*, eta
        
        do k = 1, numzones
            
            do j = 1, numspecies
          
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
                        (v_esc(k)/vbar)**2.D0 / (2.D0*eta*Apara2) * &
                        ( &
                            (Aplus*Aminus - 0.5D0) * &
                            (chiFunc(-eta,eta) - chiFunc(Aminus,Aplus)) + &
                            0.5_dp * Aplus*EXP(-Aminus**2.D0) - &
                            0.5_dp * Aminus*EXP(-Aplus**2.D0) - &
                            eta * EXP(-eta2) &
                        )

            end do

            diffCap(k) = SUM(diffCap_spec)

        end do

        captureRate = 0.D0          

        do k = 1, numzones-1

            fupper = diffCap(k)/density(k)
            flower = diffCap(k+1)/density(k+1)

            intInc = (mass(k) - mass(k+1)) * (flower + fupper) / 2.D0
            captureRate = captureRate + intInc

        end do

    end subroutine set_vars


    real(dp) function chiFunc(chiFuncInputa,chiFuncInputb)
        real(dp), intent(in) :: chiFuncInputa, chiFuncInputb

        chiFunc = (SQRT(pi)/2)* (ERF(chiFuncInputb) - ERF(chiFuncInputa))

    end function chiFunc

    
end module setup