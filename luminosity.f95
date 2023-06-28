module luminosity
    
    use star_lib
    use const_def, only: dp
    use setup
    use energyFunction
    use nchiFunction

    contains

        real(dp) function lumin(id, ierr, mchiIn, mspeciesIn, TchiIn)

            type (star_info), pointer :: starptr
            integer, intent(in) :: id
            integer, intent(out) :: ierr  

            real(dp), intent(in) :: mchiin, mspeciesIn, TchiIn
            real(dp) :: mlowerbound, rholowerbound, nchilowerbound, nspecieslowerbound, Tlowerbound
            real(dp) :: mupperbound, rhoupperbound, nchiupperbound, nspeciesupperbound, Tupperbound
            real(dp) :: flowerbound, fupperbound
            real(dp) :: integralIncrement
            real(dp) :: mH, mHe, nucmassconv
            integer :: index

            call star_ptr(id, starptr, ierr)
            if (ierr /= 0) return

            !atomic masses in grams
            mH = 1.674D-24
            mHe = 6.646477D-24

            !Convert atomic masses to GeV
            nucmassconv = 5.61D23
            mH = mH*nucmassconv
            mHe = mHe*nucmassconv


            lumin = 0._dp
            index = 1

            do index = 1, starptr% nz - 1
                mlowerbound = starptr% m(index+1)/(1.988D33)
                rholowerbound = starptr% rho(index+1)
                nchilowerbound = nchiFunc(id, ierr, starptr% r(index+1)/(6.957D10), TchiIn, mchiIn)
                Tlowerbound = starptr% T(index+1)
                if (ABS(mspeciesIn - mH)<=0.001D0) then
                    nspecieslowerbound = starptr% X(index+1)*nucmassconv/mH
                else if (ABS(mspeciesIn - mHe)<=0.001D0) then
                    nspecieslowerbound = starptr% Y(index+1)*nucmassconv/mHe
                end if
                
                mupperbound = starptr% m(index)/(1.988D33)
                rhoupperbound = starptr% rho(index)
                nchiupperbound = nchiFunc(id, ierr, starptr% r(index+1)/(6.957D10), TchiIn, mchiIn)
                Tupperbound = starptr% T(index)
                if (ABS(mspeciesIn - mH)<=0.001D0) then
                    nspeciesupperbound = starptr% X(index) * starptr% rho(index) * nucmassconv/mH
                else if (ABS(mspeciesIn - mHe)<=0.001D0) then
                    nspeciesupperbound = starptr% Y(index) * starptr% rho(index) * nucmassconv/mHe
                end if

                


                flowerbound = energyFunc(mchiIn, mspeciesIn, &
                                nchilowerbound, &
                                nspecieslowerbound, &
                                rholowerbound, &
                                Tlowerbound, &
                                TchiIn)

                fupperbound = energyFunc(mchiIn, mspeciesIn, &
                                nchiupperbound, &
                                nspeciesupperbound, &
                                rhoupperbound, &
                                Tupperbound, &
                                TchiIn)


                !print*, "flowerbound", flowerbound
                !print*, "fupperbound", fupperbound

                integralIncrement = (mupperbound - mlowerbound) * (flowerbound + fupperbound)/2.D0
                lumin = lumin + integralIncrement

                !print*, "lumin", lumin
            end do
        end function lumin

end module luminosity