module capture
    
    use const_def, only: dp
    use setup
    use escapeVelocity
    use diffCapture
    
    contains

        real(dp) function captureSpecies(id,ierr,mchiIn,mspeciesIn,nucDensArrayIn)
    
            type (star_info), pointer :: starptr
            integer, intent(in) :: id
            integer, intent(out) :: ierr

            real(dp), dimension(:), intent(in) :: nucDensArrayIn
            real(dp) :: mchiIn, mspeciesIn
            real(dp) :: mupper, rhoupper, mlower, rholower
            real(dp) :: flower, fupper
            real(dp) :: intInc, captureVal
            real(dp) :: escapeVelupper, escapeVellower
            real(dp) :: dCdVupper, dCdVlower
            real(dp) :: nucDensupper, nucDenslower
            integer :: index
            
            call star_ptr(id, starptr, ierr)
            if (ierr /= 0) return

            captureval = 0.D0
            index = 1            

            do while (index < starptr% nz)
                mlower = starptr% m(index+1)/(1.988D33)
                rholower = starptr% rho(index+1)
                escapeVellower = vescape(id, ierr, starptr% r(index+1)/(6.957D10))
                !Incoming array is mass fraction; convert to number density
                nucDenslower = nucDensArrayIn(index+1)*rholower/mspeciesIn
                !Incoming mspeciesIn is in grams; convert to GeV
                dCdVlower = diffCap(mchiIn, mspeciesIn*5.61D23, escapeVellower, nucDenslower)
    
                mupper = starptr% m(index)/(1.988D33)
                rhoupper = starptr% rho(index)
                escapeVelupper = vescape(id, ierr, starptr% r(index)/(6.957D10))
                !Incoming array is mass fraction; convert to number density
                nucDensupper = nucDensArrayIn(index)*rhoupper/mspeciesIn
                !Incoming mspeciesIn is in grams; convert to GeV
                dCdVupper = diffCap(mchiIn, mspeciesIn*5.61D23, escapeVelupper, nucDensupper)
                
                flower = dCdVlower/rholower
                fupper = dCdVupper/rhoupper

                intInc = (mupper - mlower) * (flower + fupper) / 2.D0
                            
                captureval = captureval + intInc

                index = index + 1
            end do
            
            captureSpecies = captureval
        
        end function captureSpecies

end module capture