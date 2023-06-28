module escapeVelocity
    
    use star_lib
    use const_def, only: dp
    use setup

    contains

        real(dp) function vescape(id,ierr,rval)
    
            type (star_info), pointer :: starptr
            integer, intent(in) :: id
            integer, intent(out) :: ierr
            
            real(dp), intent(in) :: rval
            real(dp) :: mR, rR, surfaceEscVel2
            real(dp) :: rIndex
            real(dp) :: rlowerbound, mlowerbound, rholowerbound
            real(dp) :: rupperbound, mupperbound, rhoupperbound
            real(dp) :: flowerbound, fupperbound
            real(dp) :: integralIncrement
            integer :: index
            !Value of G so that escape velocity is in km/s
            !assuming M in solar masses, density in g/cc, and radius in solar radii
            real(dp), parameter :: G1 = 1.126D6
            !Value of G so that escape velocity is in km/s
            !assuming M in solar masses, and radius in solar radii
            real(dp), parameter :: G2 = 1.9086D5

            call star_ptr(id, starptr, ierr)
            if (ierr /= 0) return

            !print*, "mR", starptr% r(1)/(6.957D10)

            mR = starptr% m(1)/(1.988D33)
            rR = starptr% r(1)/(6.957D10)

            surfaceEscVel2 = 2.D0 * G2 * mR/rR
            
            index  = starptr% nz
            rIndex = starptr% r(index)/(6.957D10)
            do while (rIndex < rval)
                rIndex = starptr% r(index)/(6.957D10)
                index = index - 1
            end do
            
            !print*, index, rIndex, rval
            
            if (index == 0) then
                index = 1
            end if

            vescape = 0.D0
            !print*, "vescape check1", vescape
            
            do while (index < starptr% nz)
                rlowerbound = starptr% r(index+1)/(6.957D10)
                mlowerbound = starptr% m(index+1)/(1.988D33)
                rholowerbound = starptr% rho(index+1)
            
                rupperbound = starptr% r(index)/(6.957D10)
                mupperbound = starptr% m(index)/(1.988D33)
                rhoupperbound = starptr% rho(index)

                flowerbound = mlowerbound/(rlowerbound**4._dp * rholowerbound)
                fupperbound = mupperbound/(rupperbound**4._dp * rhoupperbound)

                integralIncrement = (mupperbound - mlowerbound)/2._dp * (flowerbound + fupperbound)
                vescape = vescape + integralIncrement

                index = index + 1
            end do

            vescape = SQRT(G1*vescape/(2._dp*pi))
            print*, "vescape check2", vescape

         end function vescape

end module escapeVelocity