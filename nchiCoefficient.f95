module nchiCoefficient

    use const_def, only: dp
    use setup
    use nchiFunction
    
    contains

        real(dp) function nchiCoeff(id, ierr, mchiIn, TchiIn)
    
            type (star_info), pointer :: starptr
            integer, intent(in) :: id
            integer, intent(out) :: ierr

            real(dp), intent(in) :: mchiIn, TchiIn
            real(dp) :: mlowerbound,mupperbound, rlowerbound,rupperbound
            real(dp) :: rholowerbound,rhoupperbound,flowerbound,fupperbound
            real(dp) :: integralIncrement
            integer :: index

            call star_ptr(id, starptr, ierr)
            if (ierr /= 0) return

            nchiCoeff = 0._dp
            index = 1
            
            do while (index < starptr% nz)
                !print*, "nchicoeffIndex", index
                mlowerbound = starptr% m(index+1)/(1.988D33)
                rlowerbound = starptr% r(index+1)/(6.957D10)
                rholowerbound = starptr% rho(index+1)
            
                mupperbound = starptr% m(index)/(1.988D33)
                rupperbound = starptr% r(index)/(6.957D10)
                rhoupperbound = starptr% rho(index)

                flowerbound = nchiFunc(id, ierr, rlowerbound, TchiIn, mchiIn)/rholowerbound
                fupperbound = nchiFunc(id, ierr, rupperbound, TchiIn, mchiIn)/rhoupperbound
            
                integralIncrement = (mupperbound - mlowerbound)/2._dp * (flowerbound + fupperbound)
                nchiCoeff = nchiCoeff + integralIncrement

                index = index + 1
            end do
            
        end function nchiCoeff

end module nchiCoefficient