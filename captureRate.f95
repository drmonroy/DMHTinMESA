module captureRate
    
    use const_def, only: dp
    use setup
    use escapeVelocity
    use diffCapture
    use capture

    contains
        real(dp) function captureRateFunc(id, ierr, mchiIn)

            type (star_info), pointer :: starptr
            integer, intent(in) :: id
            integer, intent(out) :: ierr

            real(dp), intent(in) :: mchiIn
            real(dp) :: mH, mHe

            call star_ptr(id, starptr, ierr)
            if (ierr /= 0) return

            ! print*, "masscheck", starptr% m(1)/(1.988D33)

            mH = 1.674D-24
            mHe = 6.646477D-24

            captureRateFunc = &
                captureSpecies(id,ierr, mchiIn, mH, starptr% X) + &
                captureSpecies(id,ierr, mchiIn, mHe, starptr% Y)

         end function captureRateFunc

end module captureRate