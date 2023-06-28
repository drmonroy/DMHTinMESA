module nchiFunction

    use const_def, only: dp
    use setup
    use escapeVelocity
    
    contains

        real(dp) function nchiFunc(id, ierr, rIn, TchiIn, mchiIn)
    
            type (star_info), pointer :: starptr
            integer, intent(in) :: id
            integer, intent(out) :: ierr

            real(dp), intent(in) :: rIn, mchiIn, TchiIn
            real(dp) :: kB = 8.67D-14
            real(dp) :: c = 3.D5
            real(dp) :: U,Uprime,vesc,vescCore,rIndex
            integer :: index
            
            call star_ptr(id, starptr, ierr)
            if (ierr /= 0) return

            !print*, "temp", starptr% T(1)

            vesc = vescape(id, ierr, rIn)
            !print*, "vesc", vesc

            !vescCore = vescape(id, ierr, starptr% r(starptr% nz)/(6.957D10) )
            
            U = 0.5D0 * mchiIn * (vesc/c)**2

            nchiFunc = EXP(-U/(kB * TchiIn))
            ! print*, "nchifunc", nchiFunc
        end function nchiFunc

end module nchiFunction