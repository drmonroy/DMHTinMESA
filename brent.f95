module brent
    
    use const_def, only: dp
    use setup
    use luminosity
    use nchiFunction
    
    contains
            real(dp) function findTchi(id, ierr, mchiIn)

            type (star_info), pointer :: starptr
            integer, intent(in) :: id
            integer, intent(out) :: ierr            

            real(dp), intent(in) :: mchiIn
            real(dp) :: Tchilower, Tchiupper, rchi, a,b,c,d,e,p,q,r,s,fa,fb,fc, xm, minOne,minTwo
            real(dp) :: kB = 8.67*(10._dp**(-14.D0))
            real(dp) :: G = 3.59_dp*(10._dp**(-7.D0))
            real(dp) :: eps = epsilon(1.D0)
            real(dp) :: tol1, tol = 1.D-8
            real(dp) :: rIndex, Tcore, rhocore, mH, mHe
            integer :: ITmax=1000, iter, index
            integer :: maxiterCheck = 0
    
            call star_ptr(id, starptr, ierr)
            if (ierr /= 0) return

            !atomic masses in grams
            mH = 1.674D-24
            mHe = 6.646477D-24

            !Convert atomic masses to GeV
            nucmassconv = 5.61D23
            mH = mH*nucmassconv
            mHe = mHe*nucmassconv
        
            Tcore = starptr% T(starptr% nz)
            rhocore = starptr% rho(starptr% nz)
            rchi = SQRT((3._dp * kB * Tcore)/(2._dp*pi*G*rhocore * mchiIn))
            
            Tchiupper = Tcore
            ! print*, "Tchiupper", Tchiupper

            index  = 1
            rIndex = starptr% r(index)/(6.957D10)

            if (rchi > starptr% r(1)/(6.957D10)) then
               index = starptr% nz
               print*, "rchi limit", rchi
            else
               do while (rIndex > 2.5D0*rchi)
                  index = index + 1
                  rIndex = starptr% r(index)/(6.957D10)
               end do
            endif
            !print*, "rIndex", rIndex
            !print*, "rsurface", starptr% r(1)/(6.957D10)
            !print*, "rcore", starptr% r(starptr% nz)/(6.957D10)
            
            Tchilower = starptr% T(index)
            ! print*, "Tchilower", Tchilower
            
            a = Tchilower
            b = Tchiupper

            fa = lumin(id, ierr, mchiIn, mH, a) + lumin(id, ierr, mchiIn, mHe, a)
            fb = lumin(id, ierr, mchiIn, mH, b) + lumin(id, ierr, mchiIn, mHe, b)
    
            if ((fa<0._dp .and. fb<0._dp) .or. (fa>0._dp .and. fb>0._dp)) then
                print*, "Root must be bracketed!!"
                print*, "fa", fa
                print*, "fb", fb

                open(2, file="temperature.dat")
                write(2,*) "radius", "temperature"
                do k = 1, starptr% nz
                    write(2,*) starptr% r(k)/(6.957D10), starptr% T(k)
                end do

            else
                c = b
                fc = fb
                
                iter = 1
                do while (iter<=ITmax)
                    if ((fb<0._dp .and. fc<0._dp) .or. (fb>0._dp .and. fc>0._dp)) then
                        c=a
                        fc = fa
                        d = b - a
                        e = d
                    end if
                    
                    if (ABS(fc)<ABS(fb)) then
                        a = b
                        b = c
                        c = a
                        fa = fb
                        fb = fc
                        fc = fa
                    end if
                    
                    tol1 = 2._dp*eps*abs(b) + 0.5_dp*tol
                    xm = 0.5_dp*(c-b)
                    if ((ABS(xm)<= tol1) .or. (fb == 0._dp)) then
                        findTchi = b
                        maxiterCheck = iter
                        iter = ITmax
                    else
                        if ((ABS(e) >= tol1) .and. (ABS(fa) > ABS(fb))) then
                            s = fb/fa
                            if (a == c) then
                                p=2._dp*xm*s
                                q=1._dp - s
                            else
                                q=fa/fc
                                r=fb/fc
                                p=s*(2._dp*xm*q*(q-r) - (b-a)*(r-1._dp))
                                q=(q-1._dp)*(r-1._dp)*(s-1._dp);
                            end if
                            
                            if (p>0._dp) then
                            q = -q
                            end if
                            
                            p = ABS(p)
                            minOne = 3._dp * xm *q - ABS(tol1*q)
                            minTwo = ABS(e*q)
                            
                            if(2._dp*p < MIN(minOne,minTwo)) then
                                e = d
                                d = p/q
                            else
                                print*,"Interpolation failed, use bisection"
                                d = xm
                                e = d
                            end if
                        else
                            print*,"bounds decreasing too slowly, use bisection"
                            d = xm
                            e = d
                        end if
                        
                        a=b
                        fa = fb
                        if (ABS(d) > tol1) then
                            b = b + d
                        else
                            b = b + SIGN(tol1,xm)
                        end if
                        
                        fb = lumin(id, ierr, mchiIn, mH, b) + lumin(id, ierr, mchiIn, mHe, b)
                    end if
                    iter = iter + 1
                end do
                if (maxiterCheck == 0) then
                    print*,"Maximum number of iterations reached"
                else
                   ! print*, "value found! iterations used:", maxiterCheck
                endif
            end if
            findTchi = b
            ! print*, "Tchival", findTchi
        end function findTchi
end module brent