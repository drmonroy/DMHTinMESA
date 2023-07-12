module brent
    
    use const_def, only: dp
    use setup
    use luminosity
    
    contains
            function findTchi(func,x1,x2,tol)          

            real(dp) :: findTchi(1:4)
            integer, parameter :: itmax = 100
            real(dp), parameter :: eps = epsilon(1.D0)
            real(dp), intent(in) :: x1, x2, tol
            real(dp), external :: func
            real(dp) :: a,b,c,d,e,p,q,r,s,fa,fb,fc, tol1, xm, minOne,minTwo
            integer :: iter
            
            a = x1
            b = x2
            fa = func(a)
            fb = func(b)
            !print*, "fa, fb", a,b

            if ((fa>0.D0 .and. fb>0.D0) .or. (fa<0.D0 .and. fb<0.D0)) then
                print*, "Root must be bracketed!!"
                print*, "fa", fa
                print*, "fb", fb

                open(2, file="temperature.dat")
                write(2,*) "radius", "temperature"
                do k = 1, 100
                    write(2,*) b/2.D0 + (k-1)*(a-b/2.D0)/99.D0, func(b/2.D0 + (k-1)*(a-b/2.D0)/99.D0)
                end do
                
                findTchi = [-1.0D0, 0.0D0, 0.0D0, 0.0D0]
                return
            end if

            c = b
            fc = fb
                
            do iter = 1, ITmax
                if ((fb>0.D0 .and. fc>0.D0) .or. (fb<0.D0 .and. fc<0.D0)) then
                    c = a
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
                    
                tol1 = 2.D0 * eps * abs(b) + 0.5D0 * tol
                xm = 0.5D0 * (c-b)

                if ((ABS(xm) <= tol1) .or. (fb == 0.D0)) then
                    ! Root (b) has been found.
                    ! Return this root plus the root from the other interval
                    ! (and both function evaluations).
                    ! This will be a or c depending on the signs of fa, fb, and fc.
                    ! The second root will be used to approximate the emoment
                    ! equation with a straight line so that a more accurate
                    ! root can be found when the extra energy is "too high".
                    if (SIGN(1.0D0,fb) == -SIGN(1.0D0,fa)) then
                        findTchi = [b, fb, a, fa]
                    else if (SIGN(1.0D0,fb) == -SIGN(1.0D0,fc)) then
                        findTchi = [b, fb, c, fc]
                    else
                        ! print*, "---***--- ALL ZBRENT ROOTS HAVE THE SAME SIGN ---***---"
                        if (ABS(fa) < ABS(fc)) then
                            findTchi = [b, fb, a, fa]
                        else
                            findTchi = [b, fb, c, fc]
                        end if
                    end if
                    
                    return
                end if

                if ((ABS(e) >= tol1) .and. (ABS(fa) > ABS(fb))) then
                    s = fb/fa
                    if (a == c) then
                        p=2.D0 * xm * s
                        q=1.D0 - s
                    else
                        q=fa/fc
                        r=fb/fc
                        p=s*(2.D0*xm*q*(q-r) - (b-a)*(r-1.D0))
                        q=(q-1.D0)*(r-1.D0)*(s-1.D0);
                    end if
                    
                    if (p>0.D0) then
                    q = -q
                    end if
                    
                    p = ABS(p)

                    minOne = 3.D0 * xm * q - ABS(tol1*q)
                    minTwo = ABS(e*q)
                    if(2.D0*p < MIN(minOne,minTwo)) then
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
                        
                fb = func(b)
                !print*, "fb2", b
                
            end do
            
            print*,"Maximum number of iterations reached"
            findTchi = [b, 0.0D0, 0.0D0, 0.0D0]

        end function findTchi

        real(dp) function calculate_Tchi()

            real(dp) :: Txhigh, Txlow, Ttmp
            real(dp), dimension(1:4) :: Tarray
            integer :: k, tries
            integer :: model_err = -1
            logical :: Tflag
            real(dp), parameter :: tol = 1.D-8
            real(dp) :: lumin_upper, lumin_lower, slope
        
            ! if ((model_err == modelNum) .and. (.not. Tflag)) then
            !     print*, 'Txlow > Txhigh or root must be bracketed'
            ! end if
        
            Txhigh = maxT*1.1D0
            Txlow = 0.25D0*Temp(rchi_Index)
            Tflag=.false.
            tries=0
        
            do while (.not. Tflag)
                tries = tries + 1
        
                !returns -1 if root not bracketed
                Tarray = findTchi(lumin, Txhigh, Txlow, tol)
                Ttmp = Tarray(1)

                !treat as root not bracketed
                !Tx close to Tmax and slope is shallow (most likely)
                if (Txlow > maxT) then
                    Ttmp = -1.0D0
                end if
        
                if (Ttmp > 0.0D0) then

                    lumin_upper = lumin(Ttmp)
                    lumin_lower = lumin(0.999D0*Ttmp)
                    slope = (lumin_upper - lumin_lower)/(0.001D0*Ttmp)
            
                    if (slope < 0.D0) then
                       Tflag = .TRUE.
                    else
                        Tflag = .FALSE.
                    end if

                    !expand the range and try again
                    Txlow = 1.05D0*Txlow
                
                else if (tries.EQ.1) then
        
                    !if root bracket error,
                    !and if this is the first try,
                    !expand range and try again
                    Txlow = Txlow/2.0
                    Txhigh = Txhigh*1.25
                
                else
                    !if root bracket error,
                    !and this is not the first try,
                    !go back a step, find root, make sure slope is negative
                    Txlow = Txlow/1.05
                    Tarray = findTchi(lumin, Txhigh, Txlow, tol)
                    Ttmp = Tarray(1)
                    
                    lumin_upper = lumin(Ttmp)
                    lumin_lower = lumin(0.999D0*Ttmp)
                    !print*, "lumin upper/lower", Ttmp
                    slope = (lumin_upper - lumin_lower)/(0.001D0*Ttmp)
            
                    if (slope < 0.D0) then
                       Tflag = .TRUE.
                    else
                        Tflag = .FALSE.
                    end if
                    
                    if (.NOT.Tflag) then
                        
                        ! If the code gets here, a suitable root cannot be found.
                        ! In most cases there are other problems with the star at this
                        ! timestep and MESA will re-do the step.
                        ! The model_err variable (below) keeps track of this and will
                        ! kill the run if the step is not re-done, since this means
                        ! a suitable DM temperature was never found.  
                        Ttmp = maxT
                        model_err = modelNum + 1

                        ! terminate the DO WHILE loop
                        exit

                    end if
                end if
            end do

            calculate_Tchi = Ttmp
            print*, "Tchi", calculate_Tchi, "lumin", lumin(calculate_Tchi)
            
        end function calculate_Tchi

end module brent