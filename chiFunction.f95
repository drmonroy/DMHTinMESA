module chiFunction
    
    use const_def, only: dp
    use setup

    contains

        real(dp) function chiFunc(chiFuncInputa,chiFuncInputb)
    
            real(dp), intent(in) :: chiFuncInputa, chiFuncInputb
            
            chiFunc = (SQRT(pi)/2)* (ERF(chiFuncInputb) - ERF(chiFuncInputa))
                        
        end function chiFunc
        

end module chiFunction