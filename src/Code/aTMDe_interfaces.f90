!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!       arTeMiDe 3.03
!
!   Contains definitions interfaces common for the rest of artemide
!   Used in (almost) each artemide module.
!
!                               A.Vladimirov (13.10.2025)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module aTMDe_interfaces
use aTMDe_numerics
implicit none

public

!!! this is interface for function of 1 variables used in integrations
abstract interface
    function func_1D(x)
        import::dp
        real(dp) :: func_1D
        real(dp), intent(in) ::x
    end function func_1D
end interface

!!! this is interface for function of 1 variables used in integrations
abstract interface
    function func_2D(x,y)
        import::dp
        real(dp) :: func_2D
        real(dp), intent(in) ::x,y
    end function func_2D
end interface

!!! this is interface for function (-5:5) in the integration
abstract interface
    function func_1D_array5(x)
        import::dp
        real(dp),dimension(-5:5) :: func_1D_array5
        real(dp), intent(in) ::x
    end function func_1D_array5
end interface

!!! this is interface for optTMD-like function (-5:5)
abstract interface
    function optTMD(x,b,h)
        import::dp
        real(dp),dimension(-5:5) :: optTMD
        real(dp), intent(in) ::x,b
        integer,intent(in)::h
    end function optTMD
end interface

!!! this is interface for structure function
!!! with process0 being last 3 numbers of the process-numeration
abstract interface
    function strFUNC(Q2,qT,x1,x2,process0)
        import::dp
        real(dp)::strFUNC
        real(dp),intent(in)::Q2,qT,x1,x2
        integer,dimension(1:3),intent(in)::process0
    end function strFUNC
end interface

end module aTMDe_interfaces
