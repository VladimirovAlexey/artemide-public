module SnowFlake_Model
use LHA_alpha_snowflake
use IO_snowflake
implicit none
private


INCLUDE 'commonVariables.f90'

real(dp),allocatable,dimension(:)::NPparam

public::SnowFlake_Model_Initialize,SetNPparameters,PrintParameters

public::alpha
!public::Tu,Td,Ts,dTu,dTd,dTs
public::SplusU,SplusD,SplusS,SminusU,SminusD,SminusS
public::Tp,Tm

contains

!!!!! sets up all main definition
subroutine SnowFlake_Model_Initialize()

!!!!!! spacify the source of the alphaS-table
call ReadInfo("MSHT20nnlo_as118","/data/LHAPDF/share/LHAPDF/")
allocate(NPparam(0:17))

NPparam=0._dp

end subroutine SnowFlake_Model_Initialize


!!!!!! updating the value of NP-parameters
subroutine SetNPparameters(lambda)
real(dp),dimension(:),intent(in)::lambda

if(size(lambda)/=size(NPparam)) then
    write(*,*) ErrorString("Size of NP array is not correct! Nothing is done.", "SnowFlake_Model")
    write(*,*) size(lambda),size(NPparam)
    return
end if

NPparam=lambda

end subroutine SetNPparameters

subroutine PrintParameters()
write(*,*) "SNOWFLAKE PARAMETERS---->"
write(*,*) NPparam
write(*,*) "<-------"

end subroutine PrintParameters

!!!! This is a ``cut function'' which supposed to be
!!!! 1 at x1=x2=x3=0,
!!!! 0 at at border
!!!! to be consitent with the rest it must be symmetric in x1<->x3
function H(x,y,a,b,c)
    real*8,intent(in)::x,y
    real*8,intent(in)::a,b,c
    real*8::H

    if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
        H=((1-x**2)*(1-(x+y)**2))**a*(1-y**2)**b/(x**2+y**2+(x+y)**2)**c
    else
        H=0.d0
    end if
end function H

!!!! function S^+_u
function SplusU(x,y)
    real*8,intent(in)::x,y
    real*8::SplusU
    real*8,parameter::pi=3.141592653589793d0

    if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
        SplusU=H(x,y,NPparam(0),NPparam(1),NPparam(2))*(NPparam(3)+NPparam(4)*x*(-x-y))
    else
        SplusU=0.d0
    end if
end function SplusU

!!!! function S^-_u
function SminusU(x,y)
    real*8,intent(in)::x,y
    real*8::SminusU
    real*8,parameter::pi=3.141592653589793d0

    if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
        SminusU=H(x,y,NPparam(0),NPparam(1),NPparam(2))*(NPparam(5)*x+NPparam(6)*(-x-y))
    else
        SminusU=0.d0
    end if
end function SminusU

!!!! function S^+_d
function SplusD(x,y)
    real*8,intent(in)::x,y
    real*8::SplusD
    real*8,parameter::pi=3.141592653589793d0

    if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
        SplusD=H(x,y,NPparam(0),NPparam(1),NPparam(2))*(NPparam(7)+NPparam(8)*x*(-x-y))
    else
        SplusD=0.d0
    end if
end function SplusD

!!!! function S^-_d
function SminusD(x,y)
    real*8,intent(in)::x,y
    real*8::SminusD
    real*8,parameter::pi=3.141592653589793d0

    if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
        SminusD=H(x,y,NPparam(0),NPparam(1),NPparam(2))*(NPparam(9)*x+NPparam(10)*(-x-y))
    else
        SminusD=0.d0
    end if
end function SminusD

!!!! function S^+_d
function SplusS(x,y)
    real*8,intent(in)::x,y
    real*8::SplusS
    real*8,parameter::pi=3.141592653589793d0

    if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
        SplusS=H(x,y,NPparam(0),NPparam(1),NPparam(2))*(NPparam(11)+NPparam(12)*x*(-x-y))
    else
        SplusS=0.d0
    end if
end function SplusS

!!!! function S^-_d
function SminusS(x,y)
    real*8,intent(in)::x,y
    real*8::SminusS
    real*8,parameter::pi=3.141592653589793d0

    if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
        SminusS=H(x,y,NPparam(0),NPparam(1),NPparam(2))*(NPparam(13)*x+NPparam(14)*(-x-y))
    else
        SminusS=0.d0
    end if
end function SminusS

! !!!! function T_u
! function Tu(x,y)
!     real*8,intent(in)::x,y
!     real*8::Tu
!     real*8,parameter::pi=3.141592653589793d0
!
!     if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
!         Tu=NPparam(2)*H(x,y,NPparam(0),NPparam(1),NPparam(2))*(1.d0+NPparam(3)*(2*x+y))
!     else
!         Tu=0.d0
!     end if
! end function Tu
!
! !!!! function T_d
! function Td(x,y)
!     real*8,intent(in)::x,y
!     real*8::Td
!     real*8,parameter::pi=3.141592653589793d0
!
!     if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
!         Td=NPparam(8)*H(x,y,NPparam(0),NPparam(1),NPparam(2))*(1.d0+NPparam(9)*(2*x+y))
!     else
!         Td=0.d0
!     end if
! end function Td
!
! !!!! function T_s
! function Ts(x,y)
!     real*8,intent(in)::x,y
!     real*8::Ts
!
!     Ts=NPparam(14)*H(x,y,NPparam(0),NPparam(1),NPparam(2))
! end function Ts
!
! !!!! function Delta T_u
! function dTu(x,y)
!     real*8,intent(in)::x,y
!     real*8::dTu
!     real*8,parameter::pi=3.141592653589793d0
!
!     if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
!         dTu=NPparam(5)*H(x,y,NPparam(0),NPparam(1),NPparam(2))*(1.d0+NPparam(6)*(2*x+y))*y
!     else
!         dTu=0.d0
!     end if
! end function dTu
!
! !!!! function Delta T_d
! function dTd(x,y)
!     real*8,intent(in)::x,y
!     real*8::dTd
!     real*8,parameter::pi=3.141592653589793d0
!
!     if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
!         dTd=NPparam(11)*H(x,y,NPparam(0),NPparam(1),NPparam(2))*(1.d0+NPparam(12)*(2*x+y))*y
!     else
!         dTd=0.d0
!     end if
! end function dTd
!
! !!!! function Delta T_s
! function dTs(x,y)
!     real*8,intent(in)::x,y
!     real*8::dTs
!
!     dTs=NPparam(14)*H(x,y,NPparam(0),NPparam(1),NPparam(2))*y
! end function dTs
!
!!!! function T_{3F}^+
function Tp(x,y)
    real*8,intent(in)::x,y
    real*8::Tp

    if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
        Tp=NPparam(15)*H(x,y,NPparam(0),NPparam(1),NPparam(2))*(2*x+y)
    else
        Tp=0.d0
    end if
end function Tp

!!!! function T_{3F}^-
function Tm(x,y)
    real*8,intent(in)::x,y
    real*8::Tm

    if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
        Tm=NPparam(16)*H(x,y,NPparam(0),NPparam(1),NPparam(2))
    else
        Tm=0.d0
    end if
end function Tm

!!!!!! Function for alpha_s of QCD.
function alpha(mu)
real*8,intent(in)::mu
real*8::alpha
alpha=AlphaS(mu)
end function alpha

end module SnowFlake_Model
