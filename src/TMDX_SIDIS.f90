!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 2.02
!
!    Evaluation of the TMD cross-section for SIDIS-like cross-sections
!
!    if you use this module please, quote ????.????
!
!    ver 1.2: release (AV, 15.12.2017)
!    ver 1.32: part of functions migrated to TMDF, rest updated (AV, 16.08.2018)
!    ver 2.02:                            (AV,16.08.2019)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDX_SIDIS
use aTMDe_Numerics
use IO_functions
use TMDF
use TMDF_KPC
use QCDinput
use EWinput
use IntegrationRoutines

implicit none
private

!!!!!! 1=accurate but slow
!!!!!! 2=fast but not accurate
#define INTEGRATION_MODE 2

!Current version of module
character (len=10),parameter :: moduleName="TMDX-SIDIS"
character (len=5),parameter :: version="v3.01"
!Last appropriate verion of constants-file
integer,parameter::inputver=31

real(dp) :: toleranceINT=0.0001d0
real(dp) :: toleranceGEN=0.0000001d0

integer::outputlevel
integer::messageTrigger

logical::started=.false.

real(dp)::hc2

!!other global parameters, which are defined upon initialization
integer:: orderH_global

!!!----------------------POWER CORRECTIONS---------------------------------
!!!======= LP factorization
!! inclusion of power corrections to the definition
!! corrQT ~ qT/Q in cross-section
!! corrM1 ~ M1/Q
!! corrM2 ~ M2/Q
!! exactX1Z1 ~ x1,z1 as at LP factorization
!! exactScales Use exact LP factorization scale =-2q^+q^-
logical:: corrQT,corrM1,corrM2,exactX1Z1, exactScales !!!if true include
!!!======= KPC factorization
logical:: useKPC !!! use the KPC formula

!!! number of sections for PT-integral by default
integer::NumPTdefault=6
real(dp)::ptMIN_global=0.00001d0

real(dp)::c2_global

integer::messageCounter
integer::GlobalCounter
integer::CallCounter


public::TMDX_SIDIS_ShowStatistic,TMDX_SIDIS_Initialize,&
    TMDX_SIDIS_IsInitialized,TMDX_SIDIS_ResetCounters,TMDX_SIDIS_SetScaleVariation

public::xSec_SIDIS,xSec_SIDIS_BINLESS,xSec_SIDIS_List,xSec_SIDIS_List_forharpy,xSec_SIDIS_BINLESS_List_forharpy



contains

function TMDX_SIDIS_IsInitialized()
    logical::TMDX_SIDIS_IsInitialized
    TMDX_SIDIS_IsInitialized=started
end function TMDX_SIDIS_IsInitialized

!! Initialization of the package
subroutine TMDX_SIDIS_Initialize(file,prefix)
    character(len=*),intent(in)::file
    character(len=*),intent(in),optional::prefix
    character(len=300)::path
    logical::initRequired
    character(len=8)::orderMain
    integer::FILEver

    if(started) return

    if(present(prefix)) then
        path=trim(adjustl(prefix))//trim(adjustr(file))
    else
        path=trim(adjustr(file))
    end if

    OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")
    !!! Search for output level
    call MoveTO(51,'*0   ')
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) FILEver
    if(FILEver<inputver) then
        write(*,*) 'artemide.'//trim(moduleName)//': const-file version is too old.'
        write(*,*) '             Update the const-file with artemide.setup'
        write(*,*) '  '
        CLOSE (51, STATUS='KEEP')
        ERROR STOP
    end if
    call MoveTO(51,'*p2  ')
    read(51,*) outputLevel    
    if(outputLevel>2) write(*,*) '--------------------------------------------- '
    if(outputLevel>2) write(*,*) 'artemide.',moduleName,version,': initialization started ... '
    call MoveTO(51,'*p3  ')
    read(51,*) messageTrigger

    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) hc2

    call MoveTO(51,'*10   ')
    call MoveTO(51,'*p1  ')
    read(51,*) initRequired
    if(.not.initRequired) then
        if(outputLevel>2) write(*,*)'artemide.',moduleName,': initialization is not required. '
        started=.false.
        return
    end if

    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) orderMain
    SELECT CASE(orderMain)
        CASE ("LO")
            orderH_global=0
        CASE ("NLO")
            orderH_global=1
        CASE ("NNLO")
            orderH_global=2
        CASE ("N2LO")
            orderH_global=2
        CASE ("NNNLO")
            orderH_global=3
        CASE ("N3LO") !! same as NNNLO
            orderH_global=3
        CASE ("N4LO")
            orderH_global=4
        CASE DEFAULT
            if(outputLevel>0) write(*,*) WarningString('try to set unknown order. Switch to NLO.',moduleName)
            orderH_global=1
    END SELECT
    if(outputLevel>1) write(*,*) '    artemide.TMDX_SIDIS: the used order is ',trim(orderMain)
    
    call MoveTO(51,'*p2   ')
    read(51,*) useKPC
    if(outputLevel>1 .and. useKPC) write(*,*) color('    artemide.TMDX_SIDIS: using TMD factorization with KPC',c_cyan)
    if(outputLevel>1 .and. useKPC) write(*,*) color('                      Please, cite [2307.13054]',c_cyan)
    if(outputLevel>1 .and. .not.(useKPC)) write(*,*) color('    artemide.TMDX_SIDIS: using TMD factorization at LP',c_cyan)

    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) toleranceGEN
    call MoveTO(51,'*p2  ')
    read(51,*) toleranceINT
    call MoveTO(51,'*p3  ')
    read(51,*) NumPTdefault
    call MoveTO(51,'*p4  ')
    read(51,*) ptMIN_global


    !!------------------------------------LP FACTORIZATION

    call MoveTO(51,'*C   ')
    !! qT correction in kinematics
    call MoveTO(51,'*p1   ')
    read(51,*) corrQT
    if(outputLevel>2 .and. corrQT .and. (.not.useKPC)) &
            write(*,*) '    artemide.TMDX_SIDIS: qT/Q corrections in kinematics are included.'
    !! Target mass corrections
    call MoveTO(51,'*p2   ')
    read(51,*) corrM1
    if(outputLevel>2 .and. corrM1) write(*,*) '    artemide.TMDX_SIDIS: target mass corrections in kinematics are included.'
    !! produced mass corrections
    call MoveTO(51,'*p3   ')
    read(51,*) corrM2
    if(outputLevel>2 .and. corrM2) write(*,*) '    artemide.TMDX_SIDIS: product mass corrections in kinematics are included.'
    !! qT correction in x1 z1
    call MoveTO(51,'*p4   ')
    read(51,*) exactX1Z1
    if(outputLevel>2 .and. exactX1Z1 .and. (.not.useKPC)) &
            write(*,*) '    artemide.TMDX_SIDIS: Exact LP values for x1,z1 are included.'
    !!exact values for scales
    call MoveTO(51,'*p5   ')
    read(51,*) exactScales
    if(outputLevel>2 .and. exactScales .and. (.not.useKPC)) &
          write(*,*) '    artemide.TMDX_SIDIS: Exact LP values of factorization scales variables are included.'

    !!------------------------------------KPC FACTORIZATION
    if(useKPC) then
      !!!! KPC must be computed with full expression for qT
      corrQT=.true.
      exactX1Z1=.true.
      exactScales=.true.


    call MoveTO(51,'*D   ')
    end if

    CLOSE (51, STATUS='KEEP')

    if(.not.EWinput_IsInitialized()) then
        if(outputLevel>1) write(*,*) '.. initializing EWinput (from ',moduleName,')'
        if(present(prefix)) then
            call EWinput_Initialize(file,prefix)
        else
            call EWinput_Initialize(file)
        end if
    end if

    if(useKPC) then
      if(.not.TMDF_KPC_IsInitialized()) then
        if(outputLevel>1) write(*,*) '.. initializing TMDF (from ',moduleName,')'
        if(present(prefix)) then
          call TMDF_KPC_Initialize(file,prefix)
        else
          call TMDF_KPC_Initialize(file)
        end if
      end if
    else
      if(.not.TMDF_IsInitialized()) then
        if(outputLevel>1) write(*,*) '.. initializing TMDF (from ',moduleName,')'
        if(present(prefix)) then
          call TMDF_Initialize(file,prefix)
        else
          call TMDF_Initialize(file)
        end if
      end if
    end if


    c2_global=1d0

    GlobalCounter=0
    CallCounter=0
    messageCounter=0

#if INTEGRATION_MODE==2
    write(*,*)  color('--------------------------------------------------------',c_red)
    write(*,*)  color('----------------------  WARNING!  ----------------------',c_red)
    write(*,*)  color('-- TMDX_SIDIS is in the approximate integration mode  --',c_red)
    write(*,*)  color('--            Faster, but lower precision.            --',c_red)
    write(*,*)  color('--    Switch to default version by changing flag      --',c_red)
    write(*,*)  color('-- INTEGRATION_MODE in TMDX_SIDIS.f90, and recompile  --',c_red)
    write(*,*)  color('--------------------------------------------------------',c_red)
#endif
    
    started=.true.
    write(*,*)  color('----- arTeMiDe.TMD_SIDIS '//trim(version)//': .... initialized',c_green)
end subroutine TMDX_SIDIS_Initialize

subroutine TMDX_SIDIS_ResetCounters()
if(outputlevel>2) call TMDX_SIDIS_ShowStatistic()
GlobalCounter=0
CallCounter=0
messageCounter=0
end subroutine TMDX_SIDIS_ResetCounters
  

subroutine TMDX_SIDIS_ShowStatistic()

write(*,'(A,ES12.3)') 'TMDX SIDIS statistics   total calls of point xSec  :  ',Real(GlobalCounter)
write(*,'(A,ES12.3)') '                              total calls of xSecF :  ',Real(CallCounter)
write(*,'(A,F12.3)')  '                                         avarage M :  ',Real(GlobalCounter)/Real(CallCounter)
end subroutine TMDX_SIDIS_ShowStatistic
  
  
!!!!Call this after TMD initializetion but before NP, and X parameters
subroutine TMDX_SIDIS_SetScaleVariation(c2_in)
real(dp),intent(in)::c2_in

if(outputLevel>1) write(*,*) 'TMDX_SIDIS: scale variation constant c2 reset:',c2_in

if(c2_in<0.1d0 .or. c2_in>10.d0) then
    if(outputLevel>0) write(*,*) WarningString('variation in c2 is enourmous. c2 is set to 2',moduleName)
    c2_global=2d0
else
    c2_global=c2_in
end if

end subroutine TMDX_SIDIS_SetScaleVariation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! PROCESS DEFINITION !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!FUNCTIONS FOR OPERATION WITH KINEMATICS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! function makes kinematic array from the given set of pT,s,Q,x,z
!!! array has 13 often appearing entries
! 1 = pT
! 2 = Q
! 3 = Q^2
! 4 = x
! 5 = z
! 6 = y     =Q^2/x(s-M^2)
! 7 = epsilon    = (1-y-gamma2 y^2/4)/(1-y+y^2/2+gamma2 y^2/4)
! 8 = gamma2    = (2 x M/Q)^2
! 9=rho2*gamma2    = (m/z/Q)^2*gamma2
! 10=rhoT2*gamma2    = (m^2+pt^2)/(z Q)^2*gamma2
! 11=sM2    = s-M^2
! 12=M2-target
! 13=M2-product
pure function kinematicArray(pT,s,z,x,Q,M2target_in,M2product_in)
real(dp),dimension(1:13)::kinematicArray
real(dp),intent(in)::pT,s,z,x,Q,M2target_in,M2product_in
real(dp)::Q2,y,varepsilon,gamma2,rho2,rhoPEPR2,sM2,M2target,M2product

Q2=Q**2
if(corrM1 .and. M2target_in>toleranceGEN) then
    M2target=M2target_in
    sM2=s-M2target
    gamma2=4d0*M2target*x**2/Q2
else
    M2target=0._dp
    gamma2=0._dp
    sM2=s
end if

y=Q2/x/sM2!YfromSXQ2(sM2,x,Q2)

varepsilon=(1d0-y-y**2*gamma2*0.25d0)/(1d0-y+y**2*(0.5d0+0.25d0*gamma2))

if(corrM2 .and. M2product_in>toleranceGEN) then
    M2product=M2product_in
    rho2=gamma2*M2product/Q2/z**2
else
    M2product=0._dp
    rho2=0._dp
end if

if(corrQT) then
    rhoPEPR2=rho2+gamma2*(pT/z/Q)**2
else
    rhoPEPR2=rho2
end if

kinematicArray=(/pT,Q,Q2,x,z,y,varepsilon,gamma2,rho2,rhoPEPR2,sM2,M2target,M2product/)

end function kinematicArray
  
!!! xy(s-M^2)=Q^2
pure function YfromSXQ2(sM2,x,Q2)
real(dp),intent(in)::sM2,x,Q2
real(dp)::YfromSXQ2
YfromSXQ2=Q2/x/sM2
end function YfromSXQ2

pure function XfromSYQ2(sM2,y,Q2)
real(dp),intent(in)::sM2,y,Q2
real(dp)::XfromSYQ2
XfromSYQ2=Q2/y/sM2
end function XfromSYQ2

pure function QfromSXY(sM2,x,y)
real(dp),intent(in)::sM2,x,y
real(dp)::QfromSXY
QfromSXY=Sqrt(sM2*x*y)
end function QfromSXY

!!!!!Evaluate the parameters of the TMD factorized Fourier-integral
!!qT     =pT/z sqrt( ( 1+gamma2) / (1-gamma2 rho2))
!!X1=x*fac1
!!Z1=z*fac1*(...)
pure subroutine CalculateX1Z1qT(x1,z1,qT,var)
real(dp),intent(out)::x1,z1,qT
real(dp),dimension(1:13),intent(in)::var
real(dp)::fac1,fac2

qT=var(1)/var(5)*Sqrt((1d0+var(8))/(1d0-var(9)))

!!fac1 and fac2 are functions that appear in definition of x1, and z1
!! fac1= -2/gamma^2(1-sqrt[1+gamma^2(1-qT^2/Q^2)])
!! fac2= (1-sqrt[1-rho2 gamma^2])/(2 (1-qT^2/Q^2))
if(exactX1Z1) then
  if(corrQT .and. corrM1 .and. var(8)>toleranceGEN) then
      fac1=-2d0/var(8)*(1d0-sqrt(1d0+var(8)*(1d0-(qT/var(2))**2)))

      fac2=(1d0+sqrt(1d0-var(9)))/(2d0*(1d0-(qT/var(2))**2))
  else if(corrQT) then
      fac1=(1d0-(qT/var(2))**2)

      fac2=1.d0/(1d0-(qT/var(2))**2)
  else
      fac1=-2d0/var(8)*(1d0-sqrt(1d0+var(8)))

      fac2=(1d0+sqrt(1d0-var(9)))*0.5d0
  end if
else
  fac1=1._dp
  fac2=1._dp
end if

x1=var(4)*fac1
z1=var(5)*fac1*fac2

!if(x1<0.0001d0) write(*,*) '>>>>>>>>>>>>>>>',x1,z1,qT,var

end subroutine CalculateX1Z1qT
  
!!!! update a given kinematic array with new value of x.
pure subroutine SetX(x,var)
real(dp) ,dimension(1:13),intent(inout)::var
real(dp),intent(in)::x
real(dp)::g2

!   var=kinematicArray(var(1),var(11)+var(12),var(5),x,var(2),var(12),var(13))
g2=var(8)!old gamma2
var(4)=x
!!!pt same
!!!Q2 same
!!!sM2 same
!!!z same
!!! masses same

!!!y
var(6)=YfromSXQ2(var(11),x,var(3))

if(corrM1 .and. g2>toleranceGEN) then
    !!!gamma
    var(8)=4d0*var(12)*x**2/var(3)

    !! rescale rho'2 with new gamma
    var(9)=var(9)*var(8)/g2
    var(10)=var(10)*var(8)/g2
else
    var(8)=0._dp
    var(9)=0._dp
    var(10)=0._dp
end if

!!epsilon
var(7)=(1d0-var(6)-var(6)**2*var(8)*0.25d0)/(1d0-var(6)+var(6)**2*(0.5d0+0.25d0*var(8)))

end subroutine SetX
  
!!!! update a given kinematic array with new value of Q.
pure subroutine setQ(Q,var)
real(dp) ,dimension(1:13),intent(inout)::var
real(dp),intent(in)::Q

!!! in the case of Q, the array is updated completely
!!! thus we just reconstruct it
var=kinematicArray(var(1),var(11)+var(12),var(5),var(4),Q,var(12),var(13))
end subroutine setQ
  
!!!! update a given kinematic array with new value of Q2.
pure subroutine setQ2(Q2,var)
real(dp) ,dimension(1:13),intent(inout)::var
real(dp),intent(in)::Q2
!!! in the case of Q, the array is updated completely
!!! thus we just reconstruct it
var=kinematicArray(var(1),var(11)+var(12),var(5),var(4),sqrt(Q2),var(12),var(13))

end subroutine setQ2
  
!!!! update a given kinematic array with new value of Z.
pure subroutine SetZ(z,var)
real(dp) ,dimension(1:13),intent(inout)::var
real(dp),intent(in)::z
real(dp)::ratio2

ratio2=(var(5)/z)**2 !!! the factors rho are proportional to 1/z^2, to updated them I rescale by (z_old/z_new)^2
var(5)=z

!!pt same
!!Q2 same
!!gamma2 same
!!x same
!!y same
!!varepsilon same
!! sM2 and masses same

!!rho2 new
var(9)=var(9)*ratio2
var(10)=var(10)*ratio2

end subroutine SetZ
  
!!!! update a given kinematic array with new value of pt.
pure subroutine SetPT(pt,var)
real(dp),dimension(1:13),intent(inout)::var
real(dp),intent(in)::pt
var(1)=pt

!!Q2 same
!!gamma2 same
!!x same
!!y same
!!varepsilon same
!! sM2 and masses same
!!! rho2 same

!!!rho perp new
if(corrQT) then
    var(10)=var(9)+var(8)*(var(1)/var(5)/var(2))**2
else
    var(10)=var(9)
end if
end subroutine SetPT
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS FOR PREFACTORS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! hard coefficeint taken from 1004.3653 up to 2-loop
!!! it takes global values of Q,order
!!! NOTE it uses Nf=3(fixed)
function HardCoefficientSIDIS(mu)
real(dp)::HardCoefficientSIDIS,mu,alpha,LQ!=Log[Q^2/mu^2]=-2Log[c1]

HardCoefficientSIDIS=1.d0
if(orderH_global>=1) then
    LQ=-2d0*LOG(c2_global)
    alpha=As(mu*c2_global)
    HardCoefficientSIDIS=HardCoefficientSIDIS+alpha*&
    (-16.946842488404727d0 + 8d0*LQ - 2.6666666666666665d0*LQ**2)
if(orderH_global>=2) then
    HardCoefficientSIDIS=HardCoefficientSIDIS+alpha**2*&
    (-116.50054911601637d0 + 46.190372772820254d0*LQ + 16.843858371984233d0*LQ**2&
    -13.333333333333334d0*LQ**3 + 3.5555555555555554d0*LQ**4)
if(orderH_global>=3) then
    HardCoefficientSIDIS=HardCoefficientSIDIS+alpha**3*&
    (-4820.715927678687 + 2492.274201933993*LQ + 44.19495641116441*LQ**2 &
    - 237.22228827339313*LQ**3 + 43.33848430014775*LQ**4 + 7.111111111111111*LQ**5 &
    -3.1604938271604937*LQ**6)
if(orderH_global>=4) then
    HardCoefficientSIDIS=HardCoefficientSIDIS+alpha**4*&
    (26391.759725461765 - 31391.21814540276*LQ + 16136.794429475773*LQ**2 &
    - 3922.584672164565*LQ**3 - 98.8284343286739*LQ**4 + 250.360398809412*LQ**5 &
    - 47.07700456533447*LQ**6 - 3.950617283950617*LQ**7 + 1.8436213991769548*LQ**8)
end if
end if
end if
end if
end function HardCoefficientSIDIS
  
!!!!! Prefactor 2 is (universal part) x H
function PreFactor2(var,process,x1,z1,qT)
real(dp),dimension(1:13),intent(in)::var
integer,dimension(1:4),intent(in)::process
real(dp),intent(in)::x1,z1,qT
real(dp)::PreFactor2,uniPart,fac1,scaleMu

!!!! universal part

!!!! If exact scales, zeta*zeta=Q^2-qT^2.
if(exactScales) then
  scaleMu=sqrt(var(3)-qT**2)
else
  scaleMu=var(2)
end if

SELECT CASE(process(1))
  case(0)
      PreFactor2=1._dp
  CASE(1)
      !!! uniPart is the prefactor for the cross-section
      !2 pi aEm^2/Q^4 y^2/(1-epsilon)/sqrt[1-g2*rho2]*z1/z
      uniPart=pix2*alphaEM(scaleMu)**2/(var(3)**2)*var(6)**2/((1d0-var(7))*sqrt(1-var(10)))*(z1/var(5))

      !! fac1 is the prefactor for the unpolarized expression
      !!! this is 1+qT^2/Q^2(e-gamma^2/2)/(1+gamma^2)
      if(corrQT) then
          fac1=1d0+(qT**2/var(3))*(var(7)-0.5d0*var(8))/(1+var(8))
      else
          fac1=1._dp
      end if

      PreFactor2=uniPart*fac1*HardCoefficientSIDIS(scaleMu)*&
      hc2*1d9!from GeV to pbarn

  CASE(2)
      !!!! same as for case(1) but with factor Q2/y

      !!! uniPart is the prefactor for the cross-section
      !2 pi aEm^2/Q^2 y/(1-epsilon)/sqrt[1-g2*rho2]*z1/z
      uniPart=pix2*alphaEM(scaleMu)**2/var(3)*var(6)/((1d0-var(7))*sqrt(1-var(10)))*(z1/var(5))

      !! fac1 is the prefactor for the unpolarized expression
      !!! this is 1+qT^2/Q^2(e-gamma^2/2)/(1+gamma^2)
      if(corrQT) then
          fac1=1d0+(qT**2/var(3))*(var(7)-0.5d0*var(8))/(1+var(8))
      else
          fac1=1._dp
      end if

      PreFactor2=uniPart*fac1*HardCoefficientSIDIS(scaleMu)*&
      hc2*1d9!from GeV to pbarn

  CASE(3)
      !!!! same as for case(1) but with factor x/y

      !!! uniPart is the prefactor for the cross-section
      !2 pi aEm^2/Q^2 y/(1-epsilon)/sqrt[1-g2*rho2]*z1/z
      uniPart=pix2*alphaEM(scaleMu)**2/(var(3)**2)*(var(6)*var(4))/((1d0-var(7))*sqrt(1-var(10)))*(z1/var(5))

      !! fac1 is the prefactor for the unpolarized expression
      !!! this is 1+qT^2/Q^2(e-gamma^2/2)/(1+gamma^2)
      if(corrQT) then
          fac1=1d0+(qT**2/var(3))*(var(7)-0.5d0*var(8))/(1+var(8))
      else
          fac1=1._dp
      end if

      PreFactor2=uniPart*fac1*HardCoefficientSIDIS(scaleMu)*&
      hc2*1d9!from GeV to pbarn

  CASE DEFAULT
      PreFactor2=0._dp
      write(*,*) ErrorString(' unknown process p1=',moduleName), process(2)
      ERROR STOP color('Evaluation stop.',c_red_bold)
END SELECT
end function PreFactor2


!!!!! Prefactor for KPC-case is (universal part) x H
function PreFactorKPC(var,process,x1,z1,qT)
real(dp),dimension(1:13),intent(in)::var
integer,dimension(1:4),intent(in)::process
real(dp),intent(in)::x1,z1,qT
real(dp)::PreFactorKPC,uniPart,scaleMu

!!!! universal part

!!!-----------------------------------------------------------------
!!!------  Prefactors KPC case (TO BE CHECKED!!)
!!!-----------------------------------------------------------------

!!!! zeta*zeta=Q^2
scaleMu=var(2)

SELECT CASE(process(1))
  case(0)
      PreFactorKPC=1._dp
  CASE(1)
      !!! uniPart is the prefactor for the cross-section
      !pi^2 aEm^2/Q^2 y^2/(1-epsilon)/sqrt[1-g2*rho2]
      uniPart=pi2*alphaEM(scaleMu)**2/(var(3))*var(6)**2/((1d0-var(7))*sqrt(1-var(10)))

      PreFactorKPC=uniPart*HardCoefficientSIDIS(scaleMu)*&
      hc2*1d9!from GeV to pbarn

  CASE(2)
      !!!! same as for case(1) but with factor Q2/y

      uniPart=pi2*alphaEM(scaleMu)**2/(var(3))*var(6)/((1d0-var(7))*sqrt(1-var(10)))

      PreFactorKPC=uniPart*HardCoefficientSIDIS(scaleMu)*&
      hc2*1d9!from GeV to pbarn

  CASE(3)
      !!!! same as for case(1) but with factor x/y

      !!! uniPart is the prefactor for the cross-section
      uniPart=pix2*alphaEM(scaleMu)**2/(var(3))*(var(6)*var(4))/((1d0-var(7))*sqrt(1-var(10)))

      PreFactorKPC=uniPart*HardCoefficientSIDIS(scaleMu)*&
      hc2*1d9!from GeV to pbarn

  CASE DEFAULT
      PreFactorKPC=0._dp
      write(*,*) ErrorString(' unknown process p1=',moduleName), process(2)
      ERROR STOP color('Evaluation stop.',c_red_bold)
END SELECT
end function PreFactorKPC
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CUTS RELATED FUNCTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!The cuts are ymin<y<ymax,Wmin<W2<Wmax

!!! checks the value of x against cut constaints from below and return the maximal allowed value
!!! argument xmin is xmin vs. which we compare cuts.
pure function xMinWithCuts(xmin,var,cutParam)
real(dp),dimension(1:13),intent(in)::var
real(dp),dimension(1:4),intent(in)::cutParam
real(dp),intent(in)::xmin
real(dp)::x1,x2,xMinWithCuts

x1=var(3)/cutParam(2)/var(11)
x2=var(3)/(var(3)+cutParam(4)-var(12))

xMinWithCuts=max(xmin,x1,x2)
end function xMinWithCuts
  
!!! checks the value of x against cut constaints from above and return the minimal allowed value
!!! argument xmax is xmin vs. which we compare cuts.
function xMaxWithCuts(xmax,var,cutParam)
real(dp),dimension(1:13),intent(in)::var
real(dp),dimension(1:4),intent(in)::cutParam
real(dp),intent(in)::xmax
real(dp)::x1,x2,xMaxWithCuts
x1=var(3)/cutParam(1)/var(11)
x2=var(3)/(var(3)+cutParam(3)-var(12))

xMaxWithCuts=min(xmax,x1,x2)
end function xMaxWithCuts
  
!!! checks the value of x against cut constaints from above and return the minimal allowed value
!!! argument Qmin is Qmin vs. which we compare cuts.
!!! argument xmin is xmin vs. which we compare cuts.
pure function QMinWithCuts(xmin,Qmin,var,cutParam)
real(dp),dimension(1:13),intent(in)::var
real(dp),dimension(1:4),intent(in)::cutParam
real(dp),intent(in)::xmin,Qmin
real(dp)::Q1,Q2,QMinWithCuts

Q1=sqrt(xmin*cutParam(1)*var(11))
Q2=sqrt(xmin*(cutParam(3)-var(12))/(1d0-xmin))

QMinWithCuts=max(Qmin,Q1,Q2)
end function QMinWithCuts
  
!!! checks the value of Q against cut constaints from below and return the maximal allowed value
!!! argument Qmax is Qmax vs. which we compare cuts.
!!! argument xmax is xmax vs. which we compare cuts.
pure function QMaxWithCuts(xmax,Qmax,var,cutParam)
real(dp),dimension(1:13),intent(in)::var
real(dp),dimension(1:4),intent(in)::cutParam
real(dp),intent(in)::xmax,Qmax
real(dp)::Q1,Q2,QMaxWithCuts

Q1=Sqrt(xmax*cutParam(2)*var(11))
Q2=sqrt(xmax*(cutParam(4)-var(12))/(1d0-xmax))

QMaxWithCuts=min(Qmax,Q1,Q2)
end function QMaxWithCuts
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS CALCULATING CROSS-SECTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------UNTIGRATED------------------------------------------------------------------

!!! this is help function which evaluate xSec at single qt (without lists) with prefactor2
!!! this is extended (and default) version of xSec, which include all parameters
function xSec(var,process)
real(dp):: xSec,FF
real(dp)::x1,z1,qT,scaleMu,scaleZeta
real(dp),dimension(1:13)::var
integer,dimension(1:4),intent(in)::process

GlobalCounter=GlobalCounter+1

!!!! KPC formula
if(useKPC) then
    if(TMDF_KPC_IsconvergenceLost()) then
      xSec=1d9
      return
    end if

    call CalculateX1Z1qT(x1,z1,qT,var)

    scaleMu=var(2)
    scaleZeta=var(3)

    !!!!! (Q^2,qT,x1,z1,mu,procc)
    FF=KPC_SIDISconv(var(3),qT,x1,z1,scaleMu*c2_global,process(2:4))
    xSec=PreFactorKPC(var,process,x1,z1,qT)*FF

else !!!! LP formula

  call CalculateX1Z1qT(x1,z1,qT,var)

  !!! setting values of scales
  !!!! If exact scales, zeta*zeta=Q^2-qT^2.
  if(exactScales) then
    scaleMu=sqrt(var(3)-qT**2)
    scaleZeta=var(3)-qT**2
  else
    scaleMu=var(2)
    scaleZeta=var(3)
  end if


  FF=TMDF_F(var(3),qT,x1,z1,scaleMu*c2_global,scaleZeta,scaleZeta,process(2:4))
  xSec=PreFactor2(var,process,x1,z1,qT)*FF

end if

!write(*,*) "{",var(3),",",x1,"},"!,z1
end function xSec
  
    
!---------------------------------INTEGRATED over Z (and Q and X)--------------------------------------------------------------

!!! the variable doZ check should the integration over Z be performed
!!! if doZ=true, the integration is done
!!! if doZ=false the single value (at xMin) is returned
function Xsec_Zint(var,process,doZ,zMin_in,zMax_in)
real(dp),dimension(1:13):: var
logical,intent(in)::doZ
real(dp),intent(in) :: zMin_in,zMax_in
integer,dimension(1:4),intent(in)::process
real(dp) :: Xsec_Zint

real(dp)::zMin,zMax

!!!!-- check th input z-values
if(zmax_in > 1d0) then
    call Warning_Raise('upper limit of z-integration is >1. It is set to 1.',messageCounter,messageTrigger,moduleName)
    zMax=1d0
  else
    zMax=zMax_in
  end if
  if(zmin_in < 0.000001d0) then
    write(*,*) ErrorString('lower limit of z-integration is < 10^{-6}. Evaluation stop.',moduleName)
    stop
  else
    zMin=zMin_in
  end if

!!!! if order inverted return 0
if(zMin>zMax) then
  Xsec_Zint=0d0
  return
end if

!! the integration over Z is required
if(doZ) then
#if INTEGRATION_MODE==1
  !!!! slower but accurate
  Xsec_Zint=Integrate_SA(integrandOverZ,zMin,zMax,toleranceINT)
#elif INTEGRATION_MODE==2
  !!!! fast but not that accurate
  Xsec_Zint=Integrate_G7(integrandOverZ,zMin,zMax)
#endif

else
  ! no integration over Z
  ! evaluation at the central point

  call SetZ((zMin+zMax)/2d0,var)
  Xsec_Zint=Xsec(var,process)
end if

contains

function integrandOverZ(z)
real(dp),intent(in)::z
real(dp)::integrandOverZ
call SetZ(z,var)
integrandOverZ=xSec(var,process)
end function integrandOverZ

end function Xsec_Zint
  
!---------------------------------INTEGRATED over X (and Z)---------------------------------------------------------------

!!!
!!! the variable doX check should the integration over x be performed
!!! if doX=true, the integration is done
!!! if doX=false the single value (at xMin) is returned
function Xsec_Zint_Xint(var,process,doZ,zMin,zMax,doX,Xmin_in,Xmax_in,doCut,Cuts)
real(dp),dimension(1:13) :: var
logical,intent(in)::doX,doZ,doCut
real(dp),dimension(1:4),intent(in)::Cuts
real(dp),intent(in):: xmin_in,xmax_in,zMin,zMax
integer,dimension(1:4)::process
real(dp) :: Xsec_Zint_Xint

real(dp)::xmin, xmax

!!! ----- check the input variables
!!! in the case process=3 the input is y, which is to be transformed to X
!!! evaluate corresponding y's
if(process(1)==3) then
  xmin=XfromSYQ2(var(11),xmin_in,var(3))
  xmax=XfromSYQ2(var(11),xmax_in,var(3))

  !!!! This is important, because the integral is done over dX (and process 3 is over dy)
  process(1)=1
else
  xmin=xmin_in
  xmax=xmax_in
end if

if(xmax > 1) then
  call Warning_Raise('upper limit of x-integration is >1. It is set to 1.',messageCounter,messageTrigger,moduleName)
  xmax=1d0
end if
if(xmin < 0.000001d0) then
  write(*,*) ErrorString('lower limit of x-integration is < 10^{-6}. Evaluation stop.',moduleName)
  stop
end if

!! in case of cut we determine the recut values
if(doCut) then
  xmin=xMinWithCuts(xmin,var,Cuts)
  xmax=xMaxWithCuts(xmax,var,Cuts)
end if

!!! if inverted limits return zero
if(xmin>xmax) then
  Xsec_Zint_Xint=0._dp
  return
end if

if(doX) then
!!!Integration is required

#if INTEGRATION_MODE==1
    !!!! slower but accurate
    Xsec_Zint_Xint=Integrate_SA(integrandOverX,xMin,xMax,toleranceINT)
#elif INTEGRATION_MODE==2
    !!!! fast but not that accurate
    Xsec_Zint_Xint=Integrate_G7(integrandOverX,xMin,xMax)
#endif

else
  call SetX((xMin+xMax)/2d0,var)
  Xsec_Zint_Xint=Xsec_Zint(var,process,doZ,zMin,zMax)
end if

contains

function integrandOverX(x)
real(dp),intent(in)::x
real(dp)::integrandOverX
call SetX(x,var)
integrandOverX=Xsec_Zint(var,process,doZ,zMin,zMax)
end function integrandOverX

end function Xsec_Zint_Xint
  

!---------------------------------INTEGRATED over Q (and X and Z)--------------------------------------------------------------

!!! the variable doQ check should the integration over x be performed
!!! if doQ=true, the integration is done
!!! if doQ=false the single value (at xMin) is returned
function Xsec_Zint_Xint_Qint(var,process,doZ,zMin,zMax,doX,xMin,xMax,doQ,Qmin_in,Qmax_in,doCut,Cuts)
real(dp),dimension(1:13):: var
logical,intent(in)::doX,doQ,doCut,doZ
real(dp),dimension(1:4),intent(in)::Cuts
real(dp),intent(in) :: Qmin_in,Qmax_in,xMin,xMax,zMin,zMax
integer,dimension(1:4)::process
real(dp) :: Xsec_Zint_Xint_Qint

real(dp) :: Qmin, Qmax

!!! check input parameters
!!! evaluate correspnding y's
!!! in the case process=2 the integral is over y
if(process(1)==2) then
  Qmin=QfromSXY(var(11),var(4), Qmin_in)
  Qmax=QfromSXY(var(11),var(4), Qmax_in)

  !!!! This is important, because the integral is done over dQ (and process 2 is over dy)
  process(1)=1
else
  Qmin=Qmin_in
  Qmax=Qmax_in
end if

!! in case of cut we determine recut values
if(doCut) then
  if(doX) then
    Qmin=QMinWithCuts(xmin,Qmin,var,Cuts)
    Qmax=QMaxWithCuts(xmax,Qmax,var,Cuts)
  else
    Qmin=QMinWithCuts((xmin+xmax)/2d0,Qmin,var,Cuts)
    Qmax=QMaxWithCuts((xmin+xmax)/2d0,Qmin,var,Cuts)+toleranceGEN!!! this is needed to resolve 0
  end if
end if

!!!! if limits inverted return 0
if(Qmin>Qmax) then
  Xsec_Zint_Xint_Qint=0._dp
  return
end if


!! the integration over Q is required
if(doQ) then
#if INTEGRATION_MODE==1
    !!!! slower but accurate
    Xsec_Zint_Xint_Qint=Integrate_SA(integrandOverQ,Qmin,Qmax,toleranceINT)
#elif INTEGRATION_MODE==2
    !!!! fast but not that accurate
    Xsec_Zint_Xint_Qint=Integrate_G7(integrandOverQ,Qmin,Qmax)
#endif

else
  call SetQ((Qmin+Qmax)/2d0,var)
  Xsec_Zint_Xint_Qint=Xsec_Zint_Xint(var,process,doZ,zMin,zMax,doX,xmin,xmax,doCut,Cuts)
end if

contains

function integrandOverQ(Q)
real(dp),intent(in)::Q
real(dp)::integrandOverQ
call SetQ(Q,var)
integrandOverQ=2*Q*Xsec_Zint_Xint(var,process,doZ,zMin,zMax,doX,xmin,xmax,doCut,Cuts)
end function integrandOverQ

end function Xsec_Zint_Xint_Qint
  
  

!---------------------------------INTEGRATED over pT (and Z and Q and X)--------------------------------------------------------------

!!! function determines the best value of PT-sections from PT-bin size, and Q
!!! it is determined by formula Q/PT< val/ (2 k) => def+2K
function NumPT_auto(dPT,Q)
real(dp),parameter::val=40d0
real(dp),intent(in)::dPT,Q
integer::NumPT_auto

real(dp)::rat
integer::i

rat=Q/dPT

if(rat>40d0) then
  NumPT_auto=NumPTdefault
  return
else
  do i=1,5
    if(rat>(40d0/2d0/i)) then
      NumPT_auto=NumPTdefault+2*i
      return
    end if
  end do
end if
if(outputlevel>1) then
write(*,*) WarningString('Fail to automatically determine number of Pt-section for a bin.',moduleName)
write(*,*) '>>  Possibly Pt-bin is too large', dPT
end if
NumPT_auto=NumPTdefault+12

end function NumPT_auto
  
!!! the variable doPT check should the integration over pT be performed
!!! if doPT=true, the integration is done
!!! if doPT=false the single value (at xMin) is returned
function Xsec_Zint_Xint_Qint_PTint(var,process,doZ,zMin,zMax,doX,xMin,xMax,doQ,Qmin,Qmax,doPT,ptMin_in,ptMax_in,doCut,Cuts,Num)
real(dp),dimension(1:13) :: var
real(dp),dimension(1:4),intent(in) :: Cuts
logical,intent(in)::doX,doQ,doZ,doPT,doCut
integer,dimension(1:4),intent(in)::process
real(dp),intent(in) :: Qmin,Qmax,xMin,xMax,zMin,zMax,ptMax_in,ptMin_in
integer,intent(in)::Num
real(dp) :: Xsec_Zint_Xint_Qint_PTint

real(dp) :: ptMin,ptMax

if(ptMin_in<ptMIN_global) then
  pTmin=ptMIN_global
else
  pTmin=ptMin_in
end if

ptMax=ptMax_in

!!! inverted limits = 0
if(ptMin>ptMax) then
  Xsec_Zint_Xint_Qint_PTint=0._dp
  return
end if

!! the integration over PT is required
if(doPT) then
  Xsec_Zint_Xint_Qint_PTint=Integrate_SN(integrandOverPT,ptMin,ptMax,Num)

else
  ! no integration over PT
  call SetPT((PTmin+PTmax)/2d0,var)
  Xsec_Zint_Xint_Qint_PTint=Xsec_Zint_Xint_Qint(var,process,doZ,zMin,zMax,doX,xMin,xMax,doQ,Qmin,Qmax,doCut,Cuts)

end if

contains

function integrandOverPT(pT)
real(dp),intent(in)::pT
real(dp)::integrandOverPT
call SetPT(pT,var)
integrandOverPT=2*pT*Xsec_Zint_Xint_Qint(var,process,doZ,zMin,zMax,doX,xMin,xMax,doQ,Qmin,Qmax,doCut,Cuts)
end function integrandOverPT

end function Xsec_Zint_Xint_Qint_PTint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!  MAIN INTERFACE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! single value interface
subroutine xSec_SIDIS(xx,process,s,pT,z,x,Q,doCut,Cuts,masses)
  integer,intent(in),dimension(1:4)::process            !the number of process
  real(dp),intent(in)::s                                !Mandelshtam s
  real(dp),intent(in),dimension(1:2)::pT                !(qtMin,qtMax)
  real(dp),intent(in),dimension(1:2)::z                 !(zmin,zmax)
  real(dp),intent(in),dimension(1:2)::x                 !(xmin,xmax)
  real(dp),intent(in),dimension(1:2)::Q                 !(Qmin,Qmax)
  logical,intent(in)::doCut                             !triger cuts
  real(dp),intent(in),dimension(1:4)::Cuts              !(ymin,yMax,W2)
  real(dp),intent(in),dimension(1:2),optional::masses   !(mass_target,mass-product)GeV
  real(dp),intent(out)::xx

  real(dp),dimension(1:13):: var
  integer::Num

if(.not.started) ERROR STOP ErrorString('The module is not initialized. Check INI-file.',moduleName)

  CallCounter=CallCounter+1
  if(PRESENT(masses)) then
  var=kinematicArray((pt(1)+pt(2))/2d0,s,(z(1)+z(2))/2d0,(x(1)+x(2))/2d0,(Q(1)+Q(2))/2d0,masses(1)**2,masses(2)**2)
  else
  var=kinematicArray((pt(1)+pt(2))/2d0,s,(z(1)+z(2))/2d0,(x(1)+x(2))/2d0,(Q(1)+Q(2))/2d0,0._dp,0._dp)
  end if


  Num=NumPT_auto(pt(2)-pt(1),var(2))

  xx=Xsec_Zint_Xint_Qint_PTint(var,process,.true.,z(1),z(2),.true.,x(1),x(2),.true.,Q(1),Q(2),.true.,pt(1),pt(2),doCut,Cuts,Num)


end subroutine xSec_SIDIS

subroutine xSec_SIDIS_List(xx,process,s,pT,z,x,Q,doCut,Cuts,masses)
  integer,intent(in),dimension(:,:)::process            !the number of process
  real(dp),intent(in),dimension(:)::s                !Mandelshtam s
  real(dp),intent(in),dimension(:,:)::pT            !(qtMin,qtMax)
  real(dp),intent(in),dimension(:,:)::z                !(zmin,zmax)
  real(dp),intent(in),dimension(:,:)::x                !(xmin,xmax)
  real(dp),intent(in),dimension(:,:)::Q                !(Qmin,Qmax)
  logical,intent(in),dimension(:)::doCut            !triger cuts
  real(dp),intent(in),dimension(:,:)::Cuts            !(ymin,yMax,W2)
  real(dp),intent(in),dimension(:,:),optional::masses        !(mass_target,mass-product)GeV
  real(dp),dimension(:),intent(out)::xx
  integer :: i,length

if(.not.started) ERROR STOP ErrorString('The module is not initialized. Check INI-file.',moduleName)

  length=size(s)
  CallCounter=CallCounter+length

  !!! cheking sizes
  if(size(xx)/=length) then
    write(*,*) ErrorString('xSec_SIDIS_List: sizes of xSec and s lists are not equal.',moduleName)
    write(*,*) 'Evaluation stop'
    stop
  end if
  if(size(process,1)/=length) then
    write(*,*) ErrorString('xSec_SIDIS_List: sizes of process and s lists are not equal.',moduleName)
    write(*,*) 'Evaluation stop'
    stop
  end if
  if(size(pT,1)/=length) then
    write(*,*) ErrorString('xSec_SIDIS_List: sizes of pT and s lists are not equal.',moduleName)
    write(*,*) 'Evaluation stop'
    stop
  end if
  if(size(x,1)/=length) then
    write(*,*) ErrorString('xSec_SIDIS_List: sizes of x and s lists are not equal.',moduleName)
    write(*,*) 'Evaluation stop'
    stop
  end if
  if(size(Q,1)/=length) then
    write(*,*) ErrorString('xSec_SIDIS_List: sizes of Q and s lists are not equal.',moduleName)
    write(*,*) 'Evaluation stop'
    stop
  end if
  if(size(z,1)/=length) then
    write(*,*) ErrorString('xSec_SIDIS_List: sizes of z and s lists are not equal.',moduleName)
    write(*,*) 'Evaluation stop'
    stop
  end if
  if(size(doCut)/=length) then
    write(*,*) ErrorString('xSec_SIDIS_List: sizes of doCut and s lists are not equal.',moduleName)
    write(*,*) 'Evaluation stop'
    stop
  end if
  if(size(Cuts,1)/=length) then
    write(*,*) ErrorString('xSec_SIDIS_List: sizes of Cuts and s lists are not equal.',moduleName)
    write(*,*) 'Evaluation stop'
    stop
  end if
  if(size(process,2)/=4) then
    write(*,*) ErrorString('xSec_SIDIS_List: process list must be (:,1:3).',moduleName)
    write(*,*) 'Evaluation stop'
    stop
  end if
  if(size(pT,2)/=2) then
    write(*,*) ErrorString('xSec_SIDIS_List: pt list must be (:,1:2).',moduleName)
    write(*,*) 'Evaluation stop'
    stop
  end if
  if(size(x,2)/=2) then
    write(*,*) ErrorString('xSec_SIDIS_List: x list must be (:,1:2).',moduleName)
    write(*,*) 'Evaluation stop'
    stop
  end if
  if(size(Q,2)/=2) then
    write(*,*) ErrorString('xSec_SIDIS_List: Q list must be (:,1:2).',moduleName)
    write(*,*) 'Evaluation stop'
    stop
  end if
  if(size(z,2)/=2) then
    write(*,*) ErrorString('xSec_SIDIS_List: z list must be (:,1:2).',moduleName)
    write(*,*) 'Evaluation stop'
    stop
  end if
  if(size(Cuts,2)/=4) then
    write(*,*) ErrorString('xSec_SIDIS_List: cuts list must be (:,1:4).',moduleName)
    write(*,*) 'Evaluation stop'
    stop
  end if

  CallCounter=CallCounter+length
  if(PRESENT(masses)) then

  if(size(masses,1)/=length) then
    write(*,*) ErrorString('xSec_SIDIS_List: sizes of masses and s lists are not equal.',moduleName)
    write(*,*) 'Evaluation stop'
    stop
  end if
      if(size(masses,2)/=2) then
    write(*,*) ErrorString('xSec_SIDIS_List: mass list must be (:,1:2).',moduleName)
    write(*,*) 'Evaluation stop'
    stop
  end if

  !$OMP PARALLEL DO DEFAULT(SHARED)
  do i=1,length
  xx(i)=xSecFULL(process(i,1:4),s(i),pt(i,1),pt(i,2),z(i,1),z(i,2),x(i,1),x(i,2),Q(i,1),Q(i,2),doCut(i),Cuts(i,1:4),&
          masses(i,1)**2,masses(i,2)**2)
  end do
  !$OMP END PARALLEL DO

  else

  !$OMP PARALLEL DO DEFAULT(SHARED)
  do i=1,length
  xx(i)=xSecFULL(process(i,1:4),s(i),pt(i,1),pt(i,2),z(i,1),z(i,2),x(i,1),x(i,2),Q(i,1),Q(i,2),doCut(i),Cuts(i,1:4),&
          0._dp,0._dp)
  end do
  !$OMP END PARALLEL DO

  end if

end subroutine xSec_SIDIS_List


!!!! problem is that f2py does not like optional arguments.. in any form
subroutine xSec_SIDIS_List_forharpy(xx,process,s,pT,z,x,Q,doCut,Cuts,masses)
  integer,intent(in),dimension(:,:)::process            !the number of process
  real(dp),intent(in),dimension(:)::s                !Mandelshtam s
  real(dp),intent(in),dimension(:,:)::pT            !(qtMin,qtMax)
  real(dp),intent(in),dimension(:,:)::z                !(zmin,zmax)
  real(dp),intent(in),dimension(:,:)::x                !(xmin,xmax)
  real(dp),intent(in),dimension(:,:)::Q                !(Qmin,Qmax)
  logical,intent(in),dimension(:)::doCut            !triger cuts
  real(dp),intent(in),dimension(:,:)::Cuts            !(ymin,yMax,W2)
  real(dp),intent(in),dimension(:,:)::masses        !(mass_target,mass-product)GeV
  real(dp),dimension(:),intent(out)::xx
  integer :: i,length

if(.not.started) ERROR STOP ErrorString('The module is not initialized. Check INI-file.',moduleName)

  length=size(s)
  CallCounter=CallCounter+length

  !!! cheking sizes
  if(size(xx)/=length) then
    ERROR STOP ErrorString('xSec_SIDIS_List: sizes of xSec and s lists are not equal.',moduleName)
  end if
  if(size(process,1)/=length) then
    ERROR STOP ErrorString('xSec_SIDIS_List: sizes of process and s lists are not equal.',moduleName)
  end if
  if(size(pT,1)/=length) then
    ERROR STOP ErrorString('xSec_SIDIS_List: sizes of pT and s lists are not equal.',moduleName)
  end if
  if(size(x,1)/=length) then
    ERROR STOP ErrorString('xSec_SIDIS_List: sizes of x and s lists are not equal.',moduleName)
  end if
  if(size(Q,1)/=length) then
    ERROR STOP ErrorString('xSec_SIDIS_List: sizes of Q and s lists are not equal.',moduleName)
  end if
  if(size(z,1)/=length) then
    ERROR STOP ErrorString('xSec_SIDIS_List: sizes of z and s lists are not equal.',moduleName)
  end if
  if(size(doCut)/=length) then
    ERROR STOP ErrorString('xSec_SIDIS_List: sizes of doCut and s lists are not equal.',moduleName)
  end if
  if(size(Cuts,1)/=length) then
    ERROR STOP ErrorString('xSec_SIDIS_List: sizes of Cuts and s lists are not equal.',moduleName)
  end if
  if(size(process,2)/=4) then
    ERROR STOP ErrorString('xSec_SIDIS_List: process list must be (:,1:3).',moduleName)
  end if
  if(size(pT,2)/=2) then
    ERROR STOP ErrorString('xSec_SIDIS_List: pt list must be (:,1:2).',moduleName)
  end if
  if(size(x,2)/=2) then
    ERROR STOP ErrorString('xSec_SIDIS_List: x list must be (:,1:2).',moduleName)
  end if
  if(size(Q,2)/=2) then
    ERROR STOP ErrorString('xSec_SIDIS_List: Q list must be (:,1:2).',moduleName)
  end if
  if(size(z,2)/=2) then
    ERROR STOP ErrorString('xSec_SIDIS_List: z list must be (:,1:2).',moduleName)
  end if
  if(size(Cuts,2)/=4) then
    ERROR STOP ErrorString('xSec_SIDIS_List: cuts list must be (:,1:4).',moduleName)
  end if

  CallCounter=CallCounter+length

  if(size(masses,1)/=length) then
    ERROR STOP ErrorString('xSec_SIDIS_List: sizes of masses and s lists are not equal.',moduleName)
  end if
  if(size(masses,2)/=2) then
    ERROR STOP ErrorString('xSec_SIDIS_List: mass list must be (:,1:2).',moduleName)
  end if

  !$OMP PARALLEL DO DEFAULT(SHARED)
  do i=1,length

  xx(i)=xSecFULL(process(i,1:4),s(i),pt(i,1),pt(i,2),z(i,1),z(i,2),x(i,1),x(i,2),Q(i,1),Q(i,2),doCut(i),Cuts(i,1:4),&
          masses(i,1)**2,masses(i,2)**2)
  end do
  !$OMP END PARALLEL DO
end subroutine xSec_SIDIS_List_forharpy


!!! helper to incapsulate PARALLEL variables
function xSecFULL(proc,s,ptmin,ptmax,zmin,zmax,xmin,xmax,Qmin,Qmax,doCut,Cuts,m1,m2)
real(dp),intent(in)::s,ptmin,ptmax,zmin,zmax,xmin,xmax,Qmin,Qmax,Cuts(1:4),m1,m2
integer,intent(in)::proc(1:4)
logical,intent(in)::doCut
real(dp):: xSecFULL

real(dp)::var(1:13)
integer::Num

if(.not.started) ERROR STOP ErrorString('The module is not initialized. Check INI-file.',moduleName)

var=kinematicArray((ptmin+ptmax)/2d0,s,(zmin+zmax)/2d0,(xmin+xmax)/2d0,(Qmin+Qmax)/2d0,m1,m2)
Num=NumPT_auto(ptmax-ptmin,var(2))

!   write(*,*) 'aTMD:1  ',var
!   write(*,*) 'aTMD:2  ',proc,m1,m2
!   write(*,*) 'aTMD:3  ',.true.,zmin,zmax,.true.,xmin,xmax,.true.,Qmin,Qmax,.true.,ptmin,ptmax

xSecFULL=Xsec_Zint_Xint_Qint_PTint(var,proc,&
          .true.,zmin,zmax,.true.,xmin,xmax,.true.,Qmin,Qmax,.true.,ptmin,ptmax,doCut,Cuts,Num)
end function xSecFULL

!!!! problem is that f2py does not like optional arguments.. in any form
subroutine xSec_SIDIS_BINLESS_List_forharpy(xx,process,s,pT,z,x,Q,masses)
  integer,intent(in),dimension(:,:)::process            !the number of process
  real(dp),intent(in),dimension(:)::s                   !Mandelshtam s
  real(dp),intent(in),dimension(:)::pT                  !(pt)
  real(dp),intent(in),dimension(:)::z                   !(z)
  real(dp),intent(in),dimension(:)::x                   !(x)
  real(dp),intent(in),dimension(:)::Q                   !(Q)
  real(dp),intent(in),dimension(:,:)::masses            !(mass_target,mass-product)GeV
  real(dp),dimension(:),intent(out)::xx
  integer :: i,length

  if(.not.started) ERROR STOP ErrorString('The module is not initialized. Check INI-file.',moduleName)

  length=size(s)
  CallCounter=CallCounter+length

  !!! cheking sizes
  if(size(xx)/=length) then
    ERROR STOP ErrorString('xSec_SIDIS_BINLESS_List: sizes of xSec and s lists are not equal.',moduleName)
  end if
  if(size(process,1)/=length) then
    ERROR STOP ErrorString('xSec_SIDIS_BINLESS_List: sizes of process and s lists are not equal.',moduleName)
  end if
  if(size(pT)/=length) then
    ERROR STOP ErrorString('xSec_SIDIS_BINLESS_List: sizes of pT and s lists are not equal.',moduleName)
  end if
  if(size(x)/=length) then
    ERROR STOP ErrorString('xSec_SIDIS_BINLESS_List: sizes of x and s lists are not equal.',moduleName)
  end if
  if(size(Q)/=length) then
    ERROR STOP ErrorString('xSec_SIDIS_BINLESS_List: sizes of Q and s lists are not equal.',moduleName)
  end if
  if(size(z)/=length) then
    ERROR STOP ErrorString('xSec_SIDIS_BINLESS_List: sizes of z and s lists are not equal.',moduleName)
  end if
  if(size(process,2)/=4) then
    ERROR STOP ErrorString('xSec_SIDIS_BINLESS_List: process list must be (:,1:3).',moduleName)
  end if

  CallCounter=CallCounter+length

  if(size(masses,1)/=length) then
    ERROR STOP ErrorString('xSec_SIDIS_BINLESS_List: sizes of masses and s lists are not equal.',moduleName)
  end if
  if(size(masses,2)/=2) then
    ERROR STOP ErrorString('xSec_SIDIS_BINLESS_List: mass list must be (:,1:2).',moduleName)
  end if

  !$OMP PARALLEL DO DEFAULT(SHARED)
  do i=1,length
  xx(i)=xSec_SIDIS_BINLESS(process(i,1:4),s(i),pt(i),z(i),x(i),Q(i),masses(i,1)**2,masses(i,2)**2)
  end do
  !$OMP END PARALLEL DO
end subroutine xSec_SIDIS_BINLESS_List_forharpy


!!! helper to incapsulate PARALLEL variables
!!! evaluate cross-section without integrations over bins, and without cuts
function xSec_SIDIS_BINLESS(proc,s,pt,z,x,Q,m1,m2)
real(dp),intent(in)::s,pt,z,x,Q,m1,m2
integer,intent(in)::proc(1:4)
real(dp)::xSec_SIDIS_BINLESS

real(dp)::var(1:13)

if(.not.started) ERROR STOP ErrorString('The module is not initialized. Check INI-file.',moduleName)

var=kinematicArray(pt,s,z,x,Q,m1,m2)

!    write(*,*) 'aTMD:1  ',var
!    write(*,*) 'aTMD:2  ',proc,m1,m2
!    write(*,*) 'aTMD:3  ',.false.,z,z,.false.,x,x,.false.,Q,Q,.false.,pt,pt

xSec_SIDIS_BINLESS=Xsec_Zint_Xint_Qint_PTint(var,proc,&
          .false.,z,z+toleranceGEN,&
          .false.,x,x+toleranceGEN,&
          .false.,Q,Q+toleranceGEN,&
          .false.,pt,pt+toleranceGEN,&
          .false.,(/0d0,1d0,0d0,1d9/),4)
end function xSec_SIDIS_BINLESS

end module TMDX_SIDIS
