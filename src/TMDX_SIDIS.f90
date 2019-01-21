!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 1.31
!
!	Evaluation of the TMD cross-section for SIDIS-like cross-sections
!	
!	if you use this module please, quote ????.????
!
!	ver 1.2: release (AV, 15.12.2017)
!	ver 1.32: part of functions migrated to TMDF, rest updated (AV, 16.08.2018)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDX_SIDIS

use TMDF
use QCDinput
use EWinput

implicit none
  private
  
   !Current version of module
 character (len=5),parameter :: version="v1.32"
  
  real*8 :: tolerance=0.0005d0
  
  !!! Target mass of target
  real*8:: M1=0.939d0
  !!! Produced mass
  real*8:: M2=0.1d0
  
  integer::outputlevel
  
  !The variables for all parameters of the model!
  real*8:: Q_global,Q2_global,x_global,y_global,z_global,zeta_global,s_global,&
	  gamma_global,rho_global,N2_global,varepsilon_global
  real*8::prefactor2_global
  !! Set of process definition, for Prefactor 1, Prefactor 2, structure function
  integer:: processP1_global,processP2_global,processP3_global
  !!other global parameters see SetXParameters  
  integer:: orderH_global
  logical:: exactX1X2  !!! expressions for x1 x2 include qT
  logical:: AccountM1,AccountM2,AccountPT !!! account masses in calculation of kinematical variables.
  
  !!! number of sections for PT-integral by default
  integer::NumPTdefault=4
  
  real*8::c2_global,muHard_global
  
  integer::GlobalCounter,CallCounter
  
  real*8::hc2
  
  logical::started=.false.
  
  public::TMDX_SIDIS_XSetup,TMDX_SIDIS_SetTargetMass,TMDX_SIDIS_SetProducedMass,&
    TMDX_SIDIS_SetNPParameters,TMDX_SIDIS_Initialize,TMDX_SIDIS_SetScaleVariations,&
    TMDX_SIDIS_setProcess,TMDX_SIDIS_ShowStatistic
  public::  CalcXsec_SIDIS,CalcXsec_SIDIS_Xint,CalcXsec_SIDIS_Zint,CalcXsec_SIDIS_Qint
  
 interface TMDX_SIDIS_SetNPParameters
  module procedure TMDX_SIDIS_SetNPParameters,TMDX_SIDIS_SetNPParameters_rep
 end interface
 
 interface CalcXsec_SIDIS
    module procedure xSecSingle,xSecList
  end interface
  
 interface TMDX_SIDIS_setProcess
    module procedure TMDX_setProcess1,TMDX_setProcess3,TMDX_setProcess30
 end interface
 
 interface CalcXsec_SIDIS_Xint
    module procedure xSecSingle_Xint ,xSecList_Xint
 end interface
 
 interface CalcXsec_SIDIS_Zint
    module procedure xSecSingle_Zint ,xSecList_Zint
 end interface
 
  interface CalcXsec_SIDIS_Qint
    module procedure xSecSingle_Qint ,xSecList_Qint
 end interface
 
contains
  
  !!Just passes the initialization to subpackages
  !! This also set orders. Orders cannot be changes afterwards
  subroutine TMDX_SIDIS_Initialize(orderMain)
  character(len=*)::orderMain
  character(256)::line
  real*8::dummy
  
    OPEN(UNIT=51, FILE='constants', ACTION="read", STATUS="old")    
    !!! Search for output level
    do
    read(51,'(A)') line    
    if(line(1:3)=='*0 ') exit
    end do    
    do
    read(51,'(A)') line
    if(line(1:3)=='*A ') exit
    end do
    read(51,'(A)') line
    read(51,*) outputLevel
  
  processP2_global=-10221191
  
  if(outputLevel>1) write(*,*) '----- arTeMiDe.TMD_SIDIS ',version,': .... initialization'
     SELECT CASE(orderMain)
      CASE ("LO")
	orderH_global=0
      CASE ("LO+")
	orderH_global=0
      CASE ("NLO")
	orderH_global=1
      CASE ("NLO+")
	orderH_global=1
      CASE ("NNLO")
	orderH_global=2
      CASE ("NNLO+")
	orderH_global=3
      CASE DEFAULT
	if(outputLevel>0) write(*,*) 'WARNING arTeMiDe.TMDX_SIDIS:try to set unknown order. Switch to NLO.'
	orderH_global=1
     END SELECT
     
    !!!! Physical constants
    do
    read(51,'(A)') line    
    if(line(1:3)=='*1 ') exit
    end do    
    do
    read(51,'(A)') line
    if(line(1:3)=='*C ') exit
    end do
    read(51,'(A)') line
    read(51,*) hc2    !!!!!!!!!!!GeV->mbarn

    
    do
    read(51,'(A)') line    
    if(line(1:3)=='*2 ') exit
    end do    
    
    do
    read(51,'(A)') line
    if(line(1:3)=='*A ') exit
    end do
    read(51,'(A)') line
    read(51,*) tolerance
    read(51,'(A)') line
    read(51,*) NumPTdefault
    
    do
    read(51,'(A)') line
    if(line(1:3)=='*5 ') exit
    end do
    do
    read(51,'(A)') line
    if(line(1:3)=='*A ') exit
    end do
    do
    read(51,'(A)') line
    if(line(1:3)=='*2)') exit
    end do
    read(51,*) exactX1X2
    
    do
    read(51,'(A)') line
    if(line(1:3)=='*B ') exit
    end do
    do
    read(51,'(A)') line
    if(line(1:3)=='*1)') exit
    end do
    read(51,*) AccountM1
    do
    read(51,'(A)') line
    if(line(1:3)=='*2)') exit
    end do
    read(51,*) AccountM2
    do
    read(51,'(A)') line
    if(line(1:3)=='*3)') exit
    end do
    read(51,*) AccountPT
    
    CLOSE (51, STATUS='KEEP')
     
     call EWinput_Initialize(orderMain)
     call TMDF_Initialize(orderMain)
     
     call TMDX_SIDIS_SetTargetMass(0.938d0)
     call TMDX_SIDIS_SetProducedMass(0.12d0)
     
     c2_global=1d0
     
     GlobalCounter=0
     CallCounter=0
     
     started=.true.
    write(*,*)  '----- arTeMiDe.TMD_SIDIS ',version,'.... initialized'
  end subroutine TMDX_SIDIS_Initialize
  
  !!!!Call this after TMD initializetion but before NP, and X parameters
  subroutine TMDX_SIDIS_SetScaleVariations(c1_in,c2_in,c3_in,c4_in)
    real*8::c1_in,c2_in,c3_in,c4_in
    
    if(outputLevel>1) write(*,*) 'TMDX_SIDIS: scales reset:',c1_in,c2_in,c3_in,c4_in
    
    call TMDF_SetScaleVariations(c1_in,c3_in,c4_in)
    
    if(c2_in<0.1d0 .or. c2_in>10.d0) then
    if(outputLevel>0) write(*,*) 'TMDX_SIDIS WARNING: variation in c2 is enourmous. c2 is set to 2'
     c2_global=2d0
    else
    c2_global=c2_in
    end if
    
    muHard_global=Q_global*c2_global
    call EvaluatePreFactor2_qTindependent()
    
  end subroutine TMDX_SIDIS_SetScaleVariations
  
  !!Just passes settting of NP parameters to subpackage
  subroutine TMDX_SIDIS_SetNPParameters(lambda)
    real*8::lambda(:)
    GlobalCounter=0
     CallCounter=0
    call TMDF_SetNPParameters(lambda)
  end subroutine TMDX_SIDIS_SetNPParameters
  
    !!Just passes settting of NP parameters to subpackage
  subroutine TMDX_SIDIS_SetNPParameters_rep(num)
    integer::num
    GlobalCounter=0
     CallCounter=0
    call TMDF_SetNPParameters(num)
  end subroutine TMDX_SIDIS_SetNPParameters_rep
  
  !!!set variables for process definition
  subroutine TMDX_setProcess30(p0)
  integer,dimension(1:3)::p0
  
  processP1_global=p0(1)
  processP2_global=p0(2)
  processP3_global=p0(3)
  end subroutine TMDX_setProcess30
  
  !!!set variables for process definition
  subroutine TMDX_setProcess3(p1,p2,p3)
  integer::p1,p2,p3
  
  processP1_global=p1
  processP2_global=p2
  processP3_global=p3
  end subroutine TMDX_setProcess3
  
  !!!set variables for process definition
  subroutine TMDX_setProcess1(p)
  integer::p
  
  SELECT CASE(p)
   case(1)
    call TMDX_setProcess3(1,1,5) !!! p + p -> Z + gamma^*   (e.g. ATLAS, CMS, LHCb)
   case(2)
    call TMDX_setProcess3(1,1,6) !!! p + pbar -> Z + gamma^*   (e.g. CDF,D0)
   case(4)
    call TMDX_setProcess3(1,2,1001)
   case(5)
    call TMDX_setProcess3(1,1,5)
   case(7)
    call TMDX_setProcess3(1,1,6)
   case default
    write(*,*) 'ERROR: arTeMiDe_SIDIS: unknown process is called. p=',p
    write(*,*) 'Evaluation stop'
    stop
   end SELECT
  end subroutine TMDX_setProcess1
  
  !!! Set the mass of target hadron M1
  subroutine TMDX_SIDIS_SetTargetMass(M)
  real*8::M
  if(AccountM1) then
   M1=M
  else
   M1=0d0
  end if
  end subroutine TMDX_SIDIS_SetTargetMass
  
  !!! Set the mass of produced hadron M2
  subroutine TMDX_SIDIS_SetProducedMass(M)
  real*8::M
  if(AccountM2) then
   M2=M
  else
   M2=0d0
  end if
  end subroutine TMDX_SIDIS_SetProducedMass
  
  subroutine TMDX_SIDIS_XSetup(s,Q,x,z)
   real*8::s,Q,x,z
    
    if(.not.started) then
    write(*,*) 'ERROR: arTeMiDe.TMDX_SIDIS is not initialized. Evaluation terminated'
    stop
    end if
    
    s_global=s
    Q_global=Q
    x_global=x    
    z_global=z
    
    Q2_global=Q**2
    y_global=Q2_global/x/s
    
    !for a moment we fix like that
    zeta_global=Q2_global
    
    gamma_global=2*M1*x/Q
    
    rho_global=M2/z/Q
    varepsilon_global=(1-y_global-(gamma_global*y_global/2d0)**2)&
	    /(1-y_global+y_global**2/2+(gamma_global*y_global/2d0)**2)
    N2_global=SQRT((1+gamma_global**2)/(1-(gamma_global*rho_global)**2))
    
    muHard_global=Q_global*c2_global
   
    call EvaluatePreFactor2_qTindependent()
    
  end subroutine TMDX_SIDIS_XSetup
  
  !!!intrinsic change the value of Q (and only related stuff)
  subroutine TMDX_SIDIS_SetQ(Q)
    real*8::Q
   
    Q_global=Q
    Q2_global=Q**2
    zeta_global=Q2_global
    
    gamma_global=2*M1*x_global/Q
    
    rho_global=M2/z_global/Q
    varepsilon_global=(1-y_global-(gamma_global*y_global/2d0)**2)&
	  /(1-y_global+y_global**2/2+(gamma_global*y_global/2d0)**2)
    N2_global=SQRT((1+gamma_global**2)/(1-(gamma_global*rho_global)**2))
    
    muHard_global=Q_global*c2_global
   
    call EvaluatePreFactor2_qTindependent()
  end subroutine TMDX_SIDIS_SetQ
  
  !!!intrinsic change the value of y (and only related stuff)
  subroutine TMDX_SIDIS_SetX(x)
    real*8::x
   
    x_global=x
    
    gamma_global=2*M1*x/Q_global
    varepsilon_global=(1-y_global-(gamma_global*y_global/2d0)**2)&
	  /(1-y_global+y_global**2/2+(gamma_global*y_global/2d0)**2)
    N2_global=SQRT((1+gamma_global**2)/(1-(gamma_global*rho_global)**2))
    call EvaluatePreFactor2_qTindependent()
  end subroutine TMDX_SIDIS_SetX
  
    !!!intrinsic change the value of y (and only related stuff)
  subroutine TMDX_SIDIS_SetZ(z)
    real*8::z
   
    z_global=z
    
    rho_global=M2/z/Q_global
    N2_global=SQRT((1+gamma_global**2)/(1-(gamma_global*rho_global)**2))
    call EvaluatePreFactor2_qTindependent()
  end subroutine TMDX_SIDIS_SetZ
  
  
  !!! hard coefficeint taken from 1004.3653 up to 2-loop
  !!! it takes global values of Q,order
  function HardCoefficientSIDIS()
    real*8::HardCoefficientSIDIS,LQ!=Log[Q^2/mu^2]=-2Log[c1]
    
    HardCoefficientSIDIS=1.d0
    if(orderH_global>=1) then
      LQ=-2d0*LOG(c2_global)
      HardCoefficientSIDIS=HardCoefficientSIDIS+As(muHard_global)*&
      (-16.946842488404727d0 + 8d0*LQ - 2.6666666666666665d0*LQ**2)
    if(orderH_global>=2) then
      HardCoefficientSIDIS=HardCoefficientSIDIS+As(muHard_global)**2*&
      (-116.50054911601637d0 + 46.190372772820254d0*LQ + 16.843858371984233d0*LQ**2&
	  -13.333333333333334d0*LQ**3 + 3.5555555555555554d0*LQ**4)
    end if  
    end if
  end function HardCoefficientSIDIS
  
    
  subroutine TMDX_SIDIS_ShowStatistic()
      call TMDF_ShowStatistic()
  
      write(*,'(A,ES12.3)') 'TMDX SIDIS statistics   total calls of point xSec  :  ',Real(GlobalCounter)
      write(*,'(A,ES12.3)') '                              total calls of xSecF :  ',Real(CallCounter)
      write(*,'(A,F12.3)')  '                                         avarage M :  ',Real(GlobalCounter)/Real(CallCounter)
  end subroutine TMDX_SIDIS_ShowStatistic
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS FOR PREFACTORS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !!!!! Prefactor 2 is (universal part) x (cuts)
  function PreFactor2(qT)
    real*8::PreFactor2,cutPrefactor,qT,rhoT2
    
    rhoT2=rho_global**2+(qT/z_global/Q_global)**2
  
    cutPrefactor=1+(varepsilon_global-gamma_global**2/2)*(rhoT2-rho_global**2)&
      /(1d0-(gamma_global*rho_global)**2)
      
    if(AccountPT) cutPrefactor=cutPrefactor/sqrt(1+rhoT2*gamma_global**2)
  
    PreFactor2=prefactor2_global*cutPrefactor
  end function Prefactor2

  
  !!!! Set global part of the prefactor 2
  subroutine EvaluatePreFactor2_qTindependent()
  
  SELECT CASE(processP2_global)
    case(-10221191)
	prefactor2_global=1d0
    CASE(1)
	!pi aEM^2/Q^2 y/(1-eps)
	prefactor2_global=3.141592653589793d0*(alphaEM(Q_global)**2)/Q2_global*&
	    y_global/(1-varepsilon_global)*&
	    HardCoefficientSIDIS()*&
	    hc2*1d9!from GeV to pb
    CASE(2)
	!pi aEM^2/Q^4 y^2/(1-eps)
	prefactor2_global=3.141592653589793d0*(alphaEM(Q_global)**2)/(Q2_global)**2*&
	    y_global**2/(1-varepsilon_global)*&
	    HardCoefficientSIDIS()*&
	    hc2*1d9!from GeV to pb
    CASE (3) 
    !pi aEM^2/Q^3 x y/(1-eps)
	prefactor2_global=3.141592653589793d0*(alphaEM(Q_global)**2)/(Q2_global)**2*&
	    x_global*y_global/(1-varepsilon_global)*&
	    HardCoefficientSIDIS()*&
	    hc2*1d9!from GeV to pb
    CASE DEFAULT 
      write(*,*) 'ERROR: arTeMiDe.TMDX_SIDIS: unknown process p2=',processP2_global,' .Evaluation stop.'
      stop
  END SELECT
  end subroutine EvaluatePreFactor2_qTindependent
  
  !!!! Set Prefactor1
  function PreFactor1()
  real*8::Prefactor1
  SELECT CASE(processP1_global)
    CASE(1)
	PreFactor1=1d0
    CASE DEFAULT 
      write(*,*) 'ERROR: arTeMiDe.TMDX_SIDIS: unknown process p1=',processP2_global,' .Evaluation stop.'
      stop
  END SELECT
  end function PreFactor1
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS CALCULATING CROSS-SECTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !---------------------------------UNTIGRATED------------------------------------------------------------------
  
   !!! this is help function which evaluate xSec at single qt (without lists) with only prefactor 2
  function xSec(qT_in)
  real*8:: xSec,qT_in,FF,kappa
  real*8::x1,z1
   GlobalCounter=GlobalCounter+1
   
   kappa=Sqrt(1-qT_in**2/Q2_global)
   
   if(exactX1X2) then
    x1=2*x_global*kappa/(1d0+sqrt(1d0+kappa*gamma_global**2))
    z1=z_global*(1d0+sqrt(1-(rho_global*gamma_global)**2))/(1d0+sqrt(1d0+kappa*gamma_global**2))
   else
    x1=x_global
    z1=z_global
   end if
   
   FF=TMD_F(Q2_global,qT_in*N2_global/z1,x1,z1,muHard_global,zeta_global,zeta_global,processP3_global)
   xSec=Prefactor2(qT_in)*FF
  end function xSec
  
  
  !---------------------------------UNTIGRATED over X---------------------------------------------------------------
  
  !!!
  function Xsec_Xint(qt,xmin_in,xmax_in)
    real*8 :: qt,Xsec_Xint
    real*8 :: xmin_in,xmax_in

   if(xmax_in>1d0 .or. xmin_in<0d0) then
    write(*,*) 'arTeMiDe_SIDIS: CRITICAL ERROR limits of x-integration outside of (0,1)'
    write(*,*) 'Evaluation stop'
    stop
   end if
    
    Xsec_Xint=integralOverXpoint_S(qt,xmin_in,xmax_in)
  end function Xsec_Xint
  !--------------Simpsons--------------------
  !!!! parameter valueMax remembers the approximate value of integral to weight the tolerance.
  !!!! evaluation is done by adaptive simpson
  !!!! First we evaluate over 5 points and estimate the integral, and then split it to 3+3 and send to adaptive
  !!!! Thus minimal number of points =9
  function integralOverXpoint_S(qt,yMin_in,yMax_in)
  real*8 ::qt,integralOverXpoint_S
   real*8 :: X1,X2,X3,X4,X5
   real*8 :: y2,y3,y4,deltay
   real*8 :: yMin_in,yMax_in
   real*8::valueMax
   
   deltay=yMax_in-yMin_in
   y2=yMin_in+deltay/4d0
   y3=yMin_in+deltay/2d0
   y4=yMax_in-deltay/4d0  
   
   call TMDX_SIDIS_SetX(yMin_in)
   X1= xSec(qt)   
   call TMDX_SIDIS_SetX(y2)
   X2= xSec(qt)   
   call TMDX_SIDIS_SetX(y3)
   X3= xSec(qt)   
   call TMDX_SIDIS_SetX(y4)
   X4= xSec(qt)   
   call TMDX_SIDIS_SetX(yMax_in)
   X5= xSec(qt)
   
  
   !!approximate integral value
   valueMax=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   
   integralOverXpoint_S=IntegralOverXpoint_S_Rec(qt,yMin_in,y3,X1,X2,X3,valueMax)+&
	  IntegralOverXpoint_S_Rec(qt,y3,yMax_in,X3,X4,X5,valueMax)
  end function integralOverXpoint_S
  
  !!!! X1,X3,X5 are cross-sections at end (X1,X5) and central (X3) points of integraitons
  recursive function integralOverXpoint_S_Rec(qt,yMin_in,yMax_in,X1,X3,X5,valueMax) result(interX)
   real*8 ::qt
   real*8 :: interX,X1,X2,X3,X4,X5
   real*8 :: value,valueAB,valueACB
   real*8 :: yMin_in,yMax_in,y2,y3,y4,deltay
   real*8::valueMax,valueMaxNew,vv
   
   
   deltay=yMax_in-yMin_in
   y2=yMin_in+deltay/4d0
   y3=yMin_in+deltay/2d0
   y4=yMax_in-deltay/4d0
   
   valueMaxNew=valueMax
   
   call TMDX_SIDIS_SetX(y2)
   X2= xSec(qt)
   
   
   call TMDX_SIDIS_SetX(y4)
   X4= xSec(qt)
   
   
   valueAB=deltay*(X1+4d0*X3+X5)/6d0
   valueACB=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   If(ABS((valueACB-valueAB)/valueMax)>tolerance) then
    interX=integralOverXpoint_S_Rec(qt,yMin_in,y3,X1,X2,X3,valueMaxNew)&
	  +integralOverXpoint_S_Rec(qt,y3,yMax_in,X3,X4,X5,valueMaxNew)
   else
    interX=valueACB
   end if
   
  end function integralOverXpoint_S_Rec
  
    !---------------------------------UNTIGRATED over Z---------------------------------------------------------------
  
  !!!
  function Xsec_Zint(qt,zmin_in,zmax_in)
    real*8 :: qt,Xsec_Zint
    real*8 :: zmin_in,zmax_in

   if(zmax_in>1d0 .or. zmin_in<0d0) then
    write(*,*) 'arTeMiDe_SIDIS: CRITICAL ERROR limits of z-integration outside of (0,1)'
    write(*,*) 'Evaluation stop'
    stop
   end if
    
    Xsec_Zint=integralOverZpoint_S(qt,zmin_in,zmax_in)
  end function Xsec_Zint
  !--------------Simpsons--------------------
  !!!! parameter valueMax remembers the approximate value of integral to weight the tolerance.
  !!!! evaluation is done by adaptive simpson
  !!!! First we evaluate over 5 points and estimate the integral, and then split it to 3+3 and send to adaptive
  !!!! Thus minimal number of points =9
    function integralOverZpoint_S(qt,yMin_in,yMax_in)
  real*8 ::qt,integralOverZpoint_S
   real*8 :: X1,X2,X3,X4,X5
   real*8 :: y2,y3,y4,deltay
   real*8 :: yMin_in,yMax_in
   real*8::valueMax
   
   deltay=yMax_in-yMin_in
   y2=yMin_in+deltay/4d0
   y3=yMin_in+deltay/2d0
   y4=yMax_in-deltay/4d0
   
   call TMDX_SIDIS_SetZ(yMin_in)
   X1= xSec(qt)   
   call TMDX_SIDIS_SetZ(y2)
   X2= xSec(qt)   
   call TMDX_SIDIS_SetZ(y3)
   X3= xSec(qt)   
   call TMDX_SIDIS_SetZ(y4)
   X4= xSec(qt)   
   call TMDX_SIDIS_SetZ(yMax_in)
   X5= xSec(qt)
   
   !!approximate integral value
   valueMax=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   integralOverZpoint_S=IntegralOverZpoint_S_Rec(qt,yMin_in,y3,X1,X2,X3,valueMax)+&
	  IntegralOverZpoint_S_Rec(qt,y3,yMax_in,X3,X4,X5,valueMax)
  end function integralOverZpoint_S
  
  !!!! X1,X3,X5 are cross-sections at end (X1,X5) and central (X3) points of integraitons
  recursive function integralOverZpoint_S_Rec(qt,yMin_in,yMax_in,X1,X3,X5,valueMax) result(interX)
   real*8 ::qt
   real*8 :: interX,X1,X2,X3,X4,X5
   real*8 :: value,valueAB,valueACB
   real*8 :: yMin_in,yMax_in,y2,y3,y4,deltay
   real*8::valueMax,valueMaxNew,vv
   
   deltay=yMax_in-yMin_in
   y2=yMin_in+deltay/4d0
   y3=yMin_in+deltay/2d0
   y4=yMax_in-deltay/4d0
   
   valueMaxNew=valueMax
   
   call TMDX_SIDIS_SetZ(y2)
   X2= xSec(qt)
   
   call TMDX_SIDIS_SetZ(y4)
   X4= xSec(qt)
   
   valueAB=deltay*(X1+4d0*X3+X5)/6d0
   valueACB=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   If(ABS((valueACB-valueAB)/valueMax)>tolerance) then
    interX=integralOverZpoint_S_Rec(qt,yMin_in,y3,X1,X2,X3,valueMaxNew)&
	  +integralOverZpoint_S_Rec(qt,y3,yMax_in,X3,X4,X5,valueMaxNew)
   else
    interX=valueACB
   end if
   
  end function integralOverZpoint_S_Rec
  
  
  
  !---------------------------------UNTIGRATED over Q---------------------------------------------------------------
  function Xsec_Qint(qt,Q_min,Q_max)
    real*8:: qt,Xsec_Qint
    real*8:: Q_min,Q_max
    
    Xsec_Qint=integralOverQpoint_S(qt,Q_min,Q_max)
  end function Xsec_Qint
  
  !--------------Simpsons--------------------
  !!!! parameter valueMax remembers the initial value of integral to weight the tolerance.
  !!!! First we evaluate over 5 points and estimate the integral, and then split it to 3+3 and send to adaptive
  !!!! Thus minimal number of points =9
  !!!! taking into account minimum calls of y-integral we have  =81 points
  function integralOverQpoint_S(qt,QMin_in,QMax_in)
  real*8 ::qt,integralOverQpoint_S
   real*8 :: X1,X2,X3,X4,X5
   real*8 :: QMin_in,QMax_in
   real*8::valueMax,Q2,Q3,Q4,deltaQ
   
   deltaQ=QMax_in-QMin_in
   Q2=QMin_in+deltaQ/4d0
   Q3=QMin_in+deltaQ/2d0
   Q4=QMax_in-deltaQ/4d0
   
   call TMDX_SIDIS_SetQ(QMin_in)
   X1=2*Q_global*Xsec(qt)
   
   call TMDX_SIDIS_SetQ(Q2)
   X2=2*Q_global*Xsec(qt)
   
   call TMDX_SIDIS_SetQ(Q3)
   X3=2*Q_global*Xsec(qt)
   
   call TMDX_SIDIS_SetQ(Q4)
   X4=2*Q_global*Xsec(qt)
   
   call TMDX_SIDIS_SetQ(QMax_in)
   X5=2*Q_global*Xsec(qt)
   
      !!approximate integral value
   valueMax=deltaQ*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   integralOverQpoint_S=IntegralOverQpoint_S_Rec(qt,QMin_in,Q3,X1,X2,X3,valueMax)+&
	IntegralOverQpoint_S_Rec(qt,Q3,QMax_in,X3,X4,X5,valueMax)
  end function integralOverQpoint_S
  
  !!!! X1,X3,X5 are cross-sections at end (X1,X5) and central (X3) points of integraitons
  recursive function integralOverQpoint_S_Rec(qt,QMin_in,QMax_in,X1,X3,X5,valueMax) result(interX)
   real*8 ::qt
   real*8 :: interX,X1,X2,X3,X4,X5
   real*8 :: valueAB,valueACB
   real*8 :: QMin_in,QMax_in,Q2,Q3,Q4,deltaQ
   real*8::valueMax,valueMaxNew
   
   deltaQ=QMax_in-QMin_in
   Q2=QMin_in+deltaQ/4d0
   Q3=QMin_in+deltaQ/2d0
   Q4=QMax_in-deltaQ/4d0
   
   valueMaxNew=valueMax
   
   call TMDX_SIDIS_SetQ(Q2)
   X2=2*Q_global*Xsec(qt)
      
   call TMDX_SIDIS_SetQ(Q4)
   X4=2*Q_global*Xsec(qt)
   
   valueAB=deltaQ*(X1+4d0*X3+X5)/6d0
   valueACB=deltaQ*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   If(ABS((valueACB-valueAB)/valueMaxNew)>tolerance) then
    interX=integralOverQpoint_S_Rec(qt,QMin_in,Q3,X1,X2,X3,valueMaxNew)&
	  +integralOverQpoint_S_Rec(qt,Q3,Qmax_in,X3,X4,X5,valueMaxNew)
   else
    interX=valueACB
   end if
  end function integralOverQpoint_S_Rec
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!INTERFACES TO CALCULATING CROSS-SECTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !---------------------------------UNTIGRATED------------------------------------------------------------------
  
  !!qt_list is the list of requred qt -point,
  !! X_list is variable to store results (should be of the same dimension as qt_list)
  subroutine xSecList(X_list,qt_List)
    real*8, intent(in) :: qt_list(:)
    real*8, intent(out) :: X_list(:)
    integer :: i,length
    length=size(qt_list)
    CallCounter=CallCounter+length
     do i=1,length
       X_List(i)=PreFactor1()*xSec(qt_List(i))
     end do
  end subroutine xSecList
  
  !!!!Evaluate differential xSec at single point
  subroutine xSecSingle(X,qT_in)
   real*8:: X,qT_in
   CallCounter=CallCounter+1
   X=PreFactor1()*xSec(qT_in)
  end subroutine xSecSingle
  
  
!---------------------------------UNTIGRATED over X---------------------------------------------------------------
    !!qt_list is the list of requred qt -point,
  !! X_list is variable to store results (should be of the same dimension as qt_list)
  subroutine xSecList_Xint(X_list,qt_List,xMin_in,xMax_in)
    real*8, intent(in) :: qt_list(:)
    real*8, intent(out) :: X_list(:)
    real*8::xMin_in,xMax_in
    integer :: i,length
    length=size(qt_list)
    CallCounter=CallCounter+length
     do i=1,length
       X_List(i)=PreFactor1()*Xsec_Xint(qt_List(i),xMin_in,xMax_in)
     end do
  end subroutine xSecList_Xint
  
  !!
  subroutine xSecSingle_Xint(X,qt,xMin_in,xMax_in)
    real*8::X,qT
    real*8::xMin_in,xMax_in
    
   CallCounter=CallCounter+1
   X=PreFactor1()*Xsec_Xint(qt,xMin_in,xMax_in)
  end subroutine xSecSingle_Xint
  
  !---------------------------------UNTIGRATED over Z---------------------------------------------------------------
    !!qt_list is the list of requred qt -point,
  !! X_list is variable to store results (should be of the same dimension as qt_list)
  subroutine xSecList_Zint(X_list,qt_List,zMin_in,zMax_in)
    real*8, intent(in) :: qt_list(:)
    real*8, intent(out) :: X_list(:)
    real*8::zMin_in,zMax_in
    integer :: i,length
    length=size(qt_list)
    CallCounter=CallCounter+length
     do i=1,length
       X_List(i)=PreFactor1()*Xsec_Zint(qt_List(i),zMin_in,zMax_in)
     end do
  end subroutine xSecList_Zint
  
  !!
  subroutine xSecSingle_Zint(X,qt,zMin_in,zMax_in)
    real*8::X,qT
    real*8::zMin_in,zMax_in
    
   CallCounter=CallCounter+1
   X=PreFactor1()*Xsec_Zint(qt,zMin_in,zMax_in)
  end subroutine xSecSingle_Zint
  
  
  !---------------------------------UNTIGRATED over Q---------------------------------------------------------------
  !!qt_list is the list of requred qt -point,
  !! X_list is variable to store results (should be of the same dimension as qt_list)
  subroutine xSecList_Qint(X_list,qt_List,Q_min,Q_max)
    real*8, intent(in) :: qt_list(:)
    real*8, intent(out) :: X_list(:)
    real*8::Q_min,Q_max
    integer :: i,length
    length=size(qt_list)
    CallCounter=CallCounter+length
     do i=1,length
       X_List(i)=PreFactor1()*Xsec_Qint(qt_List(i),Q_min,Q_max)
     end do
  end subroutine xSecList_Qint
  
  !!
  subroutine xSecSingle_Qint(X,qt,Q_min,Q_max)
    real*8::X,qT
    real*8::Q_min,Q_max
    
   CallCounter=CallCounter+1
   X=PreFactor1()*Xsec_Qint(qt,Q_min,Q_max)
  end subroutine xSecSingle_Qint
  
  
end module TMDX_SIDIS
