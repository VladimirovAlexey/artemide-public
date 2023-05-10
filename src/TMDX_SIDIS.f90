!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 2.02
!
!	Evaluation of the TMD cross-section for SIDIS-like cross-sections
!	
!	if you use this module please, quote ????.????
!
!	ver 1.2: release (AV, 15.12.2017)
!	ver 1.32: part of functions migrated to TMDF, rest updated (AV, 16.08.2018)
!	ver 2.02:							(AV,16.08.2019)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDX_SIDIS
use aTMDe_Numerics
use IO_functions
use TMDF
use QCDinput
use EWinput

implicit none
private

!Current version of module
character (len=10),parameter :: moduleName="TMDX-SIDIS"
character (len=5),parameter :: version="v2.02"
!Last appropriate verion of constants-file
integer,parameter::inputver=9


logical::started=.false.
integer::outputlevel
integer::messageTrigger

real(dp)::hc2

!The variables for all parameters of the model!
!it is used only as input if these parameters are not set by the user.
real(dp):: Q_global,x_global,z_global,s_global
!! Set of process definition, for Prefactor 1, Prefactor 2, structure function, etc
!! = (/p1,p2,p3/)
integer,dimension(1:3)::process_global
!! cut prarameters
logical::includeCuts_global
!! (yMin,yMax,Wmin,Wmax)
real(dp),dimension(1:4)::CutParameters_global
!!! Target mass of target (squared)
real(dp):: M2_target_global=0.939d0**2
!!! Produced mass (squared)
real(dp):: M2_product_global=0.1d0**2

!!other global parameters, which are defined upon initialization
integer:: orderH_global
!!inclusion of power corrections to the definition
!! corrQT ~ qT/Q in cross-section
!! corrM1 ~ M1/Q
!! corrM2 ~ M2/Q
!! corrQTinX1Z1 ~ qT/Q in x1,z1
logical:: corrQT,corrM1,corrM2,corrQTinX1Z1 !!!if true include
integer:: exactScales  !!!=1 if exact hard scales = true, =0 otherwise


real(dp) :: tolerance=0.0005d0
!!! number of sections for PT-integral by default
integer::NumPTdefault=4

real(dp)::c2_global!,muHard_global

integer::messageCounter
integer::GlobalCounter
integer::CallCounter

!!!------------------------------Global Integration definitions-------------------
!-integration over Z
real(dp) :: toleranceZ=0.0005d0
!!! SA = adaptive simpson, S5 = Simpson with 5-points, I0 = center point x size
character(len=4)::methodZ='SA'
!-integration over X
real(dp) :: toleranceX=0.0005d0
!!! SA = adaptive simpson, S5 = Simpson with 5-points, I0 = center point x size
character(len=4)::methodX='SA'


public::TMDX_SIDIS_setProcess,TMDX_SIDIS_ShowStatistic,TMDX_SIDIS_Initialize,&
    TMDX_SIDIS_XSetup,TMDX_SIDIS_SetCuts,TMDX_SIDIS_IsInitialized,TMDX_SIDIS_ResetCounters,TMDX_SIDIS_SetScaleVariation

public::CalcXsec_SIDIS,CalcXsec_SIDIS_Zint_Xint_Qint,CalcXsec_SIDIS_PTint_Zint_Xint_Qint,xSec_SIDIS,xSec_SIDIS_List,&
    xSec_SIDIS_List_forharpy,xSec_SIDIS_BINLESS_List_forharpy

interface TMDX_SIDIS_setProcess
    module procedure TMDX_setProcess1,TMDX_setProcess3,TMDX_setProcess30
end interface

interface CalcXsec_SIDIS
    module procedure CalcXsecLIST_SIDIS,CalcXsecSINGLE_SIDIS
end interface

interface CalcXsec_SIDIS_Zint_Xint_Qint
    module procedure CalcXsecLIST_SIDIS_Zint_Xint_Qint,CalcXsecSINGLE_SIDIS_Zint_Xint_Qint
end interface

interface CalcXsec_SIDIS_PTint_Zint_Xint_Qint
    module procedure CalcXsecLISTLIST_SIDIS_PTint_Zint_Xint_Qint,&
        CalcXsecLIST_SIDIS_PTint_Zint_Xint_Qint,CalcXsecSINGLE_SIDIS_PTint_Zint_Xint_Qint
end interface

contains

function TMDX_SIDIS_IsInitialized()
    logical::TMDX_SIDIS_IsInitialized
    TMDX_SIDIS_IsInitialized=started
end function TMDX_SIDIS_IsInitialized

!! Initialization of the package
subroutine TMDX_SIDIS_Initialize(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=300)::path,line
    logical::initRequired,dummyLogical
    character(len=8)::orderMain
    integer::i,FILEver
    !$ integer:: omp_get_thread_num

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
        write(*,*) '		     Update the const-file with artemide.setup'
        write(*,*) '  '
        stop
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
        CASE ("LO+")
            orderH_global=0
        CASE ("NLO")
            orderH_global=1
        CASE ("NLO+")
            orderH_global=1
        CASE ("NNLO")
            orderH_global=2
        CASE ("N2LO")
            orderH_global=2
        CASE ("NNLO+")
            orderH_global=3
        CASE ("NNNLO")
            orderH_global=3
        CASE ("N3LO") !! same as NNNLO
            orderH_global=3
        CASE ("N3LO+")
            orderH_global=3
        CASE ("N4LO")
            orderH_global=4
        CASE DEFAULT
            if(outputLevel>0) write(*,*) WarningString('try to set unknown order. Switch to NLO.',moduleName)
            orderH_global=1
    END SELECT
    if(outputLevel>1) write(*,*) '	artemide.TMDX_SIDIS: the used order is ',trim(orderMain)
    
    !! qT correction in kinematics
    call MoveTO(51,'*p2   ')
    read(51,*) corrQT
    if(outputLevel>2 .and. corrQT) write(*,*) '	artemide.TMDX_SIDIS: qT/Q corrections in kinematics are included.'
    !! Target mass corrections
    call MoveTO(51,'*p3   ')
    read(51,*) corrM1
    if(outputLevel>2 .and. corrM1) write(*,*) '	artemide.TMDX_SIDIS: target mass corrections in kinematics are included.'
    !! produced mass corrections
    call MoveTO(51,'*p4   ')
    read(51,*) corrM2
    if(outputLevel>2 .and. corrM2) write(*,*) '	artemide.TMDX_SIDIS: product mass corrections in kinematics are included.'
    !! qT correction in x1 z1
    call MoveTO(51,'*p5   ')
    read(51,*) corrQTinX1Z1
    if(outputLevel>2 .and. corrQTinX1Z1) write(*,*) '	artemide.TMDX_SIDIS: qT/Q corrections in x1,z1 are included.'
    !!exact values for scales
    call MoveTO(51,'*p6   ')
    read(51,*) dummyLogical
    if(dummyLogical) then 
        exactScales=1
    else
        exactScales=0
    end if
    if(outputLevel>2 .and. dummyLogical) write(*,*) '	artemide.TMDX_SIDIS: qT/Q correction for scales variables are included.'
    
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) tolerance
    call MoveTO(51,'*p2  ')
    read(51,*) NumPTdefault
    call MoveTO(51,'*p3  ')
    read(51,*) toleranceZ
    call MoveTO(51,'*p4  ')
    read(51,*) methodZ
    if((.not.(methodZ=='SA'.or.methodZ=='S5'.or.methodZ=='I0')).and.outputLevel>0) then
        write(*,*) ErrorString(' method for z-bin-integration is unknow: '//methodZ,moduleName)
        write(*,*) color('switching to Adaptive-Simpson',c_red)
        methodZ='SA'
    end if
    call MoveTO(51,'*p5  ')
    read(51,*) toleranceX
    call MoveTO(51,'*p6  ')
    read(51,*) methodX
    if((.not.(methodX=='SA'.or.methodX=='S5'.or.methodX=='I0')).and.outputLevel>0) then
        write(*,*) ErrorString(' method for x-bin-integration is unknow: '//methodX,moduleName)
        write(*,*) color('switching to Adaptive-Simpson',c_red)
        methodX='SA'
    end if
        
    !$    if(outputLevel>1) write(*,*) '	artemide.TMDX_SIDIS: parallel evaluation of cross-sections is to be used'
    !$    call MoveTO(51,'*C   ')
    !$    call MoveTO(51,'*p1  ')
    !$    read(51,*) i
    !$    call OMP_set_num_threads(i)
    !$    if(outputLevel>1) write(*,*) '	artemide.TMDX_SIDIS: number of threads for parallel evaluation is set to ', i	

    !$     if(outputLevel>2) write(*,*) '------TEST OF PARALLEL PROCESSING ----------'
    !$OMP PARALLEL
    !$     if(outputLevel>2) write(*,*) '   artemide.TMDX_SIDIS:thread num ',  omp_get_thread_num(), ' ready.'
    !$OMP END PARALLEL
    CLOSE (51, STATUS='KEEP')

    if(.not.EWinput_IsInitialized()) then
        if(outputLevel>1) write(*,*) '.. initializing EWinput (from ',moduleName,')'
        if(present(prefix)) then
            call EWinput_Initialize(file,prefix)
        else
            call EWinput_Initialize(file)
        end if
    end if
        
    if(.not.TMDF_IsInitialized()) then
        if(outputLevel>1) write(*,*) '.. initializing TMDF (from ',moduleName,')'
        if(present(prefix)) then
            call TMDF_Initialize(file,prefix)
        else
            call TMDF_Initialize(file)
        end if
    end if

    includeCuts_global=.false.
    c2_global=1d0
    
    GlobalCounter=0
    CallCounter=0
    messageCounter=0
    
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
    real(dp)::c2_in

    if(outputLevel>1) write(*,*) 'TMDX_SIDIS: scale variation constant c2 reset:',c2_in

    if(c2_in<0.1d0 .or. c2_in>10.d0) then
        if(outputLevel>0) write(*,*) WarningString('variation in c2 is enourmous. c2 is set to 2',moduleName)
        c2_global=2d0
    else
        c2_global=c2_in
    end if

end subroutine TMDX_SIDIS_SetScaleVariation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!! PROCESS DEFINITION
  
!!!set variables for process definition
subroutine TMDX_setProcess30(p0)
    integer,dimension(1:3)::p0
    process_global=p0
end subroutine TMDX_setProcess30
  
!!!set variables for process definition
subroutine TMDX_setProcess3(p1,p2,p3)
    integer::p1,p2,p3
    process_global=(/p1,p2,p3/)
end subroutine TMDX_setProcess3
  
!!!set variables for process definition
subroutine TMDX_setProcess1(p)
    integer::p
    call TMDX_setProcess30(processArrayFromInteger(p))
end subroutine TMDX_setProcess1
  
    
function processArrayFromInteger(p)
    integer,intent(in)::p
    integer,dimension(1:3)::processArrayFromInteger
    SELECT CASE(p)
        case default
            write(*,*) color('ERROR: arTeMiDe_SIDIS: unknown process is called. p=',c_red_bold),p
            write(*,*) 'Evaluation stop'
            stop
    end SELECT
end function processArrayFromInteger
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS FOR OPERATION WITH KINEMATICS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  
subroutine TMDX_SIDIS_XSetup(s,z,x,Q,mTARGET,mPRODUCT)
    real(dp)::s,Q,x,z
    real(dp),optional::mTARGET,mPRODUCT

    if(.not.started) then
    write(*,*) ErrorString('module is not initialized. Evaluation terminated',moduleName)
    stop
    end if

    s_global=s
    Q_global=Q
    x_global=x    
    z_global=z

    if(PRESENT(mTARGET)) then
        M2_target_global=mTARGET**2
    !       write(*,*) 'mTarget reset ->',M2_target_global
    end if

    if(PRESENT(mPRODUCT)) then
        M2_product_global=mPRODUCT**2
    !       write(*,*) 'mProduct reset ->',M2_product_global
    end if

end subroutine TMDX_SIDIS_XSetup
  
!!! function makes kinematic array from the given set of pT,s,Q,x,z
!!! array has 13 often appearing entries
! 1 = pT
! 2 = Q
! 3 = Q^2
! 4 = x
! 5 = z
! 6 = y 	=Q^2/x(s-M^2)
! 7 = epsilon	= (1-y-gamma2 y^2/4)/(1-y+y^2/2+gamma2 y^2/4)
! 8 = gamma2	= (2 x M/Q)^2
! 9=rho2*gamma2	= (m/z/Q)^2*gamma2
! 10=rhoT2*gamma2	= (m^2+pt^2)/(z Q)^2*gamma2
! 11=sM2	= s-M^2
! 12=M2-target
! 13=M2-product
pure function kinematicArray(pT,s,z,x,Q,M2target_in,M2product_in)
real(dp),dimension(1:13)::kinematicArray
real(dp),intent(in)::pT,s,z,x,Q,M2target_in,M2product_in
real(dp)::Q2,y,varepsilon,gamma2,rho2,rhoPEPR2,sM2,M2target,M2product

    Q2=Q**2
    if(corrM1) then
        M2target=M2target_in
        sM2=s-M2target
        gamma2=4d0*M2target*x**2/Q2
    else
        M2target=0d0
        gamma2=0d0
        sM2=s
    end if
    
    y=Q2/x/sM2!YfromSXQ2(sM2,x,Q2)

    varepsilon=(1d0-y-y**2*gamma2*0.25d0)/(1d0-y+y**2*(0.5d0+0.25d0*gamma2))

    if(corrM2) then
        M2product=M2product_in
        rho2=gamma2*M2product/Q2/z**2
    else
        M2product=0d0
        rho2=0d0
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
!!qT 	=pT/z sqrt( ( 1+gamma2) / (1-gamma2 rho2))
pure subroutine CalculateqT(qT,var)
    real(dp),intent(out)::qT
    real(dp),dimension(1:13),intent(inout)::var

    qT=var(1)/var(5)*Sqrt((1d0+var(8))/(1d0-var(9)))

end subroutine CalculateqT

!!!!!Evaluate the parameters of the TMD factorized Fourier-integral
!!qT 	=pT/z sqrt( ( 1+gamma2) / (1-gamma2 rho2))
!!fac1	= -2/gamma2*(1-sqrt(1+gamma2*(1-qT^2/Q^2)))
!!X1=x*fac1
!!Z1=z*fac1*(...)
pure subroutine CalculateX1Z1qT(x1,z1,qT,var)
    real(dp),intent(out)::x1,z1,qT
    real(dp),dimension(1:13),intent(inout)::var
    real(dp)::fac1

    qT=var(1)/var(5)*Sqrt((1d0+var(8))/(1d0-var(9)))

    !!fac1 appears in definition of x1, and z1
    if(corrM1 .and. var(8)>0) then
        if(corrQTinX1Z1) then
            fac1=-2d0/var(8)*(1d0-sqrt(1d0+var(8)*(1d0-(qT/var(2))**2)))
        else
            fac1=-2d0/var(8)*(1d0-sqrt(1d0+var(8)))
        end if
    else
        if(corrQTinX1Z1) then
            fac1=(1d0-(qT/var(2))**2)
        else
            fac1=1d0
        end if
    end if

    x1=var(4)*fac1

    if (corrQTinX1Z1) then
        z1=var(5)*fac1*(1d0+sqrt(1d0-var(9)))/(2d0*(1d0-(qT/var(2))**2))
    else
        z1=var(5)*fac1*(1d0+sqrt(1d0-var(9)))*0.5d0
    end if

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
    !!!Q2 same
    !!!sM2 same
    !!!gamma
    var(8)=4d0*var(12)*x**2/var(3)
    !!!y
    var(6)=YfromSXQ2(var(11),x,var(3))
    !!epsilon
    var(7)=(1d0-var(6)-var(6)**2*var(8)*0.25d0)/(1d0-var(6)+var(6)**2*(0.5d0+0.25d0*var(8)))

    !! rescale rho'2 with new gamma
    if(g2>0) then
        var(9)=var(9)*var(8)/g2
        var(10)=var(10)*var(8)/g2
    else
        var(9)=0d0
        var(10)=0d0
    end if

    !!!masses same
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
    !   var=kinematicArray(var(1),var(11)+var(12),z,var(4),var(2),var(12),var(13))
    var(5)=z
    !!Q2 same
    !!gamma2 same
    !!y same
    !!varepsilon same

    !!rho2 new
    var(9)=var(8)*var(13)/var(3)/z**2

    !!rhoPerp new
    if(corrQT) then
        var(10)=var(9)+var(8)*(var(1)/z/var(2))**2
    else
        var(10)=var(9)
    end if

end subroutine SetZ
  
!!!! update a given kinematic array with new value of pt.
pure subroutine SetPT(pt,var)
    real(dp),dimension(1:13),intent(inout)::var
    real(dp),intent(in)::pt
    !   var=kinematicArray(pt,var(11)+var(12),var(5),var(4),var(2),var(12),var(13))
    var(1)=pt

    !!! Q2 same
    !!! gamma2 same
    !!! y same
    !!! varepsilon same
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
  
!!!! Set Prefactor1
!!!! it multiplies the cross-section as a whole, does not participate in the integration
pure function PreFactor1(p1)
    real(dp)::Prefactor1
    integer,intent(in)::p1
    PreFactor1=1d0
end function PreFactor1
  
!!!!! Prefactor 2 is (universal part) x H
function PreFactor2(var,process,x1,z1,qT)
    real(dp),dimension(1:13),intent(in)::var
    integer,dimension(1:3),intent(in)::process
    real(dp)::PreFactor2,uniPart,phasePart,x1,z1,qT,fac1,scaleMu

    !!!! universal part
    
    !!!! If exact scales, zeta*zeta=Q^2-qT^2.
    scaleMu=sqrt(var(3)-exactScales*qT**2)

    SELECT CASE(process(2))
        case(-10221191)
            uniPart=1d0
        CASE(1)
            !!!!!!!!!!!!!!!!!!!!! CHECK FACTOR 2???!
        
            !2 pi aEm^2/Q^4 y^2/(1-epsilon)*z1/z(1+varepsilon...)
            ! prefactor for unpolarized expression

            !!! this is 1+qT^2/Q^2(e-gamma^2/2)/(1+gamma^2)
            if(corrQT) then
                fac1=(1d0+(qT**2/var(3))*(var(7)-0.5d0*var(8))/(1+var(8)))/sqrt(1-var(10))
            else
                fac1=1d0
            end if
            uniPart=pix2*alphaEM(scaleMu)**2/(var(3)**2)*(var(6)**2/(1d0-var(7)))*(z1/var(5))*&
            fac1*&
            HardCoefficientSIDIS(scaleMu)*&
            hc2*1d9!from GeV to pbarn

        CASE(2)
            ! prefactor for FUU,T
            !!! this is 1+qT^2/Q^2(e-gamma^2/2)/(1+gamma^2)
            if(corrQT) then
                fac1=(1d0+(qT**2/var(3))*(var(7)-0.5d0*var(8))/(1+var(8)))/sqrt(1-var(10))
            else
                fac1=1d0
            end if
            uniPart=var(4)/pi/(1d0+0.5d0*var(8)/var(4))*&
            fac1*&
            HardCoefficientSIDIS(scaleMu)
        CASE DEFAULT 
            write(*,*) ErrorString(' unknown process p2=',moduleName),&
                process(2),color(' .Evaluation stop.',c_red_bold)
            stop
    END SELECT

    SELECT CASE(process(1))
        CASE(2) 
            !!! case of dQ2 -> dy
            phasePart=var(3)/var(6)
        CASE(3) 
            !!! case of dx -> dy
            phasePart=var(4)/var(6)
        CASE DEFAULT
            phasePart=1d0
    END SELECT

    PreFactor2=phasePart*uniPart
end function PreFactor2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! CUTS RELATED FUNCTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!The cuts are ymin<y<ymax,Wmin<W2<Wmax
subroutine TMDX_SIDIS_SetCuts(inc,yMin,yMax,Wmin,Wmax)
    logical::inc
    real(dp)::yMin,yMax,Wmin,Wmax
    real(dp)::y0,y1,W0,W1
    includeCuts_global=inc

    if(.not.inc) return

    if(yMin<0d0) then 
        if(outputlevel>0) write(*,*) WarningString('SetCuts: yMin<0. Set to 0',moduleName)
        y0=0d0
    else if(ymin>1d0) then
        if(outputlevel>0) write(*,*) WarningString('SetCuts: yMin>1. Set to 1',moduleName)
        y0=1d0
    else
        y0=yMin
    end if
    if(yMax<0d0) then 
        if(outputlevel>0) write(*,*) WarningString('SetCuts: yMax<0. Set to 0',moduleName)
        y1=0d0
    else if(ymax>1d0) then
        if(outputlevel>0) write(*,*)  WarningString('SetCuts: yMax>1. Set to 1',moduleName)
        y1=1d0
    else
        y1=yMax
    end if

    if(Wmin>0) then
        W0=Wmin
    else
        if(outputlevel>0) write(*,*) WarningString('SetCuts: Wmin<0. Set to 0',moduleName)
        W0=0d0
    end if
    if(Wmax>0) then
        W1=Wmax
    else
        if(outputlevel>0) write(*,*) WarningString('SetCuts: Wmax<0. Set to 0',moduleName)
        W1=0d0
    end if

    if(y0<=y1) then
        if(W0<=W1) then
            CutParameters_global=(/y0,y1,W0,W1/)
        else
            if(outputlevel>0) write(*,*) WarningString('SetCuts: WMin>WMax. Values exchanged',moduleName)
            CutParameters_global=(/y0,y1,W1,W0/)
        end if
    else
        if(outputlevel>0) write(*,*) WarningString('SetCuts: yMin>yMax. Values exchanged',moduleName)
        if(W0<=W1) then
            CutParameters_global=(/y1,y0,W0,W1/)
        else
            if(outputlevel>0) write(*,*) WarningString('SetCuts: WMin>WMax. Values exchanged',moduleName)
            CutParameters_global=(/y1,y0,W1,W0/)
        end if
    end if
end subroutine TMDX_SIDIS_SetCuts
  
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

!!! this is help function which evaluate xSec at single qt (without lists) with only prefactor 2
!!!! this is extended (and default) version of xSec, which include all parameters
!!! note that it is calculated with respect to qT
function xSec(var,process)
    real(dp):: xSec,FF
    real(dp)::x1,z1,qT,scaleMu,scaleZeta
    real(dp),dimension(1:13)::var
    integer,dimension(1:3),intent(in)::process
    integer::OMP_get_thread_num
    GlobalCounter=GlobalCounter+1
    call CalculateX1Z1qT(x1,z1,qT,var)
    
    !!! setting values of scales
    scaleZeta=var(3)-exactScales*qT**2  !! zeta=Q2+qT^2
    scaleMu=sqrt(scaleZeta) 

    FF=TMDF_F(var(3),qT,x1,z1,scaleMu*c2_global,scaleZeta,scaleZeta,process(3))
    xSec=PreFactor2(var,process,x1,z1,qT)*FF  

    !write(*,*) "{",var(3),",",x1,"},"!,z1
end function xSec
  
    
  !---------------------------------INTEGRATED over Z (and Q and X)--------------------------------------------------------------
  
  !!! the variable doZ check should the integration over Z be performed
  !!! if doZ=true, the integration is done
  !!! if doZ=facle the single value (at xMin) is returned
  function Xsec_Zint(var,process,doZ,zMin,zMax)
    real(dp),dimension(1:13):: var
    logical::doZ
    real(dp) :: Xsec_Zint
    real(dp) :: zMin,zMax
    integer,dimension(1:3),intent(in)::process
    
    !! the integration over Z is required
    if(doZ) then
      
      if(zmax > 1d0) then
	  if(outputlevel>1) write(*,*) WarningString('upper limit of z-integration is >1. It is set to 1.',moduleName)
        zmax=1d0
      end if
      if(zmin < 0.000001d0) then
        write(*,*) ErrorString('lower limit of z-integration is < 10^{-6}. Evaluation stop.',moduleName)
        stop
      end if
      
      if(zMin<zMax) then
        if(methodZ=='I0') then !!!! central point estimation        
            call SetZ((Zmin+zMax)/2d0,var)
            Xsec_Zint=Xsec(var,process)*(zMax-zMin)
            
        else !!!! so far only Simpsons        
            Xsec_Zint=integralOverZpoint_S(var,process,zMin,zMax)
        end if
      else
        Xsec_Zint=0d0
      end if
    else    
      ! no integration over Z
      
      if(zmax > 1d0) then
        if(outputlevel>1) write(*,*) WarningString('upper limit of z-integration is >1. It is set to 1.',moduleName)
            zmin=1d0
        end if
        if(zmin < 0.000001d0) then
            write(*,*) ErrorString('lower limit of z-integration is < 10^{-6}. Evaluation stop.',moduleName)
            stop
        end if
        
        call SetZ((Zmin+zMax)/2d0,var)
        Xsec_Zint=Xsec(var,process)
    end if
    
  end function Xsec_Zint
  !--------------Simpsons--------------------
  !!!! parameter valueMax remembers the approximate value of integral to weight the tolerance.
  !!!! evaluation is done by adaptive simpson
  !!!! First we evaluate over 5 points and estimate the integral, and then split it to 3+3 and send to adaptive
  !!!! Thus minimal number of points =9
  function integralOverZpoint_S(var,process,yMin_in,yMax_in)
   real(dp),dimension(1:13)::var
   integer,dimension(1:3),intent(in)::process
   real(dp) ::integralOverZpoint_S
   real(dp) :: X1,X2,X3,X4,X5
   real(dp) :: y2,y3,y4,deltay
   real(dp) :: yMin_in,yMax_in
   real(dp)::valueMax
   
   deltay=yMax_in-yMin_in
   y2=yMin_in+deltay/4d0
   y3=yMin_in+deltay/2d0
   y4=yMax_in-deltay/4d0
   
   call SetZ(yMin_in,var)
   X1=Xsec(var,process)  
   call SetZ(y2,var)
   X2=Xsec(var,process)   
   call SetZ(y3,var)
   X3=Xsec(var,process)  
   call SetZ(y4,var)
   X4=Xsec(var,process)  
   call SetZ(yMax_in,var)
   X5=Xsec(var,process)  
   
   !!approximate integral value
   valueMax=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   if(methodZ=='SA') then
    integralOverZpoint_S=IntegralOverZpoint_S_Rec(var,process,yMin_in,y3,X1,X2,X3,valueMax)+&
	  IntegralOverZpoint_S_Rec(var,process,y3,yMax_in,X3,X4,X5,valueMax)
   else if(methodZ=='S5') then
    integralOverZpoint_S=valueMax !!!! 5-point integration
   else
    write(*,*) ErrorString(' integralOverZpoint_S incorrect method. Check methodZ variable and initialization',moduleName)
   end if
  end function integralOverZpoint_S
  
  !!!! X1,X3,X5 are cross-sections at end (X1,X5) and central (X3) points of integraitons
  recursive function integralOverZpoint_S_Rec(var,process,yMin_in,yMax_in,X1,X3,X5,valueMax)  result(interX)
   real(dp),dimension(1:13) ::var
   integer,dimension(1:3),intent(in)::process
   real(dp) :: interX,X1,X2,X3,X4,X5
   real(dp) :: value,valueAB,valueACB
   real(dp) :: yMin_in,yMax_in,y2,y3,y4,deltay
   real(dp)::valueMax,vv
   
   deltay=yMax_in-yMin_in
   y2=yMin_in+deltay/4d0
   y3=yMin_in+deltay/2d0
   y4=yMax_in-deltay/4d0
   
   call SetZ(y2,var)
   X2=Xsec(var,process)
   
   call SetZ(y4,var)
   X4=Xsec(var,process)
   
   valueAB=deltay*(X1+4d0*X3+X5)/6d0
   valueACB=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   If(ABS((valueACB-valueAB)/valueMax)>toleranceZ) then
    interX=integralOverZpoint_S_Rec(var,process,yMin_in,y3,X1,X2,X3,valueMax)&
	  +integralOverZpoint_S_Rec(var,process,y3,yMax_in,X3,X4,X5,valueMax)
   else
    interX=valueACB
   end if
   
  end function integralOverZpoint_S_Rec
  
  !---------------------------------INTEGRATED over X (and Z)---------------------------------------------------------------
  
  !!!
  !!! the variable doX check should the integration over x be performed
  !!! if doX=true, the integration is done
  !!! if doX=false the single value (at xMin) is returned
  function Xsec_Zint_Xint(var,process,doZ,zMin,zMax,doX,Xmin_in,Xmax_in,doCut,Cuts)
    real(dp),dimension(1:13) :: var
    logical::doX,doZ,doCut
    real(dp),dimension(1:4),intent(in)::Cuts
    real(dp) :: Xsec_Zint_Xint
    real(dp) :: xmin, xmax,xmin_in,xmax_in,zMin,zMax
    integer,dimension(1:3),intent(in)::process
    
    if(doX) then    
      !!!Integration is required
      
      !!! in the case process=3 the input is y, which is to be transformed to X
      !!! evaluate corresponding y's
      if(process(1)==3) then
        xmin=XfromSYQ2(var(11),xmin_in,var(3))
        xmax=XfromSYQ2(var(11),xmax_in,var(3))
      else
        xmin=xmin_in
        xmax=xmax_in
      end if
    
      if(xmax > 1d0) then
        if(outputlevel>1) write(*,*) WarningString('upper limit of x-integration is >1. It is set to 1.',moduleName)
            xmax=1d0
        end if
        if(xmin < 0.000001d0) then
            write(*,*) ErrorString('lower limit of x-integration is < 10^{-6}. Evaluation stop.',moduleName)
        stop
      end if
      
      !! in case of cut we determine recut values
      if(doCut) then
        xmin=xMinWithCuts(xmin,var,Cuts)
        xmax=xMaxWithCuts(xmax,var,Cuts)
      end if
      
      if(xmin<xmax) then
        if(methodX=='I0') then !! central point estimaton
            call SetX((xMin+xMax)/2d0,var)
            Xsec_Zint_Xint=Xsec_Zint(var,process,doZ,zMin,zMax)*(xMax-xMin)
        else !!! so far only Simpsons
            Xsec_Zint_Xint=integralOverXpoint_S(var,process,doZ,zMin,zMax,xmin,xmax)
        end if
      else!!! it is possible that cuts cut out the integration range completely
        Xsec_Zint_Xint=0d0
      end if
      
    else
      !! no integration
      !! just single point
      if(process(1)==3) then
        xmin=XfromSYQ2(var(11),(xmin_in+xmax_in)/2d0,var(3))
        xmax=XfromSYQ2(var(11),(xmin_in+xmax_in)/2d0,var(3))
      else
        xmin=xmin_in
        xmax=xmax_in
      end if
      
      if(xmin > 1d0) then
	  if(outputlevel>1) write(*,*) WarningString('upper limit of x-integration is >1. It is set to 1.',moduleName)
	  xmin=1d0
      end if
      if(xmin < 0.000001d0) then
	  write(*,*) ErrorString('lower limit of x-integration is < 10^{-6}. Evaluation stop.',moduleName)
	  stop
      end if
      
      !! in case of cut we determine recut values
      if(doCut) then
        xmin=xMinWithCuts(xmin,var,Cuts)
        xmax=xMaxWithCuts(xmin,var,Cuts)+0.000001d0
      end if
      
      if(xmin<=xmax) then 
        call SetX((xMin+xMax)/2d0,var)
        Xsec_Zint_Xint=Xsec_Zint(var,process,doZ,zMin,zMax)
      else!!! it is possible that cuts cut out the integration range completely
        Xsec_Zint_Xint=0d0
      end if
    end if
    
  end function Xsec_Zint_Xint
  !--------------Simpsons--------------------
  !!!! parameter valueMax remembers the approximate value of integral to weight the tolerance.
  !!!! evaluation is done by adaptive simpson
  !!!! First we evaluate over 5 points and estimate the integral, and then split it to 3+3 and send to adaptive
  !!!! Thus minimal number of points =9
  function integralOverXpoint_S(var,process,doZ,zMin,zMax,yMin_in,yMax_in)
   real(dp),dimension(1:13)::var
   integer,dimension(1:3),intent(in)::process
   real(dp) ::integralOverXpoint_S
   real(dp) :: X1,X2,X3,X4,X5
   real(dp) :: y2,y3,y4,deltay
   real(dp) :: yMin_in,yMax_in
   real(dp)::valueMax
   logical::doZ
   real(dp)::zmin,zMax
   
   deltay=yMax_in-yMin_in
   y2=yMin_in+deltay/4d0
   y3=yMin_in+deltay/2d0
   y4=yMax_in-deltay/4d0
   
   call SetX(yMin_in,var)
   X1= Xsec_Zint(var,process,doZ,zMin,zMax)   
   call SetX(y2,var)
   X2= Xsec_Zint(var,process,doZ,zMin,zMax)   
   call SetX(y3,var)
   X3= Xsec_Zint(var,process,doZ,zMin,zMax)   
   call SetX(y4,var)
   X4= Xsec_Zint(var,process,doZ,zMin,zMax)   
   call SetX(yMax_in,var)
   X5= Xsec_Zint(var,process,doZ,zMin,zMax)   
   
   !!approximate integral value
   valueMax=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   if(methodX=='SA') then
    integralOverXpoint_S=IntegralOverXpoint_S_Rec(var,process,doZ,zMin,zMax,yMin_in,y3,X1,X2,X3,valueMax)+&
	  IntegralOverXpoint_S_Rec(var,process,doZ,zMin,zMax,y3,yMax_in,X3,X4,X5,valueMax)
   else if(methodX=='S5') then
    integralOverXpoint_S=valueMax !!!! 5-point integration
   else
    write(*,*) ErrorString(' integralOverXpoint_S incorrect method. Check methodX variable and initialization',moduleName)
   end if
   
   
  end function integralOverXpoint_S
  
  !!!! X1,X3,X5 are cross-sections at end (X1,X5) and central (X3) points of integraitons
  recursive function integralOverXpoint_S_Rec(var,process,doZ,zMin,zMax,yMin_in,yMax_in,X1,X3,X5,valueMax) result(interX)
   real(dp),dimension(1:13)::var
   integer,dimension(1:3),intent(in)::process
   real(dp) :: interX,X1,X2,X3,X4,X5
   real(dp) :: value,valueAB,valueACB
   real(dp) :: yMin_in,yMax_in,y2,y3,y4,deltay
   real(dp)::valueMax,vv
   logical::doZ
   real(dp)::zmin,zMax
   
   deltay=yMax_in-yMin_in
   y2=yMin_in+deltay/4d0
   y3=yMin_in+deltay/2d0
   y4=yMax_in-deltay/4d0
   
   call SetX(y2,var)
   X2= Xsec_Zint(var,process,doZ,zMin,zMax)   
   
   call SetX(y4,var)
   X4= Xsec_Zint(var,process,doZ,zMin,zMax)   
   
   valueAB=deltay*(X1+4d0*X3+X5)/6d0
   valueACB=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   If(ABS((valueACB-valueAB)/valueMax)>toleranceX) then
    interX=integralOverXpoint_S_Rec(var,process,doZ,zMin,zMax,yMin_in,y3,X1,X2,X3,valueMax)&
	  +integralOverXpoint_S_Rec(var,process,doZ,zMin,zMax,y3,yMax_in,X3,X4,X5,valueMax)
   else
    interX=valueACB
   end if
   
  end function integralOverXpoint_S_Rec
  

  !---------------------------------INTEGRATED over Q (and X and Z)--------------------------------------------------------------
  
  !!! the variable doQ check should the integration over x be performed
  !!! if doQ=true, the integration is done
  !!! if doQ=false the single value (at xMin) is returned
  function Xsec_Zint_Xint_Qint(var,process,doZ,zMin,zMax,doX,xMin,xMax,doQ,Qmin_in,Qmax_in,doCut,Cuts)
    real(dp),dimension(1:13):: var
    logical::doX,doQ,doCut,doZ
    real(dp),dimension(1:4),intent(in)::Cuts
    real(dp) :: Xsec_Zint_Xint_Qint
    real(dp) :: Qmin, Qmax,Qmin_in,Qmax_in,xMin,xMax,zMin,zMax
    integer,dimension(1:3),intent(in)::process
    
    !! the integration over Q is required
    if(doQ) then
      !!! evaluate correspnding y's
      !!! in the case process=2 the integral is over y
      if(process(1)==2) then
	Qmin=QfromSXY(var(11),var(4), Qmin_in)
	Qmax=QfromSXY(var(11),var(4), Qmax_in)
      else
	Qmin=Qmin_in
	Qmax=Qmax_in
      end if
!       write(*,*) '>>>>>',Qmin,Qmax
      !! in case of cut we determine recut values
      if(doCut) then
	if(doX) then
	  Qmin=QMinWithCuts(xmin,Qmin,var,Cuts)
	  Qmax=QMaxWithCuts(xmax,Qmax,var,Cuts)
	else
	  Qmin=QMinWithCuts(xmin,Qmin,var,Cuts)
	  Qmax=QMaxWithCuts(xmin,Qmax,var,Cuts)
	 end if
      end if
      
!       write(*,*) '>>>>>',Qmin,Qmax
      
      if(Qmin<Qmax) then 
	Xsec_Zint_Xint_Qint=integralOverQpoint_S(var,process,doZ,zMin,zMax,doX,xMin,xMax,Qmin,Qmax,doCut,Cuts)
	!Xsec_Zint_Xint_Qint=integralOverQ2point_S(var,process,doZ,zMin,zMax,doX,xMin,xMax,1d0/(Qmax**2),1d0/(Qmin**2),doCut,Cuts)
      else!!! it is possible that cuts cut out the integration range completely
	Xsec_Zint_Xint_Qint=0d0
      end if
    
    else
      ! no integration over Q
      if(process(1)==2) then
        ! y depends on Q linearly, thus there is avarages of transformed is transfromed avarage
	Qmin=QfromSXY(var(11),var(4), Qmin_in)
	Qmax=QfromSXY(var(11),var(4), Qmax_in)
      else
	Qmin=Qmin_in
	Qmax=Qmax_in
      end if
      
      !! in case of cut we determine recut values
      if(doCut) then
	if(doX) then
	  Qmin=QMinWithCuts(xmin,Qmin,var,Cuts)
	  Qmax=QMaxWithCuts(xmax,Qmin,var,Cuts)
	else
	  Qmin=QMinWithCuts((xmin+xmax)/2d0,Qmin,var,Cuts)
	  Qmax=QMaxWithCuts((xmin+xmax)/2d0,Qmin,var,Cuts)+0.000001d0!!! this is needed to resolve 0
	 end if
      end if
      
      if(Qmin<=Qmax) then
	call SetQ((Qmin+Qmax)/2d0,var)
	Xsec_Zint_Xint_Qint=Xsec_Zint_Xint(var,process,doZ,zMin,zMax,doX,xmin,xmax,doCut,Cuts)
      else!!! it is possible that cuts cut out the integration range completely        
	Xsec_Zint_Xint_Qint=0d0
      end if
    end if
  end function Xsec_Zint_Xint_Qint
  !--------------Simpsons--------------------
  !!!! parameter valueMax remembers the approximate value of integral to weight the tolerance.
  !!!! evaluation is done by adaptive simpson
  !!!! First we evaluate over 5 points and estimate the integral, and then split it to 3+3 and send to adaptive
  !!!! Thus minimal number of points =9
  function integralOverQpoint_S(var,process,doZ,zMin,zMax,doX,xMin,xMax,yMin_in,yMax_in,doCut,Cuts)
   real(dp),dimension(1:13)::var
   integer,dimension(1:3),intent(in)::process
   real(dp) ::integralOverQpoint_S
   real(dp) :: X1,X2,X3,X4,X5
   real(dp) :: y2,y3,y4,deltay
   real(dp) :: yMin_in,yMax_in
   real(dp)::valueMax
   logical::doX,doCut,doZ
   real(dp),dimension(1:4),intent(in)::Cuts
   real(dp)::xMin,xMax,zMin,zMax
   
   deltay=yMax_in-yMin_in
   y2=yMin_in+deltay/4d0
   y3=yMin_in+deltay/2d0
   y4=yMax_in-deltay/4d0
   
   call SetQ(yMin_in,var)
   X1=2d0*yMin_in*Xsec_Zint_Xint(var,process,doZ,zMin,zMax,doX,xmin,xmax,doCut,Cuts)  
   call SetQ(y2,var)
   X2=2d0*y2*Xsec_Zint_Xint(var,process,doZ,zMin,zMax,doX,xmin,xmax,doCut,Cuts)  
   call SetQ(y3,var)
   X3=2d0*y3*Xsec_Zint_Xint(var,process,doZ,zMin,zMax,doX,xmin,xmax,doCut,Cuts)  
   call SetQ(y4,var)
   X4=2d0*y4*Xsec_Zint_Xint(var,process,doZ,zMin,zMax,doX,xmin,xmax,doCut,Cuts)  
   call SetQ(yMax_in,var)
   X5=2d0*yMax_in*Xsec_Zint_Xint(var,process,doZ,zMin,zMax,doX,xmin,xmax,doCut,Cuts)  
   
   !!approximate integral value
   valueMax=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   integralOverQpoint_S=IntegralOverQpoint_S_Rec(var,process,doZ,zMin,zMax,doX,xMin,xMax,doCut,Cuts,yMin_in,y3,X1,X2,X3,valueMax)+&
	  IntegralOverQpoint_S_Rec(var,process,doZ,zMin,zMax,doX,xMin,xMax,doCut,Cuts,y3,yMax_in,X3,X4,X5,valueMax)
  end function integralOverQpoint_S
  
  !!!! X1,X3,X5 are cross-sections at end (X1,X5) and central (X3) points of integraitons
  recursive function integralOverQpoint_S_Rec(var,process,doZ,zMin,zMax,doX,xMin,xMax,doCut,Cuts,yMin_in,yMax_in,X1,X3,X5,valueMax)&
												    result(interX)
   real(dp),dimension(1:13)::var
   integer,dimension(1:3),intent(in)::process
   real(dp) :: interX,X1,X2,X3,X4,X5
   real(dp) :: value,valueAB,valueACB
   real(dp) :: yMin_in,yMax_in,y2,y3,y4,deltay
   real(dp)::valueMax,vv
   logical::doX,doCut,doZ
   real(dp),dimension(1:4),intent(in)::Cuts
   real(dp)::xMin,xMax,zMin,zMax
   
   deltay=yMax_in-yMin_in
   y2=yMin_in+deltay/4d0
   y3=yMin_in+deltay/2d0
   y4=yMax_in-deltay/4d0
   
   
   call SetQ(y2,var)
   X2=2*y2*Xsec_Zint_Xint(var,process,doZ,zMin,zMax,doX,xmin,xmax,doCut,Cuts)  
   
   call SetQ(y4,var)
   X4=2*y4*Xsec_Zint_Xint(var,process,doZ,zMin,zMax,doX,xmin,xmax,doCut,Cuts)  
   
   valueAB=deltay*(X1+4d0*X3+X5)/6d0
   valueACB=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   If(ABS((valueACB-valueAB)/valueMax)>tolerance) then
    interX=integralOverQpoint_S_Rec(var,process,doZ,zMin,zMax,doX,xMin,xMax,doCut,Cuts,yMin_in,y3,X1,X2,X3,valueMax)&
	  +integralOverQpoint_S_Rec(var,process,doZ,zMin,zMax,doX,xMin,xMax,doCut,Cuts,y3,yMax_in,X3,X4,X5,valueMax)
   else
    interX=valueACB
   end if
   
  end function integralOverQpoint_S_Rec
  
  
  !--------------Simpsons--------------------
  !!! integration over y=1/Q2  (jacobian=1/2/y^{3/2})
  !!!
  !!!! parameter valueMax remembers the approximate value of integral to weight the tolerance.
  !!!! evaluation is done by adaptive simpson
  !!!! First we evaluate over 5 points and estimate the integral, and then split it to 3+3 and send to adaptive
  !!!! Thus minimal number of points =9
  function integralOverQ2point_S(var,process,doZ,zMin,zMax,doX,xMin,xMax,yMin_in,yMax_in,doCut,Cuts)
   real(dp),dimension(1:13)::var
   integer,dimension(1:3),intent(in)::process
   real(dp) ::integralOverQ2point_S
   real(dp) :: X1,X2,X3,X4,X5
   real(dp) :: y2,y3,y4,deltay
   real(dp) :: yMin_in,yMax_in
   real(dp)::valueMax
   logical::doX,doCut,doZ
   real(dp),dimension(1:4),intent(in)::Cuts
   real(dp)::xMin,xMax,zMin,zMax
   deltay=yMax_in-yMin_in
   y2=yMin_in+deltay/4d0
   y3=yMin_in+deltay/2d0
   y4=yMax_in-deltay/4d0
   
   call SetQ2(1d0/yMin_in,var)
   X1=Xsec_Zint_Xint(var,process,doZ,zMin,zMax,doX,xmin,xmax,doCut,Cuts)/yMin_in**2
   call SetQ2(1d0/y2,var)
   X2=Xsec_Zint_Xint(var,process,doZ,zMin,zMax,doX,xmin,xmax,doCut,Cuts)/y2**2
   call SetQ2(1d0/y3,var)
   X3=Xsec_Zint_Xint(var,process,doZ,zMin,zMax,doX,xmin,xmax,doCut,Cuts)/y3**2
   call SetQ2(1d0/y4,var)
   X4=Xsec_Zint_Xint(var,process,doZ,zMin,zMax,doX,xmin,xmax,doCut,Cuts)/y4**2
   call SetQ2(1d0/yMax_in,var)
   X5=Xsec_Zint_Xint(var,process,doZ,zMin,zMax,doX,xmin,xmax,doCut,Cuts)/yMax_in**2

   !!approximate integral value
   valueMax=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   integralOverQ2point_S=IntegralOverQ2point_S_Rec(var,process,doZ,zMin,zMax,doX,xMin,xMax,doCut,Cuts,yMin_in,y3,X1,X2,X3,valueMax)&
	+IntegralOverQ2point_S_Rec(var,process,doZ,zMin,zMax,doX,xMin,xMax,doCut,Cuts,y3,yMax_in,X3,X4,X5,valueMax)
  end function integralOverQ2point_S
  
  !!!! X1,X3,X5 are cross-sections at end (X1,X5) and central (X3) points of integraitons
  recursive function integralOverQ2point_S_Rec(var,process,doZ,zMin,zMax,doX,xMin,xMax,doCut,Cuts,yMin_in,yMax_in,&
			  X1,X3,X5,valueMax) result(interX)
   real(dp),dimension(1:13) ::var
   integer,dimension(1:3),intent(in)::process
   real(dp) :: interX,X1,X2,X3,X4,X5
   real(dp) :: value,valueAB,valueACB
   real(dp) :: yMin_in,yMax_in,y2,y3,y4,deltay
   real(dp)::valueMax,vv
   logical::doX,doCut,doZ
   real(dp),dimension(1:4),intent(in)::Cuts
   real(dp)::xMin,xMax,zMin,zMax
   
   deltay=yMax_in-yMin_in
   y2=yMin_in+deltay/4d0
   y3=yMin_in+deltay/2d0
   y4=yMax_in-deltay/4d0
      
   call SetQ2(1d0/y2,var)
   X2=Xsec_Zint_Xint(var,process,doZ,zMin,zMax,doX,xmin,xmax,doCut,Cuts)/y2**2
   
   call SetQ2(1d0/y4,var)
   X4=Xsec_Zint_Xint(var,process,doZ,zMin,zMax,doX,xmin,xmax,doCut,Cuts)/y4**2
   
   valueAB=deltay*(X1+4d0*X3+X5)/6d0
   valueACB=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   If(ABS((valueACB-valueAB)/valueMax)>tolerance) then
    interX=integralOverQ2point_S_Rec(var,process,doZ,zMin,zMax,doX,xMin,xMax,doCut,Cuts,yMin_in,y3,X1,X2,X3,valueMax)&
	  +integralOverQ2point_S_Rec(var,process,doZ,zMin,zMax,doX,xMin,xMax,doCut,Cuts,y3,yMax_in,X3,X4,X5,valueMax)
   else
    interX=valueACB
   end if
   
  end function integralOverQ2point_S_Rec
  

  !---------------------------------INTEGRATED over pT (and Z and Q and X)--------------------------------------------------------------
  
  !!! function determines the best value of PT-sections from PT-bin size, and Q
  !!! it is determined by formula Q/PT< val/ (2 k) => def+2K
  function NumPT_auto(dPT,Q)
    real(dp),parameter::val=40d0
    real(dp)::dPT,Q,rat
    integer::i,NumPT_auto
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
  !!! if doZ=true, the integration is done
  !!! if doZ=false the single value (at xMin) is returned
  function Xsec_Zint_Xint_Qint_PTint(var,process,doZ,zMin,zMax,doX,xMin,xMax,doQ,Qmin,Qmax,doPT,ptMin_in,ptMax_in,doCut,Cuts,Num)
    real(dp),dimension(1:13) :: var
    real(dp),dimension(1:4),intent(in) :: Cuts
    logical::doX,doQ,doZ,doPT,doCut
    real(dp) :: Xsec_Zint_Xint_Qint_PTint
    real(dp) :: Qmin,Qmax,xMin,xMax,zMin,zMax,ptMin,ptMax,pT_cur,deltaPT,inter,ptMax_in,ptMin_in
    integer,dimension(1:3),intent(in)::process
    integer::Num,i
    
      if(ptMin_in<0.001d0) then
	pTmin=0.001d0
      else
       pTmin=ptMin_in
      end if
      
      ptMax=ptMax_in
    
    !! the integration over PT is required
    if(doPT) then
      
      if(mod(num,2)>0) then 
        write(*,*) ErrorString('number of Simpson sections is odd. Evaluation stop.',moduleName)
        stop
      end if
      !!!!!!!!!!!!!!!!!!!fixed number Simpsons
      deltaPT=(PTmax-PTmin)/Num
      
      call SetPT(PTmin,var)
      inter=2d0*PTmin*Xsec_Zint_Xint_Qint(var,process,doZ,zMin,zMax,doX,xMin,xMax,doQ,Qmin,Qmax,doCut,Cuts)
      !!!! even terms
      do i=1,Num-1,2
      pT_cur=PTmin+i*deltaPT
      call SetPT(pT_cur,var)
      inter=inter+8d0*pt_cur*Xsec_Zint_Xint_Qint(var,process,doZ,zMin,zMax,doX,xMin,xMax,doQ,Qmin,Qmax,doCut,Cuts)
      
      end do
      
      if(Num>2) then
      !!!! odd terms
      do i=2,Num-2,2
      pT_cur=PTmin+i*deltaPT
      call SetPT(pT_cur,var)
      inter=inter+4d0*pt_cur*Xsec_Zint_Xint_Qint(var,process,doZ,zMin,zMax,doX,xMin,xMax,doQ,Qmin,Qmax,doCut,Cuts)
      end do
      end if
      
      call SetPT(PTmax,var)
      inter=inter+2d0*PTmax*Xsec_Zint_Xint_Qint(var,process,doZ,zMin,zMax,doX,xMin,xMax,doQ,Qmin,Qmax,doCut,Cuts)!!!! last term
      
      
      Xsec_Zint_Xint_Qint_PTint=deltaPT/3d0*inter
      
    
    else    
      ! no integration over PT
      
      call SetPT((PTmin+PTmax)/2d0,var)
      Xsec_Zint_Xint_Qint_PTint=Xsec_Zint_Xint_Qint(var,process,doZ,zMin,zMax,doX,xMin,xMax,doQ,Qmin,Qmax,doCut,Cuts)
      
!       write(*,*) (PTmin+PTmax)/2d0," --> ",Xsec_Zint_Xint_Qint_PTint
    end if
    
  end function Xsec_Zint_Xint_Qint_PTint
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! UNIVERSAL INTERFACE TO THE CROSS-SECTION  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  !----------------------------------------------- single pt point-----------------------------------------
  !!!! just a point cross-section
  subroutine CalcXsecSINGLE_SIDIS(X,pt)
    real(dp),dimension(1:13):: var
    real(dp)::X,pt
    CallCounter=CallCounter+1
    
    var=kinematicArray(pt,s_global,z_global,x_global,Q_global,M2_target_global,M2_product_global)
    
    X=PreFactor1(process_global(1))*xSec(var,process_global)
  end subroutine CalcXsecSINGLE_SIDIS
  
    !!!! just a point cross-section
  subroutine CalcXsecLIST_SIDIS(X_list,pt_list)
    real(dp),dimension(1:13):: var
    real(dp),intent(in)::pt_list(:)
    real(dp),intent(out)::X_list(:)
    integer::length,length2,i
    
    length=size(pt_list)
    
    if(size(X_list)/=length) then
      write(*,*) ErrorString('CalcXsecLIST_SIDIS: sizes of X_list and pt_list are not equal.',moduleName)
      write(*,*) 'Evaluation stop'
      stop
    end if
    
    !$OMP PARALLEL DO DEFAULT(SHARED)
    do i=1,length
      call CalcXsecSINGLE_SIDIS(X_list(i),pt_list(i))
     end do
    !$OMP END PARALLEL DO    
    
  end subroutine CalcXsecLIST_SIDIS
  
  !----------------------------------------------- XQZ-integerations -----------------------------------------
  
  !!! integrated over XQZ
  subroutine CalcXsecSINGLE_SIDIS_Zint_Xint_Qint(X,pt,zMin,zMax,xMin,xMax,Qmin,Qmax)
    real(dp)::X,pt,xMin,xMax,Qmin,Qmax,zMin,zMax
    real(dp),dimension(1:13):: var
    CallCounter=CallCounter+1
    var=kinematicArray(pt,s_global,(zMin+zMax)/2d0,(xMin+xMax)/2d0,(Qmin+Qmax)/2d0,M2_target_global,M2_product_global)
    X=PreFactor1(process_global(1))*Xsec_Zint_Xint_Qint(var,process_global,.true.,zMin,zMax,.true.,xMin,xMax,.true.,Qmin,Qmax,&
      includeCuts_global,CutParameters_global)
  
  end subroutine CalcXsecSINGLE_SIDIS_Zint_Xint_Qint
  
  !!! integrated over XQZ
  subroutine CalcXsecLIST_SIDIS_Zint_Xint_Qint(X_list,pt_list,zMin,zMax,xMin,xMax,Qmin,Qmax)
    real(dp),intent(in)::pt_list(:)
    real(dp),intent(in)::xMin,xMax,Qmin,Qmax,zMin,zMax
    real(dp),intent(out)::X_list(:)
    integer::i,length
    
    length=size(pt_list)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!sizes checks
    if(size(X_list)/=length) then
      write(*,*) ErrorString('CalcXsecLIST_SIDIS_Zint_Qint_Xint: sizes of X_list and pt_list are not equal.',moduleName)
      write(*,*) 'Evaluation stop'
      stop
    end if
    
    
    !$OMP PARALLEL DO DEFAULT(SHARED)
    do i=1,length
      call CalcXsecSINGLE_SIDIS_Zint_Xint_Qint(X_list(i),pt_list(i),zMin,zMax,xMin,xMax,Qmin,Qmax)
    end do
    !$OMP END PARALLEL DO    
  
  end subroutine CalcXsecLIST_SIDIS_Zint_Xint_Qint

  
!----------------------------------------------- PT+XQZ-integerations -----------------------------------------
    !!! this function just to help incapsulte the variables for parallel computation
  function CalcXsecHELP_PTint_Zint_Xint_Qint(ptMin,ptMax,zMin,zMax,xMin,xMax,Qmin,Qmax)
    real(dp)::pt,xMin,xMax,Qmin,Qmax,zMin,zMax,ptMin,ptMax,CalcXsecHELP_PTint_Zint_Xint_Qint
    real(dp),dimension(1:13):: var
    integer::Num
    var=kinematicArray((ptMin+ptMax)/2d0,s_global,(zMin+zMax)/2d0,(xMin+xMax)/2d0,(Qmin+Qmax)/2d0,&
		    M2_target_global,M2_product_global)
		    
    Num=NumPT_auto(ptMax-ptMin,var(2))
  
    CalcXsecHELP_PTint_Zint_Xint_Qint=PreFactor1(process_global(1))*&
	     Xsec_Zint_Xint_Qint_PTint(var,process_global,.true.,zMin,zMax,.true.,xMin,xMax,.true.,Qmin,Qmax,.true.,ptMin,ptMax,&
		includeCuts_global,CutParameters_global,Num)
  end function CalcXsecHELP_PTint_Zint_Xint_Qint

    !!! integrated over PT+XQZ
  subroutine CalcXsecSINGLE_SIDIS_PTint_Zint_Xint_Qint(X,ptMin,ptMax,zMin,zMax,xMin,xMax,Qmin,Qmax)
    real(dp)::X,pt,xMin,xMax,Qmin,Qmax,zMin,zMax,ptMin,ptMax
    
    CallCounter=CallCounter+1
    X=CalcXsecHELP_PTint_Zint_Xint_Qint(ptMin,ptMax,zMin,zMax,xMin,xMax,Qmin,Qmax)
  
  end subroutine CalcXsecSINGLE_SIDIS_PTint_Zint_Xint_Qint
  
  subroutine CalcXsecLIST_SIDIS_PTint_Zint_Xint_Qint(X_list,ptMin_list,ptMax_list,zMin,zMax,xMin,xMax,Qmin,Qmax)
    real(dp)::xMin,xMax,Qmin,Qmax,zMin,zMax
    real(dp),intent(in)::ptMax_list(:),ptMin_list(:)
    real(dp)::X_list(:)
    integer::length,i
    
    length=size(ptMin_list)
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!sizes checks
    if(size(X_list)/=length) then
      write(*,*) ErrorString('CalcXsecLIST_SIDIS_PTint_Zint_Qint_Xint: sizes of X_list and ptMin_list are not equal.',moduleName)
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(ptMax_list)/=length) then
      write(*,*) ErrorString('CalcXsecLIST_SIDIS_PTint_Zint_Qint_Xint: sizes of ptMax_list and ptMin_list are not equal.'&
		,moduleName)
      write(*,*) 'Evaluation stop'
      stop
    end if
    
    CallCounter=CallCounter+length
    !$OMP PARALLEL DO DEFAULT(SHARED)
    do i=1,length
    X_list(i)=CalcXsecHELP_PTint_Zint_Xint_Qint(ptMin_list(i),ptMax_list(i),zMin,zMax,xMin,xMax,Qmin,Qmax)
    end do
    !$OMP END PARALLEL DO
  
  end subroutine CalcXsecLIST_SIDIS_PTint_Zint_Xint_Qint

  
  subroutine CalcXsecLISTLIST_SIDIS_PTint_Zint_Xint_Qint(X_list,pt_list,zMin,zMax,xMin,xMax,Qmin,Qmax)
    real(dp)::xMin,xMax,Qmin,Qmax,zMin,zMax,ptMin,ptMax
    real(dp),intent(in)::pt_list(:)
    real(dp)::X_list(:)
    integer::length,i
    
    length=size(pt_list)-1
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!sizes checks
    if(size(X_list)/=length) then
      write(*,*) ErrorString('CalcXsecLIST_SIDIS_PTint_Zint_Qint_Xint: sizes of X_list and ptMin_list are not equal.',moduleName)
      write(*,*) 'Evaluation stop'
      stop
    end if
    
    CallCounter=CallCounter+length
    !$OMP PARALLEL DO DEFAULT(SHARED)
    do i=1,length
    X_list(i)=CalcXsecHELP_PTint_Zint_Xint_Qint(pt_list(i),pt_list(i+1),zMin,zMax,xMin,xMax,Qmin,Qmax)
    end do
    !$OMP END PARALLEL DO
  
  end subroutine CalcXsecLISTLIST_SIDIS_PTint_Zint_Xint_Qint
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!  MAIN INTERFACE !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  !!! single value interface
  subroutine xSec_SIDIS(xx,process,s,pT,z,x,Q,doCut,Cuts,masses)
    integer,intent(in),dimension(1:3)::process			!the number of process
    real(dp),intent(in)::s					!Mandelshtam s
    real(dp),intent(in),dimension(1:2)::pT			!(qtMin,qtMax)
    real(dp),intent(in),dimension(1:2)::z				!(zmin,zmax)
    real(dp),intent(in),dimension(1:2)::x				!(xmin,xmax)
    real(dp),intent(in),dimension(1:2)::Q				!(Qmin,Qmax)    
    logical,intent(in)::doCut					!triger cuts
    real(dp),intent(in),dimension(1:4)::Cuts			!(ymin,yMax,W2)
    real(dp),intent(in),dimension(1:2),optional::masses		!(mass_target,mass-product)GeV
    real(dp),intent(out)::xx
    integer :: i,length
    
    real(dp),dimension(1:13):: var
    integer::Num
    
    CallCounter=CallCounter+1
    if(PRESENT(masses)) then
    var=kinematicArray((pt(1)+pt(2))/2d0,s,(z(1)+z(2))/2d0,(x(1)+x(2))/2d0,(Q(1)+Q(2))/2d0,masses(1)**2,masses(2)**2)
    else
    var=kinematicArray((pt(1)+pt(2))/2d0,s,(z(1)+z(2))/2d0,(x(1)+x(2))/2d0,(Q(1)+Q(2))/2d0,M2_target_global,M2_product_global)
    end if
    
    
    Num=NumPT_auto(pt(2)-pt(1),var(2))
  
    xx=PreFactor1(process(1))*&
	     Xsec_Zint_Xint_Qint_PTint(var,process,.true.,z(1),z(2),.true.,x(1),x(2),.true.,Q(1),Q(2),.true.,pt(1),pt(2),&
			    doCut,Cuts,Num)
    
    
  end subroutine xSec_SIDIS
  
  subroutine xSec_SIDIS_List(xx,process,s,pT,z,x,Q,doCut,Cuts,masses)
    integer,intent(in),dimension(:,:)::process			!the number of process
    real(dp),intent(in),dimension(:)::s				!Mandelshtam s
    real(dp),intent(in),dimension(:,:)::pT			!(qtMin,qtMax)
    real(dp),intent(in),dimension(:,:)::z				!(zmin,zmax)
    real(dp),intent(in),dimension(:,:)::x				!(xmin,xmax)
    real(dp),intent(in),dimension(:,:)::Q				!(Qmin,Qmax)        
    logical,intent(in),dimension(:)::doCut			!triger cuts
    real(dp),intent(in),dimension(:,:)::Cuts			!(ymin,yMax,W2)
    real(dp),intent(in),dimension(:,:),optional::masses		!(mass_target,mass-product)GeV
    real(dp),dimension(:),intent(out)::xx
    integer :: i,length
    
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
    if(size(process,2)/=3) then
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
    xx(i)=xSecFULL(process(i,1:3),s(i),pt(i,1),pt(i,2),z(i,1),z(i,2),x(i,1),x(i,2),Q(i,1),Q(i,2),doCut(i),Cuts(i,1:4),&
		    masses(i,1)**2,masses(i,2)**2)
    end do
    !$OMP END PARALLEL DO
    
    else
    
    !$OMP PARALLEL DO DEFAULT(SHARED)
    do i=1,length
    xx(i)=xSecFULL(process(i,1:3),s(i),pt(i,1),pt(i,2),z(i,1),z(i,2),x(i,1),x(i,2),Q(i,1),Q(i,2),doCut(i),Cuts(i,1:4),&
		    M2_target_global,M2_product_global)
    end do
    !$OMP END PARALLEL DO
    
    end if
    
  end subroutine xSec_SIDIS_List
  
  
  !!!! problem is that f2py does not like optional arguments.. in any form
  subroutine xSec_SIDIS_List_forharpy(xx,process,s,pT,z,x,Q,doCut,Cuts,masses)
    integer,intent(in),dimension(:,:)::process			!the number of process
    real(dp),intent(in),dimension(:)::s				!Mandelshtam s
    real(dp),intent(in),dimension(:,:)::pT			!(qtMin,qtMax)
    real(dp),intent(in),dimension(:,:)::z				!(zmin,zmax)
    real(dp),intent(in),dimension(:,:)::x				!(xmin,xmax)
    real(dp),intent(in),dimension(:,:)::Q				!(Qmin,Qmax)        
    logical,intent(in),dimension(:)::doCut			!triger cuts
    real(dp),intent(in),dimension(:,:)::Cuts			!(ymin,yMax,W2)
    real(dp),intent(in),dimension(:,:)::masses		!(mass_target,mass-product)GeV
    real(dp),dimension(:),intent(out)::xx
    integer :: i,length
    
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
    if(size(process,2)/=3) then
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
    
    xx(i)=xSecFULL(process(i,1:3),s(i),pt(i,1),pt(i,2),z(i,1),z(i,2),x(i,1),x(i,2),Q(i,1),Q(i,2),doCut(i),Cuts(i,1:4),&
		    masses(i,1)**2,masses(i,2)**2)
    end do
    !$OMP END PARALLEL DO
  end subroutine xSec_SIDIS_List_forharpy
  
  
  !!! helper to incapsulate PARALLEL variables
  function xSecFULL(proc,s,ptmin,ptmax,zmin,zmax,xmin,xmax,Qmin,Qmax,doCut,Cuts,m1,m2)
  real(dp)::s,ptmin,ptmax,zmin,zmax,xmin,xmax,Qmin,Qmax,Cuts(1:4),xSecFULL,var(1:13),m1,m2
  integer::proc(1:3),Num
  logical::doCut
  
  var=kinematicArray((ptmin+ptmax)/2d0,s,(zmin+zmax)/2d0,(xmin+xmax)/2d0,(Qmin+Qmax)/2d0,m1,m2)
  Num=NumPT_auto(ptmax-ptmin,var(2))
  
!   write(*,*) 'aTMD:1  ',var
!   write(*,*) 'aTMD:2  ',proc,m1,m2
!   write(*,*) 'aTMD:3  ',.true.,zmin,zmax,.true.,xmin,xmax,.true.,Qmin,Qmax,.true.,ptmin,ptmax
  
  xSecFULL=PreFactor1(proc(1))*Xsec_Zint_Xint_Qint_PTint(var,proc,&
		    .true.,zmin,zmax,.true.,xmin,xmax,.true.,Qmin,Qmax,.true.,ptmin,ptmax,doCut,Cuts,Num)
  end function xSecFULL

  !!!! problem is that f2py does not like optional arguments.. in any form
  subroutine xSec_SIDIS_BINLESS_List_forharpy(xx,process,s,pT,z,x,Q,masses)
    integer,intent(in),dimension(:,:)::process			!the number of process
    real(dp),intent(in),dimension(:)::s				!Mandelshtam s
    real(dp),intent(in),dimension(:)::pT			!(pt)
    real(dp),intent(in),dimension(:)::z				!(z)
    real(dp),intent(in),dimension(:)::x				!(x)
    real(dp),intent(in),dimension(:)::Q				!(Q)        
    real(dp),intent(in),dimension(:,:)::masses		!(mass_target,mass-product)GeV
    real(dp),dimension(:),intent(out)::xx
    integer :: i,length
    
    length=size(s)
    CallCounter=CallCounter+length
    
    !!! cheking sizes
    if(size(xx)/=length) then
      write(*,*) ErrorString('xSec_SIDIS_BINLESS_List: sizes of xSec and s lists are not equal.',moduleName)
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(process,1)/=length) then
      write(*,*) ErrorString('xSec_SIDIS_BINLESS_List: sizes of process and s lists are not equal.',moduleName)
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(pT)/=length) then
      write(*,*) ErrorString('xSec_SIDIS_BINLESS_List: sizes of pT and s lists are not equal.',moduleName)
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(x)/=length) then
      write(*,*) ErrorString('xSec_SIDIS_BINLESS_List: sizes of x and s lists are not equal.',moduleName)
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(Q)/=length) then
      write(*,*) ErrorString('xSec_SIDIS_BINLESS_List: sizes of Q and s lists are not equal.',moduleName)
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(z)/=length) then
      write(*,*) ErrorString('xSec_SIDIS_BINLESS_List: sizes of z and s lists are not equal.',moduleName)
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(process,2)/=3) then
      write(*,*) ErrorString('xSec_SIDIS_BINLESS_List: process list must be (:,1:3).',moduleName)
      write(*,*) 'Evaluation stop'
      stop
    end if
    
   CallCounter=CallCounter+length
   
   if(size(masses,1)/=length) then
      write(*,*) ErrorString('xSec_SIDIS_BINLESS_List: sizes of masses and s lists are not equal.',moduleName)
      write(*,*) 'Evaluation stop'
      stop
    end if
        if(size(masses,2)/=2) then
      write(*,*) ErrorString('xSec_SIDIS_BINLESS_List: mass list must be (:,1:2).',moduleName)
      write(*,*) 'Evaluation stop'
      stop
    end if
   
   !$OMP PARALLEL DO DEFAULT(SHARED)
    do i=1,length
    xx(i)=xSecBINLESS(process(i,1:3),s(i),pt(i),z(i),x(i),Q(i),masses(i,1)**2,masses(i,2)**2)
    end do
    !$OMP END PARALLEL DO
  end subroutine xSec_SIDIS_BINLESS_List_forharpy
  

  !!! helper to incapsulate PARALLEL variables
  !!! evaluate cross-section without integrations over bins, and without cuts
  function xSecBINLESS(proc,s,pt,z,x,Q,m1,m2)
  real(dp)::s,pt,z,x,Q,xSecBINLESS,var(1:13),m1,m2
  integer::proc(1:3)
  
  var=kinematicArray(pt,s,z,x,Q,m1,m2)
  
!    write(*,*) 'aTMD:1  ',var
!    write(*,*) 'aTMD:2  ',proc,m1,m2
!    write(*,*) 'aTMD:3  ',.false.,z,z,.false.,x,x,.false.,Q,Q,.false.,pt,pt
  
  xSecBINLESS=PreFactor1(proc(1))*Xsec_Zint_Xint_Qint_PTint(var,proc,&
		    .false.,z,z,.false.,x,x,.false.,Q,Q,.false.,pt,pt,.false.,(/0d0,1d0,0d0,1d9/),4)
  end function xSecBINLESS
  
end module TMDX_SIDIS
