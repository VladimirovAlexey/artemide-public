!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 3.00
!
!    Evaluation of the TMD cross-section for DY-like cross-sections
!
!    if you use this module please, quote 1706.01473
!
!    ver 3.0: created from v.2.06 (AV, 05.09.2023)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDX_DY
use aTMDe_Numerics
use IO_functions
use TMDF
use TMDF_KPC
use LeptonCutsDY
use QCDinput
use EWinput
use IntegrationRoutines

implicit none
private

!!!!!! 1=accurate but slow
!!!!!! 2=fast but not accurate
#define INTEGRATION_MODE 2

!Current version of module
character (len=7),parameter :: moduleName="TMDX-DY"
character (len=5),parameter :: version="v3.00"
!Last appropriate verion of constants-file
integer,parameter::inputver=30

real(dp) :: toleranceINT=0.0001d0
real(dp) :: toleranceGEN=0.0000001d0

integer::outputlevel
integer::messageTrigger

!!!!
!!!! in the module the kinematic is stored in the varibles "kinematic" real(dp),dimension(1:6)
!!!! which is (qT,s,Q,Q^2,x0,y,exp[y])
!!!! where x0=sqrt[(Q^2+q_T^2)/s]   (if exactX1X2) or x0=Q^2/s (otherwise)
!!!!
!!other global parameters see SetXParameters
integer:: orderH_global
logical::useKPC
logical::usePIresum
integer:: exactX1X2    !!!=1 if exact x's=true, =0 otherwise
integer:: exactScales  !!!=1 if exact hard scales = true, =0 otherwise

!!! number of sections for PT-integral by default
integer::NumPTdefault=6
!!! Maximum size of Q-bin. Larger bins are desected
real::maxQbinSize=30.
!!! Minimal qT, below this number the value is frozen
real(dp)::qTMin_ABS=0.0001d0

real(dp)::c2_global!,muHard_global

integer::GlobalCounter
integer::CallCounter
integer::messageCounter

real(dp)::hc2

logical::started=.false.

public::TMDX_DY_Initialize,TMDX_DY_SetScaleVariation,&
TMDX_DY_ShowStatistic,TMDX_DY_ResetCounters,TMDX_DY_IsInitialized
public::  xSec_DY,xSec_DY_List,xSec_DY_List_BINLESS,xSec_DY_List_APPROXIMATE

interface xSec_DY
  module procedure MainInterface_AsAAAloo
end interface

contains
  
function TMDX_DY_IsInitialized()
  logical::TMDX_DY_IsInitialized
  TMDX_DY_IsInitialized=started
end function TMDX_DY_IsInitialized

!! Initialization of the package
subroutine TMDX_DY_Initialize(file,prefix)
  character(len=*),intent(in)::file
  character(len=*),intent(in),optional::prefix
  character(len=300)::path
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

  !!! Check the file version
  call MoveTO(51,'*0   ')
  call MoveTO(51,'*A   ')
  call MoveTO(51,'*p1  ')
  read(51,*) FILEver
  if(FILEver<inputver) then
    write(*,*) 'artemide.'//trim(moduleName)//': const-file version is too old.'
    write(*,*) '             Update the const-file with artemide.setup'
    write(*,*) '  '
    CLOSE (51, STATUS='KEEP')
    stop
  end if

  !!! Fill the message system
  call MoveTO(51,'*p2  ')
  read(51,*) outputLevel
  if(outputLevel>2) write(*,*) '--------------------------------------------- '
  if(outputLevel>2) write(*,*) 'artemide.',moduleName,version,': initialization started ... '
  call MoveTO(51,'*p3  ')
  read(51,*) messageTrigger

  !!! other variables
  call MoveTO(51,'*B   ')
  call MoveTO(51,'*p1  ')
  read(51,*) hc2

  !$    if(outputLevel>1) write(*,*) '    artemide.TMDX_DY: parallel evaluation of cross-sections is to be used'
  !$    call MoveTO(51,'*C   ')
  !$    call MoveTO(51,'*p1  ')
  !$    read(51,*) i
  !$    call OMP_set_num_threads(i)
  !$    if(outputLevel>1) write(*,*) '    artemide.TMDX_DY: number of threads for parallel evaluation is set to ', i

  !!! I also need MZ

  !!! go to section TMD-DY
  call MoveTO(51,'*9   ')
  call MoveTO(51,'*p1  ')
  read(51,*) initRequired
  if(.not.initRequired) then
    if(outputLevel>2) write(*,*)'artemide.',moduleName,': initialization is not required. '
    started=.false.
    CLOSE (51, STATUS='KEEP')
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
    CASE ("N3LO") !!! same as NNNLO
    orderH_global=3
    CASE ("N4LO")
    orderH_global=4
    CASE DEFAULT
    if(outputLevel>0) write(*,*)  WarningString('try to set unknown order. Switch to NLO.',moduleName)
    orderH_global=1
  END SELECT
  if(outputLevel>1) write(*,*) '    artemide.TMDX_DY: the used order is ',trim(orderMain)

  call MoveTO(51,'*p2   ')
  read(51,*) useKPC
  if(outputLevel>1 .and. useKPC) write(*,*) color('    artemide.TMDX_DY: using TMD factorization with KPC',c_cyan)
  if(outputLevel>1 .and. useKPC) write(*,*) color('                      Please, cite [2307.13054]',c_cyan)
  if(outputLevel>1 .and. .not.(useKPC)) write(*,*) color('    artemide.TMDX_DY: using TMD factorization at LP',c_cyan)

  !! pi2 resummation
  call MoveTO(51,'*p3   ')
  read(51,*) usePIresum
  if(outputLevel>2 .and. usePIresum) write(*,*) '    artemide.TMDX_DY: pi-resummation in coef.function included.'

  !!!------ parameters of numerics
  call MoveTO(51,'*B   ')
  call MoveTO(51,'*p1  ')
  read(51,*) toleranceGEN
  call MoveTO(51,'*p2  ')
  read(51,*) toleranceINT
  call MoveTO(51,'*p3  ')
  read(51,*) NumPTdefault
  call MoveTO(51,'*p4  ')
  read(51,*) maxQbinSize
  call MoveTO(51,'*p5  ')
  read(51,*) qTMin_ABS

   if(.not.(useKPC)) then
    !!!------ parameters of LP factorization
    call MoveTO(51,'*C   ')
    !!exact values of x1x2
    call MoveTO(51,'*p1   ')
    read(51,*) dummyLogical
    if(dummyLogical) then
      exactX1X2=1
    else
      exactX1X2=0
    end if
    if(outputLevel>2 .and. dummyLogical) write(*,*) '    artemide.TMDX_DY: qT/Q corrections for x1 and x2 variables are included.'
    !!exact values for scales
    call MoveTO(51,'*p2   ')
    read(51,*) dummyLogical
    if(dummyLogical) then
      exactScales=1
    else
      exactScales=0
    end if
    if(outputLevel>2 .and. dummyLogical) write(*,*) '    artemide.TMDX_DY: qT/Q correction for scales variables are included.'

  else
    !!!------ parameters of KPC factorization
    exactX1X2=1
    exactScales=0

    call MoveTO(51,'*D   ')
    !!!!! nothing is here yet

  end if

  CLOSE (51, STATUS='KEEP')


  !$     if(outputLevel>2) write(*,*) '------TEST OF PARALLEL PROCESSING ----------'
  !$OMP PARALLEL
  !$     if(outputLevel>2) write(*,*) '   artemide.TMDX_DY:thread num ',  omp_get_thread_num(), ' ready.'
  !$OMP END PARALLEL

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

  !!!!! initializing Lepton Cut module
  call InitializeLeptonCutDY(toleranceINT,toleranceGEN)

  c2_global=1d0

  GlobalCounter=0
  CallCounter=0
  messageCounter=0

  started=.true.

#if INTEGRATION_MODE==2
    write(*,*)  color('--------------------------------------------------------',c_red)
    write(*,*)  color('----------------------  WARNING!  ----------------------',c_red)
    write(*,*)  color('--   TMDX_DY is in the approximate integration mode   --',c_red)
    write(*,*)  color('--            Faster, but lower precision.            --',c_red)
    write(*,*)  color('--    Switch to default version by changing flag      --',c_red)
    write(*,*)  color('--   INTEGRATION_MODE in TMDX_DY.f90, and recompile   --',c_red)
    write(*,*)  color('--------------------------------------------------------',c_red)
#endif

  write(*,*)  color('----- arTeMiDe.TMD_DY '//trim(version)//': .... initialized',c_green)
end subroutine TMDX_DY_Initialize

subroutine TMDX_DY_ResetCounters()
  if(outputlevel>2) call TMDX_DY_ShowStatistic()
  GlobalCounter=0
  CallCounter=0
  messageCounter=0
end subroutine TMDX_DY_ResetCounters
  
subroutine TMDX_DY_ShowStatistic()

    write(*,'(A,ES12.3)') 'TMDX DY statistics      total calls of point xSec  :  ',Real(GlobalCounter)
    write(*,'(A,ES12.3)') '                              total calls of xSecF :  ',Real(CallCounter)
    write(*,'(A,F12.3)')  '                                         avarage M :  ',Real(GlobalCounter)/Real(CallCounter)
end subroutine TMDX_DY_ShowStatistic
  
  
  
!!!!Call this after TMD initializetion but before NP, and X parameters
subroutine TMDX_DY_SetScaleVariation(c2_in)
  real(dp),intent(in)::c2_in

  if(outputLevel>1) write(*,*) 'TMDX_DY: c2 scale reset:',c2_in

  if(c2_in<0.1d0 .or. c2_in>10.d0) then
    if(outputLevel>0) write(*,*) WarningString('variation in c2 is enourmous. c2 is set to 2',moduleName)
      c2_global=2d0
  else
    c2_global=c2_in
  end if

end subroutine TMDX_DY_SetScaleVariation

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS FOR OPERATION WITH KINEMATICS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!!! function makes kinematic array from the given set of qT,s,Q,y
!!! array has 7 often appearing entries
function kinematicArray(qT,s,Q,y)
  real(dp),dimension(1:7)::kinematicArray
  real(dp)::qT,s,Q,y

  kinematicArray=(/qT,s,Q,Q**2,sqrt((Q**2+exactX1X2*qT**2)/s),y,exp(y)/)

end function kinematicArray
  
!!!intrinsic change the value of Q within kinematic array var
subroutine SetQ(Q,var)
  real(dp),dimension(1:7)::var
  real(dp)::Q

  var(3)=Q
  var(4)=Q**2
  var(5)=sqrt((var(4)+exactX1X2*var(1)**2)/var(2))

end subroutine SetQ
  
!!!intrinsic change the value of y within kinematic array var
subroutine SetY(y,var)
  real(dp),dimension(1:7)::var
  real(dp)::y

  var(6)=y
  var(7)=exp(y)

end subroutine SetY
  
!!!intrinsic change the value of qT within kinematic array var
subroutine SetQT(qT_in,var)
  real(dp),dimension(1:7)::var
  real(dp)::qT_in

  var(1)=qT_in
  var(5)=sqrt((var(4)+exactX1X2*var(1)**2)/var(2))

end subroutine SetQT


!!! Computes y-variable from x_F
pure function yFromXF(xF,var)
  real(dp),intent(in),dimension(1:7)::var
  real(dp),intent(in):: xF
  real(dp):: yFromXF
  yFromXF=asinh(xF/2d0/var(5))
end function yFromXF


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS FOR PREFACTORS & PROCESSES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! PROCESS DEFINITION:
!!! (i, h1, h2, n1, n2, ...) -- integers
!!! i = the integration type and global prefactor
!!! h1,h2 = the hadron types. h>0
!!! n = the integrand combination
!!! several n's implies summation of corresponding process
!!! n=0 is skept in the summation
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! for easier read coeff-functions are split into separate file
INCLUDE '/Code/TMDX/DYcoeff-func.f90'

!!!!!----------------------------------------------------------------
!!!!!  Prefactors and cuts for LP computation
!!!!!----------------------------------------------------------------

!!!!! Prefactor 2 is (universal part) x (cuts) x H
!!!!! This depends only on the integration type of the process (i)
function PreFactor2(kin,process)
  real(dp),dimension(1:7),intent(in)::kin
  integer,intent(in)::process
  real(dp)::PreFactor2

  real(dp)::uniPart,scaleMu



  !!!! universal part
  scaleMu=sqrt(kin(4)+exactScales*kin(1)**2)

  SELECT CASE(process)
  case(0)
    uniPart=1d0
  CASE(1)
    !4 pi aEm^2/3 /Nc/Q^2/s
    uniPart=pix4/9d0*(alphaEM(scaleMu)**2)/(kin(2)*kin(4))*&
        HardCoefficientDY(scaleMu)*&
        hc2*1d9!from GeV to pb
  CASE(2)
    !4 pi aEm^2/3 /Nc/Q^2/s
    ! the process=2 is for the xF-integration. It has extra weigth 2sqrt[(Q^2+q_T^2)/s] Cosh[y]
    uniPart=pix4/9d0*(alphaEM(scaleMu)**2)/(kin(2)*kin(4))*&
        HardCoefficientDY(scaleMu)*&
        hc2*1d9*&!from GeV to pb
        2._dp*kin(5)*cosh(kin(6))

  CASE (3) !Zboson in the narrow-width approximation
    !4 pi^2 aem/Nc/s Br(z->ee+mumu)
    uniPart=pi2x4/3d0*alphaEM(scaleMu)/kin(2)*&
        HardCoefficientDY(scaleMu)*&
        hc2*1d9*&!from GeV to pb
        0.03645d0!Br from PDG, ee+mumu
  CASE (4) !Wboson in the narrow-width approximation
    !4 pi^2 aem/Nc/s Br(z->ee+mumu)
    uniPart=pi2x4/3d0**alphaEM(scaleMu)/kin(2)*&
        HardCoefficientDY(scaleMu)*&
        hc2*1d9*&!from GeV to pb
        0.1086d0!Br from PDG, ee+mumu
  CASE (5) !exclusive HIGGSboson production
    ! (2\pi) *pi Mh^2 as(mu)/36/s/vev^2 * H*cT^2
    ! (1.033)^2 is correction for mT mass in Ct at LO.
    uniPart=(1d0/18d0)*MH2*(As(c2_global*scaleMu)/VEVH)**2/kin(2)*&
        HardCoefficientHIGGS(scaleMu)*(EffCouplingHFF(scaleMu)**2)*1.0677023627519822d0*&
        hc2*1d9!from GeV to pb
  CASE DEFAULT
    write(*,*) ErrorString('unknown process p='//numToStr(process)//' .Evaluation stop.',moduleName)
    stop
  END SELECT

  PreFactor2=uniPart

end function PreFactor2

!!!! Fiducial cut for leptons.
!!!! this factor depends on the type of process and computed in the corresponding module
!!!! This function perform the preliminary selection depending on the process
function LeptonCutFactorLP(kin,proc1, includeCuts_in,CutParam)
  real(dp),dimension(1:7),intent(in)::kin
  logical,intent(in)::includeCuts_in
  real(dp),dimension(1:4),intent(in)::CutParam
  integer,intent(in)::proc1
  real(dp)::LeptonCutFactorLP

  !!!!! lepton-cut prefactor
  if(includeCuts_in) then
    !!! here include cuts on the lepton tensor
    LeptonCutFactorLP=CutFactor(qT_in=kin(1),Q_in=kin(3),y_in=kin(6),CutParameters=CutParam,Cut_Type=-2)
    !LeptonCutFactorLP=CutFactor(qT_in=kin(1),Q_in=kin(3),y_in=kin(6),CutParameters=CutParam,Cut_Type=-1)
  else
    !!! this is uncut lepton tensor
    LeptonCutFactorLP=(1+0.5d0*(kin(1)/kin(3))**2)
  end if

end function LeptonCutFactorLP

!!!!!----------------------------------------------------------------
!!!!!  Prefactors and cuts for KPC computation
!!!!!----------------------------------------------------------------

!!!!! Prefactor 2 is (universal part)  x H
!!!!! This prefactor depends only on the general class of process
!!!!! proc1 is the first number of process array
!!!!! Lepton cut factor is defined in a separate function and depend on the process
function PreFactorKPC(kin,proc1)
  real(dp),dimension(1:7),intent(in)::kin
  integer,intent(in)::proc1
  real(dp)::PreFactorKPC

  real(dp)::scaleMu

  !!!! universal part
  scaleMu=kin(3)

  SELECT CASE(proc1)
  case(0)
    PreFactorKPC=1d0
  CASE(1)
!     !4 pi aEm^2/3 /Nc/s
    PreFactorKPC=pi2x2/9*(alphaEM(scaleMu)**2)/kin(2)*&
        HardCoefficientDY(scaleMu)*&
        hc2*1d9!from GeV to pb

  CASE(2)
    !4 pi aEm^2/3 /Nc/s
    ! the process=2 is for the xF-integration. It has extra weigth 2sqrt[(Q^2+q_T^2)/s] Cosh[y]
    PreFactorKPC=pi2x2/9*(alphaEM(scaleMu)**2)/kin(2)*&
        HardCoefficientDY(scaleMu)*&
        hc2*1d9*&!from GeV to pb
        2._dp*kin(5)*cosh(kin(6))

  CASE (5) !exclusive HIGGSboson production
    ! (2\pi) *pi Mh^2 as(mu)/36/s/vev^2 * H*cT^2
    ! (1.033)^2 is correction for mT mass in Ct at LO.
    ! * pi*Q2/2 (because of KPC ??? this factor is assumption)
    PreFactorKPC=(pi/36d0)*MH2*(As(c2_global*scaleMu)/VEVH)**2*kin(4)/kin(2)*&
        HardCoefficientHIGGS(scaleMu)*(EffCouplingHFF(scaleMu)**2)*1.0677023627519822d0*&
        hc2*1d9!from GeV to pb

  CASE DEFAULT
    write(*,*) ErrorString('unknown process p='//numToStr(proc1)//' .Evaluation stop.',moduleName)
    stop
  END SELECT

end function PreFactorKPC

!!!! Fiducial cut for leptons.
!!!! this factor depends on the type of process and computed in the corresponding module
!!!! This function perform the preliminary selection depending on the process
function LeptonCutFactorKPC(kin,proc1, includeCuts_in,CutParam)
  real(dp),dimension(1:7),intent(in)::kin
  logical,intent(in)::includeCuts_in
  real(dp),dimension(1:4),intent(in)::CutParam
  integer,intent(in)::proc1
  real(dp)::LeptonCutFactorKPC

  if(includeCuts_in) then
      !!! here include cuts onf lepton tensor
  !!!!!!!!!!!!!!!
  SELECT CASE(proc1)
  CASE(1,2,101,102,103) !!! UU
      LeptonCutFactorKPC=CutFactor(qT_in=kin(1),Q_in=kin(3),y_in=kin(6),CutParameters=CutParam,Cut_Type=-1)

  CASE(20,30) !!! Angular coefficients A0
      LeptonCutFactorKPC=CutFactor(qT_in=kin(1),Q_in=kin(3),y_in=kin(6),CutParameters=CutParam,Cut_Type=0)

  CASE(21,31) !!! Angular coefficients  A1
      LeptonCutFactorKPC=CutFactor(qT_in=kin(1),Q_in=kin(3),y_in=kin(6),CutParameters=CutParam,Cut_Type=1)

  CASE(22,32) !!! Angular coefficients A2
      LeptonCutFactorKPC=CutFactor(qT_in=kin(1),Q_in=kin(3),y_in=kin(6),CutParameters=CutParam,Cut_Type=2)

  CASE(23) !!! Angular coefficients A3
      LeptonCutFactorKPC=CutFactor(qT_in=kin(1),Q_in=kin(3),y_in=kin(6),CutParameters=CutParam,Cut_Type=3)

  CASE(24) !!! Angular coefficients A4
      LeptonCutFactorKPC=CutFactor(qT_in=kin(1),Q_in=kin(3),y_in=kin(6),CutParameters=CutParam,Cut_Type=4)

  CASE(35) !!! Angular coefficients A5
      LeptonCutFactorKPC=CutFactor(qT_in=kin(1),Q_in=kin(3),y_in=kin(6),CutParameters=CutParam,Cut_Type=5)

  CASE(36) !!! Angular coefficients A6
      LeptonCutFactorKPC=CutFactor(qT_in=kin(1),Q_in=kin(3),y_in=kin(6),CutParameters=CutParam,Cut_Type=6)

  CASE DEFAULT
      !!! default is spherical cut-factor
      LeptonCutFactorKPC=CutFactor(qT_in=kin(1),Q_in=kin(3),y_in=kin(6),CutParameters=CutParam,Cut_Type=-1)

  END SELECT


  else
    !!! this is uncut lepton tensor
    LeptonCutFactorKPC=1._dp
  end if

end function LeptonCutFactorKPC

!!! check is the process y-symmetric
pure function IsySymmetric(p2)
  logical::IsySymmetric
  integer,intent(in)::p2
  if(p2==1 .or. p2==3 ) then
    IsySymmetric=.true.
  else
    IsySymmetric=.false.
  end if
end function IsySymmetric


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS CALCULATING CROSS-SECTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! function determines the best value of PT-sections from PT-bin size, and Q
!!! it is determined by formula Q/PT< val/ (2 k) => def+2K
function NumPT_auto(dPT,Q)
  real,parameter::val=40.
  real(dp)::dPT,Q,rat
  integer::i,NumPT_auto
  rat=Q/dPT

  if(rat>40.) then
      NumPT_auto=NumPTdefault
      return
  else
      do i=1,25
          if(rat>(40d0/2d0/i)) then
              NumPT_auto=NumPTdefault+2*i
              return
          end if
      end do
  end if
  if(outputlevel>1) then
    write(*,*) WarningString('Fail to automatically determine number of Pt-section for a bin.',moduleName)
    write(*,*) '  ... Possibly Pt-bin is too large', dPT
  end if
  NumPT_auto=NumPTdefault+12
end function NumPT_auto

!---------------------------------UNINTEGRATED------------------------------------------------------------------

!!! this is help function which evaluate xSec at single qt (without lists) with only prefactor 2
!!!! this is extended (and default) version of xSec, which include all parameters
function xSec(var,process,incCut,CutParam)
  real(dp):: xSec,FF,LC
  real(dp)::x1,x2,scaleMu,scaleZeta
  real(dp),dimension(1:7),intent(in)::var
  logical,intent(in)::incCut
  real(dp),dimension(1:4),intent(in)::CutParam
  integer,dimension(1:4),intent(in)::process

  GlobalCounter=GlobalCounter+1

  !!!!!!!!!!!!!!!!!!!! COMPUTATION WITH KPC
  if(useKPC) then
    if(TMDF_KPC_IsconvergenceLost()) then
      xSec=1d9
      return
    end if

    !!! setting values of X
    x1=var(5)*var(7)
    x2=var(5)/var(7)

    !!! setting values of scales
    scaleZeta=var(4)  !! zeta=Q2
    scaleMu=sqrt(scaleZeta)

    FF=KPC_DYconv(var(4),var(1),x1,x2,scaleMu*c2_global,process(2:4))
    LC=LeptonCutFactorKPC(var,process(4),incCut,CutParam)
    xSec=PreFactorKPC(var,process(1))*FF*LC

  !!!!!!!!!!!!!!!!!!!! COMPUTATION WITHOUT KPC
  else
    if(TMDF_IsconvergenceLost()) then
      xSec=1d9
      return
    end if

    !!! setting values of X
    x1=var(5)*var(7)
    x2=var(5)/var(7)

    !!! setting values of scales
    scaleZeta=var(4)+exactScales*var(1)**2  !! zeta=Q2+qT^2
    scaleMu=sqrt(scaleZeta)

    !!!!! compute cross-section for each process
    FF=TMDF_F(var(4),var(1),x1,x2,scaleMu*c2_global,scaleZeta,scaleZeta,process(2:4))
    LC=LeptonCutFactorLP(var,process(4),incCut,CutParam)

    xSec=PreFactor2(var,process(1))*FF*LC
  end if

  !write(*,*) "{",var(4),",",x1,"},"!,z1
end function xSec
  
  !---------------------------------INTEGRATED over Y---------------------------------------------------------------
  

!!!function for integration over Y
function Xsec_Yint(var,process,incCut,CutParam,ymin_in,ymax_in)
  real(dp),dimension(1:7) :: var
  logical,intent(in)::incCut
  real(dp),dimension(1:4),intent(in)::CutParam
  integer,dimension(1:4)::process
  real(dp) :: Xsec_Yint

  real(dp) :: ymin, ymax,ymin_in,ymax_in
  real(dp) :: ymin_Check,ymax_Check


  if(TMDF_IsconvergenceLost()) then
    Xsec_Yint=1d9
    return
  end if

  !!! evaluate correspnding y's
  !!! in the case process=2 the integral is over xF
  if(process(1)==2) then
    ymin=yFromXF(ymin_in,var)
    ymax=yFromXF(ymax_in,var)

    process(1)=1 !!!! this is important because the actual integration is over y, and process=2 contains Jacobian.
  else
    ymin=ymin_in
    ymax=ymax_in
  end if

  ymin_Check=log(var(5))+0.000000001d0
  ymax_Check=-log(var(5))-0.000000001d0

  if(IsySymmetric(process(1)) .and. (ABS(ymax+ymin)<toleranceGEN)) then!!! symetric integral
    if(ymax > ymax_check) then
        ymax=ymax_Check
    end if!!!!! else case: automatically taken into account

#if INTEGRATION_MODE==1
    !!!! slower but accurate
    Xsec_Yint=2*Integrate_SA(integrandOverY,0._dp,ymax,toleranceINT)
#elif INTEGRATION_MODE==2
   !!!! fast but not that accurate
    Xsec_Yint=2*Integrate_G7(integrandOverY,0._dp,ymax)
#endif

  else !!!non-symmetric integral!!!!!!!!
    if(ymax<ymin_check .or. ymin>ymax_check) then !!! the case then y is outside physicsl region
      Xsec_Yint=0._dp
    else
    if(ymax > ymax_check) then
      ymax=yMax_check
    end if!!!!! else case: automatically taken into account
    if(ymin < ymin_check) then
      ymin=ymin_check
    end if!!!!! else case: automatically taken into account

#if INTEGRATION_MODE==1
    !!!! slower but accurate
    Xsec_Yint=Integrate_SA(integrandOverY,ymin,ymax,toleranceINT)
#elif INTEGRATION_MODE==2
    !!!! fast but not that accurate
    Xsec_Yint=Integrate_G7(integrandOverY,ymin,ymax)
#endif

    end if
end if

contains

function integrandOverY(y)
real(dp),intent(in)::y
real(dp)::integrandOverY
call SetY(y,var)
integrandOverY=xSec(var,process,incCut,CutParam)
end function integrandOverY

end function Xsec_Yint
  
!---------------------------------INTEGRATED over Y over Q---------------------------------------------------------------
!!!! No need for check over Y they take a place within y-integration for each value of Q(!)

!!!! to integrate over Q I use adaptive Simpson. (defined below)
!!!! before the integration over Q I check the size of Q-bin,
!!!! if Q-bin is large I split desect the integration range to smaller
function Xsec_Qint_Yint(var,process,incCut,CutParam,Qmin_in,Qmax_in,ymin_in,ymax_in)
  real(dp),dimension(1:7)::var
  logical,intent(in)::incCut
  real(dp),dimension(1:4),intent(in)::CutParam
  integer,dimension(1:4),intent(in)::process
  real(dp),intent(in) :: yMin_in,yMax_in,QMin_in,QMax_in
  real(dp):: Xsec_Qint_Yint
  integer::numSec,i
  real(dp)::dQ

  if(TMDF_IsconvergenceLost()) then
    Xsec_Qint_Yint=1d9
    return
  end if

!!!! slower but accurate
#if INTEGRATION_MODE==1
  !!! check how many maxQbins is inside the integration range (+1)
  numSec=INT((Qmax_in-Qmin_in)/maxQbinSize)+1

  !!! if the bin is smaller than maxQbinSize, integrate as is
  if(numSec==1) then

    Xsec_Qint_Yint=Integrate_SA(integrandOverQ,Qmin_in,Qmax_in,toleranceINT)
  else
    !!! else divide to smaler bins and sum the integrals
    dQ=(Qmax_in-Qmin_in)/numSec !!! size of new bins

    Xsec_Qint_Yint=0._dp
    do i=0,numSec-1
      Xsec_Qint_Yint=Xsec_Qint_Yint + &
          Integrate_SA(integrandOverQ,Qmin_in+i*dQ,Qmin_in+(i+1)*dQ,toleranceINT)
    end do
  end if
!!!! slower but accurate
#elif INTEGRATION_MODE==2
!!!! in this case I only check for the Z-boson peak
!!!! if it is in the region, I integrate over it specially
  if(Qmin_in<MZ-3 .and. MZ+3<QMax_in) then
    Xsec_Qint_Yint=Integrate_G7(integrandOverQ,Qmin_in,MZ-3)
    Xsec_Qint_Yint=Xsec_Qint_Yint+Integrate_SA(integrandOverQ,MZ-3,MZ+3,toleranceINT)
    Xsec_Qint_Yint=Xsec_Qint_Yint+Integrate_G7(integrandOverQ,MZ+3,Qmax_in)
  else
    Xsec_Qint_Yint=Integrate_G7(integrandOverQ,Qmin_in,Qmax_in)
  end if
#endif

contains

function integrandOverQ(Q)
real(dp),intent(in)::Q
real(dp)::integrandOverQ
call SetQ(Q,var)
integrandOverQ=2*Q*Xsec_Yint(var,process,incCut,CutParam,yMin_in,yMax_in)
end function integrandOverQ

end function Xsec_Qint_Yint

!---------------------------------INTEGRATED over Y over Q over pT-------------------------------------------------------------
!!!!! In fact, this is the main evaluator.
!!!integration over PT is made by Num-sections
!!!N even
function Xsec_PTint_Qint_Yint(process,incCut,CutParam,s_in,qt_min_in,qt_max_in,Q_min_in,Q_max_in,ymin_in,ymax_in,Num)
real(dp),dimension(1:7)::var
logical,intent(in)::incCut
real(dp),dimension(1:4),intent(in)::CutParam
integer,dimension(1:4),intent(in)::process
real(dp):: Xsec_PTint_Qint_Yint
real(dp),intent(in):: ymin_in,ymax_in,Q_min_in,Q_max_in,qt_min_in,qt_max_in,s_in
real(dp):: Q_min,Q_max,qt_min,qt_max,s
integer,intent(in) :: Num

if(TMDF_IsconvergenceLost()) then
  Xsec_PTint_Qint_Yint=1d9
  return
end if

!!!------------------------- checking Q----------
if(Q_min_in<0.9d0) then
  call Warning_Raise('Attempt to compute xSec with Q<0.9.',messageCounter,messageTrigger,moduleName)
  write(*,*) "Qmin =",Q_min_in," (Qmin set to 1.GeV)"
  Q_min=1._dp
else
  Q_min=Q_min_in
end if

if(Q_max_in<Q_min) then
  call Warning_Raise('Attempt to compute xSec with Qmax<Qmin. RESULT 0',messageCounter,messageTrigger,moduleName)
  Xsec_PTint_Qint_Yint=0._dp
  return
end if
Q_max=Q_max_in

!!!------------------------- checking S----------
if(s_in<0.9d0) then
  call Warning_Raise('Attempt to compute xSec with s<0.9.',messageCounter,messageTrigger,moduleName)
  write(*,*) "s =",s_in," (s set to Qmin)"
  s=Q_min
else
  s=s_in
end if

!!!------------------------- checking PT----------
if(qT_min_in<0.0d0) then
  call Warning_Raise('Attempt to compute xSec with qT<0.',messageCounter,messageTrigger,moduleName)
  write(*,*) "qTmin =",qT_min_in," (qTmin set to 0.GeV)"
  qT_min=1._dp
else
  qT_min=qT_min_in
end if

if(qT_max_in<qT_min) then
  call Warning_Raise('Attempt to compute xSec with qTmax<qTmin. RESULT 0',messageCounter,messageTrigger,moduleName)
  Xsec_PTint_Qint_Yint=0._dp
  return
end if
qT_max=qT_max_in

!!!------------------------- checking Y----------

if(ymax_in<ymin_in) then
  call Warning_Raise('Attempt to compute xSec with Ymax<Ymin. RESULT 0',messageCounter,messageTrigger,moduleName)
  Xsec_PTint_Qint_Yint=0._dp
  return
end if


if(qt_min<qTMin_ABS) then
  var=kinematicArray(qTMin_ABS,s_in,(Q_min+Q_max)/2,(ymin_in+ymax_in)/2)
else
  var=kinematicArray(qt_min,s_in,(Q_min+Q_max)/2,(ymin_in+ymax_in)/2)
end if

Xsec_PTint_Qint_Yint=Integrate_SN(integrandOverQT,qt_min,qt_max,Num)

contains

function integrandOverQT(qT)
real(dp),intent(in)::qT
real(dp)::integrandOverQT
call SetQT(qT,var)
integrandOverQT=2*qT*Xsec_Qint_Yint(var,process,incCut,CutParam,Q_min,Q_max,ymin_in,ymax_in)
end function integrandOverQT

end function Xsec_PTint_Qint_Yint

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!APROXIMATE VERSION OF THE SAME ROUTINES!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !---------------------------------INTEGRATED over Y---------------------------------------------------------------


!!!function for integration over Y
function Xsec_Yint_APROX(var,process,incCut,CutParam,ymin_in,ymax_in)
  real(dp),dimension(1:7) :: var
  logical,intent(in)::incCut
  real(dp),dimension(1:4),intent(in)::CutParam
  real(dp) :: Xsec_Yint_APROX
  real(dp) :: ymin, ymax,ymin_in,ymax_in
  real(dp) :: ymin_Check,ymax_Check
  integer,dimension(1:4),intent(in)::process

  if(TMDF_IsconvergenceLost()) then
    Xsec_Yint_APROX=1d9
    return
  end if

  !!! evaluate correspnding y's
  !!! in the case process=2 the integral is over xF
  if(process(1)==2) then
    ymin=yFromXF(ymin_in,var)
    ymax=yFromXF(ymax_in,var)
  else
    ymin=ymin_in
    ymax=ymax_in
  end if

  ymin_Check=log(var(5))+0.000000001d0
  ymax_Check=-log(var(5))-0.000000001d0

  if(IsySymmetric(process(1)) .and. (ABS(ymax+ymin)<toleranceGEN)) then!!! symetric integral
    if(ymax > ymax_check) then
        ymax=ymax_Check
    end if!!!!! else case: automatically taken into account

    Xsec_Yint_APROX=2*Integrate_G3(integrandOverY,0._dp,ymax)

  else !!!non-symmetric integral!!!!!!!!
    if(ymax<ymin_check .or. ymin>ymax_check) then !!! the case then y is outside physicsl region
      Xsec_Yint_APROX=0._dp
    else
    if(ymax > ymax_check) then
      ymax=yMax_check
    end if!!!!! else case: automatically taken into account
    if(ymin < ymin_check) then
      ymin=ymin_check
    end if!!!!! else case: automatically taken into account

    Xsec_Yint_APROX=Integrate_G3(integrandOverY,ymin,ymax)

    end if
end if

contains

function integrandOverY(y)
real(dp),intent(in)::y
real(dp)::integrandOverY
call SetY(y,var)
integrandOverY=xSec(var,process,incCut,CutParam)
end function integrandOverY

end function Xsec_Yint_APROX

!---------------------------------INTEGRATED over Y over Q---------------------------------------------------------------
!!!! No need for check over Y they take a place within y-integration for each value of Q(!)

!!!! to integrate over Q I use adaptive Simpson. (defined below)
!!!! before the integration over Q I check the size of Q-bin,
!!!! if Q-bin is large I split desect the integration range to smaller
function Xsec_Qint_Yint_APROX(var,process,incCut,CutParam,Qmin_in,Qmax_in,ymin_in,ymax_in)
  real(dp),dimension(1:7)::var
  logical,intent(in)::incCut
  real(dp),dimension(1:4),intent(in)::CutParam
  integer,dimension(1:4),intent(in)::process
  real(dp),intent(in) :: yMin_in,yMax_in,QMin_in,QMax_in
  real(dp):: Xsec_Qint_Yint_APROX
  integer::numSec,i
  real(dp)::dQ

  if(TMDF_IsconvergenceLost()) then
    Xsec_Qint_Yint_APROX=1d9
    return
  end if

!!!! in this case I only check for the Z-boson peak
!!!! if it is in the region, I integrate over it specially
  if(Qmin_in<MZ-3 .and. MZ+3<QMax_in) then
    Xsec_Qint_Yint_APROX=Integrate_G3(integrandOverQ,Qmin_in,MZ-3)
    Xsec_Qint_Yint_APROX=Xsec_Qint_Yint_APROX+Integrate_G7(integrandOverQ,MZ-3,MZ+3)
    Xsec_Qint_Yint_APROX=Xsec_Qint_Yint_APROX+Integrate_G3(integrandOverQ,MZ+3,Qmax_in)
  else
    Xsec_Qint_Yint_APROX=Integrate_G3(integrandOverQ,Qmin_in,Qmax_in)
  end if

contains

function integrandOverQ(Q)
real(dp),intent(in)::Q
real(dp)::integrandOverQ
call SetQ(Q,var)
integrandOverQ=2*Q*Xsec_Yint_APROX(var,process,incCut,CutParam,yMin_in,yMax_in)
end function integrandOverQ

end function Xsec_Qint_Yint_APROX

!---------------------------------INTEGRATED over Y over Q over pT-------------------------------------------------------------
!!!!! In fact, this is the main evaluator.
!!!integration over PT is made by Num-sections
!!!N even
function Xsec_PTint_Qint_Yint_APROX(process,incCut,CutParam,s_in,qt_min_in,qt_max_in,Q_min_in,Q_max_in,ymin_in,ymax_in)
  real(dp),dimension(1:7)::var
  logical,intent(in)::incCut
  real(dp),dimension(1:4),intent(in)::CutParam
  integer,dimension(1:4),intent(in)::process
  real(dp):: Xsec_PTint_Qint_Yint_APROX
  real(dp),intent(in):: ymin_in,ymax_in,Q_min_in,Q_max_in,qt_min_in,qt_max_in,s_in
  real(dp):: Q_min,Q_max,qt_min,qt_max,s

  if(TMDF_IsconvergenceLost()) then
    Xsec_PTint_Qint_Yint_APROX=1d9
    return
  end if

  !!!------------------------- checking Q----------
  if(Q_min_in<0.9d0) then
    call Warning_Raise('Attempt to compute xSec with Q<0.9.',messageCounter,messageTrigger,moduleName)
    write(*,*) "Qmin =",Q_min_in," (Qmin set to 1.GeV)"
    Q_min=1._dp
  else
    Q_min=Q_min_in
  end if

  if(Q_max_in<Q_min) then
    call Warning_Raise('Attempt to compute xSec with Qmax<Qmin. RESULT 0',messageCounter,messageTrigger,moduleName)
    Xsec_PTint_Qint_Yint_APROX=0._dp
    return
  end if
  Q_max=Q_max_in

  !!!------------------------- checking S----------
  if(s_in<0.9d0) then
    call Warning_Raise('Attempt to compute xSec with s<0.9.',messageCounter,messageTrigger,moduleName)
    write(*,*) "s =",s_in," (s set to Qmin)"
    s=Q_min
  else
    s=s_in
  end if

  !!!------------------------- checking PT----------
  if(qT_min_in<0.0d0) then
    call Warning_Raise('Attempt to compute xSec with qT<0.',messageCounter,messageTrigger,moduleName)
    write(*,*) "qTmin =",qT_min_in," (qTmin set to 0.GeV)"
    qT_min=1._dp
  else
    qT_min=qT_min_in
  end if

  if(qT_max_in<qT_min) then
    call Warning_Raise('Attempt to compute xSec with qTmax<qTmin. RESULT 0',messageCounter,messageTrigger,moduleName)
    Xsec_PTint_Qint_Yint_APROX=0._dp
    return
  end if
  qT_max=qT_max_in

  !!!------------------------- checking Y----------

  if(ymax_in<ymin_in) then
    call Warning_Raise('Attempt to compute xSec with Ymax<Ymin. RESULT 0',messageCounter,messageTrigger,moduleName)
    Xsec_PTint_Qint_Yint_APROX=0._dp
    return
  end if


  if(qt_min<qTMin_ABS) then
    var=kinematicArray(qTMin_ABS,s_in,(Q_min+Q_max)/2,(ymin_in+ymax_in)/2)
  else
    var=kinematicArray(qt_min,s_in,(Q_min+Q_max)/2,(ymin_in+ymax_in)/2)
  end if

  if(qt_max-qt_min<0.25d0) then
    !!!! for small bins just compute at the center
    call SetQT((qT_min+qT_max)/2,var)
    Xsec_PTint_Qint_Yint_APROX=(qT_max**2-qT_min**2)&
        *Xsec_Qint_Yint_APROX(var,process,incCut,CutParam,Q_min,Q_max,ymin_in,ymax_in)
  else
    Xsec_PTint_Qint_Yint_APROX=Integrate_G3(integrandOverQT,qt_min,qt_max)
  end if

contains

function integrandOverQT(qT)
real(dp),intent(in)::qT
real(dp)::integrandOverQT
call SetQT(qT,var)
integrandOverQT=2*qT*Xsec_Qint_Yint_APROX(var,process,incCut,CutParam,Q_min,Q_max,ymin_in,ymax_in)
end function integrandOverQT

end function Xsec_PTint_Qint_Yint_APROX

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!THE MAIN INTERFACE TO CROSS-SECTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! interface for array,s,array,array,array,logical,optional, optional
subroutine MainInterface_AsAAAloo(X,process,s,qT,Q,y,includeCuts,CutParameters,Num)
!   function xSec_DY(process,s,qT,Q,y,includeCuts,CutParameters,Num)
  integer,intent(in),dimension(1:4)::process            !the number of process
  real(dp),intent(in)::s                          !Mandelshtam s
  real(dp),intent(in),dimension(1:2)::qT           !(qtMin,qtMax)
  real(dp),intent(in),dimension(1:2)::Q            !(Qmin,Qmax)
  real(dp),intent(in),dimension(1:2)::y             !(ymin,ymax)
  logical,intent(in)::includeCuts                !include cuts
  real(dp),intent(in),dimension(1:4),optional::CutParameters    !(p1,p2,eta1,eta2)
  integer,intent(in),optional::Num                !number of sections

  real(dp)::X

  integer::nn
  real(dp),dimension(1:4)::CutParam

if(.not.started) ERROR STOP ErrorString('The module is not initialized. Check INI-file.',moduleName)

  !! determine number of sections
  if(present(Num)) then
    nn=Num
  else
    nn=NumPT_auto(qT(2)-qT(1),(Q(2)+Q(1))/2d0)
  end if

  !!! determine cut parameters
  if(includeCuts) then
    if(present(CutParameters)) then
      CutParam=CutParameters
    else
      write(*,*) ErrorString('called includeCuts=true, while CutParameters are undefined',moduleName)
      write(*,*) ErrorString('Evaluation stop',moduleName)
      stop
    end if
  else
    CutParam=(/0d0,0d0,0d0,0d0/)
  end if

  !!!! evaluation
  CallCounter=CallCounter+1
  X=Xsec_PTint_Qint_Yint(process,includeCuts,CutParameters,&
                  s,qT(1),qT(2),Q(1),Q(2),y(1),y(2),nn)

end subroutine MainInterface_AsAAAloo
  
subroutine xSec_DY_List(X,process,s,qT,Q,y,includeCuts,CutParameters,Num)
  integer,intent(in),dimension(:,:)::process            !the number of process
  real(dp),intent(in),dimension(:)::s                !Mandelshtam s
  real(dp),intent(in),dimension(:,:)::qT            !(qtMin,qtMax)
  real(dp),intent(in),dimension(:,:)::Q                !(Qmin,Qmax)
  real(dp),intent(in),dimension(:,:)::y                !(ymin,ymax)
  logical,intent(in),dimension(:)::includeCuts        !include cuts
  real(dp),intent(in),dimension(:,:)::CutParameters            !(p1,p2,eta1,eta2)
  integer,intent(in),dimension(:),optional::Num        !number of sections
  real(dp),dimension(:),intent(out)::X
  integer :: i,length
  integer,allocatable::nn(:)

if(.not.started) ERROR STOP ErrorString('The module is not initialized. Check INI-file.',moduleName)

  length=size(s)

!!! cheking sizes
if(size(X)/=length) then
  write(*,*) ErrorString('xSec_DY_List: sizes of xSec and s lists are not equal.',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(process,1)/=length) then
  write(*,*) ErrorString('xSec_DY_List: sizes of process and s lists are not equal.',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(qT,1)/=length) then
  write(*,*) ErrorString('xSec_DY_List: sizes of qT and s lists are not equal.',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(y,1)/=length) then
  write(*,*) ErrorString('xSec_DY_List: sizes of y and s lists are not equal.',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(Q,1)/=length) then
  write(*,*) ErrorString('xSec_DY_List: sizes of Q and s lists are not equal.',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(includeCuts)/=length) then
  write(*,*) ErrorString('xSec_DY_List: sizes of includeCuts and s lists are not equal.',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(CutParameters,1)/=length) then
  write(*,*) ErrorString('xSec_DY_List: sizes of CutParameters and s lists are not equal.',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(process,2)/=4) then
  write(*,*) ErrorString('xSec_DY_List: process list must be (:,1:4).',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(qT,2)/=2) then
  write(*,*) ErrorString('xSec_DY_List: qt list must be (:,1:2).',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(y,2)/=2) then
  write(*,*) ErrorString('xSec_DY_List: y list must be (:,1:2).',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(Q,2)/=2) then
  write(*,*) ErrorString('xSec_DY_List: Q list must be (:,1:2).',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if


  CallCounter=CallCounter+length

  allocate(nn(1:length))
  if(present(Num)) then
      if(size(Num,1)/=length) then
    write(*,*) 'ERROR: arTeMiDe_DY: xSec_DY_List: sizes of Num and s lists are not equal.'
    write(*,*) 'Evaluation stop'
    stop
  end if
  nn=Num
  else
      do i=1,length
          nn=NumPT_auto(qT(i,2)-qT(i,1),(Q(i,2)+Q(i,1))/2d0)
      end do
  end if
  !$OMP PARALLEL DO DEFAULT(SHARED)

    do i=1,length
      X(i)=Xsec_PTint_Qint_Yint(process(i,1:4),includeCuts(i),CutParameters(i,1:4),&
              s(i),qT(i,1),qT(i,2),Q(i,1),Q(i,2),y(i,1),y(i,2),nn(i))
    end do
  !$OMP END PARALLEL DO
  deallocate(nn)
end subroutine xSec_DY_List
  
subroutine xSec_DY_List_BINLESS(X,process,s,qT,Q,y,includeCuts,CutParameters)
  integer,intent(in),dimension(:,:)::process            !the number of process
  real(dp),intent(in),dimension(:)::s                !Mandelshtam s
  real(dp),intent(in),dimension(:)::qT            !(qtMin,qtMax)
  real(dp),intent(in),dimension(:)::Q                !(Qmin,Qmax)
  real(dp),intent(in),dimension(:)::y                !(ymin,ymax)
  logical,intent(in),dimension(:)::includeCuts        !include cuts
  real(dp),intent(in),dimension(:,:)::CutParameters            !(p1,p2,eta1,eta2)
  real(dp),dimension(:),intent(out)::X

  real(dp),allocatable,dimension(:,:)::vv
  integer :: i,length

if(.not.started) ERROR STOP ErrorString('The module is not initialized. Check INI-file.',moduleName)

  length=size(s)

  !!! cheking sizes
if(size(X)/=length) then
  write(*,*) ErrorString('xSec_DY_List_BINLESS: sizes of xSec and s lists are not equal.',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(process,1)/=length) then
  write(*,*) ErrorString('xSec_DY_List_BINLESS: sizes of process and s lists are not equal.',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(qT)/=length) then
  write(*,*) ErrorString('xSec_DY_List_BINLESS: sizes of qT and s lists are not equal.',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(y)/=length) then
  write(*,*) ErrorString('xSec_DY_List_BINLESS: sizes of y and s lists are not equal.',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(Q)/=length) then
  write(*,*) ErrorString('xSec_DY_List_BINLESS: sizes of Q and s lists are not equal.',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(includeCuts)/=length) then
  write(*,*) ErrorString('xSec_DY_List_BINLESS: sizes of includeCuts and s lists are not equal.',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(CutParameters,1)/=length) then
  write(*,*) ErrorString('xSec_DY_List_BINLESS: sizes of CutParameters and s lists are not equal.',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(process,2)/=4) then
  write(*,*) ErrorString('xSec_DY_List_BINLESS: process list must be (:,1:4).',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if

  CallCounter=CallCounter+length

  allocate(vv(1:length,1:7))

  !$OMP PARALLEL DO DEFAULT(SHARED)

    do i=1,length
      vv(i,1:7)=kinematicArray(qt(i),s(i),Q(i),y(i))
      X(i)=xSec(vv(i,1:7),process(i,1:4),includeCuts(i),CutParameters(i,1:4))

    end do
  !$OMP END PARALLEL DO
  deallocate(vv)
end subroutine xSec_DY_List_BINLESS


!!!! This subroutine simplifies all integrations to few points
!!!! It is very inaccurate, but very fast
subroutine xSec_DY_List_APPROXIMATE(X,process,s,qT,Q,y,includeCuts,CutParameters)
  integer,intent(in),dimension(:,:)::process            !the number of process
  real(dp),intent(in),dimension(:)::s                !Mandelshtam s
  real(dp),intent(in),dimension(:,:)::qT            !(qtMin,qtMax)
  real(dp),intent(in),dimension(:,:)::Q                !(Qmin,Qmax)
  real(dp),intent(in),dimension(:,:)::y                !(ymin,ymax)
  logical,intent(in),dimension(:)::includeCuts        !include cuts
  real(dp),intent(in),dimension(:,:)::CutParameters            !(p1,p2,eta1,eta2)
  real(dp),dimension(:),intent(out)::X
  integer :: i,length

if(.not.started) ERROR STOP ErrorString('The module is not initialized. Check INI-file.',moduleName)

  length=size(s)

!!! cheking sizes
if(size(X)/=length) then
  write(*,*) ErrorString('xSec_DY_List: sizes of xSec and s lists are not equal.',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(process,1)/=length) then
  write(*,*) ErrorString('xSec_DY_List: sizes of process and s lists are not equal.',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(qT,1)/=length) then
  write(*,*) ErrorString('xSec_DY_List: sizes of qT and s lists are not equal.',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(y,1)/=length) then
  write(*,*) ErrorString('xSec_DY_List: sizes of y and s lists are not equal.',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(Q,1)/=length) then
  write(*,*) ErrorString('xSec_DY_List: sizes of Q and s lists are not equal.',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(includeCuts)/=length) then
  write(*,*) ErrorString('xSec_DY_List: sizes of includeCuts and s lists are not equal.',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(CutParameters,1)/=length) then
  write(*,*) ErrorString('xSec_DY_List: sizes of CutParameters and s lists are not equal.',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(process,2)/=4) then
  write(*,*) ErrorString('xSec_DY_List: process list must be (:,1:4).',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(qT,2)/=2) then
  write(*,*) ErrorString('xSec_DY_List: qt list must be (:,1:2).',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(y,2)/=2) then
  write(*,*) ErrorString('xSec_DY_List: y list must be (:,1:2).',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if
if(size(Q,2)/=2) then
  write(*,*) ErrorString('xSec_DY_List: Q list must be (:,1:2).',moduleName)
  write(*,*) ErrorString('Evaluation stop',moduleName)
  stop
end if

  CallCounter=CallCounter+length

  !$OMP PARALLEL DO DEFAULT(SHARED)

    do i=1,length
      X(i)=Xsec_PTint_Qint_Yint_APROX(process(i,1:4),includeCuts(i),CutParameters(i,1:4),&
              s(i),qT(i,1),qT(i,2),Q(i,1),Q(i,2),y(i,1),y(i,2))
    end do
  !$OMP END PARALLEL DO
end subroutine xSec_DY_List_APPROXIMATE
  
end module TMDX_DY
