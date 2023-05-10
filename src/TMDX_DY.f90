!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 2.01
!
!	Evaluation of the TMD cross-section for DY-like cross-sections
!	
!	if you use this module please, quote 1706.01473
!
!	ver 1.0: release (AV, 10.05.2017)
!	ver 1.1: multiple updates (AV, 5.10.2017)
!	ver 1.2: module is renamed, and multiple renaming of functions (AV, 15.10.2017)
!	ver 1.31: part of functions migrated to TMDF, rest updated (AV, 1.06.2018)
!	ver 1.4: encapsulation of cuts, and process,+ multiple updates (AV, 18.01.2019)
!	ver 2.01: Added Higgs xSec, piresum option, and coefficient function moved to separate file (AV, 17.06.2019)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDX_DY
use aTMDe_Numerics
use IO_functions
use TMDF
use LeptonCutsDY
use QCDinput
use EWinput

implicit none
  private
  
   !Current version of module
 character (len=7),parameter :: moduleName="TMDX-DY"
 character (len=5),parameter :: version="v2.06"
 !Last appropriate verion of constants-file
  integer,parameter::inputver=21
  
  real(dp) :: tolerance=0.0005d0
  
  integer::outputlevel
  integer::messageTrigger
  
  !!!!
  !!!! in the module the kinematic is stored in the varibles "kinematic" real(dp),dimension(1:6)
  !!!! which is (qT,s,Q,Q^2,x0,y,exp[y])
  !!!! where x0=sqrt[(Q^2+q_T^2)/s]   (if exactX1X2) or x0=Q^2/s (otherwise)
  !!!!
  
  !The variables for all parameters of the model!
  real(dp):: s_global,Q_global,y_global
  !! Set of process definition, for Prefactor 1, Prefactor 2, structure function, etc
  !! = (/p1,p2,p3/)
  integer,dimension(1:3)::process_global
  !!cut parameters
  !!!!! this variable = (pT1,pT2,etaMIN,etaMAX)
  real(dp),dimension(1:4)::CutParameters_global
  
  !!other global parameters see SetXParameters  
  integer:: orderH_global
  logical:: includeCuts_global
  logical::usePIresum
  integer:: exactX1X2    !!!=1 if exact x's=true, =0 otherwise
  integer:: exactScales  !!!=1 if exact hard scales = true, =0 otherwise
  
  !!! number of sections for PT-integral by default
  integer::NumPTdefault=4
  !!! Maximum size of Q-bin. Larger bins are desected
  real::maxQbinSize=30.
  
  real(dp)::c2_global!,muHard_global
  

  
  integer::GlobalCounter
  integer::CallCounter
  integer::messageCounter
  
  real(dp)::hc2
  
  logical::started=.false.
  
  public::TMDX_DY_XSetup,TMDX_DY_Initialize,TMDX_DY_SetCuts,TMDX_DY_SetScaleVariation,&
    TMDX_DY_setProcess,TMDX_DY_ShowStatistic,TMDX_DY_ResetCounters,TMDX_DY_IsInitialized
  public::  CalcXsec_DY,CalcXsec_DY_Yint,CalcXsec_DY_Qint_Yint,CalcXsec_DY_PTint_Qint_Yint,CalcXsec_DY_Qint,xSec_DY,xSec_DY_List,&
    xSec_DY_List_BINLESS
  
 interface TMDX_DY_SetCuts
    module procedure SetCuts_sym,SetCuts_asym
 end interface
  
 interface CalcXsec_DY
    module procedure xSecSingle,xSecList
  end interface
  
 interface TMDX_DY_setProcess
    module procedure TMDX_setProcess1,TMDX_setProcess3,TMDX_setProcess30
 end interface
 
 interface CalcXsec_DY_Yint
    module procedure xSecSingle_Yint ,xSecList_Yint,xSecSingle_Ycomplete ,xSecList_Ycomplete
 end interface
 
  interface CalcXsec_DY_Qint
    module procedure xSecSingle_Qint ,xSecList_Qint
 end interface
 
 interface CalcXsec_DY_Qint_Yint
    module procedure xSecSingle_Qint_Yint, xSecList_Qint_Yint, xSecSingle_Qint_Ycomplete, xSecList_Qint_Ycomplete
 end interface
 

 interface CalcXsec_DY_PTint_Qint_Yint
    module procedure xSecSingle_PTint_Qint_Yint, xSecList_PTint_Qint_Yint, &
	    xSecSingle_PTint_Qint_Ycomplete, xSecList_PTint_Qint_Ycomplete,&
	    xSecSingle_PTintN_Qint_Yint, xSecList_PTintN_Qint_Yint,&
	    xSecSingle_PTintN_Qint_Ycomplete, xSecList_PTintN_Qint_Ycomplete,&
	    xSecListList_PTint_Qint_Yint,xSecListList_PTint_Qint_Ycomplete,&
	    xSecListList_PTintN_Qint_Yint,xSecListList_PTintN_Qint_Ycomplete,&
	    xSecListPY_PTint_Qint_Yint,xSecListPY_PTintN_Qint_Yint
    
 end interface
 
 interface xSec_DY
    module procedure MainInterface_AsAAAloo,MainInterface_isAAAloo
 end interface
 
contains
  
  function TMDX_DY_IsInitialized()
  logical::TMDX_DY_IsInitialized
  TMDX_DY_IsInitialized=started
 end function TMDX_DY_IsInitialized

   !! Initialization of the package
  subroutine TMDX_DY_Initialize(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=300)::path,line
    logical::initRequired,dummyLogical,TMDPDFgrid
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
    
    !!! go to section TMD-DY
    call MoveTO(51,'*9   ')
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
	orderH_global=2
      CASE ("NNNLO")
	orderH_global=3
      CASE ("N3LO") !!! same as NNNLO
	orderH_global=3
      CASE ("N3LO+")
	orderH_global=3
      CASE ("N4LO")
	orderH_global=4
      CASE DEFAULT
	if(outputLevel>0) write(*,*)  WarningString('try to set unknown order. Switch to NLO.',moduleName)
	orderH_global=1
     END SELECT
     if(outputLevel>1) write(*,*) '	artemide.TMDX_DY: the used order is ',trim(orderMain)
     
     !!exact values of x1x2
     call MoveTO(51,'*p2   ')
     read(51,*) dummyLogical
     if(dummyLogical) then 
      exactX1X2=1
     else
      exactX1X2=0
     end if
     if(outputLevel>2 .and. dummyLogical) write(*,*) '	artemide.TMDX_DY: qT/Q corrections for x1 and x2 variables are included.'
     !! pi2 resummation
     call MoveTO(51,'*p3   ')
     read(51,*) usePIresum
     if(outputLevel>2 .and. usePIresum) write(*,*) '	artemide.TMDX_DY: pi-resummation in coef.function included.'
     !!exact values for scales
     call MoveTO(51,'*p4   ')
     read(51,*) dummyLogical
     if(dummyLogical) then 
      exactScales=1
     else
      exactScales=0
     end if
     if(outputLevel>2 .and. dummyLogical) write(*,*) '	artemide.TMDX_DY: qT/Q correction for scales variables are included.'
     
     call MoveTO(51,'*B   ')
     call MoveTO(51,'*p1  ')
     read(51,*) tolerance
     call MoveTO(51,'*p2  ')
     read(51,*) NumPTdefault
     call MoveTO(51,'*p3  ')
     read(51,*) maxQbinSize
!$    if(outputLevel>1) write(*,*) '	artemide.TMDX_DY: parallel evaluation of cross-sections is to be used'
!$    call MoveTO(51,'*C   ')
!$    call MoveTO(51,'*p1  ')
!$    read(51,*) i
!$    call OMP_set_num_threads(i)
!$    if(outputLevel>1) write(*,*) '	artemide.TMDX_DY: number of threads for parallel evaluation is set to ', i	

!$     if(outputLevel>2) write(*,*) '------TEST OF PARALLEL PROCESSING ----------'
!$OMP PARALLEL
!$     if(outputLevel>2) write(*,*) '   artemide.TMDX_DY:thread num ',  omp_get_thread_num(), ' ready.'
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
    write(*,*)  color('----- arTeMiDe.TMD_DY '//trim(version)//': .... initialized',c_green)
  end subroutine TMDX_DY_Initialize

  subroutine TMDX_DY_ResetCounters()
  if(outputlevel>2) call TMDX_DY_ShowStatistic
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
    real(dp)::c1_in,c2_in,c3_in,c4_in
    
    if(outputLevel>1) write(*,*) 'TMDX_DY: c2 scale reset:',c2_in
    
    if(c2_in<0.1d0 .or. c2_in>10.d0) then
    if(outputLevel>0) write(*,*) WarningString('variation in c2 is enourmous. c2 is set to 2',moduleName)
     c2_global=2d0
    else
    c2_global=c2_in
    end if
    
  end subroutine TMDX_DY_SetScaleVariation

  !!sets the cuts
  !! argument includeCuts_global_in=logical, if .true. will add to calculation the evaluation of cut leptonic tensor
  !! call BEFORE SetXParameters
  subroutine SetCuts_sym(include_arg,pT_arg,eta_min,eta_max)
  logical:: include_arg
  real(dp):: pT_arg,eta_max,eta_min
  
  includeCuts_global=include_arg
  
  CutParameters_global=(/pT_arg,pT_arg,eta_min,eta_max/)
  
  end subroutine SetCuts_sym
  
  !!sets the cuts (asymetric)
  !! argument includeCuts_global_in=logical, if .true. will add to calculation the evaluation of cut leptonic tensor
  !! call BEFORE SetXParameters
  subroutine SetCuts_asym(include_arg,pT1_arg,pT2_arg,eta_min,eta_max)
  logical:: include_arg
  real(dp):: pT1_arg,pT2_arg,eta_max,eta_min
  
  includeCuts_global=include_arg
  
  CutParameters_global=(/pT1_arg,pT2_arg,eta_min,eta_max/)
  
  end subroutine SetCuts_asym
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!! PROCESS DEFINITION
  
  function processArrayFromInteger(p)
    integer,intent(in)::p
    integer,dimension(1:3)::processArrayFromInteger
    SELECT CASE(p)
      case(1)
	processArrayFromInteger=(/1,1,5/) !!! p + p -> Z + gamma^*   (e.g. ATLAS, CMS, LHCb)
      case(2)
	processArrayFromInteger=(/1,1,6/) !!! p + pbar -> Z + gamma^*   (e.g. CDF,D0)
      case(4)
	processArrayFromInteger=(/1,2,1001/)
      case(5)
	processArrayFromInteger=(/1,1,5/)
      case(7)
	processArrayFromInteger=(/1,1,6/)
      case default
	write(*,*) 'ERROR: arTeMiDe_DY: unknown process is called. p=',p
	write(*,*) 'Evaluation stop'
	stop
      end SELECT
  end function processArrayFromInteger
  
    !!!set variables for process definition
  subroutine TMDX_setProcess30(p0)
  integer,dimension(1:3)::p0
  
  process_global(1)=p0(1)
  process_global(2)=p0(2)
  process_global(3)=p0(3)
  end subroutine TMDX_setProcess30
  
  !!!set variables for process definition
  subroutine TMDX_setProcess3(p1,p2,p3)
  integer::p1,p2,p3
  
  process_global(1)=p1
  process_global(2)=p2
  process_global(3)=p3
  end subroutine TMDX_setProcess3
  
  !!!set variables for process definition
  subroutine TMDX_setProcess1(p)
  integer::p
  
  call TMDX_setProcess30(processArrayFromInteger(p))
  end subroutine TMDX_setProcess1
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS FOR OPERATION WITH KINEMATICS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  
  !sets the main parameters of cross-section (x,zeta,etc)
  !the variables process defines the type of process
  subroutine TMDX_DY_XSetup(s,Q,y)
    real(dp)::s,Q,y
    
    if(.not.started) then
    write(*,*) 'ERROR: arTeMiDe.TMDX_DY is not initialized. Evaluation terminated'
    stop
    end if
    
    s_global=s
    Q_global=Q
    y_global=y
        
!     if(includeCuts_global) then
!       call SetCutParameters(pT1_global,pT2_global,eta_min_global,eta_max_global)
!     end if
    
  end subroutine TMDX_DY_XSetup
  
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


  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS FOR PREFACTORS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !!!! for easier read coeff-functions are split into separate file
  INCLUDE 'DYcoeff-func.f90'
  
  !!!!! Prefactor 2 is (universal part) x (cuts) x H
  function PreFactor2(kin,process, includeCuts_in,CutParam)
    real(dp),dimension(1:7),intent(in)::kin
    logical,intent(in)::includeCuts_in
    real(dp)::PreFactor2,cutPrefactor,uniPart,scaleMu
    real(dp),dimension(1:4),intent(in)::CutParam
    integer,dimension(1:3),intent(in)::process
  
  !!!!! cut part
    if(includeCuts_in) then
       !!! here include cuts onf lepton tensor
       cutPrefactor=CutFactor4(qT=kin(1),Q_in=kin(3),y_in=kin(6),CutParameters=CutParam)
    else
	!!! this is uncut lepton tensor
       cutPrefactor=(1+0.5d0*(kin(1)/kin(3))**2)
    end if  
  
    
   !!!! universal part
    scaleMu=sqrt(kin(4)+exactScales*kin(1)**2)
    
  SELECT CASE(process(2))
    case(-10221191)
	uniPart=1d0
	cutPrefactor=1d0
    CASE(1)
	!4 pi aEm^2/3 /Nc/Q^2/s
	uniPart=pix4/9d0*(alphaEM(scaleMu)**2)/(kin(2)*kin(4))*&
	    HardCoefficientDY(scaleMu)*&
	    hc2*1d9!from GeV to pb
	!!! IsySymmetric=.true.  !!! state is IsySymmetric-function
    CASE(2)
	!4 pi aEm^2/3 /Nc/Q^2/s
	uniPart=pix4/9d0*(alphaEM(scaleMu)**2)/(kin(2)*kin(4))*&
	    HardCoefficientDY(scaleMu)*&
	    hc2*1d9!from GeV to pb
	!!!IsySymmetric=.false.	!!! state is IsySymmetric-function
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
	    
	cutPrefactor=1d0 !!! cut-prefactor is different in this case! 
    CASE DEFAULT 
      write(*,*) 'ERROR: arTeMiDe.TMDX_DY: unknown process p2=',process(2),' .Evaluation stop.'
      stop
  END SELECT
  
  !!! this is case of xF integration the weight is 2sqrt[(Q^2+q_T^2)/s] Cosh[y]
  if(process(1)==2) then
    uniPart=2d0*kin(5)*cosh(kin(6))*uniPart
  end if
  
  PreFactor2=uniPart*cutPrefactor
  
  
  end function PreFactor2
  
  
  !!!! Set Prefactor1
  !!!! it multiplies the cross-section as a whole, does not participate in the integration
  function PreFactor1(p1)
  real(dp)::Prefactor1
  integer::p1
  SELECT CASE(p1)
    CASE(1)
	PreFactor1=1d0
    CASE(2)
	PreFactor1=1d0
    CASE DEFAULT 
      write(*,*) 'ERROR: arTeMiDe.TMDX_DY: unknown process p1=',p1,' .Evaluation stop.'
      stop
  END SELECT
  end function PreFactor1
  
  !!! check is the process y-symmetric
  function IsySymmetric(p2)
  logical::IsySymmetric
  integer::p2
  if(p2==1 .or. p2==3 ) then 
    IsySymmetric=.true.
  else
    IsySymmetric=.false.
  end if
  end function IsySymmetric
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!FUNCTIONS CALCULATING CROSS-SECTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !---------------------------------UNINTEGRATED------------------------------------------------------------------
  
  !!! this is help function which evaluate xSec at single qt (without lists) with only prefactor 2
  !!!! this is extended (and default) version of xSec, which include all parameters
  function xSec(var,process,incCut,CutParam)
    real(dp):: xSec,FF
    real(dp)::x1,x2,scaleMu,scaleZeta
    real(dp),dimension(1:7),intent(in)::var
    logical,intent(in)::incCut
    real(dp),dimension(1:4),intent(in)::CutParam
    integer,dimension(1:3),intent(in)::process
    GlobalCounter=GlobalCounter+1
    
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
  
    FF=TMDF_F(var(4),var(1),x1,x2,scaleMu*c2_global,scaleZeta,scaleZeta,process(3))
    
    xSec=PreFactor2(var,process,incCut,CutParam)*FF    
    
    !write(*,*) "{",var(4),",",x1,"},"!,z1
  end function xSec
  
  !---------------------------------INTEGRATED over Y---------------------------------------------------------------
  
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
      
  function yFromXF(xF,var)
  real(dp),dimension(1:7)::var
  real(dp):: yFromXF,xF
!     yFromXF=asinh(Sqrt(s_global/(Q2_global+qt**2))*xF/2d0)
    yFromXF=asinh(xF/2d0/var(5))
  end function yFromXF
  
  !!!
  function Xsec_Yint(var,process,incCut,CutParam,ymin_in,ymax_in)
    real(dp),dimension(1:7) :: var
    logical,intent(in)::incCut
    real(dp),dimension(1:4),intent(in)::CutParam
    real(dp) :: Xsec_Yint
    real(dp) :: ymin, ymax,ymin_in,ymax_in
    real(dp) :: ymin_Check,ymax_Check
    integer,dimension(1:3),intent(in)::process
    
    if(TMDF_IsconvergenceLost()) then 
      Xsec_Yint=1d9
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
    
    if(IsySymmetric(process(2)) .and. (ABS(ymax+ymin)<tolerance)) then!!! symetric integral
    if(ymax > ymax_check) then
        ymax=ymax_Check
    end if!!!!! else case: automatically taken into account
    
    Xsec_Yint=2d0*integralOverYpoint_S(var,process,incCut,CutParam,0d0,ymax)!!! 2 since symmetric
      
    else !!!non-symmetric integral!!!!!!!!
      if(ymax<ymin_check .or. ymin>ymax_check) then !!! the case then y is outside physicsl region
	Xsec_Yint=0d0
      else
	if(ymax > ymax_check) then
	  ymax=yMax_check
	end if!!!!! else case: automatically taken into account
	if(ymin < ymin_check) then
	  ymin=ymin_check
	end if!!!!! else case: automatically taken into account
	
      Xsec_Yint=integralOverYpoint_S(var,process,incCut,CutParam,ymin,ymax)
      
      end if
  end if
  
  
  end function Xsec_Yint
  !--------------Simpsons--------------------
  !!!! parameter valueMax remembers the approximate value of integral to weight the tolerance.
  !!!! evaluation is done by adaptive simpson
  !!!! First we evaluate over 5 points and estimate the integral, and then split it to 3+3 and send to adaptive
  !!!! Thus minimal number of points =9
  function integralOverYpoint_S(var,process,incCut,CutParam,yMin_in,yMax_in)
   real(dp),dimension(1:7)::var
   logical,intent(in)::incCut
   real(dp),dimension(1:4),intent(in)::CutParam
   integer,dimension(1:3),intent(in)::process
   real(dp) ::integralOverYpoint_S
   real(dp) :: X1,X2,X3,X4,X5
   real(dp) :: y2,y3,y4,deltay
   real(dp) :: yMin_in,yMax_in
   real(dp)::valueMax
   
   deltay=yMax_in-yMin_in
   y2=yMin_in+deltay/4d0
   y3=yMin_in+deltay/2d0
   y4=yMax_in-deltay/4d0
      
   call SetY(yMin_in,var)
   X1= xSec(var,process,incCut,CutParam)   
   call SetY(y2,var)
   X2= xSec(var,process,incCut,CutParam)   
   call SetY(y3,var)
   X3= xSec(var,process,incCut,CutParam)   
   call SetY(y4,var)
   X4= xSec(var,process,incCut,CutParam)   
   call SetY(yMax_in,var)
   X5= xSec(var,process,incCut,CutParam)
   
   !!approximate integral value
   valueMax=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   integralOverYpoint_S=IntegralOverYpoint_S_Rec(var,process,incCut,CutParam,yMin_in,y3,X1,X2,X3,valueMax)+&
	  IntegralOverYpoint_S_Rec(var,process,incCut,CutParam,y3,yMax_in,X3,X4,X5,valueMax)
  end function integralOverYpoint_S
  
  !!!! X1,X3,X5 are cross-sections at end (X1,X5) and central (X3) points of integraitons
  recursive function integralOverYpoint_S_Rec(var,process,incCut,CutParam,yMin_in,yMax_in,X1,X3,X5,valueMax) result(interX)
   real(dp),dimension(1:7) ::var
   logical,intent(in)::incCut
   real(dp),dimension(1:4),intent(in)::CutParam
   integer,dimension(1:3),intent(in)::process
   real(dp) :: interX,X1,X2,X3,X4,X5
   real(dp) :: valueAB,valueACB
   real(dp) :: yMin_in,yMax_in,y2,y3,y4,deltay
   real(dp),intent(in)::valueMax
   
   deltay=yMax_in-yMin_in
   y2=yMin_in+deltay/4d0
   y3=yMin_in+deltay/2d0
   y4=yMax_in-deltay/4d0
   
   
   call SetY(y2,var)
   X2= xSec(var,process,incCut,CutParam)
   
   call SetY(y4,var)
   X4= xSec(var,process,incCut,CutParam)
   
   valueAB=deltay*(X1+4d0*X3+X5)/6d0
   valueACB=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   !write(*,*) yMin_in,yMax_in,(valueACB-valueAB)/valueMax,ABS((valueACB-valueAB)/valueMax)>tolerance
   !if(GlobalCounter>200) stop
!    if(abs(valueMax)<1d-10) write(*,*) valueMax
   
   If(ABS((valueACB-valueAB)/valueMax)>tolerance) then
    interX=integralOverYpoint_S_Rec(var,process,incCut,CutParam,yMin_in,y3,X1,X2,X3,valueMax)&
	  +integralOverYpoint_S_Rec(var,process,incCut,CutParam,y3,yMax_in,X3,X4,X5,valueMax)
   else
    interX=valueACB
   end if
   
  end function integralOverYpoint_S_Rec
  
  !---------------------------------INTEGRATED over Q---------------------------------------------------------------
  function Xsec_Qint(var,process,incCut,CutParam,Qmin_in,Qmax_in)
    real(dp),dimension(1:7)::var
    logical,intent(in)::incCut
    real(dp),dimension(1:4),intent(in)::CutParam
    integer,dimension(1:3),intent(in)::process
    real(dp),intent(in):: Qmin_in,Qmax_in
    real(dp):: Xsec_Qint
    integer::numSec,i
    real(dp)::dQ
    
    if(TMDF_IsconvergenceLost()) then 
      Xsec_Qint=1d9
      return
    end if

    !!! check how many maxQbins is inside the integration range (+1)
    numSec=INT((Qmax_in-Qmin_in)/maxQbinSize)+1

    if(numSec==1) then
      !!! if the bin is smaller than maxQbinSize, integrate as is
      Xsec_Qint=integralOverQpoint_S(var,process,incCut,CutParam,Qmin_in,Qmax_in)
    else
      !!! else divide to smaler bins and sum the integrals
      dQ=(Qmax_in-Qmin_in)/numSec !!! size of new bins

      Xsec_Qint=0d0
      do i=0,numSec-1
        Xsec_Qint=Xsec_Qint + &
            integralOverQpoint_S(var,process,incCut,CutParam,Qmin_in+i*dQ,Qmin_in+(i+1)*dQ)
      end do
    end if
    

  end function Xsec_Qint
  
  !--------------Simpsons--------------------
  !!!! parameter valueMax remembers the initial value of integral to weight the tolerance.
  !!!! First we evaluate over 5 points and estimate the integral, and then split it to 3+3 and send to adaptive
  !!!! Thus minimal number of points =9
  !!!! taking into account minimum calls of y-integral we have  =81 points
  function integralOverQpoint_S(var,process,incCut,CutParam,QMin_in,QMax_in)
  real(dp),dimension(1:7)::var
  logical,intent(in)::incCut
   real(dp),dimension(1:4),intent(in)::CutParam
   integer,dimension(1:3),intent(in)::process
    real(dp) ::integralOverQpoint_S
   real(dp) :: X1,X2,X3,X4,X5
   real(dp) :: QMin_in,QMax_in
   real(dp)::valueMax,Q2,Q3,Q4,deltaQ
   
   deltaQ=QMax_in-QMin_in
   Q2=QMin_in+deltaQ/4d0
   Q3=QMin_in+deltaQ/2d0
   Q4=QMax_in-deltaQ/4d0
   
   call SetQ(QMin_in,var)
   X1=2*QMin_in*xSec(var,process,incCut,CutParam)
   
   call SetQ(Q2,var)
   X2=2*Q2*xSec(var,process,incCut,CutParam)
   
   call SetQ(Q3,var)
   X3=2*Q3*xSec(var,process,incCut,CutParam)
   
   call SetQ(Q4,var)
   X4=2*Q4*xSec(var,process,incCut,CutParam)
   
   call SetQ(QMax_in,var)
   X5=2*Qmax_in*xSec(var,process,incCut,CutParam)
   
      !!approximate integral value
   valueMax=deltaQ*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   integralOverQpoint_S=IntegralOverQpoint_S_Rec(var,process,incCut,CutParam,QMin_in,Q3,X1,X2,X3,valueMax)+&
	IntegralOverQpoint_S_Rec(var,process,incCut,CutParam,Q3,QMax_in,X3,X4,X5,valueMax)
  end function integralOverQpoint_S
  
  !!!! X1,X3,X5 are cross-sections at end (X1,X5) and central (X3) points of integraitons
  recursive function integralOverQpoint_S_Rec(var,process,incCut,CutParam,QMin_in,QMax_in,X1,X3,X5,valueMax) result(interX)
   real(dp),dimension(1:7)::var
   logical,intent(in)::incCut
   real(dp),dimension(1:4),intent(in)::CutParam
   integer,dimension(1:3),intent(in)::process
   real(dp) :: interX,X1,X2,X3,X4,X5
   real(dp) :: valueAB,valueACB
   real(dp) :: QMin_in,QMax_in,Q2,Q3,Q4,deltaQ
   real(dp),intent(in)::valueMax
   
   deltaQ=QMax_in-QMin_in
   Q2=QMin_in+deltaQ/4d0
   Q3=QMin_in+deltaQ/2d0
   Q4=QMax_in-deltaQ/4d0
   
   call SetQ(Q2,var)
   X2=2*Q2*xSec(var,process,incCut,CutParam)
      
   call SetQ(Q4,var)
   X4=2*Q4*xSec(var,process,incCut,CutParam)
   
   valueAB=deltaQ*(X1+4d0*X3+X5)/6d0
   valueACB=deltaQ*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   If(ABS((valueACB-valueAB)/valueMax)>tolerance) then
    interX=integralOverQpoint_S_Rec(var,process,incCut,CutParam,QMin_in,Q3,X1,X2,X3,valueMax)&
	  +integralOverQpoint_S_Rec(var,process,incCut,CutParam,Q3,Qmax_in,X3,X4,X5,valueMax)
   else
    interX=valueACB
   end if
  end function integralOverQpoint_S_Rec
  
  
  !---------------------------------INTEGRATED over Y over Q---------------------------------------------------------------
  !!!! No need for check over Y they take a place within y-integration

  !!!! to integrate over Q I use adaptive Simpson. (defined below)
  !!!! before the integration over Q I check the size of Q-bin,
  !!!! if Q-bin is large I split desect the integration range to smaller
  function Xsec_Qint_Yint(var,process,incCut,CutParam,Qmin_in,Qmax_in,ymin_in,ymax_in)
  real(dp),dimension(1:7)::var
  logical,intent(in)::incCut
   real(dp),dimension(1:4),intent(in)::CutParam
   integer,dimension(1:3),intent(in)::process
   real(dp),intent(in) :: yMin_in,yMax_in,QMin_in,QMax_in
   real(dp):: Xsec_Qint_Yint
   integer::numSec,i
   real(dp)::dQ

    if(TMDF_IsconvergenceLost()) then
      Xsec_Qint_Yint=1d9
      return
    end if

    !!! check how many maxQbins is inside the integration range (+1)
    numSec=INT((Qmax_in-Qmin_in)/maxQbinSize)+1

    if(numSec==1) then
      !!! if the bin is smaller than maxQbinSize, integrate as is
      Xsec_Qint_Yint=Xsec_Qint_Yint_in(var,process,incCut,CutParam,Qmin_in,Qmax_in,ymin_in,ymax_in)
    else
      !!! else divide to smaler bins and sum the integrals
      dQ=(Qmax_in-Qmin_in)/numSec !!! size of new bins

      Xsec_Qint_Yint=0d0
      do i=0,numSec-1
        Xsec_Qint_Yint=Xsec_Qint_Yint + &
            Xsec_Qint_Yint_in(var,process,incCut,CutParam,Qmin_in+i*dQ,Qmin_in+(i+1)*dQ,ymin_in,ymax_in)
      end do
    end if
  end function Xsec_Qint_Yint

  !--------------Simpsons--------------------
  !!!! parameter valueMax remembers the initial value of integral to weight the tolerance.
  !!!! First we evaluate over 5 points and estimate the integral, and then split it to 3+3 and send to adaptive
  !!!! Thus minimal number of points =9
  !!!! taking into account minimum calls of y-integral we have  =81 points
  function Xsec_Qint_Yint_in(var,process,incCut,CutParam,Qmin_in,Qmax_in,ymin_in,ymax_in)
  real(dp),dimension(1:7)::var
  logical,intent(in)::incCut
   real(dp),dimension(1:4),intent(in)::CutParam
   integer,dimension(1:3),intent(in)::process
   real(dp),intent(in) :: yMin_in,yMax_in,QMin_in,QMax_in
   real(dp):: Xsec_Qint_Yint_in
   real(dp) :: X1,X2,X3,X4,X5
   real(dp)::valueMax,Q2,Q3,Q4,deltaQ
   
    deltaQ=QMax_in-QMin_in
   Q2=QMin_in+deltaQ/4d0
   Q3=QMin_in+deltaQ/2d0
   Q4=QMax_in-deltaQ/4d0
   
   call SetQ(QMin_in,var)
   X1=2*QMin_in*Xsec_Yint(var,process,incCut,CutParam,yMin_in,yMax_in)
   
   call SetQ(Q2,var)
   X2=2*Q2*Xsec_Yint(var,process,incCut,CutParam,yMin_in,yMax_in)
   
   call SetQ(Q3,var)
   X3=2*Q3*Xsec_Yint(var,process,incCut,CutParam,yMin_in,yMax_in)
   
   call SetQ(Q4,var)
   X4=2*Q4*Xsec_Yint(var,process,incCut,CutParam,yMin_in,yMax_in)
   
   call SetQ(QMax_in,var)
   X5=2*QMax_in*Xsec_Yint(var,process,incCut,CutParam,yMin_in,yMax_in)
   
      !!approximate integral value
   valueMax=deltaQ*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   Xsec_Qint_Yint_in=IntegralOverQYpoint_S_Rec(var,process,incCut,CutParam,QMin_in,Q3,yMin_in,yMax_in,X1,X2,X3,valueMax)+&
	IntegralOverQYpoint_S_Rec(var,process,incCut,CutParam,Q3,QMax_in,yMin_in,yMax_in,X3,X4,X5,valueMax)
  end function Xsec_Qint_Yint_in
  
  !!!! X1,X3,X5 are cross-sections at end (X1,X5) and central (X3) points of integraitons
  recursive function integralOverQYpoint_S_Rec(var,process,incCut,CutParam,&
			      QMin_in,QMax_in,yMin_in,yMax_in,X1,X3,X5,valueMax) result(interX)
   real(dp),dimension(1:7)::var
   logical,intent(in)::incCut
   real(dp),dimension(1:4),intent(in)::CutParam
   integer,dimension(1:3),intent(in)::process
   real(dp) :: interX,X1,X2,X3,X4,X5
   real(dp) :: valueAB,valueACB
   real(dp) :: yMin_in,yMax_in,QMin_in,QMax_in,Q2,Q3,Q4,deltaQ
   real(dp),intent(in)::valueMax
   
   deltaQ=QMax_in-QMin_in
   Q2=QMin_in+deltaQ/4d0
   Q3=QMin_in+deltaQ/2d0
   Q4=QMax_in-deltaQ/4d0
   
   call SetQ(Q2,var)
   X2=2*Q2*Xsec_Yint(var,process,incCut,CutParam,yMin_in,yMax_in)
      
   call SetQ(Q4,var)
   X4=2*Q4*Xsec_Yint(var,process,incCut,CutParam,yMin_in,yMax_in)
   
   valueAB=deltaQ*(X1+4d0*X3+X5)/6d0
   valueACB=deltaQ*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
   If(ABS((valueACB-valueAB)/valueMax)>tolerance) then
    interX=integralOverQYpoint_S_Rec(var,process,incCut,CutParam,QMin_in,Q3,yMin_in,yMax_in,X1,X2,X3,valueMax)&
	  +integralOverQYpoint_S_Rec(var,process,incCut,CutParam,Q3,Qmax_in,yMin_in,yMax_in,X3,X4,X5,valueMax)
   else
    interX=valueACB
   end if
  end function integralOverQYpoint_S_Rec
  
  !---------------------------------INTEGRATED over Y over Q over pT-------------------------------------------------------------
  !!!integration over PT is made by Num-sections
  !!!N even  
  function Xsec_PTint_Qint_Yint(process,incCut,CutParam,s_in,qt_min,qt_max,Q_min,Q_max,ymin_in,ymax_in,Num)
    real(dp),dimension(1:7)::var
    logical,intent(in)::incCut
    real(dp),dimension(1:4),intent(in)::CutParam
    integer,dimension(1:3),intent(in)::process
    real(dp):: Xsec_PTint_Qint_Yint,X0,Xfin
    real(dp):: ymin_in,ymax_in,Q_min,Q_max,qt_min,qt_max,s_in
    integer :: Num
    
    if(TMDF_IsconvergenceLost()) then 
      Xsec_PTint_Qint_Yint=1d9
      return
    end if
    
    if(qt_min<1d-3) then
      var=kinematicArray(1d-3,s_in,(Q_min+Q_max)/2d0,(ymin_in+ymax_in)/2d0)
    else
      var=kinematicArray(qt_min,s_in,(Q_min+Q_max)/2d0,(ymin_in+ymax_in)/2d0)
    end if
    
    X0=2d0*qt_min*Xsec_Qint_Yint(var,process,incCut,CutParam,Q_min,Q_max,ymin_in,ymax_in)
    
    call Xsec_PTint_Qint_Yint_0(process,incCut,CutParam,s_in,qt_min,qt_max,Q_min,Q_max,ymin_in,ymax_in,Num,Xfin,X0)
    Xsec_PTint_Qint_Yint=Xfin
    
  end function Xsec_PTint_Qint_Yint
  
  
  !!!integration over PT is made by Num-sections
  !!!N even
  !!! X0 is value of the function at qt_min input
  !!! X0 is value of the function at qt_max output
  !!! !!! Xfin is value of the cross-section
  subroutine Xsec_PTint_Qint_Yint_0(process,incCut,CutParam,s_in,qt_min,qt_max,Q_min,Q_max,ymin_in,ymax_in,Num,Xfin,X0)
    real(dp),dimension(1:7)::var
    logical,intent(in)::incCut
    real(dp),dimension(1:4),intent(in)::CutParam
    integer,dimension(1:3),intent(in)::process
    real(dp):: Xfin,X0
    real(dp):: ymin_in,ymax_in,Q_min,Q_max,qt_min,qt_max,s_in
    integer :: i,Num
    
    real(dp)::deltaQT,qT_cur,inter
    
    if(mod(num,2)>0) then 
      write(*,*) 'ERROR: arTeMiDe_DY: number of Simpson sections is odd. Evaluation stop.'
      stop
    end if
    
    deltaQT=(qt_max-qt_min)/Num
    inter=X0!!!first term is calculated eqarlier
    
    var=kinematicArray(qt_min,s_in,(Q_min+Q_max)/2d0,(ymin_in+ymax_in)/2d0)
    
    !!!! even terms
    do i=1,Num-1,2
    qT_cur=qt_min+i*deltaQT
    call SetQT(qT_cur,var)
    inter=inter+8d0*qt_cur*Xsec_Qint_Yint(var,process,incCut,CutParam,Q_min,Q_max,ymin_in,ymax_in)
    end do
    
    if(Num>2) then
    !!!! odd terms
    do i=2,Num-2,2
    qT_cur=qt_min+i*deltaQT
    call SetQT(qT_cur,var)
    inter=inter+4d0*qt_cur*Xsec_Qint_Yint(var,process,incCut,CutParam,Q_min,Q_max,ymin_in,ymax_in)
    end do
    end if
    
    call SetQT(qT_max,var)
    X0=2d0*qt_max*Xsec_Qint_Yint(var,process,incCut,CutParam,Q_min,Q_max,ymin_in,ymax_in)!!!! last term
    inter=inter+X0  
    
    Xfin=deltaQT/3d0*inter
    
  end subroutine Xsec_PTint_Qint_Yint_0
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!INTERFACES TO CALCULATING CROSS-SECTIONS!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !---------------------------------UN-INTEGRATED------------------------------------------------------------------
  
  !!qt_list is the list of requred qt -point,
  !! X_list is variable to store results (should be of the same dimension as qt_list)
  subroutine xSecList(X_list,qt_List)
    real(dp), intent(in) :: qt_list(:)
    real(dp), intent(out) :: X_list(:)
    integer :: i,length
    real(dp),dimension(1:7)::var
    length=size(qt_list)
    CallCounter=CallCounter+length
     do i=1,length
       var=kinematicArray(qt_List(i),s_global,Q_global,y_global)
       X_List(i)=PreFactor1(process_global(1))*xSec(var,process_global,includeCuts_global,CutParameters_global)
     end do
  end subroutine xSecList
  
  !!!!Evaluate differential xSec at single point
  subroutine xSecSingle(X,qT_in)
   real(dp),dimension(1:7)::var
   real(dp):: X,qT_in
   CallCounter=CallCounter+1
   var=kinematicArray(qt_in,s_global,Q_global,y_global)
   
   X=PreFactor1(process_global(1))*xSec(var,process_global,includeCuts_global,CutParameters_global)
  end subroutine xSecSingle
  
  
!---------------------------------INTEGRATED over Y---------------------------------------------------------------
    !!qt_list is the list of requred qt -point,
  !! X_list is variable to store results (should be of the same dimension as qt_list)
  subroutine xSecList_Yint(X_list,qt_List,yMin_in,yMax_in)
    real(dp), intent(in) :: qt_list(:)
    real(dp), intent(out) :: X_list(:)
    real(dp)::yMin_in,yMax_in
    real(dp),dimension(1:7)::var
    integer :: i,length
    length=size(qt_list)
    CallCounter=CallCounter+length
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(var)
     do i=1,length
	var=kinematicArray(qt_List(i),s_global,Q_global,y_global)
       X_List(i)=PreFactor1(process_global(1))*Xsec_Yint(var,process_global,includeCuts_global,CutParameters_global,yMin_in,yMax_in)
     end do
     !$OMP END PARALLEL DO
  end subroutine xSecList_Yint
  
  !!
  subroutine xSecSingle_Yint(X,qt,yMin_in,yMax_in)
    real(dp)::X,qT
    real(dp)::yMin_in,yMax_in
    real(dp),dimension(1:7)::var
    
   CallCounter=CallCounter+1
   var=kinematicArray(qt,s_global,Q_global,y_global)
   X=PreFactor1(process_global(1))*Xsec_Yint(var,process_global,includeCuts_global,CutParameters_global,yMin_in,yMax_in)
  end subroutine xSecSingle_Yint
  
      !!qt_list is the list of requred qt -point,
  !! X_list is variable to store results (should be of the same dimension as qt_list)
  subroutine xSecList_Ycomplete(X_list,qt_List)
    real(dp), intent(in) :: qt_list(:)
    real(dp), intent(out) :: X_list(:)
    integer :: i,length
    real(dp),dimension(1:7)::var
    length=size(qt_list)
    CallCounter=CallCounter+length
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(var)
     do i=1,length
       var=kinematicArray(qt_List(i),s_global,Q_global,y_global)
       X_List(i)=PreFactor1(process_global(1))*Xsec_Yint(var,process_global,includeCuts_global,CutParameters_global,&
		    log(var(5)),-log(var(5)))
     end do
     !$OMP END PARALLEL DO
  end subroutine xSecList_Ycomplete
  
  subroutine xSecSingle_Ycomplete(X,qt)
    real(dp)::X,qT
    real(dp),dimension(1:7)::var
    
   CallCounter=CallCounter+1
   var=kinematicArray(qt,s_global,Q_global,y_global)
   X=PreFactor1(process_global(1))*Xsec_Yint(var,process_global,includeCuts_global,CutParameters_global,log(var(5)),-log(var(5)))
  end subroutine xSecSingle_Ycomplete
  
 
  !---------------------------------INTEGRATED over Q---------------------------------------------------------------
  !!qt_list is the list of requred qt -point,
  !! X_list is variable to store results (should be of the same dimension as qt_list)
  subroutine xSecList_Qint(X_list,qt_List,Q_min,Q_max)
    real(dp), intent(in) :: qt_list(:)
    real(dp), intent(out) :: X_list(:)
    real(dp)::Q_min,Q_max
    integer :: i,length
    real(dp),dimension(1:7)::var
    length=size(qt_list)
    CallCounter=CallCounter+length
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(var)
     do i=1,length
       var=kinematicArray(qt_List(i),s_global,Q_global,y_global)
       X_List(i)=PreFactor1(process_global(1))*Xsec_Qint(var,process_global,includeCuts_global,CutParameters_global,Q_min,Q_max)
     end do
     !$OMP END PARALLEL DO
  end subroutine xSecList_Qint
  
  !!
  subroutine xSecSingle_Qint(X,qt,Q_min,Q_max)
    real(dp)::X,qT
    real(dp)::Q_min,Q_max
    real(dp),dimension(1:7)::var
    
   CallCounter=CallCounter+1
   var=kinematicArray(qt,s_global,Q_global,y_global)
   X=PreFactor1(process_global(1))*Xsec_Qint(var,process_global,includeCuts_global,CutParameters_global,Q_min,Q_max)
  end subroutine xSecSingle_Qint
  
  
!---------------------------------INTEGRATED over Y over Q---------------------------------------------------------------
  !!qt_list is the list of requred qt -point,
  !! X_list is variable to store results (should be of the same dimension as qt_list)
  subroutine xSecList_Qint_Yint(X_list,qt_List,Q_min,Q_max,yMin_in,yMax_in)
    real(dp), intent(in) :: qt_list(:)
    real(dp), intent(out) :: X_list(:)
    real(dp)::yMin_in,yMax_in,Q_min,Q_max
    integer :: i,length
    real(dp),dimension(1:7)::var
    length=size(qt_list)
    CallCounter=CallCounter+length
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(var)
     do i=1,length
       var=kinematicArray(qt_List(i),s_global,Q_global,y_global)
       X_List(i)=PreFactor1(process_global(1))*Xsec_Qint_Yint(var,process_global,includeCuts_global,CutParameters_global,&
		Q_min,Q_max,yMin_in,yMax_in)
     end do
     !$OMP END PARALLEL DO
  end subroutine xSecList_Qint_Yint
  
  !!
  subroutine xSecSingle_Qint_Yint(X,qt,Q_min,Q_max,yMin_in,yMax_in)
    real(dp)::X,qT
    real(dp)::yMin_in,yMax_in,Q_min,Q_max
    real(dp),dimension(1:7)::var
   CallCounter=CallCounter+1
   var=kinematicArray(qt,s_global,Q_global,y_global)
   X=PreFactor1(process_global(1))*Xsec_Qint_Yint(var,process_global,includeCuts_global,CutParameters_global,&
	Q_min,Q_max,yMin_in,yMax_in)
  end subroutine xSecSingle_Qint_Yint
  
      !!qt_list is the list of requred qt -point,
  !! X_list is variable to store results (should be of the same dimension as qt_list)
  subroutine xSecList_Qint_Ycomplete(X_list,qt_List,Q_min,Q_max)
    real(dp), intent(in) :: qt_list(:)
    real(dp), intent(out) :: X_list(:)
    real(dp)::Q_min,Q_max
    integer :: i,length
    real(dp),dimension(1:7)::var
    length=size(qt_list)
    CallCounter=CallCounter+length
    !$OMP PARALLEL DO DEFAULT(SHARED) PRIVATE(var)
     do i=1,length
       var=kinematicArray(qt_List(i),s_global,Q_global,y_global)
       X_List(i)=PreFactor1(process_global(1))*Xsec_Qint_Yint(var,process_global,includeCuts_global,CutParameters_global,&
			      Q_min,Q_max,log(var(5)),-log(var(5)))
     end do
     !$OMP END PARALLEL DO
  end subroutine xSecList_Qint_Ycomplete
  
  subroutine xSecSingle_Qint_Ycomplete(X,qt,Q_min,Q_max)
      real(dp)::X,qT
      real(dp)::Q_min,Q_max
      real(dp),dimension(1:7)::var
    
   CallCounter=CallCounter+1
   var=kinematicArray(qt,s_global,Q_global,y_global)
   X=PreFactor1(process_global(1))*Xsec_Qint_Yint(var,process_global,includeCuts_global,CutParameters_global,&
	      Q_min,Q_max,log(var(5)),-log(var(5)))
  end subroutine xSecSingle_Qint_Ycomplete
  
  !---------------------------------INTEGRATED over Y over Q  over PT----------------------------------------------------------
  !!qt_list is the list of requred qt -point,
  !! X_list is variable to store results (should be of the same dimension as qt_list)
  subroutine xSecList_PTintN_Qint_Yint(X_list,qtMIN_List,qtMAX_list,Q_min,Q_max,yMin_in,yMax_in,num)
    real(dp), intent(in) :: qtMIN_list(:),qtMAX_list(:)
    real(dp), intent(out) :: X_list(:)
    real(dp)::yMin_in,yMax_in,Q_min,Q_max
    integer :: i,length,num
    length=size(qtMIN_list)
    CallCounter=CallCounter+length
    !$OMP PARALLEL DO DEFAULT(SHARED)
     do i=1,length
       X_List(i)=PreFactor1(process_global(1))*Xsec_PTint_Qint_Yint(process_global,includeCuts_global,CutParameters_global,&
				s_global,qtMIN_List(i),qtMAX_list(i),Q_min,Q_max,yMin_in,yMax_in,num)
     end do
     !$OMP END PARALLEL DO
  end subroutine xSecList_PTintN_Qint_Yint
  
  subroutine xSecList_PTint_Qint_Yint(X_list,qtMIN_List,qtMAX_list,Q_min,Q_max,yMin_in,yMax_in)
    real(dp), intent(in) :: qtMIN_list(:),qtMAX_list(:)
    real(dp), intent(out) :: X_list(:)
    real(dp)::yMin_in,yMax_in,Q_min,Q_max
    
    call xSecList_PTintN_Qint_Yint(X_list,qtMIN_List,qtMAX_list,Q_min,Q_max,yMin_in,yMax_in,NumPTdefault)

  end subroutine xSecList_PTint_Qint_Yint
  
  subroutine xSecListList_PTintN_Qint_Yint(X_list,qt_List,Q_min,Q_max,yMin_in,yMax_in,num)
    real(dp), intent(in) :: qt_list(:)
    real(dp), intent(out) :: X_list(:)
    real(dp)::Q_min,Q_max,X0,Xfin,yMin_in,yMax_in
    integer :: i,length,length2,num
    real(dp),dimension(1:7)::var
    length2=size(qt_list)
    length=size(X_list)
    if( (length2-length) .ne. 1) then    
      write(*,*) 'ERROR: arTeMiDe_DY: pt integration : sizes of lists (for cross-seciont and pt-bins) are inconsistent'
      write(*,*) 'Evaluation stop'
      stop
    end if 
    
    !-----------------parallel version-------------
    !$ call xSecList_PTintN_Qint_Yint(X_list,qt_List(1:length2-1),qt_list(2:length2),Q_min,Q_max,yMin_in,yMax_in,num)
    !$ return
    
    !-----------------single-thread version-------------
    CallCounter=CallCounter+length
    
    var=kinematicArray(qt_list(1),s_global,Q_global,y_global)
    X0=2d0*qt_list(1)*Xsec_Qint_Yint(var,process_global,includeCuts_global,CutParameters_global,Q_min,Q_max,yMin_in,yMax_in)
    do i=1,length
       call Xsec_PTint_Qint_Yint_0(process_global,includeCuts_global,CutParameters_global,&
		    s_global,qt_list(i),qt_list(i+1),Q_min,Q_max,yMin_in,yMax_in,Num,Xfin,X0)
       X_List(i)=PreFactor1(process_global(1))*Xfin
     end do
  end subroutine xSecListList_PTintN_Qint_Yint
  
  subroutine xSecListList_PTint_Qint_Yint(X_list,qt_List,Q_min,Q_max,yMin_in,yMax_in)
    real(dp), intent(in) :: qt_list(:)
    real(dp), intent(out) :: X_list(:)
    real(dp)::Q_min,Q_max,yMin_in,yMax_in
    
    call xSecListList_PTintN_Qint_Yint(X_list,qt_List,Q_min,Q_max,yMin_in,yMax_in,NumPTdefault)
    
  end subroutine xSecListList_PTint_Qint_Yint
  
  !!
  subroutine xSecSingle_PTintN_Qint_Yint(X,qt_Min,qt_Max,Q_min,Q_max,yMin_in,yMax_in,num)
    real(dp)::X,qt_Min,qt_Max
    real(dp)::yMin_in,yMax_in,Q_min,Q_max
    integer::num
    
   CallCounter=CallCounter+1
   X=PreFactor1(process_global(1))*Xsec_PTint_Qint_Yint(process_global,includeCuts_global,CutParameters_global,&
				  s_global,qt_Min,qt_Max,Q_min,Q_max,yMin_in,yMax_in,num)
  end subroutine xSecSingle_PTintN_Qint_Yint
  
  subroutine xSecSingle_PTint_Qint_Yint(X,qt_Min,qt_Max,Q_min,Q_max,yMin_in,yMax_in)
    real(dp)::X,qt_Min,qt_Max
    real(dp)::yMin_in,yMax_in,Q_min,Q_max
   X=PreFactor1(process_global(1))*Xsec_PTint_Qint_Yint(process_global,includeCuts_global,CutParameters_global,&
			    s_global, qt_Min,qt_Max,Q_min,Q_max,yMin_in,yMax_in,NumPTdefault)
  end subroutine xSecSingle_PTint_Qint_Yint
  
  subroutine xSecListPY_PTint_Qint_Yint(X_list,qtMIN_List,qtMAX_list,Q_min,Q_max,yMin_List,yMax_List,num)
    real(dp), intent(in) :: qtMIN_list(:),qtMAX_list(:),yMin_List(:),yMax_List(:)
    real(dp), intent(out) :: X_list(:)
    real(dp)::Q_min,Q_max
    integer :: i,length,num
    length=size(qtMIN_list)
    if(size(qtMAX_list)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY: pt integration (2) : sizes of lists (pt min-max bins) are inconsistent'
      write(*,*) 'Evaluation stop'
      stop
    end if
    
    if(size(yMin_List)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY: pt integration (3) : sizes of lists (pt-bins vs yMin) are inconsistent'
      write(*,*) 'Evaluation stop'
      stop
    end if
    
    if(size(yMax_List)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY: pt integration (4) : sizes of lists (pt-bins vs yMax) are inconsistent'
      write(*,*) 'Evaluation stop'
      stop
    end if
    
    
    CallCounter=CallCounter+length
    !$OMP PARALLEL DO DEFAULT(SHARED)
     do i=1,length
       X_List(i)=PreFactor1(process_global(1))*Xsec_PTint_Qint_Yint(process_global,includeCuts_global,CutParameters_global,&
			  s_global,qtMIN_List(i),qtMAX_list(i),Q_min,Q_max,yMin_List(i),yMax_List(i),num)
     end do
     !$OMP END PARALLEL DO
  end subroutine xSecListPY_PTint_Qint_Yint
  
  subroutine xSecListPY_PTintN_Qint_Yint(X_list,qtMIN_List,qtMAX_list,Q_min,Q_max,yMin_List,yMax_List)
    real(dp), intent(in) :: qtMIN_list(:),qtMAX_list(:),yMin_List(:),yMax_List(:)
    real(dp), intent(out) :: X_list(:)
    real(dp)::Q_min,Q_max
    call xSecListPY_PTint_Qint_Yint(X_list,qtMIN_List,qtMAX_list,Q_min,Q_max,yMin_List,yMax_List,NumPTdefault)
  end subroutine xSecListPY_PTintN_Qint_Yint
  
  !---------------------------------INTEGRATED over Y (complete) over Q  over PT---------------------------------------------------!
  
  !!qt_list is the list of requred qt -point,
  !! X_list is variable to store results (should be of the same dimension as qt_list)
  !! I set y in (-1000,1000) since the check is made in the integration routine
  subroutine xSecList_PTintN_Qint_Ycomplete(X_list,qtMIN_List,qtMAX_list,Q_min,Q_max,num)
    real(dp), intent(in) :: qtMIN_list(:),qtMAX_list(:)
    real(dp), intent(out) :: X_list(:)
    real(dp)::Q_min,Q_max
    integer :: i,length,num
    length=size(qtMIN_list)
    CallCounter=CallCounter+length
    !$OMP PARALLEL DO DEFAULT(SHARED)
     do i=1,length
       X_List(i)=PreFactor1(process_global(1))*Xsec_PTint_Qint_Yint(process_global,includeCuts_global,CutParameters_global,&
				s_global,qtMIN_list(i),qtMAX_list(i),Q_min,Q_max,-1000d0,1000d0,num)
     end do
    !$OMP END PARALLEL DO
  end subroutine xSecList_PTintN_Qint_Ycomplete
  
  subroutine xSecSingle_PTintN_Qint_Ycomplete(X,qt_min,qt_max,Q_min,Q_max,num)
      real(dp)::X,qT_min,qT_max
      real(dp)::Q_min,Q_max
      integer::num
    
   CallCounter=CallCounter+1
   X=PreFactor1(process_global(1))*Xsec_PTint_Qint_Yint(process_global,includeCuts_global,CutParameters_global,&
			    s_global,qt_min,qt_max,Q_min,Q_max,-1000d0,1000d0,num)
  end subroutine xSecSingle_PTintN_Qint_Ycomplete
  
  subroutine xSecListList_PTintN_Qint_Ycomplete(X_list,qt_List,Q_min,Q_max,num)
    real(dp), intent(in) :: qt_list(:)
    real(dp), intent(out) :: X_list(:)
    real(dp)::Q_min,Q_max,X0,Xfin
    integer :: i,length,length2,num
    real(dp),dimension(1:7)::var
    
    length2=size(qt_list)
    length=size(X_list)
    if( (length2-length) .ne. 1) then    
      write(*,*) 'ERROR: arTeMiDe_DY: pt integration : sizes of lists (for cross-seciont and pt-bins) are inconsistent'
      write(*,*) 'Evaluation stop'
      stop
    end if  
    
    !-----------------parallel version-------------
    !$ call xSecList_PTintN_Qint_Ycomplete(X_list,qt_List(1:length2-1),qt_list(2:length2),Q_min,Q_max,num)
    !$ return
    
    !-----------------single-thread version-------------
    CallCounter=CallCounter+length
    
    var=kinematicArray(qt_list(1),s_global,Q_global,y_global)
    X0=2d0*qt_list(1)*Xsec_Qint_Yint(var,process_global,includeCuts_global,CutParameters_global,Q_min,Q_max,-1000d0,1000d0)
    do i=1,length
       call Xsec_PTint_Qint_Yint_0(process_global,includeCuts_global,CutParameters_global,&
		      s_global, qt_list(i),qt_list(i+1),Q_min,Q_max,-1000d0,1000d0,Num,Xfin,X0)
       X_List(i)=PreFactor1(process_global(1))*Xfin
     end do
  end subroutine xSecListList_PTintN_Qint_Ycomplete
  
  subroutine xSecList_PTint_Qint_Ycomplete(X_list,qtMIN_List,qtMAX_list,Q_min,Q_max)
    real(dp), intent(in) :: qtMIN_list(:),qtMAX_list(:)
    real(dp), intent(out) :: X_list(:)
    real(dp)::Q_min,Q_max
    
    call xSecList_PTintN_Qint_Ycomplete(X_list,qtMIN_List,qtMAX_list,Q_min,Q_max,NumPTdefault)
  end subroutine xSecList_PTint_Qint_Ycomplete
  
  subroutine xSecSingle_PTint_Qint_Ycomplete(X,qt_min,qt_max,Q_min,Q_max)
      real(dp)::X,qT_min,qT_max
      real(dp)::Q_min,Q_max
    
   CallCounter=CallCounter+1
   X=PreFactor1(process_global(1))*Xsec_PTint_Qint_Yint(process_global,includeCuts_global,CutParameters_global,&
			    s_global,qt_min,qt_max,Q_min,Q_max,-1000d0,1000d0,NumPTdefault)
  end subroutine xSecSingle_PTint_Qint_Ycomplete
  
  subroutine xSecListList_PTint_Qint_Ycomplete(X_list,qt_List,Q_min,Q_max)
    real(dp), intent(in) :: qt_list(:)
    real(dp), intent(out) :: X_list(:)
    real(dp)::Q_min,Q_max
    
    call xSecListList_PTintN_Qint_Ycomplete(X_list,qt_List,Q_min,Q_max,NumPTdefault)
    
  end subroutine xSecListList_PTint_Qint_Ycomplete
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!THE MAIN INTERFACE TO CROSS-SECTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!! interface for integer,s,array,array,array,logical,optional, optional
  subroutine MainInterface_isAAAloo(X,process,s,qT,Q,y,includeCuts,CutParameters,Num)
    integer,intent(in)::process					!the number of process
    real(dp),intent(in)::s					!Mandelshtam s
    real(dp),intent(in),dimension(1:2)::qT			!(qtMin,qtMax)
    real(dp),intent(in),dimension(1:2)::Q				!(Qmin,Qmax)
    real(dp),intent(in),dimension(1:2)::y				!(ymin,ymax)
    logical,intent(in)::includeCuts				!include cuts
    real(dp),intent(in),dimension(1:4),optional::CutParameters	!(p1,p2,eta1,eta2)
    integer,intent(in),optional::Num				!number of sections
    
    real(dp)::X
  
  integer::nn
  real(dp),dimension(1:4)::CutParam
  integer,dimension(1:3)::ppp
  
  !!! determine number of sections
  if(present(Num)) then
    nn=Num
  else
    nn=NumPTdefault
  end if
  
  !!! determine cut parameters
  if(includeCuts) then
    if(present(CutParameters)) then
      CutParam=CutParameters
    else
      write(*,*) 'ERROR: arTeMiDe_DY: called includeCuts=true, while CutParameters are undefined'
      write(*,*) 'Evaluation stop'
      stop
    end if
  else
    CutParam=(/0d0,0d0,0d0,0d0/)
  end if
  
  ppp=processArrayFromInteger(process)
  
  !!!! evaluation
  CallCounter=CallCounter+1
  X=PreFactor1(ppp(1))*Xsec_PTint_Qint_Yint(ppp,includeCuts,CutParameters,&
				  s,qT(1),qT(2),Q(1),Q(2),y(1),y(2),nn)
  
  end subroutine MainInterface_isAAAloo
  
  !!!! interface for array,s,array,array,array,logical,optional, optional
  subroutine MainInterface_AsAAAloo(X,process,s,qT,Q,y,includeCuts,CutParameters,Num)
!   function xSec_DY(process,s,qT,Q,y,includeCuts,CutParameters,Num)
    integer,intent(in),dimension(1:3)::process			!the number of process
    real(dp),intent(in)::s					!Mandelshtam s
    real(dp),intent(in),dimension(1:2)::qT			!(qtMin,qtMax)
    real(dp),intent(in),dimension(1:2)::Q				!(Qmin,Qmax)
    real(dp),intent(in),dimension(1:2)::y				!(ymin,ymax)
    logical,intent(in)::includeCuts				!include cuts
    real(dp),intent(in),dimension(1:4),optional::CutParameters	!(p1,p2,eta1,eta2)
    integer,intent(in),optional::Num				!number of sections
    
    real(dp)::X
  
  integer::nn
  real(dp),dimension(1:4)::CutParam
  
  
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
      write(*,*) 'ERROR: arTeMiDe_DY: called includeCuts=true, while CutParameters are undefined'
      write(*,*) 'Evaluation stop'
      stop
    end if
  else
    CutParam=(/0d0,0d0,0d0,0d0/)
  end if
  
  !!!! evaluation
  CallCounter=CallCounter+1
  X=PreFactor1(process(1))*Xsec_PTint_Qint_Yint(process,includeCuts,CutParameters,&
				  s,qT(1),qT(2),Q(1),Q(2),y(1),y(2),nn)
  
  end subroutine MainInterface_AsAAAloo
  
  subroutine xSec_DY_List(X,process,s,qT,Q,y,includeCuts,CutParameters,Num)
    integer,intent(in),dimension(:,:)::process			!the number of process
    real(dp),intent(in),dimension(:)::s				!Mandelshtam s
    real(dp),intent(in),dimension(:,:)::qT			!(qtMin,qtMax)
    real(dp),intent(in),dimension(:,:)::Q				!(Qmin,Qmax)
    real(dp),intent(in),dimension(:,:)::y				!(ymin,ymax)
    logical,intent(in),dimension(:)::includeCuts		!include cuts
    real(dp),intent(in),dimension(:,:)::CutParameters	        !(p1,p2,eta1,eta2)
    integer,intent(in),dimension(:),optional::Num		!number of sections
    real(dp),dimension(:),intent(out)::X
    integer :: i,length
    integer,allocatable::nn(:)
    
    length=size(s)
    
    !!! cheking sizes
    if(size(X)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY: xSec_DY_List: sizes of xSec and s lists are not equal.'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(process,1)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY: xSec_DY_List: sizes of process and s lists are not equal.'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(qT,1)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY: xSec_DY_List: sizes of qT and s lists are not equal.'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(y,1)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY: xSec_DY_List: sizes of y and s lists are not equal.'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(Q,1)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY: xSec_DY_List: sizes of Q and s lists are not equal.'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(includeCuts)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY: xSec_DY_List: sizes of includeCuts and s lists are not equal.'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(CutParameters,1)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY: xSec_DY_List: sizes of CutParameters and s lists are not equal.'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(process,2)/=3) then
      write(*,*) 'ERROR: arTeMiDe_DY: xSec_DY_List: process list must be (:,1:3).'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(qT,2)/=2) then
      write(*,*) 'ERROR: arTeMiDe_DY: xSec_DY_List: qt list must be (:,1:3).'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(y,2)/=2) then
      write(*,*) 'ERROR: arTeMiDe_DY: xSec_DY_List: y list must be (:,1:3).'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(Q,2)/=2) then
      write(*,*) 'ERROR: arTeMiDe_DY: xSec_DY_List: Q list must be (:,1:3).'
      write(*,*) 'Evaluation stop'
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
       X(i)=PreFactor1(process(i,1))*Xsec_PTint_Qint_Yint(process(i,1:3),includeCuts(i),CutParameters(i,1:4),&
				s(i),qT(i,1),qT(i,2),Q(i,1),Q(i,2),y(i,1),y(i,2),nn(i))
     end do
    !$OMP END PARALLEL DO
    deallocate(nn)
  end subroutine xSec_DY_List
  
  subroutine xSec_DY_List_BINLESS(X,process,s,qT,Q,y,includeCuts,CutParameters)
    integer,intent(in),dimension(:,:)::process			!the number of process
    real(dp),intent(in),dimension(:)::s				!Mandelshtam s
    real(dp),intent(in),dimension(:)::qT			!(qtMin,qtMax)
    real(dp),intent(in),dimension(:)::Q				!(Qmin,Qmax)
    real(dp),intent(in),dimension(:)::y				!(ymin,ymax)
    logical,intent(in),dimension(:)::includeCuts		!include cuts
    real(dp),intent(in),dimension(:,:)::CutParameters	        !(p1,p2,eta1,eta2)
    real(dp),dimension(:),intent(out)::X
    
    real(dp),allocatable,dimension(:,:)::vv
    integer :: i,length
    
    length=size(s)
    
    !!! cheking sizes
    if(size(X)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY_BINLESS: xSec_DY_List: sizes of xSec and s lists are not equal.'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(process,1)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY_BINLESS: xSec_DY_List: sizes of process and s lists are not equal.'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(qT)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY_BINLESS: xSec_DY_List: sizes of qT and s lists are not equal.'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(y)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY_BINLESS: xSec_DY_List: sizes of y and s lists are not equal.'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(Q)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY_BINLESS: xSec_DY_List: sizes of Q and s lists are not equal.'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(includeCuts)/=length) then
      write(*,*) 'ERROR: arTeMiDe_DY_BINLESS: xSec_DY_List: sizes of includeCuts and s lists are not equal.'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(CutParameters,1)/=length) then
      write(*,*) 'ERROR: arTeMiDe_D_BINLESSY: xSec_DY_List: sizes of CutParameters and s lists are not equal.'
      write(*,*) 'Evaluation stop'
      stop
    end if
    if(size(process,2)/=3) then
      write(*,*) 'ERROR: arTeMiDe_DY_BINLESS: xSec_DY_List: process list must be (:,1:3).'
      write(*,*) 'Evaluation stop'
      stop
    end if    
    
    CallCounter=CallCounter+length
    
    allocate(vv(1:length,1:7))
    
    !$OMP PARALLEL DO DEFAULT(SHARED)
    
     do i=1,length
       vv(i,1:7)=kinematicArray(qt(i),s(i),Q(i),y(i))
       X(i)=PreFactor1(process(i,1))*xSec(vv(i,1:7),process(i,1:3),includeCuts(i),CutParameters(i,1:4))
       
     end do
    !$OMP END PARALLEL DO
    deallocate(vv)
  end subroutine xSec_DY_List_BINLESS
  
end module TMDX_DY
