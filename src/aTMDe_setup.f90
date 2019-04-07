!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe-setup 2.0
!
!	Module that read, create and modify constants-files
!	
!	if you use this module please, quote 1803.11089
!
!				A.Vladimirov (27.05.2019)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module aTMDe_Setup

implicit none

private

character (len=5),parameter :: version="v2.00"
!! actual version of input file
integer,parameter::inputVer=1

!detalization of output: 0 = no output except critical, 1 = + WARNINGS, 2 = + states of initialization,sets,etc, 3 = + details
integer::outputLevel
integer::messageTrigger


!-----------------------physic parameters
logical::include_EWinput
!mass, and width parameters
real*8::mCHARM,mBOTTOM,mZ,mW,GammaZ,GammaW
!other physical parameters
real*8::hc2,alphaQED_MZ,sW2
!CKM matrix
real*8::Vckm_UD,Vckm_US,Vckm_CD,Vckm_CS,Vckm_CB,Vckm_UB

!--------------------- uPDF parameters
integer::number_of_uPDFs
integer,allocatable::enumeration_of_uPDFs(:),replicas_of_uPDFs(:)
character(len=100),allocatable::sets_of_uPDFs(:)

!--------------------- uFF parameters
integer::number_of_uFFs
integer,allocatable::enumeration_of_uFFs(:),replicas_of_uFFs(:)
character(len=100),allocatable::sets_of_uFFs(:)

!-------------------- TMDR options
logical::include_TMDR
character*8::TMDR_order
integer::TMDR_evolutionType,TMDR_lambdaLength
real*8::TMDR_tolerance

!-------------------- uTMDPDF parameters
logical::include_uTMDPDF
integer::uTMDPDF_lambdaLength
real*8::uTMDPDF_tolerance
integer::uTMDPDF_maxIteration
character*8::uTMDPDF_order
logical::uTMDPDF_makeGrid,uTMDPDF_withGluon
real*8::uTMDPDF_grid_bMax,uTMDPDF_grid_xMin,uTMDPDF_grid_slope
integer::uTMDPDF_grid_SizeX,uTMDPDF_grid_SizeB

!-------------------- uTMDFF parameters
logical::include_uTMDFF
integer::uTMDFF_lambdaLength
real*8::uTMDFF_tolerance
integer::uTMDFF_maxIteration
character*8::uTMDFF_order
logical::uTMDFF_makeGrid,uTMDFF_withGluon
real*8::uTMDFF_grid_bMax,uTMDFF_grid_xMin,uTMDFF_grid_slope
integer::uTMDFF_grid_SizeX,uTMDFF_grid_SizeB

!-------------------- TMDs parameters
logical::include_TMDs

!-------------------- TMDF parameters
logical::include_TMDF
real*8::TMDF_OGATAh,TMDF_tolerance

!-------------------- TMDs-inKT parameters
logical::include_TMDs_inKT
real*8::TMDs_inKT_OGATAh,TMDs_inKT_tolerance

!-------------------- TMDX-DY parameters
logical::include_TMDX_DY
character*8::TMDX_DY_order
real*8::TMDX_DY_tolerance
integer::TMDX_DY_ptSECTION
logical::TMDX_DY_exactX1X2
integer::TMDX_DY_numProc

!-------------------- TMDX-SIDIS parameters
logical::include_TMDX_SIDIS
character*8::TMDX_SIDIS_order
real*8::TMDX_SIDIS_tolerance
integer::TMDX_SIDIS_ptSECTION
logical::TMDX_SIDIS_qTcorr,TMDX_SIDIS_M1corr,TMDX_SIDIS_M2corr
integer::TMDX_SIDIS_numProc


!---------------------------------------------------
public::artemide_Setup_Default,artemide_Setup_fromFile,artemide_include,CreateConstantsFile,ReadConstantsFile
public::Set_uPDF,Set_uFF,Set_quarkMasses,Set_EWparameters
public::Set_TMDR_order,Set_TMDR_evolutionType,Set_TMDR_lengthNParray
public::Set_uTMDPDF,Set_uTMDPDF_order,Set_uTMDPDF_gridEvaluation,Set_uTMDPDF_lengthNParray
public::Set_uTMDFF,Set_uTMDFF_order,Set_uTMDFF_gridEvaluation,Set_uTMDFF_lengthNParray
contains
  
  !!!
  subroutine artemide_Setup_Default(order)
    character(len=*),intent(in)::order
    
    outputLevel=2
    call SetupDefault(trim(order))
  
  end subroutine artemide_Setup_Default
  
  subroutine artemide_Setup_fromFile(file,prefix,order)
  character(len=*)::file
  character(len=*),optional::prefix
  character(len=*),optional::order
  
  !! first we set up the default definitions for LO
  outputLevel=2
  if(present(order)) then
  call SetupDefault(trim(order))
  else
  call SetupDefault('LO')
  end if
    
  if(present(prefix)) then
    call ReadConstantsFile(file,prefix)
  else
    call ReadConstantsFile(file)
  end if
  
  end subroutine artemide_Setup_fromFile
  
  !!! This subroutine set all parameters to their defaul value according to given order
  !!! It is called in any case, in order to make sure that all parameters are defined
  subroutine SetupDefault(order)
  
    character(len=*),intent(in)::order
    
    !---------------------Global definitions
    outputLevel=2
    messageTrigger=6
    
    include_EWinput=.true.
    !-----------------------physic parameters
    mCHARM=1.400d0	!threashold mass for charm quark
    mBOTTOM=4.750d0	!threashold mass for bottom quark
    
    mZ=91.1876d0	!pole mass for Z-boson
    GammaZ=2.4952d0	!width of Z-boson
    
    mW=80.379d0		!pole mass for W-boson
    GammaW=2.085d0	!width of W-boson
    
    hc2=0.389379338d0	!transformation constant (hc)^2   GeV->mbarn
    
    alphaQED_MZ=127.955d0	!inverse alpha_QED at Z-boson mass
    
    sW2=0.23122d0	!sin^2 theta_Winberg
    
    !CKM matrix
    Vckm_UD=0.97420d0
    Vckm_US=0.2243d0
    Vckm_CD=0.218d0
    Vckm_CS=0.997d0
    Vckm_CB=0.0422d0
    Vckm_UB=0.0394d0
    
    !-------------------Parameters for uPDFs evaluation
    number_of_uPDFs=1
    if(allocated(enumeration_of_uPDFs)) deallocate(enumeration_of_uPDFs)
    if(allocated(replicas_of_uPDFs)) deallocate(replicas_of_uPDFs)
    if(allocated(sets_of_uPDFs)) deallocate(sets_of_uPDFs)
    
    allocate(enumeration_of_uPDFs(1:1))
    allocate(replicas_of_uPDFs(1:1))
    allocate(sets_of_uPDFs(1:1))
    enumeration_of_uPDFs=(/1/)
    replicas_of_uPDFs=(/0/)
    select case(order)
      case('LO','LO+')
	sets_of_uPDFs=(/trim('MMHT2014lo68cl')/)
      case('NLO','NLO+')
	sets_of_uPDFs=(/trim('MMHT2014nlo68cl')/)
      case('NNLO','NNLO+')
	sets_of_uPDFs=(/trim('MMHT2014nnlo68cl')/)
    end select
  
    !-------------------Parameters for uFFs evaluation
    ! by definition we do not initiate any FFs
    number_of_uFFs=0
    if(allocated(enumeration_of_uFFs)) deallocate(enumeration_of_uFFs)
    if(allocated(replicas_of_uFFs)) deallocate(replicas_of_uFFs)
    if(allocated(sets_of_uFFs)) deallocate(sets_of_uFFs)
    allocate(enumeration_of_uFFs(1:1))
    allocate(replicas_of_uFFs(1:1))
    allocate(sets_of_uFFs(1:1))
    enumeration_of_uFFs=(/-1/)
    replicas_of_uFFs=(/0/)
    sets_of_uFFs=(/trim('ABSENT')/)
    
    !-------------------- parameters for TMDR
    include_TMDR=.true.
    TMDR_order=trim(order)
    TMDR_evolutionType=3 !1 = improved D solution ,2 = improved gamma solution, 3 = fixed mu
    TMDR_lambdaLength=2
    TMDR_tolerance=0.0001d0	!tolerance of searching for saddle point, numerical integration etc.
    
   !-------------------- parameters for UTMDPDF
    include_uTMDPDF=.true.
    uTMDPDF_order=trim(order)
    uTMDPDF_makeGrid=.true.
    uTMDPDF_withGluon=.false.
    uTMDPDF_lambdaLength=2
    uTMDPDF_tolerance=0.0001d0	!tolerance (i.e. relative integration tolerance -- in convolution integrals)
    uTMDPDF_maxIteration=10000	!maxIteration for adaptive integration
    uTMDPDF_grid_xMin=0.00001d0
    uTMDPDF_grid_bMax=100d0
    uTMDPDF_grid_SizeX=250
    uTMDPDF_grid_SizeB=750
    uTMDPDF_grid_slope=10d0
    
    !-------------------- parameters for UTMDFF
    include_uTMDFF=.false.!!! we do not initialize TMDFF by definition
    uTMDFF_order=trim(order)
    uTMDFF_makeGrid=.true.
    uTMDFF_withGluon=.false.
    uTMDFF_lambdaLength=0
    uTMDFF_tolerance=0.0001d0	!tolerance (i.e. relative integration tolerance -- in convolution integrals)
    uTMDFF_maxIteration=10000	!maxIteration for adaptive integration
    uTMDFF_grid_xMin=0.05d0
    uTMDFF_grid_bMax=50d0
    uTMDFF_grid_SizeX=250
    uTMDFF_grid_SizeB=400
    uTMDFF_grid_slope=10d0
    
    
    !------------------ parameters for TMDs
    include_TMDs=.true.
    
    !------------------ parameters for TMDF
    include_TMDF=.true.
    TMDF_tolerance=0.0001d0	!tolerance (i.e. relative integration tolerance)
    TMDF_OGATAh=0.001d0		!Ogata quadrature integration step 
    
    !------------------ parameters for TMDs-inKT
    include_TMDs_inKT=.false.
    TMDs_inKT_tolerance=0.0001d0	!tolerance (i.e. relative integration tolerance)
    TMDs_inKT_OGATAh=0.001d0		!Ogata quadrature integration step 
    
    !------------------ parameters for TMDX-DY
    include_TMDX_DY=.true.
    TMDX_DY_tolerance=0.001d0	!tolerance (i.e. relative integration tolerance -- in kinematic integrals;)
    TMDX_DY_ptSECTION=4		!default number of sections for pt-bin integration
    TMDX_DY_order=trim(order)
    TMDX_DY_exactX1X2=.true.
    TMDX_DY_numProc=8
    
    !------------------ parameters for TMDX-DY
    include_TMDX_SIDIS=.false.
    TMDX_SIDIS_tolerance=0.001d0	!tolerance (i.e. relative integration tolerance -- in kinematic integrals;)
    TMDX_SIDIS_ptSECTION=4		!default number of sections for pt-bin integration
    TMDX_SIDIS_order=trim(order)
    TMDX_SIDIS_qTcorr=.true.
    TMDX_SIDIS_M1corr=.true.
    TMDX_SIDIS_M2corr=.true.
    TMDX_SIDIS_numProc=8
  
  end subroutine SetupDefault	

  !-------------------------------------------------------CHANGE INITILIZATION PARAMETERS-------------------------------------
  subroutine Set_outputLevel(level,numMessages)
  integer::level
  integer,optional::numMessages
  outputLevel=level
  if(present(numMessages)) messageTrigger=numMessages
  end subroutine Set_outputLevel
  
  subroutine artemide_include(a1,a2,a3,a4,a5,a6,a7,a8,a9,a10)
   character(len=*),optional::a1,a2,a3,a4,a5,a6,a7,a8,a9,a10
   if(present(a1)) call SwitchModule(a1)
   if(present(a2)) call SwitchModule(a2)
   if(present(a3)) call SwitchModule(a3)
   if(present(a4)) call SwitchModule(a4)
   if(present(a5)) call SwitchModule(a5)
   if(present(a6)) call SwitchModule(a6)
   if(present(a7)) call SwitchModule(a7)
   if(present(a8)) call SwitchModule(a8)
   if(present(a9)) call SwitchModule(a9)
   if(present(a10)) call SwitchModule(a10)
  end subroutine artemide_include
  
  subroutine SwitchModule(name)
    character(len=*)::name
    select case(trim(name))
      case('QCDinput')
	if(outputLevel>1) write(*,*) 'artemide_setup: Module QCDinput is included by default'
      case('EWinput')
	include_EWinput=.true.
	if(outputLevel>1) write(*,*) 'artemide_setup: Module EWinput is included'
      case('uTMDPDF')
	include_uTMDPDF=.true.
	if(outputLevel>1) write(*,*) 'artemide_setup: Module uTMDPDF is included'
      case('uTMDFF')
	include_uTMDFF=.true.
	if(outputLevel>1) write(*,*) 'artemide_setup: Module uTMDFF is included'
      case('TMDR')
	include_TMDR=.true.
	if(outputLevel>1) write(*,*) 'artemide_setup: Module TMDR is included'
      case('TMDs')
	include_TMDs=.true.
	if(outputLevel>1) write(*,*) 'artemide_setup: Module TMDs is included'
      case('TMDF')
	include_TMDF=.true.
	if(outputLevel>1) write(*,*) 'artemide_setup: Module TMDF is included'
      case('TMDs_inKT')
	include_TMDs_inKT=.true.
	if(outputLevel>1) write(*,*) 'artemide_setup: Module TMDs_inKT is included'
      case('TMDX_DY')
	include_TMDX_DY=.true.
	if(outputLevel>1) write(*,*) 'artemide_setup: Module TMDX_DY is included'
      case('TMDX_SIDIS')
	include_TMDX_SIDIS=.true.
	if(outputLevel>1) write(*,*) 'artemide_setup: Module TMDX_SIDIS is included'
      end select
  end subroutine SwitchModule
  
  !!! Adds a uPDF to initialization list
  subroutine Set_uPDF(hadron,setName,replica)
  integer,intent(in)::hadron
  character(len=*),intent(in)::setName
  integer,intent(in),optional::replica
  integer::pos,i
  integer,allocatable::enum_old(:),rep_old(:)
  character(len=100),allocatable::sets_old(:)
  
  !!! If there is no hadron so far in the list, we make it new
  if(number_of_uPDFs==0) then
    if(outputLevel>2) write(*,"('artemide_setup: uPDF initialization list is made')")
    number_of_uPDFs=1
    pos=1 !!! by default all allocated with single element, we have to just rewrite it.
    
    
  else if(ANY(enumeration_of_uPDFs.eq.hadron)) then!!!! hadron already in the grid. We ust redefine it
    if(outputLevel>2) write(*,"('artemide_setup: uPDF for hadron',I3,' is redefined')") hadron
    do i=1,number_of_uPDFs
      if(enumeration_of_uPDFs(i)==hadron) exit
    end do
    pos=i
    
    
  else	!!!! this hadron is NEW
    if(outputLevel>2) write(*,"('artemide_setup: uPDF for hadron',I3,' is added')") hadron
    
    !!save old arrays
    allocate(enum_old(1:number_of_uPDFs))
    allocate(rep_old(1:number_of_uPDFs))
    allocate(sets_old(1:number_of_uPDFs))
    do i=1,number_of_uPDFs
      enum_old(i)=enumeration_of_uPDFs(i)
      rep_old(i)=replicas_of_uPDFs(i)
      sets_old(i)=sets_of_uPDFs(i)
    end do
    
    !! reallocating arrays
    deallocate(enumeration_of_uPDFs,replicas_of_uPDFs,sets_of_uPDFs)
    
    number_of_uPDFs=number_of_uPDFs+1
    allocate(enumeration_of_uPDFs(1:number_of_uPDFs))
    allocate(replicas_of_uPDFs(1:number_of_uPDFs))
    allocate(sets_of_uPDFs(1:number_of_uPDFs))
    !! copy information
    do i=1,number_of_uPDFs-1
      enumeration_of_uPDFs(i)=enum_old(i)
      replicas_of_uPDFs(i)=rep_old(i)
      sets_of_uPDFs(i)=sets_old(i)
    end do
    pos=number_of_uPDFs
    
    deallocate(enum_old,sets_old,rep_old)
  end if
  
  enumeration_of_uPDFs(pos)=hadron
  sets_of_uPDFs(pos)=trim(setName)
  if(present(replica)) then
    replicas_of_uPDFs(pos)=replica
  else
    replicas_of_uPDFs(pos)=0
  end if
  
  if(outputLevel>1) write(*,"('artemide_setup: uPDF ',A,' for hadron ',I3,' added to initializaton list')") trim(setName),hadron
  
  end subroutine Set_uPDF
  
  !!! Adds a uFF to initialization list
  subroutine Set_uFF(hadron,setName,replica)
  integer,intent(in)::hadron
  character(len=*),intent(in)::setName
  integer,intent(in),optional::replica
  integer::pos,i
  integer,allocatable::enum_old(:),rep_old(:)
  character(len=100),allocatable::sets_old(:)
  
  !!! If there is no hadron so far in the list, we make it new
  if(number_of_uFFs==0) then
    if(outputLevel>2) write(*,"('artemide_setup: uFF initialization list is made')")
    
    number_of_uFFs=1
    pos=1 !!! by default all allocated with single element, we have to just rewrite it.
    
    
  else if(ANY(enumeration_of_uFFs.eq.hadron)) then!!!! hadron already in the grid. We ust redefine it
    if(outputLevel>2) write(*,"('artemide_setup: uFF for hadron',I3,' is redefined')") hadron
    do i=1,number_of_uFFs
      if(enumeration_of_uFFs(i)==hadron) exit
    end do
    pos=i
    
    
  else	!!!! this hadron is NEW
    if(outputLevel>2) write(*,"('artemide_setup: uFF for hadron',I3,' is added')") hadron
    
    !!save old arrays
    allocate(enum_old(1:number_of_uFFs))
    allocate(rep_old(1:number_of_uFFs))
    allocate(sets_old(1:number_of_uFFs))
    do i=1,number_of_uFFs
      enum_old(i)=enumeration_of_uFFs(i)
      rep_old(i)=replicas_of_uFFs(i)
      sets_old(i)=sets_of_uFFs(i)
    end do
    
    !! reallocating arrays
    deallocate(enumeration_of_uFFs,replicas_of_uFFs,sets_of_uFFs)
    
    number_of_uFFs=number_of_uFFs+1
    allocate(enumeration_of_uFFs(1:number_of_uFFs))
    allocate(replicas_of_uFFs(1:number_of_uFFs))
    allocate(sets_of_uFFs(1:number_of_uFFs))
    !! copy information
    do i=1,number_of_uFFs-1
      enumeration_of_uFFs(i)=enum_old(i)
      replicas_of_uFFs(i)=rep_old(i)
      sets_of_uFFs(i)=sets_old(i)
    end do
    pos=number_of_uFFs
    
    deallocate(enum_old,sets_old,rep_old)
  end if
  
  enumeration_of_uFFs(pos)=hadron
  sets_of_uFFs(pos)=trim(setName)
  if(present(replica)) then
    replicas_of_uFFs(pos)=replica
  else
    replicas_of_uFFs(pos)=0
  end if
  
  if(outputLevel>1) write(*,"('artemide_setup: uFF ',A,' for hadron ',I3,' added to initializaton list')") trim(setName),hadron
  
  end subroutine Set_uFF
  
  subroutine Set_quarkMasses(mC,mB)
  real,optional::mC,mB
  if(present(mC)) mCHARM=mC
  if(present(mB)) mCHARM=mC
  
  if(outputLevel>1) write(*,"('artemide_setup: quark masses reset (mCHARM,mBOTTOM)=(',F6.3,',',F6.3,')')") mCHARM,mBOTTOM
  end subroutine Set_quarkMasses
  
  subroutine Set_EWparameters(alphaInv,massZ,massW,widthZ,widthW,sin2ThetaW,UD,US,UB,CD,CS,CB)
  real*8,optional::alphaInv,massZ,massW,widthZ,widthW,sin2ThetaW,UD,US,UB,CD,CS,CB
   if(present(alphaInv)) alphaQED_MZ=alphaInv
   if(present(massZ)) MZ=massZ
   if(present(massW)) MW=massW
   if(present(widthZ)) GammaZ=widthZ
   if(present(widthW)) GammaW=widthW
   if(present(sin2ThetaW)) sW2=sin2ThetaW
   if(present(UD)) Vckm_UD=UD
   if(present(US)) Vckm_US=US
   if(present(UB)) Vckm_UB=UB
   if(present(CD)) Vckm_CD=CD
   if(present(CS)) Vckm_CS=CS
   if(present(CB)) Vckm_CB=CB
  
  if(outputLevel>1) write(*,"('artemide_setup: EW parameters reset')")
  end subroutine Set_EWparameters
  
  !-------------------------
  subroutine Set_TMDR_order(order)
  character(len=*)::order
  
  if(len(order)<8) then
    TMDR_order=trim(order)
  else
    TMDR_order=trim(order(1:8))
  end if
  
  if(outputLevel>1) write(*,"('artemide_setup: TMDR order is changed to ',A)") TMDR_order
  
  end subroutine Set_TMDR_order
  
  subroutine Set_TMDR_evolutionType(num)
  integer::num
  
  TMDR_evolutionType=num
  
  if(outputLevel>1) write(*,"('artemide_setup: TMDR evolution type is changed to ',I3)") TMDR_evolutionType
  
  end subroutine Set_TMDR_evolutionType
  
  subroutine Set_TMDR_lengthNParray(num)
  integer::num
  
  TMDR_lambdaLength=num
  
  if(outputLevel>1) write(*,"('artemide_setup: TMDR length of NP array is set to ',I3)") TMDR_lambdaLength
  
  end subroutine Set_TMDR_lengthNParray
  
  !-------------------------
  subroutine Set_uTMDPDF(hadron,setName,replica)
  integer,intent(in)::hadron
  character(len=*),intent(in)::setName
  integer,intent(in),optional::replica
    
    if(.not.include_uTMDPDF) then
      include_uTMDPDF=.true.
      if(outputLevel>1) write(*,"('artemide_setup: uTMDPDF module included into initialisaion list')") 
    end if
    
    if(outputLevel>1) write(*,"('artemide_setup: uTMDPDF ',A,' for hadron ',I3,' added to grid list')") trim(setName),hadron
    call Set_uPDF(hadron,setName,replica)
  
  end subroutine Set_uTMDPDF
  
  subroutine Set_uTMDPDF_order(order)
  character(len=*)::order
  
  if(len(order)<8) then
    uTMDPDF_order=trim(order)
  else
    uTMDPDF_order=trim(order(1:8))
  end if
  
  if(outputLevel>1) write(*,"('artemide_setup: uTMDPDF order is changed to ',A)") uTMDPDF_order
  
  end subroutine Set_uTMDPDF_order
  
  subroutine Set_uTMDPDF_gridEvaluation(prepareGrid,includeGluon)
  logical::prepareGrid
  logical,optional::includeGluon
  
  uTMDPDF_makeGrid=prepareGrid
  if(present(includeGluon)) uTMDPDF_withGluon=includeGluon
  
  if(outputLevel>1) write(*,"('artemide_setup: uTMDPDF grid evaluation is changed to (',L2,',',L2,')')") prepareGrid,includeGluon
  
  end subroutine Set_uTMDPDF_gridEvaluation
  
  subroutine Set_uTMDPDF_lengthNParray(num)
  integer::num
  
  uTMDPDF_lambdaLength=num
  
  if(outputLevel>1) write(*,"('artemide_setup: uTMDPDF length of NP array is set to ',I3)") uTMDPDF_lambdaLength
  
  end subroutine Set_uTMDPDF_lengthNParray
  
  !-------------------------
  
  subroutine Set_uTMDFF(hadron,setName,replica)
  
  integer,intent(in)::hadron
  character(len=*),intent(in)::setName
  integer,intent(in),optional::replica
    
    if(.not.include_uTMDFF) then
      include_uTMDFF=.true.
      if(outputLevel>1) write(*,"('artemide_setup: uTMDFF module included into initialisaion list')") 
    end if
    
    if(outputLevel>1) write(*,"('artemide_setup: uTMDFF ',A,' for hadron ',I3,' added to grid list')") trim(setName),hadron
    call Set_uFF(hadron,setName,replica)
  
  end subroutine Set_uTMDFF
  
  subroutine Set_uTMDFF_order(order)
  character(len=*)::order
  
  if(len(order)<8) then
    uTMDFF_order=trim(order)
  else
    uTMDFF_order=trim(order(1:8))
  end if
  
  if(outputLevel>1) write(*,"('artemide_setup: uTMDFF order is changed to ',A)") uTMDFF_order
  
  end subroutine Set_uTMDFF_order
  
  subroutine Set_uTMDFF_gridEvaluation(prepareGrid,includeGluon)
  logical::prepareGrid
  logical,optional::includeGluon
  
  uTMDFF_makeGrid=prepareGrid
  if(present(includeGluon)) uTMDFF_withGluon=includeGluon
  
  if(outputLevel>1) write(*,"('artemide_setup: uTMDFF grid evaluation is changed to (',L2,',',L2,')')") prepareGrid,includeGluon
  
  end subroutine Set_uTMDFF_gridEvaluation
  
  subroutine Set_uTMDFF_lengthNParray(num)
  integer::num
  
  uTMDFF_lambdaLength=num
  
  if(outputLevel>1) write(*,"('artemide_setup: uTMDFF length of NP array is set to ',I3)") uTMDFF_lambdaLength
  
  end subroutine Set_uTMDFF_lengthNParray
  
  
  
  !-------------------------------------------------------READ AND WRITE CONSTANTS FILE---------------------------------------
  subroutine CreateConstantsFile(file,prefix)
  character(len=*)::file
  character(len=*),optional::prefix
  character(len=300)::path
  
  integer::values(1:8),i
  
  if(outputLevel>2) write(*,*) '-----------------------------------'
  if(outputLevel>1) write(*,*) 'aTMDe_setup: Creating setup file ...'
  if(present(prefix)) then
    path=trim(adjustl(prefix))//trim(adjustr(file))
  else
    path=trim(adjustr(file))
  end if
  if(outputLevel>2) write(*,*) '   path for constants-file:',path
  call date_and_time(values=values)
  OPEN(UNIT=51, FILE=path, ACTION="write", STATUS="replace")
    write(51,"('# Setup file for artemide')")
    write(51,"('# File generated by aTMDe_setup ',A)") version
    write(51,"('# File created at ',I3,':',I2,' (',I3,'/',I2,'/',I4,')')")values(5),values(6),values(3),values(2),values(1)
    write(51,"('# ')")
    write(51,"('# FORMAT: Lines started from * determine input parameter (*????), following line determine parameter value')")
    write(51,"(' ')")
    write(51,"(' ')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('# ----                           GLOBAL PARAMETERS                      -----')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('*0   :')")
    write(51,"(' ')")
    write(51,"('*A   : ---- Global I/O ----')")
    write(51,"('*p1  : Input file Version (used to check version compatibility)')")
    write(51,*) inputVer
    write(51,"('*p2  : Output level (0 = no output except critical, 1 = + WARNINGS, 2 = + info,  3 = + details)')")
    write(51,*) outputLevel
    write(51,"('*p3  : Message trigger (the number of continious non-critical warning before stop showing them)')")
    write(51,*) messageTrigger
    write(51,"(' ')")
    write(51,"('*B   : ---- Universal physic parameters ----')")
    write(51,"('*p1  : Unit transformation constant (hc)^2 [GeV->mbarn]')")
    write(51,*) hc2
    write(51,"(' ')")
    
    
    write(51,"(' ')")
    write(51,"(' ')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('# ----                           PARAMETERS OF QCDinput                 -----')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('*1   :')")
    write(51,"(' ')")
    write(51,"('*A   : ---- quarks threashold masses ----')")
    write(51,"('*p1  : mass of charm-quark')")
    write(51,*) mCHARM
    write(51,"('*p2  : mass of bottom-quark')")
    write(51,*) mBOTTOM
    write(51,"(' ')")
    write(51,"('*B   : ---- uPDF sets----')")
    write(51,"('*p1  : total number of PDFs to initialize (0= initialization is skipped)')")
    write(51,*) number_of_uPDFs
    write(51,"('*p2  : reference number for hadrons')")
    write(51,*) enumeration_of_uPDFs
    write(51,"('*p3  : LHAPDF set names for hadrons (line-by-line corresponding to reference number')")
    do i=1,number_of_uPDFs
      write(51,*) trim(sets_of_uPDFs(i))
    end do
    write(51,"('*p4  : list of initialization replicas')")
    write(51,*) replicas_of_uPDFs
    
    write(51,"(' ')")
    write(51,"('*C   : ---- uFF sets----')")
    write(51,"('*p1  : total number of FFs to initialize (0= initialization is skipped)')")
    write(51,*) number_of_uFFs
    write(51,"('*p2  : reference number for hadrons')")
    write(51,*) enumeration_of_uFFs
    write(51,"('*p3  : LHAPDF set names for hadrons (line-by-line corresponding to reference number')")
    do i=1,number_of_uFFs
      write(51,*) trim(sets_of_uFFs(i))
    end do
    write(51,"('*p4  : list of initialization replicas')")
    write(51,*) replicas_of_uFFs
    
    
    write(51,"(' ')")
    write(51,"(' ')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('# ----                           PARAMETERS OF EWinput                  -----')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('*2   :')")
    write(51,"('*p1  : initialize EWinput module')")
    write(51,*) include_EWinput
    write(51,"(' ')")
    write(51,"('*A   : ---- EW coupling parameters ----')")
    write(51,"('*p1  : value of (alphaQED)^{-1} at MZ')")
    write(51,*) alphaQED_MZ
    write(51,"('*p2  : value of sin^2(theta_W)')")
    write(51,*) sW2
    write(51,"('*p3  : values of CKM matrix (1-line:UD, US, UB, 2-line: CD, CS, CB)')")
    write(51,*) Vckm_UD,Vckm_US,Vckm_UB
    write(51,*) Vckm_CD,Vckm_CS,Vckm_CB    
    write(51,"('*B   : ---- Z-boson ----')")
    write(51,"('*p1  : mass of Z-boson [GeV]')")
    write(51,*) mZ
    write(51,"('*p2  : width of Z-boson [GeV]')")
    write(51,*) GammaZ
    write(51,"('*C   : ---- W-boson ----')")
    write(51,"('*p1  : mass of W-boson [GeV]')")
    write(51,*) mW
    write(51,"('*p2  : width of W-boson [GeV]')")
    write(51,*) GammaW
    
    write(51,"(' ')")
    write(51,"(' ')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('# ----                            PARAMETERS OF TMDR                    -----')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('*3   :')")
    write(51,"('*p1  : initialize TMDR module')")
    write(51,*) include_TMDR
    write(51,"(' ')")
    write(51,"('*A   : ---- Main definitions ----')")
    write(51,"('*p1  : Order of evolution')")
    write(51,*) trim(TMDR_order)
    write(51,"('*p2  : Type of evolution solution (1 = improved D solution ,2 = improved gamma solution, 3 = fixed mu)')")
    write(51,*) TMDR_evolutionType
    write(51,"('*B   : ---- Parameters of NP model ----')")
    write(51,"('*p1  : Length of lambdaNP')")
    write(51,*) TMDR_lambdaLength
    write(51,"('*C   : ---- Numerical evaluation parameters ----')")
    write(51,"('*p1  : Tolerance (relative tolerance of convolution integral)')")
    write(51,*) TMDR_tolerance
    
    write(51,"(' ')")
    write(51,"(' ')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('# ----                           PARAMETERS OF uTMDPDF                  -----')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('*4   :')")
    write(51,"('*p1  : initialize uTMDPDF module')")
    write(51,*) include_uTMDPDF
    write(51,"(' ')")
    write(51,"('*A   : ---- Main definitions ----')")
    write(51,"('*p1  : Order of coefficient function')")
    write(51,*) trim(uTMDPDF_order)
    write(51,"('*B   : ---- Parameters of NP model ----')")
    write(51,"('*p1  : Length of lambdaNP')")
    write(51,*) uTMDPDF_lambdaLength
    write(51,"('*C   : ---- Numerical evaluation parameters ----')")
    write(51,"('*p1  : Tolerance (relative tolerance of convolution integral)')")
    write(51,*) uTMDPDF_tolerance
    write(51,"('*p2  : Maximum number of iterations (for adaptive integration))')")
    write(51,*) uTMDPDF_maxIteration
    write(51,"('*D   : ---- Grid preparation options ----')")
    write(51,"('*p1  : Prepare grid')")
    write(51,*) uTMDPDF_makeGrid
    write(51,"('*p2  : Include gluon TMDs into the grid')")
    write(51,*) uTMDPDF_withGluon
    write(51,"('*p3  : total number of PDFs added to the grid (by default it coincides with number of initialized PDFs)')")
    write(51,*) number_of_uPDFs
    write(51,"('*p4  : reference numbers for hadrons (by default it coincides with references for PDFs)')")
    write(51,*) enumeration_of_uPDFs
    write(51,"('*E   : ---- Parameters of grid ----')")
    write(51,"('*p1  : xGrid_Min the minimal value of x in grid (max=1), make sure that it is enough)')")
    write(51,*) uTMDPDF_grid_xMin
    write(51,"('*p2  : the maximum bT in grid (min=0), for larger approximate extrapolation is done')")
    write(51,*) uTMDPDF_grid_bMax
    write(51,"('*p3  : GridSizeX (250 is enough for 10^-8 presicion, for grid up to 10^-5)')")
    write(51,*) uTMDPDF_grid_SizeX
    write(51,"('*p4  : GridSizeB (750 is enough for ~10^-6 presicion, 500 for 10^-5 presicion)')")
    write(51,*) uTMDPDF_grid_SizeB
    write(51,"('*p5  : slope scale of griding at smaller b (better not to change it :) )')")
    write(51,*) uTMDPDF_grid_slope
    
    
    write(51,"(' ')")
    write(51,"(' ')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('# ----                           PARAMETERS OF uTMDFF                   -----')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('*5   :')")
    write(51,"('*p1  : initialize uTMDFF module')")
    write(51,*) include_uTMDFF
    write(51,"(' ')")
    write(51,"('*A   : ---- Main definitions ----')")
    write(51,"('*p1  : Order of coefficient function')")
    write(51,*) trim(uTMDFF_order)
    write(51,"('*B   : ---- Parameters of NP model ----')")
    write(51,"('*p1  : Length of lambdaNP')")
    write(51,*) uTMDFF_lambdaLength
    write(51,"('*C   : ---- Numerical evaluation parameters ----')")
    write(51,"('*p1  : Tolerance (relative tolerance of convolution integral)')")
    write(51,*) uTMDFF_tolerance
    write(51,"('*p2  : Maximum number of iterations (for adaptive integration))')")
    write(51,*) uTMDFF_maxIteration
    write(51,"('*D   : ---- Grid preparation options ----')")
    write(51,"('*p1  : Prepare grid')")
    write(51,*) uTMDFF_makeGrid
    write(51,"('*p2  : Include gluon TMDs into the grid')")
    write(51,*) uTMDFF_withGluon
    write(51,"('*p3  : total number of FFs added to the grid (by default it coincides with number of initialized FFs)')")
    write(51,*) number_of_uFFs
    write(51,"('*p4  : reference numbers for hadrons (by default it coincides with references for FFs)')")
    write(51,*) enumeration_of_uFFs
    write(51,"('*E   : ---- Parameters of grid ----')")
    write(51,"('*p1  : xGrid_Min the minimal value of x in grid (max=1), make sure that it is enough)')")
    write(51,*) uTMDFF_grid_xMin
    write(51,"('*p2  : the maximum bT in grid (min=0), for larger approximate extrapolation is done')")
    write(51,*) uTMDFF_grid_bMax
    write(51,"('*p3  : GridSizeX (250 is enough for 10^-8 presicion, for grid up to 10^-5)')")
    write(51,*) uTMDFF_grid_SizeX
    write(51,"('*p4  : GridSizeB (750 is enough for ~10^-6 presicion, 500 for 10^-5 presicion)')")
    write(51,*) uTMDFF_grid_SizeB
    write(51,"('*p5  : slope scale of griding at smaller b (better not to change it :) )')")
    write(51,*) uTMDFF_grid_slope
    
    
    write(51,"(' ')")
    write(51,"(' ')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('# ----                            PARAMETERS OF TMDs                    -----')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('*6   :')")
    write(51,"('*p1  : initialize TMDs module')")
    write(51,*) include_TMDs
    
    
    write(51,"(' ')")
    write(51,"(' ')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('# ----                            PARAMETERS OF TMDF                    -----')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('*7   :')")
    write(51,"('*p1  : initialize TMDF module')")
    write(51,*) include_TMDF
    write(51,"('*A   : ---- Numerical evaluation parameters ----')")
    write(51,"('*p1  : Tolerance (relative tolerance of summation in OGATA quadrature)')")
    write(51,*) TMDF_tolerance
    write(51,"('*p2  : Ogata quadrature integration step')")
    write(51,*) TMDF_OGATAh
    
    
    write(51,"(' ')")
    write(51,"(' ')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('# ----                          PARAMETERS OF TMDs-inKT                 -----')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('*8   :')")
    write(51,"('*p1  : initialize TMDs-inKT module')")
    write(51,*) include_TMDs_inKT
    write(51,"('*A   : ---- Numerical evaluation parameters ----')")
    write(51,"('*p1  : Tolerance (relative tolerance of summation in OGATA quadrature)')")
    write(51,*) TMDs_inKT_tolerance
    write(51,"('*p2  : Ogata quadrature integration step')")
    write(51,*) TMDs_inKT_OGATAh
    
    
    write(51,"(' ')")
    write(51,"(' ')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('# ----                           PARAMETERS OF TMDX-DY                  -----')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('*9   :')")
    write(51,"('*p1  : initialize TMDX-DY module')")
    write(51,*) include_TMDX_DY
    write(51,"('*A   : ---- Main definitions ----')")
    write(51,"('*p1  : Order of coefficient function')")
    write(51,*) trim(TMDX_DY_order)
    write(51,"('*p2  : Use the exact values of x1 and x2 (include qT/Q correction)')")
    write(51,*) TMDX_DY_exactX1X2
    write(51,"('*B   : ---- Numerical evaluation parameters ----')")
    write(51,"('*p1  : Tolerance (relative tolerance for bin-integration routines, except pt-integration))')")
    write(51,*) TMDX_DY_tolerance
    write(51,"('*p2  : Minimal number of sections for pt-integration')")
    write(51,*) TMDX_DY_ptSECTION
    write(51,"('*C   : ---- Parameters for parrallel evaluation (used only if compiled with openMP) ----')")
    write(51,"('*p1  : Maximum number of processors to use')")
    write(51,*) TMDX_DY_numProc
    
    
    write(51,"(' ')")
    write(51,"(' ')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('# ----                          PARAMETERS OF TMDX-SIDIS                -----')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('*10  :')")
    write(51,"('*p1  : initialize TMDX-SIDIS module')")
    write(51,*) include_TMDX_SIDIS
    write(51,"('*A   : ---- Main definitions ----')")
    write(51,"('*p1  : Order of coefficient function')")
    write(51,*) trim(TMDX_SIDIS_order)
    write(51,"('*p2  : Use transverse momentum corrections in kinematics')")
    write(51,*) TMDX_SIDIS_qTcorr
    write(51,"('*p3  : Use target mass corrections in kinematics')")
    write(51,*) TMDX_SIDIS_M1corr
    write(51,"('*p4  : Use product mass corrections in kinematics')")
    write(51,*) TMDX_SIDIS_M2corr
    write(51,"('*B   : ---- Numerical evaluation parameters ----')")
    write(51,"('*p1  : Tolerance (relative tolerance for bin-integration routines, except pt-integration))')")
    write(51,*) TMDX_SIDIS_tolerance
    write(51,"('*p2  : Minimal number of sections for pt-integration')")
    write(51,*) TMDX_SIDIS_ptSECTION
    write(51,"('*C   : ---- Parameters for parrallel evaluation (used only if compiled with openMP) ----')")
    write(51,"('*p1  : Maximum number of processors to use')")
    write(51,*) TMDX_SIDIS_numProc
    
  CLOSE (51, STATUS='KEEP') 
    
    if(outputLevel>1) write(*,*) 'aTMDe_setup: Constans file is written.'
  
  end subroutine CreateConstantsFile
  
 !!! move CURRET in streem to the next line that starts from pos (5 char)
 subroutine MoveTO(streem,pos)
 integer,intent(in)::streem
 character(len=5)::pos
 character(len=300)::line
    !write(*,*) 'SEARCH FOR ====> ',pos
    do
    read(streem,'(A)') line    
    if(line(1:5)==pos) exit
    end do
 end subroutine MoveTO
  
  subroutine ReadConstantsFile(file,prefix)
  character(len=*)::file
  character(len=*),optional::prefix
  character(len=300)::path
  !!!! this is version of input file. it is read first and then result is compared with the current ID
  !!!! It suppose to make compatibility betwen versions
  integer::FILEversion
  integer::i
  
  if(outputLevel>2) write(*,*) '-----------------------------------'
  if(outputLevel>1) write(*,*) 'aTMDe_setup: Reading setup file ...'
  if(present(prefix)) then
    path=trim(adjustl(prefix))//trim(adjustr(file))
  else
    path=trim(adjustr(file))
  end if
  if(outputLevel>2) write(*,*) 'path:',path
  
  OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")
    !# ----                           GLOBAL PARAMETERS                      -----
    call MoveTO(51,'*0   ')
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) FILEversion
    
    if(FILEversion>inputVer) then
      if(outputLevel>0) write(*,*) 'aTMDe_setup: Version of input file newer then the current version of code'
      if(outputLevel>0) write(*,*) 'aTMDe_setup: UPDATE ARTEMIDE!'
      if(outputLevel>0) write(*,*) 'aTMDe_setup: Attempt to load...'
    else if(FILEversion<inputVer) then
      if(outputLevel>0) write(*,*) 'aTMDe_setup: Version of input file older then the current version of code'
      if(outputLevel>0) write(*,*) 'aTMDe_setup: Attempt to load...'
    end if
    
    call MoveTO(51,'*p2  ')
    read(51,*) outputLevel
    call MoveTO(51,'*p3  ')
    read(51,*) messageTrigger
    
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) hc2
    
    !# ----                           PARAMETERS OF QCDinput                 -----
    call MoveTO(51,'*1   ')
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) mCHARM
    call MoveTO(51,'*p2  ')
    read(51,*) mBOTTOM
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    
    read(51,*) number_of_uPDFs
    if(allocated(enumeration_of_uPDFs)) deallocate(enumeration_of_uPDFs)
    if(allocated(replicas_of_uPDFs)) deallocate(replicas_of_uPDFs)
    if(allocated(sets_of_uPDFs)) deallocate(sets_of_uPDFs)
    allocate(enumeration_of_uPDFs(1:number_of_uPDFs))
    allocate(sets_of_uPDFs(1:number_of_uPDFs))
    allocate(replicas_of_uPDFs(1:number_of_uPDFs))    
    
    call MoveTO(51,'*p2  ')
    read(51,*) enumeration_of_uPDFs
    call MoveTO(51,'*p3  ')
    do i=1,number_of_uPDFs
      read(51,*) sets_of_uPDFs(i)
    end do
    call MoveTO(51,'*p4  ')
    read(51,*) replicas_of_uPDFs
    
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) number_of_uFFs
    
    if(allocated(enumeration_of_uFFs)) deallocate(enumeration_of_uFFs)
    if(allocated(replicas_of_uFFs)) deallocate(replicas_of_uFFs)
    if(allocated(sets_of_uFFs)) deallocate(sets_of_uFFs)
    allocate(enumeration_of_uFFs(1:number_of_uFFs))
    allocate(sets_of_uFFs(1:number_of_uFFs))
    allocate(replicas_of_uFFs(1:number_of_uFFs))    
    
    call MoveTO(51,'*p2  ')
    read(51,*) enumeration_of_uFFs
    call MoveTO(51,'*p3  ')
    do i=1,number_of_uFFs
      read(51,*) sets_of_uFFs(i)
    end do
    call MoveTO(51,'*p4  ')
    read(51,*) replicas_of_uFFs
    
    
    !# ----                           PARAMETERS OF EWinput                  -----
    call MoveTO(51,'*2   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_EWinput
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) alphaQED_MZ
    call MoveTO(51,'*p2  ')
    read(51,*) sW2
    call MoveTO(51,'*p3  ')
    read(51,*) Vckm_UD,Vckm_US,Vckm_UB
    read(51,*) Vckm_CD,Vckm_CS,Vckm_CB
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) mZ
    call MoveTO(51,'*p2  ')
    read(51,*) GammaZ
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) mW
    call MoveTO(51,'*p2  ')
    read(51,*) GammaW
    
    
    !# ----                            PARAMETERS OF TMDR                    -----
    call MoveTO(51,'*3   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_TMDR
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDR_order
    call MoveTO(51,'*p2  ')
    read(51,*) TMDR_evolutionType
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDR_lambdaLength
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDR_tolerance
    
    
    !# ----                           PARAMETERS OF uTMDPDF                  -----
    call MoveTO(51,'*4   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_uTMDPDF
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) uTMDPDF_order
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) uTMDPDF_lambdaLength
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) uTMDPDF_tolerance
    call MoveTO(51,'*p2  ')
    read(51,*) uTMDPDF_maxIteration
    call MoveTO(51,'*D   ')
    call MoveTO(51,'*p1  ')
    read(51,*) uTMDPDF_makeGrid
    call MoveTO(51,'*p2  ')
    read(51,*) uTMDPDF_withGluon
    call MoveTO(51,'*p3  ')
    read(51,*) i
    if(i/=number_of_uPDFs) then
      if(outputLevel>0) write(*,*) 'ESSENTIAL INCONSITENCY: number of uPDFs unequal uTMDPDFs'
      if(outputLevel>0) write(*,*) '                        information on number of uTMDPDFs ignored'
    end if
    call MoveTO(51,'*E   ')
    call MoveTO(51,'*p1  ')
    read(51,*) uTMDPDF_grid_xMin
    call MoveTO(51,'*p2  ')
    read(51,*) uTMDPDF_grid_bMax
    call MoveTO(51,'*p3  ')
    read(51,*) uTMDPDF_grid_SizeX
    call MoveTO(51,'*p4  ')
    read(51,*) uTMDPDF_grid_SizeB
    call MoveTO(51,'*p5  ')
    read(51,*) uTMDPDF_grid_slope
    
    
    !# ----                           PARAMETERS OF uTMDFF                   -----
    call MoveTO(51,'*5   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_uTMDFF
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) uTMDFF_order
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) uTMDFF_lambdaLength
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) uTMDFF_tolerance
    call MoveTO(51,'*p2  ')
    read(51,*) uTMDFF_maxIteration
    call MoveTO(51,'*D   ')
    call MoveTO(51,'*p1  ')
    read(51,*) uTMDFF_makeGrid
    call MoveTO(51,'*p2  ')
    read(51,*) uTMDFF_withGluon
    call MoveTO(51,'*p3  ')
    read(51,*) i
    if(i/=number_of_uFFs) then
      if(outputLevel>0) write(*,*) 'ESSENTIAL INCONSITENCY: number of uFFs unequal uTMDFFs'
      if(outputLevel>0) write(*,*) '                        information on number of uTMDFFs ignored'
    end if
    call MoveTO(51,'*E   ')
    call MoveTO(51,'*p1  ')
    read(51,*) uTMDFF_grid_xMin
    call MoveTO(51,'*p2  ')
    read(51,*) uTMDFF_grid_bMax
    call MoveTO(51,'*p3  ')
    read(51,*) uTMDFF_grid_SizeX
    call MoveTO(51,'*p4  ')
    read(51,*) uTMDFF_grid_SizeB
    call MoveTO(51,'*p5  ')
    read(51,*) uTMDFF_grid_slope
    
    
    !# ----                            PARAMETERS OF TMDs                    -----
    call MoveTO(51,'*6   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_TMDs
    
    
    !# ----                            PARAMETERS OF TMDF                    -----
    call MoveTO(51,'*7   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_TMDF
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDF_tolerance
    call MoveTO(51,'*p2  ')
    read(51,*) TMDF_OGATAh
    
    
    !# ----                          PARAMETERS OF TMDs-inKT                 -----
    call MoveTO(51,'*8   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_TMDs_inKT
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDs_inKT_tolerance
    call MoveTO(51,'*p2  ')
    read(51,*) TMDs_inKT_OGATAh
    
    
    !# ----                           PARAMETERS OF TMDX-DY                  -----
    call MoveTO(51,'*9   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_TMDX_DY
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDX_DY_order
    call MoveTO(51,'*p2  ')
    read(51,*) TMDX_DY_exactX1X2
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDX_DY_tolerance
    call MoveTO(51,'*p2  ')
    read(51,*) TMDX_DY_ptSECTION
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDX_DY_numProc
    
    
    !# ----                          PARAMETERS OF TMDX-SIDIS                -----
    call MoveTO(51,'*10  ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_TMDX_SIDIS
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDX_SIDIS_order
    call MoveTO(51,'*p2  ')
    read(51,*) TMDX_SIDIS_qTcorr
    call MoveTO(51,'*p3  ')
    read(51,*) TMDX_SIDIS_M1corr
    call MoveTO(51,'*p4  ')
    read(51,*) TMDX_SIDIS_M2corr
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDX_SIDIS_tolerance
    call MoveTO(51,'*p2  ')
    read(51,*) TMDX_SIDIS_ptSECTION
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDX_SIDIS_numProc
    
  CLOSE (51, STATUS='KEEP') 
    
    if(outputLevel>1) write(*,*) 'aTMDe_setup: constants-file readed loaded sucessfully.'
  
  end subroutine ReadConstantsFile
  
  
end module aTMDe_Setup