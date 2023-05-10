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
use aTMDe_Numerics
use IO_functions
implicit none

private

character (len=5),parameter :: version="v2.06"
character (len=11),parameter :: moduleName="aTMDe-setup"
!! actual version of input file
integer,parameter::inputVer=21

!detalization of output: 0 = no output except critical, 1 = + WARNINGS, 2 = + states of initialization,sets,etc, 3 = + details
integer::outputLevel
integer::messageTrigger

!! if .true. uses the given NP arrays to initialize from aTMDe_control
logical::initialize_NParrays

!-----------------------physic parameters
logical::include_EWinput
!mass, and width parameters
real(dp)::mCHARM,mBOTTOM,mTOP,mZ,mW,GammaZ,GammaW,mH,GammaH,vevH
real(dp)::mELECTRON,mMUON,mTAU
!other physical parameters
real(dp)::hc2,alphaQED_MZ,sW2,alphaQED_MTAU
!CKM matrix
real(dp)::Vckm_UD,Vckm_US,Vckm_CD,Vckm_CS,Vckm_CB,Vckm_UB

!--------------------- uPDF parameters
integer::number_of_uPDFs
integer,allocatable::enumeration_of_uPDFs(:),replicas_of_uPDFs(:)
character(len=100),allocatable::sets_of_uPDFs(:)

!--------------------- uFF parameters
integer::number_of_uFFs
integer,allocatable::enumeration_of_uFFs(:),replicas_of_uFFs(:)
character(len=100),allocatable::sets_of_uFFs(:)

!--------------------- lpPDF parameters
integer::number_of_lpPDFs
integer,allocatable::enumeration_of_lpPDFs(:),replicas_of_lpPDFs(:)
character(len=100),allocatable::sets_of_lpPDFs(:)

!--------------------- hPDF parameters
integer::number_of_hPDFs
integer,allocatable::enumeration_of_hPDFs(:),replicas_of_hPDFs(:)
character(len=100),allocatable::sets_of_hPDFs(:)

!-------------------- TMDR options
logical::include_TMDR
character*8::TMDR_order
integer::TMDR_evolutionType,TMDR_lambdaLength
real(dp)::TMDR_tolerance
real(dp),allocatable::TMDR_lambdaNP_init(:)

!-------------------- uTMDPDF parameters
logical::include_uTMDPDF
integer::uTMDPDF_lambdaLength
real(dp)::uTMDPDF_tolerance
integer::uTMDPDF_maxIteration
character*8::uTMDPDF_order
logical::uTMDPDF_makeGrid,uTMDPDF_withGluon
real(dp)::uTMDPDF_grid_bMax,uTMDPDF_grid_xMin,uTMDPDF_grid_slope
integer::uTMDPDF_grid_SizeX,uTMDPDF_grid_SizeB
real(dp),allocatable::uTMDPDF_lambdaNP_init(:)
logical::uTMDPDF_IsComposite

!-------------------- uTMDFF parameters
logical::include_uTMDFF
integer::uTMDFF_lambdaLength
real(dp)::uTMDFF_tolerance
integer::uTMDFF_maxIteration
character*8::uTMDFF_order
logical::uTMDFF_makeGrid,uTMDFF_withGluon
real(dp)::uTMDFF_grid_bMax,uTMDFF_grid_xMin,uTMDFF_grid_slope
integer::uTMDFF_grid_SizeX,uTMDFF_grid_SizeB
real(dp),allocatable::uTMDFF_lambdaNP_init(:)
logical::uTMDFF_IsComposite

!-------------------- lpTMDPDF parameters
logical::include_lpTMDPDF
integer::lpTMDPDF_lambdaLength
real(dp)::lpTMDPDF_tolerance
integer::lpTMDPDF_maxIteration
character*8::lpTMDPDF_order
logical::lpTMDPDF_makeGrid,lpTMDPDF_withGluon
real(dp)::lpTMDPDF_grid_bMax,lpTMDPDF_grid_xMin,lpTMDPDF_grid_slope
integer::lpTMDPDF_grid_SizeX,lpTMDPDF_grid_SizeB
real(dp),allocatable::lpTMDPDF_lambdaNP_init(:)
logical::lpTMDPDF_IsComposite

!-------------------- SiversTMDPDF parameters
logical::include_SiversTMDPDF
integer::SiversTMDPDF_lambdaLength
real(dp)::SiversTMDPDF_tolerance
integer::SiversTMDPDF_maxIteration
character*8::SiversTMDPDF_order
logical::SiversTMDPDF_makeGrid,SiversTMDPDF_withGluon
real(dp)::SiversTMDPDF_grid_bMax,SiversTMDPDF_grid_xMin,SiversTMDPDF_grid_slope
integer::SiversTMDPDF_grid_SizeX,SiversTMDPDF_grid_SizeB
real(dp),allocatable::SiversTMDPDF_lambdaNP_init(:)
logical::SiversTMDPDF_IsComposite
integer::number_of_SiversTMDPDFs
integer,allocatable::enumeration_of_SiversTMDPDFs(:)

!-------------------- wgtTMDPDF parameters
logical::include_wgtTMDPDF
integer::wgtTMDPDF_lambdaLength
real(dp)::wgtTMDPDF_tolerance
integer::wgtTMDPDF_maxIteration
character*8::wgtTMDPDF_order
logical::wgtTMDPDF_makeGrid,wgtTMDPDF_withGluon
real(dp)::wgtTMDPDF_grid_bMax,wgtTMDPDF_grid_xMin,wgtTMDPDF_grid_slope
integer::wgtTMDPDF_grid_SizeX,wgtTMDPDF_grid_SizeB
real(dp),allocatable::wgtTMDPDF_lambdaNP_init(:)
logical::wgtTMDPDF_IsComposite
integer::number_of_wgtTMDPDFs
integer,allocatable::enumeration_of_wgtTMDPDFs(:)

!-------------------- TMDs parameters
logical::include_TMDs

!-------------------- TMDF parameters
logical::include_TMDF
real(dp)::TMDF_OGATAh,TMDF_tolerance
real(dp)::TMDF_mass

!-------------------- TMDs-inKT parameters
logical::include_TMDs_inKT
real(dp)::TMDs_inKT_OGATAh,TMDs_inKT_tolerance

!-------------------- TMDX-DY parameters
logical::include_TMDX_DY
character*8::TMDX_DY_order
real(dp)::TMDX_DY_tolerance
integer::TMDX_DY_ptSECTION
logical::TMDX_DY_exactX1X2,TMDX_DY_piResum,TMDX_DY_exactScale
integer::TMDX_DY_numProc
real::TMDX_DY_maxQbinSize

!-------------------- TMDX-SIDIS parameters
logical::include_TMDX_SIDIS
character*8::TMDX_SIDIS_order
real(dp)::TMDX_SIDIS_tolerance
integer::TMDX_SIDIS_ptSECTION
logical::TMDX_SIDIS_qTcorr,TMDX_SIDIS_M1corr,TMDX_SIDIS_M2corr,TMDX_SIDIS_qTinX1Z1corr,TMDX_SIDIS_exactScale
integer::TMDX_SIDIS_numProc
real(dp)::TMDX_SIDIS_toleranceZ,TMDX_SIDIS_toleranceX
character(len=4)::TMDX_SIDIS_methodZ,TMDX_SIDIS_methodX

!---------------------------------------------------
public::artemide_Setup_Default,artemide_Setup_fromFile,artemide_include,CreateConstantsFile,ReadConstantsFile,CheckConstantsFile
public::Set_uPDF,Set_uFF,Set_lpPDF,Set_quarkMasses,Set_EWparameters
public::Set_TMDR_order,Set_TMDR_evolutionType,Set_TMDR_lengthNParray
public::Set_uTMDPDF,Set_uTMDPDF_order,Set_uTMDPDF_gridEvaluation,Set_uTMDPDF_lengthNParray
public::Set_uTMDFF,Set_uTMDFF_order,Set_uTMDFF_gridEvaluation,Set_uTMDFF_lengthNParray
public::Set_lpTMDPDF,Set_lpTMDPDF_order,Set_lpTMDPDF_gridEvaluation,Set_lpTMDPDF_lengthNParray
!public::Set_hTMDPDF,Set_hTMDPDF_order,Set_hTMDPDF_gridEvaluation,Set_hTMDPDF_lengthNParray
contains

INCLUDE 'Code/aTMDe_setup/const-modification.f90'
  
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

    initialize_NParrays=.false.

    include_EWinput=.true.
    !-----------------------physic parameters
    mCHARM=1.2700d0	!threashold mass for charm quark
    mBOTTOM=4.180d0	!threashold mass for bottom quark
    mTOP=172.90d0	!threashold mass for top quark

    mZ=91.1876d0	!pole mass for Z-boson
    GammaZ=2.4952d0	!width of Z-boson

    mW=80.379d0		!pole mass for W-boson
    GammaW=2.085d0	!width of W-boson

    mH=125.1d0		!pole mass for Higgs-boson
    GammaH=0.0042d0	!width of Higgs-boson [GeV]
    vevH=246.d0		!Higgs Vev

    mELECTRON=0.0005109989461d0	!mass of electron[GeV]
    mMUON=0.1056583745d0	!mass of muon[GeV]
    mTAU=1.77686d0		!mass of tau

    hc2=0.389379338d0	!transformation constant (hc)^2   GeV->mbarn

    alphaQED_MZ=127.955d0	!inverse alpha_QED at Z-boson mass
    alphaQED_MTAU=133.476d0	!inverse alpha_QED at TAU-lepton mass

    sW2=0.23122d0	!sin^2 theta_Weinberg

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

    !-------------------Parameters for lpPDFs evaluation
    ! by definition we do not initiate any lpPDF
    number_of_lpPDFs=0
    if(allocated(enumeration_of_lpPDFs)) deallocate(enumeration_of_lpPDFs)
    if(allocated(replicas_of_lpPDFs)) deallocate(replicas_of_lpPDFs)
    if(allocated(sets_of_lpPDFs)) deallocate(sets_of_lpPDFs)
    allocate(enumeration_of_lpPDFs(1:1))
    allocate(replicas_of_lpPDFs(1:1))
    allocate(sets_of_lpPDFs(1:1))
    enumeration_of_lpPDFs=(/-1/)
    replicas_of_lpPDFs=(/0/)
    sets_of_lpPDFs=(/trim('ABSENT')/)
    
    !-------------------Parameters for hPDFs evaluation
    ! by definition we do not initiate any hPDF
    number_of_hPDFs=0
    if(allocated(enumeration_of_hPDFs)) deallocate(enumeration_of_hPDFs)
    if(allocated(replicas_of_hPDFs)) deallocate(replicas_of_hPDFs)
    if(allocated(sets_of_hPDFs)) deallocate(sets_of_hPDFs)
    allocate(enumeration_of_hPDFs(1:1))
    allocate(replicas_of_hPDFs(1:1))
    allocate(sets_of_hPDFs(1:1))
    enumeration_of_hPDFs=(/-1/)
    replicas_of_hPDFs=(/0/)
    sets_of_hPDFs=(/trim('ABSENT')/)

    !-------------------- parameters for TMDR
    include_TMDR=.true.
    TMDR_order=trim(order)
    TMDR_evolutionType=3 !1 = improved D solution ,2 = improved gamma solution, 3 = fixed mu
    TMDR_lambdaLength=2
    call ReNewInitializationArray(TMDR_lambdaNP_init,TMDR_lambdaLength)! initialization values of parameters
    TMDR_tolerance=0.0001d0	!tolerance of searching for saddle point, numerical integration etc.

    !-------------------- parameters for UTMDPDF
    include_uTMDPDF=.true.
    uTMDPDF_order=trim(order)
    uTMDPDF_IsComposite=.false.
    uTMDPDF_makeGrid=.true.
    uTMDPDF_withGluon=.false.
    uTMDPDF_lambdaLength=2
    call ReNewInitializationArray(uTMDPDF_lambdaNP_init,uTMDPDF_lambdaLength)! initialization values of parameters
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
    uTMDFF_IsComposite=.false.
    uTMDFF_makeGrid=.true.
    uTMDFF_withGluon=.false.
    uTMDFF_lambdaLength=1
    call ReNewInitializationArray(uTMDFF_lambdaNP_init,uTMDFF_lambdaLength)! initialization values of parameters
    uTMDFF_tolerance=0.0001d0	!tolerance (i.e. relative integration tolerance -- in convolution integrals)
    uTMDFF_maxIteration=10000	!maxIteration for adaptive integration
    uTMDFF_grid_xMin=0.05d0
    uTMDFF_grid_bMax=50d0
    uTMDFF_grid_SizeX=250
    uTMDFF_grid_SizeB=400
    uTMDFF_grid_slope=10d0

    !-------------------- parameters for lpTMDPDF
    include_lpTMDPDF=.false. !!! we do not initialize lpTMDPDF by definition
    lpTMDPDF_order=trim(order)
    lpTMDPDF_IsComposite=.false.
    lpTMDPDF_makeGrid=.true.
    lpTMDPDF_withGluon=.true.
    lpTMDPDF_lambdaLength=1
    call ReNewInitializationArray(lpTMDPDF_lambdaNP_init,lpTMDPDF_lambdaLength)! initialization values of parameters
    lpTMDPDF_tolerance=0.0001d0	!tolerance (i.e. relative integration tolerance -- in convolution integrals)
    lpTMDPDF_maxIteration=10000	!maxIteration for adaptive integration
    lpTMDPDF_grid_xMin=0.00001d0
    lpTMDPDF_grid_bMax=100d0
    lpTMDPDF_grid_SizeX=250
    lpTMDPDF_grid_SizeB=750
    lpTMDPDF_grid_slope=10d0

    !-------------------- parameters for SiversTMDPDF
    include_SiversTMDPDF=.false. !!! we do not initialize SiversTMDPDF by definition
    SiversTMDPDF_order=trim('LO') !!! by definition Sivers is tree-order
    SiversTMDPDF_IsComposite=.false.
    SiversTMDPDF_makeGrid=.false.   !!! no need to make grid
    SiversTMDPDF_withGluon=.false.
    SiversTMDPDF_lambdaLength=1
    call ReNewInitializationArray(SiversTMDPDF_lambdaNP_init,SiversTMDPDF_lambdaLength)! initialization values of parameters
    SiversTMDPDF_tolerance=0.0001d0	!tolerance (i.e. relative integration tolerance -- in convolution integrals)
    SiversTMDPDF_maxIteration=10000	!maxIteration for adaptive integration
    SiversTMDPDF_grid_xMin=0.00001d0
    SiversTMDPDF_grid_bMax=100d0
    SiversTMDPDF_grid_SizeX=250
    SiversTMDPDF_grid_SizeB=750
    SiversTMDPDF_grid_slope=10d0
    number_of_SiversTMDPDFs=0
    if(allocated(enumeration_of_SiversTMDPDFs)) deallocate(enumeration_of_SiversTMDPDFs)
    allocate(enumeration_of_SiversTMDPDFs(1:1))
    enumeration_of_SiversTMDPDFs=(/-1/)
    
    !-------------------- parameters for wgtTMDPDF
    include_wgtTMDPDF=.false.
    wgtTMDPDF_order=trim('LO')
    wgtTMDPDF_IsComposite=.false.
    wgtTMDPDF_makeGrid=.true.
    wgtTMDPDF_withGluon=.false.
    wgtTMDPDF_lambdaLength=1
    call ReNewInitializationArray(wgtTMDPDF_lambdaNP_init,wgtTMDPDF_lambdaLength)! initialization values of parameters
    wgtTMDPDF_tolerance=0.0001d0	!tolerance (i.e. relative integration tolerance -- in convolution integrals)
    wgtTMDPDF_maxIteration=10000	!maxIteration for adaptive integration
    wgtTMDPDF_grid_xMin=0.001d0
    wgtTMDPDF_grid_bMax=100d0
    wgtTMDPDF_grid_SizeX=250
    wgtTMDPDF_grid_SizeB=500
    wgtTMDPDF_grid_slope=10d0
    number_of_wgtTMDPDFs=0
    if(allocated(enumeration_of_wgtTMDPDFs)) deallocate(enumeration_of_wgtTMDPDFs)
    allocate(enumeration_of_wgtTMDPDFs(1:1))
    enumeration_of_wgtTMDPDFs=(/-1/)
    
    

    !------------------ parameters for TMDs
    include_TMDs=.true.

    !------------------ parameters for TMDF
    include_TMDF=.true.
    TMDF_tolerance=0.0001d0	!tolerance (i.e. relative integration tolerance)
    TMDF_OGATAh=0.001d0		!Ogata quadrature integration step 
    TMDF_mass=0.938272      !mass parameter that is used as reference dimension

    !------------------ parameters for TMDs-inKT
    include_TMDs_inKT=.false.
    TMDs_inKT_tolerance=0.0001d0	!tolerance (i.e. relative integration tolerance)
    TMDs_inKT_OGATAh=0.001d0		!Ogata quadrature integration step 

    !------------------ parameters for TMDX-DY
    include_TMDX_DY=.true.
    TMDX_DY_tolerance=0.001d0	!tolerance (i.e. relative integration tolerance -- in kinematic integrals;)
    TMDX_DY_ptSECTION=4		!default number of sections for pt-bin integration
    TMDX_DY_maxQbinSize=30. !default maximum size of the Q-bin integration
    TMDX_DY_order=trim(order)
    TMDX_DY_exactX1X2=.true.
    TMDX_DY_piResum=.false.
    TMDX_DY_exactScale=.false.
    TMDX_DY_numProc=8

    !------------------ parameters for TMDX-SIDIS
    include_TMDX_SIDIS=.false.
    TMDX_SIDIS_tolerance=0.001d0	!tolerance (i.e. relative integration tolerance -- in kinematic integrals;)
    TMDX_SIDIS_ptSECTION=4		!default number of sections for pt-bin integration
    TMDX_SIDIS_order=trim(order)
    TMDX_SIDIS_qTcorr=.true.
    TMDX_SIDIS_M1corr=.true.
    TMDX_SIDIS_M2corr=.true.
    TMDX_SIDIS_qTinX1Z1corr=.true.
    TMDX_SIDIS_exactScale=.false.
    TMDX_SIDIS_numProc=8
    TMDX_SIDIS_toleranceZ=TMDX_SIDIS_tolerance
    TMDX_SIDIS_methodZ='SA'		!SA=Simpson adaptive, S5=Simpson 5-point
    TMDX_SIDIS_toleranceX=TMDX_SIDIS_tolerance
    TMDX_SIDIS_methodX='SA'		!SA=Simpson adaptive, S5=Simpson 5-point

end subroutine SetupDefault	

!!! theck the file for compatibility with the current version
!!! true= ver.FILE>=current version
!!! false= ver.FILE<current version
function CheckConstantsFile(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=516)::path
    logical::CheckConstantsFile
    integer::FILEversion

    if(present(prefix)) then
        path=trim(adjustl(prefix))//trim(adjustr(file))
    else
        path=trim(adjustr(file))
    end if

    !!!read the version of the file
    OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")
    call MoveTO(51,'*0   ')
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) FILEversion
    CLOSE (51, STATUS='KEEP') 
        
    if(FILEversion<inputVer) then
        CheckConstantsFile=.false.
    else
        CheckConstantsFile=.true.
    end if

end function CheckConstantsFile

  
  !-------------------------------------------------------READ AND WRITE CONSTANTS FILE---------------------------------------
subroutine CreateConstantsFile(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=516)::path

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
    write(51,"('*C   : ---- aTMDe-control parameters ----')")
    write(51,"('*p1  : initialize the modules by NP arrays')")
    write(51,*) initialize_NParrays
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
    write(51,"('*p3  : mass of top-quark')")
    write(51,*) mTOP
    write(51,"(' ')")
    write(51,"('*B   : ---- uPDF sets----')")
    write(51,"('*p1  : total number of PDFs to initialize (0= initialization is skipped)')")
    write(51,*) number_of_uPDFs
    write(51,"('*p2  : reference number for hadrons')")
    call writeShortIntegerList(51,enumeration_of_uPDFs)    
    write(51,"('*p3  : LHAPDF set names for hadrons (line-by-line corresponding to reference number')")
    do i=1,number_of_uPDFs
        write(51,*) trim(sets_of_uPDFs(i))
    end do
    write(51,"('*p4  : list of initialization replicas')")
    call writeShortIntegerList(51,replicas_of_uPDFs)

    write(51,"(' ')")
    write(51,"('*C   : ---- uFF sets----')")
    write(51,"('*p1  : total number of FFs to initialize (0= initialization is skipped)')")
    write(51,*) number_of_uFFs
    write(51,"('*p2  : reference number for hadrons')")
    call writeShortIntegerList(51,enumeration_of_uFFs)
    write(51,"('*p3  : LHAPDF set names for hadrons (line-by-line corresponding to reference number')")
    do i=1,number_of_uFFs
        write(51,*) trim(sets_of_uFFs(i))
    end do
    write(51,"('*p4  : list of initialization replicas')")
    call writeShortIntegerList(51,replicas_of_uFFs)

    write(51,"(' ')")
    write(51,"('*D   : ----lpPDF sets----')")
    write(51,"('*p1  : total number of PDFs to initialize (0= initialization is skipped)')")
    write(51,*) number_of_lpPDFs
    write(51,"('*p2  : reference number for hadrons')")
    call writeShortIntegerList(51,enumeration_of_lpPDFs)
    write(51,"('*p3  : LHAPDF set names for hadrons (line-by-line corresponding to reference number')")
    do i=1,number_of_lpPDFs
        write(51,*) trim(sets_of_lpPDFs(i))
    end do
    write(51,"('*p4  : list of initialization replicas')")
    call writeShortIntegerList(51,replicas_of_lpPDFs)
    
    write(51,"(' ')")
    write(51,"('*E   : ----hPDF sets----')")
    write(51,"('*p1  : total number of PDFs to initialize (0= initialization is skipped)')")
    write(51,*) number_of_hPDFs
    write(51,"('*p2  : reference number for hadrons')")
    call writeShortIntegerList(51,enumeration_of_hPDFs)
    write(51,"('*p3  : LHAPDF set names for hadrons (line-by-line corresponding to reference number')")
    do i=1,number_of_hPDFs
        write(51,*) trim(sets_of_hPDFs(i))
    end do
    write(51,"('*p4  : list of initialization replicas')")
    call writeShortIntegerList(51,replicas_of_hPDFs)


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
    write(51,"('*p4  : value of (alphaQED)^{-1} at MTAU')")
    write(51,*) alphaQED_MTAU
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
    write(51,"('*D   : ---- Higgs-boson ----')")
    write(51,"('*p1  : mass of Higgs-boson [GeV]')")
    write(51,*) mH
    write(51,"('*p2  : width of H-boson [GeV]')")
    write(51,*) GammaH
    write(51,"('*p3  : Vacuum expectation value (VEV) for Higgs potential [GeV]')")
    write(51,*) vevH
    write(51,"('*E   : ---- Leptons ----')")
    write(51,"('*p1  : mass of electron [GeV]')")
    write(51,*) mELECTRON
    write(51,"('*p2  : mass of muon [GeV]')")
    write(51,*) mMUON
    write(51,"('*p3  : mass of tau-lepton [GeV]')")
    write(51,*) mTAU

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
    write(51,"('*p2  : Initialization parameters (in column)')")
    do i=1,TMDR_lambdaLength
        write(51,*) TMDR_lambdaNP_init(i)
    end do
    write(51,*)
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
    write(51,"('*p2  : Use composite TMD-function definition')")
    write(51,*) uTMDPDF_IsComposite
    write(51,"('*B   : ---- Parameters of NP model ----')")
    write(51,"('*p1  : Length of lambdaNP')")
    write(51,*) uTMDPDF_lambdaLength
    write(51,"('*p2  : Initialization parameters (in column)')")
    do i=1,uTMDPDF_lambdaLength
        write(51,*) uTMDPDF_lambdaNP_init(i)
    end do
    write(51,"('*C   : ---- Numerical evaluation parameters ----')")
    write(51,"('*p1  : Tolerance (relative tolerance of convolution integral)')")
    write(51,*) uTMDPDF_tolerance
    write(51,"('*p2  : Maximum number of iterations (for adaptive integration)')")
    write(51,*) uTMDPDF_maxIteration
    write(51,"('*D   : ---- Grid preparation options ----')")
    write(51,"('*p1  : Prepare grid')")
    write(51,*) uTMDPDF_makeGrid
    write(51,"('*p2  : Include gluon TMDs into the grid')")
    write(51,*) uTMDPDF_withGluon
    write(51,"('*p3  : total number of PDFs added to the grid (by default it coincides with number of initialized PDFs)')")
    write(51,*) number_of_uPDFs
    write(51,"('*p4  : reference numbers for hadrons (by default it coincides with references for PDFs)')")
    call writeShortIntegerList(51,enumeration_of_uPDFs)
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
    write(51,"('*p2  : Use composite TMD-function definition')")
    write(51,*) uTMDFF_IsComposite
    write(51,"('*B   : ---- Parameters of NP model ----')")
    write(51,"('*p1  : Length of lambdaNP')")
    write(51,*) uTMDFF_lambdaLength
    write(51,"('*p2  : Initialization parameters (in column)')")
    do i=1,uTMDFF_lambdaLength
        write(51,*) uTMDFF_lambdaNP_init(i)
    end do
    write(51,"('*C   : ---- Numerical evaluation parameters ----')")
    write(51,"('*p1  : Tolerance (relative tolerance of convolution integral)')")
    write(51,*) uTMDFF_tolerance
    write(51,"('*p2  : Maximum number of iterations (for adaptive integration)')")
    write(51,*) uTMDFF_maxIteration
    write(51,"('*D   : ---- Grid preparation options ----')")
    write(51,"('*p1  : Prepare grid')")
    write(51,*) uTMDFF_makeGrid
    write(51,"('*p2  : Include gluon TMDs into the grid')")
    write(51,*) uTMDFF_withGluon
    write(51,"('*p3  : total number of FFs added to the grid (by default it coincides with number of initialized FFs)')")
    write(51,*) number_of_uFFs
    write(51,"('*p4  : reference numbers for hadrons (by default it coincides with references for FFs)')")
    call writeShortIntegerList(51,enumeration_of_uFFs)
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
    write(51,"('*B   : ---- Global garameters of structure functions----')")
    write(51,"('*p1  : Mass parameter used in the structure function (mass of hadron)')")
    write(51,*) TMDF_mass


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
    write(51,"('*p3  : Use resummation of pi^2-corrections in hard coefficient')")
    write(51,*) TMDX_DY_piResum
    write(51,"('*p4  : Use the exact values for factorization scales (include qT/Q correction)')")
    write(51,*) TMDX_DY_exactScale
    write(51,"('*B   : ---- Numerical evaluation parameters ----')")
    write(51,"('*p1  : Tolerance (relative tolerance for bin-integration routines, except pt-integration)')")
    write(51,*) TMDX_DY_tolerance
    write(51,"('*p2  : Minimal number of sections for pt-integration')")
    write(51,*) TMDX_DY_ptSECTION
    write(51,"('*p3  : Maximum size of the Q-bin integration (larger bins are desected)')")
    write(51,*) TMDX_DY_maxQbinSize
    write(51,"('*C   : ---- Parameters for parallel evaluation (used only if compiled with openMP) ----')")
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
    write(51,"('*p5  : Use transverse momentum corrections in x1 and z1')")
    write(51,*) TMDX_SIDIS_qTinX1Z1corr
    write(51,"('*p6  : Use the exact values for factorization scales (include qT/Q correction)')")
    write(51,*) TMDX_SIDIS_exactScale
    write(51,"('*B   : ---- Numerical evaluation parameters ----')")
    write(51,"('*p1  : Tolerance (relative tolerance for bin-integration routines, except pt-integration)')")
    write(51,*) TMDX_SIDIS_tolerance
    write(51,"('*p2  : Minimal number of sections for pt-integration')")
    write(51,*) TMDX_SIDIS_ptSECTION
    write(51,"('*p3  : Tolerance for Z-integration (relative tolerance for Z-bin-integration routines)')")
    write(51,*) TMDX_SIDIS_toleranceZ
    write(51,"('*p4  : Method for Z-bin integration (see manual)')")
    write(51,*) TMDX_SIDIS_methodZ
    write(51,"('*p5  : Tolerance for X-integration (relative tolerance for X-bin-integration routines)')")
    write(51,*) TMDX_SIDIS_toleranceX
    write(51,"('*p6  : Method for X-bin integration (see manual)')")
    write(51,*) TMDX_SIDIS_methodX
    write(51,"('*C   : ---- Parameters for parallel evaluation (used only if compiled with openMP) ----')")
    write(51,"('*p1  : Maximum number of processors to use')")
    write(51,*) TMDX_SIDIS_numProc


    write(51,"(' ')")
    write(51,"(' ')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('# ----                           PARAMETERS OF lpTMDPDF                 -----')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('*11  :')")
    write(51,"('*p1  : initialize lpTMDPDF module')")
    write(51,*) include_lpTMDPDF
    write(51,"(' ')")
    write(51,"('*A   : ---- Main definitions ----')")
    write(51,"('*p1  : Order of coefficient function')")
    write(51,*) trim(lpTMDPDF_order)
    write(51,"('*p2  : Use composite TMD-function definition')")
    write(51,*) lpTMDPDF_IsComposite
    write(51,"('*B   : ---- Parameters of NP model ----')")
    write(51,"('*p1  : Length of lambdaNP')")
    write(51,*) lpTMDPDF_lambdaLength
    write(51,"('*p2  : Initialization parameters (in column)')")
    do i=1,lpTMDPDF_lambdaLength
        write(51,*) lpTMDPDF_lambdaNP_init(i)
    end do
    write(51,"('*C   : ---- Numerical evaluation parameters ----')")
    write(51,"('*p1  : Tolerance (relative tolerance of convolution integral)')")
    write(51,*) lpTMDPDF_tolerance
    write(51,"('*p2  : Maximum number of iterations (for adaptive integration)')")
    write(51,*) lpTMDPDF_maxIteration
    write(51,"('*D   : ---- Grid preparation options ----')")
    write(51,"('*p1  : Prepare grid')")
    write(51,*) lpTMDPDF_makeGrid
    write(51,"('*p2  : Include gluon TMDs into the grid')")
    write(51,*) lpTMDPDF_withGluon
    write(51,"('*p3  : total number of PDFs added to the grid (by default it coincides with number of initialized PDFs)')")
    write(51,*) number_of_lpPDFs
    write(51,"('*p4  : reference numbers for hadrons (by default it coincides with references for PDFs)')")
    call writeShortIntegerList(51,enumeration_of_lpPDFs)
    write(51,"('*E   : ---- Parameters of grid ----')")
    write(51,"('*p1  : xGrid_Min the minimal value of x in grid (max=1), make sure that it is enough)')")
    write(51,*) lpTMDPDF_grid_xMin
    write(51,"('*p2  : the maximum bT in grid (min=0), for larger approximate extrapolation is done')")
    write(51,*) lpTMDPDF_grid_bMax
    write(51,"('*p3  : GridSizeX (250 is enough for 10^-8 presicion, for grid up to 10^-5)')")
    write(51,*) lpTMDPDF_grid_SizeX
    write(51,"('*p4  : GridSizeB (750 is enough for ~10^-6 presicion, 500 for 10^-5 presicion)')")
    write(51,*) lpTMDPDF_grid_SizeB
    write(51,"('*p5  : slope scale of griding at smaller b (better not to change it :) )')")
    write(51,*) lpTMDPDF_grid_slope
    
    write(51,"(' ')")
    write(51,"(' ')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('# ----                         PARAMETERS OF SiversTMDPDF               -----')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('*12  :')")
    write(51,"('*p1  : initialize SiversTMDPDF module')")
    write(51,*) include_SiversTMDPDF
    write(51,"(' ')")
    write(51,"('*A   : ---- Main definitions ----')")
    write(51,"('*p1  : Order of coefficient function')")
    write(51,*) trim(SiversTMDPDF_order)
    write(51,"('*p2  : Use composite TMD-function definition')")
    write(51,*) SiversTMDPDF_IsComposite
    write(51,"('*B   : ---- Parameters of NP model ----')")
    write(51,"('*p1  : Length of lambdaNP')")
    write(51,*) SiversTMDPDF_lambdaLength
    write(51,"('*p2  : Initialization parameters (in column)')")
    do i=1,SiversTMDPDF_lambdaLength
        write(51,*) SiversTMDPDF_lambdaNP_init(i)
    end do
    write(51,"('*C   : ---- Numerical evaluation parameters ----')")
    write(51,"('*p1  : Tolerance (relative tolerance of convolution integral)')")
    write(51,*) SiversTMDPDF_tolerance
    write(51,"('*p2  : Maximum number of iterations (for adaptive integration)')")
    write(51,*) SiversTMDPDF_maxIteration
    write(51,"('*D   : ---- Grid preparation options ----')")
    write(51,"('*p1  : Prepare grid')")
    write(51,*) SiversTMDPDF_makeGrid
    write(51,"('*p2  : Include gluon TMDs into the grid')")
    write(51,*) SiversTMDPDF_withGluon
    write(51,"('*p3  : total number of PDFs added to the grid')")
    write(51,*) number_of_SiversTMDPDFs
    write(51,"('*p4  : reference numbers for hadrons')")
    call writeShortIntegerList(51,enumeration_of_SiversTMDPDFs)
    write(51,"('*E   : ---- Parameters of grid ----')")
    write(51,"('*p1  : xGrid_Min the minimal value of x in grid (max=1), make sure that it is enough)')")
    write(51,*) SiversTMDPDF_grid_xMin
    write(51,"('*p2  : the maximum bT in grid (min=0), for larger approximate extrapolation is done')")
    write(51,*) SiversTMDPDF_grid_bMax
    write(51,"('*p3  : GridSizeX (250 is enough for 10^-8 presicion, for grid up to 10^-5)')")
    write(51,*) SiversTMDPDF_grid_SizeX
    write(51,"('*p4  : GridSizeB (750 is enough for ~10^-6 presicion, 500 for 10^-5 presicion)')")
    write(51,*) SiversTMDPDF_grid_SizeB
    write(51,"('*p5  : slope scale of griding at smaller b (better not to change it :) )')")
    write(51,*) SiversTMDPDF_grid_slope
    
    write(51,"(' ')")
    write(51,"(' ')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('# ----                         PARAMETERS OF wgtTMDPDF                  -----')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('*13  :')")
    write(51,"('*p1  : initialize wgtTMDPDF module')")
    write(51,*) include_wgtTMDPDF
    write(51,"(' ')")
    write(51,"('*A   : ---- Main definitions ----')")
    write(51,"('*p1  : Order of coefficient function (for tw-2 part only)')")
    write(51,*) trim(wgtTMDPDF_order)
    write(51,"('*p2  : Use composite TMD-function definition')")
    write(51,*) wgtTMDPDF_IsComposite
    write(51,"('*B   : ---- Parameters of NP model ----')")
    write(51,"('*p1  : Length of lambdaNP (used in fNP and tw-3 part)')")
    write(51,*) wgtTMDPDF_lambdaLength
    write(51,"('*p2  : Initialization parameters (in column)')")
    do i=1,wgtTMDPDF_lambdaLength
        write(51,*) wgtTMDPDF_lambdaNP_init(i)
    end do
    write(51,"('*C   : ---- Numerical evaluation parameters ----')")
    write(51,"('*p1  : Tolerance (relative tolerance of convolution integral)')")
    write(51,*) wgtTMDPDF_tolerance
    write(51,"('*p2  : Maximum number of iterations (for adaptive integration)')")
    write(51,*) wgtTMDPDF_maxIteration
    write(51,"('*D   : ---- Grid preparation options ----')")
    write(51,"('*p1  : Prepare grid')")
    write(51,*) wgtTMDPDF_makeGrid
    write(51,"('*p2  : Include gluon TMDs into the grid')")
    write(51,*) wgtTMDPDF_withGluon
    write(51,"('*p3  : total number of PDFs added to the grid')")
    write(51,*) number_of_wgtTMDPDFs
    write(51,"('*p4  : reference numbers for hadrons')")
    call writeShortIntegerList(51,enumeration_of_wgtTMDPDFs)
    write(51,"('*E   : ---- Parameters of grid ----')")
    write(51,"('*p1  : xGrid_Min the minimal value of x in grid (max=1), make sure that it is enough)')")
    write(51,*) wgtTMDPDF_grid_xMin
    write(51,"('*p2  : the maximum bT in grid (min=0), for larger approximate extrapolation is done')")
    write(51,*) wgtTMDPDF_grid_bMax
    write(51,"('*p3  : GridSizeX (250 is enough for 10^-8 presicion, for grid up to 10^-5)')")
    write(51,*) wgtTMDPDF_grid_SizeX
    write(51,"('*p4  : GridSizeB (750 is enough for ~10^-6 presicion, 500 for 10^-5 presicion)')")
    write(51,*) wgtTMDPDF_grid_SizeB
    write(51,"('*p5  : slope scale of griding at smaller b (better not to change it :) )')")
    write(51,*) wgtTMDPDF_grid_slope

    CLOSE (51, STATUS='KEEP') 

    if(outputLevel>1) write(*,*) 'aTMDe_setup: Constans file is made.'

end subroutine CreateConstantsFile
  
  
subroutine ReadConstantsFile(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=516)::path
    !!!! this is version of input file. it is read first and then result is compared with the current ID
    !!!! It suppose to make compatibility betwen versions
    integer::FILEversion
    integer::i   
    logical::file_exists
    
    

    if(outputLevel>2) write(*,*) '-----------------------------------'
    if(outputLevel>1) write(*,*) 'aTMDe_setup: Reading setup file ...'
    if(present(prefix)) then
        path=trim(adjustl(prefix))//trim(adjustr(file))
    else
        path=trim(adjustr(file))
    end if
    if(outputLevel>1) write(*,*)'       path:', color(trim(path),c_yellow)
    
    INQUIRE(FILE=path, EXIST=file_exists)
    if(.not.file_exists) then
        write(*,*) color('Constans file is not found.',c_red_bold)
        return
    end if

    OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")
    !# ----                           GLOBAL PARAMETERS                      -----
    call MoveTO(51,'*0   ')
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) FILEversion

    if(FILEversion>inputVer) then
        if(outputLevel>0) write(*,*) color('aTMDe_setup: Version of input file(',c_red_bold), &
                                        color(int4ToStr(FILEversion),c_red_bold),&
                                        color(') newer than the current version of code(',c_red_bold),&
                                        color(int4ToStr(inputVer),c_red_bold),&
                                        color(')',c_red_bold)
        if(outputLevel>0) write(*,*) color('aTMDe_setup: UPDATE ARTEMIDE!',c_red_bold)
        if(outputLevel>0) write(*,*) color('aTMDe_setup: UPDATE ARTEMIDE!',c_red_bold)
        if(outputLevel>0) write(*,*) 'aTMDe_setup: Attempt to load...'
    else if(FILEversion<inputVer) then
        if(outputLevel>0) write(*,*) color('aTMDe_setup: Version of input file(',c_red_bold), &
                                        color(int4ToStr(FILEversion),c_red_bold),&
                                        color(') older than the current version of code(',c_red_bold),&
                                        color(int4ToStr(inputVer),c_red_bold),&
                                        color(')',c_red_bold)
        if(outputLevel>0) write(*,*) 'aTMDe_setup: Attempt to load...'
    end if

    call MoveTO(51,'*p2  ')
    read(51,*) outputLevel
    call MoveTO(51,'*p3  ')
    read(51,*) messageTrigger

    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) hc2

    if(FILEversion>=11) then
        call MoveTO(51,'*C   ')
        call MoveTO(51,'*p1  ')
        read(51,*) initialize_NParrays
    end if

    !# ----                           PARAMETERS OF QCDinput                 -----
    call MoveTO(51,'*1   ')
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) mCHARM
    call MoveTO(51,'*p2  ')
    read(51,*) mBOTTOM
    if(FILEversion>3) then
        call MoveTO(51,'*p3  ')
        read(51,*) mTOP
    end if
    !-------PDF for uTMDPDF
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

    !-------FF for uTMDFF
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

    if(FILEversion>=3) then
        !-------PDF for lpTMDPDF
        call MoveTO(51,'*D   ')
        call MoveTO(51,'*p1  ')

        read(51,*) number_of_lpPDFs
        if(allocated(enumeration_of_lpPDFs)) deallocate(enumeration_of_lpPDFs)
        if(allocated(replicas_of_lpPDFs)) deallocate(replicas_of_lpPDFs)
        if(allocated(sets_of_lpPDFs)) deallocate(sets_of_lpPDFs)
        allocate(enumeration_of_lpPDFs(1:number_of_lpPDFs))
        allocate(sets_of_lpPDFs(1:number_of_lpPDFs))
        allocate(replicas_of_lpPDFs(1:number_of_lpPDFs))    

        call MoveTO(51,'*p2  ')
        read(51,*) enumeration_of_lpPDFs
        call MoveTO(51,'*p3  ')
        do i=1,number_of_lpPDFs
            read(51,*) sets_of_lpPDFs(i)
        end do
        call MoveTO(51,'*p4  ')
        read(51,*) replicas_of_lpPDFs
    end if
    
    if(FILEversion>=18) then
        !-------helicity PDF
        call MoveTO(51,'*E   ')
        call MoveTO(51,'*p1  ')

        read(51,*) number_of_hPDFs
        if(allocated(enumeration_of_hPDFs)) deallocate(enumeration_of_hPDFs)
        if(allocated(replicas_of_hPDFs)) deallocate(replicas_of_hPDFs)
        if(allocated(sets_of_hPDFs)) deallocate(sets_of_hPDFs)
        allocate(enumeration_of_hPDFs(1:number_of_hPDFs))
        allocate(sets_of_hPDFs(1:number_of_hPDFs))
        allocate(replicas_of_hPDFs(1:number_of_hPDFs))    

        call MoveTO(51,'*p2  ')
        read(51,*) enumeration_of_hPDFs
        call MoveTO(51,'*p3  ')
        do i=1,number_of_hPDFs
            read(51,*) sets_of_hPDFs(i)
        end do
        call MoveTO(51,'*p4  ')
        read(51,*) replicas_of_hPDFs
    end if

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
    if(FILEversion>=14) then
        call MoveTO(51,'*p4  ')
        read(51,*) alphaQED_MTAU
    end if
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
    if(FILEversion>=5) then
        call MoveTO(51,'*D   ')
        call MoveTO(51,'*p1  ')
        read(51,*) mH
        call MoveTO(51,'*p2  ')
        read(51,*) GammaH
        call MoveTO(51,'*p3  ')
        read(51,*) vevH
    end if
    if(FILEversion>=13) then
        call MoveTO(51,'*E   ')
        call MoveTO(51,'*p1  ')
        read(51,*) mELECTRON
        call MoveTO(51,'*p2  ')
        read(51,*) mMUON
        call MoveTO(51,'*p3  ')
        read(51,*) mTAU
    end if


    !# ----                            PARAMETERS OF TMDR                    -----
    call MoveTO(51,'*3   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_TMDR
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDR_order
    call MoveTO(51,'*p2  ')
    read(51,*) TMDR_evolutionType
    if((inputVer<2 .and. TMDR_evolutionType>=4) .or. (inputVer>=14 .and. TMDR_evolutionType>=4)) then
        if(outputLevel>0) write(*,*) color('aTMDe_setup: mismatch in TMDR_evolutionType (3-A-p2). Set to 3',c_red)
        TMDR_evolutionType=3
    end if

    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDR_lambdaLength
    call ReNewInitializationArray(TMDR_lambdaNP_init,TMDR_lambdaLength)!!!!we need to redo (and reallocate) the initialization of the NP-array
    if(FILEversion>=10) then
        call MoveTO(51,'*p2  ')      
        do i=1,TMDR_lambdaLength
            read(51,*) TMDR_lambdaNP_init(i)
        end do
    end if
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
    if(FILEversion>=12) then
        call MoveTO(51,'*p2  ')
        read(51,*) uTMDPDF_IsComposite
    end if
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) uTMDPDF_lambdaLength
    call ReNewInitializationArray(uTMDPDF_lambdaNP_init,uTMDPDF_lambdaLength)!!!!we need to redo (and reallocate) the initialization of the NP-array
    if(FILEversion>=10) then
        call MoveTO(51,'*p2  ')      
        do i=1,uTMDPDF_lambdaLength
    read(51,*) uTMDPDF_lambdaNP_init(i)
        end do
    end if
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
        if(outputLevel>0) write(*,*) ' '
        if(outputLevel>0) write(*,*) color('ESSENTIAL INCONSITENCY: the number of uPDFs is unequal to the number of uTMDPDFs',c_red)
        if(outputLevel>0) write(*,*) color('                        information on the number of uTMDPDFs is ignored',c_red)
        if(outputLevel>0) write(*,*) ' '
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
    if(FILEversion>=12) then
        call MoveTO(51,'*p2  ')
        read(51,*) uTMDFF_IsComposite
    end if
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) uTMDFF_lambdaLength
    call ReNewInitializationArray(uTMDFF_lambdaNP_init,uTMDFF_lambdaLength)!!!!we need to redo (and reallocate) the initialization of the NP-array
    if(FILEversion>=10) then
        call MoveTO(51,'*p2  ')      
        do i=1,uTMDFF_lambdaLength
    read(51,*) uTMDFF_lambdaNP_init(i)
        end do
    end if
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
        if(outputLevel>0) write(*,*) ' '
        if(outputLevel>0) write(*,*) color('ESSENTIAL INCONSITENCY: number of uFFs unequal uTMDFFs',c_red)
        if(outputLevel>0) write(*,*) color('                        information on number of uTMDFFs ignored',c_red)
        if(outputLevel>0) write(*,*) ' '
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
    if(FILEversion>=16) then
        call MoveTO(51,'*B   ')
        call MoveTO(51,'*p1  ')
        read(51,*) TMDF_mass
    end if


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
    if(FILEversion>=6) then
        call MoveTO(51,'*p3  ')
        read(51,*) TMDX_DY_piResum
    end if
    if(FILEversion>=17) then
        call MoveTO(51,'*p4  ')
        read(51,*) TMDX_DY_exactScale
    end if
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDX_DY_tolerance
    call MoveTO(51,'*p2  ')
    read(51,*) TMDX_DY_ptSECTION
    if(FILEversion>=21) then
        call MoveTO(51,'*p3  ')
        read(51,*) TMDX_DY_maxQbinSize
    end if
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
    if(FILEversion>=9) then
        call MoveTO(51,'*p5  ')
        read(51,*) TMDX_SIDIS_qTinX1Z1corr
    end if
    if(FILEversion>=20) then
        call MoveTO(51,'*p6  ')
        read(51,*) TMDX_SIDIS_exactScale
    end if
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDX_SIDIS_tolerance
    call MoveTO(51,'*p2  ')
    read(51,*) TMDX_SIDIS_ptSECTION
    if(FILEversion>=7) then
        call MoveTO(51,'*p3  ')
        read(51,*) TMDX_SIDIS_toleranceZ
        call MoveTO(51,'*p4  ')
        read(51,*) TMDX_SIDIS_methodZ
    end if
    if(FILEversion>=8) then
        call MoveTO(51,'*p5  ')
        read(51,*) TMDX_SIDIS_toleranceX
        call MoveTO(51,'*p6  ')
        read(51,*) TMDX_SIDIS_methodX
    end if
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDX_SIDIS_numProc

    if(FILEversion>=3) then
    !# ----                           PARAMETERS OF lpTMDPDF                  -----
        call MoveTO(51,'*11  ')
        call MoveTO(51,'*p1  ')
        read(51,*) include_lpTMDPDF
        call MoveTO(51,'*A   ')
        call MoveTO(51,'*p1  ')
        read(51,*) lpTMDPDF_order
        if(FILEversion>=12) then
            call MoveTO(51,'*p2  ')
            read(51,*) lpTMDPDF_IsComposite
        end if
        call MoveTO(51,'*B   ')
        call MoveTO(51,'*p1  ')
        read(51,*) lpTMDPDF_lambdaLength
        call ReNewInitializationArray(lpTMDPDF_lambdaNP_init,lpTMDPDF_lambdaLength)!!!!we need to redo (and reallocate) the initialization of the NP-array
        if(FILEversion>=10) then
            call MoveTO(51,'*p2  ')      
            do i=1,lpTMDPDF_lambdaLength
        read(51,*) lpTMDPDF_lambdaNP_init(i)
            end do
        end if
        call MoveTO(51,'*C   ')
        call MoveTO(51,'*p1  ')
        read(51,*) lpTMDPDF_tolerance
        call MoveTO(51,'*p2  ')
        read(51,*) lpTMDPDF_maxIteration
        call MoveTO(51,'*D   ')
        call MoveTO(51,'*p1  ')
        read(51,*) lpTMDPDF_makeGrid
        call MoveTO(51,'*p2  ')
        read(51,*) lpTMDPDF_withGluon
        call MoveTO(51,'*p3  ')
        read(51,*) i
        if(i/=number_of_lpPDFs) then
            if(outputLevel>0) write(*,*) ' '
            if(outputLevel>0) write(*,*) color(&
                                        'ESSENTIAL INCONSITENCY: the number of lpPDFs is unequal to the number of lpTMDPDFs',c_red)
            if(outputLevel>0) write(*,*) color('                        information on the number of lpTMDPDFs is ignored',c_red)
            if(outputLevel>0) write(*,*) ' '
        end if
        call MoveTO(51,'*E   ')
        call MoveTO(51,'*p1  ')
        read(51,*) lpTMDPDF_grid_xMin
        call MoveTO(51,'*p2  ')
        read(51,*) lpTMDPDF_grid_bMax
        call MoveTO(51,'*p3  ')
        read(51,*) lpTMDPDF_grid_SizeX
        call MoveTO(51,'*p4  ')
        read(51,*) lpTMDPDF_grid_SizeB
        call MoveTO(51,'*p5  ')
        read(51,*) lpTMDPDF_grid_slope
    else
        write(*,*) 'aTMDe_setup: parameters of lin.pol.gluons set default. (const-ver. < 3)'
    end if
    
    if(FILEversion>=15) then
    !# ----                           PARAMETERS OF SiversTMDPDF                  -----
        call MoveTO(51,'*12  ')
        call MoveTO(51,'*p1  ')
        read(51,*) include_SiversTMDPDF
        call MoveTO(51,'*A   ')
        call MoveTO(51,'*p1  ')
        read(51,*) SiversTMDPDF_order        
        call MoveTO(51,'*p2  ')
        read(51,*) SiversTMDPDF_IsComposite
        call MoveTO(51,'*B   ')
        call MoveTO(51,'*p1  ')
        read(51,*) SiversTMDPDF_lambdaLength
        call ReNewInitializationArray(SiversTMDPDF_lambdaNP_init,SiversTMDPDF_lambdaLength)!!!!we need to redo (and reallocate) the initialization of the NP-array        
        call MoveTO(51,'*p2  ')      
        do i=1,SiversTMDPDF_lambdaLength
            read(51,*) SiversTMDPDF_lambdaNP_init(i)
        end do
        call MoveTO(51,'*C   ')
        call MoveTO(51,'*p1  ')
        read(51,*) SiversTMDPDF_tolerance
        call MoveTO(51,'*p2  ')
        read(51,*) SiversTMDPDF_maxIteration
        call MoveTO(51,'*D   ')
        call MoveTO(51,'*p1  ')
        read(51,*) SiversTMDPDF_makeGrid
        call MoveTO(51,'*p2  ')
        read(51,*) SiversTMDPDF_withGluon
        call MoveTO(51,'*p3  ')
        read(51,*) number_of_SiversTMDPDFs        
        if(allocated(enumeration_of_SiversTMDPDFs)) deallocate(enumeration_of_SiversTMDPDFs)
        allocate(enumeration_of_SiversTMDPDFs(1:number_of_SiversTMDPDFs))
        call MoveTO(51,'*p4  ')
        read(51,*) enumeration_of_SiversTMDPDFs
        call MoveTO(51,'*E   ')
        call MoveTO(51,'*p1  ')
        read(51,*) SiversTMDPDF_grid_xMin
        call MoveTO(51,'*p2  ')
        read(51,*) SiversTMDPDF_grid_bMax
        call MoveTO(51,'*p3  ')
        read(51,*) SiversTMDPDF_grid_SizeX
        call MoveTO(51,'*p4  ')
        read(51,*) SiversTMDPDF_grid_SizeB
        call MoveTO(51,'*p5  ')
        read(51,*) SiversTMDPDF_grid_slope
    else
        write(*,*) 'aTMDe_setup: parameters of Sivers functions set default. (const-ver. < 15)'
    end if
    
    if(FILEversion>=19) then
    !# ----                           PARAMETERS OF wgtTMDPDF                  -----
        call MoveTO(51,'*13  ')
        call MoveTO(51,'*p1  ')
        read(51,*) include_wgtTMDPDF
        call MoveTO(51,'*A   ')
        call MoveTO(51,'*p1  ')
        read(51,*) wgtTMDPDF_order        
        call MoveTO(51,'*p2  ')
        read(51,*) wgtTMDPDF_IsComposite
        call MoveTO(51,'*B   ')
        call MoveTO(51,'*p1  ')
        read(51,*) wgtTMDPDF_lambdaLength
        call ReNewInitializationArray(wgtTMDPDF_lambdaNP_init,wgtTMDPDF_lambdaLength)!!!!we need to redo (and reallocate) the initialization of the NP-array        
        call MoveTO(51,'*p2  ')      
        do i=1,wgtTMDPDF_lambdaLength
            read(51,*) wgtTMDPDF_lambdaNP_init(i)
        end do
        call MoveTO(51,'*C   ')
        call MoveTO(51,'*p1  ')
        read(51,*) wgtTMDPDF_tolerance
        call MoveTO(51,'*p2  ')
        read(51,*) wgtTMDPDF_maxIteration
        call MoveTO(51,'*D   ')
        call MoveTO(51,'*p1  ')
        read(51,*) wgtTMDPDF_makeGrid
        call MoveTO(51,'*p2  ')
        read(51,*) wgtTMDPDF_withGluon
        call MoveTO(51,'*p3  ')
        read(51,*) number_of_wgtTMDPDFs        
        if(allocated(enumeration_of_wgtTMDPDFs)) deallocate(enumeration_of_wgtTMDPDFs)
        allocate(enumeration_of_wgtTMDPDFs(1:number_of_wgtTMDPDFs))
        call MoveTO(51,'*p4  ')
        read(51,*) enumeration_of_wgtTMDPDFs     
        call MoveTO(51,'*E   ')
        call MoveTO(51,'*p1  ')
        read(51,*) wgtTMDPDF_grid_xMin
        call MoveTO(51,'*p2  ')
        read(51,*) wgtTMDPDF_grid_bMax
        call MoveTO(51,'*p3  ')
        read(51,*) wgtTMDPDF_grid_SizeX
        call MoveTO(51,'*p4  ')
        read(51,*) wgtTMDPDF_grid_SizeB
        call MoveTO(51,'*p5  ')
        read(51,*) wgtTMDPDF_grid_slope
    else
        write(*,*) 'aTMDe_setup: parameters of worm-geat T functions set default. (const-ver. < 19)'
    end if


    CLOSE (51, STATUS='KEEP') 

    if(outputLevel>1) write(*,*) color('aTMDe_setup: constants-file loaded sucessfully.',c_green_bold)

end subroutine ReadConstantsFile
  
  !--------------------------------------------------------------
  !!!! this subtroutine kill re allocate the array, and sets its values as (2,0,0,0,...)
  !!!! used for empty initialization arrays
  subroutine ReNewInitializationArray(arr,n)
    real(dp),allocatable,intent(out)::arr(:)
    integer::n,i
    if(allocated(arr)) deallocate(arr)
    allocate(arr(1:n))
    arr(1)=2d0
    do i=2,n
	arr(i)=0d0
    end do
  end subroutine ReNewInitializationArray
end module aTMDe_Setup
