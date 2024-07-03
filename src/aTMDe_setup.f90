!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe-setup 2.0
!
!    Module that read, create and modify constants-files
!
!    if you use this module please, quote 1803.11089
!
!                A.Vladimirov (27.05.2019)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module aTMDe_Setup
use aTMDe_Numerics
use IO_functions
implicit none

private

character (len=5),parameter :: version="v3.00"
character (len=11),parameter :: moduleName="aTMDe-setup"
!! actual version of input file
integer,parameter::inputVer=30

!detalization of output: 0 = no output except critical, 1 = + WARNINGS, 2 = + states of initialization,sets,etc, 3 = + details
integer::outputLevel
integer::messageTrigger
integer::number_of_processors !!!! number of processors which could be used by openMP

!-----------------------physic parameters
logical::include_EWinput
!mass, and width parameters
real(dp)::mCHARM,mBOTTOM,mTOP,mZ,mW,GammaZ,GammaW,mH,GammaH,vevH
real(dp)::mELECTRON,mMUON,mTAU
!other physical parameters
real(dp)::hc2,Mnorm,alphaQED_MZ,sW2,alphaQED_MTAU
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
logical::TMDR_override
character*8::TMDR_orderGAMMA,TMDR_orderV,TMDR_orderD,TMDR_orderZETA
integer::TMDR_lambdaLength
real(dp)::TMDR_tolerance
real(dp)::TMDR_bFREEZE
real(dp)::TMDR_smooth


!-------------------- uTMDPDF parameters
logical::include_uTMDPDF
character*8::uTMDPDF_order
logical::uTMDPDF_makeGrid,uTMDPDF_withGluon,uTMDPDF_runGridTest
integer::uTMDPDF_lambdaLength
integer::uTMDPDF_numHadron
real(dp)::uTMDPDF_BMAX_ABS
real(dp)::uTMDPDF_toleranceINT
real(dp)::uTMDPDF_toleranceGEN
integer::uTMDPDF_maxIteration
integer::uTMDPDF_numSubGridsX,uTMDPDF_numSubGridsB
real(dp),allocatable::uTMDPDF_subGridsX(:),uTMDPDF_subGridsB(:)
integer::uTMDPDF_grid_SizeX,uTMDPDF_grid_SizeB
logical::uTMDPDF_makeGrid_inKT,uTMDPDF_runGridTest_inKT
integer::uTMDPDF_numSubGridsX_inKT,uTMDPDF_numSubGridsKT_inKT,uTMDPDF_numSubGridsB_inKT
real(dp),allocatable::uTMDPDF_subGridsX_inKT(:),uTMDPDF_subGridsB_inKT(:),uTMDPDF_subGridsKT_inKT(:)
integer::uTMDPDF_grid_SizeX_inKT,uTMDPDF_grid_SizeB_inKT,uTMDPDF_grid_SizeKT_inKT
real(dp)::uTMDPDF_minQ_inKT,uTMDPDF_maxQ_inKT
integer::uTMDPDF_grid_SizeQ_inKT
real(dp)::uTMDPDF_hOGATA_TMM,uTMDPDF_toleranceOGATA_TMM,uTMDPDF_muMIN_TMM


!-------------------- uTMDFF parameters
logical::include_uTMDFF
character*8::uTMDFF_order
logical::uTMDFF_makeGrid,uTMDFF_withGluon,uTMDFF_runGridTest
integer::uTMDFF_lambdaLength
integer::uTMDFF_numHadron
real(dp)::uTMDFF_BMAX_ABS
real(dp)::uTMDFF_toleranceINT
real(dp)::uTMDFF_toleranceGEN
integer::uTMDFF_maxIteration
integer::uTMDFF_numSubGridsX,uTMDFF_numSubGridsB
real(dp),allocatable::uTMDFF_subGridsX(:),uTMDFF_subGridsB(:)
integer::uTMDFF_grid_SizeX,uTMDFF_grid_SizeB
real(dp)::uTMDFF_hOGATA,uTMDFF_toleranceOGATA,uTMDFF_KT_FREEZE
real(dp)::uTMDFF_hOGATA_TMM,uTMDFF_toleranceOGATA_TMM,uTMDFF_muMIN_TMM

!-------------------- lpTMDPDF parameters
logical::include_lpTMDPDF
character*8::lpTMDPDF_order
logical::lpTMDPDF_makeGrid,lpTMDPDF_withGluon,lpTMDPDF_runGridTest
integer::lpTMDPDF_lambdaLength
integer::lpTMDPDF_numHadron
real(dp)::lpTMDPDF_BMAX_ABS
real(dp)::lpTMDPDF_toleranceINT
real(dp)::lpTMDPDF_toleranceGEN
integer::lpTMDPDF_maxIteration
integer::lpTMDPDF_numSubGridsX,lpTMDPDF_numSubGridsB
real(dp),allocatable::lpTMDPDF_subGridsX(:),lpTMDPDF_subGridsB(:)
integer::lpTMDPDF_grid_SizeX,lpTMDPDF_grid_SizeB
real(dp)::lpTMDPDF_hOGATA,lpTMDPDF_toleranceOGATA,lpTMDPDF_KT_FREEZE
real(dp)::lpTMDPDF_hOGATA_TMM,lpTMDPDF_toleranceOGATA_TMM,lpTMDPDF_muMIN_TMM

!-------------------- SiversTMDPDF parameters
logical::include_SiversTMDPDF
character*8::SiversTMDPDF_order
logical::SiversTMDPDF_makeGrid,SiversTMDPDF_withGluon,SiversTMDPDF_runGridTest
integer::SiversTMDPDF_lambdaLength
integer::SiversTMDPDF_numHadron
real(dp)::SiversTMDPDF_BMAX_ABS
real(dp)::SiversTMDPDF_toleranceINT
real(dp)::SiversTMDPDF_toleranceGEN
integer::SiversTMDPDF_maxIteration
real(dp)::SiversTMDPDF_hOGATA,SiversTMDPDF_toleranceOGATA,SiversTMDPDF_KT_FREEZE
real(dp)::SiversTMDPDF_hOGATA_TMM,SiversTMDPDF_toleranceOGATA_TMM,SiversTMDPDF_muMIN_TMM

!-------------------- wgtTMDPDF parameters
logical::include_wgtTMDPDF
character*8::wgtTMDPDF_order
logical::wgtTMDPDF_makeGrid,wgtTMDPDF_withGluon,wgtTMDPDF_runGridTest
integer::wgtTMDPDF_numHadron
character*8::wgtTMDPDF_order_tw3
logical::wgtTMDPDF_makeGrid_tw3,wgtTMDPDF_withGluon_tw3,wgtTMDPDF_runGridTest_tw3
integer::wgtTMDPDF_lambdaLength
real(dp)::wgtTMDPDF_BMAX_ABS
real(dp)::wgtTMDPDF_toleranceINT
real(dp)::wgtTMDPDF_toleranceGEN
integer::wgtTMDPDF_maxIteration
integer::wgtTMDPDF_numSubGridsX,wgtTMDPDF_numSubGridsB
real(dp),allocatable::wgtTMDPDF_subGridsX(:),wgtTMDPDF_subGridsB(:)
integer::wgtTMDPDF_grid_SizeX,wgtTMDPDF_grid_SizeB
real(dp)::wgtTMDPDF_hOGATA,wgtTMDPDF_toleranceOGATA,wgtTMDPDF_KT_FREEZE
real(dp)::wgtTMDPDF_hOGATA_TMM,wgtTMDPDF_toleranceOGATA_TMM,wgtTMDPDF_muMIN_TMM

!-------------------- BoerMuldersTMDPDF parameters
logical::include_BoerMuldersTMDPDF
character*8::BoerMuldersTMDPDF_order
logical::BoerMuldersTMDPDF_makeGrid,BoerMuldersTMDPDF_withGluon,BoerMuldersTMDPDF_runGridTest
integer::BoerMuldersTMDPDF_lambdaLength
integer::BoerMuldersTMDPDF_numHadron
real(dp)::BoerMuldersTMDPDF_BMAX_ABS
real(dp)::BoerMuldersTMDPDF_toleranceINT
real(dp)::BoerMuldersTMDPDF_toleranceGEN
integer::BoerMuldersTMDPDF_maxIteration
logical::BoerMuldersTMDPDF_makeGrid_inKT,BoerMuldersTMDPDF_runGridTest_inKT
integer::BoerMuldersTMDPDF_numSubGridsX_inKT,BoerMuldersTMDPDF_numSubGridsKT_inKT,BoerMuldersTMDPDF_numSubGridsB_inKT
real(dp),allocatable::BoerMuldersTMDPDF_subGridsX_inKT(:),BoerMuldersTMDPDF_subGridsB_inKT(:),BoerMuldersTMDPDF_subGridsKT_inKT(:)
integer::BoerMuldersTMDPDF_grid_SizeX_inKT,BoerMuldersTMDPDF_grid_SizeB_inKT,BoerMuldersTMDPDF_grid_SizeKT_inKT
real(dp)::BoerMuldersTMDPDF_minQ_inKT,BoerMuldersTMDPDF_maxQ_inKT
integer::BoerMuldersTMDPDF_grid_SizeQ_inKT
real(dp)::BoerMuldersTMDPDF_hOGATA_TMM,BoerMuldersTMDPDF_toleranceOGATA_TMM,BoerMuldersTMDPDF_muMIN_TMM

!-------------------- TMDF parameters
logical::include_TMDF
real(dp)::TMDF_OGATAh,TMDF_tolerance,TMDF_qTMIN
real(dp)::TMDF_mass

!-------------------- TMDF-KPC-DY parameters
logical::include_TMDF_KPC
real(dp)::TMDF_KPC_toleranceGEN,TMDF_KPC_toleranceINT,TMDF_KPC_qTMIN

!-------------------- TMDX-DY parameters
logical::include_TMDX_DY
character*8::TMDX_DY_order
real(dp)::TMDX_DY_toleranceINT, TMDX_DY_toleranceGEN
integer::TMDX_DY_ptSECTION
logical::TMDX_DY_exactX1X2,TMDX_DY_piResum,TMDX_DY_exactScale
real::TMDX_DY_maxQbinSize,TMDX_DY_minqTabs
logical::TMDX_DY_useKPC

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
public::artemide_Setup_fromFile,CreateConstantsFile,ReadConstantsFile,CheckConstantsFile

contains

  
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

    number_of_processors=8

    include_EWinput=.true.
    !-----------------------physic parameters
    mCHARM=1.2700d0    !threashold mass for charm quark
    mBOTTOM=4.180d0    !threashold mass for bottom quark
    mTOP=172.90d0    !threashold mass for top quark

    mZ=91.1876d0    !pole mass for Z-boson
    GammaZ=2.4952d0    !width of Z-boson

    mW=80.379d0        !pole mass for W-boson
    GammaW=2.085d0    !width of W-boson

    mH=125.1d0        !pole mass for Higgs-boson
    GammaH=0.0042d0    !width of Higgs-boson [GeV]
    vevH=246.d0        !Higgs Vev

    mELECTRON=0.0005109989461d0    !mass of electron[GeV]
    mMUON=0.1056583745d0    !mass of muon[GeV]
    mTAU=1.77686d0        !mass of tau

    hc2=0.389379338d0    !transformation constant (hc)^2   GeV->mbarn
    Mnorm=1.d0    !Mass parameter used for the normalization of TMDs

    alphaQED_MZ=127.955d0    !inverse alpha_QED at Z-boson mass
    alphaQED_MTAU=133.476d0    !inverse alpha_QED at TAU-lepton mass

    sW2=0.23122d0    !sin^2 theta_Weinberg

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
    TMDR_override=.false.
    TMDR_orderGAMMA=trim(order)
    TMDR_orderV=trim(order)
    TMDR_orderD=trim(order)
    TMDR_orderZETA=trim(order)
    TMDR_lambdaLength=2
    TMDR_tolerance=0.000001d0    !tolerance of searching for saddle point, numerical integration etc.
    TMDR_bFREEZE=0.000001d0
    TMDR_smooth=0.01d0

    !-------------------- parameters for UTMDPDF
    include_uTMDPDF=.true.
    uTMDPDF_withGluon=.false.
    uTMDPDF_numHadron=1
    uTMDPDF_order=trim("NLO")
    uTMDPDF_makeGrid=.true.
    uTMDPDF_runGridTest=.false.
    uTMDPDF_lambdaLength=2
    uTMDPDF_BMAX_ABS=100.d0
    uTMDPDF_toleranceINT=1.d-6!tolerance (i.e. relative integration tolerance)
    uTMDPDF_toleranceGEN=1.d-6!general tolerance
    uTMDPDF_maxIteration=10000    !maxIteration for adaptive integration
    uTMDPDF_numSubGridsX=4
    allocate(uTMDPDF_subGridsX(0:uTMDPDF_numSubGridsX))
    uTMDPDF_subGridsX=(/0.00001d0,0.001d0,0.1d0,0.7d0,1.d0/)
    uTMDPDF_grid_SizeX=8
    uTMDPDF_numSubGridsB=4
    allocate(uTMDPDF_subGridsB(0:uTMDPDF_numSubGridsB))
    uTMDPDF_subGridsB=(/0.00001d0,0.01d0,0.2d0,2.d0,8.d0,25.d0/)
    uTMDPDF_grid_SizeB=8
    uTMDPDF_makeGrid_inKT=.false.
    uTMDPDF_runGridTest_inKT=.false.
    uTMDPDF_numSubGridsX_inKT=4
    allocate(uTMDPDF_subGridsX_inKT(0:uTMDPDF_numSubGridsX_inKT))
    uTMDPDF_subGridsX_inKT=(/0.00001d0,0.001d0,0.1d0,0.7d0,1.d0/)
    uTMDPDF_grid_SizeX_inKT=16
    uTMDPDF_numSubGridsKT_inKT=5
    allocate(uTMDPDF_subGridsKT_inKT(0:uTMDPDF_numSubGridsKT_inKT))
    uTMDPDF_subGridsKT_inKT=(/0.01d0,1.d0,5.d0,15.d0,50.d0,200.d0/)
    uTMDPDF_grid_SizeKT_inKT=16
    uTMDPDF_minQ_inKT=1.d0
    uTMDPDF_maxQ_inKT=200.d0
    uTMDPDF_grid_SizeQ_inKT=40
    uTMDPDF_numSubGridsB_inKT=5
    allocate(uTMDPDF_subGridsB_inKT(0:uTMDPDF_numSubGridsB_inKT))
    uTMDPDF_subGridsB_inKT=(/0.00001d0,0.01d0,0.2d0,2.d0,8.d0,25.d0/)
    uTMDPDF_grid_SizeB_inKT=16
    uTMDPDF_toleranceOGATA_TMM=1.d-4    !!! OGATA tolerance (for TMM)
    uTMDPDF_hOGATA_TMM=1.d-3            !!! OGATA integration step(for TMM)
    uTMDPDF_muMIN_TMM=0.8d0         !!! min value of mu for TMM


    !-------------------- parameters for UTMDFF
    include_uTMDFF=.false.!!! we do not initialize TMDFF by definition
    uTMDFF_order=trim("NLO")
    uTMDFF_withGluon=.false.
    uTMDFF_numHadron=0
    uTMDFF_makeGrid=.true.
    uTMDFF_runGridTest=.false.
    uTMDFF_lambdaLength=2
    uTMDFF_BMAX_ABS=100.d0
    uTMDFF_toleranceINT=1.d-6!tolerance (i.e. relative integration tolerance)
    uTMDFF_toleranceGEN=1.d-6!general tolerance
    uTMDFF_maxIteration=10000    !maxIteration for adaptive integration
    uTMDFF_numSubGridsX=3
    allocate(uTMDFF_subGridsX(0:uTMDFF_numSubGridsX))
    uTMDFF_subGridsX=(/0.001d0,0.1d0,0.7d0,1.d0/)
    uTMDFF_grid_SizeX=8
    uTMDFF_numSubGridsB=4
    allocate(uTMDFF_subGridsB(0:uTMDFF_numSubGridsB))
    uTMDFF_subGridsB=(/0.00001d0,0.01d0,0.2d0,2.d0,25.d0/)
    uTMDFF_grid_SizeB=8
    uTMDFF_toleranceOGATA=1.d-4    !!! OGATA tolerance
    uTMDFF_hOGATA=1.d-3            !!! OGATA integration step
    uTMDFF_KT_FREEZE=1.d-4         !!! min value of kT
    uTMDFF_toleranceOGATA_TMM=1.d-4    !!! OGATA tolerance (for TMM)
    uTMDFF_hOGATA_TMM=1.d-3            !!! OGATA integration step(for TMM)
    uTMDFF_muMIN_TMM=0.8d0         !!! min value of mu for TMM

    !-------------------- parameters for lpTMDPDF
    include_lpTMDPDF=.false.!!! we do not initialize TMDFF by definition
    lpTMDPDF_order=trim("NLO")
    lpTMDPDF_withGluon=.true. !!! this is true by default
    lpTMDPDF_makeGrid=.true.
    lpTMDPDF_numHadron=0
    lpTMDPDF_runGridTest=.false.
    lpTMDPDF_lambdaLength=2
    lpTMDPDF_BMAX_ABS=100.d0
    lpTMDPDF_toleranceINT=1.d-6!tolerance (i.e. relative integration tolerance)
    lpTMDPDF_toleranceGEN=1.d-6!general tolerance
    lpTMDPDF_maxIteration=10000    !maxIteration for adaptive integration
    lpTMDPDF_numSubGridsX=4
    allocate(lpTMDPDF_subGridsX(0:lpTMDPDF_numSubGridsX))
    lpTMDPDF_subGridsX=(/0.00001d0,0.001d0,0.1d0,0.7d0,1.d0/)
    lpTMDPDF_grid_SizeX=8
    lpTMDPDF_numSubGridsB=4
    allocate(lpTMDPDF_subGridsB(0:lpTMDPDF_numSubGridsB))
    lpTMDPDF_subGridsB=(/0.00001d0,0.01d0,0.2d0,2.d0,25.d0/)
    lpTMDPDF_grid_SizeB=8
    lpTMDPDF_toleranceOGATA=1.d-4    !!! OGATA tolerance
    lpTMDPDF_hOGATA=1.d-3            !!! OGATA integration step
    lpTMDPDF_KT_FREEZE=1.d-4         !!! min value of kT
    lpTMDPDF_toleranceOGATA_TMM=1.d-4    !!! OGATA tolerance (for TMM)
    lpTMDPDF_hOGATA_TMM=1.d-3            !!! OGATA integration step(for TMM)
    lpTMDPDF_muMIN_TMM=0.8d0         !!! min value of mu for TMM

    !-------------------- parameters for SiversTMDPDF
    include_SiversTMDPDF=.false. !!! we do not initialize SiversTMDPDF by definition
    SiversTMDPDF_order=trim('NA') !!! by definition Sivers is tree-order
    SiversTMDPDF_makeGrid=.false.   !!! no need to make grid
    SiversTMDPDF_withGluon=.false.
    SiversTMDPDF_numHadron=1
    SiversTMDPDF_runGridTest=.false.
    SiversTMDPDF_lambdaLength=1
    SiversTMDPDF_BMAX_ABS=100.d0
    SiversTMDPDF_toleranceINT=1.d-6!tolerance (i.e. relative integration tolerance)
    SiversTMDPDF_toleranceGEN=1.d-6!general tolerance
    SiversTMDPDF_maxIteration=10000    !maxIteration for adaptive integration
    SiversTMDPDF_toleranceOGATA=1.d-4    !!! OGATA tolerance
    SiversTMDPDF_hOGATA=1.d-3            !!! OGATA integration step
    SiversTMDPDF_KT_FREEZE=1.d-4         !!! min value of kT
    SiversTMDPDF_toleranceOGATA_TMM=1.d-4    !!! OGATA tolerance (for TMM)
    SiversTMDPDF_hOGATA_TMM=1.d-3            !!! OGATA integration step(for TMM)
    SiversTMDPDF_muMIN_TMM=0.8d0         !!! min value of mu for TMM
    
    !-------------------- parameters for wgtTMDPDF
    include_wgtTMDPDF=.false.
    wgtTMDPDF_order=trim("LO")
    wgtTMDPDF_makeGrid=.true.
    wgtTMDPDF_withGluon=.false.
    wgtTMDPDF_numHadron=1
    wgtTMDPDF_runGridTest=.false.
    wgtTMDPDF_lambdaLength=2
    wgtTMDPDF_BMAX_ABS=100.d0
    wgtTMDPDF_toleranceINT=1.d-6!tolerance (i.e. relative integration tolerance)
    wgtTMDPDF_toleranceGEN=1.d-6!general tolerance
    wgtTMDPDF_maxIteration=10000    !maxIteration for adaptive integration
    wgtTMDPDF_numSubGridsX=3
    allocate(wgtTMDPDF_subGridsX(0:wgtTMDPDF_numSubGridsX))
    wgtTMDPDF_subGridsX=(/0.001d0,0.1d0,0.7d0,1.d0/)
    wgtTMDPDF_grid_SizeX=8
    wgtTMDPDF_numSubGridsB=4
    allocate(wgtTMDPDF_subGridsB(0:wgtTMDPDF_numSubGridsB))
    wgtTMDPDF_subGridsB=(/0.00001d0,0.01d0,0.2d0,2.d0,25.d0/)
    wgtTMDPDF_grid_SizeB=8
    wgtTMDPDF_order_tw3=trim("NA")
    wgtTMDPDF_makeGrid_tw3=.false.
    wgtTMDPDF_withGluon_tw3=.false. !!! this is true by default
    wgtTMDPDF_runGridTest_tw3=.false.
    wgtTMDPDF_toleranceOGATA=1.d-4    !!! OGATA tolerance
    wgtTMDPDF_hOGATA=1.d-3            !!! OGATA integration step
    wgtTMDPDF_KT_FREEZE=1.d-4         !!! min value of kT
    wgtTMDPDF_toleranceOGATA_TMM=1.d-4    !!! OGATA tolerance (for TMM)
    wgtTMDPDF_hOGATA_TMM=1.d-3            !!! OGATA integration step(for TMM)
    wgtTMDPDF_muMIN_TMM=0.8d0         !!! min value of mu for TMM
    
    !-------------------- parameters for BoerMuldersTMDPDF
    include_BoerMuldersTMDPDF=.false. !!! we do not initialize BoerMuldersTMDPDF by definition
    BoerMuldersTMDPDF_order=trim('NA') !!! by definition BoerMulders is tree-order
    BoerMuldersTMDPDF_makeGrid=.false.   !!! no need to make grid
    BoerMuldersTMDPDF_withGluon=.false.
    BoerMuldersTMDPDF_numHadron=1
    BoerMuldersTMDPDF_runGridTest=.false.
    BoerMuldersTMDPDF_lambdaLength=1
    BoerMuldersTMDPDF_BMAX_ABS=100.d0
    BoerMuldersTMDPDF_toleranceINT=1.d-6!tolerance (i.e. relative integration tolerance)
    BoerMuldersTMDPDF_toleranceGEN=1.d-6!general tolerance
    BoerMuldersTMDPDF_maxIteration=10000    !maxIteration for adaptive integration
    BoerMuldersTMDPDF_makeGrid_inKT=.false.
    BoerMuldersTMDPDF_runGridTest_inKT=.false.
    BoerMuldersTMDPDF_numSubGridsX_inKT=4
    allocate(BoerMuldersTMDPDF_subGridsX_inKT(0:BoerMuldersTMDPDF_numSubGridsX_inKT))
    BoerMuldersTMDPDF_subGridsX_inKT=(/0.00001d0,0.001d0,0.1d0,0.7d0,1.d0/)
    BoerMuldersTMDPDF_grid_SizeX_inKT=16
    BoerMuldersTMDPDF_numSubGridsKT_inKT=5
    allocate(BoerMuldersTMDPDF_subGridsKT_inKT(0:BoerMuldersTMDPDF_numSubGridsKT_inKT))
    BoerMuldersTMDPDF_subGridsKT_inKT=(/0.01d0,1.d0,5.d0,15.d0,50.d0,200.d0/)
    BoerMuldersTMDPDF_grid_SizeKT_inKT=16
    BoerMuldersTMDPDF_minQ_inKT=1.d0
    BoerMuldersTMDPDF_maxQ_inKT=200.d0
    BoerMuldersTMDPDF_grid_SizeQ_inKT=40
    BoerMuldersTMDPDF_numSubGridsB_inKT=5
    allocate(BoerMuldersTMDPDF_subGridsB_inKT(0:BoerMuldersTMDPDF_numSubGridsB_inKT))
    BoerMuldersTMDPDF_subGridsB_inKT=(/0.00001d0,0.01d0,0.2d0,2.d0,8.d0,25.d0/)
    BoerMuldersTMDPDF_grid_SizeB_inKT=16
    BoerMuldersTMDPDF_toleranceOGATA_TMM=1.d-4    !!! OGATA tolerance (for TMM)
    BoerMuldersTMDPDF_hOGATA_TMM=1.d-3            !!! OGATA integration step(for TMM)
    BoerMuldersTMDPDF_muMIN_TMM=0.8d0         !!! min value of mu for TMM

    !------------------ parameters for TMDF
    include_TMDF=.true.
    TMDF_tolerance=0.0001d0    !tolerance (i.e. relative integration tolerance)
    TMDF_OGATAh=0.001d0        !Ogata quadrature integration step
    TMDF_qTMIN=0.001d0          !minimal qT, below values is contant
    TMDF_mass=0.938272      !mass parameter that is used as reference dimension

    !------------------ parameters for TMDX-DY
    include_TMDX_DY=.true.
    TMDX_DY_toleranceINT=0.0001d0    !tolerance for integration (i.e. relative integration tolerance -- in kinematic integrals;)
    TMDX_DY_toleranceGEN=0.000001d0    !tolerance general (i.e. for comparison of variables)
    TMDX_DY_ptSECTION=4        !default number of sections for pt-bin integration
    TMDX_DY_maxQbinSize=30. !default maximum size of the Q-bin integration
    TMDX_DY_minqTabs=0.0001 !default minimum value of qT
    TMDX_DY_order=trim(order)
    TMDX_DY_exactX1X2=.true.
    TMDX_DY_piResum=.false.
    TMDX_DY_exactScale=.false.
    TMDX_DY_useKPC=.false.

    !------------------ parameters for TMDX-SIDIS
    include_TMDX_SIDIS=.false.
    TMDX_SIDIS_tolerance=0.001d0    !tolerance (i.e. relative integration tolerance -- in kinematic integrals;)
    TMDX_SIDIS_ptSECTION=4        !default number of sections for pt-bin integration
    TMDX_SIDIS_order=trim(order)
    TMDX_SIDIS_qTcorr=.true.
    TMDX_SIDIS_M1corr=.true.
    TMDX_SIDIS_M2corr=.true.
    TMDX_SIDIS_qTinX1Z1corr=.true.
    TMDX_SIDIS_exactScale=.false.
    TMDX_SIDIS_numProc=8
    TMDX_SIDIS_toleranceZ=TMDX_SIDIS_tolerance
    TMDX_SIDIS_methodZ='SA'        !SA=Simpson adaptive, S5=Simpson 5-point
    TMDX_SIDIS_toleranceX=TMDX_SIDIS_tolerance
    TMDX_SIDIS_methodX='SA'        !SA=Simpson adaptive, S5=Simpson 5-point

    !------------------ parameters for TMDF-KPC-DY
    include_TMDF_KPC=.true.
    TMDF_KPC_toleranceGEN=0.000001d0  !tolerance general (i.e. comparison etc)
    TMDF_KPC_toleranceINT=0.0001d0    !tolerance integration (i.e. relative integration tolerance)
    TMDF_KPC_qTMIN=0.001d0            !minimal value of qT (below the value is frozen)

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

    if(FILEversion<30) then
        write(*,*) color('aTMDe_setup: present version of setup-file is for artemide v2.'&
        , c_red_bold)
        write(*,*) color('..           It is incompatible with the preent version.'&
        , c_red)
        write(*,*) color('..           Please, use actual file or update it manually.'&
        ,c_red)
        stop
    end if

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
    write(51,"('*p2  : Global mass for normalization of TMDs [GeV]')")
    write(51,*) Mnorm
    write(51,"(' ')")
    write(51,"('*C   : ---- aTMDe-control parameters ----')")
    write(51,"('*p1  : Maximum number of processors to use (used only if compiled with openMP)')")
    write(51,*) number_of_processors
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
    write(51,&
    "('*p2  : Override the general definition of the pertrubative order (in this case orders are defined by the list below))')")
    write(51,*) TMDR_override
    write(51,"('*p3  : Order of Gamma-cusp (LO=1-loop, NLO=2-loop, ...)')")
    write(51,*) trim(TMDR_orderGAMMA)
    write(51,"('*p4  : Order of Gamma-V (LO=0-loop=0, NLO=1-loop, ...)n')")
    write(51,*) trim(TMDR_orderV)
    write(51,"('*p5  : Order of perturbative RAD (LO=0-loop=0, NLO=1-loop, ...)')")
    write(51,*) trim(TMDR_orderD)
    write(51,"('*p6  : Order of (exact) zeta-line (LO=a^1.., NLO=a^2..., ...)')")
    write(51,*) trim(TMDR_orderZETA)
    write(51,"(' ')")
    write(51,"('*B   : ---- Parameters of NP model ----')")
    write(51,"('*p1  : Length of lambdaNP')")
    write(51,*) TMDR_lambdaLength
    write(51,"(' ')")
    write(51,*)
    write(51,"('*C   : ---- Numerical evaluation parameters ----')")
    write(51,"('*p1  : Tolerance (tolerance used for comparisons)')")
    write(51,*) TMDR_tolerance
    write(51,"('*p2  : The smallest value of b (for smaller OPE is frozen)')")
    write(51,*) TMDR_bFREEZE
    write(51,"('*p3  : The smoothing parameters for small-b stabilization')")
    write(51,*) TMDR_smooth


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
    write(51,"('*p1  : Include gluon TMDPDFs')")
    write(51,*) uTMDPDF_withGluon
    write(51,"('*p2  : Number of hadrons (in order starting with 1)')")
    write(51,*) uTMDPDF_numHadron
    write(51,"(' ')")
    write(51,"('*B   : ---- OPE main definitions ----')")
    write(51,"('*p1  : Order of coefficient function')")
    write(51,*) trim(uTMDPDF_order)
    write(51,"('*p2  : Prepare grid')")
    write(51,*) uTMDPDF_makeGrid
    write(51,"('*p3  : run the test of the grid (takes some time)')")
    write(51,*) uTMDPDF_runGridTest
    write(51,"(' ')")
    write(51,"('*C   : ---- Parameters of NP model ----')")
    write(51,"('*p1  : Length of lambdaNP')")
    write(51,*) uTMDPDF_lambdaLength
    write(51,"('*p2  : Absolute maximum b (for larger b, TMD=0)')")
    write(51,*) uTMDPDF_BMAX_ABS
    write(51,"(' ')")
    write(51,"('*D   : ---- Numerical evaluation parameters ----')")
    write(51,"('*p1  : Tolerance (relative tolerance of convolution integral)')")
    write(51,*) uTMDPDF_toleranceINT
    write(51,"('*p2  : Tolerance general (used for various comparisons)')")
    write(51,*) uTMDPDF_toleranceGEN
    write(51,"('*p3  : Maximum number of iterations (for adaptive integration)')")
    write(51,*) uTMDPDF_maxIteration
    write(51,"(' ')")
    write(51,"('*E   : ---- (OPE) Parameters of grid ----')")
    write(51,"('*p1  : Number of subgrids in X (required to read the next line)')")
    write(51,*) uTMDPDF_numSubGridsX
    write(51,"('*p2  : Intervals for subgrids in X (must include 1., as the last point)')")
    write(51,*) uTMDPDF_subGridsX
    write(51,"('*p3  : Number of nodes in the X-subgrid')")
    write(51,*) uTMDPDF_grid_SizeX
    write(51,"('*p4  : Number of subgrids in B (required to read the next line)')")
    write(51,*) uTMDPDF_numSubGridsB
    write(51,"('*p5  : Intervals for subgrids in B (below and above ultimate points the value is frozen)')")
    write(51,*) uTMDPDF_subGridsB
    write(51,"('*p6  : Number of nodes in the B-subgrid ')")
    write(51,*) uTMDPDF_grid_SizeB
    write(51,"(' ')")
    write(51,"('*F   : ---- Transform and grid in KT-space ----')")
    write(51,"('*p1  : Prepare grid')")
    write(51,*) uTMDPDF_makeGrid_inKT
    write(51,"('*p2  : run the test of the grid (takes some time)')")
    write(51,*) uTMDPDF_runGridTest_inKT
    write(51,"('*p3  : Number of subgrids in X (required to read the next line)')")
    write(51,*) uTMDPDF_numSubGridsX_inKT
    write(51,"('*p4  : Intervals for subgrids in X (must include 1., as the last point)')")
    write(51,*) uTMDPDF_subGridsX_inKT
    write(51,"('*p5  : Number of nodes in the X-subgrid')")
    write(51,*) uTMDPDF_grid_SizeX_inKT
    write(51,"('*p6  : Number of subgrids in K (required to read the next line)')")
    write(51,*) uTMDPDF_numSubGridsKT_inKT
    write(51,"('*p7  : Intervals for subgrids in KT (below and above ultimate points the value is frozen)')")
    write(51,*) uTMDPDF_subGridsKT_inKT
    write(51,"('*p8  : Number of nodes in the KT-subgrid')")
    write(51,*) uTMDPDF_grid_SizeKT_inKT
    write(51,"('*p9  : Minimal Q in the grid')")
    write(51,*) uTMDPDF_minQ_inKT
    write(51,"('*p10 : Maximal Q in the grid')")
    write(51,*) uTMDPDF_maxQ_inKT
    write(51,"('*p11 : Number of nodes in the Q-grid')")
    write(51,*) uTMDPDF_grid_SizeQ_inKT
    write(51,"('*p12  : Number of subgrids in B (required to read the next line)')")
    write(51,*) uTMDPDF_numSubGridsB_inKT
    write(51,"('*p13  : Intervals for subgrids in B (below and above ultimate points the value is frozen)')")
    write(51,*) uTMDPDF_subGridsB_inKT
    write(51,"('*p14  : Number of nodes in the B-subgrid')")
    write(51,*) uTMDPDF_grid_SizeB_inKT
    write(51,"(' ')")
    write(51,"('*G   : ---- Computation of Transverse Momentum Moments (TMM) ----')")
    write(51,"('*p1  : Tolerance (relative tolerance of summation in OGATA quadrature)')")
    write(51,*) uTMDPDF_toleranceOGATA_TMM
    write(51,"('*p2  : Ogata quadrature integration step')")
    write(51,*) uTMDPDF_hOGATA_TMM
    write(51,"('*p3  : Minimum value of mu [GeV] (below that value the computation is terminated)')")
    write(51,*) uTMDPDF_muMIN_TMM


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
    write(51,"('*p1  : Include gluon TMDPDFs')")
    write(51,*) uTMDFF_withGluon
    write(51,"('*p2  : Number of hadrons (in order starting with 1)')")
    write(51,*) uTMDFF_numHadron
    write(51,"(' ')")
    write(51,"('*B   : ---- OPE main definitions ----')")
    write(51,"('*p1  : Order of coefficient function')")
    write(51,*) trim(uTMDFF_order)
    write(51,"('*p2  : Prepare grid')")
    write(51,*) uTMDFF_makeGrid
    write(51,"('*p3  : run the test of the grid (takes some time)')")
    write(51,*) uTMDFF_runGridTest
    write(51,"(' ')")
    write(51,"('*C   : ---- Parameters of NP model ----')")
    write(51,"('*p1  : Length of lambdaNP')")
    write(51,*) uTMDFF_lambdaLength
    write(51,"('*p2  : Absolute maximum b (for larger b, TMD=0)')")
    write(51,*) uTMDFF_BMAX_ABS
    write(51,"(' ')")
    write(51,"('*D   : ---- Numerical evaluation parameters ----')")
    write(51,"('*p1  : Tolerance (relative tolerance of convolution integral)')")
    write(51,*) uTMDFF_toleranceINT
    write(51,"('*p2  : Tolerance general (used for various comparisons)')")
    write(51,*) uTMDFF_toleranceGEN
    write(51,"('*p3  : Maximum number of iterations (for adaptive integration)')")
    write(51,*) uTMDFF_maxIteration
    write(51,"(' ')")
    write(51,"('*E   : ---- (OPE) Parameters of grid ----')")
    write(51,"('*p1  : Number of subgrids in X (required to read the next line)')")
    write(51,*) uTMDFF_numSubGridsX
    write(51,"('*p2  : Intervals for subgrids in X (must include 1., as the last point)')")
    write(51,*) uTMDFF_subGridsX
    write(51,"('*p3  : Number of nodes in the X-subgrid')")
    write(51,*) uTMDFF_grid_SizeX
    write(51,"('*p4  : Number of subgrids in B (required to read the next line)')")
    write(51,*) uTMDFF_numSubGridsB
    write(51,"('*p5  : Intervals for subgrids in B (below and above ultimate points the value is frozen)')")
    write(51,*) uTMDFF_subGridsB
    write(51,"('*p6  : Number of nodes in the B-subgrid ')")
    write(51,*) uTMDFF_grid_SizeB
    write(51,"(' ')")
    write(51,"('*F   : ---- Transformation to KT-space ----')")
    write(51,"('*p1  : Tolerance (relative tolerance of summation in OGATA quadrature)')")
    write(51,*) uTMDFF_toleranceOGATA
    write(51,"('*p2  : Ogata quadrature integration step')")
    write(51,*) uTMDFF_hOGATA
    write(51,"('*p3  : Minimum value of kT (below that value function is constant)')")
    write(51,*) uTMDFF_KT_FREEZE
    write(51,"(' ')")
    write(51,"('*G   : ---- Computation of Transverse Momentum Moments (TMM) ----')")
    write(51,"('*p1  : Tolerance (relative tolerance of summation in OGATA quadrature)')")
    write(51,*) uTMDFF_toleranceOGATA_TMM
    write(51,"('*p2  : Ogata quadrature integration step')")
    write(51,*) uTMDFF_hOGATA_TMM
    write(51,"('*p3  : Minimum value of mu [GeV] (below that value the computation is terminated)')")
    write(51,*) uTMDFF_muMIN_TMM

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
    write(51,"('*p3  : Minimal qT (for smaller values the expression is constant)')")
    write(51,*) TMDF_qTMIN
    write(51,"('*B   : ---- Global garameters of structure functions----')")
    write(51,"('*p1  : Mass parameter used in the structure function (mass of hadron)')")
    write(51,*) TMDF_mass

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
    write(51,"('*p2  : Utilize KPC formula (T=KPC formula, F=LP formula)')")
    write(51,*) TMDX_DY_useKPC
    write(51,"('*p3  : Use resummation of pi^2-corrections in hard coefficient')")
    write(51,*) TMDX_DY_piResum
    write(51,"('*B   : ---- Numerical evaluation parameters ----')")
    write(51,"('*p1  : Tolerance general (variable comparision, etc.)')")
    write(51,*) TMDX_DY_toleranceGEN
    write(51,"('*p2  : Tolerance (relative tolerance for bin-integration routines, except pt-integration)')")
    write(51,*) TMDX_DY_toleranceINT
    write(51,"('*p3  : Minimal number of sections for pt-integration')")
    write(51,*) TMDX_DY_ptSECTION
    write(51,"('*p4  : Maximum size of the Q-bin integration (larger bins are desected)')")
    write(51,*) TMDX_DY_maxQbinSize
    write(51,"('*p5  : Minimal value of qT (lower values are fixed to this number)')")
    write(51,*) TMDX_DY_minqTabs
    write(51,"('*C   : ---- Definition of LP TMD factorization ----')")
    write(51,"('*p1  : Use the exact values of x1 and x2 (include qT/Q correction)')")
    write(51,*) TMDX_DY_exactX1X2
    write(51,"('*p2  : Use the exact values for factorization scales (include qT/Q correction)')")
    write(51,*) TMDX_DY_exactScale
    write(51,"('*D   : ---- Definition of TMD factorization with KPC ----')")
    write(51,"('*p1  : Include terms induced by fiducial cuts')")
    write(51,*) .true.



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
    write(51,"('*p1  : Include gluon TMDPDFs [=T, because it is GLUON TMDPDF]')")
    write(51,*) lpTMDPDF_withGluon
    write(51,"('*p2  : Number of hadrons (in order starting with 1)')")
    write(51,*) lpTMDPDF_numHadron
    write(51,"(' ')")
    write(51,"('*B   : ---- OPE main definitions ----')")
    write(51,"('*p1  : Order of coefficient function')")
    write(51,*) trim(lpTMDPDF_order)
    write(51,"('*p2  : Prepare grid')")
    write(51,*) lpTMDPDF_makeGrid
    write(51,"('*p3  : run the test of the grid (takes some time)')")
    write(51,*) lpTMDPDF_runGridTest
    write(51,"(' ')")
    write(51,"('*C   : ---- Parameters of NP model ----')")
    write(51,"('*p1  : Length of lambdaNP')")
    write(51,*) lpTMDPDF_lambdaLength
    write(51,"('*p2  : Absolute maximum b (for larger b, TMD=0)')")
    write(51,*) lpTMDPDF_BMAX_ABS
    write(51,"(' ')")
    write(51,"('*D   : ---- Numerical evaluation parameters ----')")
    write(51,"('*p1  : Tolerance (relative tolerance of convolution integral)')")
    write(51,*) lpTMDPDF_toleranceINT
    write(51,"('*p2  : Tolerance general (used for various comparisons)')")
    write(51,*) lpTMDPDF_toleranceGEN
    write(51,"('*p3  : Maximum number of iterations (for adaptive integration)')")
    write(51,*) lpTMDPDF_maxIteration
    write(51,"(' ')")
    write(51,"('*E   : ---- (OPE) Parameters of grid ----')")
    write(51,"('*p1  : Number of subgrids in X (required to read the next line)')")
    write(51,*) lpTMDPDF_numSubGridsX
    write(51,"('*p2  : Intervals for subgrids in X (must include 1., as the last point)')")
    write(51,*) lpTMDPDF_subGridsX
    write(51,"('*p3  : Number of nodes in the X-subgrid')")
    write(51,*) lpTMDPDF_grid_SizeX
    write(51,"('*p4  : Number of subgrids in B (required to read the next line)')")
    write(51,*) lpTMDPDF_numSubGridsB
    write(51,"('*p5  : Intervals for subgrids in B (below and above ultimate points the value is frozen)')")
    write(51,*) lpTMDPDF_subGridsB
    write(51,"('*p6  : Number of nodes in the B-subgrid ')")
    write(51,*) lpTMDPDF_grid_SizeB
    write(51,"(' ')")
    write(51,"('*F   : ---- Transformation to KT-space ----')")
    write(51,"('*p1  : Tolerance (relative tolerance of summation in OGATA quadrature)')")
    write(51,*) lpTMDPDF_toleranceOGATA
    write(51,"('*p2  : Ogata quadrature integration step')")
    write(51,*) lpTMDPDF_hOGATA
    write(51,"('*p3  : Minimum value of kT (below that value function is constant)')")
    write(51,*) lpTMDPDF_KT_FREEZE
    write(51,"(' ')")
    write(51,"('*G   : ---- Computation of Transverse Momentum Moments (TMM) ----')")
    write(51,"('*p1  : Tolerance (relative tolerance of summation in OGATA quadrature)')")
    write(51,*) lpTMDPDF_toleranceOGATA_TMM
    write(51,"('*p2  : Ogata quadrature integration step')")
    write(51,*) lpTMDPDF_hOGATA_TMM
    write(51,"('*p3  : Minimum value of mu [GeV] (below that value the computation is terminated)')")
    write(51,*) lpTMDPDF_muMIN_TMM
    
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
    write(51,"('*p1  : Include gluon TMDPDFs')")
    write(51,*) SiversTMDPDF_withGluon
    write(51,"('*p2  : Number of hadrons (in order starting with 1)')")
    write(51,*) SiversTMDPDF_numHadron
    write(51,"(' ')")
    write(51,"('*B   : ---- OPE[tw3] main definitions ----')")
    write(51,"('*p1  : Order of coefficient function')")
    write(51,*) trim(SiversTMDPDF_order)
    write(51,"('*p2  : Prepare grid')")
    write(51,*) SiversTMDPDF_makeGrid
    write(51,"('*p3  : run the test of the grid (takes some time)')")
    write(51,*) SiversTMDPDF_runGridTest
    write(51,"(' ')")
    write(51,"('*C   : ---- Parameters of NP model ----')")
    write(51,"('*p1  : Length of lambdaNP')")
    write(51,*) SiversTMDPDF_lambdaLength
    write(51,"('*p2  : Absolute maximum b (for larger b, TMD=0)')")
    write(51,*) SiversTMDPDF_BMAX_ABS
    write(51,"(' ')")
    write(51,"('*D   : ---- Numerical evaluation parameters ----')")
    write(51,"('*p1  : Tolerance (relative tolerance of convolution integral)')")
    write(51,*) SiversTMDPDF_toleranceINT
    write(51,"('*p2  : Tolerance general (used for various comparisons)')")
    write(51,*) SiversTMDPDF_toleranceGEN
    write(51,"('*p3  : Maximum number of iterations (for adaptive integration)')")
    write(51,*) SiversTMDPDF_maxIteration
    write(51,"(' ')")
    write(51,"('*E   : ---- (OPE) Parameters of grid ----')")
    write(51,"(' ')")
    write(51,"('*F   : ---- Transformation to KT-space ----')")
    write(51,"('*p1  : Tolerance (relative tolerance of summation in OGATA quadrature)')")
    write(51,*) SiversTMDPDF_toleranceOGATA
    write(51,"('*p2  : Ogata quadrature integration step')")
    write(51,*) SiversTMDPDF_hOGATA
    write(51,"('*p3  : Minimum value of kT (below that value function is constant)')")
    write(51,*) SiversTMDPDF_KT_FREEZE
    write(51,"(' ')")
    write(51,"('*G   : ---- Computation of Transverse Momentum Moments (TMM) ----')")
    write(51,"('*p1  : Tolerance (relative tolerance of summation in OGATA quadrature)')")
    write(51,*) SiversTMDPDF_toleranceOGATA_TMM
    write(51,"('*p2  : Ogata quadrature integration step')")
    write(51,*) SiversTMDPDF_hOGATA_TMM
    write(51,"('*p3  : Minimum value of mu [GeV] (below that value the computation is terminated)')")
    write(51,*) SiversTMDPDF_muMIN_TMM
    
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
    write(51,"('*p1  : Include gluon TMDPDFs')")
    write(51,*) wgtTMDPDF_withGluon
    write(51,"('*p2  : Number of hadrons (in order starting with 1)')")
    write(51,*) wgtTMDPDF_numHadron
    write(51,"(' ')")
    write(51,"('*B   : ---- OPE[tw2] main definitions ----')")
    write(51,"('*p1  : Order of coefficient function')")
    write(51,*) trim(wgtTMDPDF_order)
    write(51,"('*p2  : Prepare grid')")
    write(51,*) wgtTMDPDF_makeGrid
    write(51,"('*p3  : run the test of the grid (takes some time)')")
    write(51,*) wgtTMDPDF_runGridTest
    write(51,"(' ')")
    write(51,"('*C   : ---- Parameters of NP model ----')")
    write(51,"('*p1  : Length of lambdaNP')")
    write(51,*) wgtTMDPDF_lambdaLength
    write(51,"('*p2  : Absolute maximum b (for larger b, TMD=0)')")
    write(51,*) wgtTMDPDF_BMAX_ABS
    write(51,"(' ')")
    write(51,"('*D   : ---- Numerical evaluation parameters ----')")
    write(51,"('*p1  : Tolerance (relative tolerance of convolution integral)')")
    write(51,*) wgtTMDPDF_toleranceINT
    write(51,"('*p2  : Tolerance general (used for various comparisons)')")
    write(51,*) wgtTMDPDF_toleranceGEN
    write(51,"('*p3  : Maximum number of iterations (for adaptive integration)')")
    write(51,*) wgtTMDPDF_maxIteration
    write(51,"(' ')")
    write(51,"('*E   : ---- OPE[tw2] Parameters of grid ----')")
    write(51,"('*p1  : Number of subgrids in X (required to read the next line)')")
    write(51,*) wgtTMDPDF_numSubGridsX
    write(51,"('*p2  : Intervals for subgrids in X (must include 1., as the last point)')")
    write(51,*) wgtTMDPDF_subGridsX
    write(51,"('*p3  : Number of nodes in the X-subgrid')")
    write(51,*) wgtTMDPDF_grid_SizeX
    write(51,"('*p4  : Number of subgrids in B (required to read the next line)')")
    write(51,*) wgtTMDPDF_numSubGridsB
    write(51,"('*p5  : Intervals for subgrids in B (below and above ultimate points the value is frozen)')")
    write(51,*) wgtTMDPDF_subGridsB
    write(51,"('*p6  : Number of nodes in the B-subgrid ')")
    write(51,*) wgtTMDPDF_grid_SizeB
    write(51,"(' ')")
    write(51,"('*F   : ---- OPE[tw3] main definitions ----')")
    write(51,"('*p1  : Order of coefficient function')")
    write(51,*) trim(wgtTMDPDF_order_tw3)
    write(51,"('*p2  : Prepare grid')")
    write(51,*) wgtTMDPDF_makeGrid_tw3
    write(51,"('*p3  : run the test of the grid (takes some time)')")
    write(51,*) wgtTMDPDF_runGridTest_tw3
    write(51,"(' ')")
    write(51,"('*G   : ---- OPE[tw3] Parameters of grid ----')")
    write(51,"(' ')")
    write(51,"('*H   : ---- Transformation to KT-space ----')")
    write(51,"('*p1  : Tolerance (relative tolerance of summation in OGATA quadrature)')")
    write(51,*) wgtTMDPDF_toleranceOGATA
    write(51,"('*p2  : Ogata quadrature integration step')")
    write(51,*) wgtTMDPDF_hOGATA
    write(51,"('*p3  : Minimum value of kT (below that value function is constant)')")
    write(51,*) wgtTMDPDF_KT_FREEZE
    write(51,"(' ')")
    write(51,"('*I   : ---- Computation of Transverse Momentum Moments (TMM) ----')")
    write(51,"('*p1  : Tolerance (relative tolerance of summation in OGATA quadrature)')")
    write(51,*) wgtTMDPDF_toleranceOGATA_TMM
    write(51,"('*p2  : Ogata quadrature integration step')")
    write(51,*) wgtTMDPDF_hOGATA_TMM
    write(51,"('*p3  : Minimum value of mu [GeV] (below that value the computation is terminated)')")
    write(51,*) wgtTMDPDF_muMIN_TMM

    write(51,"(' ')")
    write(51,"(' ')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('# ----                         PARAMETERS OF BoerMuldersTMDPDF               -----')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('*14  :')")
    write(51,"('*p1  : initialize BoerMuldersTMDPDF module')")
    write(51,*) include_BoerMuldersTMDPDF
    write(51,"(' ')")
    write(51,"('*A   : ---- Main definitions ----')")
    write(51,"('*p1  : Include gluon TMDPDFs')")
    write(51,*) BoerMuldersTMDPDF_withGluon
    write(51,"('*p2  : Number of hadrons (in order starting with 1)')")
    write(51,*) BoerMuldersTMDPDF_numHadron
    write(51,"(' ')")
    write(51,"('*B   : ---- OPE[tw3] main definitions ----')")
    write(51,"('*p1  : Order of coefficient function')")
    write(51,*) trim(BoerMuldersTMDPDF_order)
    write(51,"('*p2  : Prepare grid')")
    write(51,*) BoerMuldersTMDPDF_makeGrid
    write(51,"('*p3  : run the test of the grid (takes some time)')")
    write(51,*) BoerMuldersTMDPDF_runGridTest
    write(51,"(' ')")
    write(51,"('*C   : ---- Parameters of NP model ----')")
    write(51,"('*p1  : Length of lambdaNP')")
    write(51,*) BoerMuldersTMDPDF_lambdaLength
    write(51,"('*p2  : Absolute maximum b (for larger b, TMD=0)')")
    write(51,*) BoerMuldersTMDPDF_BMAX_ABS
    write(51,"(' ')")
    write(51,"('*D   : ---- Numerical evaluation parameters ----')")
    write(51,"('*p1  : Tolerance (relative tolerance of convolution integral)')")
    write(51,*) BoerMuldersTMDPDF_toleranceINT
    write(51,"('*p2  : Tolerance general (used for various comparisons)')")
    write(51,*) BoerMuldersTMDPDF_toleranceGEN
    write(51,"('*p3  : Maximum number of iterations (for adaptive integration)')")
    write(51,*) BoerMuldersTMDPDF_maxIteration
    write(51,"(' ')")
    write(51,"('*E   : ---- (OPE) Parameters of grid ----')")
    write(51,"(' ')")
    write(51,"('*F   : ---- Transform and grid in KT-space ----')")
    write(51,"('*p1  : Prepare grid')")
    write(51,*) BoerMuldersTMDPDF_makeGrid_inKT
    write(51,"('*p2  : run the test of the grid (takes some time)')")
    write(51,*) BoerMuldersTMDPDF_runGridTest_inKT
    write(51,"('*p3  : Number of subgrids in X (required to read the next line)')")
    write(51,*) BoerMuldersTMDPDF_numSubGridsX_inKT
    write(51,"('*p4  : Intervals for subgrids in X (must include 1., as the last point)')")
    write(51,*) BoerMuldersTMDPDF_subGridsX_inKT
    write(51,"('*p5  : Number of nodes in the X-subgrid')")
    write(51,*) BoerMuldersTMDPDF_grid_SizeX_inKT
    write(51,"('*p6  : Number of subgrids in K (required to read the next line)')")
    write(51,*) BoerMuldersTMDPDF_numSubGridsKT_inKT
    write(51,"('*p7  : Intervals for subgrids in KT (below and above ultimate points the value is frozen)')")
    write(51,*) BoerMuldersTMDPDF_subGridsKT_inKT
    write(51,"('*p8  : Number of nodes in the KT-subgrid')")
    write(51,*) BoerMuldersTMDPDF_grid_SizeKT_inKT
    write(51,"('*p9  : Minimal Q in the grid')")
    write(51,*) BoerMuldersTMDPDF_minQ_inKT
    write(51,"('*p10 : Maximal Q in the grid')")
    write(51,*) BoerMuldersTMDPDF_maxQ_inKT
    write(51,"('*p11 : Number of nodes in the Q-grid')")
    write(51,*) BoerMuldersTMDPDF_grid_SizeQ_inKT
    write(51,"('*p12  : Number of subgrids in B (required to read the next line)')")
    write(51,*) BoerMuldersTMDPDF_numSubGridsB_inKT
    write(51,"('*p13  : Intervals for subgrids in B (below and above ultimate points the value is frozen)')")
    write(51,*) BoerMuldersTMDPDF_subGridsB_inKT
    write(51,"('*p14  : Number of nodes in the B-subgrid')")
    write(51,*) BoerMuldersTMDPDF_grid_SizeB_inKT
    write(51,"(' ')")
    write(51,"('*G   : ---- Computation of Transverse Momentum Moments (TMM) ----')")
    write(51,"('*p1  : Tolerance (relative tolerance of summation in OGATA quadrature)')")
    write(51,*) BoerMuldersTMDPDF_toleranceOGATA_TMM
    write(51,"('*p2  : Ogata quadrature integration step')")
    write(51,*) BoerMuldersTMDPDF_hOGATA_TMM
    write(51,"('*p3  : Minimum value of mu [GeV] (below that value the computation is terminated)')")
    write(51,*) BoerMuldersTMDPDF_muMIN_TMM

    write(51,"(' ')")
    write(51,"(' ')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('# ----                      PARAMETERS OF TMDF-KPC                      -----')")
    write(51,"('# ---------------------------------------------------------------------------')")
    write(51,"('*15  :')")
    write(51,"('*p1  : initialize TMDF-KPC-DY module')")
    write(51,*) include_TMDF_KPC
    write(51,"('*A   : ---- Main definitions ----')")
    write(51,"('*p1  : NOT YET')")

    write(51,"('*B   : ---- Numerical evaluation parameters ----')")
    write(51,"('*p1  : Tolerance general (variable comparision, etc.)')")
    write(51,*) TMDF_KPC_toleranceGEN
    write(51,"('*p2  : Tolerance integral (Integration tollerance)')")
    write(51,*) TMDF_KPC_toleranceINT
    write(51,"('*p3  : Minimum qT value (below the qT-value is frozen)')")
    write(51,*) TMDF_KPC_qTMIN

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
    else if(FILEversion<30) then
        write(*,*) color('aTMDe_setup: suggested setup-file is for artemide v2.'&
        , c_red_bold)
        write(*,*) color('..           It is incompatible with the preent version.'&
        , c_red)
        write(*,*) color('..           Please, use actual file or update it manually.'&
        ,c_red)
        stop
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
    call MoveTO(51,'*p2  ')
    read(51,*) Mnorm

    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) number_of_processors

    !# ----                           PARAMETERS OF QCDinput                 -----
    call MoveTO(51,'*1   ')
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) mCHARM
    call MoveTO(51,'*p2  ')
    read(51,*) mBOTTOM
    call MoveTO(51,'*p3  ')
    read(51,*) mTOP

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

    if(number_of_uFFs>0) then
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
    end if


    !-------PDF for lpTMDPDF
    call MoveTO(51,'*D   ')
    call MoveTO(51,'*p1  ')

    read(51,*) number_of_lpPDFs

    if(number_of_lpPDFs>0) then
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
    

    !-------helicity PDF
    call MoveTO(51,'*E   ')
    call MoveTO(51,'*p1  ')

    read(51,*) number_of_hPDFs
    if(number_of_hPDFs>0) then
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
    call MoveTO(51,'*p4  ')
    read(51,*) alphaQED_MTAU
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
    call MoveTO(51,'*D   ')
    call MoveTO(51,'*p1  ')
    read(51,*) mH
    call MoveTO(51,'*p2  ')
    read(51,*) GammaH
    call MoveTO(51,'*p3  ')
    read(51,*) vevH
    call MoveTO(51,'*E   ')
    call MoveTO(51,'*p1  ')
    read(51,*) mELECTRON
    call MoveTO(51,'*p2  ')
    read(51,*) mMUON
    call MoveTO(51,'*p3  ')
    read(51,*) mTAU

    !# ----                            PARAMETERS OF TMDR                    -----
    call MoveTO(51,'*3   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_TMDR
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDR_order
    call MoveTO(51,'*p2  ')
    read(51,*) TMDR_override
    call MoveTO(51,'*p3  ')
    read(51,*) TMDR_orderGAMMA
    call MoveTO(51,'*p4  ')
    read(51,*) TMDR_orderV
    call MoveTO(51,'*p5  ')
    read(51,*) TMDR_orderD
    call MoveTO(51,'*p6  ')
    read(51,*) TMDR_orderZETA

    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDR_lambdaLength

    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDR_tolerance
    call MoveTO(51,'*p2  ')
    read(51,*) TMDR_bFREEZE
    call MoveTO(51,'*p3  ')
    read(51,*) TMDR_smooth


    !# ----                           PARAMETERS OF uTMDPDF                  -----
    call MoveTO(51,'*4   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_uTMDPDF
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) uTMDPDF_withGluon
    call MoveTO(51,'*p2  ')
    read(51,*) uTMDPDF_numHadron
    if(uTMDPDF_numHadron/=number_of_uPDFs) then
        if(outputLevel>0) write(*,*) ' '
        if(outputLevel>0) write(*,*) color('ESSENTIAL INCONSITENCY: the number of uPDFs is unequal to the number of uTMDPDFs',c_red)
        if(outputLevel>0) write(*,*) color('                        it can lead to mistakes or crash',c_red)
        if(outputLevel>0) write(*,*) ' '
    end if
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) uTMDPDF_order
    call MoveTO(51,'*p2  ')
    read(51,*) uTMDPDF_makeGrid
    call MoveTO(51,'*p3  ')
    read(51,*) uTMDPDF_runGridTest
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) uTMDPDF_lambdaLength
    call MoveTO(51,'*p2  ')
    read(51,*) uTMDPDF_BMAX_ABS
    call MoveTO(51,'*D   ')
    call MoveTO(51,'*p1  ')
    read(51,*) uTMDPDF_toleranceINT
    call MoveTO(51,'*p2  ')
    read(51,*) uTMDPDF_toleranceGEN
    call MoveTO(51,'*p3  ')
    read(51,*) uTMDPDF_maxIteration
    call MoveTO(51,'*E   ')
    call MoveTO(51,'*p1  ')
    read(51,*) uTMDPDF_numSubGridsX
    deallocate(uTMDPDF_subGridsX)
    allocate(uTMDPDF_subGridsX(0:uTMDPDF_numSubGridsX))
    call MoveTO(51,'*p2  ')
    read(51,*) uTMDPDF_subGridsX
    call MoveTO(51,'*p3  ')
    read(51,*) uTMDPDF_grid_SizeX
    call MoveTO(51,'*p4  ')
    read(51,*) uTMDPDF_numSubGridsB
    deallocate(uTMDPDF_subGridsB)
    allocate(uTMDPDF_subGridsB(0:uTMDPDF_numSubGridsB))
    call MoveTO(51,'*p5  ')
    read(51,*) uTMDPDF_subGridsB
    call MoveTO(51,'*p6  ')
    read(51,*) uTMDPDF_grid_SizeB
    call MoveTO(51,'*F   ')
    call MoveTO(51,'*p1  ')
    read(51,*) uTMDPDF_makeGrid_inKT
    call MoveTO(51,'*p2  ')
    read(51,*) uTMDPDF_runGridTest_inKT
    call MoveTO(51,'*p3  ')
    read(51,*) uTMDPDF_numSubGridsX_inKT
    deallocate(uTMDPDF_subGridsX_inKT)
    allocate(uTMDPDF_subGridsX_inKT(0:uTMDPDF_numSubGridsX_inKT))
    call MoveTO(51,'*p4  ')
    read(51,*) uTMDPDF_subGridsX_inKT
    call MoveTO(51,'*p5  ')
    read(51,*) uTMDPDF_grid_SizeX_inKT
    call MoveTO(51,'*p6  ')
    read(51,*) uTMDPDF_numSubGridsKT_inKT
    deallocate(uTMDPDF_subGridsKT_inKT)
    allocate(uTMDPDF_subGridsKT_inKT(0:uTMDPDF_numSubGridsKT_inKT))
    call MoveTO(51,'*p7  ')
    read(51,*) uTMDPDF_subGridsKT_inKT
    call MoveTO(51,'*p8  ')
    read(51,*) uTMDPDF_grid_SizeKT_inKT
    call MoveTO(51,'*p9  ')
    read(51,*) uTMDPDF_minQ_inKT
    call MoveTO(51,'*p10 ')
    read(51,*) uTMDPDF_maxQ_inKT
    call MoveTO(51,'*p11 ')
    read(51,*) uTMDPDF_grid_SizeQ_inKT
    call MoveTO(51,'*p12 ')
    read(51,*) uTMDPDF_numSubGridsB_inKT
    deallocate(uTMDPDF_subGridsB_inKT)
    allocate(uTMDPDF_subGridsB_inKT(0:uTMDPDF_numSubGridsB_inKT))
    call MoveTO(51,'*p13 ')
    read(51,*) uTMDPDF_subGridsB_inKT
    call MoveTO(51,'*p14 ')
    read(51,*) uTMDPDF_grid_SizeB_inKT
    call MoveTO(51,'*G   ')
    call MoveTO(51,'*p1  ')
    read(51,*) uTMDPDF_toleranceOGATA_TMM
    call MoveTO(51,'*p2  ')
    read(51,*) uTMDPDF_hOGATA_TMM
    call MoveTO(51,'*p3  ')
    read(51,*) uTMDPDF_muMIN_TMM

    !# ----                           PARAMETERS OF uTMDFF                   -----
    call MoveTO(51,'*5   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_uTMDFF
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) uTMDFF_withGluon
    call MoveTO(51,'*p2  ')
    read(51,*) uTMDFF_numHadron
    if(uTMDFF_numHadron/=number_of_uFFs) then
        if(outputLevel>0) write(*,*) ' '
        if(outputLevel>0) write(*,*) color('ESSENTIAL INCONSITENCY: the number of uFFs is unequal to the number of uTMDFFs',c_red)
        if(outputLevel>0) write(*,*) color('                        it can lead to mistakes or crash',c_red)
        if(outputLevel>0) write(*,*) ' '
    end if
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) uTMDFF_order
    call MoveTO(51,'*p2  ')
    read(51,*) uTMDFF_makeGrid
    call MoveTO(51,'*p3  ')
    read(51,*) uTMDFF_runGridTest
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) uTMDFF_lambdaLength
    call MoveTO(51,'*p2  ')
    read(51,*) uTMDFF_BMAX_ABS
    call MoveTO(51,'*D   ')
    call MoveTO(51,'*p1  ')
    read(51,*) uTMDFF_toleranceINT
    call MoveTO(51,'*p2  ')
    read(51,*) uTMDFF_toleranceGEN
    call MoveTO(51,'*p3  ')
    read(51,*) uTMDFF_maxIteration
    call MoveTO(51,'*E   ')
    call MoveTO(51,'*p1  ')
    read(51,*) uTMDFF_numSubGridsX
    deallocate(uTMDFF_subGridsX)
    allocate(uTMDFF_subGridsX(0:uTMDFF_numSubGridsX))
    call MoveTO(51,'*p2  ')
    read(51,*) uTMDFF_subGridsX
    call MoveTO(51,'*p3  ')
    read(51,*) uTMDFF_grid_SizeX
    call MoveTO(51,'*p4  ')
    read(51,*) uTMDFF_numSubGridsB
    deallocate(uTMDFF_subGridsB)
    allocate(uTMDFF_subGridsB(0:uTMDFF_numSubGridsB))
    call MoveTO(51,'*p5  ')
    read(51,*) uTMDFF_subGridsB
    call MoveTO(51,'*p6  ')
    read(51,*) uTMDFF_grid_SizeB
    call MoveTO(51,'*G   ')
    call MoveTO(51,'*p1  ')
    read(51,*) uTMDFF_toleranceOGATA_TMM
    call MoveTO(51,'*p2  ')
    read(51,*) uTMDFF_hOGATA_TMM
    call MoveTO(51,'*p3  ')
    read(51,*) uTMDFF_muMIN_TMM

    !# ----                            PARAMETERS OF TMDF                    -----
    call MoveTO(51,'*7   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_TMDF
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDF_tolerance
    call MoveTO(51,'*p2  ')
    read(51,*) TMDF_OGATAh
    call MoveTO(51,'*p3  ')
    read(51,*) TMDF_qTMIN
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDF_mass


    !# ----                           PARAMETERS OF TMDX-DY                  -----
    call MoveTO(51,'*9   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_TMDX_DY
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDX_DY_order
    call MoveTO(51,'*p2  ')
    read(51,*) TMDX_DY_useKPC
    call MoveTO(51,'*p3  ')
    read(51,*) TMDX_DY_piResum
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDX_DY_toleranceGEN
    call MoveTO(51,'*p2  ')
    read(51,*) TMDX_DY_toleranceINT
    call MoveTO(51,'*p3  ')
    read(51,*) TMDX_DY_ptSECTION
    call MoveTO(51,'*p4  ')
    read(51,*) TMDX_DY_maxQbinSize
    call MoveTO(51,'*p5  ')
    read(51,*) TMDX_DY_minqTabs
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDX_DY_exactX1X2
    call MoveTO(51,'*p2  ')
    read(51,*) TMDX_DY_exactScale
    call MoveTO(51,'*D   ')


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
    call MoveTO(51,'*p5  ')
    read(51,*) TMDX_SIDIS_qTinX1Z1corr
    call MoveTO(51,'*p6  ')
    read(51,*) TMDX_SIDIS_exactScale
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDX_SIDIS_tolerance
    call MoveTO(51,'*p2  ')
    read(51,*) TMDX_SIDIS_ptSECTION
    call MoveTO(51,'*p3  ')
    read(51,*) TMDX_SIDIS_toleranceZ
    call MoveTO(51,'*p4  ')
    read(51,*) TMDX_SIDIS_methodZ
    call MoveTO(51,'*p5  ')
    read(51,*) TMDX_SIDIS_toleranceX
    call MoveTO(51,'*p6  ')
    read(51,*) TMDX_SIDIS_methodX
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDX_SIDIS_numProc

    !# ----                           PARAMETERS OF lpTMDPDF                  -----
    call MoveTO(51,'*11   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_lpTMDPDF
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) lpTMDPDF_withGluon
    if(.not. lpTMDPDF_withGluon) then
        if(outputLevel>0) write(*,*) color('MINOR INCONSITENCY: gluon must be included in lpTMDPDFs.',c_yellow)
    end if
    call MoveTO(51,'*p2  ')
    read(51,*) lpTMDPDF_numHadron
    if(lpTMDPDF_numHadron/=number_of_lpPDFs) then
        if(outputLevel>0) write(*,*) ' '
        if(outputLevel>0) write(*,*) &
            color('ESSENTIAL INCONSITENCY: the number of lpPDFs is unequal to the number of lpTMDPDFs',c_red)
        if(outputLevel>0) write(*,*) color('                        it can lead to mistakes or crash',c_red)
        if(outputLevel>0) write(*,*) ' '
    end if
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) lpTMDPDF_order
    call MoveTO(51,'*p2  ')
    read(51,*) lpTMDPDF_makeGrid
    call MoveTO(51,'*p3  ')
    read(51,*) lpTMDPDF_runGridTest
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) lpTMDPDF_lambdaLength
    call MoveTO(51,'*p2  ')
    read(51,*) lpTMDPDF_BMAX_ABS
    call MoveTO(51,'*D   ')
    call MoveTO(51,'*p1  ')
    read(51,*) lpTMDPDF_toleranceINT
    call MoveTO(51,'*p2  ')
    read(51,*) lpTMDPDF_toleranceGEN
    call MoveTO(51,'*p3  ')
    read(51,*) lpTMDPDF_maxIteration
    call MoveTO(51,'*E   ')
    call MoveTO(51,'*p1  ')
    read(51,*) lpTMDPDF_numSubGridsX
    deallocate(lpTMDPDF_subGridsX)
    allocate(lpTMDPDF_subGridsX(0:lpTMDPDF_numSubGridsX))
    call MoveTO(51,'*p2  ')
    read(51,*) lpTMDPDF_subGridsX
    call MoveTO(51,'*p3  ')
    read(51,*) lpTMDPDF_grid_SizeX
    call MoveTO(51,'*p4  ')
    read(51,*) lpTMDPDF_numSubGridsB
    deallocate(lpTMDPDF_subGridsB)
    allocate(lpTMDPDF_subGridsB(0:lpTMDPDF_numSubGridsB))
    call MoveTO(51,'*p5  ')
    read(51,*) lpTMDPDF_subGridsB
    call MoveTO(51,'*p6  ')
    read(51,*) lpTMDPDF_grid_SizeB
    call MoveTO(51,'*G   ')
    call MoveTO(51,'*p1  ')
    read(51,*) lpTMDPDF_toleranceOGATA_TMM
    call MoveTO(51,'*p2  ')
    read(51,*) lpTMDPDF_hOGATA_TMM
    call MoveTO(51,'*p3  ')
    read(51,*) lpTMDPDF_muMIN_TMM
    
    !# ----                           PARAMETERS OF SiversTMDPDF                  -----
    call MoveTO(51,'*12   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_SiversTMDPDF
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) SiversTMDPDF_withGluon
    call MoveTO(51,'*p2  ')
    read(51,*) SiversTMDPDF_numHadron
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) SiversTMDPDF_order
    call MoveTO(51,'*p2  ')
    read(51,*) SiversTMDPDF_makeGrid
    call MoveTO(51,'*p3  ')
    read(51,*) SiversTMDPDF_runGridTest
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) SiversTMDPDF_lambdaLength
    call MoveTO(51,'*p2  ')
    read(51,*) SiversTMDPDF_BMAX_ABS
    call MoveTO(51,'*D   ')
    call MoveTO(51,'*p1  ')
    read(51,*) SiversTMDPDF_toleranceINT
    call MoveTO(51,'*p2  ')
    read(51,*) SiversTMDPDF_toleranceGEN
    call MoveTO(51,'*p3  ')
    read(51,*) SiversTMDPDF_maxIteration
    call MoveTO(51,'*G   ')
    call MoveTO(51,'*p1  ')
    read(51,*) SiversTMDPDF_toleranceOGATA_TMM
    call MoveTO(51,'*p2  ')
    read(51,*) SiversTMDPDF_hOGATA_TMM
    call MoveTO(51,'*p3  ')
    read(51,*) SiversTMDPDF_muMIN_TMM
    
    !# ----                           PARAMETERS OF wgtTMDPDF                  -----
    call MoveTO(51,'*13   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_wgtTMDPDF
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) wgtTMDPDF_withGluon
    call MoveTO(51,'*p2  ')
    read(51,*) wgtTMDPDF_numHadron
    if(wgtTMDPDF_numHadron/=number_of_hPDFs) then
        if(outputLevel>0) write(*,*) ' '
        if(outputLevel>0) write(*,*) &
        color('ESSENTIAL INCONSITENCY: the number of hPDFs is unequal to the number of wgtTMDPDFs',c_red)
        if(outputLevel>0) write(*,*) color('                        it can lead to mistakes or crash',c_red)
        if(outputLevel>0) write(*,*) ' '
    end if
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) wgtTMDPDF_order
    call MoveTO(51,'*p2  ')
    read(51,*) wgtTMDPDF_makeGrid
    call MoveTO(51,'*p3  ')
    read(51,*) wgtTMDPDF_runGridTest
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) wgtTMDPDF_lambdaLength
    call MoveTO(51,'*p2  ')
    read(51,*) wgtTMDPDF_BMAX_ABS
    call MoveTO(51,'*D   ')
    call MoveTO(51,'*p1  ')
    read(51,*) wgtTMDPDF_toleranceINT
    call MoveTO(51,'*p2  ')
    read(51,*) wgtTMDPDF_toleranceGEN
    call MoveTO(51,'*p3  ')
    read(51,*) wgtTMDPDF_maxIteration
    call MoveTO(51,'*E   ')
    call MoveTO(51,'*p1  ')
    read(51,*) wgtTMDPDF_numSubGridsX
    deallocate(wgtTMDPDF_subGridsX)
    allocate(wgtTMDPDF_subGridsX(0:wgtTMDPDF_numSubGridsX))
    call MoveTO(51,'*p2  ')
    read(51,*) wgtTMDPDF_subGridsX
    call MoveTO(51,'*p3  ')
    read(51,*) wgtTMDPDF_grid_SizeX
    call MoveTO(51,'*p4  ')
    read(51,*) wgtTMDPDF_numSubGridsB
    deallocate(wgtTMDPDF_subGridsB)
    allocate(wgtTMDPDF_subGridsB(0:wgtTMDPDF_numSubGridsB))
    call MoveTO(51,'*p5  ')
    read(51,*) wgtTMDPDF_subGridsB
    call MoveTO(51,'*p6  ')
    read(51,*) wgtTMDPDF_grid_SizeB
    call MoveTO(51,'*F   ')
    call MoveTO(51,'*p1  ')
    read(51,*) wgtTMDPDF_order_tw3
    call MoveTO(51,'*p2  ')
    read(51,*) wgtTMDPDF_makeGrid_tw3
    call MoveTO(51,'*p3  ')
    read(51,*) wgtTMDPDF_runGridTest_tw3
    call MoveTO(51,'*I   ')
    call MoveTO(51,'*p1  ')
    read(51,*) wgtTMDPDF_toleranceOGATA_TMM
    call MoveTO(51,'*p2  ')
    read(51,*) wgtTMDPDF_hOGATA_TMM
    call MoveTO(51,'*p3  ')
    read(51,*) wgtTMDPDF_muMIN_TMM


    !# ----                           PARAMETERS OF BoerMuldersTMDPDF                  -----
    call MoveTO(51,'*14   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_BoerMuldersTMDPDF
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) BoerMuldersTMDPDF_withGluon
    call MoveTO(51,'*p2  ')
    read(51,*) BoerMuldersTMDPDF_numHadron
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) BoerMuldersTMDPDF_order
    call MoveTO(51,'*p2  ')
    read(51,*) BoerMuldersTMDPDF_makeGrid
    call MoveTO(51,'*p3  ')
    read(51,*) BoerMuldersTMDPDF_runGridTest
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) BoerMuldersTMDPDF_lambdaLength
    call MoveTO(51,'*p2  ')
    read(51,*) BoerMuldersTMDPDF_BMAX_ABS
    call MoveTO(51,'*D   ')
    call MoveTO(51,'*p1  ')
    read(51,*) BoerMuldersTMDPDF_toleranceINT
    call MoveTO(51,'*p2  ')
    read(51,*) BoerMuldersTMDPDF_toleranceGEN
    call MoveTO(51,'*p3  ')
    read(51,*) BoerMuldersTMDPDF_maxIteration
    call MoveTO(51,'*F   ')
    call MoveTO(51,'*p1  ')
    read(51,*) BoerMuldersTMDPDF_makeGrid_inKT
    call MoveTO(51,'*p2  ')
    read(51,*) BoerMuldersTMDPDF_runGridTest_inKT
    call MoveTO(51,'*p3  ')
    read(51,*) BoerMuldersTMDPDF_numSubGridsX_inKT
    deallocate(BoerMuldersTMDPDF_subGridsX_inKT)
    allocate(BoerMuldersTMDPDF_subGridsX_inKT(0:BoerMuldersTMDPDF_numSubGridsX_inKT))
    call MoveTO(51,'*p4  ')
    read(51,*) BoerMuldersTMDPDF_subGridsX_inKT
    call MoveTO(51,'*p5  ')
    read(51,*) BoerMuldersTMDPDF_grid_SizeX_inKT
    call MoveTO(51,'*p6  ')
    read(51,*) BoerMuldersTMDPDF_numSubGridsKT_inKT
    deallocate(BoerMuldersTMDPDF_subGridsKT_inKT)
    allocate(BoerMuldersTMDPDF_subGridsKT_inKT(0:BoerMuldersTMDPDF_numSubGridsKT_inKT))
    call MoveTO(51,'*p7  ')
    read(51,*) BoerMuldersTMDPDF_subGridsKT_inKT
    call MoveTO(51,'*p8  ')
    read(51,*) BoerMuldersTMDPDF_grid_SizeKT_inKT
    call MoveTO(51,'*p9  ')
    read(51,*) BoerMuldersTMDPDF_minQ_inKT
    call MoveTO(51,'*p10 ')
    read(51,*) BoerMuldersTMDPDF_maxQ_inKT
    call MoveTO(51,'*p11 ')
    read(51,*) BoerMuldersTMDPDF_grid_SizeQ_inKT
    call MoveTO(51,'*p12 ')
    read(51,*) BoerMuldersTMDPDF_numSubGridsB_inKT
    deallocate(BoerMuldersTMDPDF_subGridsB_inKT)
    allocate(BoerMuldersTMDPDF_subGridsB_inKT(0:BoerMuldersTMDPDF_numSubGridsB_inKT))
    call MoveTO(51,'*p13 ')
    read(51,*) BoerMuldersTMDPDF_subGridsB_inKT
    call MoveTO(51,'*p14 ')
    read(51,*) BoerMuldersTMDPDF_grid_SizeB_inKT
    call MoveTO(51,'*G   ')
    call MoveTO(51,'*p1  ')
    read(51,*) BoerMuldersTMDPDF_toleranceOGATA_TMM
    call MoveTO(51,'*p2  ')
    read(51,*) BoerMuldersTMDPDF_hOGATA_TMM
    call MoveTO(51,'*p3  ')
    read(51,*) BoerMuldersTMDPDF_muMIN_TMM

    !# ----                            PARAMETERS OF TMDF-KPC              -----
    call MoveTO(51,'*15   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_TMDF_KPC
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    !read(51,*) TMDF_tolerance
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) TMDF_KPC_toleranceGEN
    call MoveTO(51,'*p2  ')
    read(51,*) TMDF_KPC_toleranceINT
    call MoveTO(51,'*p3  ')
    read(51,*) TMDF_KPC_qTMIN

    CLOSE (51, STATUS='KEEP') 

    if(outputLevel>1) write(*,*) color('aTMDe_setup: constants-file loaded sucessfully.',c_green_bold)

end subroutine ReadConstantsFile

end module aTMDe_Setup
