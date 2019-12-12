!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 2.01
!
!	Evaluation of the linearly-polarized gluon TMD PDF at low normalization point in zeta-prescription.
!	
!	if you use this module please, quote ????.????
!
!
!				A.Vladimirov (13.06.2016)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module lpTMDPDF
use IO_functions
use QCDinput
use lpTMDPDF_model

implicit none
!------------------------LOCALs -----------------------------------------------

 private 
 
  !Current version of module
 character (len=8),parameter :: moduleName="lpTMDPDF"
 character (len=5),parameter :: version="v2.02" 
 !Last appropriate verion of constants-file
  integer,parameter::inputver=3
 
  INCLUDE 'Tables/NumConst.f90'
  INCLUDE 'Tables/G7K15.f90'

!--------------------------------Working variables-----------------------------------------------!--- general
  logical:: started=.false.
!! Level of output
!! 0=only critical
!! 1=initialization details
!! 2=WARNINGS
  integer::outputLevel=2
  !! variable that count number of WRNING mesagges. In order not to spam too much
  integer::messageTrigger=6
  
  !!! The global order which is used in programm
  ! it is set by SetOrderForC subroutine
  ! 0=LO, 1=NLO, 2=NNLO
  integer :: order_global
!   
  integer::lambdaNPlength
  real*8,dimension(:),allocatable::lambdaNP
  real*8,dimension(:),allocatable::lambdaNP_grid !!! this is the value of lambda on which the grid is build

  real*8::c4_global!!!this is the variation constant for mu_OPE
  
    !!!Parameteris of numerics
  real*8 :: tolerance=0.0001d0!!! relative tolerance of the integration
  integer :: maxIteration=5000
  
!------------------------------Variables for coefficient function etc-------------------------------
  
  !!!!!Coefficient lists
  integer,parameter::parametrizationLength=14
  !! { 1/x, log[x]/x  !exact
  !! Log[x], log[x]^2 !exact
  !! 1 (exact), (1-x), x(1-x)
  !! (1-x)Log[1-x]/x, x Log[x], x^2 Log[x]
  !! (1-x) Log[1-x], (1-x)^2 Log[1-x]
  !! Log[x]^2Log[1-x],Log[x]Log[1-x]^2
  !! The Lmu and Nf parts are exact the later parts are fitted
  real*8,dimension(1:parametrizationLength) :: Coeff_q_q, Coeff_q_g, Coeff_g_q, Coeff_g_g, Coeff_q_qb, Coeff_q_qp
  !! This is list of coefficeints for the encoding the singular at x->1
  !! { 1/(1-x), (Log[1-x]/(1-x))_+}
  !! they are zero!!!
  real*8, dimension(1:2) :: CoeffSing1_q_q,CoeffSing1_g_g
  
  integer :: counter,messageCounter
  
  INCLUDE 'CommonCode/Twist2Convolution-VAR.f90'
  INCLUDE 'CommonCode/Twist2Grid-VAR.f90'
  
!!--------------------------------- variables for the griding the TMD.---------------------------------------------
  logical :: gridReady!!!!indicator that grid is ready to use. If it is .true., the TMD calculated from the grid
  logical :: prepareGrid!!!idicator that grid must be prepared
  logical,parameter :: withGluon=.true.!!!indicator the gluon is needed in the grid
  
  integer::numberOfHadrons
  integer,dimension(:),allocatable::hadronsInGRID

!!-----------------------------------------------Public interface---------------------------------------------------
    
  public::lpTMDPDF_Initialize,lpTMDPDF_SetLambdaNP,lpTMDPDF_SetScaleVariation,lpTMDPDF_resetGrid,lpTMDPDF_SetPDFreplica
  public::lpTMDPDF_IsInitialized,lpTMDPDF_CurrentNPparameters
  public::lpTMDPDF_lowScale50
!   public::CheckCoefficient  
!   public::mu_OPE
  
  interface lpTMDPDF_SetLambdaNP
    module procedure lpTMDPDF_SetLambdaNP_usual,lpTMDPDF_SetReplica_optional
  end interface
  
  contains
  
  INCLUDE 'CommonCode/Twist2Convolution.f90'
  INCLUDE 'CommonCode/Twist2Grid.f90'
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Interface subroutines!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  function lpTMDPDF_IsInitialized()
  logical::lpTMDPDF_IsInitialized
  lpTMDPDF_IsInitialized=started
  end function lpTMDPDF_IsInitialized

   !! Initialization of the package
  subroutine lpTMDPDF_Initialize(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=300)::path,line
    logical::initRequared
    character(len=8)::orderMain
    logical::bSTAR_lambdaDependent,dummyLogical
    integer::i,FILEver
    
    if(started) return
    
    if(.not.QCDinput_IsInitialized()) then
      if(outputLevel>1) write(*,*) '.. initializing QCDinput (from ',moduleName,')'
      if(present(prefix)) then
      	call QCDinput_Initialize(file,prefix)
      else
	call QCDinput_Initialize(file)
      end if
    end if
  
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
    
    
    call MoveTO(51,'*11  ')
    call MoveTO(51,'*p1  ')
    read(51,*) initRequared
    if(.not.initRequared) then
      if(outputLevel>2) write(*,*)'artemide.',moduleName,': initialization is not requared. '
      started=.false.
      return
    end if
    
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) orderMain
    
    SELECT CASE(trim(orderMain))
      CASE ("LO")
	if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: LO'
	order_global=0
	if(outputLevel>0)write(*,*) 'WARNING: arTeMiDe.lpTMDPDF_Initialize: lin.pol.gluon at LO are zero!'
      CASE ("LO+")
	if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: LO+'
	order_global=0
	if(outputLevel>0)write(*,*) 'WARNING: arTeMiDe.lpTMDPDF_Initialize: lin.pol.gluon at LO are zero!'
      CASE ("NLO")
	if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NLO'
	order_global=1
      CASE ("NLO+")
	if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NLO+'
	order_global=1
      CASE ("NNLO")
	if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NNLO'
	order_global=2
      CASE ("NNLO+")
	if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NNLO+'
	order_global=2
      CASE DEFAULT
	if(outputLevel>0)write(*,*) 'WARNING: arTeMiDe.lpTMDPDF_Initialize: unknown order for coefficient function. Switch to NLO.'
	if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NLO'
	order_global=1
     END SELECT
    
     if(outputLevel>2) then
      write(*,'(A,I1)') ' |  Coef.func.    =as^',order_global
     end if
    
    !-------------parameters of NP model
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) lambdaNPlength
    
    if(outputLevel>2) write(*,'(A,I3)') ' Number of NP parameters =',lambdaNPlength
    
    if(lambdaNPlength<=0) then
    write(*,*) 'ERROR: arTeMiDe.lpTMDPDF_Initialize: number of non-pertrubative &
		    parameters should be >=1. Check the constants-file. Evaluation STOP'
    stop
    end if
    
    allocate(lambdaNP(1:lambdaNPlength))
    call MoveTO(51,'*p2  ')
    do i=1,lambdaNPlength
      read(51,*) lambdaNP(i)
    end do
    
    !-------------Numeric parameters
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) tolerance
    call MoveTO(51,'*p2  ')
    read(51,*) maxIteration
    
    if(outputLevel>2) then
      write(*,'(A,ES10.3)') ' |  tolerance     =',tolerance
      write(*,'(A,ES10.3)') ' |  max iteration =',REAL(maxIteration)
     end if
     
    !-------------Make grid options
    call MoveTO(51,'*D   ')
    call MoveTO(51,'*p1  ')
    read(51,*) prepareGrid
    call MoveTO(51,'*p2  ')
    read(51,*) dummyLogical
    if(.not.dummyLogical) then     
      if(outputLevel>0) write(*,*) 'WARNING: arTeMiDe.lpTMDPDF_Initialize:&
		  gluon option is switched ON (otherwise there is no reason to evaluate it)'
    end if
    call MoveTO(51,'*p3  ')
    read(51,*) numberOfHadrons
    allocate(hadronsInGRID(1:numberOfHadrons))
    call MoveTO(51,'*p4  ')
    read(51,*) hadronsInGRID
    
    !-------------Parameters of grid
    call MoveTO(51,'*E   ')
    call MoveTO(51,'*p1  ')
    read(51,*) xGrid_Min
    call MoveTO(51,'*p2  ')
    read(51,*) bGrid_Max
    call MoveTO(51,'*p3  ')
    read(51,*) GridSizeX
    call MoveTO(51,'*p4  ')
    read(51,*) GridSizeB
    call MoveTO(51,'*p5  ')
    read(51,*) slope
    
     if(outputLevel>2) then
      write(*,*) 'Grid options:'
      write(*,'(A,ES10.3)') ' |  xGrid_Min                 =',xGrid_Min
      write(*,'(A,ES10.3)') ' |  bGrid_Max                 =',bGrid_Max 
      write(*,'(A,I6,A,I6,A)') ' |  (GridSizeX,GridSizeB)     =(',GridSizeX,',',GridSizeB,')'
      write(*,'(A,F6.3)') ' |  slope                     =',slope 
      write(*,'(A,I3)')   ' |  hadrons to grid           =',numberOfHadrons
      write(*,*)   ' | list of hadrons in grid    =(',hadronsInGRID,')'
     end if
    
    CLOSE (51, STATUS='KEEP') 
    
    allocate(gridMain(0:GridSizeX,0:GridSizeB,-5:5,1:numberOfHadrons))
    allocate(boundaryValues(0:GridSizeX,-5:5,1:numberOfHadrons))
    
    allocate(lambdaNP_grid(1:lambdaNPlength))
    
    c4_global=1d0
    
    call ModelInit(outputLevel,lambdaNPlength)
    
    !!!!!!!Checking the x-dependance of muOPE
    IsMuXdependent=testMU()
    
    if(IsMuXdependent) then
    if(outputLevel>2) write(*,*) 'arTeMiDe.lpTMDPDF: mu OPE is dependent on x'
    else
    if(outputLevel>2) write(*,*) 'arTeMiDe.lpTMDPDF: mu OPE is independent on x'
    end if
    
        !!!!!!!Checking the lambda-dependance of bSTAR
    bSTAR_lambdaDependent=testbSTAR(lambdaNPlength)
    
    if(bSTAR_lambdaDependent) then
    if(outputLevel>2) write(*,*) 'arTeMiDe.lpTMDPDF: bSTAR is dependent on lambda'
    else
    if(outputLevel>2) write(*,*) 'arTeMiDe.lpTMDPDF: bSTAR is independent on lambda'
    end if
    
    !!!!!!!Checking the x-dependance of FNP
    IsFnpZdependent=TestFNP(hadronsInGRID,lambdaNPlength)
    
    gridReady=.false.
    
    if(IsFnpZdependent) then
    if(outputLevel>2) write(*,*) 'arTeMiDe.lpTMDPDF: FNP is dependent on z'
    else
    if(outputLevel>2) write(*,*) 'arTeMiDe.lpTMDPDF: FNP is independent on z'
    end if
    
    !!! if fnp depende on z or bSTAR depeds on lambda
    !!! grid must be recalculate ech time. It canbe saved to single IsFnpZdependent
    if(IsFnpZdependent .or. bSTAR_lambdaDependent) then
      IsFnpZdependent=.true.
      if(outputLevel>2) write(*,*) 'arTeMiDe.lpTMDPDF: ............. convolution is lambda sensitive.'
    end if
    
    started=.true.
    messageCounter=0
    
    if(outputLevel>0) write(*,*) '----- arTeMiDe.lpTMDPDF ',version,': .... initialized'
    if(outputLevel>1) write(*,*) ' '
  end subroutine lpTMDPDF_Initialize

    !! call for parameters from the model
    !! gluonRequared option is ignored
  subroutine lpTMDPDF_SetReplica_optional(num,buildGrid, gluonRequared)
  integer:: num
  logical,optional:: buildGrid,gluonRequared
  
  if(present(buildGrid)) then
      call lpTMDPDF_SetLambdaNP_usual(GiveReplicaParameters(num),buildGrid=buildGrid)
  else
      call lpTMDPDF_SetLambdaNP_usual(GiveReplicaParameters(num))
   end if
  
  end subroutine lpTMDPDF_SetReplica_optional
  
  !! call QCDinput to change the PDF replica number
  !! unset the grid, since it should be recalculated fro different PDF replica.
  subroutine lpTMDPDF_SetPDFreplica(rep)
  integer:: rep
  
  call QCDinput_SetPDFreplica(rep)
  gridReady=.false.  
  call lpTMDPDF_resetGrid()
  end subroutine lpTMDPDF_SetPDFreplica
  
    !!!Sets the non-pertrubative parameters lambda
    !!! carries additionl option to build the grid
    !!! if need to build grid, specify the gluon requared directive.
    !!! gluon gluonRequared option is ignored
  subroutine lpTMDPDF_SetLambdaNP_usual(lambdaIN,buildGrid, gluonRequared)
    real*8,intent(in)::lambdaIN(:)
    logical,optional :: buildGrid,gluonRequared
    real*8,dimension(1:lambdaNPlength)::lambdaOLD
    logical::IsNewValues
    integer::i,ll
    messageCounter=0
    
    if(present(buildGrid)) prepareGrid=buildGrid
    
    ll=size(lambdaIN)
    if(ll<lambdaNPlength) then 
      if(outputLevel>0) write(*,"('arTeMiDe.',A,'SetLambdaNP: WARNING length of lambdaNP(,',I3,') is less then requred (',I3,')')")&
	      moduleName,ll,lambdaNPlength
      if(outputLevel>0) write(*,*)'                Reset parameters are replaced by zeros!'
      lambdaNP=0d0*lambdaNP
      lambdaNP(1:ll)=lambdaIN(1:ll)
    else if (ll>lambdaNPlength) then
      if(outputLevel>0) write(*,"('arTeMiDe.',A,'SetLambdaNP: WARNING&
	      length of lambdaNP(,',I3,') is greater then requred (',I3,')')") moduleName,ll,lambdaNPlength
      if(outputLevel>0) write(*,*)'                Array is truncated!'
      lambdaNP(1:lambdaNPlength)=lambdaIN(1:lambdaNPlength)
     else
      lambdaOLD=lambdaNP
      lambdaNP=lambdaIN
     end if
    
    
    lambdaOLD=lambdaNP
    lambdaNP=lambdaIN
    IsNewValues=.false.
    do i=1,lambdaNPlength
     if(ABS(lambdaNP(i)-lambdaOLD(i))>10d-10) then
      IsNewValues=.true.
      exit
     end if
    end do    
    
    
    if(IsNewValues.and.(outputLevel>2)) write(*,*) 'arTeMiDe.',moduleName,': NPparameters reset = (',lambdaNP,')'
    
    !! further if's are only for griding
    
    if(prepareGrid) then !!!grid is requred
    
    if(IsNewValues) then  !! values are new
    
      if(gridReady) then  !!! grid is already build
      
	if(IsFnpZdependent) then !!! check the z-dependance of FNP
	  !! if it is z- dependent, rebuild the grid
	  gridReady=.false.
	  call MakeGrid()
	  gridReady=.true.
	  
	else !!! if z-Independent just do nothing.
	 if(outputLevel>2) write(*,*) 'arTeMiDe.',moduleName,': the values are to be restored from the initial grid'
	end if
      
      else !!! grid is not ready (how comes?)
	call MakeGrid()
	gridReady=.true.
      end if
	
    else  !! values are old
    
      if(gridReady) then  !!! grid is already build
	  !!!nothing to do
      else!!rare option then parameters are not new but grit is not build
        if(outputLevel>2) write(*,*) 'arTeMiDe.',moduleName,': parameters are not reset. But grid is not ready.'
	call MakeGrid()
	gridReady=.true.
      end if
      
    end if
    
    else
      gridReady=.false.
    end if 
    
  end subroutine lpTMDPDF_SetLambdaNP_usual
  
  !!! returns current value of NP parameters
  subroutine lpTMDPDF_CurrentNPparameters(var)
    real*8,dimension(1:lambdaNPlength)::var
    var=lambdaNP
 end subroutine lpTMDPDF_CurrentNPparameters
  
  !!! This subroutine ask for the grid reconstruction (or destruction)
  !!! gluon gluonRequared option is ignored
  subroutine lpTMDPDF_resetGrid(buildGrid,gluonRequared)
    logical,optional::buildGrid,gluonRequared
    logical::previousState
    
    if(present(buildGrid)) prepareGrid=buildGrid
    
    previousState=gridReady
    gridReady=.false.
    
    !! we recalculate grid only if it was already calculated!
      if(prepareGrid .and. previousState) then
        if(outputLevel>1) write(*,*) 'arTeMiDe ',moduleName,':  Grid Reset. with c4=',c4_global
	call MakeGrid()
	gridReady=.true.
      end if
  end subroutine lpTMDPDF_resetGrid
  
  
  !!!! this routine set the variations of scales
  !!!! it is used for the estimation of errors
  subroutine lpTMDPDF_SetScaleVariation(c4_in)
    real*8::c4_in
    if(c4_in<0.1d0 .or. c4_in>10.d0) then
      if(outputLevel>0) write(*,*) 'WARNING: arTeMiDe.lpTMDPDF: variation in c4 is enourmous. c4 is set to 2'
      c4_global=2d0
    else
      c4_global=c4_in
    end if
    c4_global=c4_in
    if(outputLevel>1) write(*,*) 'lpTMDPDF: set scale variations constant c4 as:',c4_global
    call lpTMDPDF_resetGrid()
  end subroutine lpTMDPDF_SetScaleVariation

!-------------------------------------------------
 !!!!array of x times PDF(x,Q) for hadron 'hadron'
 !!!! array is (-5:5) (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
 function xf(x,Q,hadron)
      real*8 :: x,Q
      integer:: hadron
      real*8, dimension(-5:5):: xf
      
      xf=x_lp_PDF(x,Q,hadron)
      
  end function xf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!COEFFICIENT FUNCTIONS!!
  
  !!! the function which contains the functions of parameterizations
  function parametrizationString(z)
  real*8::z,lz,llz,zz
  real*8,dimension(1:parametrizationLength)::parametrizationString
      zz=1d0-z
      lz=Log(z)
      llz=Log(zz)
      parametrizationString=(/&
	    1d0/z,lz/z,&  
	    lz,lz**2,&
	    1d0,zz,z*zz,& 
	    zz*llz/z,&
	    z*lz,z**2*lz,&
	    zz*llz,zz**2*llz,&
	    lz**2*llz,lz*llz**2/)
  
  end function parametrizationString
  
    !!! the function which contains 
    !!! int_z^1 parameterization at values of z -> 1
    !!! it is used to estimate integration error at z~1
  function parametrizationStringAt1(z)
  real*8::z
  real*8,dimension(1:parametrizationLength)::parametrizationStringAt1
  
  parametrizationStringAt1=(/1d0-z, 0d0,0d0,0d0,1d0-z,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
  
  end function parametrizationStringAt1

  
    !!!!Each coefficient is split to delta, sing x->1, regular
    
  !!!!!coefficient function q<-q delta-part
  !!!! NO QUARK HERE!
  function C_q_q_delta(alpha,Nf,Lmu)
  real*8::C_q_q_delta,Nf,alpha,Lmu
  
    C_q_q_delta=0d0
  end function C_q_q_delta
  
  !!!!!coefficient function g<-g delta-part
  !!!! NO DELTA-function!
  function C_g_g_delta(alpha,Nf,Lmu)
  real*8::C_g_g_delta,Nf,alpha,Lmu
  
  C_g_g_delta=0d0
  end function C_g_g_delta
  
  !!!!!coefficient function q<-q singular-part  (1/(1-x)_+,(Log(1-x)/(1-x))_+)
   !!!! NO QUARK HERE!
  subroutine Set_CoeffSing1_q_q(alpha,Nf,Lmu)
  real*8::Nf,alpha,LLambda,Lmu
  
  CoeffSing1_q_q=(/0d0,0d0/)
  
  end subroutine Set_CoeffSing1_q_q
  
  !!!!!coefficient function g<-g singular-part  (1/(1-x)_+,(Log(1-x)/(1-x))_+)
  !!!! NO SINGULAR PART HERE!
  subroutine Set_CoeffSing1_g_g(alpha,Nf,Lmu)
  real*8::Nf,alpha,Lmu
  
  CoeffSing1_g_g=(/0d0,0d0/)
  
  end subroutine Set_CoeffSing1_g_g
  
  !!!!!coefficient function q<-q regular-part
  !!!! NO QUARK HERE!
  subroutine Set_Coeff_q_q(alpha,Nf,Lmu)  
  real*8::alpha,Nf,Lmu
  
  Coeff_q_q=(/0d0, 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
  end subroutine Set_Coeff_q_q
 
  !!!!!coefficient function q<-g regular-part  
  !!!! NO QUARK HERE!
  subroutine Set_Coeff_q_g(alpha,Nf,Lmu)  
  real*8::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_q_g=(/0d0, 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
  end subroutine Set_Coeff_q_g
  
    !!!!!coefficient function g<-q regular-part  
  subroutine Set_Coeff_g_q(alpha,Nf,Lmu)  
  real*8::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_g_q=(/0d0, 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
  
  if(order_global>=1) then
    Coeff_g_q=Coeff_g_q+alpha*(/&
      -16d0/3d0,0d0,&	!1/x, log[x]/x
      0d0,0d0,&		!log[x], log[x]^2
      16d0/3d0,0d0,0d0,&	!1 ,1-x, x(1-x)
      0d0,&		!(1-x)log[1-x]/x
      0d0,0d0,&		!xlog[x], x^2 log[x]
      0d0,0d0,&		!(1-x)Log[1-x],(1-x)^2Log[1-x]
      0d0,0d0/)		!log[1-x]Log[x]^2,log[1-x]^2log[x]
  
  !------The kernels are calculated in mathematica
    if(order_global>=2) then
     Coeff_g_q=Coeff_g_q+alpha*alpha*(/&
      352d0/9d0+128d0/27d0*Nf+Lmu*(32d0/9d0*Nf-512d0/3d0)+80d0*zeta2-128d0*zeta3,& !1/x,
      64d0*(1d0-Lmu),&	! log[x]/x
      1120d0/9d0-448d0/9d0*Lmu,&	!log[x]
      -224d0/9d0,&			!log[x]^2
      -352d0/9d0-128d0/27d0*Nf+Lmu*(-32d0/9d0*Nf+512d0/3d0)-80d0*zeta2+128d0*zeta3,&	!1
      96.93661014992554d0 - 80d0*Lmu/9d0,&	!1-x
      -12.03616313305288d0, &	!x(1-x)
      -70.94828373998075d0+256d0*Lmu/9d0 + 64d0*Nf/9d0,&		!(1-x)log[1-x]/x
      66.13401288059151d0, -4.368597177618905d0,&		!xlog[x], x^2 log[x]
      29.77672764648373d0, -7.277284712894577d0,&		!(1-x)Log[1-x],(1-x)^2Log[1-x]
      -0.21013073998329043d0, 17.922444586392995d0/)		!log[1-x]Log[x]^2,log[1-x]^2log[x]
    end if
    
    
  end if
  end subroutine Set_Coeff_g_q
  
      !!!!!coefficient function g<-g regular-part  
  subroutine Set_Coeff_g_g(alpha,Nf,Lmu)  
  real*8::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_g_g=(/0d0, 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
  if(order_global>=1) then
    Coeff_g_g=Coeff_g_g+alpha*(/&
      -12d0,0d0,&	!1/x, log[x]/x
      0d0,0d0,&		!log[x], log[x]^2
      12d0,0d0,0d0,&	!1 ,1-x, x(1-x)
      0d0,&		!(1-x)log[1-x]/x
      0d0,0d0,&		!xlog[x], x^2 log[x]
      0d0,0d0,&		!(1-x)Log[1-x],(1-x)^2Log[1-x]
      0d0,0d0/)		!log[1-x]Log[x]^2,log[1-x]^2log[x]
  
  !------The kernels are calculated in mathematica
    if(order_global>=2) then
     Coeff_g_g=Coeff_g_g+alpha*alpha*(/&
      220d0+12d0*Nf+Lmu*(64d0/9d0*Nf-408d0)+180d0*zeta2-288d0*zeta3,& !1/x,
      144d0*(1d0-Lmu),&	! log[x]/x
      228d0+8d0*Nf+Lmu*(32d0/3d0*Nf-144d0),&	!log[x]
      -72d0+16d0*Nf/3d0,&			!log[x]^2
      -208d0-16d0*Nf+Lmu*(-64d0/9d0*Nf+408d0)-180d0*zeta2+288d0*zeta3,&	!1
      -1069d0/4d0+80d0/3d0*Nf+Lmu*(64d0/9d0*Nf-48d0)+288d0*zeta2,&	!1-x
      -28d0/3d0*Nf+Lmu*(-32d0/9d0*Nf+24d0)-47.26606760498457d0, &	!x(1-x)
      -0.04368111639126644d0+144d0*Lmu,&		!(1-x)log[1-x]/x
      143.74118698578818d0, -3.936195409223967d0,&		!xlog[x], x^2 log[x]
      18.440046947186577d0, -31.890749138002473d0,&		!(1-x)Log[1-x],(1-x)^2Log[1-x]
      1.080993168990866d0, -1.8412301915911546d0/)		!log[1-x]Log[x]^2,log[1-x]^2log[x]
    end if
    
  end if
  end subroutine Set_Coeff_g_g

  !!!!!coefficient function q<-qb regular-part  
  !!!! NO QUARK HERE!
  subroutine Set_Coeff_q_qb(alpha,Nf,Lmu)  
  real*8::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_q_qb=(/0d0, 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
  end subroutine Set_Coeff_q_qb
  
  !!!!!coefficient function q<-qp regular-part
  !!!! NO QUARK HERE!
  subroutine Set_Coeff_q_qp(alpha,Nf,Lmu)  
  real*8::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_q_qp=(/0d0, 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
  end subroutine Set_Coeff_q_qp
  
  !!! This function has been used during debuging
 subroutine CheckCoefficient(as,Nf,Lmu,z)
 real*8::Lmu,as,z,Nf
 real*8, dimension(1:parametrizationLength)::func
 real*8, dimension(1:2)::func1
 
  func=parametrizationString(z)
     
  func1=(/1d0/(1d0-z),Log(1d0-z)/(1d0-z)/)
 
 !!Q->Q
!   call Set_CoeffSing1_q_q(as,Nf,Lmu) 
!   call Set_Coeff_q_q(as,Nf,Lmu)  
!   write(*,*) SUM(Coeff_q_q*func)+SUM(CoeffSing1_q_q*func1)
  
!   !!Q->G
!   call Set_Coeff_q_g(as,Nf,Lmu)  
!   write(*,*) SUM(Coeff_q_g*func)

!   !!Q->Q'
!   call Set_Coeff_q_qp(as,Nf,Lmu)  
!   write(*,*) SUM(Coeff_q_qp*func)

  !!Q->Qbar
!   call Set_Coeff_q_qb(as,Nf,Lmu)  
!   write(*,*) SUM(Coeff_q_qb*func)
 
!  !! G->Q
!   call Set_Coeff_g_q(as,Nf,Lmu)  
!   write(*,*) SUM(Coeff_g_q*func)

!	!!G->G

  call Set_CoeffSing1_g_g(as,Nf,Lmu) 
  call Set_Coeff_g_g(as,Nf,Lmu)  
!   
!   write(*,*) Coeff_g_g
!   write(*,*) '---------'
!   write(*,*) func
  
  write(*,*) SUM(Coeff_g_g*func)+SUM(CoeffSing1_g_g*func1)
 end subroutine CheckCoefficient

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Convolutions!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------NO SUCH FUNCTION IN HERE
! function lpTMDPDF_lowScale5(x,bT,hadron)
!   real*8,dimension(-5:5)::lpTMDPDF_lowScale5
!   real*8 :: x, bT
!   integer::hadron
!   
!   lpTMDPDF_lowScale5=0d0
!  
!  end function lpTMDPDF_lowScale5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !---------------------------------------------------------------------
!- This is the TMD function evaluated for all quarks simultaniously (-5..5) at x,bT,mu the GLUON INCLUDED
!--    f =  -5, -4, -3,  -2,  -1,0,1,2,3, 4 ,5
!--      = bbar ,cbar sbar,ubar,dbar,g,d,u,s, c ,b
!- The order is accumulative pertrubative order of coefficient =0,1,2 (LO,NLO,NNLO)
!---------------------------------------------------------------------
!---------------------------------------------------------------------
function lpTMDPDF_lowScale50(x,bT,hadron)
  real*8,dimension(-5:5)::lpTMDPDF_lowScale50
  real*8 :: x, bT
  integer::hadron
  
  !!! variables for restoration only
  real*8,dimension(-5:5) :: fNP_grid,fNP_current
  integer::j
  
  if(gridReady .and. ANY(hadronsInGRID.eq.hadron)) then !!! in the case the greed has been calculated
   lpTMDPDF_lowScale50=ExtractFromGrid(x,bT,hadron)
   
   !!!!!!!!!!This is procedure of restoration of function from the initial grid
   !!! if fNP is x-independent then the value can be obtained by TMDPDF(initial) fNP(current)/fNP(initial)
   if(.not.IsFnpZdependent) then
    fNP_grid=FNP(x,0d0,bT,hadron,lambdaNP_grid)
    fNP_current=FNP(x,0d0,bT,hadron,lambdaNP)
    
    do j=-5,5
      if(fNP_grid(j)==0) then
      if(lpTMDPDF_lowScale50(j)/=0.and.j/=0) then
      if(outputlevel>0 .and. messageCounter<messageTrigger) then
	write(*,*)  'WARNING: arTeMiDe',moduleName,' error in restoration: original value is zero. TMDPDF set to zero. b=',bT
	messageCounter=messageCounter+1
	  if(messageCounter>messageTrigger) &
		    write(*,*) 'WARNING: arTeMiDe',moduleName,' number of WARNINGS more then 5. Futher WARNING suppresed'
	end if
      end if
      lpTMDPDF_lowScale50(j)=0
       else
      lpTMDPDF_lowScale50(j)=lpTMDPDF_lowScale50(j)*fNP_current(j)/fNP_grid(j)
      end if
    end do
  end if
   
  else!!!! calculation
   lpTMDPDF_lowScale50=Common_lowScale50(x,bT,hadron)
 end if
  end function lpTMDPDF_lowScale50
  
end module lpTMDPDF