!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 1.31
!
!	Evaluation of the unpolarized TMD FF at low normalization point in zeta-prescription.
!	
!	if you use this module please, quote 1706.01473 (for arTeMiDe), + ????.???? (for UTMDFF)
!
!
!				A.Vladimirov (19.04.2018)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module uTMDFF
use QCDinput
use uTMDFF_model
implicit none

private

  !Current version of module
 character (len=5),parameter :: version="v1.40"
 character (len=7),parameter :: moduleName="u6TMDFF"
 
  INCLUDE 'Tables/NumConst.f90'  
  INCLUDE 'Tables/G7K15.f90'

!-------------------Physical constants-----------------------------------------
  
  !!!Threashold parameters
  real*8 :: mCHARM=1.4d0
  real*8 :: mBOTTOM=4.75d0

!--------------------------------Working variables-----------------------------------------------
  
  !! Level of output
!! 0=only critical
!! 1=initialization details
!! 2=WARNINGS
  integer::outputLevel=2
    !! variable that count number of WRNING mesagges. In order not to spam too much
  integer::messageTriger=0
  
  !!! The global order which is used in programm
  ! 0=LO, 1=NLO, 2=NNLO
  integer :: order_global
  
  !!parameters of non-perturbative input
  integer::lambdaNPlength
  real*8,dimension(:),allocatable::lambdaNP,lambdaNP_grid
  
  real*8::c4_global!!!this is the variation constant for mu_OPE
  
  !!parameters of numerics
  real*8 :: tolerance=0.0001d0!!! relative tolerance of the integration
  integer :: maxIteration=5000
  
!------------------------------Variables for coefficient function etc-------------------------------
    integer,parameter::parametrizationLength=25
  !!!!!Coefficient lists
  !! { Log[1-z], log[1-z]^2, log[1-z]^3  !exact
  !!   1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
  !! Log[z], log[z]^2, Log[z]^3 !exact
  !! 1 (exact), z, z^2
  !! zLog[z]/(1-z), z Log[z], z^2 Log[z]
  !! z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
  !! (Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  zLog[z]Log[1-z],
  !! (1-z)/z Log[1-z], (1-z)Log[1-z], (1-z) Log[1-z]^2 }
  !! The Lmu^2 part is exact the later parts are fitted, but exact if posible (e.g. Lmu and Nf parts for q->q)
  real*8,dimension(1:parametrizationLength) :: Coeff_q_q, Coeff_q_g, Coeff_g_q, Coeff_g_g, Coeff_q_qb, Coeff_q_qp
  
  !! This is list of coefficeints for the encoding the singular at x->1
  !! { 1/(1-x), (Log[1-x]/(1-x))_+}
  real*8, dimension(1:2) :: CoeffSing1_q_q,CoeffSing1_g_g
  
  !!! trigers for calculation
  logical:: IsMuXdependent, IsFnpZdependent
  integer:: counter
  
  INCLUDE 'CommonCode/Twist2Convolution-VAR.f90'
  INCLUDE 'CommonCode/Twist2Grid-VAR.f90'
  

  
!!--------------------------------- variables for the griding the TMD.---------------------------------------------
  logical :: gridReady!!!!indicator that grid is ready to use. If it is .true., the TMD calculated form the grid
  logical :: withGluon!!!indicator the gluon is needed in the grid
  integer::numberOfHadrons
  integer,dimension(:),allocatable::hadronsInGRID
  
!!-----------------------------------------------Public interface---------------------------------------------------
  public::uTMDFF_Initialize,uTMDFF_SetLambdaNP,uTMDFF_resetGrid,uTMDFF_SetScaleVariation
  public::uTMDFF_lowScale5,uTMDFF_lowScale50
  
  interface uTMDFF_SetLambdaNP
    module procedure uTMDFF_SetLambdaNP_usual,uTMDFF_SetLambdaNP_optional,uTMDFF_SetReplica,uTMDFF_SetReplica_optional
  end interface
  
  contains
  
  INCLUDE 'CommonCode/Twist2Convolution.f90'
  INCLUDE 'CommonCode/Twist2Grid.f90'
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Interface subroutines!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Initialization of the package
  subroutine uTMDFF_Initialize(orderMain)
    character(len=*)::orderMain
    character(256)::line
    integer:: orderIntemidiate,i,j,h
    real*8:: dummy1,dummy2
    real*8,dimension(-5:5):: dummyList1,dummyList2
    
    
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
    
    
    if(outputLevel>1) write(*,*) '----- arTeMiDe.uTMDFF ',version,': .... initialization'
    
    c4_global=1d0
    
    SELECT CASE(orderMain)
      CASE ("LO")
	if(outputLevel>1) write(*,*) 'Order set: LO'
	order_global=0
      CASE ("LO+")
	if(outputLevel>1) write(*,*) 'Order set: LO+'
	order_global=0
      CASE ("NLO")
	if(outputLevel>1) write(*,*) 'Order set: NLO'
	order_global=1
      CASE ("NLO+")
	if(outputLevel>1) write(*,*) 'Order set: NLO+'
	order_global=1
      CASE ("NNLO")
	if(outputLevel>1) write(*,*) 'Order set: NNLO'
	order_global=2
      CASE ("NNLO+")
	if(outputLevel>1) write(*,*) 'Order set: NNLO+'
	order_global=2
      CASE DEFAULT
	if(outputLevel>0)write(*,*) 'WARNING: arTeMiDe.uTMDFF_Initialize: unknown order for coefficient function. Switch to NLO.'
	if(outputLevel>1) write(*,*) 'Order set: NLO'
	order_global=1
     END SELECT
    
     if(outputLevel>2) then
      write(*,'(A,I1)') ' |  Coef.func     =as^',order_global
     end if
     
    !!!!search for physic constants entry
    do
    read(51,'(A)') line    
    if(line(1:3)=='*1 ') exit
    end do    
    !!!!search for masses entry entry
    do
    read(51,'(A)') line
    if(line(1:3)=='*A ') exit
    end do
    read(51,'(A)') line
    read(51,*) mCHARM
    read(51,'(A)') line
    read(51,*) mBOTTOM    
    
    if(outputLevel>2) then
      write(*,*) 'Threshold masses'
      write(*,'(A,F6.3)') ' |  mass Charm    =',mCHARM
      write(*,'(A,F6.3)') ' |  mass Bottom   =',mBOTTOM
     end if
    
    !!!!search for numeric constants entry
    do
    read(51,'(A)') line    
    if(line(1:3)=='*2 ') exit
    end do    
    !!!!search for uTMDFF entry
    do
    read(51,'(A)') line
    if(line(1:3)=='*E ') exit
    end do
    read(51,'(A)') line
    read(51,*) tolerance
    read(51,'(A)') line
    read(51,*) maxIteration
    if(outputLevel>2) then
      write(*,'(A,ES10.3)') ' |  tolerance     =',tolerance
      write(*,'(A,ES10.3)') ' |  max iteration =',REAL(maxIteration)
    end if
    
    !!!search for grid constants entry
    do
    read(51,'(A)') line    
    if(line(1:3)=='*3 ') exit
    end do    
    !!!!search for uTMDFF entry
    do
    read(51,'(A)') line
    if(line(1:3)=='*E ') exit
    end do
    read(51,'(A)') line
    read(51,*) XGrid_Min
    read(51,'(A)') line
    read(51,*) BGrid_Max 
    read(51,'(A)') line
    read(51,*) GridSizeX
    read(51,'(A)') line
    read(51,*) GridSizeB
    read(51,'(A)') line
    read(51,*) slope
    read(51,'(A)') line
    read(51,*) asymptoticBehavior
    read(51,'(A)') line
    read(51,'(A)') line!!! reading hadrons
    !!!! this mini code parse the input to line of integers
    j=1
    do i=1, len(line)    
    if(line(i:i)==',') j=j+1
    end do
        
    !!!! list of hadrons which are to grid.
    !!!! prior to grid building the number is compared and the corresponding grid is build.
    numberOfHadrons=j
    allocate(hadronsInGRID(1:numberOfHadrons))
    read(line,*) hadronsInGRID
    
     if(outputLevel>2) then
      write(*,*) 'Grid options:'
      write(*,'(A,ES10.3)') ' |  zGrid_Min                 =',xGrid_Min
      write(*,'(A,ES10.3)') ' |  bGrid_Max                 =',bGrid_Max 
      write(*,'(A,I6,A,I6,A)') ' |  (GridSizeZ,GridSizeB)     =(',GridSizeX,',',GridSizeB,')'
      write(*,'(A,F6.3)') ' |  Bslope                     =',slope 
      write(*,'(A,I1)') ' |  asymptoticBehavior        =',asymptoticBehavior
      write(*,'(A,I3)')   ' |  hadrons to grid           =',numberOfHadrons
      write(*,*)   ' | list of hadrons in grid    =(',hadronsInGRID,')'
     end if

    allocate(gridMain(0:GridSizeX,0:GridSizeB,-5:5,1:numberOfHadrons))
    allocate(extrapolateParameters(0:GridSizeX,1:2,-5:5,1:numberOfHadrons))
    allocate(taleSings(0:GridSizeX,-5:5,1:numberOfHadrons))
    
    !!!!! number of NP parameters!!!!!!!!!!
    !!!!! NP-parameter section
    do
    read(51,'(A)') line
    if(line(1:3)=='*4 ') exit
    end do
    !!!!! number NP-parameter section
    do
    read(51,'(A)') line
    if(line(1:3)=='*A ') exit
    end do
    !!!!! FF entry
    do
    read(51,'(A)') line
    if(line(1:3)=='*2)') exit
    end do
      read(51,*) lambdaNPlength
      
    if(outputLevel>2) write(*,'(A,I3)') ' Number of NP parameters =',lambdaNPlength
      
    if(lambdaNPlength<=0) then
    write(*,*) 'WARNING: arTeMiDe.uTMDF_Initialize: number of non-pertrubative &
		    parameters should be >=1. Check the constants-file. Evaluation STOP'
    stop
    end if
    
    allocate(lambdaNP(1:lambdaNPlength))
    allocate(lambdaNP_grid(1:lambdaNPlength))
    
    do i=1,lambdaNPlength
      lambdaNP(i)=1.123d0!!!!any number but non-zero, otherwice there is a possible exaption situation during initial grid-evaluation
    end do
    
    CLOSE (51, STATUS='KEEP') 
    
    if(QCDinput_IsInitialized()) then
      if(outputLevel>1) write(*,*) 'QCDinput is already initalized'
    else
      if(outputLevel>1) write(*,*) '.. initializing QCDinput'
      call QCDinput_Initialize(orderMain)
      if(outputLevel>1) write(*,*) 'QCDinput is initalized'
    end if
    
    if(outputLevel>2) write(*,*) 'Model initialization..'
    call ModelInit(outputlevel,lambdaNPlength)
    
    !!!!!!!Checking the z-dependance of muOPE
    IsMuXdependent=testMU()
    
    if(IsMuXdependent) then
    if(outputLevel>2) write(*,*) 'arTeMiDe.uTMDFF: mu OPE is dependent on z'
    else
    if(outputLevel>2) write(*,*) 'arTeMiDe.uTMDFF: mu OPE is independent on z'
    end if
    
    !!!!!!!Checking the z-dependance of FNP
    IsFnpZdependent=TestFNP(hadronsInGRID,lambdaNPlength)
    
    gridReady=.false.
    messageTriger=0
    
    if(IsFnpZdependent) then
    if(outputLevel>2) write(*,*) 'arTeMiDe.uTMDFF: FNP is dependent on z'
    else
    if(outputLevel>2) write(*,*) 'arTeMiDe.uTMDFF: FNP is independent on z'
    end if
     if(outputLevel>0) write(*,*) '----- arTeMiDe.uTMDFF ',version,': .... initialized'
     if(outputLevel>1) write(*,*) ' '
  end subroutine uTMDFF_Initialize
  
    !! call for parameters from the model
  subroutine uTMDFF_SetReplica(num)
  integer:: num
  
  call uTMDFF_SetLambdaNP_usual(GiveReplicaParameters(num))
  
  end subroutine uTMDFF_SetReplica
  
  !! call for parameters from the model
  subroutine uTMDFF_SetReplica_optional(num,buildGrid, gluonRequared)
  integer:: num
  logical::buildGrid, gluonRequared
  
  call uTMDFF_SetLambdaNP_optional(GiveReplicaParameters(num),buildGrid, gluonRequared)
  
  end subroutine uTMDFF_SetReplica_optional

 
  !!!Sets the non-pertrubative parameters lambda
  subroutine uTMDFF_SetLambdaNP_usual(lambdaIN)
    real*8,dimension(1:lambdaNPlength)::lambdaIN
    real*8,dimension(1:lambdaNPlength)::lambdaOLD
    logical::IsNewValues
    integer::i
    lambdaOLD=lambdaNP
    lambdaNP=lambdaIN
    IsNewValues=.false.
    do i=1,lambdaNPlength
     if(ABS(lambdaNP(i)-lambdaOLD(i))>10d-10) then
      IsNewValues=.true.
      exit
     end if
    end do
    !!! if values are new then grid is out dated, if not do nothing
    if(IsNewValues) then 
     if(outputLevel>2) write(*,*) 'arTeMiDe.TMDFF: NPparameters reset = (',lambdaNP,')'
     gridReady=.false.
    end if
  end subroutine uTMDFF_SetLambdaNP_usual
  
    !!!Sets the non-pertrubative parameters lambda
    !!! carries additionl option to build the grid
    !!! if need to build grid, specify the gluon requared directive.
  subroutine uTMDFF_SetLambdaNP_optional(lambdaIN,buildGrid, gluonRequared)
    real*8,dimension(1:lambdaNPlength)::lambdaIN  
    logical :: buildGrid,gluonRequared    
    real*8,dimension(1:lambdaNPlength)::lambdaOLD
    logical::IsNewValues
    integer::i
        messageTriger=0
    
    lambdaOLD=lambdaNP
    lambdaNP=lambdaIN
    IsNewValues=.false.
    do i=1,lambdaNPlength
     if(ABS(lambdaNP(i)-lambdaOLD(i))>10d-10) then
      IsNewValues=.true.
      exit
     end if
    end do    
    
    
    if(IsNewValues.and.(outputLevel>2)) write(*,*) 'arTeMiDe.TMDFF: NPparameters reset = (',lambdaNP,')'
    
    !! further if's are only for griding
    
    if(buildGrid) then !!!grid is requred
    
    if(IsNewValues) then  !! values are new
    
      if(gridReady) then  !!! grid is already build
      
	if(IsFnpZdependent) then !!! check the z-dependance of FNP
	  !! if it is z- dependent, rebuild the grid
	  gridReady=.false.
	  call MakeGrid()
	  gridReady=.true.
	  
	else !!! if z-Independent just do nothing.
	 if(outputLevel>2) write(*,*) 'arTeMiDe_uTMDFF: the values are to be restored from the initial grid'
	end if
      
      else !!! grid is not ready (how comes?)
	call MakeGrid()
	gridReady=.true.
      end if
	
    else  !! values are old
    
      if(gridReady) then  !!! grid is already build
	  !!!nothing to do
      else!!rare option then parameters are not new but grit is not build
        if(outputLevel>2) write(*,*) 'arTeMiDe_uTMDFF: parameters are not reset. But grid is not ready.'
	call MakeGrid()
	gridReady=.true.
      end if
      
    end if
    
    else
      !!! if grid is not requared nothing to do!    
    end if 
    
  end subroutine uTMDFF_SetLambdaNP_optional
  
  !!! This subroutine ask for the grid reconstraction (or distruction)
  subroutine uTMDFF_resetGrid(buildGrid,gluonRequared)
    logical::buildGrid,gluonRequared
    integer::i
    gridReady=.false.
    !call uTMDFF_SetLambdaNP_optional(lambdaNP,buildGrid,gluonRequared)  
    
!     call uTMDPDF_SetLambdaNP_optional(lambdaNP,buildGrid,gluonRequared)  
      if(buildGrid) then
	if(outputLevel>1) write(*,*) 'arTeMiDe.uTMDFF:  Grid Reset. with c4=',c4_global
	withGluon=gluonRequared
! 	call MakeGrid()
! 	gridReady=.true.
      end if
  end subroutine uTMDFF_resetGrid
  
  !!!! this routine set the variations of scales
  !!!! it is used for the estimation of errors
  subroutine uTMDFF_SetScaleVariation(c4_in)
    real*8::c4_in
    c4_global=c4_in    
  end subroutine uTMDFF_SetScaleVariation
 
  function xf(x,Q,hadron)
      real*8 :: x,Q
      integer:: hadron
      real*8, dimension(-5:5):: xf
      
      xf=xFF(x,Q,hadron)
      
  end function xf
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!COEFFICIENT FUNCTIONS!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    
  !!! the function which contains the functions of parameterizations
  function parametrizationString(z)
  real*8::z,lx,l1x
  real*8,dimension(1:parametrizationLength)::parametrizationString
      lx=Log(z)
      l1x=Log(1d0-z)
      parametrizationString=(/ l1x,l1x**2, l1x**3,&
	1d0/z, lx/z, lx**2/z, lx**3/z, &
	lx, lx**2, lx**3,&
	1d0, z, z**2,&
	z*lx/(1d0-z), z*lx, (z**2)*lx,&
	z*(lx**2)/(1d0-z), z*lx**2,  z*lx**3, &
	(lx/(1d0-z)+1d0)*l1x, lx*l1x,  z*lx*l1x,&
	(1d0-z)/z*l1x, (1d0-z)*l1x, (1d0-z)*l1x**2/)/z**2
  
  end function parametrizationString
  
    !!! the function which contains the functions of parameterizations
    !!! at values of z -> 1
  function parametrizationStringAt1(z0)
  real*8::z0
  real*8,dimension(1:parametrizationLength)::parametrizationStringAt1
    parametrizationStringAt1=(/(1d0-z0)*(-1d0+Log(1d0-z0)),(1d0-z0)*(2d0-2d0*Log(1d0-z0)+Log(1d0-z0)**2),&
	      (1d0-z0)*(-6d0+6d0*Log(1d0-z0)-3d0*Log(1d0-z0)**2+Log(1d0-z0)**3),&
	      0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,&
	      (1d0-z0)**2*Log(1d0-z0),0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0/)
  
  end function parametrizationStringAt1

  


!!!!Each coefficient is split to delta, sing x->1, regular
  
  !!!!!coefficient function q<-q delta-part
  function C_q_q_delta(alpha,Nf,Lmu)
  real*8::C_q_q_delta,Nf,alpha,Lmu
  
 C_q_q_delta=1d0
 
  if(order_global>=1) then
      C_q_q_delta=C_q_q_delta+alpha*(-4d0/3d0*zeta2-4d0*Lmu)
  end if
  if(order_global>=2) then
     C_q_q_delta=C_q_q_delta+alpha*alpha*(&
     -2416d0/81d0 + (Lmu**2)*(-14d0 + 4d0*Nf/3d0 - 128d0*zeta2/9d0) - 134d0*zeta2/3d0 + 448d0*zeta3/9d0 &
     + Nf*(352d0/243d0 + 20d0*zeta2/9d0 + 56d0*zeta3/27d0) +  Lmu*(-14d0 - 140d0*zeta2/3d0 + Nf*(4d0/9d0&
     + 40d0*zeta2/9d0) + 16d0*zeta3/3d0) + 2360d0*zeta4/9d0 )
  end if
  
  end function C_q_q_delta
  
  !!!!!coefficient function g<-g delta-part
  function C_g_g_delta(alpha,Nf,Lmu)
  real*8::C_g_g_delta,Nf,alpha,Lmu
  
  C_g_g_delta=1d0
 
  if(order_global>=1) then
   
      C_g_g_delta=C_g_g_delta+alpha*(-3d0*zeta2+(-11d0+2d0/3d0*Nf)*Lmu)
  end if
  if(order_global>=2) then
     C_g_g_delta=C_g_g_delta+alpha*alpha*(&
     -112d0 - 56d0*(Nf**2)/81d0 - 201d0*zeta2/2d0 - 72d0*(Lmu**2)*zeta2 + Lmu*(-96d0 + 32d0*Nf/3d0 &
     - 108d0*zeta3) + Nf*(548d0/27d0 + 5d0*zeta2 - 28d0*zeta3/3d0) + 154d0*zeta3 + 2385d0*zeta4/4d0)
  end if
  end function C_g_g_delta
  
  !!!!!coefficient function q<-q singular-part  (1/(1-x)_+,(Log(1-x)/(1-x))_+)
  subroutine Set_CoeffSing1_q_q(alpha,Nf,Lmu)
  real*8::Nf,alpha,Lmu,s1,s2
    
  s1=0d0!!!coeff 1/(1-x)
  s2=0d0!!!coeff log(1-x)/(1-x)
  

  if(order_global>=1) then    
     
    s1=s1+alpha*(-16d0/3d0*Lmu)
  end if
  if(order_global>=2) then
    s1=s1+alpha*alpha*(&
    -3232d0/27d0 + 448d0*Nf/81d0 + (Lmu**2)*(-8d0 + 16d0*Nf/9d0) + &
    Lmu*(-1072d0/9d0 + 160d0*Nf/27d0 + 352d0*zeta2/9d0) + 112d0*zeta3)
    
     s2=s2+alpha*alpha*(256d0/9d0*Lmu**2)
  end if
  
  CoeffSing1_q_q=(/s1,s2/)
  
  end subroutine Set_CoeffSing1_q_q
  
  !!!!!coefficient function g<-g singular-part  (1/(1-x)_+,(Log(1-x)/(1-x))_+)
  subroutine Set_CoeffSing1_g_g(alpha,Nf,Lmu)
  real*8::Nf,alpha,Lmu,s1,s2
    
  s1=0d0!!!coeff 1/(1-x)
  s2=0d0!!!coeff log(1-x)/(1-x)
 
  if(order_global>=1) then  
   
    s1=s1+alpha*(-12d0)*Lmu
  end if
  if(order_global>=2) then
    s1=s1+alpha*alpha*(&
    -808d0/3d0 + (Lmu**2)*(66d0 - 4d0*Nf) + 112d0*Nf/9d0 +&
      Lmu*(-268d0 + 40d0*Nf/3d0 + 108d0*zeta2) + 252d0*zeta3)
!     
     s2=s2+alpha*alpha*144d0*(Lmu**2)
  end if
  
  
  CoeffSing1_g_g=(/s1,s2/)
  
  end subroutine Set_CoeffSing1_g_g
  
  !!!!!coefficient function q->q
  subroutine Set_Coeff_q_q(alpha,Nf,Lmu)  
  real*8::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_q_q=(/&
  0d0,0d0,0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
  0d0,0d0,0d0,0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
  0d0,0d0,0d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
  0d0,0d0,0d0,&		!1 (exact), z, z^2
  0d0,0d0,0d0,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
  0d0,0d0,0d0,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
  0d0,0d0,0d0,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
  0d0,0d0,0d0/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
  if(order_global>=1) then
    Coeff_q_q=Coeff_q_q+alpha*(/&
      0d0,0d0,0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
      0d0,0d0,0d0,0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
      16d0/3d0,0d0,0d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
      8d0/3d0*(1d0+Lmu),8d0/3d0*(Lmu-1d0), 0d0,&		!1 (exact), z, z^2
      32d0/3d0,-16d0/3d0,0d0,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
      0d0,0d0,0d0,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
      0d0,0d0,0d0,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
      0d0,0d0,0d0/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
    

  !------The kernels are calculated in mathematica
  if(order_global>=2) then
     Coeff_q_q=Coeff_q_q+alpha*alpha*(/&
      -200d0/9d0 + 256d0*Lmu/3d0 - 256d0*(Lmu**2)/9d0, 64d0/9d0, 0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
      0d0,0d0,0d0,0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
      1496d0/9d0-32d0*(Lmu**2)/9d0+Lmu*(248d0/9d0-16d0*Nf/9d0)-8d0*Nf-560d0*zeta2/9d0,&
	  -130d0/9d0+88d0*Lmu/9d0+4d0*Nf/9d0,-140d0/27d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
      (Lmu**2)*(-28d0/9d0 - 8d0*Nf/9d0) - 296d0*Nf/81d0 + Lmu*(64d0/3d0 - 32d0*Nf/27d0 - 176d0*zeta2/9d0),&
	  (Lmu**2)*(100d0/9d0 - 8d0*Nf/9d0) - 152d0*Nf/81d0 + Lmu*(880d0/9d0 - 128d0*Nf/27d0 - 176d0*zeta2/9d0), 0d0,&		!1 (exact), z, z^2
      -128d0*(Lmu**2)/9d0 + Lmu*(-16d0/3d0 - 32d0*Nf/9d0) - 80d0*Nf/9d0, &
	  32d0*(Lmu**2)/3d0 + 8d0*Nf + Lmu*(-24d0 + 16d0*Nf/9d0), 0d0,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
      368d0*Lmu/9d0 + 8d0*Nf/9d0, -280d0*Lmu/9d0 - 4d0*Nf/9d0, 0d0,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
      -256d0*Lmu/3d0, 128d0*Lmu/3d0, 128d0*Lmu/3d0,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
      0d0, -128d0*Lmu/9d0 + 128d0*(Lmu**2)/9d0, 0d0/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
      
      !! adding approximate part
      Coeff_q_q=Coeff_q_q+alpha*alpha*(/&
      0d0,0d0,0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
      0d0,0d0,0d0,0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
      0d0,0d0,0d0,&
      -301.03247439776976d0, -989.2167272286393d0, 82.59056818996338d0,&
      -1063.98482846164d0, 206.28577290245227d0, -18.651271690975136d0,&
      83.00296625888389d0, -70.52319745715631d0, -4.911975080877064d0,&
      -1105.7693500845382d0, 327.376932797376d0, -109.45363459015105d0,&
      -174.6471655693391d0, -112.83673919345797d0, -3.5294575557396084d0/)
      
  end if
    
   !write(*,*) 'regularPart=', regularPart/x
  end if
  end subroutine Set_Coeff_q_q
  
   !!!!!coefficient function q->g
  subroutine Set_Coeff_q_g(alpha,Nf,Lmu)  
  real*8::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_q_g=(/&
  0d0,0d0,0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
  0d0,0d0,0d0,0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
  0d0,0d0,0d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
  0d0,0d0,0d0,&		!1 (exact), z,z^2
  0d0,0d0,0d0,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
  0d0,0d0,0d0,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
  0d0,0d0,0d0,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
  0d0,0d0,0d0/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
  if(order_global>=1) then
    Coeff_q_g=Coeff_q_g+alpha*(/&
      0d0,0d0,0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
      -16d0/3d0*Lmu,32d0/3d0,0d0,0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
      -32d0/3d0,0d0,0d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
      16d0/3d0*Lmu,8d0/3d0*(1d0-Lmu), 0d0,&		!1 (exact), z, z^2
      0d0,16d0/3d0,0d0,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
      0d0,0d0,0d0,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
      0d0,0d0,0d0,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
      0d0,0d0,0d0/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
    

  !------The kernels are calculated in mathematica
    if(order_global>=2) then
      Coeff_q_g=Coeff_q_g+alpha*alpha*(/&
      -40d0/9d0 - 128d0*Lmu/9d0 + 208d0*(Lmu**2)/9d0 - 80d0*zeta2/3d0, -40d0/9d0 + 80d0*Lmu/9d0, 40d0/27d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
      12512d0/27d0 - 248d0*(Lmu**2)/3d0 + Lmu*(112d0/3d0 - 64d0*zeta2) - 352d0*zeta2/3d0 - 352d0*zeta3,&
	  400d0/3d0 + 992d0*Lmu/3d0 - 32d0*Lmu**2, -848d0/3d0 + 128d0*Lmu, -320d0/3d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
      1144d0/3d0 - 64d0*Lmu - 224d0*(Lmu**2)/9d0, 64d0/3d0 + 224d0*Lmu/3d0, -1232d0/27d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
      640d0*(Lmu**2)/9d0, 56d0*(Lmu**2)/9d0, 32d0*(Lmu**2)/3d0,&		!1 (exact), z, z^2
      0d0,-320d0/9d0*Lmu**2, 0d0,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
      0d0, 0d0, 0d0,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
      0d0, 0d0, 0d0,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
      416d0*(Lmu**2)/9d0, -208d0*(Lmu**2)/9d0, 0d0/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
      
      !! adding approximate part
      Coeff_q_g=Coeff_q_g+alpha*alpha*(/&
      0d0,0d0,0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
      0d0,0d0,0d0,0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
      0d0,0d0,0d0,&
      -11519.346897372414d0 - 740.2350001396422d0*Lmu, -34241.466106099186d0 - 1895.8089140328416d0*Lmu,&
		5326.770104932414d0 + 117.62109028406978d0*Lmu,&
      -40601.87518176106d0 - 2457.115336720303d0*Lmu, 4178.463030904903d0 + 197.03426376642062d0*Lmu, &
		-1705.6350033291087d0 - 90.35072880713321d0*Lmu,&
      -966.5754106411847d0 - 65.94510382678313d0*Lmu, 1144.6267544753136d0 + 195.6628827312015d0*Lmu, &
		-84.66732541780037d0 + 0.3444214250305629d0*Lmu,&
      -11393.035115581788d0 - 812.6612648812184d0*Lmu, -25857.49295712562d0 - 1326.6218862880226d0*Lmu, &
		-10601.55795204891d0 - 581.4253894551872d0*Lmu,&
      -12368.214954397781d0 - 628.2350001396422d0*Lmu, -29991.60756399795d0 - 1608.4296229988213d0*Lmu, &
		-5.282535747460972d0 + 9.108108570274817d0*Lmu/)
    end if
    
   !write(*,*) 'regularPart=', regularPart/x
  end if
  end subroutine Set_Coeff_q_g
  
   !!!!!coefficient function g->q
  subroutine Set_Coeff_g_q(alpha,Nf,Lmu)  
  real*8::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_g_q=(/&
  0d0,0d0,0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
  0d0,0d0,0d0,0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
  0d0,0d0,0d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
  0d0,0d0,0d0,&		!1 (exact), z, z^2
  0d0,0d0,0d0,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
  0d0,0d0,0d0,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
  0d0,0d0,0d0,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
  0d0,0d0,0d0/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
  if(order_global>=1) then
    Coeff_g_q=Coeff_g_q+alpha*(/&
      0d0,0d0,0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
      0d0,0d0,0d0,0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
      2d0,0d0,0d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
      -Lmu,2d0*(1d0+Lmu), -2d0*(1d0+Lmu),&		!1 (exact), z,  z^2
      0d0,-4d0,4d0,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
      0d0,0d0,0d0,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
      0d0,0d0,0d0,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
      0d0,0d0,0d0/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
    

  !------The kernels are calculated in mathematica
    if(order_global>=2) then
    Coeff_g_q=Coeff_g_q+alpha*alpha*(/&
      -44d0/3d0 + 26d0*(Lmu**2)/3d0 + 10d0*Nf/9d0 + Lmu*(-14d0 + 4d0*Nf/3d0) + 10d0*zeta2,&
	    -7d0/2d0 - 10d0*Lmu/3d0 + Nf/3d0, -5d0/9d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
      -148d0/9d0 - 4d0*Lmu + 4d0*Lmu**2 + 8d0*zeta2, 8d0 - 16d0*Lmu, 16d0, 0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
      -157d0/3d0 + 10d0*Lmu + 14d0*(Lmu**2)/3d0 - 10d0*Nf/3d0 - 100d0*zeta2/3d0, -73d0/3d0 - 14d0*Lmu + Nf/3d0, 77d0/9d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
      -5042.153939120427d0 - 1531.8031611907763d0*Lmu + 120.72027508333485d0*Lmu**2 &
      - 66.62886623235278d0*Nf + 20.43559788824016d0*Lmu*Nf, &
      303.86339171965324d0 - 12985.247663082215d0*Lmu + 913.3142650671598d0*Lmu**2 &
      - 317.3462434009907d0*Nf + 133.2962117132032d0*Lmu*Nf, &
      3703.409480057778d0 - 412.2734474693585d0*Lmu - 17.358284535835764d0*Lmu**2 &
      + 21.346846003738534d0*Nf + 5.209836558024565d0*Lmu*Nf,&		!1 (exact), z, z^2
      -1059.0639607339292d0 - 14942.360693014281d0*Lmu + 1018.6762556146588*Lmu**2 &
      - 361.41246956998253d0*Nf + 156.7194239372457d0*Lmu*Nf, &
      -8193.44131612091d0 + 6625.5091321556465d0*Lmu - 408.80693588618846d0*Lmu**2 &
      + 116.14537890481989d0*Nf - 66.99593883060648d0*Lmu*Nf, &
      -1046.4745476064809d0 + 142.78460715142586d0*Lmu + 4.0942238266959485d0*Lmu**2 &
      - 18.841069748305667d0*Nf + 1.4503934130147407d0*Lmu*Nf,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
      -1458.9834262434276d0 + 407.23365449967133d0*Lmu - 21.2552050492032d0*Lmu**2 &
      + 5.214301688311585d0*Nf - 3.2700315397596995d0*Lmu*Nf,& 
      1590.9730212100676d0 - 528.3020502418173d0*Lmu + 22.83469849939765d0*Lmu**2 &
      - 6.201059346425027d0*Nf + 3.513030532048041d0*Lmu*Nf, &
      73.07631646309086d0 - 1.3686522702588848d0*Lmu + 0.07440056656025083d0*Lmu**2 &
      - 0.014948523947447324d0*Nf + 0.011446240977035602d0*Lmu*Nf,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
      10420.770382128236d0 - 10241.895319402418d0*Lmu + 638.9211015085228d0*Lmu**2 &
      - 200.7638366136303d0*Nf + 98.29555403916862d0*Lmu*Nf, &
      -20512.88777168424d0 + 2160.2971291961044d0*Lmu - 67.7001647874525d0*Lmu**2 &
      - 49.38978687382705d0*Nf - 10.415409899902482d0*Lmu*Nf, &
      -4683.056024325436d0 - 1584.4478029479442d0*Lmu + 158.259017038288d0*Lmu**2 &
      - 58.30447555159878d0*Nf + 24.347541101116374d0*Lmu*Nf,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
      -4929.425576232434d0 - 1575.4433993058212d0*Lmu + 118.38694175000153d0*Lmu**2 &
      - 65.41307217273041d0*Nf + 18.21337566601794d0*Lmu*Nf, &
      -15198.304463044708d0 - 2938.3075818393627d0*Lmu + 274.35887368634087d0*Lmu**2 &
      - 144.13788267277235d0*Nf + 42.20905754379859d0*Lmu*Nf, &
      -11.384648952243019d0 + 6.516317328981653d0*Lmu + 0.004805761454612681d0*Lmu**2 &
      - 0.6611076391799198d0*Nf + 0.0007393479176680426d0*Lmu*Nf/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
      
    end if
    
   !write(*,*) 'regularPart=', regularPart/x
  end if
  end subroutine Set_Coeff_g_q
  
     !!!!!coefficient function g->g
  subroutine Set_Coeff_g_g(alpha,Nf,Lmu)  
  real*8::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_g_g=(/&
  0d0,0d0,0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
  0d0,0d0,0d0,0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
  0d0,0d0,0d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
  0d0,0d0,0d0,&		!1 (exact), z, z^2
  0d0,0d0,0d0,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
  0d0,0d0,0d0,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
  0d0,0d0,0d0,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
  0d0,0d0,0d0/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
  if(order_global>=1) then
    Coeff_g_g=Coeff_g_g+alpha*(/&
      0d0,0d0,0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
      -12d0*Lmu,24d0,0d0,0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
      -24d0,0d0,0d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
      24d0*Lmu,-12d0*Lmu, 12d0*Lmu,&		!1 (exact), z,  z^2
      24d0,24d0,-24d0,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
      0d0,0d0,0d0,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
      0d0,0d0,0d0,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
      0d0,0d0,0d0/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
    

  !------The kernels are calculated in mathematica
    if(order_global>=2) then
      Coeff_g_g=Coeff_g_g+alpha*alpha*(/&
      -6d0 + 432d0*Lmu - 144d0*Lmu**2 + 2d0*Nf, 36d0, 0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
      3268d0/3d0 + Lmu**2*(-198d0 - 4d0*Nf/9d0) + 278d0*Nf/81d0 - 264d0*zeta2 + Lmu*(412d0*Nf/27d0 + 36d0*zeta2) - 792d0*zeta3,& 
	  536d0 - 72d0*Lmu**2 - 244d0*Nf/9d0 + Lmu*(792d0 + 16d0*Nf/9d0) - 360d0*zeta2,&
	  -660d0 + 288d0*Lmu - 8d0*Nf/9d0, -240d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
      819d0 + Lmu*(-12d0 - 24d0*Nf) - 50d0*Nf + Lmu**2*(-72d0 + 16d0*Nf/3d0) + 360d0*zeta2, &
	  33d0 + Lmu*(216d0 - 16d0*Nf) + 22d0*Nf/3d0, -132d0 + 88d0*Nf/9d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
      -4536.42238065816d0 + 17316.05106770077d0*Lmu - 899.5223166525981d0*Lmu**2 &
	  - 163.95224892156722d0*Nf - 17.333333050953353d0*Lmu*Nf + 10.666666751525197d0*Lmu**2*Nf, &
      18080.936974539778d0 + 118008.31936555808d0*Lmu - 7515.995501919461d0*Lmu**2 &
	  - 625.9193713983447d0*Nf - 6.666666026541869d0*Lmu*Nf - 6.666666410951244d0*Lmu**2*Nf, &
      14139.810571485807d0 - 2400.2117629303043d0*Lmu + 84.6688450985717d0*Lmu**2 &
	  - 184.19497231902477d0*Nf - 4.5925927828634645d0*Lmu*Nf + 0.4444444032182875d0*Lmu**2*Nf,&		!1 (exact), z, z^2
      27937.822349639453d0 + 132629.0291759547d0*Lmu - 8534.848973473487d0*Lmu**2 &
	  - 975.5233827623935d0*Nf + 16.0000007322339d0*Lmu*Nf + 2.993477955391831d-7*Lmu**2*Nf, &
      -35293.57801079519d0 - 52643.48365504077d0*Lmu + 3401.7806883476474d0*Lmu**2 &
	  + 498.8898228158924d0*Nf + 18.66666678376387d0*Lmu*Nf + 5.333333298547101d0*Lmu**2*Nf, &
      -5939.483006875214d0 + 1055.004893054739d0*Lmu - 6.321251200921099d0*Lmu**2 &
	  + 81.15085632871367d0*Nf - 15.999999931023279d0*Lmu*Nf + 1.4357365397679303d-8*Lmu**2*Nf,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
      -6775.737946774227d0 - 4375.194862408091d0*Lmu + 176.58169722761855d0*Lmu**2 &
	  + 30.343871838863492d0*Nf + 5.784855685667105d-8*Lmu*Nf + 5.784855685667105d-9*Lmu**2*Nf, &
      7209.990192354363d0 + 5748.449232000771d0*Lmu - 189.7036422130964d0*Lmu**2 &
	  - 51.50813149633906d0*Nf - 16.00000006398323d0*Lmu*Nf - 6.668502320598985d-9*Lmu**2*Nf, &
      -720.1068548927861d0 + 17.402351899478568d0*Lmu - 0.6180969851274909d0*Lmu**2 &
	  + 9.691339026103627d0*Nf - 2.540110346075181d-10*Lmu*Nf,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
      62244.278348304506d0 + 90595.62066115087d0*Lmu - 5307.959912413799d0*Lmu**2 &
	  - 725.7040599397907d0*Nf - 8.373370761775243d-8*Lmu*Nf + 9.071151658589846d-8*Lmu**2*Nf, &
      -73593.14246325135d0 - 13266.231345168342d0*Lmu + 562.4320363891787d0*Lmu**2 &
	  + 340.92281843959506d0*Nf + 9.023935916755379d-7*Lmu*Nf + 1.7810399835701405d-7*Lmu**2*Nf, &
      -19274.85766344364d0 + 13911.862997386164d0*Lmu - 1314.7672551699827d0*Lmu**2 &
	  + 2.9241023778057107d0*Nf + 3.285034924749133d-7*Lmu*Nf + 7.961647144010051d-8*Lmu**2*Nf,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
      -6776.731748013898d0 + 17595.79207895306d0*Lmu - 839.5223166525981d0*Lmu**2 &
	  - 73.730026699345d0*Nf + 2.823799767667272d-7*Lmu*Nf + 8.485853019485028d-8*Lmu**2*Nf, &
      -54745.95834168829d0 + 28143.140496953594d0*Lmu - 2279.2892078350114d0*Lmu**2 &
	  + 54.675824792462755d0*Nf + 9.417084318872613d-7*Lmu*Nf + 2.2242386654313947d-7*Lmu**2*Nf, &
      -17.510572839038723d0 + 0.9308218455759855d0*Lmu - 0.039924795486343194d0*Lmu**2 - 0.003977517619150036d0*Nf/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
    end if
    
   !write(*,*) 'regularPart=', regularPart/x
  end if
  end subroutine Set_Coeff_g_g
  
     !!!!!coefficient function q->qb
  subroutine Set_Coeff_q_qb(alpha,Nf,Lmu)  
  real*8::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_q_qb=(/&
  0d0,0d0,0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
  0d0,0d0,0d0,0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
  0d0,0d0,0d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
  0d0,0d0,0d0,&		!1 (exact), z, z^2
  0d0,0d0,0d0,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
  0d0,0d0,0d0,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
  0d0,0d0,0d0,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
  0d0,0d0,0d0/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
  if(order_global>=2) then
    Coeff_q_qb=Coeff_q_qb+alpha*alpha*(/&
      0d0, 0d0, 0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
      0d0, 0d0, 0d0, 0d0, &	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
      -76d0/9d0 + 16d0*Lmu/9d0, -32d0/9d0 + 8d0*Lmu/9d0, -4d0/3d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
      -448.1996646927744d0 - 244.21167904318d0*Lmu,&
      -5647.490349524448d0 - 1621.5639861911413d0*Lmu,&
      -257.33492226817515d0 + 202.68655426265354d0*Lmu,&		!1 (exact), z, z^2
      -6353.024936485396d0 - 1663.0891109716677d0*Lmu,&
      3394.3630072013384d0 + 551.6934658836014d0*Lmu,&
      84.62666948172131d0 - 66.45557368691802d0*Lmu,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
      483.57763340201967d0 + 101.50807540427206d0*Lmu,&
      -518.3011843768235d0 - 111.7376149098951d0*Lmu,&
      -0.3446236997186985d0 - 0.44982527645318066d0*Lmu,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
      -6484.394475423013d0 - 1245.2499053204756d0*Lmu,&
      3796.7810643424796d0 + 204.97152039767298d0*Lmu,&
      339.26036997443003d0 - 143.991684454412d0*Lmu,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
      -436.86534933603923d0 - 244.8429073687832d0*Lmu,&
      1330.4397468097895d0 - 317.06568392666475d0*Lmu,&
      -0.02218996661410778d0 - 0.021493311323627588d0*Lmu/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
    end if
    
  end subroutine Set_Coeff_q_qb
  
  !!!!!coefficient function q->qp
  subroutine Set_Coeff_q_qp(alpha,Nf,Lmu)  
  real*8::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_q_qp=(/&
  0d0,0d0,0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
  0d0,0d0,0d0,0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
  0d0,0d0,0d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
  0d0,0d0,0d0,&		!1 (exact), z, z^2
  0d0,0d0,0d0,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
  0d0,0d0,0d0,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
  0d0,0d0,0d0,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
  0d0,0d0,0d0/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
  if(order_global>=2) then
    Coeff_q_qp=Coeff_q_qp+alpha*alpha*(/&
      0d0, 0d0, 0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
      -592d0/81d0 - 16d0*Lmu/9d0 + 16d0*Lmu**2/9d0 + 32d0*zeta2/9d0, 32d0/9d0 - 64d0*Lmu/9d0, 64d0/9d0, 0d0, &	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
      -48d0 + 8d0*Lmu/3d0 + 8d0*Lmu**2/3d0, -22d0/3d0 - 8d0*Lmu, 44d0/9d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
      56d0*Lmu/3d0 + 4d0*Lmu**2/3d0-28.813571629909163d0,&
      -40d0*Lmu/3d0 - 4d0*Lmu**2/3d0+206.17553030550255d0,&
      -32d0*Lmu/9d0 - 16d0*Lmu**2/9d0+76.0334481394452d0,&		!1 (exact), z, z^2
      251.93541929963473d0, 40d0*Lmu/3d0 + 8d0*Lmu**2/3d0-169.05906218222225d0,&
      64d0*Lmu/9d0-32.29013619101719d0,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
      -10.685799944808078d0,-8d0*Lmu+3.7898590887852626d0,4.909581801691148d0,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
      188.58291934528984d0, -90.34188300607897d0, -1.4634823045099683d0,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
      18.62607672661471d0, -16.127886782439663d0, 0.0009797967055855182d0/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
    end if
    
  end subroutine Set_Coeff_q_qp
  
 subroutine CheckCoefficient(as,Nf,Lmu,z)
 real*8::Lmu,as,z,Nf
 real*8, dimension(1:25)::func
 real*8, dimension(1:2)::func1
 
  func=(/ Log(1d0-z),log(1d0-z)**2, log(1d0-z)**3,&
     1d0/z, log(z)/z, Log(z)**2/z, Log(z)**3/z, &
     log(z), log(z)**2, Log(z)**3,&
     1d0, z, z**2,&
     z*Log(z)/(1d0-z), z*Log(z), (z**2)*Log(z),&
     z*(Log(z)**2)/(1d0-z), z*Log(z)**2,  z*Log(z)**3, &
     (Log(z)/(1d0-z)+1d0)*Log(1d0-z), Log(z)*Log(1d0-z),  z*Log(z)*Log(1d0-z),&
     (1d0-z)/z*Log(1d0-z), (1d0-z)*Log(1d0-z), (1d0-z)*Log(1d0-z)**2/)
     
  func1=(/1d0/(1d0-z),Log(1d0-z)/(1d0-z)/)
 
  call Set_Coeff_q_g(as,Nf,Lmu)
  
  write(*,*) SUM(Coeff_q_g*func)
 
 end subroutine CheckCoefficient


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Functions for calculation of convolution!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! We evaluate the integral II=\int_z^1 dy/y  c(y) fNP(y) d(z/y)
!! The coefficient function has 3 parts c(z)=C(z)/z^2=( delta(1-z) Cd+(Cs(z))_++Cr(z) )/z^2,    this last z^2 appears only in FF due to definition of operator
!! The d(z) is given by tables in combination dd(z)=z*d(z).
!! Thus the integral reads II=1/z  \int_z^1 dy/y^2  C(y) fNP(y) dd(z/y)
!! Insrting explicit form we obtain
!! II=1/z  [ Cd fd(z)+\int_z^1 dy/y^2  Cr(y) fNP(y) dd(z/y)+ \int_z^1 dy  Cs(y) (fNP(y) dd(z/y)/y^2-fd(z))-fd(z)\int_0^z Cs(y) dy ]
!! where fd(z)=fNP(1)*dd(z)

!---------------------------------------------------------------------
!- This is the TMD function evaluated for all quarks simultaniously (-5..5) at x,bT,mu
!--    f =  -5, -4, -3,  -2,  -1,0,1,2,3, 4 ,5
!--      = bbar ,cbar sbar,ubar,dbar,g,d,u,s, c ,b
!!---Gluon contribution is undefined
!---------------------------------------------------------------------
!---------------------------------------------------------------------
function uTMDFF_lowScale5(z,bT,hadron)
real*8,dimension(-5:5)::uTMDFF_lowScale5
  real*8 :: z, bT
  integer::hadron,j
  
    !!! variables for restoration only
  real*8,dimension(-5:5) :: fNP_grid,fNP_current
  
  if(gridReady.and. ANY(hadronsInGRID.eq.hadron)) then !!! in the case the greed has been calculated
   uTMDFF_lowScale5=ExtractFromGrid(z,bT,hadron)
   
   !!!!!!!!!!This is procedure of restoration of function from the initial grid
   !!! if fNP is z-independent then the value can be obtained by TMDFF(initial) fNP(current)/fNP(initial)
   if(.not.IsFnpZdependent) then
    fNP_grid=FNP(z,0d0,bT,hadron,lambdaNP_grid)
    fNP_current=FNP(z,0d0,bT,hadron,lambdaNP)
    
    do j=-5,5
      if(fNP_grid(j)==0) then
      if(fNP_current(j)/=0 .and. ((j/=0).or.(.not.withGluon))) then
       if(outputLevel>0 .and. messageTriger<6) then
	  write(*,*) 'WARNING: arTeMiDe_uTMDPDF error in restoration: original value is zero. TMDFF set to zero. b=',bT
	messageTriger=messageTriger+1
	if(messageTriger>5) write(*,*) 'WARNING: arTeMiDe_uTMDFF number of WARNINGS more then 5. Futher WARNING suppresed'
      end if
      end if
      uTMDFF_lowScale5(j)=0
       else
      uTMDFF_lowScale5(j)=uTMDFF_lowScale5(j)*fNP_current(j)/fNP_grid(j)
      end if
    end do
    
  end if
  
  else!!!! calculation
   
  
  uTMDFF_lowScale5=Common_lowScale5(z,bT,hadron)
 end if
 end function uTMDFF_lowScale5

!---------------------------------------------------------------------
!- This is the TMD function evaluated for all quarks simultaniously (-5..5) at x,bT,mu
!--    f =  -5, -4, -3,  -2,  -1,0,1,2,3, 4 ,5
!--      = bbar ,cbar sbar,ubar,dbar,g,d,u,s, c ,b
!---------------------------------------------------------------------
!---------------------------------------------------------------------
function uTMDFF_lowScale50(z,bT,hadron)
real*8,dimension(-5:5)::uTMDFF_lowScale50
  real*8 :: z, bT
  integer::hadron,j
  
  !!! variables for restoration only
  real*8,dimension(-5:5) :: fNP_grid,fNP_current
  
  if(gridReady .and. ANY(hadronsInGRID.eq.hadron)) then !!! in the case the greed has been calculated
  uTMDFF_lowScale50=ExtractFromGrid(z,bT,hadron)
  
  if(.not.IsFnpZdependent) then
    fNP_grid=FNP(z,1d0,bT,hadron,lambdaNP_grid)
    fNP_current=FNP(z,1d0,bT,hadron,lambdaNP)
    
    do j=-5,5
      if(fNP_grid(j)==0) then
      if(uTMDFF_lowScale50(j)/=0) then
      if(outputLevel>0 .and. messageTriger<6) then
	write(*,*)  'WARNING: arTeMiDe_uTMDFF error in restoration: original value is zero. TMDPDF set to zero. b=',bT
	messageTriger=messageTriger+1
	if(messageTriger>5) write(*,*) 'WARNING: arTeMiDe_uTMDFF number of WARNINGS more then 5. Futher WARNING suppresed'
      end if
      end if
      uTMDFF_lowScale50(j)=0
       else
      uTMDFF_lowScale50(j)=uTMDFF_lowScale50(j)*fNP_current(j)/fNP_grid(j)
      end if
    end do
    
  end if
   
   
  else
  uTMDFF_lowScale50=Common_lowScale50(z,bT,hadron)
 end if
  end function uTMDFF_lowScale50

end module