!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 1.4
!
!	Evaluation of the unpolarized TMD PDF at low normalization point in zeta-prescription.
!	
!	if you use this module please, quote 1706.01473
!
!	08.11.2018  essential error in griding with x-dependance FIXED.
!
!				A.Vladimirov (19.04.2018)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module uTMDPDF
use QCDinput
use uTMDPDF_model

implicit none
!------------------------LOCALs -----------------------------------------------

 private 
 
  !Current version of module
 character (len=5),parameter :: version="v1.40"
 character (len=7),parameter :: moduleName="uTMDPDF"
 
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
  ! it is set by SetOrderForC subroutine
  ! 0=LO, 1=NLO, 2=NNLO
  integer :: order_global
  
  integer::lambdaNPlength
  real*8,dimension(:),allocatable::lambdaNP
  real*8,dimension(:),allocatable::lambdaNP_grid !!! this is the value of lambda on which the grid is build

  real*8::c4_global!!!this is the variation constant for mu_OPE
  
    !!!Parameteris of numerics
  real*8 :: tolerance=0.0001d0!!! relative tolerance of the integration
  integer :: maxIteration=5000
  
!------------------------------Variables for coefficient function etc-------------------------------
  
  !!!!!Coefficient lists
  integer,parameter::parametrizationLength=23
  !! { Log[1-x], log[1-x]^2, log[1-x]^3  !exact
  !!   1/x, log[x]/x  !exact
  !! Log[x], log[x]^2, Log[x]^3 !exact
  !! 1 (exact), x, x^2
  !! xLog[x]/(1-x), x Log[x], x^2 Log[x]
  !! x Log[x]^2/(1-x), x Log[x]^2
  !! (Log[x]/(1-x)+1)Log[1-x], Log[x]Log[1-x],  xLog[x]Log[1-x],
  !! (1-x)/x Log[1-x], (1-x)Log[1-x], (1-x)^2 Log[1-x], (1-x)Log^2[1-x] }
  !! The Lmu^2 part is exact the later parts are fitted, but exact if posible (e.g. Lmu and Nf parts for q->q)
  real*8,dimension(1:parametrizationLength) :: Coeff_q_q, Coeff_q_g, Coeff_g_q, Coeff_g_g, Coeff_q_qb, Coeff_q_qp
  !! This is list of coefficeints for the encoding the singular at x->1
  !! { 1/(1-x), (Log[1-x]/(1-x))_+}
  real*8, dimension(1:2) :: CoeffSing1_q_q,CoeffSing1_g_g
  
    !!! trigers for calculation
  logical:: IsMuXdependent, IsFnpZdependent
  
  integer :: counter
  
  INCLUDE 'CommonCode/Twist2Convolution-VAR.f90'
  INCLUDE 'CommonCode/Twist2Grid-VAR.f90'
  
!!--------------------------------- variables for the griding the TMD.---------------------------------------------
  logical :: gridReady!!!!indicator that grid is ready to use. If it is .true., the TMD calculated form the grid
  logical :: withGluon!!!indicator the gluon is needed in the grid
  
  integer::numberOfHadrons
  integer,dimension(:),allocatable::hadronsInGRID

!!-----------------------------------------------Public interface---------------------------------------------------
    
  public::uTMDPDF_Initialize,uTMDPDF_SetLambdaNP,uTMDPDF_SetScaleVariation,uTMDPDF_resetGrid,uTMDPDF_SetPDFreplica
  public::uTMDPDF_lowScale5,uTMDPDF_lowScale50	
!   public::CheckCoefficient  
!   public::mu_OPE
  
  interface uTMDPDF_SetLambdaNP
    module procedure uTMDPDF_SetLambdaNP_usual,uTMDPDF_SetLambdaNP_optional,uTMDPDF_SetReplica,uTMDPDF_SetReplica_optional
  end interface
  
  contains
  
  INCLUDE 'CommonCode/Twist2Convolution.f90'
  INCLUDE 'CommonCode/Twist2Grid.f90'
  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Interface subroutines!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !! Initialization of the package
  subroutine uTMDPDF_Initialize(orderMain)
    character(len=*)::orderMain
    character(256)::line
    integer:: orderIntemidiate,i,j,h
    real*8::dummy1,dummy2
    real*8,dimension(-5:5)::dummyList1,dummyList2
    
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
    
    messageTriger=0
    
    if(outputLevel>1) write(*,*) '----- arTeMiDe.uTMDPDF ',version,': .... initialization'
    
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
	if(outputLevel>0)write(*,*) 'WARNING: arTeMiDe.uTMDPDF_Initialize: unknown order for coefficient function. Switch to NLO.'
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
    !!!!search for uTMDPDF entry
    do
    read(51,'(A)') line
    if(line(1:3)=='*D ') exit
    end do
    read(51,'(A)') line
    read(51,*) tolerance
    read(51,'(A)') line
    read(51,*) maxIteration 
     if(outputLevel>2) then
      write(*,'(A,ES10.3)') ' |  tolerance     =',tolerance
      write(*,'(A,ES10.3)') ' |  max iteration =',REAL(maxIteration)
     end if
    
    !!!!search for grid constants entry
    do
    read(51,'(A)') line    
    if(line(1:3)=='*3 ') exit
    end do    
    !!!!search for uTMDPDF entry
    do
    read(51,'(A)') line
    if(line(1:3)=='*D ') exit
    end do
    read(51,'(A)') line
    read(51,*) xGrid_Min
    read(51,'(A)') line
    read(51,*) bGrid_Max 
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
      write(*,'(A,ES10.3)') ' |  xGrid_Min                 =',xGrid_Min
      write(*,'(A,ES10.3)') ' |  bGrid_Max                 =',bGrid_Max 
      write(*,'(A,I6,A,I6,A)') ' |  (GridSizeX,GridSizeB)     =(',GridSizeX,',',GridSizeB,')'
      write(*,'(A,F6.3)') ' |  slope                     =',slope 
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
    !!!!! PDF entry
    do
    read(51,'(A)') line
    if(line(1:3)=='*1)') exit
    end do
      read(51,*) lambdaNPlength
    
    if(outputLevel>2) write(*,'(A,I3)') ' Number of NP parameters =',lambdaNPlength
    
    if(lambdaNPlength<=0) then
    write(*,*) 'ERROR: arTeMiDe.uTMDPDF_Initialize: number of non-pertrubative &
		    parameters should be >=1. Check the constants-file. Evaluation STOP'
    stop
    end if
    
    allocate(lambdaNP(1:lambdaNPlength))
    allocate(lambdaNP_grid(1:lambdaNPlength))
    
    do i=1,lambdaNPlength
      lambdaNP(i)=1d0
    end do
    
    CLOSE (51, STATUS='KEEP') 
    
    if(QCDinput_IsInitialized()) then
      if(outputLevel>1) write(*,*) 'QCDinput is already initalized'
    else
      if(outputLevel>1) write(*,*) '.. initializing QCDinput'
      call QCDinput_Initialize(orderMain)
      if(outputLevel>1) write(*,*) 'QCDinput is initalized'
    end if
    
    call ModelInit(outputLevel,lambdaNPlength)
    
    !!!!!!!Checking the x-dependance of muOPE
    IsMuXdependent=testMU()
    
    if(IsMuXdependent) then
    if(outputLevel>2) write(*,*) 'arTeMiDe.uTMDPDF: mu OPE is dependent on x'
    else
    if(outputLevel>2) write(*,*) 'arTeMiDe.uTMDPDF: mu OPE is independent on x'
    end if
    
    !!!!!!!Checking the x-dependance of FNP
    IsFnpZdependent=TestFNP(hadronsInGRID,lambdaNPlength)
    
    gridReady=.false.
    
    if(IsFnpZdependent) then
    if(outputLevel>2) write(*,*) 'arTeMiDe.uTMDPDF: FNP is dependent on z'
    else
    if(outputLevel>2) write(*,*) 'arTeMiDe.uTMDPDF: FNP is independent on z'
    end if
    
    if(outputLevel>0) write(*,*) '----- arTeMiDe.uTMDPDF ',version,': .... initialized'
    if(outputLevel>1) write(*,*) ' '
  end subroutine uTMDPDF_Initialize
  
  
  !! call for parameters from the model
  subroutine uTMDPDF_SetReplica(num)
  integer:: num
  
  call uTMDPDF_SetLambdaNP_usual(GiveReplicaParameters(num))
  
  end subroutine uTMDPDF_SetReplica
  
    !! call for parameters from the model
  subroutine uTMDPDF_SetReplica_optional(num,buildGrid, gluonRequared)
  integer:: num
  logical :: buildGrid,gluonRequared
  
  call uTMDPDF_SetLambdaNP_optional(GiveReplicaParameters(num),buildGrid, gluonRequared)
  
  end subroutine uTMDPDF_SetReplica_optional
  
  !! call QCDinput to change the PDF replica number
  !! unset the grid, since it should be recalculated fro different PDF replica.
  subroutine uTMDPDF_SetPDFreplica(rep)
  integer:: rep
  
  call QCDinput_SetPDFreplica(rep)
  gridReady=.false.  
  
  end subroutine uTMDPDF_SetPDFreplica
  
  !!!Sets the non-pertrubative parameters lambda
  subroutine uTMDPDF_SetLambdaNP_usual(lambdaIN)
    real*8,dimension(1:lambdaNPlength)::lambdaIN  
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
    !!! if values are new then grid is out dated, if not do nothing
    if(IsNewValues) then 
     if(outputLevel>2) write(*,*) 'arTeMiDe.',moduleName,': NPparameters reset = (',lambdaNP,')'
     gridReady=.false.
    end if
    
  end subroutine uTMDPDF_SetLambdaNP_usual
  
    !!!Sets the non-pertrubative parameters lambda
    !!! carries additionl option to build the grid
    !!! if need to build grid, specify the gluon requared directive.
  subroutine uTMDPDF_SetLambdaNP_optional(lambdaIN,buildGrid, gluonRequared)
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
    
    
    if(IsNewValues.and.(outputLevel>2)) write(*,*) 'arTeMiDe.',moduleName,': NPparameters reset = (',lambdaNP,')'
    
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
      !!! if grid is not requared nothing to do!    
    end if 
    
  end subroutine uTMDPDF_SetLambdaNP_optional
  
  !!! This subroutine ask for the grid reconstruction (or destruction)
  subroutine uTMDPDF_resetGrid(buildGrid,gluonRequared)
    logical::buildGrid,gluonRequared
    gridReady=.false.
    
!     call uTMDPDF_SetLambdaNP_optional(lambdaNP,buildGrid,gluonRequared)  
      if(buildGrid) then
        if(outputLevel>1) write(*,*) 'arTeMiDe ',moduleName,':  Grid Reset. with c4=',c4_global
	withGluon=gluonRequared
	call MakeGrid()
	gridReady=.true.
      end if
!     write(*,*) 'CHECK'
  end subroutine uTMDPDF_resetGrid
  
  
  !!!! this routine set the variations of scales
  !!!! it is used for the estimation of errors
  subroutine uTMDPDF_SetScaleVariation(c4_in)
    real*8::c4_in
    c4_global=c4_in    
  end subroutine uTMDPDF_SetScaleVariation
 
!-------------------------------------------------
 !!!!array of x times PDF(x,Q) for hadron 'hadron'
 !!!! array is (-5:5) (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
 function xf(x,Q,hadron)
      real*8 :: x,Q
      integer:: hadron
      real*8, dimension(-5:5):: xf
      
      xf=xPDF(x,Q,hadron)
      
  end function xf

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!COEFFICIENT FUNCTIONS!!
  
  !!! the function which contains the functions of parameterizations
  function parametrizationString(z)
  real*8::z,lz,l1z
  real*8,dimension(1:parametrizationLength)::parametrizationString
      lz=Log(z)
      l1z=Log(1d0-z)
      parametrizationString=(/l1z,l1z**2,l1z**3,&
	    1d0/z,lz/z,&  
	    lz,lz**2,lz**3,&
	    1d0,z,z**2,& 
	    z/(1d0-z)*lz, z*lz, (z**2)*lz,& 
	    z/(1d0-z)*lz**2,z*lz**2,& 
	    (lz/(1d0-z)+1d0)*l1z,lz*l1z,z*lz*l1z,& 
	    (1d0-z)/z*l1z, (1d0-z)*l1z, ((1d0-z)**2)*l1z,(1d0-z)*l1z**2/)
  
  end function parametrizationString
  
    !!! the function which contains the functions of parameterizations
    !!! at values of z -> 1
  function parametrizationStringAt1(z)
  real*8::z,l1z
  real*8,dimension(1:parametrizationLength)::parametrizationStringAt1
    l1z=Log(1d0-z)
     parametrizationStringAt1=(/(1d0-z)*(l1z-1),(1d0-z)*(2d0-2d0*l1z+l1z**2),&
	      (1d0-z)*(-6d0+6d0*l1z-3d0*l1z**2+l1z**3),&
	      0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,&
	      ((1d0-z)**2)*l1z,&
	      0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
  
  end function parametrizationStringAt1

  
    !!!!Each coefficient is split to delta, sing x->1, regular
    
  !!!!!coefficient function q<-q delta-part
  function C_q_q_delta(alpha,Nf,Lmu)
  real*8::C_q_q_delta,Nf,alpha,Lmu
  
 C_q_q_delta=1d0
 
  if(order_global>=1) then
      !!!checked (3.12.17) AV!!!
      C_q_q_delta=C_q_q_delta+alpha*(-4d0/3d0*zeta2-4d0*Lmu)
  end if
  if(order_global>=2) then
      !!!checked (3.12.17) AV!!!
     C_q_q_delta=C_q_q_delta+alpha*alpha*(&
     -2416d0/81d0-134d0/3d0*zeta2+448d0/9d0*zeta3+200d0/9d0*zeta4+&
     Nf*(352d0/243d0+20d0/9d0*zeta2+56d0/27d0*zeta3)+&
     Lmu*(-14d0-140d0/3d0*zeta2+16d0/3d0*zeta3+Nf*(4d0/9d0+40d0/9d0*zeta2))+&
     Lmu**2*(-14d0+4d0/3d0*Nf-128d0/9d0*zeta2))     
     !-19.413533328856175d0 + 7.5971237666274805d0*Nf&
     !+ Lmu**2*(-37.394617839619215d0 + 1.3333333333333333d0*Nf)&
     !+ Lmu*(-84.35261963606607d0 + 7.75526251932545d0*Nf))
  end if  
  end function C_q_q_delta
  
  !!!!!coefficient function g<-g delta-part
  function C_g_g_delta(alpha,Nf,Lmu)
  real*8::C_g_g_delta,Nf,alpha,Lmu
  
  C_g_g_delta=1d0
 
  if(order_global>=1) then
  !!!checked (20.11.17) AV!!!
      C_g_g_delta=C_g_g_delta+alpha*(-3d0*zeta2-(11d0-2d0/3d0*Nf)*Lmu)
  end if
  if(order_global>=2) then
     C_g_g_delta=C_g_g_delta+alpha*alpha*(&
    -112d0 - 56d0*Nf**2/81d0 - 201d0*zeta2/2d0 - 72d0*Lmu**2*zeta2 +&
    Lmu*(-96d0 + 32d0*Nf/3d0 - 108d0*zeta3) + Nf*(548d0/27 + 5d0*zeta2 - 28d0*zeta3/3d0) + 154d0*zeta3 + 225d0*zeta4/4d0)
    
    ! C_g_g_delta=C_g_g_delta+alpha*alpha*(&
    !-31.318428735417683d0 - 118.4352528130723d0*Lmu**2 + 17.30176886771455d0*Nf - 0.691358024691358d0*Nf**2&
    !+ Lmu*(-225.82214554123618d0 + 10.666666666666666d0*Nf))
  end if
  
  end function C_g_g_delta
  
  !!!!!coefficient function q<-q singular-part  (1/(1-x)_+,(Log(1-x)/(1-x))_+)
  subroutine Set_CoeffSing1_q_q(alpha,Nf,Lmu)
  real*8::Nf,alpha,LLambda,Lmu,s1,s2
    
  s1=0d0!!!coeff 1/(1-x)
  s2=0d0!!!coeff log(1-x)/(1-x)
  

  if(order_global>=1) then    
    !!!checked (20.11.17) AV!!!
    s1=s1+alpha*(-16d0/3d0)*Lmu
  end if
  if(order_global>=2) then
    !!!checked (21.11.17) AV!!!
    s1=s1+alpha*alpha*(&
       -3232d0/27d0+112d0*zeta3+448d0/81d0*Nf+&
       Lmu*(-1072d0/9d0+352d0/9d0*zeta2+160d0/27d0*Nf)+&
       Lmu**2*(-8d0+16d0/9d0*Nf))
       
!      14.926669450170849d0 + 5.530864197530864d0*Nf&
!      + Lmu**2*(-8d0 + 1.7777777777777777d0*Nf) &
!      + Lmu*(-54.77591205215826d0 + 5.925925925925926d0*Nf))

     s2=s2+alpha*alpha*(256d0/9d0*(Lmu**2))
  end if
  
  CoeffSing1_q_q=(/s1,s2/)
  
  end subroutine Set_CoeffSing1_q_q
  
  !!!!!coefficient function g<-g singular-part  (1/(1-x)_+,(Log(1-x)/(1-x))_+)
  subroutine Set_CoeffSing1_g_g(alpha,Nf,Lmu)
  real*8::Nf,alpha,Lmu,s1,s2
    
  s1=0d0!!!coeff 1/(1-x)
  s2=0d0!!!coeff log(1-x)/(1-x)
 
  if(order_global>=1) then  
  !!!checked (20.11.17) AV!!!
    s1=s1+alpha*(-12d0)*Lmu
  end if
  if(order_global>=2) then
    s1=s1+alpha*alpha*(&
    -808d0/3d0 + Lmu**2*(66d0 - 4d0*Nf) + 112d0*Nf/9d0 + Lmu*(-268d0 + 40d0*Nf/3d0 + 108d0*zeta2) + 252d0*zeta3)
    
!     s1=s1+alpha*alpha*(&
!     33.585006262884406d0 + Lmu**2*(66d0 - 4d0*Nf) + 12.444444444444445d0*Nf&
!     + Lmu*(-90.34712078039155d0 + 13.333333333333334d0*Nf))
    
     s2=s2+alpha*alpha*144d0*(Lmu**2)
  end if
  
  
  CoeffSing1_g_g=(/s1,s2/)
  
  end subroutine Set_CoeffSing1_g_g
  
  !!!!!coefficient function q<-q regular-part  
  subroutine Set_Coeff_q_q(alpha,Nf,Lmu)  
  real*8::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_q_q=(/&
  0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
  0d0,0d0,&     !1/x, Log(x)/x
  0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
  0d0,0d0,0d0,& !1 ,x , x^2
  0d0,0d0,0d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
  0d0,0d0,&     !x/(1-x)Log^2(x), x Log^2(x)
  0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
  0d0,0d0,0d0,0d0/) !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
  if(order_global>=1) then
          
      Coeff_q_q=Coeff_q_q+alpha*(/&
	0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
	0d0,0d0,&     !1/x, Log(x)/x
	0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
	8d0/3d0*(1d0+Lmu),8d0/3d0*(Lmu-1d0),0d0,& !1 ,x , x^2
	0d0,0d0,0d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
	0d0,0d0,&     !x/(1-x)Log^2(x), x Log^2(x)
	0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
	0d0,0d0,0d0,0d0/) !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
  !------The kernels are calculated in mathematica
    if(order_global>=2) then
    
     Coeff_q_q=Coeff_q_q+alpha*alpha*((/&!!!exact part
      200d0/9d0 - 256d0*Lmu/9d0 - 256d0*Lmu**2/9d0, -64d0/9d0, 0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
      0d0,0d0,&     !1/x, Log(x)/x
      -8d0 - 32d0*Lmu**2/9d0 + 40d0*Nf/27d0 + Lmu*(-152d0/9d0 + 16d0*Nf/9d0), &
	-2d0 - 40d0*Lmu/9d0 + 4d0*Nf/9d0, -20d0/27d0,& !Log(x),Log^2(x),Log^3(x)
      Lmu**2*(-28d0/9d0 - 8d0*Nf/9d0) - 296d0*Nf/81d0 + Lmu*(64d0/3d0 - 32d0*Nf/27d0 - 176d0*zeta2/9d0), &
	Lmu**2*(100d0/9d0 - 8d0*Nf/9d0) - 152d0*Nf/81d0 + Lmu*(880d0/9d0 - 128d0*Nf/27d0 - 176d0*zeta2/9d0), 0d0,& !1 ,x , x^2
      -128d0*Lmu**2/9d0 + 80d0*Nf/27d0 + Lmu*(-112d0/3d0 + 32d0*Nf/9d0), &
	32d0*Lmu**2/3d0 + Lmu*(184d0/9d0 - 16d0*Nf/9d0) - 40d0*Nf/27d0, 0d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
      -16d0*Lmu + 8d0*Nf/9d0, 104d0*Lmu/9d0 - 4d0*Nf/9d0,&     !x/(1-x)Log^2(x), x Log^2(x)
      256d0*Lmu/9d0, -128d0*Lmu/9d0, -128d0*Lmu/9d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
      0d0, -128d0*Lmu/9d0 + 128d0*Lmu**2/9d0, 0d0, 0d0/)& !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
    +(/&!!aprroximate 
      0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
      0d0,0d0,&     !1/x, Log(x)/x
      0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
      1076.6744297226016d0, 7792.719665777814d0, 111.49810429898287d0,& !1 ,x , x^2
      8980.334190376141d0, -3795.008745809993d0, 82.30795871692112d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
      -201.0129463471822d0, 206.75145891009598d0,&     !x/(1-x)Log^2(x), x Log^2(x)
      5603.371344939401d0, -526.1352578350587d0, 1382.8610999663256d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
      1092.9256332669593d0, 2547.784733022028d0, -147.17479558391307d0, 3.564983084988843d0/)) !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
    end if
    
!  write(*,*) 'regularPart=', regularPart/x
  end if
  end subroutine Set_Coeff_q_q
 
  !!!!!coefficient function q<-g regular-part  
  subroutine Set_Coeff_q_g(alpha,Nf,Lmu)  
  real*8::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_q_g=(/&
  0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
  0d0,0d0,&     !1/x, Log(x)/x
  0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
  0d0,0d0,0d0,& !1 ,x , x^2
  0d0,0d0,0d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
  0d0,0d0,&     !x/(1-x)Log^2(x), x Log^2(x)
  0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
  0d0,0d0,0d0,0d0/) !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
  if(order_global>=1) then
    Coeff_q_g=Coeff_q_g+alpha*(/ &
      0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
      0d0,0d0,&     !1/x, Log(x)/x
      0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
      -Lmu,2d0+2d0*Lmu,-2d0-2d0*Lmu,& !1 ,x , x^2
      0d0,0d0,0d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
      0d0,0d0,&     !x/(1-x)Log^2(x), x Log^2(x)
      0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
      0d0,0d0,0d0,0d0/) !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
  
    if(order_global>=2) then
     Coeff_q_g=Coeff_q_g+alpha*alpha*((/&!!!Exact part
      -5d0/3d0 + 26d0*Lmu**2/3d0, 10d0*Lmu/3d0, 5d0/9d0,& !log(1-x),log(1-x)^2,log(1-x)^3
      172d0/9d0 - 52d0*Lmu/3d0 + 4d0*Lmu**2 - 8d0*zeta2, 0d0,&     !1/x, Log(x)/x
      34d0/3d0 - 22d0*Lmu/3d0 + 14d0*Lmu**2/3d0, -7d0/6d0 + 14d0*Lmu/3d0, 7d0/9d0,& !Log(x),Log^2(x),Log^3(x)
      7d0*Lmu**2/3d0, 80d0*Lmu**2/3d0, -31d0*Lmu**2,& !1 ,x , x^2
      0d0, 80d0*Lmu**2/3d0, -16d0*Lmu**2/3d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
      0d0,0d0,&     !x/(1-x)Log^2(x), x Log^2(x)
      0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
      0d0,-52d0*Lmu**2/3d0, 52d0*Lmu**2/3d0,0d0/)& !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
      +(/&!!!Approximate part
	0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
	0d0,0d0,&     !1/x, Log(x)/x
	0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
	11804.917158263002d0 - 241.5977366488703d0*Lmu, &
	  43420.91122559322d0 - 1282.2433019065484d0*Lmu, &
	  -3972.3493167408856d0 - 200.38789519756068d0*Lmu,& !1 ,x , x^2
	51274.90217179903d0 - 1729.992223130878d0*Lmu, &
	  -11322.606629412014d0 + 831.0967099095751d0*Lmu, &
	  1490.8267572457225d0 - 43.62058264210284d0*Lmu,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
	365.38010578193774d0 + 17.18487164919196d0*Lmu, &
	  -402.6162297779727d0 - 2.712778270327616d0*Lmu,&     !x/(1-x)Log^2(x), x Log^2(x)
	21027.101170409886d0 - 1172.7172087577621d0*Lmu, &
	  19269.1080641264d0 + 259.31259625889305d0*Lmu, &
	  11357.814388977074d0 - 185.71439659057268d0*Lmu,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
	11798.406220601846d0 - 263.2306304278585d0*Lmu, &
	  29476.38189858057d0 - 264.4846500672005d0*Lmu, &
	  -227.4273932698467d0 + 82.74615814331906d0*Lmu, &
	  15.408922558612357d0 - 6.5920606591307065d0*Lmu/)) !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
    end if
    
  end if
  end subroutine Set_Coeff_q_g
  
    !!!!!coefficient function g<-q regular-part  
  subroutine Set_Coeff_g_q(alpha,Nf,Lmu)  
  real*8::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_g_q=(/&
    0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
    0d0,0d0,&     !1/x, Log(x)/x
    0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
    0d0,0d0,0d0,& !1 ,x , x^2
    0d0,0d0,0d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
    0d0,0d0,&     !x/(1-x)Log^2(x), x Log^2(x)
    0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
    0d0,0d0,0d0,0d0/) !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
  if(order_global>=1) then
    Coeff_g_q=Coeff_g_q+alpha*(/&
      0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
      -16d0/3d0*Lmu,0d0,&     !1/x, Log(x)/x
      0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
      16d0/3d0*Lmu,8d0/3d0*(1d0-Lmu),0d0,& !1 ,x , x^2
      0d0,0d0,0d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
      0d0,0d0,&     !x/(1-x)Log^2(x), x Log^2(x)
      0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
      0d0,0d0,0d0,0d0/) !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
  
  !------The kernels are calculated in mathematica
    if(order_global>=2) then
     Coeff_g_q=Coeff_g_q+alpha*alpha*((/ &!!Exact part
      -184d0/9d0 + 208d0*Lmu**2/9d0 + 32d0*Nf/27d0 + Lmu*(-208d0/3d0 + 32d0*Nf/9d0), &
	 -44d0/9d0 - 80d0*Lmu/9d0 + 8d0*Nf/9d0, -40d0/27d0,& !log(1-x),log(1-x)^2,log(1-x)^3
      -12640d0/27d0 - 248d0*Lmu**2/3d0 + 896d0*Nf/81d0 + 352d0*zeta2/3d0 + Lmu*(-16d0 + 320d0*Nf/27d0 + 16d0*zeta2)&
	  + 192d0*zeta3, -32d0*Lmu**2,&     !1/x, Log(x)/x
      -200d0/3d0 + 448d0*Lmu/9d0 - 224d0*Lmu**2/9d0, 112d0/9d0 - 224d0*Lmu/9d0, -112d0/27d0,& !Log(x),Log^2(x),Log^3(x)
      640d0*Lmu**2/9d0 - 320d0*Lmu*Nf/27d0, 56d0*Lmu**2/9d0 + 208d0*Lmu*Nf/27d0, 32d0*Lmu**2/3d0,& !1 ,x , x^2
      0d0,-320d0*Lmu**2/9d0,0d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
      0d0,0d0,&     !x/(1-x)Log^2(x), x Log^2(x)
      0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
      416d0*Lmu**2/9d0 + 64d0*Lmu*Nf/9d0, -208d0*Lmu**2/9d0 - 32d0*Lmu*Nf/9d0,0d0,0d0/)&
      +(/&!!approximate part
      0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
      0d0,0d0,&     !1/x, Log(x)/x
      0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
      25387.766684267783d0 + 995.0373392459464d0*Lmu - 155.5078210719438d0*Nf, &
	128641.80191589122d0 + 3632.526455640396d0*Lmu - 529.7740230354465d0*Nf, &
	-4691.90206636766d0 - 626.2634197396076d0*Lmu + 33.2493657381115d0*Nf,& !1 ,x , x^2
      149304.47081486747d0 + 4104.585347056393d0*Lmu - 643.5386512087849d0*Nf, &
	-49319.02516115282d0 - 501.0681979811703d0*Lmu + 160.4570922525814d0*Nf, &
	3077.798040737399d0 + 226.19867367258644d0*Lmu - 15.139140368040819d0*Nf,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
      -1297.381684493681d0 - 5.587857511419472d0*Lmu + 2.266359003377411d0*Nf, &
	1379.187692679483d0 - 13.897150833036367d0*Lmu - 2.2865882004537834d0*Nf,&     !x/(1-x)Log^2(x), x Log^2(x)
      79157.51026314059d0 + 1713.027825100937d0*Lmu - 313.45913506070985d0*Nf, &
	20294.80059647656d0 + 1827.284781273226d0*Lmu - 169.33020796238313d0*Nf, &
	27727.02916158908d0 + 870.6474261030947d0*Lmu - 115.18065110064754d0*Nf,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
      25431.533693585556d0 + 1011.076509038249d0*Lmu - 138.52016675095615d0*Nf, &
	62304.61943144707d0 + 2503.327007358015d0*Lmu - 297.70931938249043d0*Nf,&
	-1802.796656932475d0 - 3.5567918590974164d0*Lmu + 4.76862827955143d0*Nf, &
	5.352770139788647d0 - 9.08059951934387d0*Lmu + 0.9080131334622416d0*Nf/))
    end if
    
  end if
  end subroutine Set_Coeff_g_q
  
      !!!!!coefficient function g<-g regular-part  
  subroutine Set_Coeff_g_g(alpha,Nf,Lmu)  
  real*8::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_g_g=(/&
    0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
    0d0,0d0,&     !1/x, Log(x)/x
    0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
    0d0,0d0,0d0,& !1 ,x , x^2
    0d0,0d0,0d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
    0d0,0d0,&     !x/(1-x)Log^2(x), x Log^2(x)
    0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
    0d0,0d0,0d0,0d0/)
  if(order_global>=1) then
    Coeff_g_g=Coeff_g_g+alpha*Lmu*(/&
      0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
      -12d0,0d0,&     !1/x, Log(x)/x
      0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
      24d0,-12d0,12d0,& !1 ,x , x^2
      0d0,0d0,0d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
      0d0,0d0,&     !x/(1-x)Log^2(x), x Log^2(x)
      0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
      0d0,0d0,0d0,0d0/)
  
  !------The kernels are calculated in mathematica
    if(order_global>=2) then
     Coeff_g_g=Coeff_g_g+alpha*alpha*((/ &!!Exact part
      6d0 - 144d0*Lmu - 144d0*Lmu**2 - 2d0*Nf, -36d0, 0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
      -3160d0/3d0 + Lmu**2*(-198d0 - 4d0*Nf/9d0) + 226d0*Nf/9d0 + 264d0*zeta2 + Lmu*(244d0*Nf/9d0 + 36d0*zeta2)&
	+ 432d0*zeta3,-72d0*Lmu**2,&     !1/x, Log(x)/x
      -293d0 + 74d0*Nf/3d0 + Lmu**2*(-72d0 + 16d0*Nf/3d0) + Lmu*(12d0 + 24d0*Nf), &
	3d0 + 6d0*Nf + Lmu*(-72d0 + 16d0*Nf/3d0), -12d0 + 8d0*Nf/9d0,& !Log(x),Log^2(x),Log^3(x)
      -4d0*Lmu*Nf/3d0 + Lmu**2*(84d0 + 32d0*Nf/3d0), Lmu**2*(-150d0 - 20d0*Nf/3d0) - 4d0*Lmu*Nf/3d0, &
	Lmu**2*(198d0 + 4d0*Nf/9d0) - 340d0*Lmu*Nf/9d0,& !1 ,x , x^2
      -72d0*Lmu**2, 24d0*Lmu*Nf + Lmu**2*(-216d0 + 16d0*Nf/3d0), 72d0*Lmu**2,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
      0d0,16d0*Lmu*Nf/3d0,&     !x/(1-x)Log^2(x), x Log^2(x)
      0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
      144d0*Lmu**2, 144d0*Lmu**2, -144d0*Lmu**2,0d0/)&
      +(/&!!approximate part
      0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
      0d0,0d0,&     !1/x, Log(x)/x
      0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
      28350.905925013587d0 - 11826.30667028584d0*Lmu - 1165.597999487027d0*Nf, &
	204106.38400427837d0 - 83819.2024142788d0*Lmu - 9112.918558788755d0*Nf, &
	-3654.296204682914d0 - 4654.703280638828d0*Lmu - 202.82515955533333d0*Nf,& !1 ,x , x^2
      228606.3613205384d0 - 100331.34185957731d0*Lmu - 10439.786162275559d0*Nf, &
	-92515.7827039022d0+ 44024.02341566855d0*Lmu + 4572.107887712472d0*Nf, &
	4598.901824247715d0 + 129.71334826178415d0*Lmu - 82.57482307730146d0*Nf,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
      -5961.667327556318d0 + 539.0169565557716d0*Lmu + 242.87579888233213d0*Nf, &
	6430.443013955799d0 - 686.8969602769131d0*Lmu - 258.29777371318767d0*Nf,&     !x/(1-x)Log^2(x), x Log^2(x)
      144229.2381370356d0 - 56326.8843090263d0*Lmu - 6707.548773597168d0*Nf, &
	-11622.452033932963d0 + 209.23681132222114d0*Lmu + 958.7180437906424d0*Nf, &
	35791.91490687339d0 - 18745.132737698765d0*Lmu - 1542.1541360841081d0*Nf,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
      28788.30003557479d0 - 11546.565659033551d0*Lmu - 1171.3757772648048d0*Nf, &
	67492.96763028382d0 - 35083.678877258644d0*Lmu - 2764.3635575309518d0*Nf, &
	-3881.2814494298837d0 + 1923.6823271680485d0*Lmu + 176.34274530829245d0*Nf, &
	18.352295609756112d0 - 0.22756238102609413d0*Lmu - 0.04234942977167983d0*Nf/))
    end if
    
  end if
  end subroutine Set_Coeff_g_g

   !!!!!coefficient function q<-qb regular-part  
  subroutine Set_Coeff_q_qb(alpha,Nf,Lmu)  
  real*8::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_q_qb=(/&
    0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
    0d0,0d0,&     !1/x, Log(x)/x
    0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
    0d0,0d0,0d0,& !1 ,x , x^2
    0d0,0d0,0d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
    0d0,0d0,&     !x/(1-x)Log^2(x), x Log^2(x)
    0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
    0d0,0d0,0d0,0d0/) !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
  !------The kernels are calculated in mathematica
    if(order_global>=2) then
     Coeff_q_qb=Coeff_q_qb+alpha*alpha*(/&
      0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
      0d0,0d0,&     !1/x, Log(x)/x
      -4d0/3d0 + 16d0*Lmu/9d0, 8d0*Lmu/9d0, 4d0/27d0,& !Log(x),Log^2(x),Log^3(x)
      540.6779037661292d0 + 488.677462854424d0*Lmu, &
	4561.881499265996d0 + 3844.6673748808585d0*Lmu, &
	330.3186826846845d0 + 284.40816788543395d0*Lmu,& !1 ,x , x^2
      5432.878085716809d0 + 4617.753005620716d0*Lmu, &
	-2563.28806902233d0 - 2127.204749929736d0*Lmu, &
	-17.78991469587654d0 - 7.338839370492167d0*Lmu,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
      -76.36755190995123d0 - 27.78047292742369d0*Lmu, &
	78.87763768206408d0 + 27.108672170009157d0*Lmu,&     !x/(1-x)Log^2(x), x Log^2(x)
      3443.1424448947796d0 + 2682.892435111155d0*Lmu, &
	-599.7983485750569d0 - 193.78350958785984d0*Lmu, &
	839.4963238597323d0 + 837.4027040063779d0*Lmu,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
      544.0265746128981d0 + 488.0462345288209d0*Lmu, &
	1417.3624064262308d0 + 1497.1299554721447d0*Lmu, &
	-113.40180701234924d0 - 106.96326421760017d0*Lmu, &
	0.009338040302161874d0 + 0.00862378569953127d0*Lmu/) !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
    end if
  end subroutine Set_Coeff_q_qb
  
     !!!!!coefficient function q<-qp regular-part  
  subroutine Set_Coeff_q_qp(alpha,Nf,Lmu)  
  real*8::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_q_qp=(/&
    0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
    0d0,0d0,&     !1/x, Log(x)/x
    0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
    0d0,0d0,0d0,& !1 ,x , x^2
    0d0,0d0,0d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
    0d0,0d0,&     !x/(1-x)Log^2(x), x Log^2(x)
    0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
    0d0,0d0,0d0,0d0/) !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
  !------The kernels are calculated in mathematica
    if(order_global>=2) then
      Coeff_q_qp=Coeff_q_qp+alpha*alpha*((/ &!!Exact
	0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
	688d0/81d0 - 208d0*Lmu/27d0 + 16d0*Lmu**2/9d0 - 32d0*zeta2/9d0, 0d0,&     !1/x, Log(x)/x
	8d0/3d0 - 8d0*Lmu/3d0 + 8d0*Lmu**2/3d0, -2d0/3d0 + 8d0*Lmu/3d0, 4d0/9d0,& !Log(x),Log^2(x),Log^3(x)
	32d0*Lmu/3d0 + 4d0*Lmu**2/3d0, -16d0*Lmu - 4d0*Lmu**2/3d0, 352d0*Lmu/27d0 - 16d0*Lmu**2/9d0,& !1 ,x , x^2
	0d0, -8d0*Lmu + 8d0*Lmu**2/3d0, -64d0*Lmu/9d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
	0d0, 8d0*Lmu/3d0,&     !x/(1-x)Log^2(x), x Log^2(x)
	0d0,0d0,0d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
	0d0,0d0,0d0,0d0/)& !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
	+(/&!Approximate
	0d0,0d0,0d0,& !log(1-x),log(1-x)^2,log(1-x)^3
	0d0,0d0,&     !1/x, Log(x)/x
	0d0,0d0,0d0,& !Log(x),Log^2(x),Log^3(x)
	-603.924227035499d0, -4636.485211211568d0, -49.76555465398209d0,& !1 ,x , x^2
	-5287.52982020046d0, 2269.612280502602d0, -58.06494427244658d0,& !x/(1-x)Log(x), x Log(x), x^2 Log(x)
	119.75596348356056d0, -129.79997791500733d0,&     !x/(1-x)Log^2(x), x Log^2(x)
	-3369.724995744234d0, 427.8946185110229d0, -812.4665998422748d0,& !(Log(x)/(1-x)+1)Log(1-x) , Log(x)Log(1-x), xLog(x)Log(1-x)
	-600.6972087253562d0, -1469.0061919285804d0, 92.73615445019001d0, -0.021528986881156446d0/)) !(1-x)/xLog(1-x), (1-x)Log(1-x), (1-x)^2 Log(1-x), (1-x)Log^2(1-x)
    end if
  end subroutine Set_Coeff_q_qp
  
  !!! This function has been used during debuging
 subroutine CheckCoefficient(as,Nf,Lmu,z)
 real*8::Lmu,as,z,Nf,lz,l1z
 real*8, dimension(1:23)::func
 real*8, dimension(1:2)::func1
 
  lz=Log(z)
  l1z=Log(1d0-z)
  func=(/l1z,l1z**2,l1z**3,&
	    1d0/z,lz/z,&  
	    lz,lz**2,lz**3,&
	    1d0,z,z**2,& 
	    z/(1d0-z)*lz, z*lz, (z**2)*lz,& 
	    z/(1d0-z)*lz**2,z*lz**2,& 
	    (lz/(1d0-z)+1d0)*l1z,lz*l1z,z*lz*l1z,& 
	    (1d0-z)/z*l1z, (1d0-z)*l1z, ((1d0-z)**2)*l1z,(1d0-z)*l1z**2/)
     
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
  write(*,*) SUM(Coeff_g_g*func)+SUM(CoeffSing1_g_g*func1)
 end subroutine CheckCoefficient

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Convolutions!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
!- This is the TMD function evaluated for all quarks simultaniously (-5..5) at x,bT,mu
!--    f =  -5, -4, -3,  -2,  -1,0,1,2,3, 4 ,5
!--      = bbar ,cbar sbar,ubar,dbar,??,d,u,s, c ,b
!!---Gluon contribution is undefined
!- The order is accumulative pertrubative order of coefficient =0,1,2 (LO,NLO,NNLO)
!---------------------------------------------------------------------
!---------------------------------------------------------------------
function uTMDPDF_lowScale5(x,bT,hadron)
  real*8,dimension(-5:5)::uTMDPDF_lowScale5
  real*8 :: x, bT
  integer::hadron
  
  !!! variables for restoration only
  real*8,dimension(-5:5) :: fNP_grid,fNP_current
  integer::j
  
  
  if(gridReady .and. ANY(hadronsInGRID.eq.hadron)) then !!! in the case the greed has been calculated AND the hadron is in the grid
    
   uTMDPDF_lowScale5=ExtractFromGrid(x,bT,hadron)
   
   !!!!!!!!!!This is procedure of restoration of function from the initial grid
   !!! if fNP is x-independent then the value can be obtained by TMDPDF(initial) fNP(current)/fNP(initial)
   if(.not.IsFnpZdependent) then
    fNP_grid=FNP(x,0d0,bT,hadron,lambdaNP_grid)
    fNP_current=FNP(x,0d0,bT,hadron,lambdaNP)
    
    do j=-5,5
     if(fNP_grid(j)==0) then
      if(fNP_current(j)/=0 .and. ((j/=0).or.(.not.withGluon))) then
      
       if(outputLevel>0 .and. messageTriger<6) then
	  write(*,*) 'WARNING: arTeMiDe.',moduleName,' error in restoration: original value is zero. TMDPDF set to zero. b=',bT
	  messageTriger=messageTriger+1
	  if(messageTriger>5) write(*,*) 'WARNING: arTeMiDe',moduleName,' number of WARNINGS more then 5. Futher WARNING suppresed'
       end if
       end if
      uTMDPDF_lowScale5(j)=0!!!! this is case then 0/0
       else
      uTMDPDF_lowScale5(j)=uTMDPDF_lowScale5(j)*fNP_current(j)/fNP_grid(j)
      end if
    end do
    
  end if
   
  else!!!! calculation
   
  uTMDPDF_lowScale5=Common_lowScale5(x,bT,hadron)
  
  end if
 
 end function uTMDPDF_lowScale5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !---------------------------------------------------------------------
!- This is the TMD function evaluated for all quarks simultaniously (-5..5) at x,bT,mu the GLUON INCLUDED
!--    f =  -5, -4, -3,  -2,  -1,0,1,2,3, 4 ,5
!--      = bbar ,cbar sbar,ubar,dbar,g,d,u,s, c ,b
!- The order is accumulative pertrubative order of coefficient =0,1,2 (LO,NLO,NNLO)
!---------------------------------------------------------------------
!---------------------------------------------------------------------
function uTMDPDF_lowScale50(x,bT,hadron)
  real*8,dimension(-5:5)::uTMDPDF_lowScale50
  real*8 :: x, bT
  integer::hadron
  
  !!! variables for restoration only
  real*8,dimension(-5:5) :: fNP_grid,fNP_current
  integer::j
  
  if(gridReady .and. ANY(hadronsInGRID.eq.hadron)) then !!! in the case the greed has been calculated
   uTMDPDF_lowScale50=ExtractFromGrid(x,bT,hadron)
   
   !!!!!!!!!!This is procedure of restoration of function from the initial grid
   !!! if fNP is x-independent then the value can be obtained by TMDPDF(initial) fNP(current)/fNP(initial)
   if(.not.IsFnpZdependent) then
    fNP_grid=FNP(x,0d0,bT,hadron,lambdaNP_grid)
    fNP_current=FNP(x,0d0,bT,hadron,lambdaNP)
    
    do j=-5,5
      if(fNP_grid(j)==0) then
      if(uTMDPDF_lowScale50(j)/=0.and.j/=0) then
      if(outputlevel>0 .and. messageTriger<6) then
	write(*,*)  'WARNING: arTeMiDe',moduleName,' error in restoration: original value is zero. TMDPDF set to zero. b=',bT
	messageTriger=messageTriger+1
	  if(messageTriger>5) write(*,*) 'WARNING: arTeMiDe',moduleName,' number of WARNINGS more then 5. Futher WARNING suppresed'
	end if
      end if
      uTMDPDF_lowScale50(j)=0
       else
      uTMDPDF_lowScale50(j)=uTMDPDF_lowScale50(j)*fNP_current(j)/fNP_grid(j)
      end if
    end do
  end if
   
  else!!!! calculation
   uTMDPDF_lowScale50=Common_lowScale50(x,bT,hadron)
 end if
  end function uTMDPDF_lowScale50
  
end module uTMDPDF