!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 2.00
!
!	Evaluation of the TMD evolution kernel
!	Here we use the improved gamma-solution, and the universal TMD definition.
!	
!	if you use this module please, quote 1803.11089
!
!				A.Vladimirov (17.04.2018)
!			v1.32   A.Vladimirov (30.08.2018)
!				b-freeze at 1d-6 A.Vladimirov (16.09.2018)
!			v1.41   transpose-issue fixed A.Vladimirov (11.03.2019)
!				29.03.2019  Update to version 2.00 (AV).
!			v2.01 	Added zeta-line with non-pertrubative D A.Vladimirov (06.06.2019)
!				Added gluon evolution A.Vladimirov (12.06.2019)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDR
use QCDinput

implicit none

  private
!   public
 
 !Current version of module
 character (len=5),parameter :: version="v2.01"
 character (len=7),parameter :: moduleName="TMDR"
 !Last appropriate verion of constants-file
  integer,parameter::inputver=2
 
!------------------------------------------Physical and mathematical constants------------------------------------------
  
  INCLUDE 'Tables/NumConst.f90'
  
!-------------------------------------------Anomalous dimensions--------------------------------------------------------
  
! These are coefficients of the anomalous dimension D.
! D=dnk(n,k,Nf) as^n Lmu^k 
! So, n is power of as, Nf is power of Log, Nf is number of flavours (=3,4,5)
! Here For quark (last check, AV 16.04.2018)
  real*8,parameter, dimension(1:3,0:3,3:5) :: dnk_q=reshape((/&
  	0.0d0, -15.759631021381722d0, -260.5342569653856d0,& 			!Nf=3,k=0
	2.6666666666666665d0, 24.347721597095045d0, 25.439076311068938d0,&	!Nf=3,k=1
	0.0d0, 12.0d0, 304.4628277071887d0,&					!Nf=3,k=2
	0.0d0,  0.0d0, 72.0d0&							!Nf=3,k=3
	,&
	0.0d0, -18.525063120147152d0, -305.8365169792602d0,&			!Nf=4,k=0
	2.6666666666666665d0, 21.384758634132083d0, -93.99776462634762d0,&	!Nf=4,k=1
	0.0d0, 11.11111111111111d0, 246.65076639554516d0,&			!Nf=4,k=2
	0.0d0,  0.0d0, 61.72839506172839d0&					!Nf=4,k=3
	,&
	0.0d0, -21.290495218912582d0, -342.0455323786133d0,&			!Nf=5,k=0
	2.6666666666666665d0, 18.42179567116912d0, -206.85024342384636d0,&	!Nf=5,k=1
	0.0d0, 10.222222222222221d0, 192.78932236785215d0,&			!Nf=5,k=2
	0.0d0,  0.0d0, 52.24691358024691d0&					!Nf=5,k=3
	/), shape(dnk_q))

! These are coefficients of the anomalous dimension D.
! D=dnk(n,k,Nf) as^n Lmu^k 
! So, n is power of as, Nf is power of Log, Nf is number of flavours (=3,4,5)
! Here For gluon (last check, AV 16.04.2018)
  real*8,parameter, dimension(1:3,0:3,3:5) :: dnk_g=reshape((/&
  	0.d0, -35.459169798108874d0, -586.2020781721177d0,& 			!Nf=3,k=0
	6.d0, 54.78237359346385d0, 57.237921699905115d0,&			!Nf=3,k=1
	0.0d0, 27.d0, 685.0413623411746d0,&					!Nf=3,k=2
	0.0d0,  0.0d0, 162.0d0&							!Nf=3,k=3
	,&
	0.0d0, -41.68139202033109d0, -688.1321632033355d0,&			!Nf=4,k=0
	6.d0, 48.11570692679719d0, -211.49497040928216d0,&			!Nf=4,k=1
	0.0d0, 25.d0, 554.9642243899766d0,&					!Nf=4,k=2
	0.0d0,  0.0d0, 138.88888888888889d0&					!Nf=4,k=3
	,&
	0.0d0, -47.90361424255332d0, -769.60244785188d0,&			!Nf=5,k=0
	6.d0, 41.44904026013052d0, -465.41304770365434d0,&			!Nf=5,k=1
	0.0d0, 23.d0, 433.7759753276673d0,&					!Nf=5,k=2
	0.0d0,  0.0d0, 117.55555555555556d0&					!Nf=5,k=3
	/), shape(dnk_g))
	
! These are beta-function coefficeints at Nf=3,4,5
! (last check, AV 30.08.2018)
  real*8,parameter, dimension(3:5) :: beta0=(/9d0,25d0/3d0,23d0/3d0 /)
  real*8,parameter, dimension(3:5) :: beta1=(/64d0,154d0/3d0,116d0/3d0 /)
  real*8,parameter, dimension(3:5) :: beta2=(/3863d0/6d0,21943d0/54d0,9769d0/54d0/)
  real*8,parameter, dimension(3:5) :: beta3=(/12090.3781d0,8035.18642d0,4826.15633d0/)
	
  !These are gammaV coefficients
  ! gammaV_q(n,Nf) is for as^n, Nf=3,4,5 is number of flavours	
  ! Here for quark (last check, AV 16.04.2018)
  real*8,parameter, dimension(1:3,3:5) :: gammaV_q= reshape((/&
	-8.0d0, -8.0d0, -8.0d0, &
	-29.243530284415503d0, -14.050795508138547d0, 1.1419392681384102d0, &
	-738.2562930508085d0, -491.96573445169145d0, -249.38756710544408d0 /), shape(gammaV_q),ORDER=[2,1])
	
  !These are gammaV coefficients
  ! gammaV_q(n,Nf) is for as^n, Nf=3,4,5 is number of flavours	
  ! Here for gluon (last check, AV 16.04.2018)
  real*8,parameter, dimension(1:3,3:5) :: gammaV_g= reshape((/&
	-18.0d0, -16.666666666666668d0, -15.333333333333334d0, &
	-227.8995118764504d0, -200.7014703660655d0, -173.50342885568062d0, &
	-3957.378545284555d0, -3206.1850399832547d0, -2485.5387880396356d0 /), shape(gammaV_g),ORDER=[2,1])

  !These are gammaCUSP coefficients
  ! gammaCUSP_q(n,Nf) is for as^n, Nf=3,4,5 is number of flavours	
  ! Here for quark (last check, AV 16.04.2018)
  ! The 4-loop gamma cusp is taken from (1808.08981) (3.8)
  real*8,parameter, dimension(1:4,3:5) :: GammaCusp_q= reshape((/&
	5.333333333333333d0, 48.69544319419009d0, 618.2248693918798d0, 7034.85d0, &
	5.333333333333333d0, 42.76951726826417d0, 429.5065747522099d0, 3353.07d0, &
	5.333333333333333d0, 36.84359134233824d0, 239.20803319895987d0, 140.973d0/), shape(GammaCusp_q))
  
  ! These are gammaCUSP coefficients
  ! gammaCUSP_q(n,Nf) is for as^n, Nf=3,4,5 is number of flavours	
  ! Here for gluon (last check, AV 16.04.2018)
  ! The 4-loop gamma cusp is taken from (1808.08981) (3.9) (!!!! it is incomplete)
  real*8,parameter, dimension(1:4,3:5) :: GammaCusp_g= reshape((/&
	12.d0, 109.5647471869277d0, 1391.0059561317298d0, 9882.32d0,&
	12.d0, 96.23141385359438d0, 966.3897931924722d0, 1524.02d0, &
	12.d0, 82.89808052026105d0, 538.2180746976597d0,-5786.52d0 /), shape(GammaCusp_g))
	
	
!------------------------------------------Working variables------------------------------------------------------------

!! Orders of anomalous dimensions as^n
  integer::orderCusp,orderV,orderD,orderDresum,orderZETA
  
!! Level of output
!! 0=only critical
!! 1=initialization details
!! 2=WARNINGS
  integer::outputLevel=2
  integer::messageTrigger=5
  logical::started=.false.
  
!! Precision tolerance used in various routines
  real*8::tolerance=0.001d0

!! Evolution type 
  integer:: EvolutionType=4
  
  integer::counter,messageCounter
  
  
!------------------------------------------Non-pertrubative parameters--------------------------------------------------
  
  !!Number of non-pertrubative parameters
  integer::NPlength=0
  !! array of non-pertrubative parameters
  real*8,allocatable,dimension(:):: NPparam

  public::TMDR_R,TMDR_Rzeta,TMDR_R_toSL
  public:: TMDR_Initialize,TMDR_setNPparameters,LowestQ,TMDR_IsInitialized,TMDR_CurrentNPparameters
  
  public::DNP!,GammaCusp,gammaV
  
  public::zetaMUpert,zetaMuResum,zetaSL,zetaNP
  
  
  interface TMDR_setNPparameters
   module procedure TMDR_setNPparameters, TMDR_SetReplica
  end interface 
  
  interface TMDR_R
    module procedure TMDR_R_type1,TMDR_R_typeZ
  end interface
  
  interface TMDR_Rzeta
    module procedure TMDR_Rzeta_type1,TMDR_Rzeta_typeZ3,TMDR_Rzeta_type3
  end interface
  
 contains
  
   INCLUDE 'Model/TMDR_model.f90'

  function TMDR_IsInitialized()
  logical::TMDR_IsInitialized
  TMDR_IsInitialized=started
  end function TMDR_IsInitialized
   
 !!! move CURRET in streem to the next line that starts from pos (5 char)
 subroutine MoveTO(streem,pos)
 integer,intent(in)::streem
 character(len=5)::pos
 character(len=300)::line
    do
    read(streem,'(A)') line    
    if(line(1:5)==pos) exit
    end do
 end subroutine MoveTO
   
!!! Initializing routing
!!! Filles the prebuiled arrays
!!! orderAD, is order of anomalous dimension for evolution
!!!! order zeta is the order of pertrubative zeta expression, used only in the "universal TMD"
 subroutine TMDR_Initialize(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=300)::path,line
    logical::initRequared
    character(len=8)::orderMain
    integer::i,FILEver
    
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
    
    
    call MoveTO(51,'*3   ')
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
    
    SELECT CASE(orderMain)
      CASE ("LO")
	if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: LO'
	orderCusp=1
	orderV=0
	orderD=0
	orderDresum=0
	orderZETA=0
      CASE ("LO+")
	if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: LO+'
	orderCusp=1
	orderV=1
	orderD=1
	orderDresum=0
	orderZETA=0
      CASE ("NLO")
	if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NLO'
	orderCusp=2
	orderV=1
	orderD=1
	orderDresum=1
	orderZETA=1
      CASE ("NLO+")
	if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NLO+'
	orderCusp=2
	orderV=2
	orderD=2
	orderDresum=1
	orderZETA=1
      CASE ("NNLO")
	if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NNLO'
	orderCusp=3
	orderV=2
	orderD=2
	orderDresum=2
	orderZETA=2
      CASE ("NNLO+")
	if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NNLO+'
	orderCusp=4!3
	orderV=3
	orderD=3
	orderDresum=3!2
	orderZETA=3!2
      CASE ("NNNLO")
	if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NNNLO'
	orderCusp=4
	orderV=3
	orderD=3
	orderDresum=3
	orderZETA=3
      CASE DEFAULT
	if(outputLevel>0) write(*,*) 'arTeMiDe.TMDR_Initialize:try to set unknown ADs orders. Switch to NLO.'
	if(outputLevel>1) write(*,*)  trim(moduleName)//' Order set: NLO'
	orderCusp=2
	orderV=1
	orderD=1
	orderDresum=1
	orderZETA=1
     END SELECT
      
     if(outputLevel>2) then
      write(*,'(A,I1)') ' |  GammaCusp     =as^',orderCusp
      write(*,'(A,I1)') ' |  gammaV        =as^',orderV
      write(*,'(A,I1)') ' |  D             =as^',orderD
      write(*,'(A,I1)') ' |  Dresum        =as^',orderDresum
      write(*,'(A,I1)') ' |  zeta_mu       =as^',orderZETA
     end if

    call MoveTO(51,'*p2  ')
    read(51,*) EvolutionType
    
    if(outputLevel>2) write(*,'(A,I3)') ' Evolution type =',EvolutionType
    
    !-------------------------
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) NPlength
    
    if(outputLevel>2) write(*,'(A,I3)') ' Number of NP parameters =',NPlength
    
    !!allocating number of NP input
    allocate(NPparam(1:NPlength))
    do i=1,NPlength
    NPparam(i)=0d0
    end do
    
    !--------------------------
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) tolerance
    
    if(outputLevel>2) write(*,'(A,ES10.3)') ' |   tolerance=',tolerance
    
    CLOSE (51, STATUS='KEEP') 
    
    
    if(.not.QCDinput_IsInitialized()) then
      if(outputLevel>1) write(*,*) '.. initializing QCDinput (from ',moduleName,')'
      if(present(prefix)) then
      	call QCDinput_Initialize(file,prefix)
      else
	call QCDinput_Initialize(file)
      end if
    end if
    
    if(outputLevel>2) write(*,*) 'Model initialization..'
    call ModelInitialization()
    
    started=.true.
    counter=0
    messageCounter=0
    if(outputLevel>0) write(*,*) '----- arTeMiDe.TMDR ',version,': .... initialized'
    if(outputLevel>1) write(*,*) ' '
  end subroutine TMDR_Initialize

    !Reset the NP parameters
 subroutine TMDR_setNPparameters(lambdaIN)
    real*8,intent(in)::lambdaIN(:)
    integer::ll
    
    ll=size(lambdaIN)
    if(ll<NPlength) then 
      if(outputLevel>0) write(*,"('arTeMiDe.',A,'SetLambdaNP: WARNING length of lambdaNP(,',I3,') is less then requred (',I3,')')")&
	      moduleName,ll,NPlength
      if(outputLevel>0) write(*,*)'                Reset parameters are replaced by zeros!'
      NPparam=0d0*NPparam
      NPparam(1:ll)=lambdaIN(1:ll)
    else if (ll>NPlength) then
      if(outputLevel>0) write(*,"('arTeMiDe.',A,'SetLambdaNP: WARNING&
	      length of lambdaNP(,',I3,') is greater then requred (',I3,')')") moduleName,ll,NPlength
      if(outputLevel>0) write(*,*)'                Array is truncated!'
      NPparam(1:NPlength)=lambdaIN(1:NPlength)
     else
      NPparam=lambdaIN
     end if
    
    if(outputLevel>2) write(*,*) 'arTeMiDe.TMDR: NPparameters reset = (',NPparam,')'
    
 end subroutine TMDR_setNPparameters

 subroutine TMDR_SetReplica(rep)
 integer::rep,i
 real*8,dimension(1:NPlength):: InitVar
   InitVar(1:2)=ReplicaParameters(rep)
    if(NPlength>2) then
     do i=3,NPlength
      InitVar(i)=0d0
     end do
    end if
    
    !!! we also initialize the variables
    call TMDR_setNPparameters(InitVar)
    
 end subroutine TMDR_SetReplica 
 
 subroutine TMDR_CurrentNPparameters(var)
 real*8,dimension(1:NPlength)::var
 var=NPparam
 
 end subroutine TMDR_CurrentNPparameters
!-------------------Anomalous dimensions------------------------------
  
  ! The TMD rapidity anomalous dimension D, pertrubative expression
  function Dpert(mu,bT,f)
    real*8 :: bT,mu,Dpert
    integer:: f,Nf,n,k
    
     ! calculate the boundary value
    real*8 ::LL,astrong,inter
    
    LL=2*LOG(bt*mu*C0_inv_const)
    astrong=As(mu)
    Dpert=0d0
    
    if(mu>=mBOTTOM) then
      Nf=5
    else if(mu>=mCHARM) then
      Nf=4
    else
      Nf=3
    end if
    
    if(f==0) then!gluon
      do n=1,orderD
	inter=0d0
	do k=0,orderD
	 inter=inter+dnk_g(n,k,Nf)*(LL**k)
	end do
	Dpert=Dpert+inter*(astrong**n)
      end do
     else!quark
      do n=1,orderD
	inter=0d0
	do k=0,orderD
	 inter=inter+dnk_q(n,k,Nf)*(LL**k)
	end do
	Dpert=Dpert+inter*(astrong**n)	
      end do
    end if

  end function Dpert
  
  ! The TMD UV anomalous dimension gammaV
  function gammaV(mu,f)
    real*8 :: gammaV,mu
    integer:: f,Nf
    
     ! calculate the boundary value
    real*8 ::astrong
    
    astrong=As(mu)
    
    if(mu>=mBOTTOM) then
      Nf=5
    else if(mu>=mCHARM) then
      Nf=4
    else
      Nf=3
    end if
   
   if(f==0) then!gluon case
    if(orderV>=1) then
    gammaV=astrong*gammaV_g(1,Nf)
     if(orderV>=2) then
      gammaV=gammaV+gammaV_g(2,Nf)*astrong**2     
      if(orderV>=3) then
	gammaV=gammaV+gammaV_g(3,Nf)*astrong**3
      end if
     end if
    else
    gammaV=0
    end if
   else !! quark case
    if(orderV>=1) then
    gammaV=astrong*gammaV_q(1,Nf)
     if(orderV>=2) then
      gammaV=gammaV+gammaV_q(2,Nf)*astrong**2     
      if(orderV>=3) then
	gammaV=gammaV+gammaV_q(3,Nf)*astrong**3
      end if
     end if
    else
    gammaV=0
    end if
   end if
   
  end function gammaV

  ! The Gamma Cusp anomalous dimension gammaV
  function gammaCUSP(mu,f)
    real*8 :: gammaCUSP,mu
    integer:: f,Nf
    
     ! calculate the boundary value
    real*8 ::astrong
    
    astrong=As(mu)
    
    if(mu>=mBOTTOM) then
      Nf=5
    else if(mu>=mCHARM) then
      Nf=4
    else
      Nf=3
    end if
   
   if(f==0) then!gluon case
    gammaCUSP=astrong*GammaCusp_g(1,Nf)
     if(orderCusp>=2) then
      gammaCUSP=gammaCUSP+GammaCusp_g(2,Nf)*astrong**2     
      if(orderCusp>=3) then
	gammaCUSP=gammaCUSP+GammaCusp_g(3,Nf)*astrong**3
	if(orderCusp>=4) then
	 write(*,*) 'ERROR: arTeMiDe.TMDR. Gamma CUSP at 4-loops for gluon is unknown'
         write(*,*) 'Evaluation stop'
	 stop
	 gammaCUSP=gammaCUSP+GammaCusp_g(4,Nf)*astrong**4
	 end if
      end if
     end if
   else !! quark case
    gammaCUSP=astrong*GammaCusp_q(1,Nf)
     if(orderCusp>=2) then
      gammaCUSP=gammaCUSP+GammaCusp_q(2,Nf)*astrong**2     
      if(orderCusp>=3) then
	gammaCUSP=gammaCUSP+GammaCusp_q(3,Nf)*astrong**3
	if(orderCusp>=4) then
	  gammaCUSP=gammaCUSP+GammaCusp_q(4,Nf)*astrong**4
	end if
      end if
     end if
   end if
   
   
  end function gammaCUSP
  
  !the resummed version of rapidity anomalous dimension taken from 1208.1281 (0,1,2-loop)
! 3-loops calculated in Math file.
! IT has natural counting
 function Dresum(mu,b,f)
    real*8::Dresum,mu,b
    integer::f
    real*8:: X,alpha,lX
    
    alpha=As(mu)
    If(mu<mCHARM) then !! nf=3
    !X=1-b0*Lmu
    X=18d0*alpha*Log(mu*b*C0_inv_const)
    lX=Log(1d0-X)
    
    !last check (30.09.2018)+(12.06.2019)
    Dresum=-8d0/27d0*lX
    
    if(orderDresum>=1) then
    !last check (30.09.2018)+(12.06.2019)
    Dresum=Dresum+alpha/(1d0-X)*(0.5983065149035645d0*X - 512d0*lX/243d0)    
      if(orderDresum>=2) then
	!last check (30.09.2018)+(12.06.2019)
        Dresum=Dresum+((alpha/(1d0-X))**2)*(&
        -15.75963102138173d0 - 19.23770595326028d0*lX +16384d0*lX**2/2187d0&
        + 15.108120124066376d0*X - 4.447521424630353d0*X**2)
        !!!!!4-loop
        if(orderDresum>=3) then
	  !last check (30.09.2018) error?
! 	  Dresum=Dresum+((alpha/(1d0-X))**3)*(&
! 	  -260.5342569653865d0 - 20.10001091244974d0*lX + 190.07464445795534d0*lX**2 &
! 	  -35.51545326762519d0*lX**3 + 146.58790833901065d0*X - 44.18188284306254d0*lX*X &
! 	  +11.703452376448753d0*X**2 - 37.75565322877677d0*X**3)

	  !(12.06.2019) difference in polynomial part
	  Dresum=Dresum+((alpha/(1d0-X))**3)*(&
	  -260.5342569653858d0 - 20.10001091244908d0*lX + 190.07464445795534d0*lX**2 &
	  - 35.51545326762519d0*lX**3 + 191.81857011678846d0*X - 44.18188284306254d0*lX*X &
	  - 33.52720940132923d0*X**2 - 22.678765969517485d0*X**3)
	end if
      end if
    end if
    
    else if(mu<mBOTTOM) then !! nf=4
    X=50d0/3d0*alpha*Log(mu*b*C0_inv_const)
    lX=Log(1d0-X)
    
    !last check (30.09.2018)+(12.06.2019)
    Dresum=-8d0/25d0*lX
    if(orderDresum>=1) then
        !last check (30.09.2018)+(12.06.2019)
        Dresum=Dresum+alpha/(1d0-X)*(0.5949710360958496d0*X - 1232d0*lX/625d0)
    if(orderDresum>=2) then
	!last check (30.09.2018)+(12.06.2019)
        Dresum=Dresum+((alpha/(1d0-X))**2)*(&
        -18.52506312014714d0 - 15.807613582350436d0*lX + 94864d0*lX**2/15625d0 &
        + 9.962780902782157d0*X - 3.2507308958355225d0*X**2)
        !!!!!4-loop
        if(orderDresum>=3) then
	  !last check (30.09.2018) error?
! 	  Dresum=Dresum+((alpha/(1d0-X))**3)*(&
! 	  -305.836516979261d0 + 69.4831476117966d0*lX + 134.7740830272787d0*lX**2 &
! 	  -24.93278890666667d0*lX**3 + 42.43828293158347d0*X - 21.321725724444445d0*lX*X &
! 	  +80.87368989899981d0*X**2 - 58.80945335817513d0*X**3)

	  !(12.06.2019) difference in polynomial part
	  Dresum=Dresum+((alpha/(1d0-X))**3)*(&
	  -305.83651697926075d0 + 69.48314761179631d0*lX + 134.77408302727866d0*lX**2 &
	  - 24.932788906666666d0*lX**3 + 100.03436997158337d0*X - 21.321725724444445d0*lX*X &
	  + 23.277602858999842d0*X**2 - 39.61075767817519d0*X**3)
	end if
    end if
    end if
    else !!! nf =5 
    X=46d0/3d0*alpha*Log(mu*b*C0_inv_const)
    lX=Log(1d0-X)
    
    !last check (12.06.2019)
    Dresum=-8d0/23d0*lX
    if(orderDresum>=1) then
	!last check (30.09.2018)+(12.06.2019)
        Dresum=Dresum+alpha/(1d0-X)*(0.6485896055022099d0*X - 928d0*lX/529d0)
     if(orderDresum>=2) then
	!last check (30.09.2018)+(12.06.2019)
        Dresum=Dresum+((alpha/(1d0-X))**2)*(&
        -21.2904952189126d0 - 12.118685999181197d0*lX + 53824d0*lX**2/12167d0&
        + 3.4818379050987946d0*X - 2.0609284500060885d0*X**2)
        !!!!4-loop
        if(orderDresum>=3) then
	  !last check (30.09.2018) error??
! 	  Dresum=Dresum+((alpha/(1d0-X))**3)*(&
! 	  -342.0455323786137d0 + 136.0753964300541d0*lX + 83.43151323800433d0*lX**2 &
! 	  -14.874122567219718d0*lX**3 - 69.48701186506418d0*X + 3.227921887389228d0*lX*X &
! 	  +155.46041408237463d0*X**2 - 81.95172616030685d0*X**3)

	  !(12.06.2019) difference in polynomial part
	  Dresum=Dresum+((alpha/(1d0-X))**3)*(&
	  -342.04553237861387d0 + 136.07539643005407d0*lX + 83.4315132380043d0*lX**2 &
	  - 14.874122567219718d0*lX**3 + 22.601705526240153d0*X + 3.2279218873892277d0*lX*X &
	  + 63.371696691070326d0*X**2 - 51.255487029872086d0*X**3)
	end if
        
    end if
    end if
    end if
    
    !! gluon case is (quark case)*9/4 (at all orders)
    if(f==0) then
     Dresum=Dresum*9d0/4d0
    end if
    
    if(ISNAN(Dresum)) then
        write(*,*) 'ERROR: arTeMiDe.TMDR: Dresum is NaN. At mu=',mu,'b=',b
        write(*,*) 'Lmu=',2*Log(mu*b*C0_inv_const), 'X=',X,'log(1-x)=',lX
        write(*,*) 'Evaluation STOP'
        stop
    end if
    

!     write(*,*) X,Dresum
 end function Dresum

!-------------------zeta-lines -------------------------------------

   !!! the value of zeta_mu in the pertrubation theory with ri=0
 function zetaMUpert(mu,bt,f)
  real*8:: zetaMUpert,mu,bt,alpha,Lmu,val
  integer::f
  
   val=0d0
      if(f==0) then !gluon
      
      if(mu>mBOTTOM) then !Nf=5	
	val=-23d0/18d0
	if(orderZETA>=1) then !!NLO+
	 Lmu=2*LOG(mu*bt*C0_inv_const)
	 alpha=As(mu)
	 val=val+alpha*(23d0/36d0*Lmu**2-9.62347649908429d0)      
	 if(orderZETA>=2) then !!NNLO+	
	  val=val+alpha**2*(529d0/216d0*Lmu**3+7.63577743510649d0*Lmu**2-62.39423231081216d0*Lmu-147.47086756154818d0)
	 end if
	end if
      
      else if(mu>mCHARM) then !Nf=4
	val=-25d0/18d0
	if(orderZETA>=1) then !!NLO+
	 Lmu=2*LOG(mu*bt*C0_inv_const)
	 alpha=As(mu)
	 val=val+alpha*(25d0/36d0*Lmu**2-9.06067674358926d0)
	 if(orderZETA>=2) then !!NNLO+	
	  val=val+alpha**2*(625d0/216d0*Lmu**3+9.846725338749675d0*Lmu**2-61.873995702646845d0*Lmu-140.01561863166336d0)
	 end if
	end if
      
      else !Nf=3
      	val=-1.5d0
	if(orderZETA>=1) then !!NLO+
	 Lmu=2*LOG(mu*bt*C0_inv_const)
	 alpha=As(mu)
	 val=val+alpha*(3d0/4d0*Lmu**2-8.2509634078473d0)
	 if(orderZETA>=2) then !!NNLO+	
	  val=val+alpha**2*(27d0/8d0*Lmu**3+12.181130032516315d0*Lmu**2-59.29131645913094d0*Lmu-129.42141411914173d0)
	 end if
	end if
      end if
      
      else !quark
	
	val=-1.5d0
	
	if(orderZETA>=1) then
	 Lmu=2*LOG(mu*bt*C0_inv_const)
	 alpha=As(mu)
		
	 if(mu>mBOTTOM) then !Nf=5
	  val=val+alpha*(23d0/36d0*Lmu**2+6.58440582426247d0)
	  if(orderZETA>=2) then !!NNLO+	
	    val=val+alpha**2*(529d0/216d0*Lmu**3+7.63577743510649d0*Lmu**2-0.26401673798290304d0*Lmu-89.10266383253594d0)
	  end if
	  
	else if(mu>mCHARM) then !Nf=4
	 val=val+alpha*(25d0/36d0*Lmu**2+5.920953238895728d0)
	  if(orderZETA>=2) then !!NNLO+	
	   val=val+alpha**2*(625d0/216d0*Lmu**3+9.846725338749675d0*Lmu**2+0.5494625577072547d0*Lmu-76.27100645590781d0)
	  end if
	else !Nf=3
	 val=val+alpha*(3d0/4d0*Lmu**2+5.25750065352898d0)
	  if(orderZETA>=2) then !!NNLO+	
	   val=val+alpha**2*(27d0/8d0*Lmu**3+12.181130032516315d0*Lmu**2 + 1.4967718170623598d0*Lmu-61.40054443982149d0)	
	  end if
	end if
      end if
      end if
      
      zetaMUpert=mu*C0_const/bT*EXP(-val)
 end function zetaMUpert

 !!!! the value of zeta-line resummed over log[b]
 function zetaMUresum(mu,b,f)
  real*8::zetaMUresum,mu,b
  integer::f
  real*8::X,val,lX,alpha,XlX

  !val=g_0+a g_1+...
  
   alpha=As(mu)
   
   !----------------------------------------------------------QUARK----------------------------------------
    If(mu<mCHARM) then !! ----------------------------------------nf=3
    X=18d0*alpha*Log(mu*b*C0_inv_const)
    lX=Log(1d0-X)
    XlX=X/((1d0-X)*lX)
    
      !last check (30.09.2018)
    val=(X+lX)/lX/9d0
    
    if(orderZETA>=1) then
    
    !last check (30.09.2018)
    val=val+alpha*(32d0/81d0*lX-64d0/81d0*XlX-2.5144883998789602d0+0.2243649430888367d0*X*XlX/lX)
    
    end if
    if(orderZETA>=2) then
!     !last check (30.09.2018) ?Error?
!     val=val+(alpha**2)*(&
!     (-2.809327846364883d0*lX)/(X-1d0)&
!     +(0.7344199383092869d0 + X*(0.5427459498748931d0 +1.5321619581807067d0*X))/(X-1d0)**2&
!     +(2.4745769670392037d0 + lX*(3.8610071741320366d0 + (-5.9098616330181475d0+0.4533704191432824d0*X)/X**2)&
!       - 5.9098616330181475d0/X - 1.6678205342363737d0*X + 0.4530566491853122d0*X/lX)*XlX**2)
      
    !(12.06.2019)
    val=val+(alpha**2)*(&
    (-2.809327846364883*lX)/(-1 + X) + (0.734419938309286 + 0.54274594987489*X + 1.5321619581807064*X**2)/(-1 + X)**2 &
    + ((5.909861633018146 - 0.45337041914327586*X - 3.86100717413203*X**2)*XlX)/((-1 + X)*X) &
    + (2.474576967039215 - 5.909861633018147/X - 1.6678205342363828*X + (0.4530566491853128*X)/lX)*XlX**2)
    
    end if
    if(orderZETA>=3) then
    !last check (30.09.2018) ?Error?
!     val=val+(alpha**3)*(&
!     (-9.988721231519586d0*lX**2)/(-1d0 + X)**2 &
!     +(lX*(-25.650274604347043d0 + 38.96856957970649d0*X))/(-1d0 + X)**3 &
!     +(-116.72083321043493d0 + 157.99440284791794d0*X - 104.75296305816975*X**2 &
!     +32.15628667503184d0*X**3)/(-1d0 + X)**3 &
!     +(13.21548593880104d0 + lX*(35.90538079987103d0 - 25.582672835245578d0/X**2 - 16.566763069323986d0/X &
!     -14.158369960791271d0*X) - 23.8673838411903d0/X - 6.7356082664341095d0*X &
!     +(0.9148502637943364d0*X)/lX + (lX**2*(-1.7152889940552896d0 + 138.17793148235785d0*X &
!     -113.52541590924972d0*X**2 + 26.948941036555656d0*X**3))/X**3)*XlX**3)

    !(12.06.2019)
    val=val+(alpha**3)*(&
    (-9.988721231519586*lX**2)/(-1 + X)**2 - (0.9148502637943364*X**4)/(lX**4*(-1 + X)**3) &
    + (lX*(-25.650274604347043 + 38.96856957970649*X))/(-1 + X)**3 &
    + (6.7356082664341095*X**2*(3.543463767055726 - 1.9620330363715506*X + X**2))/(lX**3*(-1 + X)**3) &
    + (-116.72083321043495 + 157.9944028479179*X - 104.75296305816973*X**2 + 32.156286675031836*X**3)/(-1 + X)**3 &
    + ((-1.7152889940552796 + 155.13942964902444*X - 138.9676631592496*X**2 + 35.42969011988898*X**3)*XlX)/((-1 + X)**2*X) &
    + ((25.582672835245585 - 0.3947350973427885*X - 18.94388263320436*X**2 + 8.504537238569082*X**3)*XlX**2)/((-1 + X)*X))
    
    end if
    
    else if(mu<mBOTTOM) then !! --------------------------------------- nf=4
    X=50d0/3d0*alpha*Log(mu*b*C0_inv_const)
    lX=Log(1d0-X)
    XlX=X/((1d0-X)*lX)
    
    !last check (30.09.2018)
    val=3d0*(X+lX)/lX/25d0
    
    if(orderZETA>=1) then
    
    !last check (30.09.2018)
    val=val+alpha*(-2.4623141385359437d0+231d0*lX/625d0-462d0*XlX/625d0+0.22311413853594375d0*X*XlX/lX)
    
    
    
    end if
    if(orderZETA>=2) then
    
	!last check (30.09.2018) ?error?
!         val=val+(alpha**2)*(&
!         (-2.276736d0*lX)/(X-1d0)&
!         + (-0.07638842671474691d0 + (1.7423213067387915d0 + 0.6108031199759604d0*X)*X)/(X-1d0)**2 &
! 	+ (0.9872766517804755d0 + lX*(7.523963912264815d0 + (-6.946898670055185d0 - 1.9514483355910248d0*X)/X**2)&
! 	  - 6.946898670055185d0/X - 1.219024085938312d0*X + 0.41483265678863507d0*X/lX)*XlX**2)

	!(12.06.2019)
	val=val+(alpha**2)*(&
	(-2.276736*lX)/(-1 + X) + (-0.07638842671474769 + 1.742321306738789*X + 0.6108031199759592*X**2)/(-1 + X)**2 &
	+ ((6.946898670055183 + 1.9514483355910373*X - 7.523963912264807*X**2)*XlX)/((-1 + X)*X) &
	+ (0.9872766517804829 - 6.946898670055182/X - 1.2190240859383212*X + (0.4148326567886355*X)/lX)*XlX**2)
    end if
    
    if(orderZETA>=3) then
    !last check (30.09.2018) ?error?
!     val=val+(alpha**3)*(&
!     (-7.01234688d0*lX**2)/(-1d0 + X)**2 &
!      +(lX*(-18.257793687614754d0 + 27.607589527614756d0*X))/(-1d0 + X)**3 &
!      +(-138.8607743573032d0 + 184.52451857767127d0*X - 100.23304496672974d0*X**2 &
!      +32.07840713113214d0*X**3)/(-1d0 + X)**3 &
!      +(6.226625493462215d0 + lX*(59.335199839722506d0 - 42.01916318769151d0/X**2 - 27.831158601204876d0/X &
!      -22.05354500931568d0*X) - 25.83252187109758d0/X - 4.533025146478245d0*X &
!      +(0.7712919237997766d0*X)/lX + (lX**2*(-16.186641316593924d0 + 177.02234836509135d0*X &
!      -111.05224763708378d0*X**2 + 21.725060882031794d0*X**3))/X**3)*XlX**3)

      !(12.06.2019)
      val=val+(alpha**3)*(&
      (-7.01234688*lX**2)/(-1 + X)**2 - (0.7712919237997766*X**4)/(lX**4*(-1 + X)**3) &
      + (lX*(-18.257793687614754 + 27.607589527614756*X))/(-1 + X)**3 &
      + (4.533025146478245*X**2*(5.698737826585219 - 1.3736137109893911*X + X**2))/(lX**3*(-1 + X)**3) &
      + (-138.86077435730317 + 184.5245185776713*X - 100.23304496672974*X**2 + 32.07840713113214*X**3)/(-1 + X)**3 &
      + ((-16.186641316593928 + 198.62088100509132*X - 143.45004659708377*X**2 + 32.524327202031785*X**3)*XlX)/((-1 + X)**2*X) &
      + ((42.01916318769152 + 6.232625961204816*X - 37.736667199722454*X**2 + 14.854034129315712*X**3)*XlX**2)/((-1 + X)*X))
    
    end if
    
    
    else !!! --------------------------------------- nf =5 
    
    X=46d0/3d0*alpha*Log(mu*b*C0_inv_const)
    lX=Log(1d0-X)
    XlX=X/((1d0-X)*lX)
    
    !last check (30.09.2018)
    val=3d0*(X+lX)/lX/23d0
    
    if(orderZETA>=1) then
!         
      !last check (30.09.2018)
    val=val+alpha*(-2.401066092611533d0 + 174d0*lX/529d0 - 348d0*XlX/529d0 + 0.24322110206332886d0*X*XlX/lX)
    end if
    if(orderZETA>=2) then
	
! 	!last check (30.09.2018) ?error?
! 	val=val+(alpha**2)*(&
! 	(-1.658913454425906d0*lX)/(X-1d0)&
! 	+ (-1.4666874639336553d0 + (3.9789482118726536d0 - 0.8533472935130869d0*X)*X)/(X-1d0)**2&
! 	+ (-1.147671467270234d0 + lX*(12.41758546925291d0 + (-7.983935707092223d0-5.6603301030017965d0*X)/X**2)&
! 	  - 7.983935707092223d0/X - 0.7728481687522705d0*X + (0.4535332010815679d0*X)/lX)*XlX**2)
	  
	!(12.06.2019)
	val=val+(alpha**2)*(&
	(-1.658913454425906*lX)/(-1 + X) + (-1.4666874639336558 + 3.978948211872648*X - 0.8533472935130866*X**2)/(-1 + X)**2 &
	+ ((7.983935707092219 + 5.660330103001812*X - 12.417585469252895*X**2)*XlX)/((-1 + X)*X) &
	+ (-1.1476714672702268 - 7.983935707092221/X - 0.7728481687522828*X + (0.45353320108156825*X)/lX)*XlX**2)

    end if
    
     if(orderZETA>=3) then
    !last check (30.09.2018) ?error?
!     val=val+(alpha**3)*(&
!      (-4.183346972030546d0*lX**2)/(-1d0 + X)**2 &
!      +(lX*(-11.460061760095265d0 + 17.03785772280266d0*X))/(-1d0 + X)**3 &
!      +(-153.02822534041636d0 + 201.76779077125593d0*X - 95.24211018739625d0*X**2 &
!      +31.949115180427246d0*X**3)/(-1d0 + X)**3 &
!      +(-1.9927299189294243d0 + lX*(89.24835213177958d0 - 62.621068210206815d0/X**2 - 45.88349889123617d0/X &
!      -30.731897310115087d0*X) - 29.775211836053224d0/X - 2.8822524110838543d0*X &
!      +(0.8457011449184908d0*X)/lX + (lX**2*(-32.845856374153605d0 + 200.62114361461371d0*X &
!      - 92.74911411722383d0*X**2 + 12.506166030132125d0*X**3))/X**3)*XlX**3)
    
    	!(12.06.2019)
	val=val+(alpha**3)*(&
	(-4.183346972030546*lX**2)/(-1 + X)**2 - (0.8457011449184908*X**4)/(lX**4*(-1 + X)**3) &
	+ (lX*(-11.460061760095265 + 17.03785772280266*X))/(-1 + X)**3 &
	+ (2.8822524110838543*X**2*(10.330535841188325 + 0.6913793917792482*X + X**2))/(lX**3*(-1 + X)**3) &
	+ (-153.02822534041638 + 201.76779077125593*X - 95.24211018739621*X**2 + 31.949115180427235*X**3)/(-1 + X)**3 &
	+ ((-32.845856374153605 + 235.15441263635282*X - 144.5490176498325*X**2 + 29.772800541001683*X**3)*XlX)/((-1 + X)**2*X) &
	+ ((62.621068210206815 + 11.350229869497092*X - 54.715083110040446*X**2 + 19.22080763620205*X**3)*XlX**2)/((-1 + X)*X))
    
    end if
    end if
   
   !!!! in the gluon case we add difference
   if(f==0) then
    if(mu<mCHARM) then !nf=3
    
    !if(orderZETA>=1) val=val+alpha*() !!! +0
    if(orderZETA>=2) val=val+alpha**2*13.508464061376284d0*XlX
    if(orderZETA>=3) val=val+alpha**3*((-96.06018888089804*(lX**2 + &
	    lX*(0.007931353012306063 + 0.49603432349384596*X)*X - 0.2839618810968089*X**2)*XlX**2)/X**2)
    
    else if(mu<mBOTTOM) then  !nf=4
    if(orderZETA>=1) val=val+alpha/9d0
    if(orderZETA>=2) val=val+alpha**2*14.981629982484977d0*XlX
    if(orderZETA>=3) val=val+alpha**3*((-92.28684069210746*(lX**2 + &
	    lX*(0.007445441655734091 + 0.4962772791721318*X)*X - 0.3018318973700538*X**2)*XlX**2)/X**2)
    
    else  !nf=5
    if(orderZETA>=1) val=val+alpha*2d0/9d0
    if(orderZETA>=2) val=val+alpha**2*16.20788232334676*XlX
    if(orderZETA>=3) val=val+alpha**3*((-81.74410215253148*(lX**2 + &
	    lX*(-0.08375968910657726 + 0.5418798445532871*X)*X - 0.3697240315847729*X**2)*XlX**2)/X**2)
    end if
   
   end if

    
    if(ISNAN(val)) then
        write(*,*) 'ERROR: arTeMiDe.TMDR: zetaMuResum is NaN. At mu=',mu,'b=',b
        write(*,*) 'Lmu=',2*Log(mu*b*C0_inv_const), '1-X=',X
        write(*,*) 'Evaluation STOP'
        stop
    end if
    
    zetaMUresum=mu**2*Exp(-val/alpha)

 end function zetaMUresum

 !!!!!!!!!! exact value of zeta-line at given b,mu expanded over as
 function zetaSL(mu,b,f)
  real*8::zetaSL,b,mu
  integer::f
  real*8::dd,alpha,GD
  
  if(b<1d-6) b=1d-6
  
  alpha=As(mu)
  dd=DNP(mu,b,f)
  if(f==0) then
    GD=valueOfGD_type4_G(dd,alpha,mu)
  else
    GD=valueOfGD_type4_Q(dd,alpha,mu)
  end if
  
  !write(*,*) '..........',b,dd,alpha,GD
  
  zetaSL=mu**2*exp(-GD/dd)
 
 end function zetaSL
!-------------------- R-kernels   full----------------------------------

!-----------------------Improved D picture------------------------------
!!! Evolution exponent in the improved D-picture
 function TMDR_R_type1(b,muf,zetaf,mui,zetai,mu0,f)
  real*8::TMDR_R_type1,b,muf,zetaf,mui,zetai,mu0
  integer::f
  
  counter=0
  
  if(b<1d-6) b=1d-6
  
  TMDR_R_type1=EXP(IntegralG1(muf,mui,zetaf,f)-(IntegralG2(mui,mu0,f)+DNP(mu0,b,f))*Log(zetaf/zetai))
  
  !write(*,*) 'TMDR_R_type1: number of AD calls ',counter

 end function TMDR_R_type1
 
 !!!This is the integral   \int_{mu_i}^{muf} dmu/mu (Gamma(mu)Log[mu^2/zeta]-gammaV)
!!! evaluation by adaptive simpson
 function IntegralG1(muf,mui,zeta,f)
  real*8::IntegralG1,muf,mui,zeta,lnz
  integer::f
  
  real*8::y1,y3,y5,X1,X3,X5,maxV
  
  lnz=Log(zeta)
  
  y1=Log(muf)
  y5=Log(mui)
  if(y1<y5) then
   y3=y1
   y1=y5
   y5=y3
  end if
  y3=(y1+y5)/2d0
  
  X1=gammaCUSP(Exp(y1),f)*(2*y1-lnz)-gammaV(Exp(y1),f)
  X3=gammaCUSP(Exp(y3),f)*(2*y3-lnz)-gammaV(Exp(y3),f)
  X5=gammaCUSP(Exp(y5),f)*(2*y5-lnz)-gammaV(Exp(y5),f)
  
  counter=counter+3
  
  !maxV is approximate value of the integral, just to normalize the errors
  maxV=ABS((y1-y5)*(X1+4d0*X3+X5)/6d0)

  !!!! swithc sign if the integral is opposite
  if(muf>mui) then
  IntegralG1=integral1_S(y1,y5,lnz,f,X1,X3,X5,maxV)
  else
  IntegralG1=-integral1_S(y1,y5,lnz,f,X1,X3,X5,maxV)
  end if
  
 end function IntegralG1

!!! expected that y1>y2 (logarithm scale)
 recursive function integral1_S(y1,y5,lnz,f,X1,X3,X5,maxV) result(interX)
   real*8 ::y1,y5,lnz
   integer::f
   real*8 :: interX,X1,X2,X3,X4,X5
   real*8 :: value,valueAB,valueACB
   real*8 :: y2,y3,y4,deltay,maxV
   
   deltay=y1-y5
   y2=y5+deltay/4d0
   y3=y5+deltay/2d0
   y4=y1-deltay/4d0
   
   
   X2=gammaCUSP(Exp(y2),f)*(2*y2-lnz)-gammaV(Exp(y2),f)
   X4=gammaCUSP(Exp(y4),f)*(2*y4-lnz)-gammaV(Exp(y4),f)
   
   counter=counter+2
   
   valueAB=deltay*(X1+4d0*X3+X5)/6d0
   valueACB=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
!    write(*,*) y1,y3,y5,valueACB-valueAB
   
   If(ABS(valueACB-valueAB)/15>tolerance*maxV) then
    interX=integral1_S(y1,y3,lnz,f,X1,X2,X3,maxV)&
	  +integral1_S(y3,y5,lnz,f,X3,X4,X5,maxV)
   else
    interX=valueACB
   end if
   
 end function integral1_S


 !!!This is the integral   \int_{mu_i}^{muf} dmu/mu Gamma(mu)
!!! evaluation by adaptive simpson
 function IntegralG2(muf,mui,f)
  real*8::IntegralG2,muf,mui
  integer::f
  
  real*8::y1,y3,y5,X1,X3,X5,maxV
  
  if(Abs(muf-mui)<tolerance) then
    IntegralG2=0d0
    return
  end if
  
  y1=Log(muf)
  y5=Log(mui)
  if(y1<y5) then
   y3=y1
   y1=y5
   y5=y3
  end if
  y3=(y1+y5)/2d0
  
  X1=gammaCUSP(Exp(y1),f)
  X3=gammaCUSP(Exp(y3),f)
  X5=gammaCUSP(Exp(y5),f)
  
  counter=counter+3
  
  !maxV is approximate value of the integral, just to normalize the errors
  maxV=ABS((y1-y5)*(X1+4d0*X3+X5)/6d0)

  !!!! swithc sign if the integral is opposite
  if(muf>mui) then
  IntegralG2=integral2_S(y1,y5,f,X1,X3,X5,maxV)
  else
  IntegralG2=-integral2_S(y1,y5,f,X1,X3,X5,maxV)
  end if
  
 end function IntegralG2

!!! expected that y1>y2 (logarithm scale)
 recursive function integral2_S(y1,y5,f,X1,X3,X5,maxV) result(interX)
   real*8 ::y1,y5
   integer::f
   real*8 :: interX,X1,X2,X3,X4,X5
   real*8 :: value,valueAB,valueACB
   real*8 :: y2,y3,y4,deltay,maxV
   
   deltay=y1-y5
   y2=y5+deltay/4d0
   y3=y5+deltay/2d0
   y4=y1-deltay/4d0
   
   
   X2=gammaCUSP(Exp(y2),f)
   X4=gammaCUSP(Exp(y4),f)
   
   counter=counter+2
   
   valueAB=deltay*(X1+4d0*X3+X5)/6d0
   valueACB=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
!    write(*,*) y1,y3,y5,valueACB-valueAB
   
   If(ABS(valueACB-valueAB)/15>tolerance*maxV) then
    interX=integral2_S(y1,y3,f,X1,X2,X3,maxV)&
	  +integral2_S(y3,y5,f,X3,X4,X5,maxV)
   else
    interX=valueACB
   end if
   
 end function integral2_S

!-----------------------Improved gamma picture------------------------------
!!! Evolution exponent in the improved gamma-picture
 function TMDR_R_type2(b,muf,zetaf,mui,zetai,f)
  real*8::TMDR_R_type2,b,muf,zetaf,mui,zetai
  integer::f
  
  counter=0
  
  if(b<1d-6) b=1d-6
  
  TMDR_R_type2=EXP(-IntegralG3(muf,mui,b,f)+DNP(muf,b,f)*Log(muf**2/zetaf)-DNP(mui,b,f)*Log(mui**2/zetai))
  
  !write(*,*) 'TMDR_R_type2: number of AD calls ',counter
  
  if(TMDR_R_type2>1d6) then
    write(*,*) 'arTeMiDe.TMDR: ERROR -- Evolution factor is TOOO HUGE check the formula'
    write(*,*) 'NP parameters =', NPparam
    write(*,*) 'b=',b,'zetaf=',zetaf,'muf=',muf,'zetai=',zetai,'mui=',mui
    write(*,*) 'int=',IntegralG3(muf,mui,b,f),'t1=',DNP(muf,b,f)*Log(muf**2/zetaf),'t2=',DNP(mui,b,f)*Log(mui**2/zetai)
    write(*,*) 'Evaluation contirune with R=10^6'
    TMDR_R_type2=1d6
  end if

 end function TMDR_R_type2
 
 !!!This is the integral   \int_{mu_i}^{muf} dmu/mu (2 D(mu,b)+gammaV)
!!! evaluation by adaptive simpson
 function IntegralG3(muf,mui,b,f)
  real*8::IntegralG3,muf,mui,b
  integer::f
  
  real*8::y1,y3,y5,X1,X3,X5,maxV
  
  y1=Log(muf)
  y5=Log(mui)
  if(y1<y5) then
   y3=y1
   y1=y5
   y5=y3
  end if
  y3=(y1+y5)/2d0
  
  X1=2*DNP(Exp(y1),b,f)+gammaV(Exp(y1),f)
  X3=2*DNP(Exp(y3),b,f)+gammaV(Exp(y3),f)
  X5=2*DNP(Exp(y5),b,f)+gammaV(Exp(y5),f)
  
  counter=counter+3
  
  !maxV is approximate value of the integral, just to normalize the errors
  maxV=ABS((y1-y5)*(X1+4d0*X3+X5)/6d0)

  !!!! swithc sign if the integral is opposite
  if(muf>mui) then
  IntegralG3=integral3_S(y1,y5,b,f,X1,X3,X5,maxV)
  else
  IntegralG3=-integral3_S(y1,y5,b,f,X1,X3,X5,maxV)
  end if
  
 end function IntegralG3

!!! expected that y1>y2 (logarithm scale)
 recursive function integral3_S(y1,y5,b,f,X1,X3,X5,maxV) result(interX)
   real*8 ::y1,y5,b
   integer::f
   real*8 :: interX,X1,X2,X3,X4,X5
   real*8 :: value,valueAB,valueACB
   real*8 :: y2,y3,y4,deltay,maxV
   
   deltay=y1-y5
   y2=y5+deltay/4d0
   y3=y5+deltay/2d0
   y4=y1-deltay/4d0
   
   
   X2=2*DNP(Exp(y2),b,f)+gammaV(Exp(y2),f)
   X4=2*DNP(Exp(y4),b,f)+gammaV(Exp(y4),f)
   
   counter=counter+2
   
   valueAB=deltay*(X1+4d0*X3+X5)/6d0
   valueACB=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
!    write(*,*) y1,y3,y5,valueACB-valueAB
   
   If(ABS(valueACB-valueAB)/15>tolerance*maxV) then
    interX=integral3_S(y1,y3,b,f,X1,X2,X3,maxV)&
	  +integral3_S(y3,y5,b,f,X3,X4,X5,maxV)
   else
    interX=valueACB
   end if
   
 end function integral3_S
  
  
 !!! general interface for evolutions to zeta-line types 2,3,4
 function TMDR_R_typeZ(b,muf,zetaf,mui,zetai,f)
  real*8::TMDR_R_typeZ,b,muf,zetaf,mui,zetai
  integer::f
  
  SELECT CASE(EvolutionType)
    CASE(2)
      TMDR_R_typeZ=TMDR_R_type2(b,muf,zetaf,mui,zetai,f)
    CASE(3)
      TMDR_R_typeZ=TMDR_Rzeta_type3(b,muf,zetaf,f)/TMDR_Rzeta_type3(b,mui,zetai,f)
    CASE(4)
      TMDR_R_typeZ=TMDR_R_type4(b,muf,zetaf,mui,zetai,f)
    CASE DEFAULT
      TMDR_R_typeZ=TMDR_R_type4(b,muf,zetaf,mui,zetai,f)
    end SELECT
 
 end function TMDR_R_typeZ
!-------------------- R-kernels   to-zeta-line----------------------------------
 

!-----------------------Improved D picture------------------------------
!!! Evolution exponent in the improved D-picture to zeta-line
 function TMDR_Rzeta_type1(b,muf,zetaf,mui,mu0,f)
  real*8::TMDR_Rzeta_type1,b,muf,zetaf,mui,mu0,zetaP
  integer::f
  
  counter=0
  
  if(b<1d-6) b=1d-6
  
  zetaP=zetaNP(mui,b,f)
  
  TMDR_Rzeta_type1=EXP(IntegralG1(muf,mui,zetaf,f)-(IntegralG2(mui,mu0,f)+DNP(mu0,b,f))*Log(zetaf/zetaP))

  !write(*,*) 'TMDR_Rzeta_type1: number of AD calls ',counter
  
 end function TMDR_Rzeta_type1
 
!-----------------------Improved gamma picture------------------------------
!!! Evolution exponent in the improved gamma-picture to zeta-line
 function TMDR_Rzeta_type2(b,muf,zetaf,mui,f)
  real*8::TMDR_Rzeta_type2,b,muf,zetaf,mui,zetaP
  integer::f
  
  counter=0
  
  if(b<1d-6) b=1d-6
  
  zetaP=zetaNP(mui,b,f)
  
  TMDR_Rzeta_type2=EXP(-IntegralG3(muf,mui,b,f)+DNP(muf,b,f)*Log(muf**2/zetaf)-DNP(mui,b,f)*Log(mui**2/zetaP))
  
  !write(*,*) 'TMDR_Rzeta_type2: number of AD calls ',counter
  
  if(TMDR_Rzeta_type2>1d6) then
    write(*,*) 'arTeMiDe.TMDR: ERROR -- Evolution factor is TOOO HUGE check the formula'
    write(*,*) 'b=',b,'zetaf=',zetaf,'muf=',muf,'zetaP=',zetaP,'mui=',mui
    write(*,*) 'int=',IntegralG3(muf,mui,b,f),'t1=',DNP(muf,b,f)*Log(muf**2/zetaf),'t2=',DNP(mui,b,f)*Log(mui**2/zetaP)
    write(*,*) 'Evaluation continue with R=10^6'
    TMDR_Rzeta_type2=1d6
  end if

 end function TMDR_Rzeta_type2
 
 !!! general interface for evolutions to zeta-line for case of 3 give scales (zetaf,muf,mui)
 function TMDR_Rzeta_typeZ3(b,muf,zetaf,mui,f)
  real*8::TMDR_Rzeta_typeZ3,b,muf,zetaf,mui
  integer::f
 
  SELECT CASE(EvolutionType)
    CASE(2)
      TMDR_Rzeta_typeZ3=TMDR_Rzeta_type2(b,muf,zetaf,mui,f)
    CASE(4)
      TMDR_Rzeta_typeZ3=TMDR_Rzeta_type4(b,muf,zetaf,mui,f)
    CASE DEFAULT
      TMDR_Rzeta_typeZ3=TMDR_Rzeta_type4(b,muf,zetaf,mui,f)
  END SELECT
 end function TMDR_Rzeta_typeZ3
 
!!! Evolution exponent in the improved gamma-picture to zeta-line (defined by zetaNP)
 function TMDR_Rzeta_type3(b,muf,zetaf,f)
  real*8::TMDR_Rzeta_type3,b,muf,zetaf,zetaP
  integer::f
  
  if(b<1d-6) b=1d-6
  
  zetaP=zetaNP(muf,b,f)
  
  TMDR_Rzeta_type3=EXP(-DNP(muf,b,f)*Log(zetaf/zetaP))
  
  !write(*,*) 'HERE'
  
  if(TMDR_Rzeta_type3>1d6) then
    write(*,*) 'arTeMiDe.TMDR: ERROR -- Evolution factor is TOOO HUGE check the formula'
    write(*,*) 'b=',b,'zetaf=',zetaf,'muf=',muf,'zetaP=',zetaP    
    write(*,*) 'DNP=',DNP(muf,b,f), 'log(zeta/zetamu)=',Log(zetaf/zetaP)
    write(*,*) 'NPparameters= (',NPparam,')'
    write(*,*) 'Evaluation continue with R=10^6'
    TMDR_Rzeta_type3=1d6
    stop
  end if

 end function TMDR_Rzeta_type3

 !!! Evolution exponent in the exact solution to zeta-line (defined by zetaNP)
 function TMDR_Rzeta_type4(b,muf,zetaf,mui,f)
  real*8::TMDR_Rzeta_type4,b,muf,zetaf,mui,zetai
  integer::f
  real*8::dd1,alpha1,GD1 !!! at point f
  real*8::dd2,alpha2,GD2 !!! at point i
  
  if(b<1d-6) b=1d-6
    
  zetai=zetaNP(mui,b,f)
  
  alpha1=As(muf)
  dd1=DNP(muf,b,f)  
  alpha2=As(mui)
  dd2=DNP(mui,b,f)  
  
  if(f==0) then
    GD1=valueOfGD_type4_G(dd1,alpha1,muf)
    GD2=valueOfGD_type4_G(dd2,alpha2,mui)  
  else
    GD1=valueOfGD_type4_Q(dd1,alpha1,muf)
    GD2=valueOfGD_type4_Q(dd2,alpha2,mui)
  end if
  
  !write(*,*) '----> ',dd1,dd2,GD1,GD2
  if(b>10d0 .and. GD1<GD2) GD2=GD1
  
  
  TMDR_Rzeta_type4=EXP(-dd1*Log(zetaf/muf**2)-GD1+dd2*Log(zetai/mui**2)+GD2)
  
  !write(*,*) b,muf,zetaf,mui,zetai,TMDR_Rzeta_type4
  
  if(TMDR_Rzeta_type4>1d6) then
    write(*,*) 'arTeMiDe.TMDR: ERROR -- Evolution factor T4 (TMD_Rzeta) is TOOO HUGE check the formula'
    write(*,*) 'b=',b,'ln zetaf=',log(zetaf) ,'muf=',muf,'ln zetaSL(f)=',log(muf**2)-GD1/dd1
    write(*,*) 'DNP(f)=',dd1, 'g*dd(f)=',GD1
    write(*,*) 'b=',b,'ln zetai=',log(zetai),'mui=',mui,'ln zetaSL(i)=',log(mui**2)-GD2/dd2
    write(*,*) 'DNP(i)=',dd2, 'g*dd(i)=',GD2
    write(*,*) 'NPparameters= (',NPparam,')'
    write(*,*) 'Evaluation continue with R=10^6'
    TMDR_Rzeta_type4=1d6
    stop
  end if

 end function TMDR_Rzeta_type4
 
 !!! Evolution exponent in the exact solution from point to point
 function TMDR_R_type4(b,muf,zetaf,mui,zetai,f)
  real*8::TMDR_R_type4,b,muf,zetaf,mui,zetai
  integer::f
  real*8::dd1,alpha1,GD1 !!! at point f
  real*8::dd2,alpha2,GD2 !!! at point i
  
  if(b<1d-6) b=1d-6
  
  alpha1=As(muf)
  dd1=DNP(muf,b,f)  
  
  alpha2=As(mui)
  dd2=DNP(mui,b,f)
  
  if(f==0) then
    GD1=valueOfGD_type4_G(dd1,alpha1,muf)
    GD2=valueOfGD_type4_G(dd2,alpha2,mui)
  else
    GD1=valueOfGD_type4_Q(dd1,alpha1,muf)
    GD2=valueOfGD_type4_Q(dd2,alpha2,mui)
  end if
  
  
  TMDR_R_type4=EXP(-dd1*Log(zetaf/muf**2)-GD1+dd2*Log(zetai/mui**2)+GD2)
  
  !write(*,*) b,muf,zetaf,mui,zetai,TMDR_R_type4
  
  if(TMDR_R_type4>1d6) then
    write(*,*) 'arTeMiDe.TMDR: ERROR -- Evolution factor T4 (TMD_R) is TOOO HUGE check the formula'
    write(*,*) 'b=',b,'zetaf=',zetaf,'muf=',muf
    write(*,*) 'DNP(f)=',dd1, 'g*dd(f)=',GD1
    write(*,*) 'b=',b,'zetai=',zetai,'mui=',mui
    write(*,*) 'DNP(i)=',dd2, 'g*dd(i)=',GD2
    write(*,*) 'NPparameters= (',NPparam,')'
    write(*,*) 'Evaluation continue with R=10^6'
    TMDR_R_type4=1d6
    stop
  end if

 end function TMDR_R_type4
 
 
 ! expression for G*D, at given D,alpha, and mu
 ! evaluated for type4 evolution for Quark
 function valueOfGD_type4_Q(dd,alpha,mu)
  real*8::valueOfGD_type4_Q,dd,alpha,mu
  real*8::val,p,ee
  
  !!p=2 beta0 d /Gamma0
  !!val=dd*g
  
    !---------Nf=3-------------
  if(mu<mCHARM) then
    p=27d0*dd/8d0
    ee=Exp(-p)
    val=(p-1d0+ee)*8d0/243d0
    
    if(orderZETA>=1) then
      val=val+alpha*(0.06647850165595162*(1d0-ee) - 0.510922946100396d0*p - 0.11705532693187014d0*p**2)
    end if
    if(orderZETA>=2) then
      val=val+alpha**2*(-1.3616080333159475d0+ 1.7886575769140867d0/ee - 0.4270495435981395d0*ee + 0.21760590764719615d0*p)
    end if
    if(orderZETA>=3) then
      val=val+alpha**3*(2.7269428180833457d0 - 9.811304303151902d0/ee**2 + 3.6118084990477106d0/ee &
      + 3.472552986020849d0*ee + 7.770743143365355d0*p)
    end if
    !---------Nf=4-------------
  else if(mu<mBOTTOM) then
    p=25d0*dd/8d0
    ee=Exp(-p)
    val=(p-1d0+ee)*24d0/625d0
    
    if(orderZETA>=1) then
      val=val+alpha*(0.071396524331502d0*(1-ee) - 0.551396524331502d0*p - 0.118272d0*p**2)
    end if
    if(orderZETA>=2) then
      val=val+alpha**2*(-2.3832241553760207d0 + 2.706938637790099d0/ee - 0.32371448241407974d0*ee - 0.02444429654871905d0*p)
    end if
    if(orderZETA>=3) then
      val=val+alpha**3*(1.344164561788757d0 - 11.814572956005321d0/ee**2 + 5.032969018668324d0/ee &
      + 5.4374393755482275d0*ee + 9.06362014286148d0*p)
    end if
   !---------Nf=5-------------
   else
    p=23d0*dd/8d0
    ee=Exp(-p)
    val=(p-1d0+ee)*24d0/529d0
    
    if(orderZETA>=1) then
      val=val+alpha*(0.08459864419594047d0*(1d0-ee) - 0.6063377746307231d0*p - 0.11440782444316594d0*p**2)
    end if
    if(orderZETA>=2) then
      val=val+alpha**2*(-3.8090080018501724d0 + 3.9989494168802584d0/ee - 0.1899414150300844d0*ee - 0.5101521613682279d0*p)
    end if
    if(orderZETA>=3) then
      val=val+alpha**3*(-0.07909327666015716d0 - 14.515461613202893d0/ee**2 + 7.456821446063274d0/ee &
      + 7.137733443799741d0*ee + 10.43484998657383d0*p)
    end if
   end if
   
   valueOfGD_type4_Q=val/alpha
 end function valueOfGD_type4_Q

 ! expression for G*D, at given D,alpha, and mu
 ! evaluated for type4 evolution for GLUON
 function valueOfGD_type4_G(dd,alpha,mu)
  real*8::valueOfGD_type4_G,dd,alpha,mu
  real*8::val,p,ee
  
  !!p=2 beta0 d /Gamma0
  !!val=dd*g
  
    !---------Nf=3-------------
  if(mu<mCHARM) then
    p=3d0*dd/2d0
    ee=Exp(-p)
    val=(p-1d0+ee)*2d0/27d0
    
    if(orderZETA>=1) then
      val=val+alpha*(0.14957662872589117d0*(1- ee) - 1.1495766287258913d0*p - 0.26337448559670784d0*p**2)
    end if
    if(orderZETA>=2) then
      val=val+alpha**2*(5.94202463262331d0 - 4.9811631595274966d0/ee - 0.9608614730958139d0*ee + 0.48961329220619126d0*p)
    end if
    if(orderZETA>=3) then
      val=val+alpha**3*(37.9017218781724d0 - 35.65658059748694d0/ee**2 - 10.058385499232372d0/ee &
      + 7.81324421854691d0*ee + 17.484172072572047d0*p)
    end if
    !---------Nf=4-------------
  else if(mu<mBOTTOM) then
    p=25d0*dd/18d0
    ee=Exp(-p)
    val=(p-1d0+ee)*54d0/625d0
    
    if(orderZETA>=1) then
      val=val+alpha*(0.1606421797458795d0*(1-ee) - 1.1606421797458795d0*p - 0.266112d0*p**2)
    end if
    if(orderZETA>=2) then
      val=val+alpha**2*(5.424519237793142d0 - 4.696161652361463d0/ee - 0.7283575854316795d0*ee - 0.05499966723461787d0*p)
    end if
    if(orderZETA>=3) then
      val=val+alpha**3*(36.000271049515874d0 - 39.50300913206475d0/ee**2 - 8.731500512434682d0/ee &
      + 12.234238594983513d0*ee + 20.393145321438325d0*p)
    end if
   !---------Nf=5-------------
   else
    p=23d0*dd/18d0
    ee=Exp(-p)
    val=(p-1d0+ee)*54d0/529d0
    
    if(orderZETA>=1) then
      val=val+alpha*(0.19034694944086608d0*(1-ee) - 1.190346949440866d0*p - 0.25741760499712335d0*p**2)
    end if
    if(orderZETA>=2) then
      val=val+alpha**2*(4.114161640195447d0 - 3.686793456377757d0/ee - 0.42736818381768993d0*ee - 1.1478423630785126d0*p)
    end if
    if(orderZETA>=3) then
      val=val+alpha**3*(34.48806902386751d0 - 43.67322352127637d0/ee**2 - 6.874745751140516d0/ee &
      + 16.059900248549415d0*ee + 23.47841246979112d0*p)
    end if
   end if
   
   valueOfGD_type4_G=val/alpha
 end function valueOfGD_type4_G

 
 
 !--------------------------------------------------------------------------------------------------------
 !------------------------EVOLUTION TO SPECIAL LINE-------------------------------------------------------
  !!! Evolution exponent in the exact solution to exact zeta-line
 function TMDR_R_toSL(b,muf,zetaf,f)
  real*8::TMDR_R_toSL,b,muf,zetaf
  integer::f
  real*8::dd,alpha,GD
  
  if(b<1d-6) b=1d-6
  
  alpha=As(muf)
  dd=DNP(muf,b,f)
  if(f==0) then
    GD=valueOfGD_type4_G(dd,alpha,muf)
  else
    GD=valueOfGD_type4_Q(dd,alpha,muf)
  end if
  
  TMDR_R_toSL=EXP(-dd*Log(zetaf/muf**2)-GD)
  
  
  if(TMDR_R_toSL>1d6) then
    write(*,*) 'arTeMiDe.TMDR: ERROR -- Evolution factor N4 is TOOO HUGE check the formula'
    write(*,*) 'b=',b,'zetaf=',zetaf,'muf=',muf
    write(*,*) 'DNP=',dd, 'g*dd=',GD
    write(*,*) 'NPparameters= (',NPparam,')'
    write(*,*) 'Evaluation continue with R=10^6'
    TMDR_R_toSL=1d6
    stop
  end if

 end function TMDR_R_toSL
 
 
 !! we search for the lowest avalible Q, under which the evolution is inverted (for quark)
  !! it is defined by Q^2<zeta_Q(b->infty)
  !! with variabtions it is Q^2<zeta_{c Q}(b->infty)
  !! the function return the values for c={0.5,1,2} at b=25 (let it be assimptotic value)
  !! The function returns 1 if Q<1
 function LowestQ()
  real*8,dimension(1:3)::LowestQ
  
  real*8,parameter::b=25d0
  integer,parameter::f=1
  real*8::Qs1,Qs2,Qs3,V1,V2,V3
  
  !c=1
  Qs1=1d0
  V1=Qs1**2-zetaNP(Qs1,b,f)
  
  Qs3=25d0
  V3=Qs3**2-zetaNP(Qs3,b,f)
  
  if(V1<0d0 .and. V3>0d0) then
    do
      Qs2=(Qs3+Qs1)/2d0
      V2=Qs2**2-zetaNP(Qs2,b,f)
      if(V2>0d0) then 
	Qs3=Qs2
	V3=V2
      else
	Qs1=Qs2
	V1=V2
      end if
      if(Qs3-Qs1<tolerance) exit
    end do
    LowestQ(2)=(Qs3+Qs1)/2d0
  else if(V1>0d0 .and. V3>0d0) then
    LowestQ(2)=1d0
  else
    write(*,*) 'arTeMiDe.TMDR: LowestQ ERROR. Negative value at large Q (c=1).'
    stop
  end if
  
    !c=0.5
  Qs1=1d0
  V1=Qs1**2-zetaNP(0.5d0*Qs1,b,f)
  
  Qs3=25d0
  V3=Qs3**2-zetaNP(0.5d0*Qs3,b,f)
  
  if(V1<0d0 .and. V3>0d0) then
    do
      Qs2=(Qs3+Qs1)/2d0
      V2=Qs2**2-zetaNP(0.5d0*Qs2,b,f)
      if(V2>0d0) then 
	Qs3=Qs2
	V3=V2
      else
	Qs1=Qs2
	V1=V2
      end if
      if(Qs3-Qs1<tolerance) exit
    end do
    LowestQ(1)=(Qs3+Qs1)/2d0
  else if(V1>0d0 .and. V3>0d0) then
    LowestQ(1)=1d0
  else
    write(*,*) 'arTeMiDe.TMDR: LowestQ ERROR. Negative value at large Q (c=1).'
    stop
  end if
  
    !c=2
  Qs1=1d0
  V1=Qs1**2-zetaNP(2d0*Qs1,b,f)
  
  Qs3=25d0
  V3=Qs3**2-zetaNP(2d0*Qs3,b,f)
  
  if(V1<0d0 .and. V3>0d0) then
    do
      Qs2=(Qs3+Qs1)/2d0
      V2=Qs2**2-zetaNP(2d0*Qs2,b,f)
      if(V2>0d0) then 
	Qs3=Qs2
	V3=V2
      else
	Qs1=Qs2
	V1=V2
      end if
      if(Qs3-Qs1<tolerance) exit
    end do
    LowestQ(3)=(Qs3+Qs1)/2d0
  else if(V1>0d0 .and. V3>0d0) then
    LowestQ(3)=1d0
  else
    write(*,*) 'arTeMiDe.TMDR: LowestQ ERROR. Negative value at large Q (c=1).'
    stop
  end if
 
 end function LowestQ


end module TMDR