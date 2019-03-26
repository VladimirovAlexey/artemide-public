!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 1.41
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
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDR
use QCDinput

implicit none

  private
!   public
 
 !Current version of module
 character (len=5),parameter :: version="v1.41"
 
!------------------------------------------Physical and mathematical constants------------------------------------------

!!!Threashold parameters
! Set in the constants file
  real*8 :: mCHARM=1.4d0
  real*8 :: mBOTTOM=4.75d0
  
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
  
!! Precision tolerance used in various routines
  real*8::tolerance=0.001d0

!! Evolutio type 
  integer:: EvolutionType=3
  
  integer::counter=0
  
!------------------------------------------Non-pertrubative parameters--------------------------------------------------
  
  !!Number of non-pertrubative parameters
  integer::NPlength=0
  !! array of non-pertrubative parameters
  real*8,allocatable,dimension(:):: NPparam

  public::TMDR_R,TMDR_Rzeta
  public:: TMDR_Initialize,TMDR_setNPparameters,LowestQ
  
    public::DNP!,GammaCusp,gammaV
  interface TMDR_setNPparameters
   module procedure TMDR_setNPparameters, TMDR_SetReplica
  end interface 
  
  interface TMDR_R
    module procedure TMDR_R_type1,TMDR_R_type2
  end interface
  
  interface TMDR_Rzeta
    module procedure TMDR_Rzeta_type1,TMDR_Rzeta_type2,TMDR_Rzeta_type3
  end interface
  
 contains
  
   INCLUDE 'Model/TMDR_model.f90'
  
  !!! Initializing routing
!!! Filles the prebuiled arrays
!!! orderAD, is order of anomalous dimension for evolution
!!!! order zeta is the order of pertrubative zeta expression, used only in the "universal TMD"
subroutine TMDR_Initialize(orderMain)
    character(len=*)::orderMain
    character (len=50) :: line
    integer ::i
    
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
    
    if(outputLevel>1) write(*,*) '----- arTeMiDe.TMDR ',version,': .... initialization'
    
      SELECT CASE(orderMain)
      CASE ("LO")
	if(outputLevel>1) write(*,*) 'Order set: LO'
	orderCusp=1
	orderV=0
	orderD=0
	orderDresum=0
	orderZETA=0
      CASE ("LO+")
	if(outputLevel>1) write(*,*) 'Order set: LO+'
	orderCusp=1
	orderV=1
	orderD=1
	orderDresum=0
	orderZETA=0
      CASE ("NLO")
	if(outputLevel>1) write(*,*) 'Order set: NLO'
	orderCusp=2
	orderV=1
	orderD=1
	orderDresum=1
	orderZETA=1
      CASE ("NLO+")
	if(outputLevel>1) write(*,*) 'Order set: NLO+'
	orderCusp=2
	orderV=2
	orderD=2
	orderDresum=1
	orderZETA=1
      CASE ("NNLO")
	if(outputLevel>1) write(*,*) 'Order set: NNLO'
	orderCusp=3
	orderV=2
	orderD=2
	orderDresum=2
	orderZETA=2
      CASE ("NNLO+")
	if(outputLevel>1) write(*,*) 'Order set: NNLO+'
	orderCusp=4!3
	orderV=3
	orderD=3
	orderDresum=3!2
	orderZETA=3!2
      CASE ("NNNLO")
	if(outputLevel>1) write(*,*) 'Order set: NNNLO'
	orderCusp=4
	orderV=3
	orderD=3
	orderDresum=3
	orderZETA=3
      CASE DEFAULT
	if(outputLevel>0) write(*,*) 'arTeMiDe.TMDR_Initialize:try to set unknown ADs orders. Switch to NLO.'
	if(outputLevel>1) write(*,*) 'Order set: NLO'
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
    
    if(QCDinput_IsInitialized()) then
      if(outputLevel>1) write(*,*) 'QCDinput is already initalized'
    else
      if(outputLevel>1) write(*,*) '.. initializing QCDinput'
      call QCDinput_Initialize(orderMain)
      if(outputLevel>1) write(*,*) 'QCDinput is initalized'
    end if

     !!!!!!! --------------Read threashod parameters-------------------------
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
!   
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
    !!!!search for TMDR entry
    do
    read(51,'(A)') line
    if(line(1:3)=='*C ') exit
    end do
    read(51,'(A)') line
    read(51,*) tolerance
    
    if(outputLevel>2) write(*,'(A,ES10.3)') ' |   tolerance=',tolerance
! 
!     
    !!!!!!! --------------Length of NP input-------------------------
    
    do
    read(51,'(A)') line    
    if(line(1:3)=='*4 ') exit
    end do  
    do
    read(51,'(A)') line
    if(line(1:3)=='*A ') exit
    end do
    do
    read(51,'(A)') line
    if(line(1:3)=='*0)') exit
    end do
    read(51,*) NPlength
    
    if(outputLevel>2) write(*,'(A,I3)') ' Number of NP parameters =',NPlength
    
    !!allocating number of NP input
    allocate(NPparam(1:NPlength))
    do i=1,NPlength
    NPparam(i)=0d0
    end do
!   
    !!!!!! --------------Evolution type------------------------------
    
    do
    read(51,'(A)') line
    if(line(1:3)=='*B ') exit
    end do
    do
    read(51,'(A)') line
    if(line(1:3)=='*1)') exit
    end do
    read(51,*) EvolutionType
    
    if(outputLevel>2) write(*,'(A,I3)') ' Evolution type =',EvolutionType
    
    CLOSE (51, STATUS='KEEP') 
    
    if(outputLevel>2) write(*,*) 'Model initialization..'
    call ModelInitialization()
!           
     if(outputLevel>0) write(*,*) '----- arTeMiDe.TMDR ',version,': .... initialized'
     if(outputLevel>1) write(*,*) ' '
     
  end subroutine TMDR_Initialize

    !Reset the NP parameters
 subroutine TMDR_setNPparameters(NPin)
    real*8,dimension(1:NPlength)::NPin
    
    NPparam=NPin
    
    if(outputLevel>2) write(*,*) 'arTeMiDe.TMDR: NPparameters reset = (',NPparam,')'
    
 end subroutine TMDR_setNPparameters



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
! so far only quarks
 function Dresum(mu,b,f)
    real*8::Dresum,mu,b
    integer::f
    real*8:: X,alpha,lX
    
    if(f==0) then
    write(*,*) 'ERROR: arTeMiDe.TMDR. Resummed version of D is made only for quark'
    write(*,*) 'Evaluation stop'
    stop
    end if
    
    alpha=As(mu)
    If(mu<mCHARM) then !! nf=3
    !X=1-b0*Lmu
    X=18d0*alpha*Log(mu*b*C0_inv_const)
    lX=Log(1d0-X)
    
    !last check (30.09.2018)
    Dresum=-8d0/27d0*lX
    
    if(orderDresum>=1) then
    !last check (30.09.2018)
    Dresum=Dresum+alpha/(1d0-X)*(0.5983065149035645d0*X - 512d0*lX/243d0)    
      if(orderDresum>=2) then
	!last check (30.09.2018)
        Dresum=Dresum+((alpha/(1d0-X))**2)*(&
        -15.75963102138173d0 - 19.23770595326028d0*lX +16384d0*lX**2/2187d0&
        + 15.108120124066376d0*X - 4.447521424630353d0*X**2)
        !!!!!4-loop
        if(orderDresum>=3) then
	  !last check (30.09.2018)
	  Dresum=Dresum+((alpha/(1d0-X))**3)*(&
	  -260.5342569653865d0 - 20.10001091244974d0*lX + 190.07464445795534d0*lX**2 &
	  -35.51545326762519d0*lX**3 + 146.58790833901065d0*X - 44.18188284306254d0*lX*X &
	  +11.703452376448753d0*X**2 - 37.75565322877677d0*X**3)
	end if
      end if
    end if
    
    else if(mu<mBOTTOM) then !! nf=4
    X=50d0/3d0*alpha*Log(mu*b*C0_inv_const)
    lX=Log(1d0-X)
    
    !last check (30.09.2018)
    Dresum=-8d0/25d0*lX
    if(orderDresum>=1) then
        !last check (30.09.2018)
        Dresum=Dresum+alpha/(1d0-X)*(0.5949710360958496d0*X - 1232d0*lX/625d0)
    if(orderDresum>=2) then
	!last check (30.09.2018)
        Dresum=Dresum+((alpha/(1d0-X))**2)*(&
        -18.52506312014714d0 - 15.807613582350436d0*lX + 94864d0*lX**2/15625d0&
        + 9.962780902782157d0*X - 3.2507308958355225d0*X**2)
        !!!!!4-loop
        if(orderDresum>=3) then
	  !last check (30.09.2018)
	  Dresum=Dresum+((alpha/(1d0-X))**3)*(&
	  -305.836516979261d0 + 69.4831476117966d0*lX + 134.7740830272787d0*lX**2 &
	  -24.93278890666667d0*lX**3 + 42.43828293158347d0*X - 21.321725724444445d0*lX*X &
	  +80.87368989899981d0*X**2 - 58.80945335817513d0*X**3)
	end if
    end if
    end if
    else !!! nf =5 
    X=46d0/3d0*alpha*Log(mu*b*C0_inv_const)
    lX=Log(1d0-X)
    
    Dresum=-8d0/23d0*lX
    if(orderDresum>=1) then
	!last check (30.09.2018)
        Dresum=Dresum+alpha/(1d0-X)*(0.6485896055022099d0*X - 928d0*lX/529d0)
     if(orderDresum>=2) then
	!last check (30.09.2018)
        Dresum=Dresum+((alpha/(1d0-X))**2)*(&
        -21.2904952189126d0 - 12.118685999181197d0*lX + 53824d0*lX**2/12167d0&
        + 3.4818379050987946d0*X - 2.0609284500060885d0*X**2)
        !!!!4-loop
        if(orderDresum>=3) then
	  !last check (30.09.2018)
	  Dresum=Dresum+((alpha/(1d0-X))**3)*(&
	  -342.0455323786137d0 + 136.0753964300541d0*lX + 83.43151323800433d0*lX**2 &
	  -14.874122567219718d0*lX**3 - 69.48701186506418d0*X + 3.227921887389228d0*lX*X &
	  +155.46041408237463d0*X**2 - 81.95172616030685d0*X**3)
	end if
        
    end if
    end if
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

    
 function zetaMUresum(mu,b,f)
  real*8::zetaMUresum,mu,b
  integer::f
  real*8::X,val,lX,alpha,XlX
  
  if(f==0) then
    write(*,*) 'ERROR: arTeMiDe.TMDR. Resummed version of zeta-line is made only for quark'
    write(*,*) 'Evaluation stop'
    stop
  end if
  
  !val=g_0+a g_1+...
  
   alpha=As(mu)
    
    If(mu<mCHARM) then !! nf=3
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
    !last check (30.09.2018)
    val=val+(alpha**2)*(&
    (-2.809327846364883d0*lX)/(X-1d0)&
    +(0.7344199383092869d0 + X*(0.5427459498748931d0 +1.5321619581807067d0*X))/(X-1d0)**2&
    +(2.4745769670392037d0 + lX*(3.8610071741320366d0 + (-5.9098616330181475d0+0.4533704191432824d0*X)/X**2)&
      - 5.9098616330181475d0/X - 1.6678205342363737d0*X + 0.4530566491853122d0*X/lX)*XlX**2)
    
    end if
    if(orderZETA>=3) then
    !last check (30.09.2018)
    val=val+(alpha**3)*(&
    (-9.988721231519586d0*lX**2)/(-1d0 + X)**2 &
    +(lX*(-25.650274604347043d0 + 38.96856957970649d0*X))/(-1d0 + X)**3 &
    +(-116.72083321043493d0 + 157.99440284791794d0*X - 104.75296305816975*X**2 &
    +32.15628667503184d0*X**3)/(-1d0 + X)**3 &
    +(13.21548593880104d0 + lX*(35.90538079987103d0 - 25.582672835245578d0/X**2 - 16.566763069323986d0/X &
    -14.158369960791271d0*X) - 23.8673838411903d0/X - 6.7356082664341095d0*X &
    +(0.9148502637943364d0*X)/lX + (lX**2*(-1.7152889940552896d0 + 138.17793148235785d0*X &
    -113.52541590924972d0*X**2 + 26.948941036555656d0*X**3))/X**3)*XlX**3)
    
    end if
    
    else if(mu<mBOTTOM) then !! nf=4
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
    
	!last check (30.09.2018)
        val=val+(alpha**2)*(&
        (-2.276736d0*lX)/(X-1d0)&
        + (-0.07638842671474691d0 + (1.7423213067387915d0 + 0.6108031199759604d0*X)*X)/(X-1d0)**2 &
	+ (0.9872766517804755d0 + lX*(7.523963912264815d0 + (-6.946898670055185d0 - 1.9514483355910248d0*X)/X**2)&
	  - 6.946898670055185d0/X - 1.219024085938312d0*X + 0.41483265678863507d0*X/lX)*XlX**2)
    end if
    
    if(orderZETA>=3) then
    !last check (30.09.2018)
    val=val+(alpha**3)*(&
    (-7.01234688d0*lX**2)/(-1d0 + X)**2 &
     +(lX*(-18.257793687614754d0 + 27.607589527614756d0*X))/(-1d0 + X)**3 &
     +(-138.8607743573032d0 + 184.52451857767127d0*X - 100.23304496672974d0*X**2 &
     +32.07840713113214d0*X**3)/(-1d0 + X)**3 &
     +(6.226625493462215d0 + lX*(59.335199839722506d0 - 42.01916318769151d0/X**2 - 27.831158601204876d0/X &
     -22.05354500931568d0*X) - 25.83252187109758d0/X - 4.533025146478245d0*X &
     +(0.7712919237997766d0*X)/lX + (lX**2*(-16.186641316593924d0 + 177.02234836509135d0*X &
     -111.05224763708378d0*X**2 + 21.725060882031794d0*X**3))/X**3)*XlX**3)
    
    end if
    
    
    else !!! nf =5 
    
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
	
	!last check (30.09.2018)
	val=val+(alpha**2)*(&
	(-1.658913454425906d0*lX)/(X-1d0)&
	+ (-1.4666874639336553d0 + (3.9789482118726536d0 - 0.8533472935130869d0*X)*X)/(X-1d0)**2&
	+ (-1.147671467270234d0 + lX*(12.41758546925291d0 + (-7.983935707092223d0-5.6603301030017965d0*X)/X**2)&
	  - 7.983935707092223d0/X - 0.7728481687522705d0*X + (0.4535332010815679d0*X)/lX)*XlX**2)

    end if
    
     if(orderZETA>=3) then
    !last check (30.09.2018)
    val=val+(alpha**3)*(&
     (-4.183346972030546d0*lX**2)/(-1d0 + X)**2 &
     +(lX*(-11.460061760095265d0 + 17.03785772280266d0*X))/(-1d0 + X)**3 &
     +(-153.02822534041636d0 + 201.76779077125593d0*X - 95.24211018739625d0*X**2 &
     +31.949115180427246d0*X**3)/(-1d0 + X)**3 &
     +(-1.9927299189294243d0 + lX*(89.24835213177958d0 - 62.621068210206815d0/X**2 - 45.88349889123617d0/X &
     -30.731897310115087d0*X) - 29.775211836053224d0/X - 2.8822524110838543d0*X &
     +(0.8457011449184908d0*X)/lX + (lX**2*(-32.845856374153605d0 + 200.62114361461371d0*X &
     - 92.74911411722383d0*X**2 + 12.506166030132125d0*X**3))/X**3)*XlX**3)
    
    end if
    
    end if
    
    if(ISNAN(val)) then
        write(*,*) 'ERROR: arTeMiDe.TMDR: zetaMuResum is NaN. At mu=',mu,'b=',b
        write(*,*) 'Lmu=',2*Log(mu*b*C0_inv_const), '1-X=',X
        write(*,*) 'Evaluation STOP'
        stop
    end if
    
    zetaMUresum=mu**2*Exp(-val/alpha)
!     write(*,*) 'alpha=',alpha

 end function zetaMUresum

 
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
    write(*,*) 'Evaluation contirune with R=10^6'
    TMDR_Rzeta_type2=1d6
  end if

 end function TMDR_Rzeta_type2
 
!!! Evolution exponent in the improved gamma-picture to zeta-line
 function TMDR_Rzeta_type3(b,muf,zetaf,f)
  real*8::TMDR_Rzeta_type3,b,muf,zetaf,zetaP
  integer::f
  
  if(b<1d-6) b=1d-6
  
  zetaP=zetaNP(muf,b,f)
  
  TMDR_Rzeta_type3=EXP(-DNP(muf,b,f)*Log(zetaf/zetaP))
  
  !write(*,*) 'TMDR_Rzeta_type2: number of AD calls ',counter
  
  if(TMDR_Rzeta_type3>1d6) then
    write(*,*) 'arTeMiDe.TMDR: ERROR -- Evolution factor is TOOO HUGE check the formula'
    write(*,*) 'b=',b,'zetaf=',zetaf,'muf=',muf,'zetaP=',zetaP
    write(*,*) 'b=',b,'zetaf=',zetaf,'muf=',muf,'zetaP=',zetaP
    write(*,*) 'DNP=',DNP(muf,b,f), 'log(zeta/zetamu)=',Log(zetaf/zetaP)
    write(*,*) 'Evaluation contirune with R=10^6'
    TMDR_Rzeta_type3=1d6
    stop
  end if

 end function TMDR_Rzeta_type3

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