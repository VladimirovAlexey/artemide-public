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
use aTMDe_Numerics
use IO_functions
use QCDinput
use TMD_AD

implicit none

  private
!   public
 
 !Current version of module
 character (len=5),parameter :: version="v2.03"
 character (len=7),parameter :: moduleName="TMDR"
 !Last appropriate verion of constants-file
  integer,parameter::inputver=10

  character(256)::replicaFILE
  logical::usereplicaFILE=.false.
  character(50)::name
  
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
  real(dp)::tolerance=0.0001d0

!! Evolution type 
  integer:: EvolutionType=4
  
  integer::counter,messageCounter
  
  
!------------------------------------------Non-pertrubative parameters--------------------------------------------------
  
  !!Number of non-pertrubative parameters
  integer::NPlength=0
  !! array of non-pertrubative parameters
  real(dp),allocatable,dimension(:):: NPparam

  public::TMDR_R,TMDR_Rzeta,TMDR_R_toSL
  public:: TMDR_Initialize,TMDR_setNPparameters,LowestQ,TMDR_IsInitialized,TMDR_CurrentNPparameters
  
  public::DNP!,GammaCusp,gammaV
  
  public::zetaMUpert,zetaSL,zetaNP!,zetaMuResum
  
  
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
	orderZETA=-1
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
	if(outputLevel>0) write(*,*) &
                WarningString(' Initialize:try to set unknown ADs orders. Switch to NLO.',modulename)
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
    call MoveTO(51,'*p2  ')
    do i=1,NPlength
      read(51,*) NPparam(i)
    end do
    
    !--------------------------
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) tolerance
    
    if(outputLevel>2) write(*,'(A,ES10.3)') ' |   tolerance=',tolerance
    
    CLOSE (51, STATUS='KEEP') 
    
    call TMD_AD_Initialize()
    
    
    if(.not.QCDinput_IsInitialized()) then
      if(outputLevel>1) write(*,*) '.. initializing QCDinput (from ',moduleName,')'
      if(present(prefix)) then
      	call QCDinput_Initialize(file,prefix)
      else
	call QCDinput_Initialize(file)
      end if
    end if
    
    if(outputLevel>2) write(*,*) 'Model initialization..'
    call ModelInitialization(NPparam)
    
    started=.true.
    counter=0
    messageCounter=0
    if(outputLevel>0) write(*,*) color('----- arTeMiDe.TMDR '//trim(version)//': .... initialized',c_green)
    if(outputLevel>1) write(*,*) ' '
  end subroutine TMDR_Initialize

    !Reset the NP parameters
 subroutine TMDR_setNPparameters(lambdaIN)
    real(dp),intent(in)::lambdaIN(:)
    integer::ll
    
    ll=size(lambdaIN)
    if(ll<NPlength) then 
      if(outputLevel>0) write(*,"(A,I3,A,I3,')')")&
	      WarningString('length of lambdaNP(',moduleName),ll,&
	      color(') is less then requred (',c_red),NPlength
      if(outputLevel>0) write(*,*)color('                Rest parameters are replaced by zeros!',c_red)
      NPparam=0d0*NPparam
      NPparam(1:ll)=lambdaIN(1:ll)
    else if (ll>NPlength) then
      if(outputLevel>0) write(*,"(A,I3,A,I3,')')")&
	      WarningString('length of lambdaNP(',moduleName),ll,&
	      color(') is greater then requred (',c_red),NPlength
      if(outputLevel>0) write(*,*)color('                Array is truncated!',c_red)
      NPparam(1:NPlength)=lambdaIN(1:NPlength)
     else
      NPparam=lambdaIN
     end if
    
    if(outputLevel>2) write(*,*) 'arTeMiDe.TMDR: NPparameters reset = (',NPparam,')'
    
 end subroutine TMDR_setNPparameters

 subroutine TMDR_SetReplica(rep)
 integer::rep,i
 real(dp),dimension(1:NPlength):: InitVar
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
 real(dp),dimension(1:NPlength)::var
 var=NPparam
 
 end subroutine TMDR_CurrentNPparameters
!-------------------Anomalous dimensions------------------------------
  
  ! The TMD rapidity anomalous dimension D, pertrubative expression
  function Dpert(mu,bT,f)
    real(dp) :: bT,mu,Dpert
    integer:: f,Nf,n,k
    real(dp) ::LL,astrong,inter
    
    LL=2d0*LOG(bt*mu*C0_inv_const)
    astrong=As(mu)
    Dpert=0d0
    
    Nf=ActiveNf(mu)
    
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
    real(dp) :: gammaV,mu
    integer:: f,Nf,i    
    real(dp) ::astrong
    
    astrong=As(mu)
    
    Nf=ActiveNf(mu)
   
    gammaV=0d0
    if(f==0) then   !!gluon case
     do i=1,orderV
        gammaV=gammaV+gammaV_g(i,Nf)*(astrong**i)
     end do
    else  !! quark case
     do i=1,orderV
        gammaV=gammaV+gammaV_q(i,Nf)*(astrong**i)
     end do
    end if

  end function gammaV

  ! The Gamma Cusp anomalous dimension gammaV
  function gammaCUSP(mu,f)
    real(dp) :: gammaCUSP,mu
    integer:: f,Nf,i
    real(dp) ::astrong
    
    astrong=As(mu)
    
    Nf=activeNf(mu)
    
    gammaCUSP=0d0
    if(f==0) then   !!gluon case
     do i=1,orderCusp
        gammaCUSP=gammaCUSP+GammaCusp_g(i-1,Nf)*(astrong**i)
     end do
    else  !! quark case
     do i=1,orderCusp
        gammaCUSP=gammaCUSP+GammaCusp_q(i-1,Nf)*(astrong**i)
     end do
    end if
   
  end function gammaCUSP

!!the resummed version of rapidity anomalous dimension
function Dresum(mu,bT,f)
    real(dp)::Dresum,mu,bT
    integer::f,Nf
    integer::n,k,l
    real(dp):: X,alpha,lX,LL,commulant,inter
    
    alpha=As(mu)
    Nf=ActiveNf(mu)
    LL=2d0*LOG(bT*mu*C0_inv_const)
    
    X=betaQCD(0,Nf)*alpha*LL
    lX=Log(1d0-X)
    
    if(f==0) then !! gluon case
        commulant=lX
        
        do n=1,orderDresum
        inter=0d0
            do k=0,n
            do l=0,n
                inter=inter+(X**k)*(lX**l)*dnkl_g(n,k,l,Nf)
            end do
            end do
        commulant=commulant+((alpha/(1d0-X))**n)*inter
        end do
        
        Dresum=-GammaCusp_g(0,Nf)/betaQCD(0,Nf)*commulant
    else !!! quark case
        commulant=lX
        do n=1,orderDresum
        inter=0d0
            do k=0,n
            do l=0,n
                inter=inter+(X**k)*(lX**l)*dnkl_q(n,k,l,Nf)
            end do
            end do
        commulant=commulant+((alpha/(1d0-X))**n)*inter
        end do
        
        Dresum=-GammaCusp_q(0,Nf)/betaQCD(0,Nf)/2d0*commulant  
    end if
    
    if(ISNAN(Dresum)) then
        write(*,*) ErrorString('Dresum is NaN.',moduleName)
        write(*,*) 'At mu=',mu,'b=',bT,'Lmu=',2d0*Log(mu*bT*C0_inv_const), 'X=',X,'log(1-x)=',lX
        write(*,*) 'Evaluation STOP'
        stop
    end if
end function Dresum

!-------------------zeta-lines -------------------------------------

!! the value of zeta_mu in the pertrubation theory with ri=0
function zetaMUpert(mu,bt,f)
  real(dp):: zetaMUpert,mu,bt
  integer::f,Nf,n,k
  real(dp)::alpha,LL,val,iter
  
  if(orderZETA<0) then
    zetaMUpert=1d0
    return
  end if
  
  LL=2d0*LOG(bt*mu*C0_inv_const)
  alpha=As(mu)
  Nf=ActiveNf(mu)
  
  if(f==0) then !!! gluon
    val=vnk_g(0,0,Nf)
    
    do n=1,orderZETA
        iter = 0d0
        do k=0,orderZETA+1            
            iter=iter+(LL**k)*vnk_g(n,k,Nf)
        end do
        val=val+(alpha**n)*iter
    end do   
  else !!!! quark
    val=vnk_q(0,0,Nf)
    
    do n=1,orderZETA
        iter = 0d0
        do k=0,orderZETA+1
            iter=iter+(LL**k)*vnk_q(n,k,Nf)
        end do
        val=val+(alpha**n)*iter
    end do
  end if
  
  zetaMUpert=mu*C0_const/bT*EXP(-val)
  
end function zetaMUpert

!!!!!!!!!! exact value of zeta-line at given b,mu expanded over as
 function zetaSL(mu,b,f)
  real(dp)::zetaSL,b,mu
  integer::f
  real(dp)::dd,alpha,GD
  
  if(b<1d-6) b=1d-6
  
  if (orderZETA<0) then
    zetaSL=1d0
    return
  end if
  
  alpha=As(mu)
  dd=DNP(mu,b,f)
  if(f==0) then
    GD=valueOfGD_type4_G(dd,alpha,mu)
  else
    GD=valueOfGD_type4_Q(dd,alpha,mu)
  end if
  
  !write(*,*) '..........',b,dd,alpha,GD
  
  zetaSL=(mu**2)*exp(-GD/dd)
 
 end function zetaSL
 
!!! expression for G*D, at given D,alpha, and mu
!!! evaluated for type4 evolution for Quark
!!!!
!!!! At NNNLO I use NNLO solution. The reason is the negative values of constants at large-d. It is unphysical
!!!! To be solved...
!!!!
 function valueOfGD_type4_Q(dd,alpha,mu)
    real(dp)::valueOfGD_type4_Q,dd,alpha,mu
    integer::Nf
    real(dp)::p,pFF,z1,zm1,zm2,val
    
    if (orderZETA<0) then
        valueOfGD_type4_Q=0d0
        return
    end if
    
    Nf=ActiveNf(mu)
    pFF=pFACTOR_Q(Nf)
    p=dd*pFF
    
    z1=zFUNC(p,1)
    
    val=z1*OMEGA_q(0,1,Nf)
        
    if(orderZETA>=1) then
        val=val+alpha*(z1*OMEGA_q(1,1,Nf)+OMEGA_q(1,2,Nf)+p*OMEGA_q(1,3,Nf))
        if(orderZETA>=2) then
            zm1=zFUNC(p,-1)
            val=val+(alpha**2)*(z1*OMEGA_q(2,1,Nf)+zm1*OMEGA_q(2,2,Nf)+OMEGA_q(2,3,Nf))
!             if(orderZETA>=3) then
!                 zm2=zFUNC(p,-2)
!                 val=val+(alpha**3)*(z1*OMEGA_q(3,1,Nf)+zm1*OMEGA_q(3,2,Nf)+zm2*OMEGA_q(3,3,Nf)+OMEGA_q(3,4,Nf))
!             end if
        end if
    end if
    
    valueOfGD_type4_Q=dd/alpha*val
        
end function valueOfGD_type4_Q

!!! expression for G*D, at given D,alpha, and mu
!!! evaluated for type4 evolution for Quark
!!!!
!!!! At NNNLO I use NNLO solution. The reason is the negative values of constants at large-d. It is unphysical
!!!! To be solved...
!!!!
 function valueOfGD_type4_G(dd,alpha,mu)
    real(dp)::valueOfGD_type4_G,dd,alpha,mu
    integer::Nf
    real(dp)::p,pFF,z1,zm1,zm2,val
    
    if (orderZETA<0) then
        valueOfGD_type4_G=0d0
        return
    end if
    
    Nf=ActiveNf(mu)
    pFF=pFACTOR_G(Nf)
    p=dd*pFF
    
    z1=zFUNC(p,1)
    
    val=z1*OMEGA_g(0,1,Nf)
        
    if(orderZETA>=1) then
        val=val+alpha*(z1*OMEGA_g(1,1,Nf)+OMEGA_g(1,2,Nf)+p*OMEGA_g(1,3,Nf))
        if(orderZETA>=2) then
            zm1=zFUNC(p,-1)
            val=val+(alpha**2)*(z1*OMEGA_g(2,1,Nf)+zm1*OMEGA_g(2,2,Nf)+OMEGA_g(2,3,Nf))
!             if(orderZETA>=3) then
!                 zm2=zFUNC(p,-2)
!                 val=val+(alpha**3)*(z1*OMEGA_g(3,1,Nf)+zm1*OMEGA_g(3,2,Nf)+zm2*OMEGA_g(3,3,Nf)+OMEGA_g(3,4,Nf))
!             end if
        end if
    end if
    
    valueOfGD_type4_G=dd/alpha*val
        
end function valueOfGD_type4_G

!!! this is z=(exp(- np)-1+np)/p
function zFUNC(p,n)
    real(dp)::zFUNC,p
    integer::n
    
    if(Abs(p)>0.00001d0) then
        zFUNC=(real(n,dp)*p-1d0+Exp(-real(n,dp)*p))/p
    else
        zFUNC=real(n**2,dp)*p/2d0-real(n**3,dp)*(p**2)/6d0
    end if

end function zFUNC
 
!-------------------- R-kernels   full----------------------------------

!-----------------------Improved D picture------------------------------
!!! Evolution exponent in the improved D-picture
 function TMDR_R_type1(b,muf,zetaf,mui,zetai,mu0,f)
  real(dp)::TMDR_R_type1,b,muf,zetaf,mui,zetai,mu0
  integer::f
  
  counter=0
  
  if(b<1d-6) b=1d-6
  
  TMDR_R_type1=EXP(IntegralG1(muf,mui,zetaf,f)-(IntegralG2(mui,mu0,f)+DNP(mu0,b,f))*Log(zetaf/zetai))
  
  !write(*,*) 'TMDR_R_type1: number of AD calls ',counter

 end function TMDR_R_type1
 
 !!!This is the integral   \int_{mu_i}^{muf} dmu/mu (Gamma(mu)Log[mu^2/zeta]-gammaV)
!!! evaluation by adaptive simpson
 function IntegralG1(muf,mui,zeta,f)
  real(dp)::IntegralG1,muf,mui,zeta,lnz
  integer::f
  
  real(dp)::y1,y3,y5,X1,X3,X5,maxV
  
  lnz=Log(zeta)
  
  y1=Log(muf)
  y5=Log(mui)
  if(y1<y5) then
   y3=y1
   y1=y5
   y5=y3
  end if
  y3=(y1+y5)/2d0
  
  X1=gammaCUSP(Exp(y1),f)*(2d0*y1-lnz)-gammaV(Exp(y1),f)
  X3=gammaCUSP(Exp(y3),f)*(2d0*y3-lnz)-gammaV(Exp(y3),f)
  X5=gammaCUSP(Exp(y5),f)*(2d0*y5-lnz)-gammaV(Exp(y5),f)
  
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
   real(dp) ::y1,y5,lnz
   integer::f
   real(dp) :: interX,X1,X2,X3,X4,X5
   real(dp) :: value,valueAB,valueACB
   real(dp) :: y2,y3,y4,deltay,maxV
   
   deltay=y1-y5
   y4=y5+deltay/4d0
   y3=y5+deltay/2d0
   y2=y1-deltay/4d0
   
   
   X2=gammaCUSP(Exp(y2),f)*(2d0*y2-lnz)-gammaV(Exp(y2),f)
   X4=gammaCUSP(Exp(y4),f)*(2d0*y4-lnz)-gammaV(Exp(y4),f)
   
   counter=counter+2
   
   valueAB=deltay*(X1+4d0*X3+X5)/6d0
   valueACB=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
!    write(*,*) y1,y3,y5,valueACB-valueAB
   
   If(ABS(valueACB-valueAB)/15d0>tolerance*maxV) then
    interX=integral1_S(y1,y3,lnz,f,X1,X2,X3,maxV)&
	  +integral1_S(y3,y5,lnz,f,X3,X4,X5,maxV)
   else
    interX=valueACB
   end if
   
 end function integral1_S


 !!!This is the integral   \int_{mu_i}^{muf} dmu/mu Gamma(mu)
!!! evaluation by adaptive simpson
 function IntegralG2(muf,mui,f)
  real(dp)::IntegralG2,muf,mui
  integer::f
  
  real(dp)::y1,y3,y5,X1,X3,X5,maxV
  
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
   real(dp) ::y1,y5
   integer::f
   real(dp) :: interX,X1,X2,X3,X4,X5
   real(dp) :: value,valueAB,valueACB
   real(dp) :: y2,y3,y4,deltay,maxV
   
   deltay=y1-y5
   y4=y5+deltay/4d0
   y3=y5+deltay/2d0
   y2=y1-deltay/4d0
   
   
   X2=gammaCUSP(Exp(y2),f)
   X4=gammaCUSP(Exp(y4),f)
   
   counter=counter+2
   
   valueAB=deltay*(X1+4d0*X3+X5)/6d0
   valueACB=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
!    write(*,*) y1,y3,y5,valueACB-valueAB
   
   If(ABS(valueACB-valueAB)/15d0>tolerance*maxV) then
    interX=integral2_S(y1,y3,f,X1,X2,X3,maxV)&
	  +integral2_S(y3,y5,f,X3,X4,X5,maxV)
   else
    interX=valueACB
   end if
   
 end function integral2_S

!-----------------------Improved gamma picture------------------------------
!!! Evolution exponent in the improved gamma-picture
 function TMDR_R_type2(b,muf,zetaf,mui,zetai,f)
  real(dp)::TMDR_R_type2,b,muf,zetaf,mui,zetai
  integer::f
  
  counter=0
  
  if(b<1d-6) b=1d-6
  
  TMDR_R_type2=EXP(-IntegralG3(muf,mui,b,f)+DNP(muf,b,f)*Log(muf**2/zetaf)-DNP(mui,b,f)*Log(mui**2/zetai))
  
  !write(*,*) 'TMDR_R_type2: number of AD calls ',counter
  
  if(TMDR_R_type2>1d6) then
    write(*,*) ErrorString('Evolution factor (type2-1) is TOO HUGE check the formula',moduleName)
    write(*,*) 'NP parameters =', NPparam
    write(*,*) 'b=',b,'zetaf=',zetaf,'muf=',muf,'zetai=',zetai,'mui=',mui
    write(*,*) 'int=',IntegralG3(muf,mui,b,f),'t1=',DNP(muf,b,f)*Log(muf**2/zetaf),'t2=',DNP(mui,b,f)*Log(mui**2/zetai)
    write(*,*) 'Evaluation continue with R=10^6'
    TMDR_R_type2=1d6
  end if

 end function TMDR_R_type2
 
 !!!This is the integral   \int_{mu_i}^{muf} dmu/mu (2 D(mu,b)+gammaV)
!!! evaluation by adaptive simpson
 function IntegralG3(muf,mui,b,f)
  real(dp)::IntegralG3,muf,mui,b
  integer::f
  
  real(dp)::y1,y3,y5,X1,X3,X5,maxV
  
  y1=Log(muf)
  y5=Log(mui)
  if(y1<y5) then
   y3=y1
   y1=y5
   y5=y3
  end if
  y3=(y1+y5)/2d0
  
  X1=2d0*DNP(Exp(y1),b,f)+gammaV(Exp(y1),f)
  X3=2d0*DNP(Exp(y3),b,f)+gammaV(Exp(y3),f)
  X5=2d0*DNP(Exp(y5),b,f)+gammaV(Exp(y5),f)
  
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
   real(dp) ::y1,y5,b
   integer::f
   real(dp) :: interX,X1,X2,X3,X4,X5
   real(dp) :: value,valueAB,valueACB
   real(dp) :: y2,y3,y4,deltay,maxV
   
   deltay=y1-y5
   y4=y5+deltay/4d0
   y3=y5+deltay/2d0
   y2=y1-deltay/4d0
   
   
   X2=2d0*DNP(Exp(y2),b,f)+gammaV(Exp(y2),f)
   X4=2d0*DNP(Exp(y4),b,f)+gammaV(Exp(y4),f)
   
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
  real(dp)::TMDR_R_typeZ,b,muf,zetaf,mui,zetai
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
  real(dp)::TMDR_Rzeta_type1,b,muf,zetaf,mui,mu0,zetaP
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
  real(dp)::TMDR_Rzeta_type2,b,muf,zetaf,mui,zetaP
  integer::f
  
  counter=0
  
  if(b<1d-6) b=1d-6
  
  zetaP=zetaNP(mui,b,f)
  
  TMDR_Rzeta_type2=EXP(-IntegralG3(muf,mui,b,f)+DNP(muf,b,f)*Log(muf**2/zetaf)-DNP(mui,b,f)*Log(mui**2/zetaP))
  
  !write(*,*) 'TMDR_Rzeta_type2: number of AD calls ',counter
  
  if(TMDR_Rzeta_type2>1d6) then
    write(*,*) ErrorString('Evolution factor (type2-2) is TOO HUGE check the formula',moduleName)
    write(*,*) 'b=',b,'zetaf=',zetaf,'muf=',muf,'zetaP=',zetaP,'mui=',mui
    write(*,*) 'int=',IntegralG3(muf,mui,b,f),'t1=',DNP(muf,b,f)*Log(muf**2/zetaf),'t2=',DNP(mui,b,f)*Log(mui**2/zetaP)
    write(*,*) 'Evaluation continue with R=10^6'
    TMDR_Rzeta_type2=1d6
  end if

 end function TMDR_Rzeta_type2
 
 !!! general interface for evolutions to zeta-line for case of 3 give scales (zetaf,muf,mui)
 function TMDR_Rzeta_typeZ3(b,muf,zetaf,mui,f)
  real(dp)::TMDR_Rzeta_typeZ3,b,muf,zetaf,mui
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
  real(dp)::TMDR_Rzeta_type3,b,muf,zetaf,zetaP
  integer::f
  
  if(b<1d-6) b=1d-6
  
  zetaP=zetaNP(muf,b,f)
  
  TMDR_Rzeta_type3=EXP(-DNP(muf,b,f)*Log(zetaf/zetaP))
  
  !write(*,*) 'HERE'
  
  if(TMDR_Rzeta_type3>1d6) then
    write(*,*) ErrorString('Evolution factor(type3) is TOO HUGE check the formula',moduleName)
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
  real(dp)::TMDR_Rzeta_type4,b,muf,zetaf,mui,zetai
  integer::f
  real(dp)::dd1,alpha1,GD1 !!! at point f
  real(dp)::dd2,alpha2,GD2 !!! at point i
  
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
    write(*,*) ErrorString('Evolution factor T4 (TMD_Rzeta) is TOO HUGE check the formula',moduleName)
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
  real(dp)::TMDR_R_type4,b,muf,zetaf,mui,zetai
  integer::f
  real(dp)::dd1,alpha1,GD1 !!! at point f
  real(dp)::dd2,alpha2,GD2 !!! at point i
  
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
    write(*,*) ErrorString('Evolution factor T4 (TMD_R) is TOO HUGE check the formula',moduleName)
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
 
  !--------------------------------------------------------------------------------------------------------
 !------------------------EVOLUTION TO SPECIAL LINE-------------------------------------------------------
  !!! Evolution exponent in the exact solution to exact zeta-line
 function TMDR_R_toSL(b,muf,zetaf,f)
  real(dp)::TMDR_R_toSL,b,muf,zetaf
  integer::f
  real(dp)::dd,alpha,GD
  
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
    write(*,*) ErrorString('Evolution factor N4 is TOO HUGE check the formula',moduleName)
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
  real(dp),dimension(1:3)::LowestQ
  
  real(dp),parameter::b=25d0
  integer,parameter::f=1
  real(dp)::Qs1,Qs2,Qs3,V1,V2,V3
  
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
    write(*,*) ErrorString('LowestQ: Negative value at large Q (c=1).',moduleName)
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
    write(*,*) ErrorString('LowestQ: Negative value at large Q (c=1).',moduleName)
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
    write(*,*) ErrorString('LowestQ: Negative value at large Q (c=1).',moduleName)
    stop
  end if
 
 end function LowestQ


end module TMDR
