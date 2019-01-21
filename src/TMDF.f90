!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 1.31
!
!	Evaluation of the TMD structure function
!	
!	if you use this module please, quote 1706.01473
!
!	ver 1.31: release (AV, 30.05.2018)
!
!				A.Vladimirov (30.05.2018)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDF
  use TMDs
  use EWinput
  implicit none

  private
!   public
 
 !Current version of module
 character (len=5),parameter :: version="v1.31"

 !------------------------Physical parameters
!  real*8::MZ2=(91.d0)**2
!  real*8::GammaZ2=(2.5d0)**2
 
!------------------------------------------Tables-----------------------------------------------------------------------
    integer,parameter::Nmax=200
    INCLUDE 'Tables/BesselZero.f90'
!------------------------------------------Working variables------------------------------------------------------------
  
  logical::started=.false.
  
  !! Level of output
  !! 0=only critical
  !! 1=initialization details
  !! 2=WARNINGS
  integer::outputLevel=2
  
  logical:: convergenceLost=.false.
  
  real*8::hOGATA,tolerance
  !!!weights of ogata quadrature
  real*8,dimension(0:3,1:Nmax)::ww
  !!!nodes of ogata quadrature
  real*8,dimension(0:3,1:Nmax)::bb
  
  !Counters of calls "Integrand", Global in between setNP
  integer::GlobalCounter
  !Counter of total calls of TMD_F (reset at setNP)
  integer::CallCounter
  !Counter of maximum number of calls per integral
  integer::MaxCounter
!-----------------------------------------Public interface--------------------------------------------------------------
  public::TMDF_SetNPParameters,TMDF_SetScaleVariations,TMDF_Initialize,TMDF_ShowStatistic
  real*8,public:: TMD_F
  public::TMDF_convergenceISlost,TMDF_IsconvergenceLost
  
  real*8,allocatable::currentNP(:)
  
  interface TMDF_SetNPParameters
    module procedure TMDF_SetNPParameters,TMDF_SetNPParameters_rep
  end interface
  
 contains
 
      !!! This subroutine can be called only ones per programm run in the very beginning
  !!! It initializes TMDRm abd TMDconvolution subpackages
  !!! It set the pertrubative orders(!), so pertrubative orders cannot be changed afterwards.
  !!! It also set the paths forPDF and As greeds according to orderPDF (NLO or NNLO).
  subroutine TMDF_Initialize(orderMain)
    character(len=*)::orderMain
    character(256)::line
    integer:: j
    real*8::X
    
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
    
    
    if(outputLevel>1) write(*,*) '----- arTeMiDe.TMDF ',version,': .... initialization'
    
    do
    read(51,'(A)') line    
    if(line(1:3)=='*2 ') exit
    end do    
    do
    read(51,'(A)') line
    if(line(1:3)=='*AA') exit
    end do
    read(51,'(A)') line
    read(51,*) hOGATA
    read(51,'(A)') line
    read(51,*) tolerance
    
    if(outputLevel>2) write(*,'(A,ES8.2)') ' | h for Ogata quadrature	: ',hOGATA
    if(outputLevel>2) write(*,'(A,ES8.2)') ' | tolerance			: ',tolerance
    
    CLOSE (51, STATUS='KEEP') 
    
    if(outputLevel>1) write(*,*) 'arTeMiDe.TMDF: preparing Ogata tables'
    call PrepareTables()
    if(outputLevel>2) write(*,'(A,I4)') ' | Maximum number of nodes	:',Nmax
    if(outputLevel>1) write(*,*) 'arTeMiDe.TMDF: Ogata tables prepared'
    
    GlobalCounter=0
    CallCounter=0
    MaxCounter=0
    
    call TMDs_Initialize(orderMain)
    call EWinput_Initialize(orderMain)
    
    if(started) then
    if(outputLevel>1) write(*,*) 'arTeMiDe.TMDF already initialized'
    else      
     started=.true.
     if(outputLevel>0) write(*,*) '----- arTeMiDe.TMDF ',version,': .... initialized'
     if(outputLevel>1) write(*,*) ' '
    end if
    
  end subroutine TMDF_Initialize

  !!!!!!!Functions which carry the trigger on convergences.... Its used in xSec, and probably in other places.
  function TMDF_IsconvergenceLost()
  logical::TMDF_IsconvergenceLost
  !!! checks TMDs trigger
  if(TMDs_IsconvergenceLost()) convergenceLost=.true.
  
  TMDF_IsconvergenceLost=convergenceLost
  end function TMDF_IsconvergenceLost
  
  subroutine TMDF_convergenceISlost()  
  convergenceLost=.true.
  if(outputLevel>1) write(*,*) 'arTeMiDe.TMDF: convergence triger set to be lost.'
  end subroutine TMDF_convergenceISlost
  
  subroutine TMDF_ShowStatistic()
    if(convergenceLost) then
      write(*,*) '         TMDF statistics: convergence has been lost.'
    else
      write(*,'(A,ES12.3)') 'TMDF statistics               total calls of TMDs  :  ',Real(2*GlobalCounter)
      write(*,'(A,ES12.3)') '                              total calls of TMD_F :  ',Real(CallCounter)
      write(*,'(A,F12.3)')  '                                         avarage M :  ',Real(GlobalCounter)/Real(CallCounter)
      write(*,'(A,I12)')    '                                     maximum calls :  ',MaxCounter
	
    end if
  end subroutine TMDF_ShowStatistic
  
  !passes the NP parameters to TMDs
  subroutine TMDF_SetNPParameters(lambda)
    real*8 :: lambda(:)
   
   !!some output statistic
   if(outputLevel>2) call TMDF_ShowStatistic()
   
   !!! save current parameters
   if(.not. allocated(currentNP)) allocate(currentNP(1:size(lambda)))   
   currentNP=lambda
   
   convergenceLost=.false.
   
   GlobalCounter=0
   CallCounter=0
   MaxCounter=0
   
   call TMDs_setNPparameters(lambda)

  end subroutine TMDF_SetNPParameters
  
   !passes the NP parameters to TMDs
  subroutine TMDF_SetNPParameters_rep(num)
   integer::num
   
   !!some output statistic
   if(outputLevel>2) call TMDF_ShowStatistic()
   
   !!! save current parameters
   if(.not. allocated(currentNP)) allocate(currentNP(1:1))   
   currentNP(1)=real(num)
   
   convergenceLost=.false.
   
   GlobalCounter=0
   CallCounter=0
   MaxCounter=0
   
   call TMDs_setNPparameters(num)

  end subroutine TMDF_SetNPParameters_rep
  
!passes the scale variation to TMDs
 subroutine TMDF_SetScaleVariations(c1_in,c3_in,c4_in)
    real*8::c1_in,c3_in,c4_in
    
    call TMDs_SetScaleVariations(c1_in,c3_in,c4_in)
    
 end subroutine TMDF_SetScaleVariations
 
 !!!Prepare tables for Ogata quadrature with given h
 subroutine PrepareTables()
  real*8,parameter::piHalf=1.5707963267948966d0
  real*8,parameter::pi=3.141592653589793d0
  integer::i,k
  real*8::t!=h*xi
  real*8::psiPart!=tanh[pi/2 Sinh[h xi]]
  
  do k=0,3
  do i=1,Nmax
    t=hOGATA*JZero(k,i)
    psiPart=Tanh(piHalf*Sinh(t))
      bb(k,i)=JZero(k,i)*psiPart
      ww(k,i)=BESSEL_JN(k,bb(k,i))/JZero(k,i)/(BESSEL_JN(k+1,JZero(k,i))**2)*(pi*t*Cosh(t)+Sinh(pi*Sinh(t)))/(1d0+Cosh(pi*Sinh(t)))
!      write(*,*) psiPart,b(k,i),w(k,i)
  end do
  end do
 
 end subroutine PrepareTables
 
 !!!This is the defining module function
 !!! It evaluates the integral 
 !!!  int_0^infty   b db/2  Jn(b qT) zff F1 F2
 !!!
 function TMD_F(Q2,qT,x1,x2,mu,zeta1,zeta2,process)
  real*8::qT,x1,x2,mu,zeta1,zeta2,Q2
  integer::process
  real*8::integral,eps
  real*8::v1,v2,v3,v4
  integer::k,n
  
  CallCounter=CallCounter+1
  integral=0d0
  
  if(qT<0.0000001d0) then  
  integral=0d0
  else
  !!!in the case of lost convergence we return huge number (divergent xSec)
  if(TMDF_IsconvergenceLost()) then	
	TMD_F=1d10		
  else
  
  v1=1d0
  v2=1d0
  v3=1d0
  v4=1d0
  
  !!Here we set the order of Bessel
  if(process<10000) then
  n=0
  else if(process<20000) then
  n=1
  else if(process<30000) then
  n=2
  else
  n=3
  end if
  
  do k=1,Nmax!!! maximum of number of bessel roots preevaluated in the head
    eps=ww(n,k)*(bb(n,k)**(n+1))*Integrand(Q2,bb(n,k)/qT,x1,x2,mu,zeta1,zeta2,process)
    
    integral=integral+eps
!     if(k>8) then
      v4=v3
      v3=v2
      v2=v1
      v1=ABS(eps)
!       if(v1+v2+v3+v4>0.7d0*ABS(integral)) write(17,*) bb(n,k)/qT,x1,x2
!       write(*,*) k, eps,ww(n,k),eps/ww(n,k)
      if(v1+v2+v3+v4<=tolerance*ABS(integral)) exit
!     end if
  end do
  if(k>=Nmax) then	
    if(outputlevel>0) WRITE(*,*) 'WARNING arTeMiDe.TMDF: OGATA quadrature diverge. TMD decaing too slow? '
      if(outputlevel>1) then
      write(*,*) 'Current set of NP parameters ------------'
      write(*,*) currentNP
      write(*,*) 'Information over the last call ----------'
      write(*,*) 'bt/qT= ',bb(n,Nmax)/qT, 'qT=',qT
      write(*,*) 'W=',Integrand(Q2,bb(n,Nmax)/qT,x1,x2,mu,zeta1,zeta2,process), 'eps/integral =', eps/integral
      write(*,*) 'v1+v2+v3+v4=',v1+v2+v3+v4, '>',tolerance*ABS(integral)
      write(*,*) '(x1,x2)=(',x1,',',x2,')'
      write(*,*) 'process =',process,' it is ',CallCounter,'call.'
      write(*,*) '------------------------------------------'
      end if
    call TMDF_convergenceISlost()
  end if
  if(k>MaxCounter) MaxCounter=k-1
!   write(*,*) 'Integral=',integral
  TMD_F=integral/(qT**(n+2)) 
  end if 
  end if
  !write(*,*) 'Last call: ',k
 end function TMD_F
 
 function Integrand(Q2,b,x1,x2,mu,zeta1,zeta2,process)
 real*8::Integrand
 real*8::b,x1,x2,mu,zeta1,zeta2,Q2
 integer::process
 real*8,dimension(-5:5)::FA,FB,FAB
   !!cross-seciton parameters
!   real*8,parameter:: paramU=0.40329064872689d0 !! ((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) for eq=2/3
!   real*8,parameter:: paramD=0.51983027428079d0 !! ((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) for eq=1/3
!   real*8,parameter:: paramS=0.51983027428079d0 !! ((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) for eq=1/3
!   real*8,parameter:: paramC=0.40329064872689d0 !! ((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) for eq=2/3
!   real*8,parameter:: paramB=0.51983027428079d0 !! ((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) for eq=1/3
 
 !increment counter 
 GlobalCounter=GlobalCounter+1
 
 SELECT CASE(process)
  !!!test case
  CASE(0,10000,20000,30000)
    Integrand=Exp(-0.2d0*b)
!--------------------------------------------------------------------------------
  CASE (1) !pp->gamma
	! e_q^2 *F_q(A)*F_qbar(B)
	if(zeta1==zeta2) then
	 FAB=uPDF_uPDF(x1,x2,b,mu,zeta1,1,1)
	else
	 FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	 FB=uTMDPDF_5(x2,b,mu,zeta2,1)
	 FAB=FA*(FB(5:-5:-1))
	end if
	
	Integrand=FAB(1)/9.d0&
	  +FAB(2)*4.d0/9.d0&
	  +FAB(3)/9.d0&
	  +FAB(4)/9.d0&
	  +FAB(5)*4d0/9d0&
	  +FAB(-1)/9.d0&
	  +FAB(-2)*4.d0/9.d0&
	  +FAB(-3)/9.d0&
	  +FAB(-4)/9.d0&
	  +FAB(-5)*4d0/9d0
	
! 	Integrand=FA(1)*FB(-1)/9.d0&
! 	  +FA(2)*FB(-2)*4.d0/9.d0&
! 	  +FA(3)*FB(-3)/9.d0&
! 	  +FA(4)*FB(-4)/9.d0&
! 	  +FA(5)*FB(-5)*4d0/9d0&
! 	  +FA(-1)*FB(1)/9.d0&
! 	  +FA(-2)*FB(2)*4.d0/9.d0&
! 	  +FA(-3)*FB(3)/9.d0&
! 	  +FA(-4)*FB(4)/9.d0&
! 	  +FA(-5)*FB(5)*4d0/9d0
!--------------------------------------------------------------------------------  
  CASE (2) !ppbar->gamma
	! e_q^2 *F_q(A)*F_q(B)
	if(zeta1==zeta2) then
	 FAB=uPDF_anti_uPDF(x1,x2,b,mu,zeta1,1,1)
	else
	 FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	 FB=uTMDPDF_5(x2,b,mu,zeta2,1)
	 FAB=FA*FB
	end if
	
	Integrand=FAB(1)/9.d0&
	  +FAB(2)*4.d0/9.d0&
	  +FAB(3)/9.d0&
	  +FAB(4)/9.d0&
	  +FAB(5)*4d0/9d0&
	  +FAB(-1)/9.d0&
	  +FAB(-2)*4.d0/9.d0&
	  +FAB(-3)/9.d0&
	  +FAB(-4)/9.d0&
	  +FAB(-5)*4d0/9d0
	
! 	Integrand=FA(1)*FB(1)/9.d0&
! 	  +FA(2)*FB(2)*4.d0/9.d0&
! 	  +FA(3)*FB(3)/9.d0&
! 	  +FA(4)*FB(4)/9.d0&
! 	  +FA(5)*FB(5)*4d0/9d0&
! 	  +FA(-1)*FB(-1)/9.d0&
! 	  +FA(-2)*FB(-2)*4.d0/9.d0&
! 	  +FA(-3)*FB(-3)/9.d0&
! 	  +FA(-4)*FB(-4)/9.d0&
! 	  +FA(-5)*FB(-5)*4d0/9d0
!--------------------------------------------------------------------------------  
  CASE (3) !pp->Z
	  !((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) *F_q(A)*F_qbar(B)
	if(zeta1==zeta2) then
	 FAB=uPDF_uPDF(x1,x2,b,mu,zeta1,1,1)
	else
	 FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	 FB=uTMDPDF_5(x2,b,mu,zeta2,1)
	 FAB=FA*(FB(5:-5:-1))
	end if
	
	Integrand=&
	  FAB(1)*paramD&
	  +FAB(2)*paramU&
	  +FAB(3)*paramS&
	  +FAB(4)*paramC&
	  +FAB(5)*paramB&
	  +FAB(-1)*paramD&
	  +FAB(-2)*paramU&
	  +FAB(-3)*paramS&
	  +FAB(-4)*paramC&
	  +FAB(-5)*paramB
	
! 	Integrand=&
! 	  FA(1)*FB(-1)*paramD&
! 	  +FA(2)*FB(-2)*paramU&
! 	  +FA(3)*FB(-3)*paramS&
! 	  +FA(4)*FB(-4)*paramC&
! 	  +FA(5)*FB(-5)*paramB&
! 	  +FA(-1)*FB(1)*paramD&
! 	  +FA(-2)*FB(2)*paramU&
! 	  +FA(-3)*FB(3)*paramS&
! 	  +FA(-4)*FB(4)*paramC&
! 	  +FA(-5)*FB(5)*paramB
!--------------------------------------------------------------------------------  
  CASE (4) !ppbar->Z
	  !((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) *F_q(A)*F_qbar(B)
	if(zeta1==zeta2) then
	 FAB=uPDF_anti_uPDF(x1,x2,b,mu,zeta1,1,1)
	else
	 FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	 FB=uTMDPDF_5(x2,b,mu,zeta2,1)
	 FAB=FA*FB
	end if
	
	Integrand=&
	  FAB(1)*paramD&
	  +FAB(2)*paramU&
	  +FAB(3)*paramS&
	  +FAB(4)*paramC&
	  +FAB(5)*paramB&
	  +FAB(-1)*paramD&
	  +FAB(-2)*paramU&
	  +FAB(-3)*paramS&
	  +FAB(-4)*paramC&
	  +FAB(-5)*paramB
	
! 	Integrand=&
! 	  FA(1)*FB(1)*paramD&
! 	  +FA(2)*FB(2)*paramU&
! 	  +FA(3)*FB(3)*paramS&
! 	  +FA(4)*FB(4)*paramC&
! 	  +FA(5)*FB(5)*paramB&
! 	  +FA(-1)*FB(-1)*paramD&
! 	  +FA(-2)*FB(-2)*paramU&
! 	  +FA(-3)*FB(-3)*paramS&
! 	  +FA(-4)*FB(-4)*paramC&
! 	  +FA(-5)*FB(-5)*paramB
!--------------------------------------------------------------------------------  
  CASE (5) !pp->Z+gamma
	if(zeta1==zeta2) then
	 FAB=uPDF_uPDF(x1,x2,b,mu,zeta1,1,1)
	else
	 FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	 FB=uTMDPDF_5(x2,b,mu,zeta2,1)
	 FAB=FA*(FB(5:-5:-1))
	end if
	
	Integrand=XIntegrandForDYwithZgamma(FAB,Q2)
!--------------------------------------------------------------------------------  
  CASE (6) !ppbar->Z+gamma
	if(zeta1==zeta2) then
	 FAB=uPDF_anti_uPDF(x1,x2,b,mu,zeta1,1,1)
	else
	 FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	 FB=uTMDPDF_5(x2,b,mu,zeta2,1)
	 FAB=FA*FB
	end if
	
	!! we invert the order of FB
	Integrand=XIntegrandForDYwithZgamma(FAB,Q2)
!--------------------------------------------------------------------------------  
  CASE (1001) !p+Cu->gamma* !!this is for E288
	FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	FB=uTMDPDF_5(x2,b,mu,zeta2,1)
	Integrand=116d0/567d0*(FA(2)*FB(-2)+FA(-2)*FB(2))+136d0/567d0*(FA(-2)*FB(1)+FA(2)*FB(-1))&
	      +34d0/567d0*(FA(-1)*FB(2)+FA(1)*FB(-2))+29d0/567d0*(FA(-1)*FB(1)+FA(1)*FB(-1))&
	      +1d0/9d0*(FA(-3)*FB(3)+FA(3)*FB(-3)+4d0*FA(-4)*FB(4)+4d0*FA(4)*FB(-4)+FA(-5)*FB(5)+FA(5)*FB(-5))
	      
  !--------------------------------------------------------------------------------  
  CASE (1002) !p+2H->gamma* !!this is for E772
	FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	FB=uTMDPDF_5(x2,b,mu,zeta2,1)
	Integrand=2d0/9d0*(FA(2)*FB(-2)+FA(-2)*FB(2))+2d0/9d0*(FA(-2)*FB(1)+FA(2)*FB(-1))&
	      +1d0/18d0*(FA(-1)*FB(2)+FA(1)*FB(-2))+1d0/18d0*(FA(-1)*FB(1)+FA(1)*FB(-1))&
	      +1d0/9d0*(FA(-3)*FB(3)+FA(3)*FB(-3)+4d0*FA(-4)*FB(4)+4d0*FA(4)*FB(-4)+FA(-5)*FB(5)+FA(5)*FB(-5))
  !----------------------------------------------------------------------------------
  !-------------------------SIDIS----------------------------------------------------
  !----------------------------------------------------------------------------------
    CASE (2001) !ppbar->gamma
	! e_q^2 *F_q(A)*F_q(B)
	FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	FB=uTMDFF_5(x2,b,mu,zeta2,1)
	
	Integrand=FA(1)*FB(1)/9.d0&
	  +FA(2)*FB(2)*4.d0/9.d0&
	  +FA(3)*FB(3)/9.d0&
	  +FA(4)*FB(4)/9.d0&
	  +FA(5)*FB(5)*4d0/9d0&
	  +FA(-1)*FB(-1)/9.d0&
	  +FA(-2)*FB(-2)*4.d0/9.d0&
	  +FA(-3)*FB(-3)/9.d0&
	  +FA(-4)*FB(-4)/9.d0&
	  +FA(-5)*FB(-5)*4d0/9d0
  CASE DEFAULT
    write(*,*) 'ERROR:arTeMiDe.TMDF: undefined process: ',process
    write(*,*) 'Evaluation stop'
    stop
 END SELECT
 
  if(ISNAN(Integrand)) then
   write(*,*) 'arTeMiDe TMDF: CRITICAL ERROR. Integrand evaluated to NaN'
   write(*,*) 'bT=',b, 'x1,x2=',x1,x2,' process=',process
   write(*,*) 'mu=',mu, 'Q2=',Q2
   write(*,*) 'Current set of NP parameters ------------'
   write(*,*) currentNP
   write(*,*) 'arTeMiDe: ConvergenceLOST trigger ON'
   call TMDF_convergenceISlost()
   Integrand=1d10
   end if
  
   if(Integrand>1d32) then
   write(*,*) 'arTeMiDe TMDF: CRITICAL ERROR. Integrand evaluated to >10^32'
   write(*,*) 'bT=',b, 'x1,x2=',x1,x2,' process=',process
   write(*,*) 'mu=',mu, 'Q2=',Q2
   write(*,*) 'Current set of NP parameters ------------'
   write(*,*) currentNP
   write(*,*) 'arTeMiDe: onvergenceLOST trigger ON'
   call TMDF_convergenceISlost()
   Integrand=1d10
   end if
 
 end function Integrand

!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!!! The hadron tensonr for the DY icludes Z + gamma, evaluated at FA and FB 
function XIntegrandForDYwithZgamma(FAB,Q2)
     real*8::XIntegrandForDYwithZgamma,Q2
    !!cross-seciton parameters
     real*8,dimension(-5:5):: FAB
     
     !!!parameters of Z boson coupling
!      real*8,parameter:: paramU=0.40329064872689d0 !! ((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) for eq=2/3
!      real*8,parameter:: paramD=0.51983027428079d0 !! ((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) for eq=1/3
!      real*8,parameter:: paramS=0.51983027428079d0  !! ((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) for eq=1/3
!      real*8,parameter:: paramC=0.40329064872689d0 !! ((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) for eq=2/3
!      real*8,parameter:: paramB=0.51983027428079d0  !! ((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) for eq=1/3
!      real*8,parameter:: paramL=0.35358707798999d0  !! ((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) for eq=1
     
     !!!parameters of Z-gamma boson coupling
!      real*8,parameter:: paramMIXU=0.1515661518957d0 !! e(T3-2e sW^2)/2sWcW for eq=+2/3,, T3=+1/2
!      real*8,parameter:: paramMIXD=0.1367184036034d0 !! e(T3-2e sW^2)/2sWcW for eq=-1/3,T3=-1/2
!      real*8,parameter:: paramMIXS=0.1367184036034d0  !! e(T3-2e sW^2)/2sWcW for eq=-1/3, T3=-1/2
!      real*8,parameter:: paramMIXC=0.1515661518957d0 !! e(T3-2e sW^2)/2sWcW for eq=+2/3,, T3=+1/2
!      real*8,parameter:: paramMIXB=0.1367184036034d0  !! e(T3-2e sW^2)/2sWcW for eq=-1/3, T3=-1/2
!      real*8,parameter:: paramMIXL=0.0445432448766d0  !! e(T3-2e sW^2)/2sWcW for eq=-1, T3=-1/2
     
     XIntegrandForDYwithZgamma=&
     (&!gamma-part
	  4d0/9d0*FAB(2)&
	  +1d0/9d0*FAB(1)&
	  +1d0/9d0*FAB(3)&
	  +4d0/9d0*FAB(4)&
	  +1d0/9d0*FAB(5)&
	  +4d0/9d0*FAB(-2)&
	  +1d0/9d0*FAB(-1)&
	  +1d0/9d0*FAB(-3)&
	  +4d0/9d0*FAB(-4)&
	  +1d0/9d0*FAB(-5))&
     +&!gamma-Z interference
     paramMIXL*(&
	  paramMIXU*FAB(2)&
	  +paramMIXD*FAB(1)&
	  +paramMIXS*FAB(3)&
	  +paramMIXC*FAB(4)&
	  +paramMIXB*FAB(5)&
	  +paramMIXU*FAB(-2)&
	  +paramMIXD*FAB(-1)&
	  +paramMIXS*FAB(-3)&
	  +paramMIXC*FAB(-4)&
	  +paramMIXB*FAB(-5))*&
	  2d0*Q2*(Q2-MZ2)/((Q2-MZ2)**2+GammaZ2*MZ2)&
     +&!ZZ-contributions
     paramL*(&
	  paramU*FAB(2)&
	  +paramD*FAB(1)&
	  +paramS*FAB(3)&
	  +paramC*FAB(4)&
	  +paramB*FAB(5)&
	  +paramU*FAB(-2)&
	  +paramD*FAB(-1)&
	  +paramS*FAB(-3)&
	  +paramC*FAB(-4)&
	  +paramB*FAB(-5))*&
	  Q2*Q2/((Q2-MZ2)**2+GammaZ2*MZ2)
     
     
!      XIntegrandForDYwithZgamma=&
!      (&!gamma-part
! 	  4d0/9d0*FA(2)*FB(-2)&
! 	  +1d0/9d0*FA(1)*FB(-1)&
! 	  +1d0/9d0*FA(3)*FB(-3)&
! 	  +4d0/9d0*FA(4)*FB(-4)&
! 	  +1d0/9d0*FA(5)*FB(-5)&
! 	  +4d0/9d0*FA(-2)*FB(2)&
! 	  +1d0/9d0*FA(-1)*FB(1)&
! 	  +1d0/9d0*FA(-3)*FB(3)&
! 	  +4d0/9d0*FA(-4)*FB(4)&
! 	  +1d0/9d0*FA(-5)*FB(5))&
!      +&!gamma-Z interference
!      paramMIXL*(&
! 	  paramMIXU*FA(2)*FB(-2)&
! 	  +paramMIXD*FA(1)*FB(-1)&
! 	  +paramMIXS*FA(3)*FB(-3)&
! 	  +paramMIXC*FA(4)*FB(-4)&
! 	  +paramMIXB*FA(5)*FB(-5)&
! 	  +paramMIXU*FA(-2)*FB(2)&
! 	  +paramMIXD*FA(-1)*FB(1)&
! 	  +paramMIXS*FA(-3)*FB(3)&
! 	  +paramMIXC*FA(-4)*FB(4)&
! 	  +paramMIXB*FA(-5)*FB(5))*&
! 	  2d0*Q2*(Q2-MZ2)/((Q2-MZ2)**2+GammaZ2*MZ2)&
!      +&!ZZ-contributions
!      paramL*(&
! 	  paramU*FA(2)*FB(-2)&
! 	  +paramD*FA(1)*FB(-1)&
! 	  +paramS*FA(3)*FB(-3)&
! 	  +paramC*FA(4)*FB(-4)&
! 	  +paramB*FA(5)*FB(-5)&
! 	  +paramU*FA(-2)*FB(2)&
! 	  +paramD*FA(-1)*FB(1)&
! 	  +paramS*FA(-3)*FB(3)&
! 	  +paramC*FA(-4)*FB(4)&
! 	  +paramB*FA(-5)*FB(5))*&
! 	  Q2*Q2/((Q2-MZ2)**2+GammaZ2*MZ2)
     
end function XIntegrandForDYwithZgamma
 
 
end module TMDF

! program example
!   use TMDF
!   implicit none
!   
!   call TMDF_Initialize("NLO")
!   
!   call TMDF_SetNPParameters((/1.67d0,0.327d0,0.112d0,0.828d0/))
!   
!   write(*,*) TMD_F(100d0,5d0,0.5d0,0.5d0,100d0,10000d0,1000d0,1)
!   write(*,*) TMD_F(100d0,5d0,0.5d0,0.5d0,100d0,10000d0,1000d0,2)
!   write(*,*) TMD_F(100d0,5d0,0.5d0,0.5d0,100d0,10000d0,1000d0,3)
!   write(*,*) TMD_F(100d0,5d0,0.5d0,0.5d0,100d0,10000d0,1000d0,4)
!   write(*,*) TMD_F(100d0,5d0,0.5d0,0.5d0,100d0,10000d0,1000d0,5)
!   write(*,*) TMD_F(100d0,5d0,0.5d0,0.5d0,100d0,10000d0,1000d0,6)
!   
!   call TMDF_ShowStatistic()
! 
! end program example