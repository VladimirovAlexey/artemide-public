!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 1.41
!
!	Evaluation of the TMD structure function
!	
!	if you use this module please, quote 1706.01473
!
!	ver 1.31: release (AV, 30.05.2018)
!	ver 1.41: fixed potential bug in the initialisation order (AV, 28.02.2019)
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
 
 character (len=7),parameter :: moduleName="TMDR"
 character (len=5),parameter :: version="v2.00"

!------------------------------------------Tables-----------------------------------------------------------------------
    integer,parameter::Nmax=200
    INCLUDE 'Tables/BesselZero.f90'
!------------------------------------------Working variables------------------------------------------------------------
  
  integer::outputLevel=2
  !! variable that count number of WRNING mesagges. In order not to spam too much
  integer::messageTrigger=6
  logical::started=.false.
  
  logical:: convergenceLost=.false.
  
  real*8::hOGATA,tolerance
  !!!weights of ogata quadrature
  real*8,dimension(0:3,1:Nmax)::ww,ww0
  !!!nodes of ogata quadrature
  real*8,dimension(0:3,1:Nmax)::bb,bb0
  
  integer::GlobalCounter
  integer::CallCounter
  integer::MaxCounter
  integer::messageCounter
!-----------------------------------------Public interface--------------------------------------------------------------
  public::TMDF_Initialize,TMDF_ShowStatistic,TMDF_ResetCounters
  public:: TMDF_F
  public::TMDF_convergenceISlost,TMDF_IsconvergenceLost,TMDF_IsInitialized
  
 contains
  function TMDF_IsInitialized()
  logical::TMDF_IsInitialized
  TMDF_IsInitialized=started
  end function TMDF_IsInitialized

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

   !! Initialization of the package
  subroutine TMDF_Initialize(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=300)::path,line
    logical::initRequared
    
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
    call MoveTO(51,'*p2  ')
    read(51,*) outputLevel    
    if(outputLevel>2) write(*,*) '--------------------------------------------- '
    if(outputLevel>2) write(*,*) 'artemide.',moduleName,version,': initialization started ... '
    call MoveTO(51,'*p3  ')
    read(51,*) messageTrigger
    
    call MoveTO(51,'*7   ')
    call MoveTO(51,'*p1  ')
    read(51,*) initRequared
    if(.not.initRequared) then
      if(outputLevel>2) write(*,*)'artemide.',moduleName,': initialization is not requared. '
      started=.false.
      return
    end if
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) tolerance
    call MoveTO(51,'*p2  ')
    read(51,*) hOGATA
      
      if(outputLevel>2) write(*,'(A,ES8.2)') ' | h for Ogata quadrature	: ',hOGATA
      if(outputLevel>2) write(*,'(A,ES8.2)') ' | tolerance			: ',tolerance
      
      CLOSE (51, STATUS='KEEP') 
      
      if(outputLevel>1) write(*,*) 'arTeMiDe.TMDF: preparing Ogata tables'
      call PrepareTables()
      if(outputLevel>2) write(*,'(A,I4)') ' | Maximum number of nodes	:',Nmax
      if(outputLevel>1) write(*,*) 'arTeMiDe.TMDF: Ogata tables prepared'
      
      convergenceLost=.false.
      GlobalCounter=0
      CallCounter=0
      MaxCounter=0
      messageCounter=0
      
      if(.not.TMDs_IsInitialized()) then
	if(outputLevel>1) write(*,*) '.. initializing TMDs (from ',moduleName,')'
	if(present(prefix)) then
	  call TMDs_Initialize(file,prefix)
	else
	  call TMDs_Initialize(file)
	end if
      end if
      
      if(.not.EWinput_IsInitialized()) then
	if(outputLevel>1) write(*,*) '.. initializing EWinput (from ',moduleName,')'
	if(present(prefix)) then
	  call EWinput_Initialize(file,prefix)
	else
	  call EWinput_Initialize(file)
	end if
      end if
      
      
      started=.true.
      if(outputLevel>0) write(*,*) '----- arTeMiDe.TMDF ',version,': .... initialized'
      if(outputLevel>1) write(*,*) ' '
    
  end subroutine TMDF_Initialize

  !!!!!!!Functions which carry the trigger on convergences.... Its used in xSec, and probably in other places.
  function TMDF_IsconvergenceLost()
  logical::TMDF_IsconvergenceLost
  !!! checks TMDs trigger
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
      write(*,'(A,ES12.3)') '                              total calls of TMDF_F :  ',Real(CallCounter)
      write(*,'(A,F12.3)')  '                                         avarage M :  ',Real(GlobalCounter)/Real(CallCounter)
      write(*,'(A,I12)')    '                                     maximum calls :  ',MaxCounter
	
    end if
  end subroutine TMDF_ShowStatistic
  
  !passes the NP parameters to TMDs
  subroutine TMDF_ResetCounters()
   if(outputLevel>2) call TMDF_ShowStatistic()
   
   convergenceLost=.false.
   
   GlobalCounter=0
   CallCounter=0
   MaxCounter=0
   messageCounter=0

  end subroutine TMDF_ResetCounters
  

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
  
  !!these are tables for the step h=h*0.05, they are used in the case of small q_T
  do k=0,3
  do i=1,Nmax
    t=hOGATA*JZero(k,i)*0.05d0
    psiPart=Tanh(piHalf*Sinh(t))
    bb0(k,i)=JZero(k,i)*psiPart
    ww0(k,i)=BESSEL_JN(k,bb0(k,i))/JZero(k,i)/(BESSEL_JN(k+1,JZero(k,i))**2)*(pi*t*Cosh(t)+Sinh(pi*Sinh(t)))/(1d0+Cosh(pi*Sinh(t)))
  end do
  end do
 
 end subroutine PrepareTables
 
 !!!This is the defining module function
 !!! It evaluates the integral 
 !!!  int_0^infty   b db/2  Jn(b qT) zff F1 F2
 !!!
 function TMDF_F(Q2,qT,x1,x2,mu,zeta1,zeta2,process)
  real*8::TMDF_F
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
	TMDF_F=1d10		
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
  
  if(qT>1d0) then
  do k=1,Nmax!!! maximum of number of bessel roots preevaluated in the head
    eps=ww(n,k)*(bb(n,k)**(n+1))*Integrand(Q2,bb(n,k)/qT,x1,x2,mu,zeta1,zeta2,process)
    
    integral=integral+eps
    
!     write(*,*) k,bb(n,k)/qT,eps,integral
    
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
      !write(*,*) 'Current set of NP parameters ------------'
      !write(*,*) currentNP
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
  
  else !!! in the case of small q_T we use smaller step for ogata
  do k=1,Nmax!!! maximum of number of bessel roots preevaluated in the head
    eps=ww0(n,k)*(bb0(n,k)**(n+1))*Integrand(Q2,bb0(n,k)/qT,x1,x2,mu,zeta1,zeta2,process)
    
    integral=integral+eps
    
!     write(*,*) k,bb0(n,k)/qT,eps,integral
    
!     if(k>8) then
      v4=v3
      v3=v2
      v2=v1
      v1=ABS(eps)
!       if(v1+v2+v3+v4>0.7d0*ABS(integral)) write(17,*) bb0(n,k)/qT,x1,x2
!       write(*,*) k, eps,ww0(n,k),eps/ww0(n,k)
      if(v1+v2+v3+v4<=tolerance*ABS(integral)) exit
!     end if
  end do
  if(k>=Nmax) then	
    if(outputlevel>0) WRITE(*,*) 'WARNING arTeMiDe.TMDF: OGATA quadrature diverge. TMD decaing too slow? '
      if(outputlevel>1) then
      !write(*,*) 'Current set of NP parameters ------------'
      !write(*,*) currentNP
      write(*,*) 'Information over the last call ----------'
      write(*,*) 'bt/qT= ',bb0(n,Nmax)/qT, 'qT=',qT
      write(*,*) 'W=',Integrand(Q2,bb0(n,Nmax)/qT,x1,x2,mu,zeta1,zeta2,process), 'eps/integral =', eps/integral
      write(*,*) 'v1+v2+v3+v4=',v1+v2+v3+v4, '>',tolerance*ABS(integral)
      write(*,*) '(x1,x2)=(',x1,',',x2,')'
      write(*,*) 'process =',process,' it is ',CallCounter,'call.'
      write(*,*) '------------------------------------------'
      end if
    call TMDF_convergenceISlost()
  end if
  end if
  
  if(k>MaxCounter) MaxCounter=k-1
!   write(*,*) 'Integral=',integral
  TMDF_F=integral/(qT**(n+2)) 
  end if 
  end if
  !write(*,*) 'Last call: ',k
 end function TMDF_F
 
 function Integrand(Q2,b,x1,x2,mu,zeta1,zeta2,process)
 real*8::Integrand
 real*8::b,x1,x2,mu,zeta1,zeta2,Q2
 integer::process,h
 real*8,dimension(-5:5)::FA,FB,FAB
 
 !increment counter 
 GlobalCounter=GlobalCounter+1
 
 if(b>1000d0) then
  Integrand=0d0
  return
 end if
 
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
	  +FAB(4)*4d0/9.d0&
	  +FAB(5)/9d0&
	  +FAB(-1)/9.d0&
	  +FAB(-2)*4.d0/9.d0&
	  +FAB(-3)/9.d0&
	  +FAB(-4)*4d0/9.d0&
	  +FAB(-5)/9d0
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
	  +FAB(4)*4d0/9.d0&
	  +FAB(5)/9d0&
	  +FAB(-1)/9.d0&
	  +FAB(-2)*4.d0/9.d0&
	  +FAB(-3)/9.d0&
	  +FAB(-4)*4d0/9.d0&
	  +FAB(-5)/9d0
	
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
  CASE (7) !pp-> W+
	 FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	 FB=uTMDPDF_5(x2,b,mu,zeta2,1)
	 
	Integrand=paramW_L*(&
	paramW_UD*(FA(2)*FB(-1)+FA(-1)*FB(2))&		!u*dbar+dbar*u
	+paramW_US*(FA(2)*FB(-3)+FA(-3)*FB(2))&		!u*sbar+sbar*u
	+paramW_UB*(FA(2)*FB(-5)+FA(-5)*FB(2))&		!u*bbar+bbar*u
	+paramW_CD*(FA(4)*FB(-1)+FA(-1)*FB(4))&		!c*dbar+dbar*c
	+paramW_CS*(FA(4)*FB(-3)+FA(-3)*FB(4))&		!c*sbar+sbar*c
	+paramW_CB*(FA(4)*FB(-5)+FA(-5)*FB(4))&		!c*bbar+bbar*c
	)*Q2*Q2/((Q2-MW2)**2+GammaW2*MW2)
!--------------------------------------------------------------------------------  
  CASE (8) !pp-> W-
	 FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	 FB=uTMDPDF_5(x2,b,mu,zeta2,1)
	 
	Integrand=paramW_L*(&
	paramW_UD*(FA(1)*FB(-2)+FA(-2)*FB(1))&		!d*ubar+ubar*d
	+paramW_US*(FA(3)*FB(-2)+FA(-2)*FB(3))&		!s*ubar+ubar*s
	+paramW_UB*(FA(5)*FB(-2)+FA(-2)*FB(5))&		!b*ubar+ubar*b
	+paramW_CD*(FA(1)*FB(-4)+FA(-4)*FB(1))&		!d*cbar+cbar*d
	+paramW_CS*(FA(3)*FB(-4)+FA(-4)*FB(3))&		!s*cbar+cbar*s
	+paramW_CB*(FA(5)*FB(-4)+FA(-4)*FB(5))&		!b*cbar+cbar*b
	)*Q2*Q2/((Q2-MW2)**2+GammaW2*MW2)
!--------------------------------------------------------------------------------  
  CASE (9) !pp-> W+ + W-
	 FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	 FB=uTMDPDF_5(x2,b,mu,zeta2,1)
	 
	Integrand=paramW_L*(&
	paramW_UD*(FA(2)*FB(-1)+FA(1)*FB(-2)+FA(-2)*FB(1)+FA(-1)*FB(2))&	!u*dbar+d*ubar+ubar*d+dbar*u
	+paramW_US*(FA(2)*FB(-3)+FA(3)*FB(-2)+FA(-2)*FB(3)+FA(-3)*FB(2))&	!u*sbar+s*ubar+ubar*s+sbar*u
	+paramW_UB*(FA(2)*FB(-5)+FA(5)*FB(-2)+FA(-2)*FB(5)+FA(-5)*FB(2))&	!u*bbar+b*ubar+ubar*b+bbar*u
	+paramW_CD*(FA(4)*FB(-1)+FA(1)*FB(-4)+FA(-4)*FB(1)+FA(-1)*FB(4))&	!c*dbar+d*cbar+cbar*d+dbar*c
	+paramW_CS*(FA(4)*FB(-3)+FA(3)*FB(-4)+FA(-4)*FB(3)+FA(-3)*FB(4))&	!c*sbar+s*cbar+cbar*s+sbar*c
	+paramW_CB*(FA(4)*FB(-5)+FA(5)*FB(-4)+FA(-4)*FB(5)+FA(-5)*FB(4))&	!c*bbar+b*cbar+cbar*b+bbar*c
	)*Q2*Q2/((Q2-MW2)**2+GammaW2*MW2)
!--------------------------------------------------------------------------------  
  CASE (10) !ppbar-> W+
	 FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	 FB=uTMDPDF_5(x2,b,mu,zeta2,1)
	 FB=FB(5:-5:-1) !! inverse the quark order
	 
	Integrand=paramW_L*(&
	paramW_UD*(FA(2)*FB(-1)+FA(-1)*FB(2))&		!u*dbar+dbar*u
	+paramW_US*(FA(2)*FB(-3)+FA(-3)*FB(2))&		!u*sbar+sbar*u
	+paramW_UB*(FA(2)*FB(-5)+FA(-5)*FB(2))&		!u*bbar+bbar*u
	+paramW_CD*(FA(4)*FB(-1)+FA(-1)*FB(4))&		!c*dbar+dbar*c
	+paramW_CS*(FA(4)*FB(-3)+FA(-3)*FB(4))&		!c*sbar+sbar*c
	+paramW_CB*(FA(4)*FB(-5)+FA(-5)*FB(4))&		!c*bbar+bbar*c
	)*Q2*Q2/((Q2-MW2)**2+GammaW2*MW2)
!--------------------------------------------------------------------------------  
  CASE (11) !ppbar-> W-
	 FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	 FB=uTMDPDF_5(x2,b,mu,zeta2,1)
	 FB=FB(5:-5:-1) !! inverse the quark order
	 
	Integrand=paramW_L*(&
	paramW_UD*(FA(1)*FB(-2)+FA(-2)*FB(1))&		!d*ubar+ubar*d
	+paramW_US*(FA(3)*FB(-2)+FA(-2)*FB(3))&		!s*ubar+ubar*s
	+paramW_UB*(FA(5)*FB(-2)+FA(-2)*FB(5))&		!b*ubar+ubar*b
	+paramW_CD*(FA(1)*FB(-4)+FA(-4)*FB(1))&		!d*cbar+cbar*d
	+paramW_CS*(FA(3)*FB(-4)+FA(-4)*FB(3))&		!s*cbar+cbar*s
	+paramW_CB*(FA(5)*FB(-4)+FA(-4)*FB(5))&		!b*cbar+cbar*b
	)*Q2*Q2/((Q2-MW2)**2+GammaW2*MW2)
!--------------------------------------------------------------------------------  
  CASE (12) !ppbar-> W+ + W-
	 FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	 FB=uTMDPDF_5(x2,b,mu,zeta2,1)
	 FB=FB(5:-5:-1) !! inverse the quark order
	 
	Integrand=paramW_L*(&
	paramW_UD*(FA(2)*FB(-1)+FA(1)*FB(-2)+FA(-2)*FB(1)+FA(-1)*FB(2))&	!u*dbar+d*ubar+ubar*d+dbar*u
	+paramW_US*(FA(2)*FB(-3)+FA(3)*FB(-2)+FA(-2)*FB(3)+FA(-3)*FB(2))&	!u*sbar+s*ubar+ubar*s+sbar*u
	+paramW_UB*(FA(2)*FB(-5)+FA(5)*FB(-2)+FA(-2)*FB(5)+FA(-5)*FB(2))&	!u*bbar+b*ubar+ubar*b+bbar*u
	+paramW_CD*(FA(4)*FB(-1)+FA(1)*FB(-4)+FA(-4)*FB(1)+FA(-1)*FB(4))&	!c*dbar+d*cbar+cbar*d+dbar*c
	+paramW_CS*(FA(4)*FB(-3)+FA(3)*FB(-4)+FA(-4)*FB(3)+FA(-3)*FB(4))&	!c*sbar+s*cbar+cbar*s+sbar*c
	+paramW_CB*(FA(4)*FB(-5)+FA(5)*FB(-4)+FA(-4)*FB(5)+FA(-5)*FB(4))&	!c*bbar+b*cbar+cbar*b+bbar*c
	)*Q2*Q2/((Q2-MW2)**2+GammaW2*MW2)
	
!--------------------------------------------------------------------------------  
  CASE (13) !pp-> W+
	 FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	 FB=uTMDPDF_5(x2,b,mu,zeta2,1)
	 
	Integrand=&
	paramW_UD*(FA(2)*FB(-1)+FA(-1)*FB(2))&		!u*dbar+dbar*u
	+paramW_US*(FA(2)*FB(-3)+FA(-3)*FB(2))&		!u*sbar+sbar*u
	+paramW_UB*(FA(2)*FB(-5)+FA(-5)*FB(2))&		!u*bbar+bbar*u
	+paramW_CD*(FA(4)*FB(-1)+FA(-1)*FB(4))&		!c*dbar+dbar*c
	+paramW_CS*(FA(4)*FB(-3)+FA(-3)*FB(4))&		!c*sbar+sbar*c
	+paramW_CB*(FA(4)*FB(-5)+FA(-5)*FB(4))		!c*bbar+bbar*c
	
!--------------------------------------------------------------------------------  
  CASE (14) !pp-> W-
	 FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	 FB=uTMDPDF_5(x2,b,mu,zeta2,1)
	 
	Integrand=&
	paramW_UD*(FA(1)*FB(-2)+FA(-2)*FB(1))&		!d*ubar+ubar*d
	+paramW_US*(FA(3)*FB(-2)+FA(-2)*FB(3))&		!s*ubar+ubar*s
	+paramW_UB*(FA(5)*FB(-2)+FA(-2)*FB(5))&		!b*ubar+ubar*b
	+paramW_CD*(FA(1)*FB(-4)+FA(-4)*FB(1))&		!d*cbar+cbar*d
	+paramW_CS*(FA(3)*FB(-4)+FA(-4)*FB(3))&		!s*cbar+cbar*s
	+paramW_CB*(FA(5)*FB(-4)+FA(-4)*FB(5))		!b*cbar+cbar*b
	
!--------------------------------------------------------------------------------  
  CASE (15) !pp-> W+ + W-
	 FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	 FB=uTMDPDF_5(x2,b,mu,zeta2,1)
	 
	Integrand=&
	paramW_UD*(FA(2)*FB(-1)+FA(1)*FB(-2)+FA(-2)*FB(1)+FA(-1)*FB(2))&	!u*dbar+d*ubar+ubar*d+dbar*u
	+paramW_US*(FA(2)*FB(-3)+FA(3)*FB(-2)+FA(-2)*FB(3)+FA(-3)*FB(2))&	!u*sbar+s*ubar+ubar*s+sbar*u
	+paramW_UB*(FA(2)*FB(-5)+FA(5)*FB(-2)+FA(-2)*FB(5)+FA(-5)*FB(2))&	!u*bbar+b*ubar+ubar*b+bbar*u
	+paramW_CD*(FA(4)*FB(-1)+FA(1)*FB(-4)+FA(-4)*FB(1)+FA(-1)*FB(4))&	!c*dbar+d*cbar+cbar*d+dbar*c
	+paramW_CS*(FA(4)*FB(-3)+FA(3)*FB(-4)+FA(-4)*FB(3)+FA(-3)*FB(4))&	!c*sbar+s*cbar+cbar*s+sbar*c
	+paramW_CB*(FA(4)*FB(-5)+FA(5)*FB(-4)+FA(-4)*FB(5)+FA(-5)*FB(4))	!c*bbar+b*cbar+cbar*b+bbar*c
	
!--------------------------------------------------------------------------------  
  CASE (16) !ppbar-> W+
	 FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	 FB=uTMDPDF_5(x2,b,mu,zeta2,1)
	 FB=FB(5:-5:-1) !! inverse the quark order
	 
	Integrand=&
	paramW_UD*(FA(2)*FB(-1)+FA(-1)*FB(2))&		!u*dbar+dbar*u
	+paramW_US*(FA(2)*FB(-3)+FA(-3)*FB(2))&		!u*sbar+sbar*u
	+paramW_UB*(FA(2)*FB(-5)+FA(-5)*FB(2))&		!u*bbar+bbar*u
	+paramW_CD*(FA(4)*FB(-1)+FA(-1)*FB(4))&		!c*dbar+dbar*c
	+paramW_CS*(FA(4)*FB(-3)+FA(-3)*FB(4))&		!c*sbar+sbar*c
	+paramW_CB*(FA(4)*FB(-5)+FA(-5)*FB(4))		!c*bbar+bbar*c
	
!--------------------------------------------------------------------------------  
  CASE (17) !ppbar-> W-
	 FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	 FB=uTMDPDF_5(x2,b,mu,zeta2,1)
	 FB=FB(5:-5:-1) !! inverse the quark order
	 
	Integrand=&
	paramW_UD*(FA(1)*FB(-2)+FA(-2)*FB(1))&		!d*ubar+ubar*d
	+paramW_US*(FA(3)*FB(-2)+FA(-2)*FB(3))&		!s*ubar+ubar*s
	+paramW_UB*(FA(5)*FB(-2)+FA(-2)*FB(5))&		!b*ubar+ubar*b
	+paramW_CD*(FA(1)*FB(-4)+FA(-4)*FB(1))&		!d*cbar+cbar*d
	+paramW_CS*(FA(3)*FB(-4)+FA(-4)*FB(3))&		!s*cbar+cbar*s
	+paramW_CB*(FA(5)*FB(-4)+FA(-4)*FB(5))		!b*cbar+cbar*b
	
!--------------------------------------------------------------------------------  
  CASE (18) !ppbar-> W+ + W-
	 FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	 FB=uTMDPDF_5(x2,b,mu,zeta2,1)
	 FB=FB(5:-5:-1) !! inverse the quark order
	 
	Integrand=&
	paramW_UD*(FA(2)*FB(-1)+FA(1)*FB(-2)+FA(-2)*FB(1)+FA(-1)*FB(2))&	!u*dbar+d*ubar+ubar*d+dbar*u
	+paramW_US*(FA(2)*FB(-3)+FA(3)*FB(-2)+FA(-2)*FB(3)+FA(-3)*FB(2))&	!u*sbar+s*ubar+ubar*s+sbar*u
	+paramW_UB*(FA(2)*FB(-5)+FA(5)*FB(-2)+FA(-2)*FB(5)+FA(-5)*FB(2))&	!u*bbar+b*ubar+ubar*b+bbar*u
	+paramW_CD*(FA(4)*FB(-1)+FA(1)*FB(-4)+FA(-4)*FB(1)+FA(-1)*FB(4))&	!c*dbar+d*cbar+cbar*d+dbar*c
	+paramW_CS*(FA(4)*FB(-3)+FA(3)*FB(-4)+FA(-4)*FB(3)+FA(-3)*FB(4))&	!c*sbar+s*cbar+cbar*s+sbar*c
	+paramW_CB*(FA(4)*FB(-5)+FA(5)*FB(-4)+FA(-4)*FB(5)+FA(-5)*FB(4))	!c*bbar+b*cbar+cbar*b+bbar*c
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
  !--------------------------------------------------------------------------------  
  CASE (1003) !pbar+W->gamma* !!this is for E537
	!Wolfram has A=183,	Z=74,	N=109
	FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	FB=uTMDPDF_5(x2,b,mu,zeta2,1)
	Integrand=296d0/1647d0*(FA(2)*FB(2)+FA(-2)*FB(-2))+436d0/1647d0*(FA(2)*FB(1)+FA(-2)*FB(-1))&
	      +109d0/1647d0*(FA(1)*FB(2)+FA(-1)*FB(-2))+74d0/1647d0*(FA(1)*FB(1)+FA(-1)*FB(-1))&
	      +1d0/9d0*(FA(3)*FB(3)+FA(-3)*FB(-3)+4d0*FA(4)*FB(4)+4d0*FA(-4)*FB(-4)+FA(5)*FB(5)+FA(-5)*FB(-5))
  !----------------------------------------------------------------------------------
  !-------------------------SIDIS----------------------------------------------------
  !----------------------------------------------------------------------------------
    CASE (2001:2009) !p->hN where n=last number
	! e_q^2 *F_q(A)*F_q(B)
	h=process-2000
	FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	FB=uTMDFF_5(x2,b,mu,zeta2,h)
	Integrand=FA(1)*FB(1)/9.d0&
	  +FA(2)*FB(2)*4.d0/9.d0&
	  +FA(3)*FB(3)/9.d0&
	  +FA(4)*FB(4)*4d0/9.d0&
	  +FA(5)*FB(5)/9d0&
	  +FA(-1)*FB(-1)/9.d0&
	  +FA(-2)*FB(-2)*4.d0/9.d0&
	  +FA(-3)*FB(-3)/9.d0&
	  +FA(-4)*FB(-4)*4d0/9.d0&
	  +FA(-5)*FB(-5)/9d0
  CASE DEFAULT
    write(*,*) 'ERROR:arTeMiDe.TMDF: undefined process: ',process
    write(*,*) 'Evaluation stop'
    stop
 END SELECT
 
  if(ISNAN(Integrand)) then
   write(*,*) 'arTeMiDe TMDF: CRITICAL ERROR. Integrand evaluated to NaN'
   write(*,*) 'bT=',b, 'x1,x2=',x1,x2,' process=',process
   write(*,*) 'mu=',mu, 'Q2=',Q2
   !write(*,*) 'Current set of NP parameters ------------'
   !write(*,*) currentNP
   write(*,*) 'arTeMiDe: ConvergenceLOST trigger ON'
   call TMDF_convergenceISlost()
   Integrand=1d10
   end if
  
   if(Integrand>1d32) then
   write(*,*) 'arTeMiDe TMDF: CRITICAL ERROR. Integrand evaluated to >10^32'
   write(*,*) 'bT=',b, 'x1,x2=',x1,x2,' process=',process
   write(*,*) 'mu=',mu, 'Q2=',Q2
   !write(*,*) 'Current set of NP parameters ------------'
   !write(*,*) currentNP
   write(*,*) 'arTeMiDe: convergenceLOST trigger ON'
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
     
end function XIntegrandForDYwithZgamma
end module TMDF