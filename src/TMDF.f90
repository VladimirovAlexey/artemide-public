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
  use aTMDe_Numerics
  use IO_functions
  use TMDs
  use EWinput
  implicit none

  private
!   public
 
 character (len=7),parameter :: moduleName="TMDF"
 character (len=5),parameter :: version="v2.03"
 !Last appropriate verion of constants-file
 integer,parameter::inputver=1

!------------------------------------------Tables-----------------------------------------------------------------------
    integer,parameter::Nmax=1000
    INCLUDE 'Tables/BesselZero1000.f90'
!------------------------------------------Working variables------------------------------------------------------------
  
  integer::outputLevel=2
  !! variable that count number of WRNING mesagges. In order not to spam too much
  integer::messageTrigger=6
  logical::started=.false.
  
  logical:: convergenceLost=.false.
  
  !!!!! I split the qT over runs qT<qTSegmentationBoundary
  !!!!! In each segment I have the ogata quadrature with h=hOGATA*hSegmentationWeight
  !!!!! It helps to convergen integrals, since h(optimal) ~ qT
  integer,parameter::hSegmentationNumber=5
  real(dp),dimension(1:hSegmentationNumber),parameter::hSegmentationWeight=(/0.001d0,0.01d0,0.1d0,1d0,5d0/)
  real(dp),dimension(1:hSegmentationNumber),parameter::qTSegmentationBoundary=(/0.001d0,0.01d0,0.1d0,1d0,50d0/)
  
  real(dp)::hOGATA,tolerance
  !!!weights of ogata quadrature
  real(dp),dimension(1:hSegmentationNumber,0:3,1:Nmax)::ww,ww0
  !!!nodes of ogata quadrature
  real(dp),dimension(1:hSegmentationNumber,0:3,1:Nmax)::bb,bb0
  
  
  
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

   !! Initialization of the package
  subroutine TMDF_Initialize(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=300)::path,line
    logical::initRequared
    integer::FILEver
    
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
      if(outputLevel>0) write(*,*) color('----- arTeMiDe.TMDF '//trim(version)//': .... initialized',c_green)
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
  if(outputLevel>1) write(*,*) WarningString('convergenceLOST trigger ON',moduleName)
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
  integer::i,k,j
  real(dp)::hS!=h*hSegmentationWeight
  real(dp)::xi,qqq
   
  do j=1,hSegmentationNumber
  do k=0,3
  do i=1,Nmax
    
    hS=hOGATA*hSegmentationWeight(j)    
    xi=JZero(k,i)
    
!     ww(j,k,i)=BESSEL_JN(k,bb(j,k,i))/xi/(BESSEL_JN(k+1,xi)**2)&
! 	    *(pi*xi*hS*Cosh(xi*hS)+Sinh(pi*Sinh(xi*hS)))/(1d0+Cosh(pi*Sinh(xi*hS)))
    
    !!! if we too far away in xI*hS, the double exponential grow rapidly.
    !!! and for >6, it generates term 10^{300} and exceed the presision

    if(xi*hS>6.d0) then
        bb(j,k,i)=xi*Tanh(piHalf*Sinh(xi*hS))
        ww(j,k,i)=BESSEL_JN(k,bb(j,k,i))/xi/(BESSEL_JN(k+1,xi)**2)
        
    else
        bb(j,k,i)=xi*Tanh(piHalf*Sinh(xi*hS))
        ww(j,k,i)=BESSEL_JN(k,bb(j,k,i))/xi/(BESSEL_JN(k+1,xi)**2)&
        *(pi*xi*hS*Cosh(xi*hS)/(2d0*Cosh(piHalf*Sinh(xi*hS))**2)+Tanh(piHalf*Sinh(xi*hS)))
    end if

  end do
  end do
  end do 
 end subroutine PrepareTables
 
 !!!This is the defining module function
 !!! It evaluates the integral 
 !!!  int_0^infty   b db/2  Jn(b qT) zff F1 F2
 !!!
 function TMDF_F(Q2,qT,x1,x2,mu,zeta1,zeta2,process)
  real(dp)::TMDF_F
  real(dp)::qT,x1,x2,mu,zeta1,zeta2,Q2
  integer::process
  real(dp)::integral,eps,delta
  real(dp)::v1,v2,v3,v4
  integer::k,n,j,Nsegment
    
  CallCounter=CallCounter+1
  integral=0d0
  
  if(qT<0.0000001d0 .or. x1>=1d0 .or. x2>=1d0) then  
  integral=0d0
  else if(TMDF_IsconvergenceLost()) then	
  !!!in the case of lost convergence we return huge number (divergent xSec)
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
  
  !!! define segment of qT
  do j=1,hSegmentationNumber
    if(qT<qTSegmentationBoundary(j)) exit
  end do
  if(j>hSegmentationNumber) then
    Nsegment=hSegmentationNumber
  else
    Nsegment=j
  end if
  
  !!! sum over OGATA nodes
  do k=1,Nmax!!! maximum of number of bessel roots preevaluated in the head
    eps=ww(Nsegment,n,k)*(bb(Nsegment,n,k)**(n+1))*Integrand(Q2,bb(Nsegment,n,k)/qT,x1,x2,mu,zeta1,zeta2,process)
    
    v4=v3
    v3=v2
    v2=v1
    v1=ABS(eps)
    
    delta=(v1+v2+v3+v4)
    integral=integral+eps
    
    !!! here we check that residual term is smaller than already collected integral
    !!! also checking the zerothness of the integral. If already collected integral is null it is null
    !!! Here is potential bug. If the first 10 points give zero (whereas some later points do not), the integral will be zero
    if((delta<tolerance*abs(integral) .or. abs(integral)<1d-32) .and. k>=10) exit
  end do
  if(k>=Nmax) then	
    if(outputlevel>0) WRITE(*,*) WarningString('OGATA quadrature diverge. TMD decaing too slow? ',moduleName)
      if(outputlevel>1) then
      write(*,*) 'Information over the last call ----------'
      write(*,*) 'bt/qT= ',bb(Nsegment,n,Nmax)/qT, 'qT=',qT, '| segmentation zone=',Nsegment,&
	      ' ogata h=',hOGATA*hSegmentationWeight(Nsegment)
      write(*,*) 'W=',Integrand(Q2,bb(Nsegment,n,Nmax)/qT,x1,x2,mu,zeta1,zeta2,process), 'eps/integral =', eps/integral
      write(*,*) 'residual term=',delta, '>',tolerance
      write(*,*) '(x1,x2)=(',x1,',',x2,')'
      write(*,*) 'process =',process,' it is ',CallCounter,'call.'
      write(*,*) '------------------------------------------'
      end if
    call TMDF_convergenceISlost()
  end if
  
  if(k>MaxCounter) MaxCounter=k-1
!   write(*,*) 'Integral=',integral
  TMDF_F=integral/(qT**(n+2)) 
  end if
  !write(*,*) 'Last call: ',k
  
!    write(*,'("{",F6.2,",",F18.16,"},")') qT,x1*x2*TMDF_F
 end function TMDF_F
 
 function Integrand(Q2,b,x1,x2,mu,zeta1,zeta2,process)
 real(dp)::Integrand
 real(dp)::b,x1,x2,mu,zeta1,zeta2,Q2
 integer::process,h
 real(dp),dimension(-5:5)::FA,FB,FAB
 
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
  CASE(9999,19999,29999,39999)
    Integrand=Exp(-mu*b)*(1d0+x1*b**2+x2*b**4)
  CASE(9998,19998,29998,39998)
    Integrand=Exp(-mu*b**2)*(1d0+x1*b**2+x2*b**4)
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
  CASE(20) !pp -> Higgs (unpol.part+lin.pol.part)
	FA=uTMDPDF_50(x1,b,mu,zeta1,1)
	FB=uTMDPDF_50(x2,b,mu,zeta2,1)
	Integrand=FA(0)*FB(0) !!!! unpolarized part
	
	FA=lpTMDPDF_50(x1,b,mu,zeta1,1)
	FB=lpTMDPDF_50(x2,b,mu,zeta2,1)
	Integrand=Integrand+FA(0)*FB(0) !!!! linearly polarized part
!--------------------------------------------------------------------------------  
  CASE(21) !pp -> Higgs (unpol.part)    
	FA=uTMDPDF_50(x1,b,mu,zeta1,1)
	FB=uTMDPDF_50(x2,b,mu,zeta2,1)
	Integrand=FA(0)*FB(0)
  
!--------------------------------------------------------------------------------  
  CASE(22) !pp -> Higgs (lin.pol.part)
	FA=lpTMDPDF_50(x1,b,mu,zeta1,1)
	FB=lpTMDPDF_50(x2,b,mu,zeta2,1)
	Integrand=FA(0)*FB(0)
	
!--------------------------------------------------------------------------------
  CASE (101) !p h->gamma
	! e_q^2 *F_q(A)*F_qbar(B)
	if(zeta1==zeta2) then
	 FAB=uPDF_uPDF(x1,x2,b,mu,zeta1,1,2)
	else
	 FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	 FB=uTMDPDF_5(x2,b,mu,zeta2,2)
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
  CASE (102) !pbar h->gamma
	! e_q^2 *F_q(A)*F_q(B)
	if(zeta1==zeta2) then
	 FAB=uPDF_anti_uPDF(x1,x2,b,mu,zeta1,1,2)
	else
	 FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	 FB=uTMDPDF_5(x2,b,mu,zeta2,2)
	 FAB=FA*FB
	end if
	!! in fact, we must revert this array, but the coefficients are symmetric
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
  CASE (103) !p hbar->gamma
	! e_q^2 *F_q(A)*F_qbar(B)
	if(zeta1==zeta2) then
	 FAB=uPDF_anti_uPDF(x1,x2,b,mu,zeta1,1,2)
	else
	 FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	 FB=uTMDPDF_5(x2,b,mu,zeta2,2)
	 FAB=FA*FB
	end if
	!! in fact, we must revert this array, but the coefficients are symmetric
	
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
  CASE (104) !pbar hbar->gamma
	! e_q^2 *F_q(A)*F_q(B)
	! e_q^2 *F_q(A)*F_qbar(B)
	if(zeta1==zeta2) then
	 FAB=uPDF_uPDF(x1,x2,b,mu,zeta1,1,2)
	else
	 FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	 FB=uTMDPDF_5(x2,b,mu,zeta2,2)
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
  !--------------------------------------------------------------------------------  
  CASE (1004) !pminus+W->gamma* !!this is for E537
	!Wolfram has A=183,	Z=74,	N=109
	FA=uTMDPDF_5(x1,b,mu,zeta1,2)
	FB=uTMDPDF_5(x2,b,mu,zeta2,1)
	Integrand=296d0/1647d0*(FA(-2)*FB(2)+FA(2)*FB(-2))+436d0/1647d0*(FA(-2)*FB(1)+FA(2)*FB(-1))&
	      +109d0/1647d0*(FA(-1)*FB(2)+FA(1)*FB(-2))+74d0/1647d0*(FA(-1)*FB(1)+FA(1)*FB(-1))&
	      +1d0/9d0*(FA(-3)*FB(3)+FA(3)*FB(-3)+4d0*FA(-4)*FB(4)+4d0*FA(4)*FB(-4)+FA(-5)*FB(5)+FA(5)*FB(-5))
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
  !--------------------------------------------------------------------------------  
  CASE (2011:2019) !d->hN where n=last number (d=deutron=(p+n)/2)
	! e_q^2 *F_q(A)*F_q(B)
	h=process-2010
	FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	FB=uTMDFF_5(x2,b,mu,zeta2,h)
	Integrand=(FA(1)+FA(2))*(FB(1)+4d0*FB(2))/18d0&
	  +FA(3)*FB(3)/9.d0&
	  +FA(4)*FB(4)*4d0/9.d0&
	  +FA(5)*FB(5)/9d0&
	  +(FA(-1)+FA(-2))*(FB(-1)+4d0*FB(-2))/18d0&
	  +FA(-3)*FB(-3)/9.d0&
	  +FA(-4)*FB(-4)*4d0/9.d0&
	  +FA(-5)*FB(-5)/9d0
  !--------------------------------------------------------------------------------  
  CASE (2021:2029) !p->bar-hN where n=last number
	! e_q^2 *F_q(A)*F_bar-q(B)
	h=process-2020
	FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	FB=uTMDFF_5(x2,b,mu,zeta2,h)
	Integrand=FA(1)*FB(-1)/9.d0&
	  +FA(2)*FB(-2)*4.d0/9.d0&
	  +FA(3)*FB(-3)/9.d0&
	  +FA(4)*FB(-4)*4d0/9.d0&
	  +FA(5)*FB(-5)/9d0&
	  +FA(-1)*FB(1)/9.d0&
	  +FA(-2)*FB(2)*4.d0/9.d0&
	  +FA(-3)*FB(3)/9.d0&
	  +FA(-4)*FB(4)*4d0/9.d0&
	  +FA(-5)*FB(5)/9d0
!--------------------------------------------------------------------------------  
  CASE (2031:2039) !d->bar-hN where n=last number (d=deutron=(p+n)/2)
	! e_q^2 *F_q(A)*F_bar-q(B)
	h=process-2030
	FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	FB=uTMDFF_5(x2,b,mu,zeta2,h)
	Integrand=(FA(1)+FA(2))*(FB(-1)+4d0*FB(-2))/18d0&
	  +FA(3)*FB(-3)/9.d0&
	  +FA(4)*FB(-4)*4d0/9.d0&
	  +FA(5)*FB(-5)/9d0&
	  +(FA(-1)+FA(-2))*(FB(1)+4d0*FB(2))/18d0&
	  +FA(-3)*FB(3)/9.d0&
	  +FA(-4)*FB(4)*4d0/9.d0&
	  +FA(-5)*FB(5)/9d0
!--------------------------------------------------------------------------------  
!--------------------------------------------------------------------------------  
   CASE (2101) !p->h? where h?=h1+h2
	! e_q^2 *F_q(A)*F_q(B)
	FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	FB=uTMDFF_5(x2,b,mu,zeta2,1)+uTMDFF_5(x2,b,mu,zeta2,2)
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
!--------------------------------------------------------------------------------  
    CASE (2102) !p->h? where h?=h1+h2+h3
	! e_q^2 *F_q(A)*F_q(B)
	FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	FB=uTMDFF_5(x2,b,mu,zeta2,1)+uTMDFF_5(x2,b,mu,zeta2,2)+uTMDFF_5(x2,b,mu,zeta2,3)
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
!--------------------------------------------------------------------------------  
    CASE (2103) !d->h? where h?=h1+h2 (d=deutron=(p+n)/2)
	! e_q^2 *F_q(A)*F_q(B)
	FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	FB=uTMDFF_5(x2,b,mu,zeta2,1)+uTMDFF_5(x2,b,mu,zeta2,2)
	Integrand=(FA(1)+FA(2))*(FB(1)+4d0*FB(2))/18d0&
	  +FA(3)*FB(3)/9.d0&
	  +FA(4)*FB(4)*4d0/9.d0&
	  +FA(5)*FB(5)/9d0&
	  +(FA(-1)+FA(-2))*(FB(-1)+4d0*FB(-2))/18d0&
	  +FA(-3)*FB(-3)/9.d0&
	  +FA(-4)*FB(-4)*4d0/9.d0&
	  +FA(-5)*FB(-5)/9d0
!--------------------------------------------------------------------------------  
    CASE (2104) !d->h? where h?=h1+h2+h3 (d=deutron=(p+n)/2)
	! e_q^2 *F_q(A)*F_q(B)
	FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	FB=uTMDFF_5(x2,b,mu,zeta2,1)+uTMDFF_5(x2,b,mu,zeta2,2)+uTMDFF_5(x2,b,mu,zeta2,3)
	Integrand=(FA(1)+FA(2))*(FB(1)+4d0*FB(2))/18d0&
	  +FA(3)*FB(3)/9.d0&
	  +FA(4)*FB(4)*4d0/9.d0&
	  +FA(5)*FB(5)/9d0&
	  +(FA(-1)+FA(-2))*(FB(-1)+4d0*FB(-2))/18d0&
	  +FA(-3)*FB(-3)/9.d0&
	  +FA(-4)*FB(-4)*4d0/9.d0&
	  +FA(-5)*FB(-5)/9d0
!------------------------------------------------------------------------------------
  CASE (2111) !p->bar h? where h?=h1+h2
	! e_q^2 *F_q(A)*F_bq(B)
	FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	FB=uTMDFF_5(x2,b,mu,zeta2,1)+uTMDFF_5(x2,b,mu,zeta2,2)
	Integrand=FA(1)*FB(-1)/9.d0&
	  +FA(2)*FB(-2)*4.d0/9.d0&
	  +FA(3)*FB(-3)/9.d0&
	  +FA(4)*FB(-4)*4d0/9.d0&
	  +FA(5)*FB(-5)/9d0&
	  +FA(-1)*FB(1)/9.d0&
	  +FA(-2)*FB(2)*4.d0/9.d0&
	  +FA(-3)*FB(3)/9.d0&
	  +FA(-4)*FB(4)*4d0/9.d0&
	  +FA(-5)*FB(5)/9d0
!--------------------------------------------------------------------------------  
    CASE (2112) !p->bar h? where h?=h1+h2+h3
	! e_q^2 *F_q(A)*F_bq(B)
	FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	FB=uTMDFF_5(x2,b,mu,zeta2,1)+uTMDFF_5(x2,b,mu,zeta2,2)+uTMDFF_5(x2,b,mu,zeta2,3)
	Integrand=FA(1)*FB(-1)/9.d0&
	  +FA(2)*FB(-2)*4.d0/9.d0&
	  +FA(3)*FB(-3)/9.d0&
	  +FA(4)*FB(-4)*4d0/9.d0&
	  +FA(5)*FB(-5)/9d0&
	  +FA(-1)*FB(1)/9.d0&
	  +FA(-2)*FB(2)*4.d0/9.d0&
	  +FA(-3)*FB(3)/9.d0&
	  +FA(-4)*FB(4)*4d0/9.d0&
	  +FA(-5)*FB(5)/9d0
!--------------------------------------------------------------------------------  
    CASE (2113) !d->bar h? where h?=h1+h2 (d=deutron=(p+n)/2)
	! e_q^2 *F_q(A)*F_bq(B)
	FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	FB=uTMDFF_5(x2,b,mu,zeta2,1)+uTMDFF_5(x2,b,mu,zeta2,2)
	Integrand=(FA(1)+FA(2))*(FB(-1)+4d0*FB(-2))/18d0&
	  +FA(3)*FB(-3)/9.d0&
	  +FA(4)*FB(-4)*4d0/9.d0&
	  +FA(5)*FB(-5)/9d0&
	  +(FA(-1)+FA(-2))*(FB(1)+4d0*FB(2))/18d0&
	  +FA(-3)*FB(3)/9.d0&
	  +FA(-4)*FB(4)*4d0/9.d0&
	  +FA(-5)*FB(5)/9d0
!--------------------------------------------------------------------------------  
    CASE (2114) !d->bar h? where h?=h1+h2+h3 (d=deutron=(p+n)/2)
	! e_q^2 *F_q(A)*F_bq(B)
	FA=uTMDPDF_5(x1,b,mu,zeta1,1)
	FB=uTMDFF_5(x2,b,mu,zeta2,1)+uTMDFF_5(x2,b,mu,zeta2,2)+uTMDFF_5(x2,b,mu,zeta2,3)
	Integrand=(FA(1)+FA(2))*(FB(-1)+4d0*FB(-2))/18d0&
	  +FA(3)*FB(-3)/9.d0&
	  +FA(4)*FB(-4)*4d0/9.d0&
	  +FA(5)*FB(-5)/9d0&
	  +(FA(-1)+FA(-2))*(FB(1)+4d0*FB(2))/18d0&
	  +FA(-3)*FB(3)/9.d0&
	  +FA(-4)*FB(4)*4d0/9.d0&
	  +FA(-5)*FB(5)/9d0
  CASE DEFAULT
    write(*,*) ErrorString('undefined process: ',moduleName),process
    write(*,*) color('Evaluation stop',c_red_bold)
    stop
 END SELECT
 
  if(ISNAN(Integrand)) then
   write(*,*) ErrorString('Integrand evaluated to NaN',moduleName)
   write(*,*) 'bT=',b, 'x1,x2=',x1,x2,' process=',process
   write(*,*) 'mu=',mu, 'Q2=',Q2
   !write(*,*) 'Current set of NP parameters ------------'
   !write(*,*) currentNP
   call TMDF_convergenceISlost()
   Integrand=1d10
   end if
  
   if(Integrand>1d32) then
   write(*,*) ErrorString('Integrand evaluated to >10^32',moduleName)
   write(*,*) 'bT=',b, 'x1,x2=',x1,x2,' process=',process
   write(*,*) 'mu=',mu, 'Q2=',Q2
   !write(*,*) 'Current set of NP parameters ------------'
   !write(*,*) currentNP
   call TMDF_convergenceISlost()
   Integrand=1d10
   end if
 
 end function Integrand

!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!!! The hadron tensonr for the DY icludes Z + gamma, evaluated at FA and FB 
function XIntegrandForDYwithZgamma(FAB,Q2)
     real(dp)::XIntegrandForDYwithZgamma,Q2
    !!cross-seciton parameters
     real(dp),dimension(-5:5):: FAB
     
     !!!parameters of Z boson coupling
!      real(dp),parameter:: paramU=0.40329064872689d0 !! ((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) for eq=2/3
!      real(dp),parameter:: paramD=0.51983027428079d0 !! ((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) for eq=1/3
!      real(dp),parameter:: paramS=0.51983027428079d0  !! ((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) for eq=1/3
!      real(dp),parameter:: paramC=0.40329064872689d0 !! ((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) for eq=2/3
!      real(dp),parameter:: paramB=0.51983027428079d0  !! ((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) for eq=1/3
!      real(dp),parameter:: paramL=0.35358707798999d0  !! ((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) for eq=1
     
     !!!parameters of Z-gamma boson coupling
!      real(dp),parameter:: paramMIXU=0.1515661518957d0 !! e(T3-2e sW^2)/2sWcW for eq=+2/3,, T3=+1/2
!      real(dp),parameter:: paramMIXD=0.1367184036034d0 !! e(T3-2e sW^2)/2sWcW for eq=-1/3,T3=-1/2
!      real(dp),parameter:: paramMIXS=0.1367184036034d0  !! e(T3-2e sW^2)/2sWcW for eq=-1/3, T3=-1/2
!      real(dp),parameter:: paramMIXC=0.1515661518957d0 !! e(T3-2e sW^2)/2sWcW for eq=+2/3,, T3=+1/2
!      real(dp),parameter:: paramMIXB=0.1367184036034d0  !! e(T3-2e sW^2)/2sWcW for eq=-1/3, T3=-1/2
!      real(dp),parameter:: paramMIXL=0.0445432448766d0  !! e(T3-2e sW^2)/2sWcW for eq=-1, T3=-1/2
     
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
