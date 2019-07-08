!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 1.4
!
!	Evaluation of the TMDs in KT space
!	
!	if you use this module please, quote 1902.08474
!
!	ver 1.0: release (AV, 23.12.2018)
!	ver 2.00: release (AV, 29.03.2019)
!
!				A.Vladimirov (23.12.2018)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDs_inKT
  use TMDs
  implicit none
  
  private
!   public
 
 character (len=10),parameter :: moduleName="TMDs-inKT"
 character (len=5),parameter :: version="v2.01"
 !Last appropriate verion of constants-file
  integer,parameter::inputver=1
 
 !------------------------------------------Tables-----------------------------------------------------------------------
    integer,parameter::Nmax=200
    INCLUDE 'Tables/BesselZero.f90'
 
   logical:: convergenceLost=.false.
 
   real*8::hOGATA,tolerance
  !!!weights of ogata quadrature
  real*8,dimension(0:3,1:Nmax)::ww
  !!!nodes of ogata quadrature
  real*8,dimension(0:3,1:Nmax)::bb
  
  integer::GlobalCounter
  integer::CallCounter
  integer::MaxCounter
  integer::messageCounter
 
!------------------------------------------Physical and mathematical constants------------------------------------------
  
!------------------------------------------Working variables------------------------------------------------------------
  
  logical::started=.false.
  
  integer::outputLevel=2
  integer::messageTrigger=5
  
  public::TMDs_inKT_Initialize,TMDs_inKT_ShowStatistic,TMDs_inKT_IsInitialized,TMDs_inKT_ResetCounters
	  
  real*8,dimension(-5:5),public::uTMDPDF_kT_50,uTMDPDF_kT_5,uTMDFF_kT_5,uTMDFF_kT_50,lpTMDPDF_kT_50
  
  interface uTMDPDF_kT_5
    module procedure uTMDPDF_kT_5_Ev,uTMDPDF_kT_5_optimal
  end interface
  
  interface uTMDPDF_kT_50
    module procedure uTMDPDF_kT_50_Ev,uTMDPDF_kT_50_optimal
  end interface
  
  interface uTMDFF_kT_5
    module procedure uTMDFF_kT_5_Ev,uTMDFF_kT_5_optimal
  end interface
  
  interface uTMDFF_kT_50
    module procedure uTMDFF_kT_50_Ev,uTMDFF_kT_50_optimal
  end interface
  
  interface lpTMDPDF_kT_50
    module procedure lpTMDPDF_kT_50_Ev,lpTMDPDF_kT_50_optimal
  end interface

contains 
 function TMDs_inKT_IsInitialized()
  logical::TMDs_inKT_IsInitialized
  TMDs_inKT_IsInitialized=started
 end function TMDs_inKT_IsInitialized

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
  subroutine TMDs_inKT_Initialize(file,prefix)
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
    
    call MoveTO(51,'*8   ')
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
      
      if(outputLevel>1) write(*,*) 'arTeMiDe.TMDs-inKT: preparing Ogata tables'
      call PrepareTables()
      if(outputLevel>2) write(*,'(A,I4)') ' | Maximum number of nodes	:',Nmax
      if(outputLevel>1) write(*,*) 'arTeMiDe.TMDs-inKT: Ogata tables prepared'
      
      GlobalCounter=0
      CallCounter=0
      MaxCounter=0
      convergenceLost=.false.
      
      if(.not.TMDs_IsInitialized()) then
	if(outputLevel>1) write(*,*) '.. initializing TMDs (from ',moduleName,')'
	if(present(prefix)) then
	  call TMDs_Initialize(file,prefix)
	else
	  call TMDs_Initialize(file)
	end if
      end if
      
      started=.true.
      if(outputLevel>0) write(*,*) '----- arTeMiDe.TMDs-inKT ',version,': .... initialized'
      if(outputLevel>1) write(*,*) ' '
    
  end subroutine TMDs_inKT_Initialize


  
  subroutine TMDs_inKT_ResetCounters()
    
    convergenceLost=.false.
    GlobalCounter=0
    CallCounter=0
    MaxCounter=0
    
  end subroutine TMDs_inKT_ResetCounters
  
    !!!!!!!Functions which carry the trigger on convergences.... Its used in xSec, and probably in other places.
  function TMDs_inKT_IsconvergenceLost()
  logical::TMDs_inKT_IsconvergenceLost
  
  TMDs_inKT_IsconvergenceLost=convergenceLost
  end function TMDs_inKT_IsconvergenceLost
  
  subroutine TMDs_inKT_convergenceISlost()  
  convergenceLost=.true.
  if(outputLevel>1) write(*,*) 'arTeMiDe.TMDs-in-kT: convergence triger set to be lost.'
  end subroutine TMDs_inKT_convergenceISlost
  
  subroutine TMDs_inKT_ShowStatistic()
    if(convergenceLost) then
      write(*,*) '         TMDs-in-kT statistics: convergence has been lost.'
    else
      write(*,'(A,ES12.3)') 'TMDs-in-kT statistics         total calls of TMDs  :  ',Real(2*GlobalCounter)
      write(*,'(A,ES12.3)') '                              total calls of TMD_F :  ',Real(CallCounter)
      write(*,'(A,F12.3)')  '                                         avarage M :  ',Real(GlobalCounter)/Real(CallCounter)
      write(*,'(A,I12)')    '                                     maximum calls :  ',MaxCounter
	
    end if
  end subroutine TMDs_inKT_ShowStatistic
  
  
   !!!Prepare tables for Ogata quadrature with given h
   !!! note that the factor 1/(2pi) is taken into ww
!  subroutine PrepareTables()
!   real*8,parameter::piHalf=1.5707963267948966d0
!   real*8,parameter::pi=3.141592653589793d0
!   integer::i
!   real*8::t!=h*xi
!   real*8::psiPart!=tanh[pi/2 Sinh[h xi]]
!   
!   do i=1,Nmax
!     t=hOGATA*JZero(0,i)
!     psiPart=Tanh(piHalf*Sinh(t))
!       bb(i)=JZero(0,i)*psiPart
!       ww(i)=BESSEL_J0(bb(i))/JZero(0,i)/(BESSEL_J1(JZero(0,i))**2)*(pi*t*Cosh(t)+Sinh(pi*Sinh(t)))/(1d0+Cosh(pi*Sinh(t)))
! !      write(*,*) psiPart,b(k,i),w(k,i)
!   end do
!  
!  end subroutine PrepareTables
 
 
  !!!Prepare tables for Ogata quadrature with given h
  !!! note that the factor 1/(2pi) is taken into ww
  !!! the difference between definition in TMDF and here is 1/pi
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
      ww(k,i)=BESSEL_JN(k,bb(k,i))/pi/JZero(k,i)/(BESSEL_JN(k+1,JZero(k,i))**2)&
			*(pi*t*Cosh(t)+Sinh(pi*Sinh(t)))/(1d0+Cosh(pi*Sinh(t)))
!      write(*,*) psiPart,b(k,i),w(k,i)
  end do
  end do
  
 end subroutine PrepareTables
 
 !--------------------------------------INTERFACES TO TMD------------------------------------------------
 
 !---------------------------------------------------uTMDPDF
 function uTMDPDF_kT_5_Ev(x,qT,mu,zeta,hadron)
 real*8::uTMDPDF_kT_5_Ev(-5:5)
 real*8::x,qT,mu,zeta
 integer::hadron
 uTMDPDF_kT_5_Ev=Fourier(x,qT,mu,zeta,1,hadron) 
 end function uTMDPDF_kT_5_Ev
 
 function uTMDPDF_kT_50_Ev(x,qT,mu,zeta,hadron)
 real*8::uTMDPDF_kT_50_Ev(-5:5)
 real*8::x,qT,mu,zeta
 integer::hadron
 uTMDPDF_kT_50_Ev=Fourier(x,qT,mu,zeta,2,hadron) 
 end function uTMDPDF_kT_50_Ev
 
  function uTMDPDF_kT_5_optimal(x,qT,hadron)
 real*8::uTMDPDF_kT_5_optimal(-5:5)
 real*8::x,qT
 integer::hadron
 uTMDPDF_kT_5_optimal=Fourier(x,qT,10d0,10d0,3,hadron) 
 end function uTMDPDF_kT_5_optimal
 
 function uTMDPDF_kT_50_optimal(x,qT,hadron)
 real*8::uTMDPDF_kT_50_optimal(-5:5)
 real*8::x,qT
 integer::hadron
 uTMDPDF_kT_50_optimal=Fourier(x,qT,10d0,10d0,4,hadron) 
 end function uTMDPDF_kT_50_optimal
 
 !---------------------------------------------------uTMDFF
 function uTMDFF_kT_5_Ev(x,qT,mu,zeta,hadron)
 real*8::uTMDFF_kT_5_Ev(-5:5)
 real*8::x,qT,mu,zeta
 integer::hadron
 uTMDFF_kT_5_Ev=Fourier(x,qT,mu,zeta,5,hadron) 
 end function uTMDFF_kT_5_Ev
 
 function uTMDFF_kT_50_Ev(x,qT,mu,zeta,hadron)
 real*8::uTMDFF_kT_50_Ev(-5:5)
 real*8::x,qT,mu,zeta
 integer::hadron
 uTMDFF_kT_50_Ev=Fourier(x,qT,mu,zeta,6,hadron) 
 end function uTMDFF_kT_50_Ev
 
  function uTMDFF_kT_5_optimal(x,qT,hadron)
 real*8::uTMDFF_kT_5_optimal(-5:5)
 real*8::x,qT
 integer::hadron
 uTMDFF_kT_5_optimal=Fourier(x,qT,10d0,10d0,7,hadron) 
 end function uTMDFF_kT_5_optimal
 
 function uTMDFF_kT_50_optimal(x,qT,hadron)
 real*8::uTMDFF_kT_50_optimal(-5:5)
 real*8::x,qT
 integer::hadron
 uTMDFF_kT_50_optimal=Fourier(x,qT,10d0,10d0,8,hadron) 
 end function uTMDFF_kT_50_optimal
 
 
 !---------------------------------------------------lpTMDPDF
 
 function lpTMDPDF_kT_50_Ev(x,qT,mu,zeta,hadron)
 real*8::lpTMDPDF_kT_50_Ev(-5:5)
 real*8::x,qT,mu,zeta
 integer::hadron
 lpTMDPDF_kT_50_Ev=Fourier(x,qT,mu,zeta,9,hadron) 
 end function lpTMDPDF_kT_50_Ev
 
 function lpTMDPDF_kT_50_optimal(x,qT,hadron)
 real*8::lpTMDPDF_kT_50_optimal(-5:5)
 real*8::x,qT
 integer::hadron
 lpTMDPDF_kT_50_optimal=Fourier(x,qT,10d0,10d0,10,hadron) 
 end function lpTMDPDF_kT_50_optimal
 
 
 !------------------------------------------FOURIER--------------------------------
  !!!This is the defining module function
 !!! It evaluates the integral 
 !!!  int_0^infty   b db/2pi  J0(b qT) F1
 !!! The function F1 is given via number.. num
 function Fourier(x,qT_in,mu,zeta,num,hadron)
  real*8::qT,x,mu,zeta,qT_in
  integer::num,hadron
  real*8::integral(-5:5),eps(-5:5)
  real*8::v1,v2,v3,v4
  integer::k,n
  real*8::Fourier(-5:5)
  
  CallCounter=CallCounter+1
  integral=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
  
  if(qT_in<0.001d0) then  
  qT=0.001d0  
  else
  qT=qT_in
  end if
  
  !!!in the case of lost convergence we return huge number (divergent xSec)
  if(TMDs_inKT_IsconvergenceLost()) then	
	Fourier=integral+1d10		
  else
  
  v1=1d0
  v2=1d0
  v3=1d0
  v4=1d0
  
  !!!! select the order of Bessel function for transform
  SELECT CASE(num)
    CASE(9,10)
      n=2
    CASE DEFAULT
      n=0
  END SELECT
  
  do k=1,Nmax!!! maximum of number of bessel roots preevaluated in the head
    eps=ww(n,k)*bb(n,k)*Integrand(bb(n,k)/qT,x,mu,zeta,num,hadron)
    
    integral=integral+eps
!     if(k>8) then
      v4=v3
      v3=v2
      v2=v1
      v1=ABS(eps(0))+ABS(eps(1))+ABS(eps(2))!!! we check by u+d+g
!       if(v1+v2+v3+v4>0.7d0*ABS(integral)) write(17,*) bb(n,k)/qT,x1,x2
!       write(*,*) k, eps,ww(n,k),eps/ww(n,k)
      if(v1+v2+v3+v4<=tolerance*(ABS(integral(0))+ABS(integral(1))+ABS(integral(2)))) exit
!     end if
  end do
  if(k>=Nmax) then	
    if(outputlevel>0) WRITE(*,*) 'WARNING arTeMiDe.TMDs-in-kT: OGATA quadrature diverge. TMD decaing too slow? '
      if(outputlevel>1) then
      write(*,*) 'Information over the last call ----------'
      write(*,*) 'bt/qT= ',bb(n,Nmax)/qT, 'qT=',qT
      write(*,*) 'W=',Integrand(bb(n,Nmax)/qT,x,mu,zeta,num,hadron), 'eps/integral =', eps/integral
      write(*,*) 'v1+v2+v3+v4=',v1+v2+v3+v4, '>',tolerance*ABS(integral)
      write(*,*) 'x=',x,'type =',num,' it is ',CallCounter,'call.'
      write(*,*) '------------------------------------------'
      end if
    call TMDs_inKT_convergenceISlost()
  end if
  if(k>MaxCounter) MaxCounter=k-1
!   write(*,*) 'Integral=',integral
  Fourier=integral/(qT**2) 
  end if 
  !write(*,*) 'Last call: ',k
 end function Fourier
 
 function Integrand(b,x,mu,zeta,num,hadron)
 real*8::b,x,mu,zeta
 integer::num,hadron
 real*8::Integrand(-5:5)
 
 !increment counter 
 GlobalCounter=GlobalCounter+1
 SELECT CASE(num)
  CASE(1) !!! uTMDPDF  quarks
   Integrand=uTMDPDF_5(x,b,mu,zeta,hadron)
   Integrand(0)=0d0
  
  CASE(2) !!! uTMDPDF  quarks+gluon
   Integrand=uTMDPDF_50(x,b,mu,zeta,hadron)
   
  CASE(3) !!! uTMDPDF  quarks OPTIMAL
   Integrand=uTMDPDF_5(x,b,hadron)
   Integrand(0)=0d0
  
  CASE(4) !!! uTMDPDF  quarks+gluon OPTIMAL
   Integrand=uTMDPDF_50(x,b,hadron)
   
  CASE(5) !!! uTMDFF  quarks
   Integrand=uTMDFF_5(x,b,mu,zeta,hadron)
   Integrand(0)=0d0
  
  CASE(6) !!! uTMDFF  quarks+gluon
   Integrand=uTMDFF_50(x,b,mu,zeta,hadron)
   
  CASE(7) !!! uTMDFF  quarks OPTIMAL
   Integrand=uTMDFF_5(x,b,hadron)
   Integrand(0)=0d0
  
  CASE(8) !!! uTMDFF  quarks+gluon OPTIMAL
   Integrand=uTMDFF_50(x,b,hadron)
   
  CASE(9) !!! lin.pol.gluon TMDPDF
  !!! minus is due to definition (see manual)
   Integrand=-lpTMDPDF_50(x,b,mu,zeta,hadron)
  
  CASE(10) !!! lim.pol.gluon TMDPDF OPTIMAL
  !!! minus is due to definition (see manual)
   Integrand=-lpTMDPDF_50(x,b,hadron)
   
 END SELECT
 
 end function Integrand
 
  !-------------------------------------------------functions for uTMDPDF---------------------------------------
end module TMDs_inKT