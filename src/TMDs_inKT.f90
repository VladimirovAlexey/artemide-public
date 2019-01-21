!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 1.4
!
!	Evaluation of the TMDs in KT space
!	
!	if you use this module please, quote 190?.?????
!
!	ver 1.0: release (AV, 23.12.2018)
!
!				A.Vladimirov (23.12.2018)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDs_inKT
  use TMDs
  implicit none
  
  private
!   public
 
 !Current version of module
 character (len=5),parameter :: version="v1.40"
 
 !------------------------------------------Tables-----------------------------------------------------------------------
    integer,parameter::Nmax=200
    INCLUDE 'Tables/BesselZero.f90'
 
   logical:: convergenceLost=.false.
 
   real*8::hOGATA,tolerance
  !!!weights of ogata quadrature
  real*8,dimension(1:Nmax)::ww
  !!!nodes of ogata quadrature
  real*8,dimension(1:Nmax)::bb
  
  !Counters of calls "Integrand", Global in between setNP
  integer::GlobalCounter
  !Counter of total calls of TMD_F (reset at setNP)
  integer::CallCounter
  !Counter of maximum number of calls per integral
  integer::MaxCounter
 
!------------------------------------------Physical and mathematical constants------------------------------------------
  
!------------------------------------------Working variables------------------------------------------------------------
  
  logical::started=.false.
  
  !! Level of output
  !! 0=only critical
  !! 1=initialization details
  !! 2=WARNINGS
  integer::outputLevel=2
  
  public::TMDs_inKT_Initialize,TMDs_inKT_SetNPParameters,TMDs_inKT_SetScaleVariations,&
	  TMDs_inKT_ShowStatistic
	  
  real*8,dimension(-5:5),public::uTMDPDF_kT_50,uTMDPDF_kT_5,uTMDFF_kT_5,uTMDFF_kT_50
  
  interface TMDs_inKT_SetNPParameters
    module procedure TMDs_inKT_SetNPParameters,TMDs_inKT_SetNPParameters_rep
  end interface
  
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

contains 
  !!! This subroutine can be called only ones per programm run in the very beginning
  !!! It initializes TMDs subpackage
  !!! It set the pertrubative orders(!), so pertrubative orders cannot be changed afterwards.
  subroutine TMDs_inKT_Initialize(orderMain)
    character(len=*)::orderMain
    character(256)::line
    integer:: j
    
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
    
    if(outputLevel>1) write(*,*) '----- arTeMiDe.TMDs-in-KT ',version,': .... initialization'
    
    do
    read(51,'(A)') line    
    if(line(1:3)=='*2 ') exit
    end do    
    do
    read(51,'(A)') line
    if(line(1:3)=='*BB') exit
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
    
    
    if(started) then
    if(outputLevel>1) write(*,*) 'arTeMiDe.TMDs-in-KT already initialized'
    else
    
    
    CLOSE (51, STATUS='KEEP') 
    
    call TMDs_Initialize(orderMain)
    
      started=.true.
     if(outputLevel>0) write(*,*) '----- arTeMiDe.TMDs-inKT ',version,': .... initialized'
     if(outputLevel>1) write(*,*) ' '
    end if
  end subroutine TMDs_inKT_Initialize

  !sets the NP parameters of TMD model (bmax,gB,etc)
  subroutine TMDs_inKT_SetNPParameters(lambda)
    real*8 :: lambda(:)
    
   call TMDs_SetNPParameters(lambda)

  end subroutine TMDs_inKT_SetNPParameters
  
    !sets the NP parameters of TMD model (in replicas)
  subroutine TMDs_inKT_SetNPParameters_rep(n)
    integer::n
    
   call TMDs_SetNPParameters(n)

  end subroutine TMDs_inKT_SetNPParameters_rep
 
  !!!! this routine set the variations of scales
  !!!! it is used for the estimation of errors
  subroutine TMDs_inKT_SetScaleVariations(c1_in,c3_in,c4_in)
    real*8::c1_in,c3_in,c4_in
    
    call TMDs_SetScaleVariations(c1_in,c3_in,c4_in)
    
  end subroutine TMDs_inKT_SetScaleVariations
  
    !!!!!!!Functions which carry the trigger on convergences.... Its used in xSec, and probably in other places.
  function TMDs_inKT_IsconvergenceLost()
  logical::TMDs_inKT_IsconvergenceLost
  !!! checks TMDs trigger
  if(TMDs_IsconvergenceLost()) convergenceLost=.true.
  
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
 subroutine PrepareTables()
  real*8,parameter::piHalf=1.5707963267948966d0
  real*8,parameter::pi=3.141592653589793d0
  integer::i
  real*8::t!=h*xi
  real*8::psiPart!=tanh[pi/2 Sinh[h xi]]
  
  do i=1,Nmax
    t=hOGATA*JZero(0,i)
    psiPart=Tanh(piHalf*Sinh(t))
      bb(i)=JZero(0,i)*psiPart
      ww(i)=BESSEL_J0(bb(i))/JZero(0,i)/pi/(BESSEL_J1(JZero(0,i))**2)*(pi*t*Cosh(t)+Sinh(pi*Sinh(t)))/(1d0+Cosh(pi*Sinh(t)))
!      write(*,*) psiPart,b(k,i),w(k,i)
  end do
 
 end subroutine PrepareTables
 
 !--------------------------------------INTERFACES TO TMD------------------------------------------------
 
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
 
  !!!This is the defining module function
 !!! It evaluates the integral 
 !!!  int_0^infty   b db/2pi  J0(b qT) F1
 !!! The function F1 is given via number.. num
 function Fourier(x,qT_in,mu,zeta,num,hadron)
  real*8::qT,x,mu,zeta,qT_in
  integer::num,hadron
  real*8::integral(-5:5),eps(-5:5)
  real*8::v1,v2,v3,v4
  integer::k
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
  
  do k=1,Nmax!!! maximum of number of bessel roots preevaluated in the head
    eps=ww(k)*bb(k)*Integrand(bb(k)/qT,x,mu,zeta,num,hadron)
    
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
      write(*,*) 'bt/qT= ',bb(Nmax)/qT, 'qT=',qT
      write(*,*) 'W=',Integrand(bb(Nmax)/qT,x,mu,zeta,num,hadron), 'eps/integral =', eps/integral
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
   
 END SELECT
 
 end function Integrand
 
  !-------------------------------------------------functions for uTMDPDF---------------------------------------
end module TMDs_inKT