!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 1.41
!
!    Evaluation of the TMD structure function
!
!    if you use this module please, quote 1706.01473
!
!    ver 1.31: release (AV, 30.05.2018)
!    ver 1.41: fixed potential bug in the initialisation order (AV, 28.02.2019)
!    ver 3.00: removal of TMDs-module (AV, 12.02.2024)
!
!                A.Vladimirov (30.05.2018)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module TMDF_ogata
INCLUDE 'Code/TMDF/Fourier_byOgata.f90'
end module TMDF_ogata

module TMDF
use aTMDe_Numerics
use IO_functions
use TMDF_ogata
use EWinput
use uTMDPDF
use uTMDFF
use lpTMDPDF
use SiversTMDPDF
use wgtTMDPDF
use BoerMuldersTMDPDF

implicit none

private
!   public

character (len=7),parameter :: moduleName="TMDF"
character (len=5),parameter :: version="v3.01"
!Last appropriate verion of constants-file
integer,parameter::inputver=31

!------------------------------------------Working variables------------------------------------------------------------
integer::outputLevel=2
!! variable that count number of WRNING mesagges. In order not to spam too much
integer::messageTrigger=6
logical::started=.false.

logical:: convergenceLost=.false.

logical::include_uTMDPDF
logical::include_uTMDFF
logical::include_lpTMDPDF
logical::include_SiversTMDPDF
logical::include_wgtTMDPDF
logical::include_BoerMuldersTMDPDF

real(dp):: global_mass_scale=0.938_dp
real(dp):: qtMIN=0.0001d0
real(dp):: HardScaleMIN=0.8d0

integer::messageCounter
!-----------------------------------------Public interface--------------------------------------------------------------
public::TMDF_Initialize,TMDF_ResetCounters
public:: TMDF_F
public::TMDF_convergenceISlost,TMDF_IsconvergenceLost,TMDF_IsInitialized

public::Integrand

contains
 
function TMDF_IsInitialized()
    logical::TMDF_IsInitialized
    TMDF_IsInitialized=started
end function TMDF_IsInitialized

   !! Initialization of the package
subroutine TMDF_Initialize(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=300)::path
    logical::initRequired
    integer::FILEver
    real(dp)::hOGATA,tolerance

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
        write(*,*) '             Update the const-file with artemide.setup'
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
    read(51,*) initRequired
    if(.not.initRequired) then
        if(outputLevel>2) write(*,*)'artemide.',moduleName,': initialization is not required. '
        started=.false.
        return
    end if
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) tolerance
    call MoveTO(51,'*p2  ')
    read(51,*) hOGATA
    call MoveTO(51,'*p3  ')
    read(51,*) qtMIN
    call MoveTO(51,'*p4  ')
    read(51,*) HardScaleMIN
    
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) global_mass_scale

    if(outputLevel>2) write(*,'(A,ES8.2)') ' | h for Ogata quadrature    : ',hOGATA
    if(outputLevel>2) write(*,'(A,ES8.2)') ' | tolerance            : ',tolerance

    CLOSE (51, STATUS='KEEP') 

    !!! then we read it again from the beginning to fill parameters of other modules
    OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")

    !! uTMDPDF
    call MoveTO(51,'*4   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_uTMDPDF

    !! uTMDFF
    call MoveTO(51,'*5   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_uTMDFF

    !! lpTMDPDF
    call MoveTO(51,'*11  ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_lpTMDPDF

    !! SiversTMDPDF
    call MoveTO(51,'*12  ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_SiversTMDPDF

    !! wgtTMDPDF
    call MoveTO(51,'*13  ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_wgtTMDPDF

    !! BoerMuldersTMDPDF
    call MoveTO(51,'*14  ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_BoerMuldersTMDPDF

    CLOSE (51, STATUS='KEEP')

    call PrepareTables(tolerance,hOGATA)

    convergenceLost=.false.
    messageCounter=0

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!! Initialization of submodules !!!!!!!!!!!!!!!!!!!!!!
    if(.not.EWinput_IsInitialized()) then
        if(outputLevel>1) write(*,*) '.. initializing EWinput (from ',moduleName,')'
        if(present(prefix)) then
            call EWinput_Initialize(file,prefix)
        else
            call EWinput_Initialize(file)
        end if
    end if

    if(include_uTMDPDF .and. (.not.uTMDPDF_IsInitialized())) then
        if(outputLevel>1) write(*,*) '.. initializing uTMDPDF (from ',moduleName,')'
        if(present(prefix)) then
            call uTMDPDF_Initialize(file,prefix)
        else
            call uTMDPDF_Initialize(file)
        end if
    end if

    if(include_uTMDFF .and. (.not.uTMDFF_IsInitialized())) then
        if(outputLevel>1) write(*,*) '.. initializing uTMDFF (from ',moduleName,')'
        if(present(prefix)) then
            call uTMDFF_Initialize(file,prefix)
        else
            call uTMDFF_Initialize(file)
        end if
    end if

    if(include_lpTMDPDF .and. (.not.lpTMDPDF_IsInitialized())) then
        if(outputLevel>1) write(*,*) '.. initializing lpTMDPDF (from ',moduleName,')'
        if(present(prefix)) then
            call lpTMDPDF_Initialize(file,prefix)
        else
            call lpTMDPDF_Initialize(file)
        end if
    end if

    if(include_SiversTMDPDF .and. (.not.SiversTMDPDF_IsInitialized())) then
        if(outputLevel>1) write(*,*) '.. initializing SiversTMDPDF (from ',moduleName,')'
        if(present(prefix)) then
            call SiversTMDPDF_Initialize(file,prefix)
        else
            call SiversTMDPDF_Initialize(file)
        end if
    end if

    if(include_wgtTMDPDF .and. (.not.wgtTMDPDF_IsInitialized())) then
        if(outputLevel>1) write(*,*) '.. initializing wgtTMDPDF (from ',moduleName,')'
        if(present(prefix)) then
            call wgtTMDPDF_Initialize(file,prefix)
        else
            call wgtTMDPDF_Initialize(file)
        end if
    end if

    if(include_BoerMuldersTMDPDF .and. (.not.BoerMuldersTMDPDF_IsInitialized())) then
        if(outputLevel>1) write(*,*) '.. initializing SiversTMDPDF (from ',moduleName,')'
        if(present(prefix)) then
            call BoerMuldersTMDPDF_Initialize(file,prefix)
        else
            call BoerMuldersTMDPDF_Initialize(file)
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
  
!passes the NP parameters to TMDs
subroutine TMDF_ResetCounters()
    convergenceLost=.false.
    messageCounter=0
end subroutine TMDF_ResetCounters

!!!This is the defining module function
!!! It evaluates the integral
!!!  int_0^infty   b^(n+1) db/2  Jn(b qT) zff F1 F2
!!!
function TMDF_F(Q2,qT,x1,x2,mu,zeta1,zeta2,process) result(integral_result)
real(dp)::integral_result
real(dp),intent(in)::qT,x1,x2,mu,zeta1,zeta2,Q2
integer,dimension(1:3),intent(in)::process

integer::n
logical::ISconvergent
if(x1>=1d0 .or. x2>=1d0) then
  integral_result=0d0
else if(Q2<HardScaleMIN .or. mu<HardScaleMIN .or. zeta1<HardScaleMIN .or. zeta2<HardScaleMIN) then
    write(*,*) ErrorString("Some hard scale is too low (< "//trim(real8ToStr(HardScaleMIN))//")",moduleName)
    write(*,*) "(Q2 ,mu, zeta1,zeta2) = ",Q2 ,mu, zeta1,zeta2
    error stop ErrorString("Evaluation stop",moduleName)
else if(TMDF_IsconvergenceLost()) then
  !!!in the case of lost convergence we return huge number (divergent xSec)
    integral_result=1d10
else

!!Here we set the order of Bessel
if(process(3)<10000) then
n=0
else if(process(3)<20000) then
n=1
else if(process(3)<30000) then
n=2
else
n=3
end if

if(qT<=qtMIN) then
    call Fourier_byOgata(n,F_toInt,qtMIN,integral_result,ISconvergent)
else
    call Fourier_byOgata(n,F_toInt,qT,integral_result,ISconvergent)
end if

if(.not.ISconvergent) call TMDF_convergenceISlost
end if

contains

function F_toInt(b)
real(dp)::F_toInt
real(dp),intent(in)::b

F_toInt=Integrand(Q2,b,x1,x2,mu,zeta1,zeta2,process)

end function F_toInt


end function TMDF_F

function Integrand(Q2,b,x1,x2,mu,zeta1,zeta2,process_array)
 real(dp)::Integrand
 real(dp),intent(in)::b,x1,x2,mu,zeta1,zeta2,Q2
 integer,dimension(1:3),intent(in)::process_array
 integer::process,h,h1,h2
 real(dp),dimension(-5:5)::FA,FB,FAB
 
 process=process_array(3) !! general process name
 h1=process_array(1)!! first hadron
 h2=process_array(2)!! second hadron
 
 if(b>1000d0) then
  Integrand=0d0
  return
 end if
  
 SELECT CASE(process)
  !!!test cases
  CASE(0,10000,20000,30000)
    !Integrand=Exp(-0.2d0*b)
    Integrand=b**4*Exp(-5.d0*b)
  CASE(9999,19999,29999,39999)
    Integrand=Exp(-mu*b)*(1d0+x1*b**2+x2*b**4)
  CASE(9998,19998,29998,39998)
    Integrand=Exp(-mu*b**2)*(1d0+x1*b**2+x2*b**4)
!--------------------------------------------------------------------------------
  CASE (1) !h1+h2 -> gamma
    ! e_q^2 *F_q(A)*F_qbar(B)
    if(zeta1==zeta2) then
      FAB=uPDF_uPDF(x1,x2,b,mu,zeta1,h1,h2)!!! -h2, to multiply quarks by anti-quarks in FAB
    else
     FA=uTMDPDF_inB(x1,b,mu,zeta1,h1)
     FB=uTMDPDF_inB(x2,b,mu,zeta2,-h2)!!! -h2, to multiply quarks by anti-quarks in FAB

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
  CASE (2) !h1+h2->Z
      !((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) *F_q(A)*F_qbar(B)
    if(zeta1==zeta2) then
      FAB=uPDF_uPDF(x1,x2,b,mu,zeta1,h1,h2)
    else
     FA=uTMDPDF_inB(x1,b,mu,zeta1,h1)
     FB=uTMDPDF_inB(x2,b,mu,zeta2,-h2)!!! -h2, to multiply quarks by anti-quarks in FAB

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
  CASE (3) !h1+h2->Z+gamma
    if(zeta1==zeta2) then
      FAB=uPDF_uPDF(x1,x2,b,mu,zeta1,h1,h2)
    else
     FA=uTMDPDF_inB(x1,b,mu,zeta1,h1)
     FB=uTMDPDF_inB(x2,b,mu,zeta2,-h2)!!! -h2, to multiply quarks by anti-quarks in FAB

     FAB=FA*FB
    end if

    Integrand=XIntegrandForDYwithZgamma(FAB,Q2)
!--------------------------------------------------------------------------------  
  CASE (4) !h1+h2-> W+
     FA=uTMDPDF_inB(x1,b,mu,zeta1,h1)
     FB=uTMDPDF_inB(x2,b,mu,zeta2,h2)

    Integrand=paramW_L*(&
    paramW_UD*(FA(2)*FB(-1)+FA(-1)*FB(2))&        !u*dbar+dbar*u
    +paramW_US*(FA(2)*FB(-3)+FA(-3)*FB(2))&        !u*sbar+sbar*u
    +paramW_UB*(FA(2)*FB(-5)+FA(-5)*FB(2))&        !u*bbar+bbar*u
    +paramW_CD*(FA(4)*FB(-1)+FA(-1)*FB(4))&        !c*dbar+dbar*c
    +paramW_CS*(FA(4)*FB(-3)+FA(-3)*FB(4))&        !c*sbar+sbar*c
    +paramW_CB*(FA(4)*FB(-5)+FA(-5)*FB(4))&        !c*bbar+bbar*c
    )*Q2*Q2/((Q2-MW2)**2+GammaW2*MW2)
!--------------------------------------------------------------------------------  
  CASE (5) !h1+h2-> W-
     FA=uTMDPDF_inB(x1,b,mu,zeta1,h1)
     FB=uTMDPDF_inB(x2,b,mu,zeta2,h2)

    Integrand=paramW_L*(&
    paramW_UD*(FA(1)*FB(-2)+FA(-2)*FB(1))&        !d*ubar+ubar*d
    +paramW_US*(FA(3)*FB(-2)+FA(-2)*FB(3))&        !s*ubar+ubar*s
    +paramW_UB*(FA(5)*FB(-2)+FA(-2)*FB(5))&        !b*ubar+ubar*b
    +paramW_CD*(FA(1)*FB(-4)+FA(-4)*FB(1))&        !d*cbar+cbar*d
    +paramW_CS*(FA(3)*FB(-4)+FA(-4)*FB(3))&        !s*cbar+cbar*s
    +paramW_CB*(FA(5)*FB(-4)+FA(-4)*FB(5))&        !b*cbar+cbar*b
    )*Q2*Q2/((Q2-MW2)**2+GammaW2*MW2)
!--------------------------------------------------------------------------------  
  CASE (6) !h1+h2-> W+ + W-
     FA=uTMDPDF_inB(x1,b,mu,zeta1,h1)
     FB=uTMDPDF_inB(x2,b,mu,zeta2,h2)

    Integrand=paramW_L*(&
    paramW_UD*(FA(2)*FB(-1)+FA(1)*FB(-2)+FA(-2)*FB(1)+FA(-1)*FB(2))&    !u*dbar+d*ubar+ubar*d+dbar*u
    +paramW_US*(FA(2)*FB(-3)+FA(3)*FB(-2)+FA(-2)*FB(3)+FA(-3)*FB(2))&    !u*sbar+s*ubar+ubar*s+sbar*u
    +paramW_UB*(FA(2)*FB(-5)+FA(5)*FB(-2)+FA(-2)*FB(5)+FA(-5)*FB(2))&    !u*bbar+b*ubar+ubar*b+bbar*u
    +paramW_CD*(FA(4)*FB(-1)+FA(1)*FB(-4)+FA(-4)*FB(1)+FA(-1)*FB(4))&    !c*dbar+d*cbar+cbar*d+dbar*c
    +paramW_CS*(FA(4)*FB(-3)+FA(3)*FB(-4)+FA(-4)*FB(3)+FA(-3)*FB(4))&    !c*sbar+s*cbar+cbar*s+sbar*c
    +paramW_CB*(FA(4)*FB(-5)+FA(5)*FB(-4)+FA(-4)*FB(5)+FA(-5)*FB(4))&    !c*bbar+b*cbar+cbar*b+bbar*c
    )*Q2*Q2/((Q2-MW2)**2+GammaW2*MW2)

!--------------------------------------------------------------------------------  
  CASE (7) !h1+h2-> W+ (for zero-width)
     FA=uTMDPDF_inB(x1,b,mu,zeta1,h1)
     FB=uTMDPDF_inB(x2,b,mu,zeta2,h2)

    Integrand=&
    paramW_UD*(FA(2)*FB(-1)+FA(-1)*FB(2))&        !u*dbar+dbar*u
    +paramW_US*(FA(2)*FB(-3)+FA(-3)*FB(2))&        !u*sbar+sbar*u
    +paramW_UB*(FA(2)*FB(-5)+FA(-5)*FB(2))&        !u*bbar+bbar*u
    +paramW_CD*(FA(4)*FB(-1)+FA(-1)*FB(4))&        !c*dbar+dbar*c
    +paramW_CS*(FA(4)*FB(-3)+FA(-3)*FB(4))&        !c*sbar+sbar*c
    +paramW_CB*(FA(4)*FB(-5)+FA(-5)*FB(4))        !c*bbar+bbar*c

!--------------------------------------------------------------------------------  
  CASE (8) !h1+h2-> W- (for zero-width)
     FA=uTMDPDF_inB(x1,b,mu,zeta1,h1)
     FB=uTMDPDF_inB(x2,b,mu,zeta2,h2)

    Integrand=&
    paramW_UD*(FA(1)*FB(-2)+FA(-2)*FB(1))&        !d*ubar+ubar*d
    +paramW_US*(FA(3)*FB(-2)+FA(-2)*FB(3))&        !s*ubar+ubar*s
    +paramW_UB*(FA(5)*FB(-2)+FA(-2)*FB(5))&        !b*ubar+ubar*b
    +paramW_CD*(FA(1)*FB(-4)+FA(-4)*FB(1))&        !d*cbar+cbar*d
    +paramW_CS*(FA(3)*FB(-4)+FA(-4)*FB(3))&        !s*cbar+cbar*s
    +paramW_CB*(FA(5)*FB(-4)+FA(-4)*FB(5))        !b*cbar+cbar*b

!--------------------------------------------------------------------------------  
  CASE (9) !h1+h2-> W+- (for zero-width)
     FA=uTMDPDF_inB(x1,b,mu,zeta1,h1)
     FB=uTMDPDF_inB(x2,b,mu,zeta2,h2)

    Integrand=&
    paramW_UD*(FA(2)*FB(-1)+FA(1)*FB(-2)+FA(-2)*FB(1)+FA(-1)*FB(2))&    !u*dbar+d*ubar+ubar*d+dbar*u
    +paramW_US*(FA(2)*FB(-3)+FA(3)*FB(-2)+FA(-2)*FB(3)+FA(-3)*FB(2))&    !u*sbar+s*ubar+ubar*s+sbar*u
    +paramW_UB*(FA(2)*FB(-5)+FA(5)*FB(-2)+FA(-2)*FB(5)+FA(-5)*FB(2))&    !u*bbar+b*ubar+ubar*b+bbar*u
    +paramW_CD*(FA(4)*FB(-1)+FA(1)*FB(-4)+FA(-4)*FB(1)+FA(-1)*FB(4))&    !c*dbar+d*cbar+cbar*d+dbar*c
    +paramW_CS*(FA(4)*FB(-3)+FA(3)*FB(-4)+FA(-4)*FB(3)+FA(-3)*FB(4))&    !c*sbar+s*cbar+cbar*s+sbar*c
    +paramW_CB*(FA(4)*FB(-5)+FA(5)*FB(-4)+FA(-4)*FB(5)+FA(-5)*FB(4))    !c*bbar+b*cbar+cbar*b+bbar*c

!--------------------------------------------------------------------------------  
  CASE(10) !h1+h2 -> Higgs (unpol.part+lin.pol.part)
    FA=uTMDPDF_inB(x1,b,mu,zeta1,h1)
    FB=uTMDPDF_inB(x2,b,mu,zeta2,h2)
    Integrand=FA(0)*FB(0) !!!! unpolarized part

    FA=lpTMDPDF_inB(x1,b,mu,zeta1,h1)
    FB=lpTMDPDF_inB(x2,b,mu,zeta2,h2)
    Integrand=Integrand+FA(0)*FB(0) !!!! linearly polarized part
!--------------------------------------------------------------------------------  
  CASE(11) !h1+h2 -> Higgs (unpol.part)
    FA=uTMDPDF_inB(x1,b,mu,zeta1,h1)
    FB=uTMDPDF_inB(x2,b,mu,zeta2,h2)
    Integrand=FA(0)*FB(0)
  
!--------------------------------------------------------------------------------  
  CASE(12) !h1+h2 -> Higgs (lin.pol.part)
    FA=lpTMDPDF_inB(x1,b,mu,zeta1,h1)
    FB=lpTMDPDF_inB(x2,b,mu,zeta2,h2)
    Integrand=FA(0)*FB(0)

!--------------------------------------------------------------------------------
  CASE(24) !h1+h2 -> A4
    if(zeta1==zeta2) then
      FAB=uPDF_uPDF(x1,x2,b,mu,zeta1,h1,h2)
    else
     FA=uTMDPDF_inB(x1,b,mu,zeta1,h1)
     FB=uTMDPDF_inB(x2,b,mu,zeta2,-h2)!!! -h2, to multiply quarks by anti-quarks in FAB

     FAB=FA*FB
    end if
    Integrand=-2*XTMD_pairZmZm_A(FAB,Q2)

!--------------------------------------------------------------------------------  
  CASE (101) !h1+Cu->gamma* !!this is for E288
    !!!! strictly hadron 1
    FA=uTMDPDF_inB(x1,b,mu,zeta1,h1)
    FB=uTMDPDF_inB(x2,b,mu,zeta2,1)
    Integrand=116d0/567d0*(FA(2)*FB(-2)+FA(-2)*FB(2))+136d0/567d0*(FA(-2)*FB(1)+FA(2)*FB(-1))&
          +34d0/567d0*(FA(-1)*FB(2)+FA(1)*FB(-2))+29d0/567d0*(FA(-1)*FB(1)+FA(1)*FB(-1))&
          +1d0/9d0*(FA(-3)*FB(3)+FA(3)*FB(-3)+4d0*FA(-4)*FB(4)+4d0*FA(4)*FB(-4)+FA(-5)*FB(5)+FA(5)*FB(-5))

  !--------------------------------------------------------------------------------  
  CASE (102) !h1+2H->gamma* !!this is for E772
    !!!! strictrly hadron 1
    FA=uTMDPDF_inB(x1,b,mu,zeta1,h1)
    FB=uTMDPDF_inB(x2,b,mu,zeta2,1)
    Integrand=2d0/9d0*(FA(2)*FB(-2)+FA(-2)*FB(2))+2d0/9d0*(FA(-2)*FB(1)+FA(2)*FB(-1))&
          +1d0/18d0*(FA(-1)*FB(2)+FA(1)*FB(-2))+1d0/18d0*(FA(-1)*FB(1)+FA(1)*FB(-1))&
          +1d0/9d0*(FA(-3)*FB(3)+FA(3)*FB(-3)+4d0*FA(-4)*FB(4)+4d0*FA(4)*FB(-4)+FA(-5)*FB(5)+FA(5)*FB(-5))
  !--------------------------------------------------------------------------------  
  CASE (103) !h1+W->gamma* !!this is for E537
    !!!! strictrly hadron 1
    !Wolfram has A=183,    Z=74,    N=109
    FA=uTMDPDF_inB(x1,b,mu,zeta1,h1)
    FB=uTMDPDF_inB(x2,b,mu,zeta2,1)
    Integrand=296d0/1647d0*(FA(-2)*FB(2)+FA(2)*FB(-2))+436d0/1647d0*(FA(-2)*FB(1)+FA(2)*FB(-1))&
        +109d0/1647d0*(FA(-1)*FB(2)+FA(1)*FB(-2))+74d0/1647d0*(FA(-1)*FB(1)+FA(1)*FB(-1))&
        +1d0/9d0*(FA(-3)*FB(3)+FA(3)*FB(-3)+4d0*FA(-4)*FB(4)+4d0*FA(4)*FB(-4)+FA(-5)*FB(5)+FA(5)*FB(-5))
  !----------------------------------------------------------------------------------
  !-------------------------SIDIS----------------------------------------------------
  !----------------------------------------------------------------------------------
  CASE (2001) !h1->h2 where !!!! unpolarized SIDIS
    ! e_q^2 *F_q(A)*F_q(B)
    FA=uTMDPDF_inB(x1,b,mu,zeta1,h1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,h2)
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
    CASE (2002) !d->h2 where d is deutron prepared from hadron 1 [i.e u->(u+d)/2, d->(u+d)/2]
    ! e_q^2 *F_q(A)*F_q(B)
    FA=uTMDPDF_inB(x1,b,mu,zeta1,h1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,h2)
    Integrand=(FA(1)+FA(2))*(FB(1)+4d0*FB(2))/18d0&
      +FA(3)*FB(3)/9.d0&
      +FA(4)*FB(4)*4d0/9.d0&
      +FA(5)*FB(5)/9d0&
      +(FA(-1)+FA(-2))*(FB(-1)+4d0*FB(-2))/18d0&
      +FA(-3)*FB(-3)/9.d0&
      +FA(-4)*FB(-4)*4d0/9.d0&
      +FA(-5)*FB(-5)/9d0
!--------------------------------------------------------------------------------  
    CASE (2003) !n->h2 where n=last number (n=neutron=p(u<->d))
    ! e_q^2 *F_q(A)*F_q(B)
    FA=uTMDPDF_inB(x1,b,mu,zeta1,h1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,h2)
    Integrand=FA(2)*FB(1)/9.d0&
      +FA(1)*FB(2)*4.d0/9.d0&
      +FA(3)*FB(3)/9.d0&
      +FA(4)*FB(4)*4d0/9.d0&
      +FA(5)*FB(5)/9d0&
      +FA(-2)*FB(-1)/9.d0&
      +FA(-1)*FB(-2)*4.d0/9.d0&
      +FA(-3)*FB(-3)/9.d0&
      +FA(-4)*FB(-4)*4d0/9.d0&
      +FA(-5)*FB(-5)/9d0
!--------------------------------------------------------------------------------  
   CASE (2101) !p->h? where h?=h1+h2
    ! e_q^2 *F_q(A)*F_q(B)
    FA=uTMDPDF_inB(x1,b,mu,zeta1,1)
    if(h2>0) then
        FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)
    else
        FB=uTMDFF_inB(x2,b,mu,zeta2,-1)+uTMDFF_inB(x2,b,mu,zeta2,-2)
    end if
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
    FA=uTMDPDF_inB(x1,b,mu,zeta1,1)
    if(h2>0) then
        FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)+uTMDFF_inB(x2,b,mu,zeta2,3)
    else
        FB=uTMDFF_inB(x2,b,mu,zeta2,-1)+uTMDFF_inB(x2,b,mu,zeta2,-2)+uTMDFF_inB(x2,b,mu,zeta2,-3)
    end if
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
    FA=uTMDPDF_inB(x1,b,mu,zeta1,1)
    if(h2>0) then
        FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)
    else
        FB=uTMDFF_inB(x2,b,mu,zeta2,-1)+uTMDFF_inB(x2,b,mu,zeta2,-2)
    end if
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
    FA=uTMDPDF_inB(x1,b,mu,zeta1,1)
    if(h2>0) then
        FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)+uTMDFF_inB(x2,b,mu,zeta2,3)
    else
        FB=uTMDFF_inB(x2,b,mu,zeta2,-1)+uTMDFF_inB(x2,b,mu,zeta2,-2)+uTMDFF_inB(x2,b,mu,zeta2,-3)
    end if
    Integrand=(FA(1)+FA(2))*(FB(1)+4d0*FB(2))/18d0&
      +FA(3)*FB(3)/9.d0&
      +FA(4)*FB(4)*4d0/9.d0&
      +FA(5)*FB(5)/9d0&
      +(FA(-1)+FA(-2))*(FB(-1)+4d0*FB(-2))/18d0&
      +FA(-3)*FB(-3)/9.d0&
      +FA(-4)*FB(-4)*4d0/9.d0&
      +FA(-5)*FB(-5)/9d0
!--------------------------------------------------------------------------------  
   CASE (2105) !n->h? where h?=h1+h2 (n=neutron=p(u<->d))
    ! e_q^2 *F_q(A)*F_q(B)
    FA=uTMDPDF_inB(x1,b,mu,zeta1,1)
    if(h2>0) then
        FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)
    else
        FB=uTMDFF_inB(x2,b,mu,zeta2,-1)+uTMDFF_inB(x2,b,mu,zeta2,-2)
    end if
    Integrand=FA(2)*FB(1)/9.d0&
      +FA(1)*FB(2)*4.d0/9.d0&
      +FA(3)*FB(3)/9.d0&
      +FA(4)*FB(4)*4d0/9.d0&
      +FA(5)*FB(5)/9d0&
      +FA(-2)*FB(-1)/9.d0&
      +FA(-1)*FB(-2)*4.d0/9.d0&
      +FA(-3)*FB(-3)/9.d0&
      +FA(-4)*FB(-4)*4d0/9.d0&
      +FA(-5)*FB(-5)/9d0
!--------------------------------------------------------------------------------  
    CASE (2106) !n->h? where h?=h1+h2+h3 (n=neutron=p(u<->d))
    ! e_q^2 *F_q(A)*F_q(B)
    FA=uTMDPDF_inB(x1,b,mu,zeta1,1)
    if(h2>0) then
        FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)+uTMDFF_inB(x2,b,mu,zeta2,3)
    else
        FB=uTMDFF_inB(x2,b,mu,zeta2,-1)+uTMDFF_inB(x2,b,mu,zeta2,-2)+uTMDFF_inB(x2,b,mu,zeta2,-3)
    end if
    Integrand=FA(2)*FB(1)/9.d0&
      +FA(1)*FB(2)*4.d0/9.d0&
      +FA(3)*FB(3)/9.d0&
      +FA(4)*FB(4)*4d0/9.d0&
      +FA(5)*FB(5)/9d0&
      +FA(-2)*FB(-1)/9.d0&
      +FA(-1)*FB(-2)*4.d0/9.d0&
      +FA(-3)*FB(-3)/9.d0&
      +FA(-4)*FB(-4)*4d0/9.d0&
      +FA(-5)*FB(-5)/9d0
!--------------------------------------------------------------------------------
   CASE (2107) !p->h? where h?=h1+h2 [from 3+4]
    ! e_q^2 *F_q(A)*F_q(B)
    FA=uTMDPDF_inB(x1,b,mu,zeta1,1)
    if(h2>0) then
        FB=uTMDFF_inB(x2,b,mu,zeta2,3)+uTMDFF_inB(x2,b,mu,zeta2,4)
    else
        FB=uTMDFF_inB(x2,b,mu,zeta2,-3)+uTMDFF_inB(x2,b,mu,zeta2,-4)
    end if
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
    CASE (2108) !d->h? where h?=h1+h2 (d=deutron=(p+n)/2) [from 3+4]
    ! e_q^2 *F_q(A)*F_q(B)
    FA=uTMDPDF_inB(x1,b,mu,zeta1,1)
    if(h2>0) then
        FB=uTMDFF_inB(x2,b,mu,zeta2,3)+uTMDFF_inB(x2,b,mu,zeta2,4)
    else
        FB=uTMDFF_inB(x2,b,mu,zeta2,-3)+uTMDFF_inB(x2,b,mu,zeta2,-4)
    end if
    Integrand=(FA(1)+FA(2))*(FB(1)+4d0*FB(2))/18d0&
      +FA(3)*FB(3)/9.d0&
      +FA(4)*FB(4)*4d0/9.d0&
      +FA(5)*FB(5)/9d0&
      +(FA(-1)+FA(-2))*(FB(-1)+4d0*FB(-2))/18d0&
      +FA(-3)*FB(-3)/9.d0&
      +FA(-4)*FB(-4)*4d0/9.d0&
      +FA(-5)*FB(-5)/9d0
!--------------------------------------------------------------------------------
   CASE (2109) !n->h? where h?=h1+h2 (n=neutron=p(u<->d))[from 3+4]
    ! e_q^2 *F_q(A)*F_q(B)
    FA=uTMDPDF_inB(x1,b,mu,zeta1,1)
    if(h2>0) then
        FB=uTMDFF_inB(x2,b,mu,zeta2,3)+uTMDFF_inB(x2,b,mu,zeta2,4)
    else
        FB=uTMDFF_inB(x2,b,mu,zeta2,-3)+uTMDFF_inB(x2,b,mu,zeta2,-4)
    end if
    Integrand=FA(2)*FB(1)/9.d0&
      +FA(1)*FB(2)*4.d0/9.d0&
      +FA(3)*FB(3)/9.d0&
      +FA(4)*FB(4)*4d0/9.d0&
      +FA(5)*FB(5)/9d0&
      +FA(-2)*FB(-1)/9.d0&
      +FA(-1)*FB(-2)*4.d0/9.d0&
      +FA(-3)*FB(-3)/9.d0&
      +FA(-4)*FB(-4)*4d0/9.d0&
      +FA(-5)*FB(-5)/9d0
!----------------------------------------------------------------------------------------------------------------------------------  
!-----------------------------------------------------Sivers asymetries------------------------------------------------------------
!----------------------------------------------------------------------------------------------------------------------------------
!--------------------------------------------------------------------------------
    CASE (10001) !pp->gamma
        ! e_q^2 *F_q(A)*F_qbar(B)
        FA=-SiversTMDPDF_inB(x1,b,mu,zeta1,h1)  !!!! -1 is due to definition of Sivers function (+1) for SIDIS (-1) for DY
        FB=uTMDPDF_inB(x2,b,mu,zeta2,-h2) !!! -h2, to multiply quarks by anti-quarks in FAB
        FAB=FA*FB
        
        Integrand=-global_mass_scale*(&
        FAB(1)/9.d0&
        +FAB(2)*4.d0/9.d0&
        +FAB(3)/9.d0&
        +FAB(4)*4d0/9.d0&
        +FAB(5)/9d0&
        +FAB(-1)/9.d0&
        +FAB(-2)*4.d0/9.d0&
        +FAB(-3)/9.d0&
        +FAB(-4)*4d0/9.d0&
        +FAB(-5)/9d0)
    !--------------------------------------------------------------------------------  
    CASE (10003) !pp->Z+gamma
        FA=-SiversTMDPDF_inB(x1,b,mu,zeta1,h1)  !!!! -1 is due to definition of Sivers function (+1) for SIDIS (-1) for DY
        FB=uTMDPDF_inB(x2,b,mu,zeta2,-h2) !!! -h2, to multiply quarks by anti-quarks in FAB
        FAB=FA*FB
            
        Integrand=-global_mass_scale*XIntegrandForDYwithZgamma(FAB,Q2)
    !--------------------------------------------------------------------------------  
    CASE (10004) !pp-> W+
        FA=-SiversTMDPDF_inB(x1,b,mu,zeta1,h1)  !!!! -1 is due to definition of Sivers function (+1) for SIDIS (-1) for DY
        FB=uTMDPDF_inB(x2,b,mu,zeta2,h2)
        
        Integrand=-global_mass_scale*paramW_L*(&
        paramW_UD*(FA(2)*FB(-1)+FA(-1)*FB(2))&        !u*dbar+dbar*u
        +paramW_US*(FA(2)*FB(-3)+FA(-3)*FB(2))&        !u*sbar+sbar*u
        +paramW_UB*(FA(2)*FB(-5)+FA(-5)*FB(2))&        !u*bbar+bbar*u
        +paramW_CD*(FA(4)*FB(-1)+FA(-1)*FB(4))&        !c*dbar+dbar*c
        +paramW_CS*(FA(4)*FB(-3)+FA(-3)*FB(4))&        !c*sbar+sbar*c
        +paramW_CB*(FA(4)*FB(-5)+FA(-5)*FB(4))&        !c*bbar+bbar*c
        )*Q2*Q2/((Q2-MW2)**2+GammaW2*MW2)
    !--------------------------------------------------------------------------------  
    CASE (10005) !pp-> W-
        FA=-SiversTMDPDF_inB(x1,b,mu,zeta1,h1)  !!!! -1 is due to definition of Sivers function (+1) for SIDIS (-1) for DY
        FB=uTMDPDF_inB(x2,b,mu,zeta2,h2)
        
        Integrand=-global_mass_scale*paramW_L*(&
        paramW_UD*(FA(1)*FB(-2)+FA(-2)*FB(1))&        !d*ubar+ubar*d
        +paramW_US*(FA(3)*FB(-2)+FA(-2)*FB(3))&        !s*ubar+ubar*s
        +paramW_UB*(FA(5)*FB(-2)+FA(-2)*FB(5))&        !b*ubar+ubar*b
        +paramW_CD*(FA(1)*FB(-4)+FA(-4)*FB(1))&        !d*cbar+cbar*d
        +paramW_CS*(FA(3)*FB(-4)+FA(-4)*FB(3))&        !s*cbar+cbar*s
        +paramW_CB*(FA(5)*FB(-4)+FA(-4)*FB(5))&        !b*cbar+cbar*b
        )*Q2*Q2/((Q2-MW2)**2+GammaW2*MW2)
    !--------------------------------------------------------------------------------  
    CASE (10011) !h+p(s)->gamma
    ! e_q^2 *F_q(A)*F_qbar(B)
    FA=uTMDPDF_inB(x1,b,mu,zeta2,h1)
    FB=-SiversTMDPDF_inB(x2,b,mu,zeta1,h2)    !!!! -1 is due to definition of Sivers function (+1) for SIDIS (-1) for DY
    FAB=FA*(FB(5:-5:-1))

    !!!! extra factor -1 (total +1) is due to definition of Sivers, h1+h2(s)=-h1(s)+h2
    Integrand=+global_mass_scale*(&
        FAB(1)/9.d0&
        +FAB(2)*4.d0/9.d0&
        +FAB(3)/9.d0&
        +FAB(4)*4d0/9.d0&
        +FAB(5)/9d0&
        +FAB(-1)/9.d0&
        +FAB(-2)*4.d0/9.d0&
        +FAB(-3)/9.d0&
        +FAB(-4)*4d0/9.d0&
        +FAB(-5)/9d0)

    !--------------------------------------------------------------------------------
    CASE (12011:12019) !Sivers asymmetry d->hN where n=last number (d=deutron=(p+n)/2)
    ! e_q^2 *F_q(A)*F_q(B)
    h=process-12010
    FA=SiversTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,h)
    Integrand=-global_mass_scale*(&
        (FA(1)+FA(2))*(FB(1)+4d0*FB(2))/18d0&
        +FA(3)*FB(3)/9.d0&
        +FA(4)*FB(4)*4d0/9.d0&
        +FA(5)*FB(5)/9d0&
        +(FA(-1)+FA(-2))*(FB(-1)+4d0*FB(-2))/18d0&
        +FA(-3)*FB(-3)/9.d0&
        +FA(-4)*FB(-4)*4d0/9.d0&
        +FA(-5)*FB(-5)/9d0)
    !--------------------------------------------------------------------------------
    CASE (12021:12029) !Sivers asymmetry p->bar-hN where n=last number
    ! e_q^2 *F_q(A)*F_bar-q(B)
    h=process-12020
    FA=SiversTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,h)
    Integrand=-global_mass_scale*(&
        FA(1)*FB(-1)/9.d0&
        +FA(2)*FB(-2)*4.d0/9.d0&
        +FA(3)*FB(-3)/9.d0&
        +FA(4)*FB(-4)*4d0/9.d0&
        +FA(5)*FB(-5)/9d0&
        +FA(-1)*FB(1)/9.d0&
        +FA(-2)*FB(2)*4.d0/9.d0&
        +FA(-3)*FB(3)/9.d0&
        +FA(-4)*FB(4)*4d0/9.d0&
        +FA(-5)*FB(5)/9d0)
    !--------------------------------------------------------------------------------
    CASE (12031:12039) !Sivers asymmetry d->bar-hN where n=last number (d=deutron=(p+n)/2)
    ! e_q^2 *F_q(A)*F_bar-q(B)
    h=process-12030
    FA=SiversTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,h)
    Integrand=-global_mass_scale*(&
        (FA(1)+FA(2))*(FB(-1)+4d0*FB(-2))/18d0&
        +FA(3)*FB(-3)/9.d0&
        +FA(4)*FB(-4)*4d0/9.d0&
        +FA(5)*FB(-5)/9d0&
        +(FA(-1)+FA(-2))*(FB(1)+4d0*FB(2))/18d0&
        +FA(-3)*FB(3)/9.d0&
        +FA(-4)*FB(4)*4d0/9.d0&
        +FA(-5)*FB(5)/9d0)
    !--------------------------------------------------------------------------------
    CASE (12041:12049) !Sivers asymmetry p->hN where n=last number  (n=neutron=p(u<->d))
    ! e_q^2 *F_q(A)*F_q(B)
    h=process-12040
    FA=SiversTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,h)
    Integrand=-global_mass_scale*(&
        FA(2)*FB(1)/9.d0&
        +FA(1)*FB(2)*4.d0/9.d0&
        +FA(3)*FB(3)/9.d0&
        +FA(4)*FB(4)*4d0/9.d0&
        +FA(5)*FB(5)/9d0&
        +FA(-2)*FB(-1)/9.d0&
        +FA(-1)*FB(-2)*4.d0/9.d0&
        +FA(-3)*FB(-3)/9.d0&
        +FA(-4)*FB(-4)*4d0/9.d0&
        +FA(-5)*FB(-5)/9d0)
    !--------------------------------------------------------------------------------
    CASE (12051:12059) !Sivers asymmetry p->bar-hN where n=last number (n=neutron=p(u<->d))
    ! e_q^2 *F_q(A)*F_bar-q(B)
    h=process-12050
    FA=SiversTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,h)
    Integrand=-global_mass_scale*(&
        FA(2)*FB(-1)/9.d0&
        +FA(1)*FB(-2)*4.d0/9.d0&
        +FA(3)*FB(-3)/9.d0&
        +FA(4)*FB(-4)*4d0/9.d0&
        +FA(5)*FB(-5)/9d0&
        +FA(-2)*FB(1)/9.d0&
        +FA(-1)*FB(2)*4.d0/9.d0&
        +FA(-3)*FB(3)/9.d0&
        +FA(-4)*FB(4)*4d0/9.d0&
        +FA(-5)*FB(5)/9d0)
    !--------------------------------------------------------------------------------  
    !--------------------------------------------------------------------------------  
    CASE (12101) !p->h? where h?=h1+h2
    ! e_q^2 *F_q(A)*F_q(B)
    FA=SiversTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)
    Integrand=-global_mass_scale*(&
        FA(1)*FB(1)/9.d0&
        +FA(2)*FB(2)*4.d0/9.d0&
        +FA(3)*FB(3)/9.d0&
        +FA(4)*FB(4)*4d0/9.d0&
        +FA(5)*FB(5)/9d0&
        +FA(-1)*FB(-1)/9.d0&
        +FA(-2)*FB(-2)*4.d0/9.d0&
        +FA(-3)*FB(-3)/9.d0&
        +FA(-4)*FB(-4)*4d0/9.d0&
        +FA(-5)*FB(-5)/9d0)
    !--------------------------------------------------------------------------------  
    CASE (12102) !p->h? where h?=h1+h2+h3
    ! e_q^2 *F_q(A)*F_q(B)
    FA=SiversTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)+uTMDFF_inB(x2,b,mu,zeta2,3)
    Integrand=-global_mass_scale*(&
        FA(1)*FB(1)/9.d0&
        +FA(2)*FB(2)*4.d0/9.d0&
        +FA(3)*FB(3)/9.d0&
        +FA(4)*FB(4)*4d0/9.d0&
        +FA(5)*FB(5)/9d0&
        +FA(-1)*FB(-1)/9.d0&
        +FA(-2)*FB(-2)*4.d0/9.d0&
        +FA(-3)*FB(-3)/9.d0&
        +FA(-4)*FB(-4)*4d0/9.d0&
        +FA(-5)*FB(-5)/9d0)
    !--------------------------------------------------------------------------------  
    CASE (12103) !d->h? where h?=h1+h2 (d=deutron=(p+n)/2)
    ! e_q^2 *F_q(A)*F_q(B)
    FA=SiversTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)
    Integrand=-global_mass_scale*(&
        (FA(1)+FA(2))*(FB(1)+4d0*FB(2))/18d0&
        +FA(3)*FB(3)/9.d0&
        +FA(4)*FB(4)*4d0/9.d0&
        +FA(5)*FB(5)/9d0&
        +(FA(-1)+FA(-2))*(FB(-1)+4d0*FB(-2))/18d0&
        +FA(-3)*FB(-3)/9.d0&
        +FA(-4)*FB(-4)*4d0/9.d0&
        +FA(-5)*FB(-5)/9d0)
    !--------------------------------------------------------------------------------  
    CASE (12104) !d->h? where h?=h1+h2+h3 (d=deutron=(p+n)/2)
    ! e_q^2 *F_q(A)*F_q(B)
    FA=SiversTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)+uTMDFF_inB(x2,b,mu,zeta2,3)
    Integrand=-global_mass_scale*(&
        (FA(1)+FA(2))*(FB(1)+4d0*FB(2))/18d0&
        +FA(3)*FB(3)/9.d0&
        +FA(4)*FB(4)*4d0/9.d0&
        +FA(5)*FB(5)/9d0&
        +(FA(-1)+FA(-2))*(FB(-1)+4d0*FB(-2))/18d0&
        +FA(-3)*FB(-3)/9.d0&
        +FA(-4)*FB(-4)*4d0/9.d0&
        +FA(-5)*FB(-5)/9d0)
    !--------------------------------------------------------------------------------  
    CASE (12105) !n->h? where h?=h1+h2 (n=neutron=p(u<->d))
    ! e_q^2 *F_q(A)*F_q(B)
    FA=SiversTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)
    Integrand=-global_mass_scale*(&
        FA(2)*FB(1)/9.d0&
        +FA(1)*FB(2)*4.d0/9.d0&
        +FA(3)*FB(3)/9.d0&
        +FA(4)*FB(4)*4d0/9.d0&
        +FA(5)*FB(5)/9d0&
        +FA(-2)*FB(-1)/9.d0&
        +FA(-1)*FB(-2)*4.d0/9.d0&
        +FA(-3)*FB(-3)/9.d0&
        +FA(-4)*FB(-4)*4d0/9.d0&
        +FA(-5)*FB(-5)/9d0)
    !--------------------------------------------------------------------------------  
    CASE (12106) !n->h? where h?=h1+h2+h3 (n=neutron=p(u<->d))
    ! e_q^2 *F_q(A)*F_q(B)
    FA=SiversTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)+uTMDFF_inB(x2,b,mu,zeta2,3)
    Integrand=-global_mass_scale*(&
        FA(2)*FB(1)/9.d0&
        +FA(1)*FB(2)*4.d0/9.d0&
        +FA(3)*FB(3)/9.d0&
        +FA(4)*FB(4)*4d0/9.d0&
        +FA(5)*FB(5)/9d0&
        +FA(-2)*FB(-1)/9.d0&
        +FA(-1)*FB(-2)*4.d0/9.d0&
        +FA(-3)*FB(-3)/9.d0&
        +FA(-4)*FB(-4)*4d0/9.d0&
        +FA(-5)*FB(-5)/9d0)
    !------------------------------------------------------------------------------------
    CASE (12111) !p->bar h? where h?=h1+h2
    ! e_q^2 *F_q(A)*F_bq(B)
    FA=SiversTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)
    Integrand=-global_mass_scale*(&
        FA(1)*FB(-1)/9.d0&
        +FA(2)*FB(-2)*4.d0/9.d0&
        +FA(3)*FB(-3)/9.d0&
        +FA(4)*FB(-4)*4d0/9.d0&
        +FA(5)*FB(-5)/9d0&
        +FA(-1)*FB(1)/9.d0&
        +FA(-2)*FB(2)*4.d0/9.d0&
        +FA(-3)*FB(3)/9.d0&
        +FA(-4)*FB(4)*4d0/9.d0&
        +FA(-5)*FB(5)/9d0)
    !--------------------------------------------------------------------------------  
    CASE (12112) !p->bar h? where h?=h1+h2+h3
    ! e_q^2 *F_q(A)*F_bq(B)
    FA=SiversTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)+uTMDFF_inB(x2,b,mu,zeta2,3)
    Integrand=-global_mass_scale*(&
        FA(1)*FB(-1)/9.d0&
        +FA(2)*FB(-2)*4.d0/9.d0&
        +FA(3)*FB(-3)/9.d0&
        +FA(4)*FB(-4)*4d0/9.d0&
        +FA(5)*FB(-5)/9d0&
        +FA(-1)*FB(1)/9.d0&
        +FA(-2)*FB(2)*4.d0/9.d0&
        +FA(-3)*FB(3)/9.d0&
        +FA(-4)*FB(4)*4d0/9.d0&
        +FA(-5)*FB(5)/9d0)
    !--------------------------------------------------------------------------------  
    CASE (12113) !d->bar h? where h?=h1+h2 (d=deutron=(p+n)/2)
    ! e_q^2 *F_q(A)*F_bq(B)
    FA=SiversTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)
    Integrand=-global_mass_scale*(&
        (FA(1)+FA(2))*(FB(-1)+4d0*FB(-2))/18d0&
        +FA(3)*FB(-3)/9.d0&
        +FA(4)*FB(-4)*4d0/9.d0&
        +FA(5)*FB(-5)/9d0&
        +(FA(-1)+FA(-2))*(FB(1)+4d0*FB(2))/18d0&
        +FA(-3)*FB(3)/9.d0&
        +FA(-4)*FB(4)*4d0/9.d0&
        +FA(-5)*FB(5)/9d0)
    !--------------------------------------------------------------------------------  
    CASE (12114) !d->bar h? where h?=h1+h2+h3 (d=deutron=(p+n)/2)
    ! e_q^2 *F_q(A)*F_bq(B)
    FA=SiversTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)+uTMDFF_inB(x2,b,mu,zeta2,3)
    Integrand=-global_mass_scale*(&
        (FA(1)+FA(2))*(FB(-1)+4d0*FB(-2))/18d0&
        +FA(3)*FB(-3)/9.d0&
        +FA(4)*FB(-4)*4d0/9.d0&
        +FA(5)*FB(-5)/9d0&
        +(FA(-1)+FA(-2))*(FB(1)+4d0*FB(2))/18d0&
        +FA(-3)*FB(3)/9.d0&
        +FA(-4)*FB(4)*4d0/9.d0&
        +FA(-5)*FB(5)/9d0)
    !------------------------------------------------------------------------------------
    CASE (12115) !n->bar h? where h?=h1+h2 (n=neutron=p(u<->d))
    ! e_q^2 *F_q(A)*F_bq(B)
    FA=SiversTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)
    Integrand=-global_mass_scale*(&
        FA(2)*FB(-1)/9.d0&
        +FA(1)*FB(-2)*4.d0/9.d0&
        +FA(3)*FB(-3)/9.d0&
        +FA(4)*FB(-4)*4d0/9.d0&
        +FA(5)*FB(-5)/9d0&
        +FA(-2)*FB(1)/9.d0&
        +FA(-1)*FB(2)*4.d0/9.d0&
        +FA(-3)*FB(3)/9.d0&
        +FA(-4)*FB(4)*4d0/9.d0&
        +FA(-5)*FB(5)/9d0)
    !--------------------------------------------------------------------------------  
    CASE (12116) !n->bar h? where h?=h1+h2+h3 (n=neutron=p(u<->d))
    ! e_q^2 *F_q(A)*F_bq(B)
    FA=SiversTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)+uTMDFF_inB(x2,b,mu,zeta2,3)
    Integrand=-global_mass_scale*(&
        FA(2)*FB(-1)/9.d0&
        +FA(1)*FB(-2)*4.d0/9.d0&
        +FA(3)*FB(-3)/9.d0&
        +FA(4)*FB(-4)*4d0/9.d0&
        +FA(5)*FB(-5)/9d0&
        +FA(-2)*FB(1)/9.d0&
        +FA(-1)*FB(2)*4.d0/9.d0&
        +FA(-3)*FB(3)/9.d0&
        +FA(-4)*FB(4)*4d0/9.d0&
        +FA(-5)*FB(5)/9d0)
    !--------------------------------------------------------------------------------
    !----------------------------------WGT-related asymetries------------------------
    !-------------------------------------------------------------------------------
    !---------------------------------SIDIS-----------------------------------------    
    CASE (13001:13009) !A_LT asymmetry p->hN where n=last number
    ! e_q^2 *F_q(A)*F_q(B)
    h=process-13000
    FA=wgtTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,h)
    Integrand=+global_mass_scale*(&
        FA(1)*FB(1)/9.d0&
        +FA(2)*FB(2)*4.d0/9.d0&
        +FA(3)*FB(3)/9.d0&
        +FA(4)*FB(4)*4d0/9.d0&
        +FA(5)*FB(5)/9d0&
        +FA(-1)*FB(-1)/9.d0&
        +FA(-2)*FB(-2)*4.d0/9.d0&
        +FA(-3)*FB(-3)/9.d0&
        +FA(-4)*FB(-4)*4d0/9.d0&
        +FA(-5)*FB(-5)/9d0)
    !--------------------------------------------------------------------------------  
    CASE (13011:13019) !A_LT asymmetry d->hN where n=last number (d=deutron=(p+n)/2)
    ! e_q^2 *F_q(A)*F_q(B)
    h=process-13010
    FA=wgtTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,h)
    Integrand=+global_mass_scale*(&
        (FA(1)+FA(2))*(FB(1)+4d0*FB(2))/18d0&
        +FA(3)*FB(3)/9.d0&
        +FA(4)*FB(4)*4d0/9.d0&
        +FA(5)*FB(5)/9d0&
        +(FA(-1)+FA(-2))*(FB(-1)+4d0*FB(-2))/18d0&
        +FA(-3)*FB(-3)/9.d0&
        +FA(-4)*FB(-4)*4d0/9.d0&
        +FA(-5)*FB(-5)/9d0)
    !--------------------------------------------------------------------------------  
    CASE (13021:13029) !A_LT asymmetry p->bar-hN where n=last number
    ! e_q^2 *F_q(A)*F_bar-q(B)
    h=process-13020
    FA=wgtTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,h)
    Integrand=+global_mass_scale*(&
        FA(1)*FB(-1)/9.d0&
        +FA(2)*FB(-2)*4.d0/9.d0&
        +FA(3)*FB(-3)/9.d0&
        +FA(4)*FB(-4)*4d0/9.d0&
        +FA(5)*FB(-5)/9d0&
        +FA(-1)*FB(1)/9.d0&
        +FA(-2)*FB(2)*4.d0/9.d0&
        +FA(-3)*FB(3)/9.d0&
        +FA(-4)*FB(4)*4d0/9.d0&
        +FA(-5)*FB(5)/9d0)
    !--------------------------------------------------------------------------------  
    CASE (13031:13039) !A_LT asymmetry d->bar-hN where n=last number (d=deutron=(p+n)/2)
    ! e_q^2 *F_q(A)*F_bar-q(B)
    h=process-13030
    FA=wgtTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,h)
    Integrand=+global_mass_scale*(&
        (FA(1)+FA(2))*(FB(-1)+4d0*FB(-2))/18d0&
        +FA(3)*FB(-3)/9.d0&
        +FA(4)*FB(-4)*4d0/9.d0&
        +FA(5)*FB(-5)/9d0&
        +(FA(-1)+FA(-2))*(FB(1)+4d0*FB(2))/18d0&
        +FA(-3)*FB(3)/9.d0&
        +FA(-4)*FB(4)*4d0/9.d0&
        +FA(-5)*FB(5)/9d0)
    !--------------------------------------------------------------------------------
    CASE (13041:13049) !A_LT asymmetry p->hN where n=last number  (n=neutron=p(u<->d))
    ! e_q^2 *F_q(A)*F_q(B)
    h=process-13040
    FA=wgtTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,h)
    Integrand=+global_mass_scale*(&
        FA(2)*FB(1)/9.d0&
        +FA(1)*FB(2)*4.d0/9.d0&
        +FA(3)*FB(3)/9.d0&
        +FA(4)*FB(4)*4d0/9.d0&
        +FA(5)*FB(5)/9d0&
        +FA(-2)*FB(-1)/9.d0&
        +FA(-1)*FB(-2)*4.d0/9.d0&
        +FA(-3)*FB(-3)/9.d0&
        +FA(-4)*FB(-4)*4d0/9.d0&
        +FA(-5)*FB(-5)/9d0)
    !--------------------------------------------------------------------------------  
    CASE (13051:13059) !A_LT asymmetry p->bar-hN where n=last number (n=neutron=p(u<->d))
    ! e_q^2 *F_q(A)*F_bar-q(B)
    h=process-13050
    FA=wgtTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,h)
    Integrand=+global_mass_scale*(&
        FA(2)*FB(-1)/9.d0&
        +FA(1)*FB(-2)*4.d0/9.d0&
        +FA(3)*FB(-3)/9.d0&
        +FA(4)*FB(-4)*4d0/9.d0&
        +FA(5)*FB(-5)/9d0&
        +FA(-2)*FB(1)/9.d0&
        +FA(-1)*FB(2)*4.d0/9.d0&
        +FA(-3)*FB(3)/9.d0&
        +FA(-4)*FB(4)*4d0/9.d0&
        +FA(-5)*FB(5)/9d0)
    !--------------------------------------------------------------------------------  
    !--------------------------------------------------------------------------------  
    CASE (13101) !p->h? where h?=h1+h2
    ! e_q^2 *F_q(A)*F_q(B)
    FA=wgtTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)
    Integrand=+global_mass_scale*(&
        FA(1)*FB(1)/9.d0&
        +FA(2)*FB(2)*4.d0/9.d0&
        +FA(3)*FB(3)/9.d0&
        +FA(4)*FB(4)*4d0/9.d0&
        +FA(5)*FB(5)/9d0&
        +FA(-1)*FB(-1)/9.d0&
        +FA(-2)*FB(-2)*4.d0/9.d0&
        +FA(-3)*FB(-3)/9.d0&
        +FA(-4)*FB(-4)*4d0/9.d0&
        +FA(-5)*FB(-5)/9d0)
    !--------------------------------------------------------------------------------  
    CASE (13102) !p->h? where h?=h1+h2+h3
    ! e_q^2 *F_q(A)*F_q(B)
    FA=wgtTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)+uTMDFF_inB(x2,b,mu,zeta2,3)
    Integrand=+global_mass_scale*(&
        FA(1)*FB(1)/9.d0&
        +FA(2)*FB(2)*4.d0/9.d0&
        +FA(3)*FB(3)/9.d0&
        +FA(4)*FB(4)*4d0/9.d0&
        +FA(5)*FB(5)/9d0&
        +FA(-1)*FB(-1)/9.d0&
        +FA(-2)*FB(-2)*4.d0/9.d0&
        +FA(-3)*FB(-3)/9.d0&
        +FA(-4)*FB(-4)*4d0/9.d0&
        +FA(-5)*FB(-5)/9d0)
    !--------------------------------------------------------------------------------  
    CASE (13103) !d->h? where h?=h1+h2 (d=deutron=(p+n)/2)
    ! e_q^2 *F_q(A)*F_q(B)
    FA=wgtTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)
    Integrand=+global_mass_scale*(&
        (FA(1)+FA(2))*(FB(1)+4d0*FB(2))/18d0&
        +FA(3)*FB(3)/9.d0&
        +FA(4)*FB(4)*4d0/9.d0&
        +FA(5)*FB(5)/9d0&
        +(FA(-1)+FA(-2))*(FB(-1)+4d0*FB(-2))/18d0&
        +FA(-3)*FB(-3)/9.d0&
        +FA(-4)*FB(-4)*4d0/9.d0&
        +FA(-5)*FB(-5)/9d0)
    !--------------------------------------------------------------------------------  
    CASE (13104) !d->h? where h?=h1+h2+h3 (d=deutron=(p+n)/2)
    ! e_q^2 *F_q(A)*F_q(B)
    FA=wgtTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)+uTMDFF_inB(x2,b,mu,zeta2,3)
    Integrand=+global_mass_scale*(&
        (FA(1)+FA(2))*(FB(1)+4d0*FB(2))/18d0&
        +FA(3)*FB(3)/9.d0&
        +FA(4)*FB(4)*4d0/9.d0&
        +FA(5)*FB(5)/9d0&
        +(FA(-1)+FA(-2))*(FB(-1)+4d0*FB(-2))/18d0&
        +FA(-3)*FB(-3)/9.d0&
        +FA(-4)*FB(-4)*4d0/9.d0&
        +FA(-5)*FB(-5)/9d0)
    !--------------------------------------------------------------------------------  
    CASE (13105) !n->h? where h?=h1+h2 (n=neutron=p(u<->d))
    ! e_q^2 *F_q(A)*F_q(B)
    FA=wgtTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)
    Integrand=+global_mass_scale*(&
        FA(2)*FB(1)/9.d0&
        +FA(1)*FB(2)*4.d0/9.d0&
        +FA(3)*FB(3)/9.d0&
        +FA(4)*FB(4)*4d0/9.d0&
        +FA(5)*FB(5)/9d0&
        +FA(-2)*FB(-1)/9.d0&
        +FA(-1)*FB(-2)*4.d0/9.d0&
        +FA(-3)*FB(-3)/9.d0&
        +FA(-4)*FB(-4)*4d0/9.d0&
        +FA(-5)*FB(-5)/9d0)
    !--------------------------------------------------------------------------------  
    CASE (13106) !n->h? where h?=h1+h2+h3 (n=neutron=p(u<->d))
    ! e_q^2 *F_q(A)*F_q(B)
    FA=wgtTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)+uTMDFF_inB(x2,b,mu,zeta2,3)
    Integrand=+global_mass_scale*(&
        FA(2)*FB(1)/9.d0&
        +FA(1)*FB(2)*4.d0/9.d0&
        +FA(3)*FB(3)/9.d0&
        +FA(4)*FB(4)*4d0/9.d0&
        +FA(5)*FB(5)/9d0&
        +FA(-2)*FB(-1)/9.d0&
        +FA(-1)*FB(-2)*4.d0/9.d0&
        +FA(-3)*FB(-3)/9.d0&
        +FA(-4)*FB(-4)*4d0/9.d0&
        +FA(-5)*FB(-5)/9d0)
    !------------------------------------------------------------------------------------
    CASE (13111) !p->bar h? where h?=h1+h2
    ! e_q^2 *F_q(A)*F_bq(B)
    FA=wgtTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)
    Integrand=+global_mass_scale*(&
        FA(1)*FB(-1)/9.d0&
        +FA(2)*FB(-2)*4.d0/9.d0&
        +FA(3)*FB(-3)/9.d0&
        +FA(4)*FB(-4)*4d0/9.d0&
        +FA(5)*FB(-5)/9d0&
        +FA(-1)*FB(1)/9.d0&
        +FA(-2)*FB(2)*4.d0/9.d0&
        +FA(-3)*FB(3)/9.d0&
        +FA(-4)*FB(4)*4d0/9.d0&
        +FA(-5)*FB(5)/9d0)
    !--------------------------------------------------------------------------------  
    CASE (13112) !p->bar h? where h?=h1+h2+h3
    ! e_q^2 *F_q(A)*F_bq(B)
    FA=wgtTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)+uTMDFF_inB(x2,b,mu,zeta2,3)
    Integrand=+global_mass_scale*(&
        FA(1)*FB(-1)/9.d0&
        +FA(2)*FB(-2)*4.d0/9.d0&
        +FA(3)*FB(-3)/9.d0&
        +FA(4)*FB(-4)*4d0/9.d0&
        +FA(5)*FB(-5)/9d0&
        +FA(-1)*FB(1)/9.d0&
        +FA(-2)*FB(2)*4.d0/9.d0&
        +FA(-3)*FB(3)/9.d0&
        +FA(-4)*FB(4)*4d0/9.d0&
        +FA(-5)*FB(5)/9d0)
    !--------------------------------------------------------------------------------  
    CASE (13113) !d->bar h? where h?=h1+h2 (d=deutron=(p+n)/2)
    ! e_q^2 *F_q(A)*F_bq(B)
    FA=wgtTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)
    Integrand=+global_mass_scale*(&
        (FA(1)+FA(2))*(FB(-1)+4d0*FB(-2))/18d0&
        +FA(3)*FB(-3)/9.d0&
        +FA(4)*FB(-4)*4d0/9.d0&
        +FA(5)*FB(-5)/9d0&
        +(FA(-1)+FA(-2))*(FB(1)+4d0*FB(2))/18d0&
        +FA(-3)*FB(3)/9.d0&
        +FA(-4)*FB(4)*4d0/9.d0&
        +FA(-5)*FB(5)/9d0)
    !--------------------------------------------------------------------------------  
    CASE (13114) !d->bar h? where h?=h1+h2+h3 (d=deutron=(p+n)/2)
    ! e_q^2 *F_q(A)*F_bq(B)
    FA=wgtTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)+uTMDFF_inB(x2,b,mu,zeta2,3)
    Integrand=+global_mass_scale*(&
        (FA(1)+FA(2))*(FB(-1)+4d0*FB(-2))/18d0&
        +FA(3)*FB(-3)/9.d0&
        +FA(4)*FB(-4)*4d0/9.d0&
        +FA(5)*FB(-5)/9d0&
        +(FA(-1)+FA(-2))*(FB(1)+4d0*FB(2))/18d0&
        +FA(-3)*FB(3)/9.d0&
        +FA(-4)*FB(4)*4d0/9.d0&
        +FA(-5)*FB(5)/9d0)
    !------------------------------------------------------------------------------------
    CASE (13115) !n->bar h? where h?=h1+h2 (n=neutron=p(u<->d))
    ! e_q^2 *F_q(A)*F_bq(B)
    FA=wgtTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)
    Integrand=+global_mass_scale*(&
        FA(2)*FB(-1)/9.d0&
        +FA(1)*FB(-2)*4.d0/9.d0&
        +FA(3)*FB(-3)/9.d0&
        +FA(4)*FB(-4)*4d0/9.d0&
        +FA(5)*FB(-5)/9d0&
        +FA(-2)*FB(1)/9.d0&
        +FA(-1)*FB(2)*4.d0/9.d0&
        +FA(-3)*FB(3)/9.d0&
        +FA(-4)*FB(4)*4d0/9.d0&
        +FA(-5)*FB(5)/9d0)
    !--------------------------------------------------------------------------------  
    CASE (13116) !n->bar h? where h?=h1+h2+h3 (n=neutron=p(u<->d))
    ! e_q^2 *F_q(A)*F_bq(B)
    FA=wgtTMDPDF_inB(x1,b,mu,zeta1,1)
    FB=uTMDFF_inB(x2,b,mu,zeta2,1)+uTMDFF_inB(x2,b,mu,zeta2,2)+uTMDFF_inB(x2,b,mu,zeta2,3)
    Integrand=+global_mass_scale*(&
        FA(2)*FB(-1)/9.d0&
        +FA(1)*FB(-2)*4.d0/9.d0&
        +FA(3)*FB(-3)/9.d0&
        +FA(4)*FB(-4)*4d0/9.d0&
        +FA(5)*FB(-5)/9d0&
        +FA(-2)*FB(1)/9.d0&
        +FA(-1)*FB(2)*4.d0/9.d0&
        +FA(-3)*FB(3)/9.d0&
        +FA(-4)*FB(4)*4d0/9.d0&
        +FA(-5)*FB(5)/9d0)
    !--------------------------------------------------------------------------------
    !-------------------------SPECIAL DY CASES---------------------------------------  
    CASE (13200) !pp->Z+gamma        
        FA=wgtTMDPDF_inB(x1,b,mu,zeta1,h1)
        FB=uTMDPDF_inB(x2,b,mu,zeta2,-h2)!!! -h2, to multiply quarks by anti-quarks in FAB
        FAB=FA*FB
            
        Integrand=global_mass_scale*XIntegrandForDYwithZgamma_GTU(FAB,Q2)
    !--------------------------------------------------------------------------------  
    CASE (13201) !pp-> W+
        FA=wgtTMDPDF_inB(x1,b,mu,zeta1,h1)
        FB=uTMDPDF_inB(x2,b,mu,zeta2,h2)
        
        Integrand=-global_mass_scale*paramW_L*(& !!! -1=is due to the -gL^2 in the coupling for lepton
        paramW_UD*(FA(2)*FB(-1)-FA(-1)*FB(2))&        !u*dbar+dbar*u
        +paramW_US*(FA(2)*FB(-3)-FA(-3)*FB(2))&        !u*sbar+sbar*u
        +paramW_UB*(FA(2)*FB(-5)-FA(-5)*FB(2))&        !u*bbar+bbar*u
        +paramW_CD*(FA(4)*FB(-1)-FA(-1)*FB(4))&        !c*dbar+dbar*c
        +paramW_CS*(FA(4)*FB(-3)-FA(-3)*FB(4))&        !c*sbar+sbar*c
        +paramW_CB*(FA(4)*FB(-5)-FA(-5)*FB(4))&        !c*bbar+bbar*c
        )*Q2*Q2/((Q2-MW2)**2+GammaW2*MW2)
    !--------------------------------------------------------------------------------  
    CASE (13202) !pp-> W-
        FA=wgtTMDPDF_inB(x1,b,mu,zeta1,h1)
        FB=uTMDPDF_inB(x2,b,mu,zeta2,h2)
        
        Integrand=-global_mass_scale*paramW_L*(& !!! -1=is due to the -gL^2 in the coupling for lepton
        paramW_UD*(FA(1)*FB(-2)-FA(-2)*FB(1))&        !d*ubar+ubar*d
        +paramW_US*(FA(3)*FB(-2)-FA(-2)*FB(3))&        !s*ubar+ubar*s
        +paramW_UB*(FA(5)*FB(-2)-FA(-2)*FB(5))&        !b*ubar+ubar*b
        +paramW_CD*(FA(1)*FB(-4)-FA(-4)*FB(1))&        !d*cbar+cbar*d
        +paramW_CS*(FA(3)*FB(-4)-FA(-4)*FB(3))&        !s*cbar+cbar*s
        +paramW_CB*(FA(5)*FB(-4)-FA(-4)*FB(5))&        !b*cbar+cbar*b
        )*Q2*Q2/((Q2-MW2)**2+GammaW2*MW2) 
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

!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!!! The hadron tensonr for the structure function GTU in DY icludes Z + gamma
function XIntegrandForDYwithZgamma_GTU(FAB,Q2)
     real(dp)::XIntegrandForDYwithZgamma_GTU,Q2
    !!cross-seciton parameters
     real(dp),dimension(-5:5):: FAB     
     
     !gamma-part=0
    XIntegrandForDYwithZgamma_GTU=&!gamma-Z interference
      paramMIXL_A*(&
      paramMIXU*FAB(2)&
      +paramMIXD*FAB(1)&
      +paramMIXS*FAB(3)&
      +paramMIXC*FAB(4)&
      +paramMIXB*FAB(5)&
      -paramMIXU*FAB(-2)&
      -paramMIXD*FAB(-1)&
      -paramMIXS*FAB(-3)&
      -paramMIXC*FAB(-4)&
      -paramMIXB*FAB(-5))*&
      2d0*Q2*(Q2-MZ2)/((Q2-MZ2)**2+GammaZ2*MZ2)&
     +&!ZZ-contributions
     paramL_A*(&
      paramU*FAB(2)&
      +paramD*FAB(1)&
      +paramS*FAB(3)&
      +paramC*FAB(4)&
      +paramB*FAB(5)&
      -paramU*FAB(-2)&
      -paramD*FAB(-1)&
      -paramS*FAB(-3)&
      -paramC*FAB(-4)&
      -paramB*FAB(-5))*&
      Q2*Q2/((Q2-MZ2)**2+GammaZ2*MZ2)
     
end function XIntegrandForDYwithZgamma_GTU

!!! Combination Delta^{GG'} z_{-l}z_{-f} {FF}_A
function XTMD_pairZmZm_A(FAB,Q2)
     real(dp)::XTMD_pairZmZm_A
     real(dp),intent(in)::Q2
    !!cross-seciton parameters
     real(dp),dimension(-5:5),intent(in):: FAB

     XTMD_pairZmZm_A=&  !zM_gg =0
     zM_gZ_L*(& !gamma-Z interference
      zM_gZ_U*FAB(2)&
      +zM_gZ_D*FAB(1)&
      +zM_gZ_S*FAB(3)&
      +zM_gZ_C*FAB(4)&
      +zM_gZ_B*FAB(5)&
      -zM_gZ_U*FAB(-2)&
      -zM_gZ_D*FAB(-1)&
      -zM_gZ_S*FAB(-3)&
      -zM_gZ_C*FAB(-4)&
      -zM_gZ_B*FAB(-5))*&
      2d0*Q2*(Q2-MZ2)/((Q2-MZ2)**2+GammaZ2*MZ2)&
     +&!ZZ-contributions
       zM_ZZ_L*(&
       zM_ZZ_U*FAB(2)&
      +zM_ZZ_D*FAB(1)&
      +zM_ZZ_S*FAB(3)&
      +zM_ZZ_C*FAB(4)&
      +zM_ZZ_B*FAB(5)&
      -zM_ZZ_U*FAB(-2)&
      -zM_ZZ_D*FAB(-1)&
      -zM_ZZ_S*FAB(-3)&
      -zM_ZZ_C*FAB(-4)&
      -zM_ZZ_B*FAB(-5))*&
      Q2*Q2/((Q2-MZ2)**2+GammaZ2*MZ2)

end function XTMD_pairZmZm_A

end module TMDF
