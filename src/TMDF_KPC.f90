!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 3.00
!
!    Evaluation of the convolution integral for DY TMD cross-section with KPC
!
!    if you use this module please, quote 2307.13054
!
!    ver 3.0: created (AV, 07.09.2023)
!    ver 3.03: optimization of integration algorithms (Sara&Oscar, 14.01.2026)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDF_KPC
use aTMDe_Numerics
use aTMDe_Integration
use aTMDe_IO
use EWinput
use uTMDPDF
use uTMDFF
use lpTMDPDF
use SiversTMDPDF
use wgtTMDPDF
use BoerMuldersTMDPDF
use wglTMDPDF
use CollinsTMDFF

implicit none

private

!!!!!! 1=accurate but slow
!!!!!! 2=fast but not accurate
#define INTEGRATION_MODE 2

character (len=8),parameter :: moduleName="TMDF-KPC"
character (len=5),parameter :: version="v3.03"
!Last appropriate verion of constants-file
integer,parameter::inputver=37

integer::outputLevel=2
type(Warning_OBJ)::Warning_Handler

logical::started=.false.
!! flag for loss of convergence
logical:: convergenceLost=.false.

logical::include_uTMDPDF
logical::include_uTMDFF
logical::include_lpTMDPDF
logical::include_SiversTMDPDF
logical::include_wgtTMDPDF
logical::include_BoerMuldersTMDPDF
logical::include_wglTMDPDF
logical::include_CollinsTMDFF

!------------------------------------------Working variables------------------------------------------------------------
!! tolerances for integration and general
real(dp)::toleranceGEN
real(dp)::toleranceINT
real(dp)::qTMIN

!increment counters
integer::GlobalCounter=0 !!!total counter of calls of TMD pairs
integer::LocalCounter=0 !!!counter of calls of TMD pairs within the current integrand

! parameter of TMD proportionality
real(dp)::M2=1._dp


!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Public declarations
public::TMDF_KPC_IsInitialized,TMDF_KPC_Initialize,TMDF_KPC_IsconvergenceLost
public::KPC_DYconv,KPC_SIDISconv

contains

INCLUDE 'Code/TMDF_KPC/TMDpairs.f90'
INCLUDE 'Code/TMDF_KPC/KERNELpairs_DY.f90'
INCLUDE 'Code/TMDF_KPC/KERNELpairs_SIDIS.f90'

function TMDF_KPC_IsInitialized()
    logical::TMDF_KPC_IsInitialized
    TMDF_KPC_IsInitialized=started
end function TMDF_KPC_IsInitialized

   !! Initialization of the package
subroutine TMDF_KPC_Initialize(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=300)::path
    logical::initRequired
    integer::FILEver,messageTrigger
    real(dp)::dummy

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
        CLOSE (51, STATUS='KEEP')
        stop
    end if
    call MoveTO(51,'*p2  ')
    read(51,*) outputLevel
    if(outputLevel>2) write(*,*) '--------------------------------------------- '
    if(outputLevel>2) write(*,*) 'artemide.',moduleName,version,': initialization started ... '
    call MoveTO(51,'*p3  ')
    read(51,*) messageTrigger

    !!! mass parameter
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p2  ')
    read(51,*) dummy
    M2=dummy**2

    call MoveTO(51,'*15   ')
    call MoveTO(51,'*p1  ')
    read(51,*) initRequired
    if(.not.initRequired) then
        if(outputLevel>2) write(*,*)'artemide.',moduleName,': initialization is not required. '
        started=.false.
        CLOSE (51, STATUS='KEEP')
        return
    end if
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) toleranceGEN
    call MoveTO(51,'*p2  ')
    read(51,*) toleranceINT
    call MoveTO(51,'*p3  ')
    read(51,*) qTMIN

    if(outputLevel>2) write(*,'(A,ES8.2)') ' | tolerance general    : ',toleranceGEN
    if(outputLevel>2) write(*,'(A,ES8.2)') ' | tolerance integral   : ',toleranceINT

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

    !! wglTMDPDF
    call MoveTO(51,'*16  ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_wglTMDPDF

    !! BoerMuldersTMDPDF
    call MoveTO(51,'*18  ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_CollinsTMDFF

    CLOSE (51, STATUS='KEEP')
    Warning_Handler=Warning_OBJ(moduleName=moduleName,messageCounter=0,messageTrigger=messageTrigger)

    convergenceLost=.false.
    GlobalCounter=0
    LocalCounter=0

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
        if(outputLevel>1) write(*,*) '.. initializing BoerMuldersTMDPDF (from ',moduleName,')'
        if(present(prefix)) then
            call BoerMuldersTMDPDF_Initialize(file,prefix)
        else
            call BoerMuldersTMDPDF_Initialize(file)
        end if
    end if

    if(include_wglTMDPDF .and. (.not.wglTMDPDF_IsInitialized())) then
    if(outputLevel>1) write(*,*) '.. initializing wglTMDPDF (from ',moduleName,')'
    if(present(prefix)) then
        call wglTMDPDF_Initialize(file,prefix)
    else
        call wglTMDPDF_Initialize(file)
    end if
    end if

    if(include_CollinsTMDFF .and. (.not.CollinsTMDFF_IsInitialized())) then
    if(outputLevel>1) write(*,*) '.. initializing CollinsTMDFF (from ',moduleName,')'
    if(present(prefix)) then
        call CollinsTMDFF_Initialize(file,prefix)
    else
        call CollinsTMDFF_Initialize(file)
    end if
    end if

    started=.true.

#if INTEGRATION_MODE==2
    write(*,*)  color('--------------------------------------------------------',c_red)
    write(*,*)  color('----------------------  WARNING!  ----------------------',c_red)
    write(*,*)  color('--   TMDF_KPC is in the approximate integration mode  --',c_red)
    write(*,*)  color('--            Faster, but lower precision.            --',c_red)
    write(*,*)  color('--    Switch to default version by changing flag      --',c_red)
    write(*,*)  color('--   INTEGRATION_MODE in TMDF_KPC.f90, and recompile  --',c_red)
    write(*,*)  color('--------------------------------------------------------',c_red)
#endif

    if(outputLevel>0) write(*,*) color('----- arTeMiDe.TMDF '//trim(version)//': .... initialized',c_green)
    if(outputLevel>1) write(*,*) ' '


end subroutine TMDF_KPC_Initialize


!!!--------------------------------------------------------------------------------------------------
!!!!!!!Functions which carry the trigger on convergences.... Its used in xSec, and probably in other places.
function TMDF_KPC_IsconvergenceLost()
    logical::TMDF_KPC_IsconvergenceLost
    !!! checks TMDs trigger
    TMDF_KPC_IsconvergenceLost=convergenceLost
end function TMDF_KPC_IsconvergenceLost

subroutine TMDF_KPC_convergenceISlost()
    convergenceLost=.true.
    if(outputLevel>1) call Warning_Handler%WarningRaise('convergenceLOST trigger ON')
end subroutine TMDF_KPC_convergenceISlost

!!!--------------------------------------------------------------------------------------------------
!!!-----------------------------------DY PART -------------------------------------------------------
!!!--------------------------------------------------------------------------------------------------


!!!--------------------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Function that computes the integral for KPC convolution in DY
!!! proc1 = (int,int,int) is the process def for TMD*TMD, and for the integral kernel
!!! Q2, qT, x1,x2, mu are usual DY variables
!!! NOTE: that Q2 is DY kinematic variable, and mu is the factorization scale
!!! THIS IS A SYMMETRIC VERSION (i.e. it should contain only cos(theta)
!!!-----
!!! The integral is 2D, over theta and alpha (which are complicated combinations)
!!! First evaluate over theta (0,pi), then over alpha (0,pi/2)
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!! ATTENTION: The integral has been updated in order to fit the KPCs.
!!! When we started working on the fit of DY data with KPCs, we realized that the integration over alpha
!!! was very slow because the integrand was not particularly smooth and featured a peak that, depending
!!! on the values of Q and qT, could become quite narrow. To address this, we implemented a change of
!!! variables (described below) and began working with a new variable, omega, which ranges from –1 to 1.
!!! In its final form, the integral is evaluated first over theta (0,pi) and then over omega (–1,1),
!!! using the K15 integration method for both integrations.
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!! NOTE: For computing observables that require this convolution integral, we recommend using
!!! the previous version with the integration over alpha, or alternatively the new version but with an
!!! adaptive method such as GK or SA, since using K15 together with the new change of variables
!!! sacrifices a bit of precision (though not too much) in order to prioritize speed.
function KPC_DYconv(Q2,qT_in,x1,x2,mu,proc1)
    real(dp),intent(in)::Q2,qT_in,x1,x2,mu
    integer,intent(in),dimension(1:3)::proc1
    real(dp)::KPC_DYconv

    real(dp)::tau2,deltaT,qT
    logical::exist

    if(qT_in<qTMIN) then
        qT=qTMIN
    else
        qT=qT_in
    end if

    LocalCounter=0

    tau2=Q2+qT**2
    deltaT=qT**2/tau2

#if INTEGRATION_MODE==1
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!!   Addaptive implementations of the integral for computing observables
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! First version, integral over theta (0,pi), then over alpha (0,pi/2)
    ! KPC_DYconv=Integrate_SA(Integrand_forAlpha,0._dp,piHalf,toleranceINT)

    ! Second version, integral over theta (0,pi), then over omega (-1,1)
    KPC_DYconv=Integrate_SA(Integrand_forOmega,-1._dp,1._dp,toleranceINT)

#elif INTEGRATION_MODE==2
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!!   Fixed-number of points implementations of the integral for fitting KPCs
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KPC_DYconv=Integrate_K15(Integrand_forOmega,-1._dp,1._dp)
#endif


contains

!!!! This is integral over Theta over integral over Alpha
function Integrand_forTheta(theta)
    real(dp)::Integrand_forTheta
    real(dp),intent(in)::theta
    real(dp)::cT

    cT=cos(theta)

    Integrand_forTheta=INT_overALPHA(Q2,tau2,deltaT,x1,x2,mu,proc1,cT)

end function Integrand_forTheta

!!!! This is integral over Alpha over integral over Theta
function Integrand_forAlpha(alpha)
    real(dp)::Integrand_forAlpha
    real(dp),intent(in)::alpha
    real(dp)::sA

    sA=sin(alpha)

    Integrand_forAlpha=INT_overTHETA(Q2,tau2,deltaT,x1,x2,mu,proc1,sA)
    !write(*,"('{',F16.12,',',F16.8,'},')") alpha,Integrand_forAlpha

end function Integrand_forAlpha

!!! Integral over omega with a special change of variables that resolves the peak observed in alpha.
!!! NOTE: The implemented change of variables is the following:
!!! omega = arctan[rho*pi*(alpha-delta)/(alpha*(pi-4*delta)+pi*delta)]/arctan[rho], where delta is the
!!! position of the peak and rho is its width.
!!! In addition, omega lies in [-1,1], and satisfies that omega = 0 <--> alpha = delta, meaning that,
!!! after the change of variables the peak ends up roughly centered.
!!! ATTENTION: The position and the width of the peak are, respectively and approximately, delta = qT/Q and rho = qT/Q,
!!! but depending on the relation between Q and qT, one of them needs to be slightly corrected.
!!! After running many tests, we identified two different cases (see below in the code).
function Integrand_forOmega(omega)
    real(dp) :: Integrand_forOmega
    real(dp), intent(in) :: omega
    real(dp) :: delta,rho,atan_rho,tan_o,jacobian
    real(dp) :: alpha,sA

    ! First case: Q is much greater than qT -> the position of the peak needs to be corrected
    if (Sqrt(Q2)/qT > 39.d0) then

        delta=qT/Sqrt(Q2)+0.5d0*(Sqrt(Q2)/qT)**0.35d0/Sqrt(Q2)
        rho=qT/Sqrt(Q2)

    ! Second case: Q is greater than qT -> the width of the peak needs to be corrected
    else

        delta=qT/Sqrt(Q2)
        rho=qT/Sqrt(Q2)+qT/(Sqrt(Q2)/qT)**(0.18d0)/2.5d0

    end if

    ! Some auxiliary variables to shorten expressions
    atan_rho=atan(rho)
    tan_o=tan(omega*atan_rho)

    ! Alpha as a function of omega
    alpha=pi*delta*(rho+tan_o)/(pi*rho+(4*delta-pi)*tan_o)
    sA=sin(alpha)

    ! Jacobian
    jacobian=2*pi*(pi-2*delta)*delta*rho*atan_rho*(1+tan_o**2)/(pi*rho+(4*delta-pi)*tan_o)**2

    Integrand_forOmega = jacobian*INT_overTHETA(Q2,tau2,deltaT,x1,x2,mu,proc1,sA)

end function Integrand_forOmega

end function KPC_DYconv


!!! the integral over alpha at given theta
function INT_overALPHA(Q2,tau2,deltaT,x1,x2,mu,proc1,cT)
    real(dp),intent(in)::Q2,tau2,deltaT,x1,x2,mu,cT
    integer,intent(in),dimension(1:3)::proc1
    real(dp)::INT_overALPHA

    INT_overALPHA=Integrate_SA(Integrand_forALPHA,0._dp,piHalf,toleranceINT)

contains

function Integrand_forALPHA(alpha)
    real(dp)::Integrand_forALPHA
    real(dp),intent(in)::alpha
    real(dp)::S,Lam,xi1,xi2,K1,K2,sinA

    sinA=sin(alpha)
    S=Sqrt(deltaT)*sinA*cT
    Lam=(1-deltaT)*(1-sinA*sinA)

    xi1=x1/2*(1+S+sqrt(Lam))
    xi2=x2/2*(1-S+sqrt(Lam))
    K1=tau2/4*((1+S)**2-Lam)
    K2=tau2/4*((1-S)**2-Lam)

    !!!! Some times K1 and K2 became too close to zero... and turns negative
    if(K1<toleranceGEN) K1=toleranceGEN
    if(K2<toleranceGEN) K2=toleranceGEN

    LocalCounter=LocalCounter+1

    !!! it is devided by 2 (instead of 4), because the integral over cos(theta) is over (0,pi).
    Integrand_forALPHA=TMD_pair(Q2,xi1,xi2,K1,K2,mu,proc1)*DY_KERNEL(Q2,tau2,tau2-Q2,S,Lam,sinA,cT,proc1(3))*sinA
end function Integrand_forALPHA

end function INT_overALPHA

!!! the integral over theta at given alpha
function INT_overTHETA(Q2,tau2,deltaT,x1,x2,mu,proc1,sA)
    real(dp),intent(in)::Q2,tau2,deltaT,x1,x2,mu,sA
    integer,intent(in),dimension(1:3)::proc1
    real(dp)::INT_overTHETA


#if INTEGRATION_MODE==1
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!!   Adaptive/Recommended for computing observables
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INT_overTHETA=Integrate_SA(Integrand_forTHETA,0._dp,pi,toleranceINT)
    !INT_overTHETA=Integrate_GK(Integrand_forTHETA,0._dp,pi,toleranceINT)

#elif INTEGRATION_MODE==2
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!!   Recommended implementations of the integral for fitting KPCs
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INT_overTHETA=Integrate_K15(Integrand_forTHETA,0._dp,pi)
#endif

contains

function Integrand_forTHETA(theta)
    real(dp)::Integrand_forTHETA
    real(dp),intent(in)::theta
    real(dp)::S,Lam,xi1,xi2,K1,K2,cosT

    cosT=cos(theta)
    S=Sqrt(deltaT)*sA*cosT
    Lam=(1-deltaT)*(1-sA*sA)

    xi1=x1/2*(1+S+sqrt(Lam))
    xi2=x2/2*(1-S+sqrt(Lam))
    K1=tau2/4*((1+S)**2-Lam)
    K2=tau2/4*((1-S)**2-Lam)

    !!!! Some times K1 and K2 became too close to zero... and turns negative
    if(K1<toleranceGEN) K1=toleranceGEN
    if(K2<toleranceGEN) K2=toleranceGEN

    LocalCounter=LocalCounter+1

    !!! it is divided by 2 (instead of 4), because the integral over cos(theta) is over (0,pi).
    Integrand_forTHETA=TMD_pair(Q2,xi1,xi2,K1,K2,mu,proc1)*DY_KERNEL(Q2,tau2,tau2-Q2,S,Lam,sA,cosT,proc1(3))*sA

end function Integrand_forTHETA

end function INT_overTHETA


!!!--------------------------------------------------------------------------------------------------
!!!----------------------------------- SIDIS PART ---------------------------------------------------
!!!--------------------------------------------------------------------------------------------------

!!!--------------------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Function that computes the integral for KPC convolution in SIDIS
!!! proc1 = (int,int,int) is the process def for TMD*TMD, and for the integral kernel
!!! Q2, qT, x1, z1, vareps and mu are the usual SIDIS variables
!!! NOTE: that Q2 is SIDIS kinematic variable, and mu is the factorization scale

!!! The integral is 2D, over theta and delta (which are complicated combinations)
!!! First evaluate over theta (0,pi), then over delta (0,1)

!!! ATTENTION: The integral has been updated in order to fit the KPCs.
!!! We realized that the integration over delta was very slow because the integrand
!!! was not particularly smooth and featured a peak that, depending on the values of Q and qT,
!!! could become quite narrow. To address this, we implemented a change of
!!! variables (described below) and began working with a new variable, omega, which ranges from –1 to 1.
!!! In its final form, the integral is evaluated first over theta (0,pi) and then over omega (–1,1),
!!! using the K21 integration method for both integrations.
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!! NOTE: For computing observables that require this convolution integral, we recommend using
!!! the previous version with the integration over delta, or alternatively the new version but with an
!!! adaptive method such as GK, since using K21 together with the new change of variables
!!! sacrifices a bit of precision (though not too much) in order to prioritize speed.

function KPC_SIDISconv(Q2, qT_in, x1, z1, vareps, mu, proc1)
    real(dp) :: KPC_SIDISconv
    real(dp), intent(in) :: Q2, x1, z1, mu, qT_in, vareps
    integer, intent(in), dimension(1:3) :: proc1
    real(dp) :: tau2, dT, qT, x

    if(qT_in < qTMIN) then
        qT = qTMIN
    else
        qT = qT_in
    end if

    LocalCounter = 0

    tau2 = Q2 - qT**2
    dT = qT**2/tau2
    x = x1*Q2/tau2

#if INTEGRATION_MODE==1
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!!   Adaptive implementations of the integral for computing observables
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    ! First version, integral over theta (0,pi), then over delta (0,1)
    KPC_SIDISconv = Integrate_GK(Integrand_forDeltat2, 0._dp, 1._dp, toleranceINT)

    ! Second version, integral over theta (0,pi), then over omega (-1,1)
!     KPC_SIDISconv=Integrate_GK(Integrand_forDelta_th,-1._dp,1._dp,toleranceINT)

#elif INTEGRATION_MODE==2
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!!   Fixed-number of points implementations of the integral for fitting KPCs
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    KPC_SIDISconv=Integrate_K21(Integrand_forDelta_th,-1._dp,1._dp)

#endif


contains

! Integral over Deltat over integral over theta
function Integrand_forDeltat2(Deltat)
    real(dp) :: Integrand_forDeltat2
    real(dp), intent(in) :: Deltat

    Integrand_forDeltat2 = INT_overTHETA_SIDIS(Q2, tau2, dT, x1, z1, vareps, mu, proc1, Deltat)

end function Integrand_forDeltat2


!!! Integral over omega with a special change of variables that resolves the peak observed in delta.
!!! NOTE: The implemented change of variables is the following:
!!! omega = arctan[a (delta-d)/(delta(1-2d)+d)]/arctan[a], where d is the position of the peak and a is its width.
!!! In addition, omega lies in [-1,1], and satisfies that omega = 0 <--> delta = d, meaning that,
!!! after the change of variables the peak ends up roughly centered.
!!! ATTENTION: The position and the width of the peak are, respectively and approximately, d = qT/Q and a = qT/Q,
!!! but depending on the vaue of x, one of them needs to be slightly corrected.
!!! After running many tests, we identified two different cases (see below in the code).

function Integrand_forDelta_th(omega)
    real(dp) :: Integrand_forDelta_th
    real(dp), intent(in) :: omega
    real(dp) :: d,a,delta,atan_a,tan_o,jacobian

    ! First case: "high" x -> the position of the peak needs to be corrected
    if (x > 0.2d0) then

        d=qT/Sqrt(Q2)+x*0.02d0*(Sqrt(Q2)/qT)**(0.3d0)
        a=Sqrt(Q2)/qT
!
    ! Second case: "low" x -> the width of the peak needs to be corrected
    else

        d=qT/Sqrt(Q2)
        a=qT/Sqrt(Q2)+5.1d0/(Sqrt(Q2)/qT)**(0.65d0)/x

    end if

    atan_a=atan(a)
    tan_o=tan(omega*atan_a)

    !!!! expression for delta as function of omega
    delta=d*(a+tan_o)/(a+(2*d-1)*tan_o)
    !!!
    jacobian=2*a*(1-d)*d*atan_a*(1+tan_o**2)/(a+(2*d-1)*tan_o)**2

    if(ISNAN(delta)) then
!     write(*,*) "---->",d,a,atan_a,tan_o,delta,jacobian

    end if

    !Deltat = sqrt(Deltat2)
    Integrand_forDelta_th = jacobian*INT_overTHETA_SIDIS(Q2, tau2, dT, x1, z1, vareps, mu, proc1, delta)
    !write(*,"('{',F16.12,',',F16.8,'},')") Deltat,Integrand_forDeltat2

end function Integrand_forDelta_th

end function KPC_SIDISconv

! Integral over theta at given Deltat
function INT_overTHETA_SIDIS(Q2, tau2, dT, x1, z1, vareps, mu, proc1, Deltat)
    real(dp), intent(in) :: Q2, tau2, dT, x1, z1, vareps, mu, Deltat
    integer, intent(in), dimension(1:3) :: proc1
    real(dp) :: INT_overTHETA_SIDIS


#if INTEGRATION_MODE==1
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!!   Adaptive/Recommended for computing observables
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INT_overTHETA_SIDIS=2*Integrate_GK(Integrand_forTHETA_SIDIS, 0.d0, pi, toleranceINT)

#elif INTEGRATION_MODE==2
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
!!!   Recommended implementations of the integral for fitting KPCs
!!! ~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
    INT_overTHETA_SIDIS=2*Integrate_K21(Integrand_forTHETA_SIDIS, 0.d0, pi)
#endif

contains

! Integrand
function Integrand_forTHETA_SIDIS(theta)
    real(dp) :: Integrand_forTHETA_SIDIS
    real(dp), intent(in) :: theta
    real(dp) :: S, Lam, xi, zeta, K, kh, cosT,RRR

    cosT = cos(theta)

    S = ((2-x1)*dT+2*Deltat*cosT*sqrt(dT*(1+dT)*(1-x1)))/x1

    Lam = (1.d0+dT)*((1+dT)+4*(2-x1)/x1**2*sqrt((1-x1)*dT*(1+dT))*Deltat*cosT&
    +4*(1-x1)/x1**2*(Deltat**2+dT*(1+(Deltat*cosT)**2)))


    xi = x1*(1-S+sqrt(Lam))/2
    zeta = 2._dp*z1/(1+S+sqrt(Lam))

    K = tau2*(Lam-(S-1._dp)**2)/4
    kh = tau2*(Lam-(S+1._dp)**2)/4

    if(K < toleranceGEN) K = toleranceGEN
    if(kh < toleranceGEN) kh = toleranceGEN

    LocalCounter = LocalCounter + 1
    RRR=TMD_pair(Q2,xi,zeta,K,kh,mu,proc1)
    Integrand_forTHETA_SIDIS = 2*(1-x1)*zeta*(1+dT)/(x1**2*sqrt(Lam))*Deltat&
    *SIDIS_KERNEL(Q2,tau2,dT*tau2, vareps, S,Lam,proc1(3))*RRR

!     write(*,*) "--------->>>", Q2, tau2, dT, x1, z1
!     write(*,*) "Deltat,cosT -->>>", Deltat,cosT
!     write(*,*) "Lam, S--->",Lam,S
!     write(*,*) "K,kh--->",K,kh
!     write(*,*) "xi,zeta--->",xi,zeta
!     write(*,*) "1--->", 2*(1-x1)*zeta*(1+dT)/(x1**2*sqrt(Lam))*Deltat
!     write(*,*) "2--->", SIDIS_KERNEL(Q2,tau2,dT*tau2,S,Lam,proc1(3))
!     write(*,*) TMD_pair(Q2,xi,zeta,K,kh,mu,proc1)
!     stop

!     write(10,*) Deltat,theta,Integrand_forTHETA_SIDIS


end function Integrand_forTHETA_SIDIS

end function INT_overTHETA_SIDIS

end module TMDF_KPC
