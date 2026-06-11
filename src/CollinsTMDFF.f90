!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 3.0
!
!	Evaluation of the Collins TMD FF at low normalization point in zeta-prescription.
!
!    ver 3.03: release (SP & AV, 27.10.2025)
!
!                S.Piloneta (27.10.2025)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module Fourier_Levin_CollinsTMDFF
INCLUDE 'Code/KTspace/Fourier_Levin.f90'
end module Fourier_Levin_CollinsTMDFF

module gridInKT_CollinsTMDFF
INCLUDE 'Code/KTspace/grid_inKT.f90'
end module gridInKT_CollinsTMDFF

module CollinsTMDFF
use aTMDe_Numerics
use aTMDe_IO
use aTMDe_Ogata
use QCDinput
use TMDR
use CollinsTMDFF_OPE
use CollinsTMDFF_model
use Fourier_Levin_CollinsTMDFF
use gridInKT_CollinsTMDFF

implicit none
!------------------------LOCALs -----------------------------------------------

private 

!Current version of module
character (len=5),parameter :: version="v3.03"
character (len=12),parameter :: moduleName="CollinsTMDFF"
!Last appropriate version of constants-file
integer,parameter::inputver=37

!--------------------------------Working variables-----------------------------------------------
!--- general
logical:: started=.false.
!! Level of output
!! 0=only critical
!! 1=initialization details
!! 2=WARNINGS
integer::outputLevel=2
type(Warning_OBJ)::Warning_Handler

!!! the length and array of NP parameters
integer::lambdaNPlength
real(dp),dimension(:),allocatable::lambdaNP
real(dp)::BMAX_ABS=100._dp !!! for large values of b returns 0
real(dp)::toleranceGEN !!! tolerance general

!!!------------------------------ General parameters----------------------------------------------
logical::includeGluon=.false.   !! gluons included/non-included
integer::numOfHadrons=1         !! total number of hadrons to compute
real(dp)::TMDmass=1._dp         !! mass parameter used as mass-scale

!!!------------------------------ Parameters of transform to KT-space -------------------------------------------

integer,parameter::TMDtypeN=1 !!!!! this is the order of Bessel-transform (IT IS STRICT FOR TMD)

type(OgataIntegrator)::Hankel
!!!------------------------------ Parameters of transform in KT space and KT-grid ------------------------------
logical::makeGrid_inKT,gridIsReady_inKT
integer::numKsubgrids,kGridSize
!!!------------------------------ Parameters of transform to TMM -------------------------------------------

real(dp)::muTMM_min=0.8_dp  !!!!! minimal mu

!!-----------------------------------------------Public interface---------------------------------------------------

public::CollinsTMDFF_Initialize,CollinsTMDFF_IsInitialized
public::CollinsTMDFF_SetScaleVariation_tw3
public::CollinsTMDFF_SetFFreplica_tw3
public::CollinsTMDFF_SetLambdaNP,CollinsTMDFF_CurrentLambdaNP
public::CollinsTMDFF_inB,CollinsTMDFF_inKT,CollinsTMDFF_TMM_G,CollinsTMDFF_TMM_X

interface CollinsTMDFF_inB
    module procedure TMD_opt,TMD_ev
end interface

interface CollinsTMDFF_inKT
    module procedure TMD_opt_inKT,TMD_ev_inKT, TMD_grid_inKT
end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Interface subroutines!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function CollinsTMDFF_IsInitialized()
    logical::CollinsTMDFF_IsInitialized
    CollinsTMDFF_IsInitialized=started
end function CollinsTMDFF_IsInitialized

!! Initialization of the package
subroutine CollinsTMDFF_Initialize(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=300)::path
    logical::initRequired
    integer::FILEver,messageTrigger
    real(dp)::hOGATA_TMM,toleranceOGATA_TMM

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
        CLOSE (51, STATUS='KEEP')
        ERROR STOP
    end if

    call MoveTO(51,'*p2  ')
    read(51,*) outputLevel    
    if(outputLevel>1) write(*,*) '--------------------------------------------- '
    if(outputLevel>1) write(*,*) 'artemide.',moduleName,version,': initialization started ... '

    call MoveTO(51,'*p3  ')
    read(51,*) messageTrigger

    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p2  ')
    read(51,*) TMDmass

        !! TMDR
    call MoveTO(51,'*3   ')
    call MoveTO(51,'*p1  ')
    read(51,*) initRequired
    if(.not.initRequired) then
        write(*,*) ErrorString('TMDR module MUST be included.',moduleName)
        write(*,*) ErrorString('Check initialization-file. Evaluation stop.',moduleName)
        CLOSE (51, STATUS='KEEP')
        ERROR STOP
    end if

    call MoveTO(51,'*18   ')
    call MoveTO(51,'*p1  ')
    read(51,*) initRequired
    if(.not.initRequired) then
        if(outputLevel>1) write(*,*)'artemide.',moduleName,': initialization is not required. '
        started=.false.
        CLOSE (51, STATUS='KEEP')
        return
    end if

    !-------------general parameters
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) includeGluon
    call MoveTO(51,'*p2  ')
    read(51,*) numOfHadrons

    !-------------parameters of NP model
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) lambdaNPlength

    call MoveTO(51,'*p2  ')
    read(51,*) BMAX_ABS


    if(lambdaNPlength<=0) then
    write(*,*) ErrorString(&
    'Initialize: number of non-perturbative parameters should be >=1. Check the constants-file. Evaluation STOP',moduleName)
            CLOSE (51, STATUS='KEEP')
    ERROR STOP
    end if

    !!!!! ---- parameters of numerical evaluation
    call MoveTO(51,'*D   ')
    call MoveTO(51,'*p2  ')
    read(51,*) toleranceGEN

    !!!!! ---- parameters kT-transform and grid
    call MoveTO(51,'*F   ')
    call MoveTO(51,'*p1  ')
    read(51,*) makeGrid_inKT
    call MoveTO(51,'*p6  ')
    read(51,*) numKsubgrids
    call MoveTO(51,'*p8  ')
    read(51,*) kGridSize

    !!!!! ---- parameters of TMM-transformation
    call MoveTO(51,'*G   ')
    call MoveTO(51,'*p1  ')
    read(51,*) toleranceOGATA_TMM
    call MoveTO(51,'*p2  ')
    read(51,*) hOGATA_TMM
    call MoveTO(51,'*p3  ')
    read(51,*) muTMM_min

    CLOSE (51, STATUS='KEEP') 
    Warning_Handler=Warning_OBJ(moduleName=moduleName,messageCounter=0,messageTrigger=messageTrigger)

    if(outputLevel>2 .and. includeGluon) write(*,'(A)') ' ... gluons are included'
    if(outputLevel>2 .and. .not.includeGluon) write(*,'(A)') ' ... gluons are not included'
    if(outputLevel>2) write(*,'(A,I3)') ' Number of hadrons to be considered =',numOfHadrons
    if(outputLevel>2) write(*,'(A,I3)') ' Number of NP parameters =',lambdaNPlength
    if(outputLevel>2) write(*,'(A,F12.2)') ' Absolute maximum b      =',BMAX_ABS

    allocate(lambdaNP(1:lambdaNPlength))

    call Initialize_Fourier_Levin(path,'*18  ','*F   ',moduleName,outputLevel,TMDtypeN)
    if(makeGrid_inKT) then
        call Initialize_GridInKT(path,'*18  ','*F   ',numOfHadrons,includeGluon,moduleName,outputLevel)
    end if

    gridIsReady_inKT=.false.

    !!!!!! TODO: fix the minimal value of KT
    Hankel=OgataIntegrator(moduleName,outputLevel,TMDtypeN, toleranceOGATA_TMM,hOGATA_TMM,TMDmass,0.0001_dp)

    if(.not.TMDR_IsInitialized()) then
        if(outputLevel>2) write(*,*) '.. initializing TMDR (from ',moduleName,')'
        if(present(prefix)) then
            call TMDR_Initialize(file,prefix)
        else
            call TMDR_Initialize(file)
        end if
    end if

    if(.not.CollinsTMDFF_OPE_IsInitialized()) then
        if(outputLevel>2) write(*,*) '.. initializing CollinsTMDFF_OPE (from ',moduleName,')'
        if(present(prefix)) then
            call CollinsTMDFF_OPE_Initialize(file,prefix)
        else
            call CollinsTMDFF_OPE_Initialize(file)
        end if
    end if

    call ModelInitialization(lambdaNPlength)
    if(outputLevel>0) write(*,*) color('----- arTeMiDe.CollinsTMDFF_model : .... initialized',c_green)

    started=.true.

    if(outputLevel>0) write(*,*) color('----- arTeMiDe.CollinsTMDFF '//trim(version)//': .... initialized',c_green)
    if(outputLevel>1) write(*,*) ' '
!!!   The grid in KT can be done only after some non-perturbative values are inset
end subroutine CollinsTMDFF_Initialize

!!!!!!!!!! ------------------------ SUPPORINTG ROUTINES --------------------------------------
!!! update FF replica
subroutine CollinsTMDFF_SetFFreplica_tw3(rep,hadron)
    integer,intent(in):: rep,hadron

    call CollinsTMDFF_OPE_tw3_SetFFreplica(rep,hadron)
end subroutine CollinsTMDFF_SetFFreplica_tw3

!!!! this routine set the variations of scales
subroutine CollinsTMDFF_SetScaleVariation_tw3(c4_in)
    real(dp),intent(in)::c4_in
    call CollinsTMDFF_OPE_tw3_SetScaleVariation(c4_in)
end subroutine CollinsTMDFF_SetScaleVariation_tw3

!!!Sets the non-pertrubative parameters lambda
!!! carries additionl option to build the grid
!!! if need to build grid, specify the gluon required directive.
subroutine CollinsTMDFF_SetLambdaNP(lambdaIN)
    real(dp),intent(in)::lambdaIN(:)
    integer::ll
    call Warning_Handler%Reset()

    ll=size(lambdaIN)
    if(ll/=lambdaNPlength) then
        if(outputLevel>0) write(*,"(A,I3,A,I3,')')")&
                WarningString('SetLambdaNP:length of lambdaNP(',moduleName),&
            ll,color(') does not match the required (',c_red),lambdaNPlength
        if(outputLevel>0) write(*,*)color('                NOTHING IS DONE!',c_red)
        return
    end if

    lambdaNP=lambdaIN
    call ModelUpdate(lambdaNP)
    gridIsReady_inKT=.false.

    if(outputLevel>2) write(*,*) 'arTeMiDe.',moduleName,': NPparameters reset = (',lambdaNP,')'

    if(makeGrid_inKT) then
        call updateGrid_inKT()
        gridIsReady_inKT=.true.
    end if
end subroutine CollinsTMDFF_SetLambdaNP

!!! returns current value of NP parameters
function CollinsTMDFF_CurrentLambdaNP()
    real(dp),dimension(1:lambdaNPlength)::CollinsTMDFF_CurrentLambdaNP
    CollinsTMDFF_CurrentLambdaNP=lambdaNP
end function CollinsTMDFF_CurrentLambdaNP

!!!!! Function that constracts the grid in KT-space
!!!!! I take the TMD and send it to Fourier, and then resend it to the grid...
subroutine updateGrid_inKT()
real(dp)::Q,x
integer::h

call PrepareGrid_inKT(toGrid)

contains

function toFourier(b)
    real(dp),dimension(-5:5) :: toFourier
    real(dp), intent(in) ::b
    real(dp):: Rkernel,RkernelG

    if(includeGluon) then
    Rkernel=TMDR_Rzeta(b,Q,Q**2,1)
    RkernelG=TMDR_Rzeta(b,Q,Q**2,0)

    toFourier=CollinsTMDFF_OPE_tw3_convolution(x,b,abs(h))*FNP(x,b,abs(h),lambdaNP)*&
        (/Rkernel,Rkernel,Rkernel,Rkernel,Rkernel,RkernelG,Rkernel,Rkernel,Rkernel,Rkernel,Rkernel/)

    else
        Rkernel=TMDR_Rzeta(b,Q,Q**2,1)
        toFourier=Rkernel*CollinsTMDFF_OPE_tw3_convolution(x,b,abs(h))*FNP(x,b,abs(h),lambdaNP)
    end if
end function toFourier

function toGrid(x_in,Q_in,h_in,arraySize1,arraySize2)
    integer,intent(in)::arraySize1,arraySize2
    real(dp),dimension(1:arraySize1,0:arraySize2,-5:5) :: toGrid
    real(dp), intent(in) ::x_in,Q_in
    integer,intent(in)::h_in

    x=x_in
    Q=Q_in
    h=h_in

    toGrid=Fourier_Levin_array(toFourier)
end function toGrid

end subroutine updateGrid_inKT

! subroutine testGrid_inKT()
!
! call TestGrid_inKT_internal(ToCompare)
!
! contains
!
! function ToCompare(x,kT,Q,h)
! real(dp),dimension(-5:5)::ToCompare
! real(dp),intent(in)::x,kT,Q
! integer,intent(in)::h
! ToCompare=TMD_ev_inKT(x,kT,Q,Q**2,h)
! end function ToCompare
!
! end subroutine testGrid_inKT

!!!!!!!--------------------------- DEFINING ROUTINES ------------------------------------------

!!!!! the names are neutral because these procedures are feed to Fourier transform. And others universal sub programs.

!!!!!!! the function that actually returns the CollinsTMDFF!
function TMD_opt(x,bT,hadron)
  real(dp),dimension(-5:5)::TMD_opt
  real(dp),intent(in) :: x, bT
  integer,intent(in)::hadron

  !!! test boundaries
    if(x>1d0) then
        call Warning_Handler%WarningRaise('Called x>1 (return 0). x='//numToStr(x))
        TMD_opt=0._dp
        return
    else if(x==1.d0) then !!! funny but sometimes FORTRAN can compare real numbers exactly
        TMD_opt=0._dp
        return
    else if(bT>BMAX_ABS) then
        TMD_opt=0._dp
        return
    else if(x<1d-12) then
        ERROR STOP ErrorString('Called x<0. x='//numToStr(x)//' . Evaluation STOP',moduleName)
    else if(bT<0d0) then
        ERROR STOP ErrorString('Called b<0. b='//numToStr(bT)//' . Evaluation STOP',moduleName)
    end if

    TMD_opt=CollinsTMDFF_OPE_tw3_convolution(x,bT,abs(hadron))*FNP(x,bT,abs(hadron),lambdaNP)

    if(hadron<0) TMD_opt=TMD_opt(5:-5:-1)

end function TMD_opt

!!!!!!!! the function that actually returns the CollinsTMDFF evolved to (mu,zeta) value
function TMD_ev(x,bt,muf,zetaf,hadron)
    real(dp)::TMD_ev(-5:5)
    real(dp),intent(in):: x,bt,muf,zetaf
    integer,intent(in)::hadron
    real(dp):: Rkernel,RkernelG

    if(includeGluon) then
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
        RkernelG=TMDR_Rzeta(bt,muf,zetaf,0)

        TMD_ev=TMD_opt(x,bT,hadron)*&
            (/Rkernel,Rkernel,Rkernel,Rkernel,Rkernel,RkernelG,Rkernel,Rkernel,Rkernel,Rkernel,Rkernel/)

    else
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
        TMD_ev=Rkernel*TMD_opt(x,bT,hadron)
    end if


    !!! forcefully set =0 below threshold
    if(muf<mBOTTOM) then
    TMD_ev(5)=0_dp
    TMD_ev(-5)=0_dp
    end if
    if(muf<mCHARM) then
    TMD_ev(4)=0_dp
    TMD_ev(-4)=0_dp
    end if

end function TMD_ev

!!!!! the names are neutral because these procedures are feed to Fourier transform. And others universal sub programs.

!!!!!!! the function that actually returns the uTMDPDF optimal value
function TMD_opt_inKT(x,kT,hadron)
  real(dp),dimension(-5:5)::TMD_opt_inKT
  real(dp),intent(in) :: x, kT
  integer,intent(in)::hadron

  !!! test boundaries
    if(x>1d0) then
        call Warning_Handler%WarningRaise('Called x>1 (return 0). x='//numToStr(x))
        TMD_opt_inKT=0._dp
        return
     else if(x==1.d0) then !!! funny but sometimes FORTRAN can compare real numbers exactly
        TMD_opt_inKT=0._dp
        return
    else if(x<1d-12) then
        ERROR STOP ErrorString('Called x<0. x='//numToStr(x)//' . Evaluation STOP',moduleName)
    else if(kT<0d0) then
        ERROR STOP ErrorString('Called kT<0. kT='//numToStr(kT)//' . Evaluation STOP',moduleName)
    end if

    TMD_opt_inKT=Fourier_Levin(toFourier,kT)

    if(hadron<0) TMD_opt_inKT=TMD_opt_inKT(5:-5:-1)

    contains

    function toFourier(b)
        real(dp),dimension(-5:5) :: toFourier
        real(dp), intent(in) ::b
        toFourier=CollinsTMDFF_OPE_tw3_convolution(x,b,abs(hadron))*FNP(x,b,abs(hadron),lambdaNP)
    end function toFourier

end function TMD_opt_inKT
!
!!!!!!!! the function that actually returns the uTMDPDF evolved to (mu,zeta) value
function TMD_ev_inKT(x,kT,muf,zetaf,hadron)
    real(dp)::TMD_ev_inKT(-5:5)
    real(dp),intent(in):: x,kT,muf,zetaf
    integer,intent(in)::hadron
    real(dp):: Rkernel,RkernelG

    TMD_ev_inKT=Fourier_Levin(toFourier,kT)

    if(hadron<0) TMD_ev_inKT=TMD_ev_inKT(5:-5:-1)

    !!! forcefully set =0 below threshold
    if(muf<mBOTTOM) then
    TMD_ev_inKT(5)=0_dp
    TMD_ev_inKT(-5)=0_dp
    end if
    if(muf<mCHARM) then
    TMD_ev_inKT(4)=0_dp
    TMD_ev_inKT(-4)=0_dp
    end if

contains

function toFourier(b)
    real(dp),dimension(-5:5) :: toFourier
    real(dp), intent(in) ::b

    if(includeGluon) then
    Rkernel=TMDR_Rzeta(b,muf,zetaf,1)
    RkernelG=TMDR_Rzeta(b,muf,zetaf,0)

    toFourier=CollinsTMDFF_OPE_tw3_convolution(x,b,abs(hadron))*FNP(x,b,abs(hadron),lambdaNP)*&
        (/Rkernel,Rkernel,Rkernel,Rkernel,Rkernel,RkernelG,Rkernel,Rkernel,Rkernel,Rkernel,Rkernel/)

    else
        Rkernel=TMDR_Rzeta(b,muf,zetaf,1)
        toFourier=Rkernel*CollinsTMDFF_OPE_tw3_convolution(x,b,abs(hadron))*FNP(x,b,abs(hadron),lambdaNP)
    end if
end function toFourier

end function TMD_ev_inKT


!!!!!!!! the function that actually returns the uTMDPDF evolved to (mu,mu^2) value
!!!!!!! This is exactly what is stored in the grid. So if grid is buid attempt to extract from it.
function TMD_grid_inKT(x,kT,mu,hadron)
    real(dp)::TMD_grid_inKT(-5:5)
    real(dp),intent(in):: x,kT,mu
    integer,intent(in)::hadron

    if(gridIsReady_inKT) then
!         if(ISNAN(kT)) then
!             write(*,*) "HERE",kT,x,mu
!             stop
!         end if

        TMD_grid_inKT=ExtractFromGrid_inKT(x,kT,mu,abs(hadron))

        if(hadron<0) TMD_grid_inKT=TMD_grid_inKT(5:-5:-1)

    else if(makeGrid_inKT) then
        ERROR STOP ErrorString("Attempt to extract TMD from grid, while grid is not ready. CHECK!",moduleName)
    else
        TMD_grid_inKT=TMD_ev_inKT(x,kT,mu,mu**2,hadron)
    end if

end function TMD_grid_inKT

!!!!!!!! TMM G_{n,n} at (x,mu)
function CollinsTMDFF_TMM_G(x,mu,hadron)
    real(dp)::CollinsTMDFF_TMM_G(-5:5)
    real(dp),intent(in):: x,mu
    integer,intent(in)::hadron

    if(mu<muTMM_min) then
        write(*,*) ErrorString("ERROR in KT-moment computation. mu<mu_min:",moduleName),muTMM_min
        error stop
    end if

    CollinsTMDFF_TMM_G=Hankel%Moment_G(F,mu)


    !!! forcefully set =0 below threshold
    if(mu<mBOTTOM) then
    CollinsTMDFF_TMM_G(5)=0_dp
    CollinsTMDFF_TMM_G(-5)=0_dp
    end if
    if(mu<mCHARM) then
    CollinsTMDFF_TMM_G(4)=0_dp
    CollinsTMDFF_TMM_G(-4)=0_dp
    end if
contains
    function F(b)
    real(dp),dimension(-5:5)::F
    real(dp),intent(in)::b
    F=TMD_opt(x,b,hadron)
    end function F

end function CollinsTMDFF_TMM_G

!!!!!!!! TMM G_{n+1,n} at (x,mu)
function CollinsTMDFF_TMM_X(x,mu,hadron)
    real(dp)::CollinsTMDFF_TMM_X(-5:5)
    real(dp),intent(in):: x,mu
    integer,intent(in)::hadron

    if(mu<muTMM_min) then
        write(*,*) ErrorString("ERROR in KT-moment computation. mu<mu_min:",moduleName),muTMM_min
        error stop
    end if

    CollinsTMDFF_TMM_X=Hankel%Moment_X(F,mu)

    !!! forcefully set =0 below threshold
    if(mu<mBOTTOM) then
    CollinsTMDFF_TMM_X(5)=0_dp
    CollinsTMDFF_TMM_X(-5)=0_dp
    end if
    if(mu<mCHARM) then
    CollinsTMDFF_TMM_X(4)=0_dp
    CollinsTMDFF_TMM_X(-4)=0_dp
    end if
contains
    function F(b)
    real(dp),dimension(-5:5)::F
    real(dp),intent(in)::b
    F=TMD_opt(x,b,hadron)
    end function F

end function CollinsTMDFF_TMM_X

end module CollinsTMDFF
