!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 3.0
!
!	Evaluation of the unpolarized TMD PDF at low normalization point in zeta-prescription.
!	
!	if you use this module please, quote 1706.01473
!
!	18.08.2023  Implementation in ver.3.0
!
!				A.Vladimirov (18.08.2023)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module eeTMDFF
use aTMDe_Numerics
use IO_functions
use QCDinput
use TMDR
use eeTMDFF_model

implicit none
!------------------------LOCALs -----------------------------------------------

private 

!Current version of module
character (len=5),parameter :: version="v3.01"
character (len=7),parameter :: moduleName="eeTMDFF"
!Last appropriate version of constants-file
integer,parameter::inputver=32

!--------------------------------Working variables-----------------------------------------------
!--- general
logical:: started=.false.
!! Level of output
!! 0=only critical
!! 1=initialization details
!! 2=WARNINGS
integer::outputLevel=2
!! variable that count number of WRNING mesagges. In order not to spam too much
integer::messageTrigger=6

!!! the length and array of NP parameters
integer::lambdaNPlength
real(dp),dimension(:),allocatable::lambdaNP
real(dp)::BMAX_ABS=100._dp !!! for large values of b returns 0
!!!! Numerical parameters
real(dp) :: toleranceINT=1d-6  !!! tolerance for numerical integration
real(dp) :: toleranceGEN=1d-6  !!! tolerance for other purposes
integer :: maxIteration=4000   !!! maximum iteration in the integrals (not used at the moment)

integer :: messageCounter

!!! ---------------------------- Perturbative variables
!!! Perturbative order
integer :: orderMain=2 !! LO=0, NLO=1,...
real(dp) :: c4_global=1_dp  !!! scale variation parameter

!!!------------------------------ General parameters----------------------------------------------
logical::includeGluon=.false.   !! gluons included/non-included
integer::numOfHadrons=1         !! total number of hadrons to compute
real(dp)::TMDmass=1._dp         !! mass parameter used as mass-scale

!!!------------------------------ Parameters of transform to KT-space -------------------------------------------

integer,parameter::TMDtypeN=0 !!!!! this is the order of Bessel-transform (IT IS STRICT FOR TMD)

!!!------------------------------ Parameters of transform in KT space and KT-grid ------------------------------
logical::makeGrid_inKT,gridIsReady_inKT

!!-----------------------------------------------Public interface---------------------------------------------------
public::eeTMDFF_Initialize,eeTMDFF_IsInitialized,eeTMDFF_SetScaleVariation
public::eeTMDFF_SetLambdaNP,eeTMDFF_CurrentLambdaNP
public::eeTMDFF_inB,eeTMDFF_inKT

interface eeTMDFF_inB
    module procedure TMD_opt,TMD_ev
end interface

interface eeTMDFF_inKT
    module procedure TMD_opt_inKT,TMD_ev_inKT
end interface

contains

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Interface subroutines!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function eeTMDFF_IsInitialized()
    logical::eeTMDFF_IsInitialized
    eeTMDFF_IsInitialized=started
end function eeTMDFF_IsInitialized

!! Initialization of the package
subroutine eeTMDFF_Initialize(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=300)::path
    character(len=8)::order_global
    logical::initRequired
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
        CLOSE (51, STATUS='KEEP')
        stop
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
        stop
    end if

    call MoveTO(51,'*17  ')
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

    !-------------general parameters
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) order_global

    SELECT CASE(trim(order_global))
        CASE ("NA")
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NA',color(" (TMD=fNP)",c_yellow)
            orderMain=-50
        CASE ("LO")
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: LO'
            orderMain=0
        CASE ("NLO")
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NLO'
            orderMain=1
        CASE ("NNLO")
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NNLO'
            orderMain=2
        CASE ("N2LO")
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NNLO'
            orderMain=2
        CASE ("NNNLO")
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: N3LO'
            orderMain=3
        CASE ("N3LO")
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: N3LO'
            orderMain=3
        CASE DEFAULT
            if(outputLevel>0)write(*,*) &
                WarningString('Initialize: unknown order for coefficient function. Switch to NLO.',moduleName)
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NLO'
            orderMain=1
    END SELECT

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
        stop
    end if

    !!!!! ---- parameters of numerical evaluation
    call MoveTO(51,'*D   ')
    call MoveTO(51,'*p1  ')
    read(51,*) toleranceINT
    call MoveTO(51,'*p2  ')
    read(51,*) toleranceGEN
    call MoveTO(51,'*p3  ')
    read(51,*) maxIteration

    !!!!! ---- parameters kT-transform and grid
    call MoveTO(51,'*F   ')
    call MoveTO(51,'*p1  ')
    read(51,*) makeGrid_inKT
    if(makeGrid_inKT) write(*,*) WarningString('The KT-transformation is not yet implemented.',moduleName)

    CLOSE (51, STATUS='KEEP') 

    c4_global=1d0

    if(outputLevel>2 .and. includeGluon) write(*,'(A)') ' ... gluons are included'
    if(outputLevel>2 .and. .not.includeGluon) write(*,'(A)') ' ... gluons are not included'
    if(outputLevel>2) write(*,'(A,I3)') ' Number of hadrons to be considered =',numOfHadrons
    if(outputLevel>2) write(*,'(A,I3)') ' Number of NP parameters =',lambdaNPlength
    if(outputLevel>2) write(*,'(A,F12.2)') ' Absolute maximum b      =',BMAX_ABS

    allocate(lambdaNP(1:lambdaNPlength))

    call ModelInitialization(lambdaNPlength)
    if(outputLevel>0) write(*,*) color('----- arTeMiDe.eeTMDFF_model : .... initialized',c_green)


    if(.not.TMDR_IsInitialized()) then
        if(outputLevel>2) write(*,*) '.. initializing TMDR (from ',moduleName,')'
        if(present(prefix)) then
            call TMDR_Initialize(file,prefix)
        else
            call TMDR_Initialize(file)
        end if
    end if

    started=.true.
    messageCounter=0

    if(outputLevel>0) write(*,*) color('----- arTeMiDe.eeTMDFF '//trim(version)//': .... initialized',c_green)
    if(outputLevel>1) write(*,*) ' '

end subroutine eeTMDFF_Initialize

!!!!!!!!!! ------------------------ SUPPORINTG ROUTINES --------------------------------------

!!!! this routine set the variations of scales
subroutine eeTMDFF_SetScaleVariation(c4_in)
    real(dp),intent(in)::c4_in
    if(c4_in<0.1d0 .or. c4_in>10.d0) then
        if(outputLevel>0) write(*,*) WarningString('variation in c4 is enourmous. c4 is set to 2',moduleName)
        c4_global=2d0
    else if(abs(c4_in-c4_global)<toleranceGEN) then
        if(outputLevel>1) write(*,*) color('uTMDFF: c4-variation is ignored. c4='//real8ToStr(c4_global),c_yellow)
    else
        c4_global=c4_in
        if(outputLevel>1) write(*,*) color('uTMDFF: set scale variations c4 as:'//real8ToStr(c4_global),c_yellow)
    end if
end subroutine eeTMDFF_SetScaleVariation

!!!Sets the non-pertrubative parameters lambda
!!! carries additionl option to build the grid
!!! if need to build grid, specify the gluon required directive.
subroutine eeTMDFF_SetLambdaNP(lambdaIN)
    real(dp),intent(in)::lambdaIN(:)
    integer::ll
    messageCounter=0

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
    if(makeGrid_inKT) then
        write(*,*) ErrorString('KT-space of Jet is not implemented . Evaluation STOP',moduleName)
        ERROR STOP
    end if

    if(outputLevel>2) write(*,*) 'arTeMiDe.',moduleName,': NPparameters reset = (',lambdaNP,')'

end subroutine eeTMDFF_SetLambdaNP

!!! returns current value of NP parameters
function eeTMDFF_CurrentLambdaNP()
    real(dp),dimension(1:lambdaNPlength)::eeTMDFF_CurrentLambdaNP
    eeTMDFF_CurrentLambdaNP=lambdaNP
end function eeTMDFF_CurrentLambdaNP

!!!!!!!--------------------------- DEFINING ROUTINES ------------------------------------------


!!!!! the names are neutral because these procedures are feed to Fourier transform. And others universal sub programs.

!!!!!!! the function that actually returns the eeTMDFF optimal value
function TMD_opt(bT,hadron)
  real(dp),dimension(-5:5)::TMD_opt
  real(dp),intent(in) :: bT
  integer,intent(in)::hadron
  real(dp),dimension(-5:5)::Jet_beam_function

  !!! test boundaries
    if(bT>BMAX_ABS) then
        TMD_opt=0._dp
        return
    else if(bT<0d0) then
        write(*,*) ErrorString('Called b<0. b='//numToStr(bT)//' . Evaluation STOP',moduleName)
        ERROR STOP
    end if

    !!!!!-----------------coefficient function to write -------------------------
    !!!! b-> bSTAR(bT)
    !!!! mu->muOPE(bt,c4_global)

    Jet_beam_function=1.d0*(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)!!!
    if(orderMain>0) then !!NLO
        Jet_beam_function=Jet_beam_function
        if(orderMain>1) then !! NNLO
            Jet_beam_function=Jet_beam_function
        end if
    end if
    !!!!!-----------------coefficient function to write -------------------------

    TMD_opt=Jet_beam_function*FNP(bT,abs(hadron),lambdaNP)

    if(hadron<0) TMD_opt=TMD_opt(5:-5:-1)

end function TMD_opt
!
!!!!!!!! the function that actually returns the eeTMDFF evolved to (mu,zeta) value
function TMD_ev(bt,muf,zetaf,hadron)
    real(dp)::TMD_ev(-5:5)
    real(dp),intent(in):: bt,muf,zetaf
    integer,intent(in)::hadron
    real(dp):: Rkernel,RkernelG

    if(includeGluon) then
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
        RkernelG=TMDR_Rzeta(bt,muf,zetaf,0)

        TMD_ev=TMD_opt(bT,hadron)*&
            (/Rkernel,Rkernel,Rkernel,Rkernel,Rkernel,RkernelG,Rkernel,Rkernel,Rkernel,Rkernel,Rkernel/)

    else
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
        TMD_ev=Rkernel*TMD_opt(bT,hadron)
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

!!!!!!! the function that actually returns the eeTMDFF optimal value
function TMD_opt_inKT(kT,hadron)
  real(dp),dimension(-5:5)::TMD_opt_inKT
  real(dp),intent(in) :: kT
  integer,intent(in)::hadron

    write(*,*) ErrorString('KT-space of Jet is not implemented . Evaluation STOP',moduleName)
    ERROR STOP

end function TMD_opt_inKT
!
!!!!!!!! the function that actually returns the eeTMDFF evolved to (mu,zeta) value
function TMD_ev_inKT(kT,muf,zetaf,hadron)
    real(dp)::TMD_ev_inKT(-5:5)
    real(dp),intent(in):: kT,muf,zetaf
    integer,intent(in)::hadron
    real(dp):: Rkernel,RkernelG

    write(*,*) ErrorString('KT-space of Jet is not implemented . Evaluation STOP',moduleName)
    ERROR STOP
end function TMD_ev_inKT

end module eeTMDFF
