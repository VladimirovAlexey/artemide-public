!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 3.00
!
!    Evaluation of the small-b OPE for BoerMuldersTMDPDF
!
!    if you use this module please, quote     2209.00962
!
!    ver 3.00: release (AV, 29.01.2024)
!
!                A.Vladimirov (29.07.2024)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!---- Concept---
! This module
! * computes the tw3-convolution C x f[x,b] at the optimal point
! * creates and saves the grids for it
! * only global variables are kept here
! * the most part of the code is universal, and shared by many such modules

module BoerMuldersTMDPDF_OPE
use aTMDe_Numerics
use IntegrationRoutines
use IO_functions
use QCDinput
use BoerMuldersTMDPDF_model
implicit none

!------------------------LOCALs -----------------------------------------------

private 

!Current version of module
character (len=5),parameter :: version="v3.00"
character (len=21),parameter :: moduleName="BoerMuldersTMDPDF_OPE"
!Last appropriate version of constants-file
integer,parameter::inputver=30

!--- general
logical:: started=.false.
!! Level of output
!! 0=only critical
!! 1=initialization details
!! 2=WARNINGS
integer::outputLevel=2
!! variable that count number of WARNING mesagges. In order not to spam too much
integer::messageTrigger=6
integer::messageCounter=0 !!! actual counter

!!!------------------------- PARAMETERS DEFINED IN THE INI-file--------------

!!! Perturbative order
integer :: orderMainTW3=-50 !! LO=0, NLO=1,...

!!!! X-Grid parameters

!!!! B-Grid parameters

!!!! Numerical parameters
real(dp) :: toleranceINT=1d-6  !!! tolerance for numerical integration
real(dp) :: toleranceGEN=1d-6  !!! tolerance for other purposes
integer :: maxIteration=4000   !!! maximum iteration in the integrals (not used at the moment)

!!!! grid preparation for tw3 part
logical :: useGridTW3=.false.  !!!idicator that grid must be prepared
logical :: withGluonTW3=.false.   !!!indicator the gluon is needed in the grid
logical :: runTestTW3=.false.   !!!trigger to run the test

!!!------------------------- HARD-CODED PARAMETERS ----------------------

!!!------------------------- DYNAMICAL-GLOBAL PARAMETERS -------------------
real(dp) :: c4_global=1_dp  !!! scale variation parameter
logical :: gridReady!!!!indicator that grid is ready to use. If it is .true., the TMD calculated from the grid

!!!------------------------- SPECIAL VARIABLES FOR GRID (used by TMDGrid-XB)------------------
integer::numberOfHadrons=1                !!!total number of hadrons to be stored

!!--------------------------------------Public interface-----------------------------------------
public::BoerMuldersTMDPDF_OPE_IsInitialized,BoerMuldersTMDPDF_OPE_Initialize
public::BoerMuldersTMDPDF_OPE_tw3_convolution,BoerMuldersTMDPDF_OPE_tw3_resetGrid,BoerMuldersTMDPDF_OPE_tw3_testGrid
public::BoerMuldersTMDPDF_OPE_tw3_SetPDFreplica,BoerMuldersTMDPDF_OPE_tw3_SetScaleVariation

!!!!!!----FOR TEST
!public::MakeGrid,ExtractFromGrid,CxF_compute,TestGrid

contains

!! Coefficient function
!INCLUDE 'Code/BoerMuldersTMDPDF/coeffFunc-new2.f90'


function BoerMuldersTMDPDF_OPE_IsInitialized()
    logical::BoerMuldersTMDPDF_OPE_IsInitialized
    BoerMuldersTMDPDF_OPE_IsInitialized=started
end function BoerMuldersTMDPDF_OPE_IsInitialized

!! Initialization of the package
subroutine BoerMuldersTMDPDF_OPE_Initialize(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=300)::path
    logical::initRequired
    character(len=8)::order_global
    integer::i,FILEver

    if(started) return

    if(.not.QCDinput_IsInitialized()) then
        if(outputLevel>1) write(*,*) '.. initializing QCDinput (from ',moduleName,')'
        if(present(prefix)) then
            call QCDinput_Initialize(file,prefix)
        else
            call QCDinput_Initialize(file)
        end if
    end if

    if(present(prefix)) then
        path=trim(adjustl(prefix))//trim(adjustr(file))
    else
        path=trim(adjustr(file))
    end if

    !----------------- reading ini-file --------------------------------------
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
    if(outputLevel>1) write(*,*) '--------------------------------------------- '
    if(outputLevel>1) write(*,*) 'artemide.',moduleName," ",version,': initialization started ... '
    call MoveTO(51,'*p3  ')
    read(51,*) messageTrigger

    !!!!!!!!!!! OPEN MP Initialization
    !$ call MoveTO(51,'*C   ')
    !$ call MoveTO(51,'*p1  ')
    !$read(51,*) i
    !$ call OMP_set_num_threads(i)

    call MoveTO(51,'*14   ')
    call MoveTO(51,'*p1  ')
    read(51,*) initRequired
    if(.not.initRequired) then
        if(outputLevel>1) write(*,*)'artemide.',moduleName,': initialization is not required. '
        started=.false.
        CLOSE (51, STATUS='KEEP')
        return
    end if

    !----- GENERALS
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) withGluonTW3
    call MoveTO(51,'*p2  ')
    read(51,*) numberOfHadrons

    !----- ORDER
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) order_global

    SELECT CASE(trim(order_global))
        CASE ("NA")
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NA',color(" (TMD=fNP)",c_yellow)
            orderMainTW3=-50
        CASE DEFAULT
            if(outputLevel>0)write(*,*) &
                WarningString('Initialize: tw3-convolution is not implemented so far. Switch to NA.',moduleName)
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NA'
            orderMainTW3=-50
        END SELECT

    if(outputLevel>2 .and. orderMainTW3>-1) write(*,'(A,I1)') ' |  Coef.func.    =as^',orderMainTW3

    call MoveTO(51,'*p2  ')
    read(51,*) useGridTW3
    call MoveTO(51,'*p3  ')
    read(51,*) runTestTW3

    !!!!! ---- parameters of numerical evaluation
    call MoveTO(51,'*D   ')
    call MoveTO(51,'*p1  ')
    read(51,*) toleranceINT
    call MoveTO(51,'*p2  ')
    read(51,*) toleranceGEN
    call MoveTO(51,'*p3  ')
    read(51,*) maxIteration

    if(outputLevel>2) then
        write(*,'(A,ES10.3)') ' |  tolerance     =',toleranceINT
        write(*,'(A,ES10.3)') ' |  max iteration =',REAL(maxIteration)
     end if

    CLOSE (51, STATUS='KEEP')
    c4_global=1d0

    !!! call initializations for Grids
!     call XGrid_Initialize()
!     call BGrid_Initialize()
!     call TMDGrid_XB_Initialize()
    
    !!! Model initialisation is called from the BoerMuldersTMDPDF-module
    gridReady=.false.

!     if(useGrid) then
!         call MakeGrid()
!         gridReady=.true.
!
!         if(runTest) call TestGrid()
!     end if

    started=.true.
    messageCounter=0

    if(outputLevel>0) write(*,*) color('----- arTeMiDe.BoerMuldersTMDPDF_OPE '//trim(version)//': .... initialized',c_green)
    if(outputLevel>1) write(*,*) ' '

end subroutine BoerMuldersTMDPDF_OPE_Initialize

!!!!!!!--------------------------- DEFINING ROUTINES ------------------------------------------

! !!!!array of x times PDF(x,Q) for hadron 'hadron'
! !!!! array is (-5:5) (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
! function xf(x,Q,hadron)
!     real(dp) :: x,Q
!     integer:: hadron
!     real(dp), dimension(-5:5):: xf
!
!     xf=x_hPDF(x,Q,hadron)
!
! end function xf

!!! this is convolution with twist3 PDF
function BoerMuldersTMDPDF_OPE_tw3_convolution(x,b,h,addGluon)
    real(dp),dimension(-5:5)::BoerMuldersTMDPDF_OPE_tw3_convolution
    real(dp),intent(in)::x,b
    integer,intent(in)::h
    logical,optional,intent(in)::addGluon

    logical::gluon

    !!! check gluonity
    if(present(addGluon)) then
        gluon=addGluon
    else
        gluon=withGluonTW3
    end if

    !!!! test for boundaries is done in BoerMuldersTMDPDF_lowScale5 (on the enty to this procedure)

    !!!! case NA
    if(orderMainTW3==-50) then
        if(gluon) then
            BoerMuldersTMDPDF_OPE_tw3_convolution=1._dp
        else
            BoerMuldersTMDPDF_OPE_tw3_convolution=(/1._dp,1._dp,1._dp,1._dp,1._dp,0._dp,1._dp,1._dp,1._dp,1._dp,1._dp/)
        end if
        return
    end if

    !!!!! perturbative convolution is not implemented yet

end function BoerMuldersTMDPDF_OPE_tw3_convolution


!!!!!!!!!! ------------------------ SUPPORINTG ROUTINES --------------------------------------

!!! This subroutine force reconstruction of the grid (if griding is ON)
subroutine BoerMuldersTMDPDF_OPE_tw3_resetGrid()
    !!!!! not implemented yet
end subroutine BoerMuldersTMDPDF_OPE_tw3_resetGrid

!!! This subroutine force reconstruction of the grid (if griding is ON)
subroutine BoerMuldersTMDPDF_OPE_tw3_testGrid()
    !!!!! not implemented yet
end subroutine BoerMuldersTMDPDF_OPE_tw3_testGrid

!! call QCDinput to change the PDF replica number
!! unset the grid, since it should be recalculated fro different PDF replica.
subroutine BoerMuldersTMDPDF_OPE_tw3_SetPDFreplica(rep,hadron)
    integer,intent(in):: rep,hadron

    !!!! not implemented yet

end subroutine BoerMuldersTMDPDF_OPE_tw3_SetPDFreplica

!!!! this routine set the variations of scales
!!!! it is used for the estimation of errors
subroutine BoerMuldersTMDPDF_OPE_tw3_SetScaleVariation(c4_in)
    real(dp),intent(in)::c4_in
    !!!! not implemented yet
end subroutine BoerMuldersTMDPDF_OPE_tw3_SetScaleVariation

end module BoerMuldersTMDPDF_OPE
