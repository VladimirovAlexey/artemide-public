!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 3.00
!
!    Evaluation of the small-b OPE for SiversTMDPDF
!
!    if you use this module please, quote 1901.04519
!
!    ver 3.00: release (AV, 21.07.2023)
!
!                A.Vladimirov (21.07.2023)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!---- Concept---
! This module
! * computes the tw3-convolution C x f[x,b] at the optimal point
! * creates and saves the grids for it
! * only global variables are kept here
! * the most part of the code is universal, and shared by many such modules

module SiversTMDPDF_OPE
use aTMDe_Numerics
use aTMDe_Integration
use aTMDe_IO
use aTMDe_optGrid
use QCDinput
use SiversTMDPDF_model
implicit none

!------------------------LOCALs -----------------------------------------------

private 

!Current version of module
character (len=5),parameter :: version="v3.00"
character (len=16),parameter :: moduleName="SiversTMDPDF_OPE"
!Last appropriate version of constants-file
integer,parameter::inputver=36

!--- general
logical:: started=.false.
!! Level of output
!! 0=only critical
!! 1=initialization details
!! 2=WARNINGS
integer::outputLevel=2
type(Warning_OBJ)::Warning_Handler

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

type(optGrid)::mainGrid

!!!------------------------- HARD-CODED PARAMETERS ----------------------

!!!------------------------- DYNAMICAL-GLOBAL PARAMETERS -------------------
real(dp) :: c4_global=1_dp  !!! scale variation parameter

!!!------------------------- SPECIAL VARIABLES FOR GRID (used by TMDGrid-XB)------------------
integer::numberOfHadrons=1                !!!total number of hadrons to be stored

!!--------------------------------------Public interface-----------------------------------------
public::SiversTMDPDF_OPE_IsInitialized,SiversTMDPDF_OPE_Initialize
public::SiversTMDPDF_OPE_tw3_convolution,SiversTMDPDF_OPE_tw3_resetGrid
public::SiversTMDPDF_OPE_tw3_SetPDFreplica,SiversTMDPDF_OPE_tw3_SetScaleVariation

!!!!!!----FOR TEST
!public::MakeGrid,ExtractFromGrid,CxF_compute,TestGrid

contains

!! Coefficient function
!INCLUDE 'Code/SiversTMDPDF/coeffFunc-new2.f90'


function SiversTMDPDF_OPE_IsInitialized()
    logical::SiversTMDPDF_OPE_IsInitialized
    SiversTMDPDF_OPE_IsInitialized=started
end function SiversTMDPDF_OPE_IsInitialized

!! Initialization of the package
subroutine SiversTMDPDF_OPE_Initialize(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=300)::path
    logical::initRequired
    character(len=8)::order_global
    integer::i,FILEver,messageTrigger

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

    call MoveTO(51,'*12   ')
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

    Warning_Handler=Warning_OBJ(moduleName=moduleName,messageCounter=0,messageTrigger=messageTrigger)

    c4_global=1d0

    mainGrid=optGrid(path,'*12   ','*E   ',numberOfHadrons,withGluonTW3,moduleName,outputLevel)
    
    !!! Model initialisation is called from the SiversTMDPDF-module

    if(useGridTW3) then
        call mainGrid%MakeGrid(tw3_convolution)
        if(runTestTW3) call mainGrid%Test(tw3_convolution)
    end if

    started=.true.

    if(outputLevel>0) write(*,*) color('----- arTeMiDe.SiversTMDPDF_OPE '//trim(version)//': .... initialized',c_green)
    if(outputLevel>1) write(*,*) ' '

end subroutine SiversTMDPDF_OPE_Initialize

!!!!!!!--------------------------- DEFINING ROUTINES ------------------------------------------

!!! This is defining procedure it computes the convolution of tw3 PDF with coefficient function
!!! this is convolution with twist3 PDF
function tw3_convolution(x,b,h)
    real(dp),dimension(-5:5)::tw3_convolution
    real(dp),intent(in)::x,b
    integer,intent(in)::h

    !!!! test for boundaries is done in SiversTMDPDF_lowScale5 (on the enty to this procedure)

    !!!! case NA
    if(orderMainTW3==-50) then
        if(withGluonTW3) then
            tw3_convolution=1._dp
        else
            tw3_convolution=(/1._dp,1._dp,1._dp,1._dp,1._dp,0._dp,1._dp,1._dp,1._dp,1._dp,1._dp/)
        end if
        return
    end if

    !!!!! perturbative convolution is not implemented yet

end function tw3_convolution

!!! this is convolution with twist3 PDF
function SiversTMDPDF_OPE_tw3_convolution(x,b,h,addGluon)
    real(dp),dimension(-5:5)::SiversTMDPDF_OPE_tw3_convolution
    real(dp),intent(in)::x,b
    integer,intent(in)::h
    logical,optional,intent(in)::addGluon

    if(useGridTW3) then
        SiversTMDPDF_OPE_tw3_convolution=mainGrid%Extract(x,b,h)
    else
        SiversTMDPDF_OPE_tw3_convolution=tw3_convolution(x,b,h)
    end if

end function SiversTMDPDF_OPE_tw3_convolution


!!!!!!!!!! ------------------------ SUPPORINTG ROUTINES --------------------------------------

!!! This subroutine force reconstruction of the grid (if griding is ON)
subroutine SiversTMDPDF_OPE_tw3_resetGrid()
    if(useGridTW3) then
        if(outputLevel>1) write(*,*) 'arTeMiDe ',moduleName,':  Grid Reset. with c4=',c4_global
        call mainGrid%MakeGrid(tw3_convolution)!
    end if
end subroutine SiversTMDPDF_OPE_tw3_resetGrid

!! call QCDinput to change the PDF replica number
!! unset the grid, since it should be recalculated fro different PDF replica.
subroutine SiversTMDPDF_OPE_tw3_SetPDFreplica(rep,hadron)
    integer,intent(in):: rep,hadron

    !!!! not implemented yet

end subroutine SiversTMDPDF_OPE_tw3_SetPDFreplica

!!!! this routine set the variations of scales
!!!! it is used for the estimation of errors
subroutine SiversTMDPDF_OPE_tw3_SetScaleVariation(c4_in)
    real(dp),intent(in)::c4_in
    !!!! not implemented yet
end subroutine SiversTMDPDF_OPE_tw3_SetScaleVariation

end module SiversTMDPDF_OPE
