!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 3.00
!
!    Evaluation of the small-b OPE for wglTMDPDF
!
!    if you use this module please, quote 2209.00962
!
!    ver 3.00: release (AV, 07.01.2024)
!
!                A.Vladimirov (07.01.2024)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!---- Concept---
! This module
! * computes the tw2-convolution C x f[x,b] at the optimal point
! * creates and saves the grids for it
! * only global variables are kept here
! * the most part of the code is universal, and shared by many such modules

module Grid_wglTMDPDF
INCLUDE 'Code/Twist2/Twist2_ChGrid.f90'
end module Grid_wglTMDPDF

module wglTMDPDF_OPE
use aTMDe_Numerics
use IntegrationRoutines
use IO_functions
use QCDinput
use TMD_AD, only : Dpert_atL
use wglTMDPDF_model
use Grid_wglTMDPDF
implicit none

!------------------------LOCALs -----------------------------------------------

private 

!Current version of module
character (len=5),parameter :: version="v3.01"
character (len=13),parameter :: moduleName="wglTMDPDF_OPE"
!Last appropriate version of constants-file
integer,parameter::inputver=34

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
integer :: orderMain=2 !! LO=0, NLO=1,...
integer :: orderMainTW3=-50 !! LO=0, NLO=1,...
!!! Order of large-X resummation
logical :: resumLargeX=.false. !!!! could not be resummation of largeX
integer :: orderLX=0 !! LO=0 [no-resummation], NLO=1,...

!!! Phase space limitations parameters
real(dp) :: xMin=0.001_dp !!! min x
real(dp) :: BMAX=25._dp !!! maximum B
real(dp) :: BMIN=1d-6 !!! minimum B

!!!! Numerical parameters
real(dp) :: toleranceINT=1d-6  !!! tolerance for numerical integration
real(dp) :: toleranceGEN=1d-6  !!! tolerance for other purposes
integer :: maxIteration=4000   !!! maximum iteration in the integrals (not used at the moment)

!!!total number of hadrons to be used
integer::numberOfHadrons=1

logical(dp) :: IsMuYdependent = .true.  !!! if mu is y independent, computation is much(!) faster

!!!! grid preparation
logical :: useGrid=.true.  !!!idicator that grid must be prepared
logical :: withGluon=.false.   !!!indicator the gluon is needed in the grid
logical :: runTest=.false.   !!!trigger to run the test

!!!! grid preparation for tw3 part
logical :: useGridTW3=.false.  !!!idicator that grid must be prepared
logical :: withGluonTW3=.false.   !!!indicator the gluon is needed in the grid
logical :: runTestTW3=.false.   !!!trigger to run the test

!!!------------------------- HARD-CODED PARAMETERS ----------------------
!!! Coefficient lists
integer,parameter::parametrizationLength=4 !!![exact]

!!!------------------------- DYNAMICAL-GLOBAL PARAMETERS -------------------
real(dp) :: c4_global=1_dp  !!! scale variation parameter
logical :: gridReady!!!!indicator that grid is ready to use. If it is .true., the TMD calculated from the grid



!!--------------------------------------Public interface-----------------------------------------
public::wglTMDPDF_OPE_IsInitialized,wglTMDPDF_OPE_Initialize,wglTMDPDF_OPE_convolution
public::wglTMDPDF_OPE_resetGrid,wglTMDPDF_OPE_SetPDFreplica,wglTMDPDF_OPE_SetScaleVariation
public::wglTMDPDF_OPE_tw3_convolution,wglTMDPDF_OPE_tw3_testGrid,wglTMDPDF_OPE_tw3_resetGrid
public::wglTMDPDF_OPE_tw3_SetPDFreplica,wglTMDPDF_OPE_tw3_SetScaleVariation

!!!!!!----FOR TEST
!public::MakeGrid,ExtractFromGrid,CxF_compute,TestGrid

contains

!! Coefficient function
INCLUDE 'Code/wglTMDPDF/coeffFunc.f90'

!! Coefficient function
INCLUDE 'Code/wglTMDPDF/coeffFunc_largeX.f90'

!! Mellin convolution routine
INCLUDE 'Code/Twist2/Twist2_WW.f90'


function wglTMDPDF_OPE_IsInitialized()
    logical::wglTMDPDF_OPE_IsInitialized
    wglTMDPDF_OPE_IsInitialized=started
end function wglTMDPDF_OPE_IsInitialized

!! Initialization of the package
subroutine wglTMDPDF_OPE_Initialize(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=300)::path
    logical::initRequired
    character(len=8)::order_global
    integer::i,FILEver
    real(dp),allocatable::subGridsX(:),subGridsB(:)

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

    call MoveTO(51,'*16   ')
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
    read(51,*) withGluon
    withGluonTW3=withGluon
    call MoveTO(51,'*p2  ')
    read(51,*) numberOfHadrons

    !----- ORDER
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
        CASE DEFAULT
            if(outputLevel>0)write(*,*) &
                WarningString('Initialize: unknown order for coefficient function. Switch to NLO.',moduleName)
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NLO'
            orderMain=1
        END SELECT

    if(outputLevel>2 .and. orderMain>-1) write(*,'(A,I1)') ' |  Coef.func.    =as^',orderMain

    call MoveTO(51,'*p2  ')
    read(51,*) useGrid
    call MoveTO(51,'*p3  ')
    read(51,*) runTest

    call MoveTO(51,'*p4  ')
    read(51,*) resumLargeX
    if(resumLargeX) then
        call MoveTO(51,'*p5  ')
        read(51,*) order_global

        SELECT CASE(trim(order_global))
            CASE ("LO")
                if(outputLevel>1) write(*,*) trim(moduleName)//' Large-X order set: LO'
                orderLX =0
            CASE ("NLO")
                if(outputLevel>1) write(*,*) trim(moduleName)//' Large-X order set: NLO'
                orderLX=1
            CASE ("NNLO")
                if(outputLevel>1) write(*,*) trim(moduleName)//' Large-X order set: NNLO'
                orderLX=2
            CASE ("N2LO")
                if(outputLevel>1) write(*,*) trim(moduleName)//' Large-X order set: NNLO'
                orderLX=2
            CASE ("NNNLO")
                if(outputLevel>1) write(*,*) trim(moduleName)//' Large-X order set: N3LO'
                orderLX=3
            CASE ("N3LO")
                if(outputLevel>1) write(*,*) trim(moduleName)//' Large-X order set: N3LO'
                orderLX=3
            CASE DEFAULT
                if(outputLevel>0) then
                    write(*,*) &
                    WarningString('Initialize: unknown order for large-X resummation of coefficient function.',moduleName)
                    write(*,*) WarningString('set to same order as the common part.',moduleName)
                end if
                orderLX=orderMain
            END SELECT

    if(outputLevel>2 .and. orderLX>-1) write(*,'(A,I1)') ' |  Large-X       =as^',orderLX
    end if

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

    !-------------Parameters of grid
    call MoveTO(51,'*E   ')
    call MoveTO(51,'*p1  ')
    read(51,*) i
    allocate(subGridsX(0:i))
    call MoveTO(51,'*p2  ')
    read(51,*) subGridsX
    call MoveTO(51,'*p4  ')
    read(51,*) i
    allocate(subGridsB(0:i))
    call MoveTO(51,'*p5  ')
    read(51,*) subGridsB

    !-------------------tw3-part is place holder for a moment
    call MoveTO(51,'*F   ')
    call MoveTO(51,'*p1  ')
    read(51,*) order_global

    SELECT CASE(trim(order_global))
        CASE ("NA")
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order[tw3] set: NA',color(" (TMD=fNP)",c_yellow)
            orderMainTW3=-50
        CASE DEFAULT
            if(outputLevel>0)write(*,*) &
                WarningString('Initialize: the perturbative tw3 convolution is not implemented yet. Switch to NA.',moduleName)
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NA'
            orderMainTW3=-50
        END SELECT

    if(outputLevel>2 .and. orderMainTW3>-1) write(*,'(A,I1)') ' |  Coef.func.    =as^',orderMainTW3

    !!!! ignored
    call MoveTO(51,'*p2  ')
    read(51,*) useGridTW3
    call MoveTO(51,'*p3  ')
    read(51,*) runTestTW3

    CLOSE (51, STATUS='KEEP')
    c4_global=1d0


    xMin=subGridsX(0)
    bMin=subGridsB(0)
    bMax=subGridsB(size(subGridsB)-1)

    if(abs(subGridsX(size(subGridsX)-1)-1)>toleranceGEN) then
        write(*,*) ErrorString("The last subgrid in X must complete by x=1. Initialization terminated",moduleName)
        stop
    end if

    call Twist2_ChGrid_Initialize(path,'*13  ','*E   ',numberOfHadrons,withGluon,moduleName,outputLevel)
    
    !!! Model initialisation is called from the wglTMDPDF-module
    
    !!!!!!!Checking the x-dependance of muOPE
    IsMuYdependent=testMU()

    if(IsMuYdependent) then
        if(outputLevel>2) write(*,*) trim(moduleName)//': mu OPE is dependent on x'
    else
        if(outputLevel>2) write(*,*) trim(moduleName)//': mu OPE is independent on x'
    end if
    gridReady=.false.

    if(useGrid) then
        if(resumLargeX) then
            call Twist2_ChGrid_MakeGrid(CxF_largeX_compute)
            gridReady=.true.
            if(runTest) call TestGrid(CxF_largeX_compute)
        else
            call Twist2_ChGrid_MakeGrid(CxF_compute)
            gridReady=.true.
            if(runTest) call TestGrid(CxF_compute)
        end if
    end if

    started=.true.
    messageCounter=0

    if(outputLevel>0) write(*,*) color('----- arTeMiDe.wglTMDPDF_OPE '//trim(version)//': .... initialized',c_green)
    if(outputLevel>1) write(*,*) ' '

end subroutine wglTMDPDF_OPE_Initialize

!!!!!!!--------------------------- DEFINING ROUTINES ------------------------------------------

!!!!array of x times PDF(x,Q) for hadron 'hadron'
!!!! array is (-5:5) (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
function xf(x,Q,hadron)
    real(dp) :: x,Q
    integer:: hadron
    real(dp), dimension(-5:5):: xf
    
    xf=x_hPDF(x,Q,hadron)
    
end function xf


function wglTMDPDF_OPE_convolution(x,b,h,addGluon)
    real(dp),dimension(-5:5)::wglTMDPDF_OPE_convolution
    real(dp),intent(in)::x,b
    integer,intent(in)::h
    logical,optional,intent(in)::addGluon

    logical::gluon

    !!! check gluonity
    if(present(addGluon)) then
        gluon=addGluon
    else
        gluon=withGluon
    end if

    !!!! test for boundaries is done in wglTMDPDF_lowScale5 (on the enty to this procedure)

    !!!! case NA
    if(orderMain==-50) then
        if(gluon) then
            wglTMDPDF_OPE_convolution=1._dp
        else
            wglTMDPDF_OPE_convolution=(/1._dp,1._dp,1._dp,1._dp,1._dp,0._dp,1._dp,1._dp,1._dp,1._dp,1._dp/)
        end if
        return
    end if

    !!!! The factor 1/x is replace by 1 because the WW part of wgl is
    !!!! -x^2\int dy/y^2 c[x/y]h[y] = x \int dy/y (-x/y c[x/y] h[y])
    !!!! so x is canceled, and the factor x/y is included in the coefficient function

    !!! computation
    if(useGrid) then
        if(gridReady) then
            wglTMDPDF_OPE_convolution=ExtractFromGrid(x,b,h)
        else
            call Warning_Raise('Called OPE_convolution while grid is not ready.',messageCounter,messageTrigger,moduleName)
            call wglTMDPDF_OPE_resetGrid()
            wglTMDPDF_OPE_convolution=ExtractFromGrid(x,b,h)
        end if
    else
        if(resumLargeX) then
            wglTMDPDF_OPE_convolution=CxF_largeX_compute(x,b,h,gluon)
        else
            wglTMDPDF_OPE_convolution=CxF_compute(x,b,h,gluon)
        end if
    end if

end function wglTMDPDF_OPE_convolution

!!! this is convolution with twist3 PDF
function wglTMDPDF_OPE_tw3_convolution(x,b,h,addGluon)
    real(dp),dimension(-5:5)::wglTMDPDF_OPE_tw3_convolution
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

    !!!! test for boundaries is done in wglTMDPDF_lowScale5 (on the enty to this procedure)

    !!!! case NA
    if(orderMainTW3==-50) then
        if(gluon) then
            wglTMDPDF_OPE_tw3_convolution=1._dp
        else
            wglTMDPDF_OPE_tw3_convolution=(/1._dp,1._dp,1._dp,1._dp,1._dp,0._dp,1._dp,1._dp,1._dp,1._dp,1._dp/)
        end if
        return
    end if

    !!!!! perturbative convolution is not implemented yet

end function wglTMDPDF_OPE_tw3_convolution


!!!!!!!!!! ------------------------ SUPPORINTG ROUTINES --------------------------------------
!!! This subroutine force reconstruction of the grid (if griding is ON)
subroutine wglTMDPDF_OPE_resetGrid()
    gridReady=.false.
    if(useGrid) then
        if(outputLevel>1) write(*,*) 'arTeMiDe ',moduleName,':  Grid Reset. with c4=',c4_global

        if(resumLargeX) then
            call Twist2_ChGrid_MakeGrid(CxF_largeX_compute)
        else
            call Twist2_ChGrid_MakeGrid(CxF_compute)
        end if

        gridReady=.true.
    end if
end subroutine wglTMDPDF_OPE_resetGrid

!! call QCDinput to change the PDF replica number
!! unset the grid, since it should be recalculated fro different PDF replica.
subroutine wglTMDPDF_OPE_SetPDFreplica(rep,hadron)
    integer,intent(in):: rep,hadron
    logical::newPDF

    call QCDinput_SethPDFreplica(rep,hadron,newPDF)
    if(newPDF) then
        gridReady=.false.
        call wglTMDPDF_OPE_resetGrid()
    else
        if(outputLevel>1) write(*,"('arTeMiDe ',A,':  replica of PDF (',I4,') is the same as the used one. Nothing is done!')") &
        moduleName, rep
    end if

end subroutine wglTMDPDF_OPE_SetPDFreplica

!!!! this routine set the variations of scales
!!!! it is used for the estimation of errors
subroutine wglTMDPDF_OPE_SetScaleVariation(c4_in)
    real(dp),intent(in)::c4_in
    if(c4_in<0.1d0 .or. c4_in>10.d0) then
        if(outputLevel>0) write(*,*) WarningString('variation in c4 is enourmous. c4 is set to 2',moduleName)
        c4_global=2d0
        call wglTMDPDF_OPE_resetGrid()
    else if(abs(c4_in-c4_global)<toleranceGEN) then
        if(outputLevel>1) write(*,*) color('wglTMDPDF: c4-variation is ignored. c4='//real8ToStr(c4_global),c_yellow)
    else
        c4_global=c4_in
        if(outputLevel>1) write(*,*) color('wglTMDPDF: set scale variations c4 as:'//real8ToStr(c4_global),c_yellow)
        call wglTMDPDF_OPE_resetGrid()
    end if
end subroutine wglTMDPDF_OPE_SetScaleVariation

!!! This subroutine force reconstruction of the grid (if griding is ON)
subroutine wglTMDPDF_OPE_tw3_resetGrid()
    !!!!! not implemented yet
end subroutine wglTMDPDF_OPE_tw3_resetGrid

!!! This subroutine force reconstruction of the grid (if griding is ON)
subroutine wglTMDPDF_OPE_tw3_testGrid()
    !!!!! not implemented yet
end subroutine wglTMDPDF_OPE_tw3_testGrid

!! call QCDinput to change the PDF replica number
!! unset the grid, since it should be recalculated fro different PDF replica.
subroutine wglTMDPDF_OPE_tw3_SetPDFreplica(rep,hadron)
    integer,intent(in):: rep,hadron

    !!!! not implemented yet

end subroutine wglTMDPDF_OPE_tw3_SetPDFreplica

!!!! this routine set the variations of scales
!!!! it is used for the estimation of errors
subroutine wglTMDPDF_OPE_tw3_SetScaleVariation(c4_in)
    real(dp),intent(in)::c4_in
    !!!! not implemented yet
end subroutine wglTMDPDF_OPE_tw3_SetScaleVariation

end module wglTMDPDF_OPE
