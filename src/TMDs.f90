!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 2.0
!
!	Evaluation of the TMDs
!	
!	if you use this module please, quote 1706.01473
!
!	ver 1.0: release (AV, 10.05.2017)
!	ver 1.4: release (AV, 23.12.2018)
!	ver 1.41:release (AV, 23.12.2018)
!	ver 2.00:release (AV, 29.03.2019)
!	ver 2.01:release (AV, 08.06.2019)
!	ver 2.05:release (AV, 22.05.2020)
!
!				A.Vladimirov (23.12.2018)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDs
use aTMDe_Numerics
use IO_functions
use QCDinput
use TMDR
use uTMDPDF
use uTMDFF
use lpTMDPDF
use SiversTMDPDF
use wgtTMDPDF
implicit none

private
!   public

character (len=7),parameter :: moduleName="TMDs"
character (len=5),parameter :: version="v2.06"
!Last appropriate verion of constants-file
integer,parameter::inputver=19

!------------------------------------------Physical and mathematical constants------------------------------------------
!------------------------------------------Working variables------------------------------------------------------------

logical::started=.false.
integer::outputLevel=2
integer::messageTrigger=5

!!!! The TMD evolution can be tritted differently
integer::EvolutionType
logical::include_uTMDPDF
logical::include_uTMDFF
logical::include_lpTMDPDF
logical::include_SiversTMDPDF
logical::include_wgtTMDPDF


!!!parameters for the uncertanty estimation
real(dp)::c1_global,c3_global

!-----------------------------------------Public interface--------------------------------------------------------------
public::TMDs_SetScaleVariations,TMDs_Initialize,TMDs_IsInitialized
real(dp),dimension(-5:5),public:: uTMDPDF_5,uTMDPDF_50
real(dp),dimension(-5:5),public:: uTMDFF_5,uTMDFF_50
real(dp),dimension(-5:5),public:: lpTMDPDF_50
real(dp),dimension(-5:5),public:: SiversTMDPDF_5,SiversTMDPDF_50
real(dp),dimension(-5:5),public:: wgtTMDPDF_5,wgtTMDPDF_50

public::uPDF_uPDF,uPDF_anti_uPDF

interface uTMDPDF_5
    module procedure uTMDPDF_5_Ev,uTMDPDF_5_optimal
end interface

interface uTMDPDF_50
    module procedure uTMDPDF_50_Ev,uTMDPDF_50_optimal
end interface

interface uTMDFF_5
    module procedure uTMDFF_5_Ev,uTMDFF_5_optimal
end interface

interface uTMDFF_50
    module procedure uTMDFF_50_Ev,uTMDFF_50_optimal
end interface

interface lpTMDPDF_50
    module procedure lpTMDPDF_50_Ev,lpTMDPDF_50_optimal
end interface

interface SiversTMDPDF_5
    module procedure SiversTMDPDF_5_Ev,SiversTMDPDF_5_optimal
end interface

interface SiversTMDPDF_50
    module procedure SiversTMDPDF_50_Ev,SiversTMDPDF_50_optimal
end interface

interface wgtTMDPDF_5
    module procedure wgtTMDPDF_5_Ev,wgtTMDPDF_5_optimal
end interface

interface wgtTMDPDF_50
    module procedure wgtTMDPDF_50_Ev,wgtTMDPDF_50_optimal
end interface

contains 

INCLUDE 'Model/TMDs_model.f90'
INCLUDE 'Code/TMDs/TMD-calls.f90'

function TMDs_IsInitialized()
    logical::TMDs_IsInitialized
    TMDs_IsInitialized=started
end function TMDs_IsInitialized

!!! Initializing routing
!!! Filles the prebuiled arrays
!!! orderAD, is order of anomalous dimension for evolution
!!!! order zeta is the order of pertrubative zeta expression, used only in the "universal TMD"
subroutine TMDs_Initialize(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=300)::path,line
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
        stop
    end if
    call MoveTO(51,'*p2  ')
    read(51,*) outputLevel    
    if(outputLevel>2) write(*,*) '--------------------------------------------- '
    if(outputLevel>2) write(*,*) 'artemide.',moduleName,version,': initialization started ... '
    call MoveTO(51,'*p3  ')
    read(51,*) messageTrigger


    call MoveTO(51,'*6   ')
    call MoveTO(51,'*p1  ')
    read(51,*) initRequired
    if(.not.initRequired) then
        if(outputLevel>2) write(*,*)'artemide.',moduleName,': initialization is not required. '
        started=.false.
        return
    end if
    CLOSE (51, STATUS='KEEP') 
    !!! then we read it again from the beginning to fill parameters of other modules
    OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")

    !! TMDR
    call MoveTO(51,'*3   ')
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p2  ')
    read(51,*) EvolutionType
    if(outputLevel>2) write(*,'(A,I3)') ' Evolution type =',EvolutionType

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

    CLOSE (51, STATUS='KEEP')

    if(.not.QCDinput_IsInitialized()) then
        if(outputLevel>1) write(*,*) '.. initializing QCDinput (from ',moduleName,')'
        if(present(prefix)) then
            call QCDinput_Initialize(file,prefix)
        else
            call QCDinput_Initialize(file)
        end if
    end if

    if(.not.TMDR_IsInitialized()) then
        if(outputLevel>1) write(*,*) '.. initializing TMDR (from ',moduleName,')'
        if(present(prefix)) then
            call TMDR_Initialize(file,prefix)
        else
            call TMDR_Initialize(file)
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

    c1_global=1d0
    c3_global=1d0
    
    
    started=.true.
    if(outputLevel>0) write(*,*) color('----- arTeMiDe.TMDs '//trim(version)//': .... initialized',c_green)
    if(outputLevel>1) write(*,*) 
end subroutine TMDs_Initialize

!!!! this routine set the variations of scales
!!!! it is used for the estimation of errors
subroutine TMDs_SetScaleVariations(c1_in,c3_in)
    real(dp),intent(in)::c1_in,c3_in

    If(EvolutionType==1) then
        if(c1_in<0.1d0 .or. c1_in>10.d0) then
            if(outputLevel>0) write(*,*) WarningString('variation in c1 is enourmous. c1 is set to 2',moduleName)
            c1_global=2d0
        else
            c1_global=c1_in
        end if
    else
        if(ABS(c1_in-1d0)>0.00001d0) then
            if(outputLevel>0) write(*,*) WarningString('variation of c1 is sensless. &
            There is no such constant within current evolution scheme. Nothing is done',moduleName)
        end if
    end if

    if(EvolutionType==3) then
        if(ABS(c3_in-1d0)>0.00001d0) then
            if(outputLevel>0) write(*,*) WarningString("variation of c3 is sensless. &
            There is no such constant within current evolution scheme. Nothing is done",moduleName)
        end if
    else
    if(c3_in<0.1d0 .or. c3_in>10.d0) then
        if(outputLevel>0) write(*,*) WarningString('variation in c3 is enourmous. c3 is set to 2',moduleName)
            c3_global=2d0
        else
            c3_global=c3_in
        end if
    end if

    if(outputLevel>1) write(*,*) 'TMDs: set scale variations constants (c1,c3) as:',c1_global,c3_global
end subroutine TMDs_SetScaleVariations 

end module TMDs
