!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe-control 2.0
!
!	Module that control the flow of artemide package.
!	Initialize and reset submodules
!	
!	if you use this module please, quote 1803.11089
!
!				A.Vladimirov (30.05.2019)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module aTMDe_control
use aTMDe_Numerics
use IO_functions
use QCDinput
use EWinput
use uTMDPDF
use lpTMDPDF
use SiversTMDPDF
use wgtTMDPDF
use uTMDFF
use TMDR
use TMDs
use TMDF
use TMDX_DY
use TMDX_SIDIS
use TMDs_inKT
use aTMDe_Setup
implicit none

private
character (len=14),parameter :: moduleName="aTMDe-control"
character (len=5),parameter :: version="v2.06"
!Last appropriate verion of constants-file
integer,parameter::inputver=11
character (len=15),parameter :: constNAME="aTMDe-temporary"

integer::outputLevel=2
integer::messageTrigger=5

!!!! indicator of the aTMDe_control intialization
logical::isStarted=.false.
!!!! indicators of modules usage  
logical::include_EWinput,include_uTMDPDF,include_uTMDFF,include_TMDR,include_TMDs,include_TMDF
logical::include_lpTMDPDF,include_SiversTMDPDF,include_wgtTMDPDF
logical::include_TMDX_DY,include_TMDX_SIDIS,include_TMDs_inKT

!!!! legths of non-perturbative arrays
integer::NPlength_total
integer::NPlength_TMDR,NPlength_uTMDPDF,NPlength_uTMDFF,NPlength_lpTMDPDF,NPlength_SiversTMDPDF,NPlength_wgtTMDPDF

!!!! non-pertrubative parameters for individual modules
real(dp),allocatable::lambdaNP_TMDR(:),lambdaNP_uTMDPDF(:),lambdaNP_uTMDFF(:),lambdaNP_lpTMDPDF(:),&
            lambdaNP_SiversTMDPDF(:),lambdaNP_wgtTMDPDF(:)

!!!! Saved values of scale-variation parameters
real(dp)::c1_saved,c2_saved,c3_saved,c4_saved

public::artemide_Initialize
public::artemide_SetNPparameters,artemide_SetNPparameters_TMDR,artemide_SetNPparameters_uTMDFF,artemide_SetNPparameters_uTMDPDF
public::artemide_SetNPparameters_lpTMDPDF,artemide_SetNPparameters_SiversTMDPDF,artemide_SetNPparameters_wgtTMDPDF
public::artemide_SetReplica_TMDR,artemide_SetReplica_uTMDFF,artemide_SetReplica_uTMDPDF
public::artemide_SetReplica_lpTMDPDF,artemide_SetReplica_SiversTMDPDF,artemide_SetReplica_wgtTMDPDF
public::artemide_SetScaleVariations
public::artemide_ShowStatistics
public::artemide_GetReplicaFromFile,artemide_NumOfReplicasInFile

contains
  
!! Module that creates the constants-file,
!! and the read it
!! intialize it-self and artemide modules
subroutine artemide_Initialize(file,prefix,order)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=*),optional::order
    character(len=300)::path
    integer::FILEver,i
    logical::initilize_NParrays
    !-----------------------------------------------------------
    if(present(prefix).and.present(order)) then
        call artemide_Setup_fromFile(file,prefix=prefix,order=order)
    else if(present(prefix)) then
        call artemide_Setup_fromFile(file,prefix=prefix)
    else if(present(order)) then
        call artemide_Setup_fromFile(file,order=order)
    else
        call artemide_Setup_fromFile(file)
    end if

    if(present(prefix)) then
        call CreateConstantsFile(constNAME,prefix)
    else
        call CreateConstantsFile(constNAME)
    end if
    !-----------------------------------------------------------

    if(present(prefix)) then
        path=trim(adjustl(prefix))//trim(adjustr(constNAME))
    else
        path=trim(adjustr(constNAME))
    end if

    OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")
    !!! Search for output level
    call MoveTO(51,'*0   ')
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) FILEver
    if(FILEver<inputver) then
        write(*,*) color('artemide.'//trim(moduleName)//': const-file version is too old.',c_red_bold)
        write(*,*) color('		     Update the const-file with artemide.setup',c_red_bold)
        write(*,*) '  '
        stop
    end if
    call MoveTO(51,'*p2  ')
    read(51,*) outputLevel    
    if(outputLevel>2) write(*,*) '--------------------------------------------- '
    if(outputLevel>2) write(*,*) 'artemide.control: initialization started ... '
    call MoveTO(51,'*p3  ')
    read(51,*) outputLevel    

    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) initilize_NParrays

    call MoveTO(51,'*2   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_EWinput

    call MoveTO(51,'*3   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_TMDR
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) NPlength_TMDR
    !allocate lambda's and read initialization NP-array
    if(include_TMDR) then
        allocate(lambdaNP_TMDR(1:NPlength_TMDR))
        call MoveTO(51,'*p2  ')
        do i=1,NPlength_TMDR
            read(51,*) lambdaNP_TMDR(i)
        end do
    end if


    call MoveTO(51,'*4   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_uTMDPDF
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) NPlength_uTMDPDF
    !allocate lambda's and read initialization NP-array
    if(include_uTMDPDF) then
        allocate(lambdaNP_uTMDPDF(1:NPlength_uTMDPDF))
        call MoveTO(51,'*p2  ')
        do i=1,NPlength_uTMDPDF
            read(51,*) lambdaNP_uTMDPDF(i)
        end do
    end if

    call MoveTO(51,'*5   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_uTMDFF
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) NPlength_uTMDFF
    !allocate lambda's and read initialization NP-array
    if(include_uTMDFF) then
        allocate(lambdaNP_uTMDFF(1:NPlength_uTMDFF))
        call MoveTO(51,'*p2  ')
        do i=1,NPlength_uTMDFF
            read(51,*) lambdaNP_uTMDFF(i)
        end do
    end if

    call MoveTO(51,'*6   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_TMDs

    call MoveTO(51,'*7   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_TMDF

    call MoveTO(51,'*8   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_TMDs_inKT

    call MoveTO(51,'*9   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_TMDX_DY

    call MoveTO(51,'*10   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_TMDX_SIDIS

    call MoveTO(51,'*11   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_lpTMDPDF
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) NPlength_lpTMDPDF
    !allocate lambda's and read initialization NP-array
    if(include_lpTMDPDF) then
        allocate(lambdaNP_lpTMDPDF(1:NPlength_lpTMDPDF))
        call MoveTO(51,'*p2  ')
        do i=1,NPlength_lpTMDPDF
            read(51,*) lambdaNP_lpTMDPDF(i)
        end do
    end if
    
    call MoveTO(51,'*12   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_SiversTMDPDF
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) NPlength_SiversTMDPDF
    !allocate lambda's and read initialization NP-array
    if(include_SiversTMDPDF) then
        allocate(lambdaNP_SiversTMDPDF(1:NPlength_SiversTMDPDF))
        call MoveTO(51,'*p2  ')
        do i=1,NPlength_SiversTMDPDF
            read(51,*) lambdaNP_SiversTMDPDF(i)
        end do
    end if
    
    call MoveTO(51,'*13   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_wgtTMDPDF
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) NPlength_wgtTMDPDF
    !allocate lambda's and read initialization NP-array
    if(include_wgtTMDPDF) then
        allocate(lambdaNP_wgtTMDPDF(1:NPlength_wgtTMDPDF))
        call MoveTO(51,'*p2  ')
        do i=1,NPlength_wgtTMDPDF
            read(51,*) lambdaNP_wgtTMDPDF(i)
        end do
    end if

    CLOSE (51, STATUS='KEEP')
    !-----------------------------------------------------------
    if(outputLevel>2) write(*,*) 'artemide.control: initialization file is read. Initialization of modules ... '

    if(present(prefix)) then
        call QCDinput_Initialize(constNAME,prefix)
    else
        call QCDinput_Initialize(constNAME)
    end if

    if(include_EWinput) then
        if(present(prefix)) then
            call EWinput_Initialize(constNAME,prefix)
        else
            call EWinput_Initialize(constNAME)
        end if
    end if

    if(include_uTMDPDF) then
        if(present(prefix)) then
            call uTMDPDF_Initialize(constNAME,prefix)
        else
            call uTMDPDF_Initialize(constNAME)
        end if
    end if

    if(include_uTMDFF) then
        if(present(prefix)) then
            call uTMDFF_Initialize(constNAME,prefix)
        else
            call uTMDFF_Initialize(constNAME)
        end if
    end if

    if(include_lpTMDPDF) then
        if(present(prefix)) then
            call lpTMDPDF_Initialize(constNAME,prefix)
        else
            call lpTMDPDF_Initialize(constNAME)
        end if
    end if
    
    if(include_SiversTMDPDF) then
        if(present(prefix)) then
            call SiversTMDPDF_Initialize(constNAME,prefix)
        else
            call SiversTMDPDF_Initialize(constNAME)
        end if
    end if
    
    if(include_wgtTMDPDF) then
        if(present(prefix)) then
            call wgtTMDPDF_Initialize(constNAME,prefix)
        else
            call wgtTMDPDF_Initialize(constNAME)
        end if
    end if

    if(include_TMDR) then
        if(present(prefix)) then
            call TMDR_Initialize(constNAME,prefix)
        else
            call TMDR_Initialize(constNAME)
        end if
    end if

    if(include_TMDs) then
        if(present(prefix)) then
            call TMDs_Initialize(constNAME,prefix)
        else
            call TMDs_Initialize(constNAME)
        end if
    end if

    if(include_TMDF) then
        if(present(prefix)) then
            call TMDF_Initialize(constNAME,prefix)
        else
            call TMDF_Initialize(constNAME)
        end if
    end if

    if(include_TMDs_inKT) then
        if(present(prefix)) then
            call TMDs_inKT_Initialize(constNAME,prefix)
        else
            call TMDs_inKT_Initialize(constNAME)
        end if
    end if

    if(include_TMDX_DY) then
        if(present(prefix)) then
            call TMDX_DY_Initialize(constNAME,prefix)
        else
            call TMDX_DY_Initialize(constNAME)
        end if
    end if

    if(include_TMDX_SIDIS) then
        if(present(prefix)) then
            call TMDX_SIDIS_Initialize(constNAME,prefix)
        else
            call TMDX_SIDIS_Initialize(constNAME)
        end if
    end if

    !------------------------------------------------------
    !!! calculate total NP lendth
    NPlength_total=0
    if(include_TMDR) NPlength_total=NPlength_total+NPlength_TMDR
    if(include_uTMDPDF) NPlength_total=NPlength_total+NPlength_uTMDPDF
    if(include_uTMDFF) NPlength_total=NPlength_total+NPlength_uTMDFF
    if(include_lpTMDPDF) NPlength_total=NPlength_total+NPlength_lpTMDPDF
    if(include_SiversTMDPDF) NPlength_total=NPlength_total+NPlength_SiversTMDPDF
    if(include_wgtTMDPDF) NPlength_total=NPlength_total+NPlength_wgtTMDPDF

    if(outputLevel>2) write(*,*) ' artemide.control: Total number of NP parameters:',NPlength_total

    c1_saved=1d0
    c2_saved=1d0
    c3_saved=1d0
    c4_saved=1d0

    call artemide_SetScaleVariations(c1_saved,c2_saved,c3_saved,c4_saved)

    isStarted=.true.

    if(outputLevel>1) write(*,*) color(' artemide.control: initialization done.',c_green_bold)

    if(initilize_NParrays) then
    if(outputLevel>1) write(*,*) 'artemide.control: setting initial NP arrays.'
    if(include_TMDR) call artemide_SetNPparameters_TMDR(lambdaNP_TMDR)
    if(include_uTMDPDF) call artemide_SetNPparameters_uTMDPDF(lambdaNP_uTMDPDF)
    if(include_uTMDFF) call artemide_SetNPparameters_uTMDFF(lambdaNP_uTMDFF)
    if(include_lpTMDPDF) call artemide_SetNPparameters_lpTMDPDF(lambdaNP_lpTMDPDF)
    if(include_SiversTMDPDF) call artemide_SetNPparameters_SiversTMDPDF(lambdaNP_SiversTMDPDF)
    if(include_wgtTMDPDF) call artemide_SetNPparameters_wgtTMDPDF(lambdaNP_wgtTMDPDF)
    if(outputLevel>1) write(*,*) color(' artemide.control: initial NP arrays set.',c_green_bold)
    end if
end subroutine artemide_Initialize
  
  !-------------------------------------------------------------- NP parameters ---------------------------
subroutine artemide_SetNPparameters(lambdaNP)
    real(dp),intent(in)::lambdaNP(:)
    real(dp),dimension(1:NPlength_total)::lambda_cur
    integer::ll,num

    ll=size(lambdaNP)

    if(ll<NPlength_total) then
        if(outputLevel>0) write(*,"(A,I4,A,I4,A)") &
            color('artemide.SetNPparameters: ERROR: the length of NP parameters array (',c_red),ll,&
            color(') is smaller then the total number of NP parameters (',c_red),NPlength_total,color(')',c_red)
        if(outputLevel>0) write(*,*) color('NOTHING IS DONE',c_red)
    return
    end if

    if(ll>NPlength_total) then
        if(outputLevel>0) write(*,"(A,I4,A,I4,A)")&
            color('artemide.SetNPparameters: ERROR: the length of NP parameters array (',c_red), ll,&
            color(') is larger then total the number of NP parameters (',c_red),NPlength_total,color(')',c_red)
        if(outputLevel>0) write(*,*) color('The array is trucated',c_red)
    end if

    lambda_cur=lambdaNP(1:NPlength_total)

    num=0
    !!! split the incoming array
    if(include_TMDR) then
        lambdaNP_TMDR=lambda_cur(num+1:num+NPlength_TMDR)
        num=num+NPlength_TMDR
    end if
    if(include_uTMDPDF) then
        lambdaNP_uTMDPDF=lambda_cur(num+1:num+NPlength_uTMDPDF)
        num=num+NPlength_uTMDPDF
    end if
    if(include_uTMDFF) then
        lambdaNP_uTMDFF=lambda_cur(num+1:num+NPlength_uTMDFF)
        num=num+NPlength_uTMDFF
    end if
    if(include_lpTMDPDF) then
        lambdaNP_lpTMDPDF=lambda_cur(num+1:num+NPlength_lpTMDPDF)
        num=num+NPlength_uTMDPDF
    end if
    if(include_SiversTMDPDF) then
        lambdaNP_SiversTMDPDF=lambda_cur(num+1:num+NPlength_SiversTMDPDF)
        num=num+NPlength_SiversTMDPDF
    end if
    if(include_wgtTMDPDF) then
        lambdaNP_wgtTMDPDF=lambda_cur(num+1:num+NPlength_wgtTMDPDF)
        num=num+NPlength_wgtTMDPDF
    end if

    !!! sending NP arrays to packages
    if(include_TMDR) call TMDR_setNPparameters(lambdaNP_TMDR)
    if(include_uTMDPDF) call uTMDPDF_SetLambdaNP(lambdaNP_uTMDPDF)
    if(include_uTMDFF) call uTMDFF_SetLambdaNP(lambdaNP_uTMDFF)
    if(include_lpTMDPDF) call lpTMDPDF_SetLambdaNP(lambdaNP_lpTMDPDF)
    if(include_SiversTMDPDF) call SiversTMDPDF_SetLambdaNP(lambdaNP_SiversTMDPDF)
    if(include_wgtTMDPDF) call wgtTMDPDF_SetLambdaNP(lambdaNP_wgtTMDPDF)

    !!! reseting other packages
    if(include_TMDF) call TMDF_ResetCounters()
    if(include_TMDX_DY) call TMDX_DY_ResetCounters()
    if(include_TMDX_SIDIS) call TMDX_SIDIS_ResetCounters()

end subroutine artemide_SetNPparameters
  
subroutine artemide_SetNPparameters_TMDR(lambdaNP)
    real(dp),intent(in)::lambdaNP(:)
    integer::ll

    if(.not.include_TMDR) then
        if(outputLevel>0) &
            write(*,*) &
            ErrorString('attempt to set NP-parameters for TMDR, while TMDR module is not included in the current setup',moduleName)
        if(outputLevel>0) write(*,*) color('NOTHING IS DONE',c_red)
        return
    end if

    ll=size(lambdaNP)


    if(ll<NPlength_TMDR) then
        if(outputLevel>0) write(*,"(A,I4,A,I4,A)")&
            color('artemide.SetNPparameters-TMDR: ERROR: the length of NP parameters array (',c_red),ll,&
            color(') is smaller then the total number of NP parameters for TMDR (',c_red),NPlength_TMDR,&
            color(')',c_red)
        if(outputLevel>0) write(*,*) color('NOTHING IS DONE',c_red)
        return
    end if

    if(ll>NPlength_TMDR) then
        if(outputLevel>0) write(*,"(A,I4,A,I4,A)")&
            color('artemide.SetNPparameters: ERROR: the length of NP parameters array (',c_red),ll,&
            color(') is larger then total the number of NP parameters for TMDR (',c_red),NPlength_TMDR,&
            color(')',c_red)
        if(outputLevel>0) write(*,*) color('The array is trucated',c_red)
    end if

    lambdaNP_TMDR=lambdaNP(1:NPlength_TMDR)

    call TMDR_setNPparameters(lambdaNP_TMDR)

    !!! reseting other packages
    if(include_TMDF) call TMDF_ResetCounters()
    if(include_TMDX_DY) call TMDX_DY_ResetCounters()
    if(include_TMDX_SIDIS) call TMDX_SIDIS_ResetCounters()

end subroutine artemide_SetNPparameters_TMDR
  
subroutine artemide_SetNPparameters_uTMDPDF(lambdaNP)
    real(dp),intent(in)::lambdaNP(:)
    integer::ll

    if(.not.include_uTMDPDF) then
        if(outputLevel>0) &
            write(*,*) &
            ErrorString('attempt to set NP-parameters for uTMDPDF, while uTMDPDF module is not included in the current setup',&
            moduleName)
        if(outputLevel>0) write(*,*) color('NOTHING IS DONE',c_red)
        return
    end if

    ll=size(lambdaNP)


    if(ll<NPlength_uTMDPDF) then
        if(outputLevel>0) write(*,"(A,I4,A,I4,A)")&
            color('artemide.SetNPparameters-uTMDPDF: ERROR: the length of NP parameters array (',c_red),ll,&
            color(') is smaller then the total number of NP parameters for uTMDPDF (',c_red),NPlength_uTMDPDF,&
            color(')',c_red)
        if(outputLevel>0) write(*,*) color('NOTHING IS DONE',c_red)
        return
    end if

    if(ll>NPlength_uTMDPDF) then
        if(outputLevel>0) write(*,"(A,I4,A,I4,A)")&
            color('artemide.SetNPparameters: ERROR: the length of NP parameters array (',c_red),ll,&
            color(') is larger then total the number of NP parameters for uTMDPDF (',c_red),NPlength_uTMDPDF,&
            color(')',c_red)
        if(outputLevel>0) write(*,*) color('The array is trucated',c_red)
    end if

    lambdaNP_uTMDPDF=lambdaNP(1:NPlength_uTMDPDF)

    call uTMDPDF_SetLambdaNP(lambdaNP_uTMDPDF)

    !!! reseting other packages
    if(include_TMDF) call TMDF_ResetCounters()
    if(include_TMDX_DY) call TMDX_DY_ResetCounters()
    if(include_TMDX_SIDIS) call TMDX_SIDIS_ResetCounters()

end subroutine artemide_SetNPparameters_uTMDPDF
  
subroutine artemide_SetNPparameters_uTMDFF(lambdaNP)
    real(dp),intent(in)::lambdaNP(:)
    integer::ll

    if(.not.include_uTMDFF) then
        if(outputLevel>0) &
        write(*,*) ErrorString('attempt to set NP-parameters for uTMDFF, while uTMDFF module is not included in the current setup',&
        moduleName)
        if(outputLevel>0) write(*,*) color('NOTHING IS DONE',c_red)
        return
    end if

    ll=size(lambdaNP)


    if(ll<NPlength_uTMDFF) then
        if(outputLevel>0) write(*,"(A,I4,A,I4,A)")&
            color('artemide.SetNPparameters-uTMDFF: ERROR: the length of NP parameters array (',c_red),ll,&
            color(') is smaller then the total number of NP parameters for uTMDFF (',c_red),NPlength_uTMDFF,&
            color(')',c_red)
        if(outputLevel>0) write(*,*) color('NOTHING IS DONE',c_red)
        return
    end if

    if(ll>NPlength_uTMDFF) then
        if(outputLevel>0) write(*,"(A,I4,A,I4,A)")&
            color('artemide.SetNPparameters: ERROR: the length of NP parameters array (',c_red),ll,&
            color(') is larger then total the number of NP parameters for uTMDFF (',c_red),NPlength_uTMDFF,&
            color(')',c_red)
        if(outputLevel>0) write(*,*) color('The array is trucated',c_red)
    end if

    lambdaNP_uTMDFF=lambdaNP(1:NPlength_uTMDFF)

    call uTMDFF_SetLambdaNP(lambdaNP_uTMDFF)

    !!! reseting other packages
    if(include_TMDF) call TMDF_ResetCounters()
    if(include_TMDX_DY) call TMDX_DY_ResetCounters()
    if(include_TMDX_SIDIS) call TMDX_SIDIS_ResetCounters()

end subroutine artemide_SetNPparameters_uTMDFF

   
subroutine artemide_SetNPparameters_lpTMDPDF(lambdaNP)
    real(dp),intent(in)::lambdaNP(:)
    integer::ll

    if(.not.include_lpTMDPDF) then
        if(outputLevel>0) &
            write(*,*) ErrorString('attempt to set NP-parameters for lpTMDPDF, &
        while lpTMDPDF module is not included in the current setup',moduleName)
        if(outputLevel>0) write(*,*) color('NOTHING IS DONE',c_red)
        return
    end if

    ll=size(lambdaNP)


    if(ll<NPlength_lpTMDPDF) then
        if(outputLevel>0) write(*,"(A,I4,A,I4,A)")&
            color('artemide.SetNPparameters-lpTMDPDF: ERROR: the length of NP parameters array (',c_red),ll,&
            color(') is smaller then the total number of NP parameters for lpTMDPDF (',c_red),NPlength_lpTMDPDF,&
            color(')',c_red)
        if(outputLevel>0) write(*,*) color('NOTHING IS DONE',c_red)
        return
    end if

    if(ll>NPlength_lpTMDPDF) then
        if(outputLevel>0) write(*,"(A,I4,A,I4,A)")&
            color('artemide.SetNPparameters: ERROR: the length of NP parameters array (',c_red),ll,&
            color(') is larger then total the number of NP parameters for lpTMDPDF (',c_red),NPlength_lpTMDPDF,&
            color(')',c_red)
        if(outputLevel>0) write(*,*) color('The array is trucated',c_red)
    end if

    lambdaNP_lpTMDPDF=lambdaNP(1:NPlength_lpTMDPDF)

    call lpTMDPDF_SetLambdaNP(lambdaNP_lpTMDPDF)

    !!! reseting other packages
    if(include_TMDF) call TMDF_ResetCounters()
    if(include_TMDX_DY) call TMDX_DY_ResetCounters()
    if(include_TMDX_SIDIS) call TMDX_SIDIS_ResetCounters()

end subroutine artemide_SetNPparameters_lpTMDPDF

subroutine artemide_SetNPparameters_SiversTMDPDF(lambdaNP)
    real(dp),intent(in)::lambdaNP(:)
    integer::ll

    if(.not.include_SiversTMDPDF) then
        if(outputLevel>0) &
            write(*,*) ErrorString('attempt to set NP-parameters for SiversTMDPDF, &
        while SiversTMDPDF module is not included in the current setup',moduleName)
        if(outputLevel>0) write(*,*) color('NOTHING IS DONE',c_red)
        return
    end if

    ll=size(lambdaNP)


    if(ll<NPlength_SiversTMDPDF) then
        if(outputLevel>0) write(*,"(A,I4,A,I4,A)")&
            color('artemide.SetNPparameters-SiversTMDPDF: ERROR: the length of NP parameters array (',c_red),ll,&
            color(') is smaller then the total number of NP parameters for SiversTMDPDF (',c_red),NPlength_SiversTMDPDF,&
            color(')',c_red)
        if(outputLevel>0) write(*,*) color('NOTHING IS DONE',c_red)
        return
    end if

    if(ll>NPlength_SiversTMDPDF) then
        if(outputLevel>0) write(*,"(A,I4,A,I4,A)")&
            color('artemide.SetNPparameters: ERROR: the length of NP parameters array (',c_red),ll,&
            color(') is larger then total the number of NP parameters for SiversTMDPDF (',c_red),NPlength_SiversTMDPDF,&
            color(')',c_red)
        if(outputLevel>0) write(*,*) color('The array is trucated',c_red)
    end if

    lambdaNP_SiversTMDPDF=lambdaNP(1:NPlength_SiversTMDPDF)

    call SiversTMDPDF_SetLambdaNP(lambdaNP_SiversTMDPDF)

    !!! reseting other packages
    if(include_TMDF) call TMDF_ResetCounters()
    if(include_TMDX_DY) call TMDX_DY_ResetCounters()
    if(include_TMDX_SIDIS) call TMDX_SIDIS_ResetCounters()

end subroutine artemide_SetNPparameters_SiversTMDPDF

subroutine artemide_SetNPparameters_wgtTMDPDF(lambdaNP)
    real(dp),intent(in)::lambdaNP(:)
    integer::ll

    if(.not.include_wgtTMDPDF) then
        if(outputLevel>0) &
            write(*,*) ErrorString('attempt to set NP-parameters for wgtTMDPDF, &
        while wgtTMDPDF module is not included in the current setup',moduleName)
        if(outputLevel>0) write(*,*) color('NOTHING IS DONE',c_red)
        return
    end if

    ll=size(lambdaNP)


    if(ll<NPlength_wgtTMDPDF) then
        if(outputLevel>0) write(*,"(A,I4,A,I4,A)")&
            color('artemide.SetNPparameters-wgtTMDPDF: ERROR: the length of NP parameters array (',c_red),ll,&
            color(') is smaller then the total number of NP parameters for wgtTMDPDF (',c_red),NPlength_wgtTMDPDF,&
            color(')',c_red)
        if(outputLevel>0) write(*,*) color('NOTHING IS DONE',c_red)
        return
    end if

    if(ll>NPlength_wgtTMDPDF) then
        if(outputLevel>0) write(*,"(A,I4,A,I4,A)")&
            color('artemide.SetNPparameters: ERROR: the length of NP parameters array (',c_red),ll,&
            color(') is larger then total the number of NP parameters for wgtTMDPDF (',c_red),NPlength_wgtTMDPDF,&
            color(')',c_red)
        if(outputLevel>0) write(*,*) color('The array is trucated',c_red)
    end if

    lambdaNP_wgtTMDPDF=lambdaNP(1:NPlength_wgtTMDPDF)

    call wgtTMDPDF_SetLambdaNP(lambdaNP_wgtTMDPDF)

    !!! reseting other packages
    if(include_TMDF) call TMDF_ResetCounters()
    if(include_TMDX_DY) call TMDX_DY_ResetCounters()
    if(include_TMDX_SIDIS) call TMDX_SIDIS_ResetCounters()

end subroutine artemide_SetNPparameters_wgtTMDPDF
 
subroutine artemide_SetReplica_TMDR(num)
    integer,intent(in)::num

    if(.not.include_TMDR) then
        if(outputLevel>0) &
        write(*,*) ErrorString('attempt to set a replica for TMDR,&
                while TMDR module is not included in the current setup',moduleName)
        if(outputLevel>0) write(*,*) color('NOTHING IS DONE',c_red)
        return
    end if

    call TMDR_setNPparameters(num)
    call TMDR_CurrentNPparameters(lambdaNP_TMDR)

    !!! reseting other packages
    if(include_TMDF) call TMDF_ResetCounters()
    if(include_TMDX_DY) call TMDX_DY_ResetCounters()
    if(include_TMDX_SIDIS) call TMDX_SIDIS_ResetCounters()

end subroutine artemide_SetReplica_TMDR
  
subroutine artemide_SetReplica_uTMDPDF(num)
    integer,intent(in)::num

    if(.not.include_uTMDPDF) then
        if(outputLevel>0) &
        write(*,*) ErrorString('attempt to set a replica for uTMDPDF,&
                while uTMDPDF module is not included in the current setup',moduleName)
        if(outputLevel>0) write(*,*) color('NOTHING IS DONE',c_red)
        return
    end if

    call uTMDPDF_SetLambdaNP(num)
    call uTMDPDF_CurrentNPparameters(lambdaNP_uTMDPDF)

    !!! reseting other packages
    if(include_TMDF) call TMDF_ResetCounters()
    if(include_TMDX_DY) call TMDX_DY_ResetCounters()
    if(include_TMDX_SIDIS) call TMDX_SIDIS_ResetCounters()

end subroutine artemide_SetReplica_uTMDPDF
  
subroutine artemide_SetReplica_uTMDFF(num)
    integer,intent(in)::num

    if(.not.include_uTMDFF) then
        if(outputLevel>0) &
        write(*,*) ErrorString('attempt to set a replica for uTMDFF,&
            while uTMDFF module is not included in the current setup',moduleName)
        if(outputLevel>0) write(*,*) color('NOTHING IS DONE',c_red)
        return
    end if

    call uTMDFF_SetLambdaNP(num)
    call uTMDFF_CurrentNPparameters(lambdaNP_uTMDFF)

    !!! reseting other packages
    if(include_TMDF) call TMDF_ResetCounters()
    if(include_TMDX_DY) call TMDX_DY_ResetCounters()
    if(include_TMDX_SIDIS) call TMDX_SIDIS_ResetCounters()

end subroutine artemide_SetReplica_uTMDFF
  
subroutine artemide_SetReplica_lpTMDPDF(num)
    integer,intent(in)::num

    if(.not.include_lpTMDPDF) then
        if(outputLevel>0) &
        write(*,*) ErrorString('attempt to set a replica for lpTMDPDF,&
            while lpTMDPDF module is not included in the current setup',ModuleName)
        if(outputLevel>0) write(*,*) 'NOTHING IS DONE'
        return
    end if

    call lpTMDPDF_SetLambdaNP(num)
    call lpTMDPDF_CurrentNPparameters(lambdaNP_lpTMDPDF)

    !!! reseting other packages
    if(include_TMDF) call TMDF_ResetCounters()
    if(include_TMDX_DY) call TMDX_DY_ResetCounters()
    if(include_TMDX_SIDIS) call TMDX_SIDIS_ResetCounters()

end subroutine artemide_SetReplica_lpTMDPDF

subroutine artemide_SetReplica_SiversTMDPDF(num)
    integer,intent(in)::num

    if(.not.include_SiversTMDPDF) then
        if(outputLevel>0) &
        write(*,*) ErrorString('attempt to set a replica for SiversTMDPDF,&
            while SiversTMDPDF module is not included in the current setup',ModuleName)
        if(outputLevel>0) write(*,*) 'NOTHING IS DONE'
        return
    end if

    call SiversTMDPDF_SetLambdaNP(num)
    call SiversTMDPDF_CurrentNPparameters(lambdaNP_SiversTMDPDF)

    !!! reseting other packages
    if(include_TMDF) call TMDF_ResetCounters()
    if(include_TMDX_DY) call TMDX_DY_ResetCounters()
    if(include_TMDX_SIDIS) call TMDX_SIDIS_ResetCounters()

end subroutine artemide_SetReplica_SiversTMDPDF

subroutine artemide_SetReplica_wgtTMDPDF(num)
    integer,intent(in)::num

    if(.not.include_wgtTMDPDF) then
        if(outputLevel>0) &
        write(*,*) ErrorString('attempt to set a replica for wgtTMDPDF,&
            while wgtTMDPDF module is not included in the current setup',ModuleName)
        if(outputLevel>0) write(*,*) 'NOTHING IS DONE'
        return
    end if

    call wgtTMDPDF_SetLambdaNP(num)
    call wgtTMDPDF_CurrentNPparameters(lambdaNP_wgtTMDPDF)

    !!! reseting other packages
    if(include_TMDF) call TMDF_ResetCounters()
    if(include_TMDX_DY) call TMDX_DY_ResetCounters()
    if(include_TMDX_SIDIS) call TMDX_SIDIS_ResetCounters()

end subroutine artemide_SetReplica_wgtTMDPDF
  
!------------------------------------------------------- Other routines ---------------------------------
subroutine artemide_ShowStatistics()
    integer::i
    write(*,*) '-------------------------  artemide statistics  -------------------------'
    write(*,*) 'Synopsis of NP input:'
    if(include_TMDR) then
        write(*,"('--- TMDR     : ',I3,' parameters')") NPlength_TMDR
        do i=1,NPlength_TMDR
            write(*,"(F10.4,' ')",advance='no') lambdaNP_TMDR(i)
        end do
        write(*,*) ' '
    end if
    if(include_uTMDPDF) then
        write(*,"('--- uTMDPDF  : ',I3,' parameters')") NPlength_uTMDPDF
        do i=1,NPlength_uTMDPDF
            write(*,"(F10.4,' ')",advance='no') lambdaNP_uTMDPDF(i)
        end do
        write(*,*) ' '
    end if
    if(include_uTMDFF) then
        write(*,"('--- uTMDFF   : ',I3,' parameters')") NPlength_uTMDFF
        do i=1,NPlength_uTMDFF
            write(*,"(F10.4,' ')",advance='no') lambdaNP_uTMDFF(i)
        end do
        write(*,*) ' '
    end if
    if(include_lpTMDPDF) then
        write(*,"('--- lpTMDPDF : ',I3,' parameters')") NPlength_lpTMDPDF
        do i=1,NPlength_lpTMDPDF
            write(*,"(F10.4,' ')",advance='no') lambdaNP_lpTMDPDF(i)
        end do
        write(*,*) ' '
    end if
    if(include_SiversTMDPDF) then
        write(*,"('--- SiversTMDPDF : ',I3,' parameters')") NPlength_SiversTMDPDF
        do i=1,NPlength_SiversTMDPDF
            write(*,"(F10.4,' ')",advance='no') lambdaNP_SiversTMDPDF(i)
        end do
        write(*,*) ' '
    end if
    if(include_wgtTMDPDF) then
        write(*,"('--- wgtTMDPDF : ',I3,' parameters')") NPlength_wgtTMDPDF
        do i=1,NPlength_wgtTMDPDF
            write(*,"(F10.4,' ')",advance='no') lambdaNP_wgtTMDPDF(i)
        end do
        write(*,*) ' '
    end if
    
    if(include_TMDF) call TMDF_ShowStatistic()
    if(include_TMDX_DY) call TMDX_DY_ShowStatistic()
    if(include_TMDX_SIDIS) call TMDX_SIDIS_ShowStatistic()


    write(*,*) '-------------------------------------------------------------------------'
end subroutine artemide_ShowStatistics

!!! changes the scale variations
subroutine artemide_SetScaleVariations(c1,c2,c3,c4)
    real(dp),intent(in)::c1,c2,c3,c4
    if(outputLevel>2) write(*,"('artemide.control: set scale variations (c1,c2,c3,c4)= (',F5.2,F5.2,F5.2,F5.2,')')") c1,c2,c3,c4

    if(include_uTMDPDF) call uTMDPDF_SetScaleVariation(c4)
    if(include_uTMDFF) call uTMDFF_SetScaleVariation(c4)
    if(include_lpTMDPDF) call lpTMDPDF_SetScaleVariation(c4)
    if(include_SiversTMDPDF) call SiversTMDPDF_SetScaleVariation(c4)
    if(include_wgtTMDPDF) call wgtTMDPDF_SetScaleVariation(c4)
    if(include_TMDs) call TMDs_SetScaleVariations(c1,c3)
    if(include_TMDX_DY) call TMDX_DY_SetScaleVariation(c2)
    if(include_TMDX_SIDIS) call TMDX_SIDIS_SetScaleVariation(c2)

end subroutine artemide_SetScaleVariations

!!!! generates some base(typical) array of NP parameters (2,0.02,0...,1,0...,etc)
function BaseNPString()
    real(dp)::BaseNPString(1:NPlength_total)
    integer::i,j
    j=1
    if(include_TMDR) then
        BaseNPString(j)=2d0
        j=j+1
        BaseNPString(j)=0.02d0
        j=j+1
        do i=3,NPlength_TMDR
            BaseNPString(j)=0d0
            j=j+1
        end do
    end if
    if(include_uTMDPDF) then
        BaseNPString(j)=1d0
        j=j+1
        do i=1,NPlength_uTMDPDF
            BaseNPString(j)=0d0
            j=j+1
        end do
    end if
    if(include_uTMDFF) then
        BaseNPString(j)=1d0
        j=j+1
        do i=1,NPlength_uTMDFF
            BaseNPString(j)=0d0
            j=j+1
        end do
    end if
    if(include_lpTMDPDF) then
        BaseNPString(j)=1d0
        j=j+1
        do i=1,NPlength_lpTMDPDF
            BaseNPString(j)=0d0
            j=j+1
        end do
    end if
    if(include_SiversTMDPDF) then
        BaseNPString(j)=1d0
        j=j+1
        do i=1,NPlength_SiversTMDPDF
            BaseNPString(j)=0d0
            j=j+1
        end do
    end if
    if(include_wgtTMDPDF) then
        BaseNPString(j)=1d0
        j=j+1
        do i=1,NPlength_wgtTMDPDF
            BaseNPString(j)=0d0
            j=j+1
        end do
    end if
    if(j-1/=NPlength_total) then
        write(*,*) ERRORstring("something went wrong on the creation of base NP parameters.",moduleName)
        stop
    end if
end function BaseNPString
  
  
!------------------------------------------------------- Working with .rep files -------------------------
!!!! this search for the string of NP parameters in the file "file" and returns rep=n
subroutine artemide_GetReplicaFromFile(file,rep,repString)
    character(len=*)::file
    integer::rep,i,k1,k2,numR,lenArray,ver
    real(dp),allocatable,intent(out)::repString(:)
    real(dp),allocatable::ParametersTOread(:)

    logical::file_exists

    INQUIRE(FILE=trim(FILE), EXIST=file_exists)
    if(.not.file_exists) then
        write(*,*) ERRORstring('replica file is not found: '//trim(FILE),moduleName)
        stop
    end if

    if(.not.isStarted .and. outputLevel>1) &
    write(*,*) WarningString(' is not initialized. Continue to read replicas without checks.',moduleName)

    OPEN(UNIT=51, FILE=trim(FILE), ACTION="read", STATUS="old")
    !!! the version
    call MoveTO(51,'*V   ')
    read(51,*) ver
    
    call MoveTO(51,'*B   ')
    !!!! check total length
    call MoveTO(51,'*0   ')
    read(51,*) lenArray
    if(isStarted) then
        if(lenArray-1/=NPlength_total) then
            write(*,*) ERRORstring('number of NP parameters in replica-file does not match the model',moduleName)
            stop
        end if
    !!! --------------------check lengths by modules

        if(include_TMDR) then
            call MoveTO(51,'*3   ')
            read(51,*) k1,k2
            if(k2-k1+1/=NPlength_TMDR) then
                write(*,*) ERRORstring('number of NP parameters in replica-file does not match the TMDR model',moduleName)
                stop
            end if
        end if

        if(include_uTMDPDF) then                    
            call MoveTO(51,'*4   ')
            read(51,*) k1,k2
            if(k2-k1+1/=NPlength_uTMDPDF) then
                write(*,*) ERRORstring('number of NP parameters in replica-file does not match the uTMDPDF model',moduleName)
                stop
            end if
        end if

        if(include_uTMDFF) then
            call MoveTO(51,'*5   ')
            read(51,*) k1,k2
            if(k2-k1+1/=NPlength_uTMDFF) then
                write(*,*) ERRORstring('number of NP parameters in replica-file does not match the uTMDFF model',moduleName)
                stop
            end if
        end if 
        
        if(ver>=3 .and. include_lpTMDPDF) then
            call MoveTO(51,'*11   ')
            read(51,*) k1,k2
            if(k2-k1+1/=NPlength_lpTMDPDF) then
                write(*,*) ERRORstring('number of NP parameters in replica-file does not match the lpTMDPDF model',moduleName)
                stop
            end if
        end if
        
        if(ver>=15 .and. include_SiversTMDPDF) then
            call MoveTO(51,'*12   ')
            read(51,*) k1,k2
            if(k2-k1+1/=NPlength_SiversTMDPDF) then
                write(*,*) ERRORstring('number of NP parameters in replica-file does not match the SiversTMDPDF model',moduleName)
                stop
            end if
        end if
        
        if(ver>=19 .and. include_wgtTMDPDF) then
            call MoveTO(51,'*13   ')
            read(51,*) k1,k2
            if(k2-k1+1/=NPlength_wgtTMDPDF) then
                write(*,*) ERRORstring('number of NP parameters in replica-file does not match the wgtTMDPDF model',moduleName)
                stop
            end if
        end if
    end if

    allocate(ParametersTOread(1:lenArray))
    if(allocated(repString)) deallocate(repString)
    allocate(repString(1:lenArray-1))


    !!---------------------Continue to reading the file
    !! indices related to the module
    call MoveTO(51,'*C   ')
    read(51,*) numR!!! full number of replicas
    !! read the replicas
    if(rep<1) then			! case of technical replicas
        call MoveTO(51,'*D   ')
        read(51,*) ParametersTOread
        if(rep==-1) then				! replica -1
            if(Int(ParametersTOread(1))==-1) then
                repString=ParametersTOread(2:lenArray)
            else
                write(*,*) ERRORstring('error in the reading the replica (-1)',moduleName)
                repString=BaseNPString()
            end if
        else if(rep==0) then			! replica 0
            read(51,*) ParametersTOread
            if(Int(ParametersTOread(1))==0) then
                repString=ParametersTOread(2:lenArray)
            else
                write(*,*) ERRORstring('error in the reading the replica (0)',moduleName)
                repString=BaseNPString()
            end if
        else					! ! replica unknown
            write(*,*) ERRORstring('there is no such replica',moduleName),rep
            repString=BaseNPString()
        end if

    else if(rep<=numR) then				! main pull of replicas
    call MoveTO(51,'*R   ')
        do i=1,numR
            read(51,*) ParametersTOread
            repString=ParametersTOread(2:lenArray)
            if(Int(ParametersTOread(1))==rep) then
                repString=ParametersTOread(2:lenArray)
                goto 101
            end if
        end do
        write(*,*) ERRORstring('there is no replica',moduleName),rep
        repString=BaseNPString()
    else
        write(*,*) ERRORstring('there is no such replica',moduleName),rep
    end if

    101 CLOSE (51, STATUS='KEEP') 
    deallocate(ParametersTOread)

end subroutine artemide_GetReplicaFromFile

!!!! this return the number of replicas stored in the file
function artemide_NumOfReplicasInFile(file)
    character(len=*)::file
    integer::numR,artemide_NumOfReplicasInFile

    logical::file_exists

    INQUIRE(FILE=trim(FILE), EXIST=file_exists)
    if(.not.file_exists) then
        write(*,*) ERRORstring('replica file is not found: '//trim(FILE),moduleName)
        stop
    end if
    OPEN(UNIT=51, FILE=trim(FILE), ACTION="read", STATUS="old")
    !! indices related to the module
    call MoveTO(51,'*C   ')
    read(51,*) numR!!! full number of replicas
    artemide_NumOfReplicasInFile=numR
    CLOSE (51, STATUS='KEEP') 

end function artemide_NumOfReplicasInFile
end module aTMDe_control
