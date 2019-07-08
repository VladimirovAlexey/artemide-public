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
use QCDinput
use EWinput
use uTMDPDF
use lpTMDPDF
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
  character (len=5),parameter :: version="v2.01"
  !Last appropriate verion of constants-file
  integer,parameter::inputver=3
  character (len=15),parameter :: constNAME="aTMDe-temporary"
  
  integer::outputLevel=2
  integer::messageTrigger=5
  
  !!!! indicators of modules usage
  logical::include_EWinput,include_uTMDPDF,include_uTMDFF,include_TMDR,include_TMDs,include_TMDF,include_lpTMDPDF
  logical::include_TMDX_DY,include_TMDX_SIDIS,include_TMDs_inKT
  
  !!!! legths of non-perturbative arrays
  integer::NPlength_total
  integer::NPlength_TMDR,NPlength_uTMDPDF,NPlength_uTMDFF,NPlength_lpTMDPDF
  
  !!!! non-pertrubative parameters for individual modules
  real*8,allocatable::lambdaNP_TMDR(:),lambdaNP_uTMDPDF(:),lambdaNP_uTMDFF(:),lambdaNP_lpTMDPDF(:)
  
  !!!! Saved values of scale-variation parameters
  real*8::c1_saved,c2_saved,c3_saved,c4_saved
  
  public::artemide_Initialize
  public::artemide_SetNPparameters,artemide_SetNPparameters_TMDR,artemide_SetNPparameters_uTMDFF,artemide_SetNPparameters_uTMDPDF
  public::artemide_SetNPparameters_lpTMDPDF
  public::artemide_SetReplica_TMDR,artemide_SetReplica_uTMDFF,artemide_SetReplica_uTMDPDF,artemide_SetReplica_lpTMDPDF
  public::artemide_SetScaleVariations
  public::artemide_ShowStatistics

contains
  
   !!! move CURRET in streem to the next line that starts from pos (5 char)
 subroutine MoveTO(streem,pos)
 integer,intent(in)::streem
 character(len=5)::pos
 character(len=300)::line
    do
    read(streem,'(A)') line    
    if(line(1:5)==pos) exit
    end do
 end subroutine MoveTO
  
  !! Module that creates the constants-file,
  !! and the read it
  !! intialize it-self and artemide modules
  subroutine artemide_Initialize(file,prefix,order)
  character(len=*)::file
  character(len=*),optional::prefix
  character(len=*),optional::order
  character(len=300)::path
  integer::FILEver
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
    if(outputLevel>2) write(*,*) 'artemide.control: initialization started ... '
    call MoveTO(51,'*p3  ')
    read(51,*) outputLevel    
    
    call MoveTO(51,'*2   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_EWinput
    
    call MoveTO(51,'*3   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_TMDR
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) NPlength_TMDR
    
    call MoveTO(51,'*4   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_uTMDPDF
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) NPlength_uTMDPDF
    
    call MoveTO(51,'*5   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_uTMDFF
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) NPlength_uTMDFF
    
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
    
    CLOSE (51, STATUS='KEEP')
   
   !-----------------------------------------------------------
    if(outputLevel>2) write(*,*) 'artemide.control: initialization file read. Initialization of modules ... '
    
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
   
   if(outputLevel>2) write(*,*) ' artemide.control: Total number of NP parameters:',NPlength_total
   
   !allocate lambda's
   if(include_TMDR) allocate(lambdaNP_TMDR(1:NPlength_TMDR))
   if(include_uTMDPDF) allocate(lambdaNP_uTMDPDF(1:NPlength_uTMDPDF))
   if(include_uTMDFF) allocate(lambdaNP_uTMDFF(1:NPlength_uTMDFF))
   if(include_lpTMDPDF) allocate(lambdaNP_lpTMDPDF(1:NPlength_lpTMDPDF))
   
   c1_saved=1d0
   c2_saved=1d0
   c3_saved=1d0
   c4_saved=1d0
   
   call artemide_SetScaleVariations(c1_saved,c2_saved,c3_saved,c4_saved)
   
   if(outputLevel>1) write(*,*) ' artemide.control: initialization done.'
  
  end subroutine artemide_Initialize
  
  !-------------------------------------------------------------- NP parameters ---------------------------
  subroutine artemide_SetNPparameters(lambdaNP)
  real*8,intent(in)::lambdaNP(:)
  real*8,dimension(1:NPlength_total)::lambda_cur
  integer::ll,num
  
  ll=size(lambdaNP)
  
  if(ll<NPlength_total) then
    if(outputLevel>0) write(*,"('artemide.SetNPparameters: ERROR: the length of NP parameters array ('&
	    ,I4,') is smaller then the total number of NP parameters (',I4,')')")ll,NPlength_total
    if(outputLevel>0) write(*,*) 'NOTHING IS DONE'
    return
  end if
  
  if(ll>NPlength_total) then
    if(outputLevel>0) write(*,"('artemide.SetNPparameters: ERROR: the length of NP parameters array ('&
	    ,I4,') is larger then total the number of NP parameters (',I4,')')")ll,NPlength_total
    if(outputLevel>0) write(*,*) 'The array is trucated'
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
  
  !!! sending NP arrays to packages
  if(include_TMDR) call TMDR_setNPparameters(lambdaNP_TMDR)
  if(include_uTMDPDF) call uTMDPDF_SetLambdaNP(lambdaNP_uTMDPDF)
  if(include_uTMDFF) call uTMDFF_SetLambdaNP(lambdaNP_uTMDFF)
  if(include_lpTMDPDF) call lpTMDPDF_SetLambdaNP(lambdaNP_lpTMDPDF)
  
  !!! reseting other packages
  if(include_TMDF) call TMDF_ResetCounters()
  if(include_TMDX_DY) call TMDX_DY_ResetCounters()
  if(include_TMDX_SIDIS) call TMDX_SIDIS_ResetCounters()
  
  end subroutine artemide_SetNPparameters
  
 subroutine artemide_SetNPparameters_TMDR(lambdaNP)
  real*8,intent(in)::lambdaNP(:)
  integer::ll

  if(.not.include_TMDR) then
    if(outputLevel>0) &
      write(*,*) &
      'artemide.control: ERROR: attempt to set NP-parameters for TMDR, while TMDR module is not included in the current setup'
    if(outputLevel>0) write(*,*) 'NOTHING IS DONE'
    return
  end if
  
  ll=size(lambdaNP)
  
  
  if(ll<NPlength_TMDR) then
    if(outputLevel>0) write(*,"('artemide.SetNPparameters-TMDR: ERROR: the length of NP parameters array ('&
	    ,I4,') is smaller then the total number of NP parameters for TMDR (',I4,')')")ll,NPlength_TMDR
    if(outputLevel>0) write(*,*) 'NOTHING IS DONE'
    return
  end if
  
  if(ll>NPlength_total) then
    if(outputLevel>0) write(*,"('artemide.SetNPparameters: ERROR: the length of NP parameters array ('&
	    ,I4,') is larger then total the number of NP parameters for TMDR (',I4,')')")ll,NPlength_TMDR
    if(outputLevel>0) write(*,*) 'The array is trucated'
  end if
  
  lambdaNP_TMDR=lambdaNP(1:NPlength_TMDR)
  
  call TMDR_setNPparameters(lambdaNP_TMDR)
  
  !!! reseting other packages
  if(include_TMDF) call TMDF_ResetCounters()
  if(include_TMDX_DY) call TMDX_DY_ResetCounters()
  if(include_TMDX_SIDIS) call TMDX_SIDIS_ResetCounters()
  
  end subroutine artemide_SetNPparameters_TMDR
  
  subroutine artemide_SetNPparameters_uTMDPDF(lambdaNP)
  real*8,intent(in)::lambdaNP(:)
  integer::ll
  
  if(.not.include_uTMDPDF) then
    if(outputLevel>0) &
      write(*,*) 'artemide.control: ERROR: attempt to set NP-parameters for uTMDPDF,',&
      ' while uTMDPDF module is not included in the current setup'
    if(outputLevel>0) write(*,*) 'NOTHING IS DONE'
    return
  end if
  
  ll=size(lambdaNP)
  
  
  if(ll<NPlength_uTMDPDF) then
    if(outputLevel>0) write(*,"('artemide.SetNPparameters-uTMDPDF: ERROR: the length of NP parameters array ('&
	    ,I4,') is smaller then the total number of NP parameters for uTMDPDF (',I4,')')")ll,NPlength_uTMDPDF
    if(outputLevel>0) write(*,*) 'NOTHING IS DONE'
    return
  end if
  
  if(ll>NPlength_total) then
    if(outputLevel>0) write(*,"('artemide.SetNPparameters: ERROR: the length of NP parameters array ('&
	    ,I4,') is larger then total the number of NP parameters for uTMDPDF (',I4,')')")ll,NPlength_uTMDPDF
    if(outputLevel>0) write(*,*) 'The array is trucated'
  end if
  
  lambdaNP_uTMDPDF=lambdaNP(1:NPlength_uTMDPDF)
  
  call uTMDPDF_SetLambdaNP(lambdaNP_uTMDPDF)
  
  !!! reseting other packages
  if(include_TMDF) call TMDF_ResetCounters()
  if(include_TMDX_DY) call TMDX_DY_ResetCounters()
  if(include_TMDX_SIDIS) call TMDX_SIDIS_ResetCounters()
  
  end subroutine artemide_SetNPparameters_uTMDPDF
  
  subroutine artemide_SetNPparameters_uTMDFF(lambdaNP)
  real*8,intent(in)::lambdaNP(:)
  integer::ll
  
  if(.not.include_uTMDFF) then
    if(outputLevel>0) &
      write(*,*) 'artemide.control: ERROR: attempt to set NP-parameters for uTMDFF,',&
      ' while uTMDFF module is not included in the current setup'
    if(outputLevel>0) write(*,*) 'NOTHING IS DONE'
    return
  end if
  
  ll=size(lambdaNP)
  
  
  if(ll<NPlength_uTMDFF) then
    if(outputLevel>0) write(*,"('artemide.SetNPparameters-uTMDFF: ERROR: the length of NP parameters array ('&
	    ,I4,') is smaller then the total number of NP parameters for uTMDFF (',I4,')')")ll,NPlength_uTMDFF
    if(outputLevel>0) write(*,*) 'NOTHING IS DONE'
    return
  end if
  
  if(ll>NPlength_total) then
    if(outputLevel>0) write(*,"('artemide.SetNPparameters: ERROR: the length of NP parameters array ('&
	    ,I4,') is larger then total the number of NP parameters for uTMDFF (',I4,')')")ll,NPlength_uTMDFF
    if(outputLevel>0) write(*,*) 'The array is trucated'
  end if
  
  lambdaNP_uTMDFF=lambdaNP(1:NPlength_uTMDFF)
  
  call uTMDFF_SetLambdaNP(lambdaNP_uTMDFF)
  
  !!! reseting other packages
  if(include_TMDF) call TMDF_ResetCounters()
  if(include_TMDX_DY) call TMDX_DY_ResetCounters()
  if(include_TMDX_SIDIS) call TMDX_SIDIS_ResetCounters()
  
  end subroutine artemide_SetNPparameters_uTMDFF

   
  subroutine artemide_SetNPparameters_lpTMDPDF(lambdaNP)
  real*8,intent(in)::lambdaNP(:)
  integer::ll
  
  if(.not.include_lpTMDPDF) then
    if(outputLevel>0) &
      write(*,*) 'artemide.control: ERROR: attempt to set NP-parameters for lpTMDPDF,',&
      ' while lpTMDPDF module is not included in the current setup'
    if(outputLevel>0) write(*,*) 'NOTHING IS DONE'
    return
  end if
  
  ll=size(lambdaNP)
  
  
  if(ll<NPlength_lpTMDPDF) then
    if(outputLevel>0) write(*,"('artemide.SetNPparameters-lpTMDPDF: ERROR: the length of NP parameters array ('&
	    ,I4,') is smaller then the total number of NP parameters for lpTMDPDF (',I4,')')")ll,NPlength_lpTMDPDF
    if(outputLevel>0) write(*,*) 'NOTHING IS DONE'
    return
  end if
  
  if(ll>NPlength_total) then
    if(outputLevel>0) write(*,"('artemide.SetNPparameters: ERROR: the length of NP parameters array ('&
	    ,I4,') is larger then total the number of NP parameters for lpTMDPDF (',I4,')')")ll,NPlength_lpTMDPDF
    if(outputLevel>0) write(*,*) 'The array is trucated'
  end if
  
  lambdaNP_lpTMDPDF=lambdaNP(1:NPlength_lpTMDPDF)
  
  call lpTMDPDF_SetLambdaNP(lambdaNP_lpTMDPDF)
  
  !!! reseting other packages
  if(include_TMDF) call TMDF_ResetCounters()
  if(include_TMDX_DY) call TMDX_DY_ResetCounters()
  if(include_TMDX_SIDIS) call TMDX_SIDIS_ResetCounters()
  
  end subroutine artemide_SetNPparameters_lpTMDPDF
 
  subroutine artemide_SetReplica_TMDR(num)
  integer::num
  
  if(.not.include_TMDR) then
    if(outputLevel>0) &
	write(*,*) 'artemide.control: ERROR: attempt to set a replica for TMDR,',&
	    ' while TMDR module is not included in the current setup'
    if(outputLevel>0) write(*,*) 'NOTHING IS DONE'
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
  integer::num
  
  if(.not.include_uTMDPDF) then
    if(outputLevel>0) &
	write(*,*) 'artemide.control: ERROR: attempt to set a replica for uTMDPDF,',&
	    ' while uTMDPDF module is not included in the current setup'
    if(outputLevel>0) write(*,*) 'NOTHING IS DONE'
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
  integer::num
  
  if(.not.include_uTMDFF) then
    if(outputLevel>0) &
	write(*,*) 'artemide.control: ERROR: attempt to set a replica for uTMDFF,',&
	    ' while uTMDFF module is not included in the current setup'
    if(outputLevel>0) write(*,*) 'NOTHING IS DONE'
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
  integer::num
  
  if(.not.include_uTMDPDF) then
    if(outputLevel>0) &
	write(*,*) 'artemide.control: ERROR: attempt to set a replica for lpTMDPDF,',&
	    ' while lpTMDPDF module is not included in the current setup'
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
  if(include_TMDF) call TMDF_ShowStatistic()
  if(include_TMDX_DY) call TMDX_DY_ShowStatistic()
  if(include_TMDX_SIDIS) call TMDX_SIDIS_ShowStatistic()
  
  
  write(*,*) '-------------------------------------------------------------------------'
  end subroutine artemide_ShowStatistics

  !!! changes the scale variations
  subroutine artemide_SetScaleVariations(c1,c2,c3,c4)
  real*8::c1,c2,c3,c4
  if(outputLevel>2) write(*,"('artemide.control: set scale variations (c1,c2,c3,c4)= (',F5.2,F5.2,F5.2,F5.2,')')") c1,c2,c3,c4
  
  if(include_uTMDPDF) call uTMDPDF_SetScaleVariation(c4)
  if(include_uTMDFF) call uTMDFF_SetScaleVariation(c4)
  if(include_lpTMDPDF) call lpTMDPDF_SetScaleVariation(c4)
  if(include_TMDs) call TMDs_SetScaleVariations(c1,c3)
  if(include_TMDX_DY) call TMDX_DY_SetScaleVariation(c2)
  if(include_TMDX_SIDIS) call TMDX_SIDIS_SetScaleVariation(c2)
  
  end subroutine artemide_SetScaleVariations
end module aTMDe_control