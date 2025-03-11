!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 2.01
!
! Interface module to the user defined alpha-s, PDF, FF, etc.
! Could be interfaced to LHAPDF
!
!    ver.2.0: 28.03.2019 AV
!    ver.2.01: 14.06.2019 AV (+lpPDF)
!    ver.2.07: 09.11.2021 AV (+g1T)
!    ver.3.01: 22.08.2024 AV (introduction of LHA)
!                A.Vladimirov
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!! insertion of LHA-modules for reading PDF-grids

!!!------------------------------uPDF-------------------
module uLHAPDF_1
use IO_functions
implicit none
character (len=6),parameter :: moduleName="uPDF_1"
#include "Code/LHA/LHA_PDF.f90"
!INCLUDE 'Code/LHA/LHA_PDF.f90'
end module uLHAPDF_1

module uLHAPDF_2
use IO_functions
implicit none

character (len=6),parameter :: moduleName="uPDF_2"
#include "Code/LHA/LHA_PDF.f90"
!INCLUDE 'Code/LHA/LHA_PDF.f90'
end module uLHAPDF_2

module uLHAPDF_3
use IO_functions
implicit none
character (len=6),parameter :: moduleName="uPDF_3"
#include "Code/LHA/LHA_PDF.f90"
!INCLUDE 'Code/LHA/LHA_PDF.f90'
end module uLHAPDF_3

!!!------------------------------uFF-------------------
module uLHAFF_1
use IO_functions
implicit none
character (len=5),parameter :: moduleName="uFF_1"
#include "Code/LHA/LHA_PDF.f90"
!INCLUDE 'Code/LHA/LHA_PDF.f90'
end module uLHAFF_1

module uLHAFF_2
use IO_functions
implicit none
character (len=5),parameter :: moduleName="uFF_2"
#include "Code/LHA/LHA_PDF.f90"
!INCLUDE 'Code/LHA/LHA_PDF.f90'
end module uLHAFF_2

module uLHAFF_3
use IO_functions
implicit none
character (len=5),parameter :: moduleName="uFF_3"
#include "Code/LHA/LHA_PDF.f90"
!INCLUDE 'Code/LHA/LHA_PDF.f90'
end module uLHAFF_3

module uLHAFF_4
use IO_functions
implicit none
character (len=5),parameter :: moduleName="uFF_4"
#include "Code/LHA/LHA_PDF.f90"
!INCLUDE 'Code/LHA/LHA_PDF.f90'
end module uLHAFF_4


module uLHAFF_5
use IO_functions
implicit none
character (len=5),parameter :: moduleName="uFF_5"
#include "Code/LHA/LHA_PDF.f90"
!INCLUDE 'Code/LHA/LHA_PDF.f90'
end module uLHAFF_5

module uLHAFF_6
use IO_functions
implicit none
character (len=5),parameter :: moduleName="uFF_6"
#include "Code/LHA/LHA_PDF.f90"
!INCLUDE 'Code/LHA/LHA_PDF.f90'
end module uLHAFF_6

!!!------------------------------lpPDF-------------------
module lpLHAPDF_1
use IO_functions
implicit none
character (len=7),parameter :: moduleName="lpPDF_1"
#include "Code/LHA/LHA_PDF.f90"
!INCLUDE 'Code/LHA/LHA_PDF.f90'
end module lpLHAPDF_1

!!!------------------------------gPDF (helicity)-------------------
module gLHAPDF_1
use IO_functions
implicit none
character (len=6),parameter :: moduleName="gPDF_1"
#include "Code/LHA/LHA_PDF.f90"
!INCLUDE 'Code/LHA/LHA_PDF.f90'
end module gLHAPDF_1

module gLHAPDF_2
use IO_functions
implicit none
character (len=6),parameter :: moduleName="gPDF_2"
#include "Code/LHA/LHA_PDF.f90"
!INCLUDE 'Code/LHA/LHA_PDF.f90'
end module gLHAPDF_2

!!!------------------------------hPDF (transversity)-------------------
module hLHAPDF_1
use IO_functions
implicit none
character (len=6),parameter :: moduleName="hPDF_1"
#include "Code/LHA/LHA_PDF.f90"
!INCLUDE 'Code/LHA/LHA_PDF.f90'
end module hLHAPDF_1

module hLHAPDF_2
use IO_functions
implicit none
character (len=6),parameter :: moduleName="hPDF_2"
#include "Code/LHA/LHA_PDF.f90"
!INCLUDE 'Code/LHA/LHA_PDF.f90'
end module hLHAPDF_2

!!!----------------------------------------------------------------------
!!!----------------------------------------------------------------------
!!!----------------------------------------------------------------------

module QCDinput
use aTMDe_Numerics
use IO_functions
use LHA_alpha, only : ReadInfo_alpha => ReadInfo, AlphaS_fromLHA => AlphaS
!!
use uLHAPDF_1, only : ReadInfo_uPDF1 => ReadInfo, SetReplica_uPDF1 => SetReplica, xPDF_uPDF1 => xPDF
use uLHAPDF_2, only : ReadInfo_uPDF2 => ReadInfo, SetReplica_uPDF2 => SetReplica, xPDF_uPDF2 => xPDF
use uLHAPDF_3, only : ReadInfo_uPDF3 => ReadInfo, SetReplica_uPDF3 => SetReplica, xPDF_uPDF3 => xPDF
!!
use uLHAFF_1, only : ReadInfo_uFF1 => ReadInfo, SetReplica_uFF1 => SetReplica, xPDF_uFF1 => xPDF
use uLHAFF_2, only : ReadInfo_uFF2 => ReadInfo, SetReplica_uFF2 => SetReplica, xPDF_uFF2 => xPDF
use uLHAFF_3, only : ReadInfo_uFF3 => ReadInfo, SetReplica_uFF3 => SetReplica, xPDF_uFF3 => xPDF
use uLHAFF_4, only : ReadInfo_uFF4 => ReadInfo, SetReplica_uFF4 => SetReplica, xPDF_uFF4 => xPDF
use uLHAFF_5, only : ReadInfo_uFF5 => ReadInfo, SetReplica_uFF5 => SetReplica, xPDF_uFF5 => xPDF
use uLHAFF_6, only : ReadInfo_uFF6 => ReadInfo, SetReplica_uFF6 => SetReplica, xPDF_uFF6 => xPDF
!!
use lpLHAPDF_1, only : ReadInfo_lpPDF1 => ReadInfo, SetReplica_lpPDF1 => SetReplica, xPDF_lpPDF1 => xPDF
!!
use gLHAPDF_1, only : ReadInfo_gPDF1 => ReadInfo, SetReplica_gPDF1 => SetReplica, xPDF_gPDF1 => xPDF
use gLHAPDF_2, only : ReadInfo_gPDF2 => ReadInfo, SetReplica_gPDF2 => SetReplica, xPDF_gPDF2 => xPDF
!!
use hLHAPDF_1, only : ReadInfo_hPDF1 => ReadInfo, SetReplica_hPDF1 => SetReplica, xPDF_hPDF1 => xPDF
use hLHAPDF_2, only : ReadInfo_hPDF2 => ReadInfo, SetReplica_hPDF2 => SetReplica, xPDF_hPDF2 => xPDF
implicit none

private



public::QCDinput_Initialize,QCDinput_IsInitialized
public::As,activeNf
public::xPDF,xFF,x_lp_PDF,x_gPDF,x_hPDF
public:: QCDinput_SetPDFreplica, QCDinput_SetFFreplica, QCDinput_SetlpPDFreplica,QCDinput_SetgPDFreplica,QCDinput_SethPDFreplica

character (len=8),parameter :: moduleName="QCDinput"
!Current version of module
character (len=5),parameter :: version="v3.01"
!Last appropriate verion of constants-file
integer,parameter::inputver=35
!--- general
logical:: started=.false.
integer::outputLevel

real(dp),public::mCHARM,mBOTTOM,mTOP

!---uPDFs
integer::num_of_uPDFs
integer,allocatable::current_replica_uPDFs(:)

!---uFFs
integer::num_of_uFFs
integer,allocatable::current_replica_uFFs(:)

!---lpPDFs
integer::num_of_lpPDFs
integer,allocatable::current_replica_lpPDFs(:)

!---gPDFs
integer::num_of_gPDFs
integer,allocatable::current_replica_gPDFs(:)

!---hPDFs
integer::num_of_hPDFs
integer,allocatable::current_replica_hPDFs(:)

contains
 
 
function QCDinput_IsInitialized()
logical::QCDinput_IsInitialized
QCDinput_IsInitialized=started
end function QCDinput_IsInitialized


 !------------------------------Functions below are to be changed by user (if needed)
 subroutine QCDinput_Initialize(file,prefix)
  character(len=*),intent(in)::file
  character(len=*),optional,intent(in)::prefix
  character(len=300)::path
  character(len=600)::lineToRead

  character(len=64),allocatable::names(:)

  character(len=:),allocatable::pathToLHA
  character(len=:),allocatable::alphaNAME
  character(len=64),allocatable::names_uPDF(:),names_uFF(:),names_lpPDF(:),names_gPDF(:),names_hPDF(:)
  integer::i,FILEver
  integer,allocatable::replicas(:)

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
      CLOSE (51, STATUS='KEEP')
      write(*,*) 'artemide.'//trim(moduleName)//': const-file version is too old.'
      write(*,*) '             Update the const-file with artemide.setup'
      ERROR STOP '  '
    end if
    call MoveTO(51,'*p2  ')
    read(51,*) outputLevel    
    if(outputLevel>1) write(*,*) '--------------------------------------------- '
    if(outputLevel>1) write(*,*) 'artemide.QCDinput: initialization started ... '
    
    !!! Search for QCDinput initialization options
    call MoveTO(51,'*1   ')

    !!!! reading the path to LHAPDF
    call MoveTO(51,'*p1  ')
    read(51,"(A)") lineToRead
    pathToLHA=trim(adjustl(lineToRead))//"/"
    call MoveTO(51,'*p2  ')
    read(51,"(A)") lineToRead
    alphaNAME=trim(adjustl(lineToRead))

    !!! Search for parameters
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) mCHARM
    call MoveTO(51,'*p2  ')
    read(51,*) mBOTTOM
    call MoveTO(51,'*p3  ')
    read(51,*) mTOP
    
    if(outputLevel>2) write(*,*) '    mass of charm : ',mCHARM
    if(outputLevel>2) write(*,*) '    mass of bottom: ',mBOTTOM
    if(outputLevel>2) write(*,*) '    mass of top   : ',mTOP
    
    
    !!!---------------------------------------- Search for uPDF initialization options
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) num_of_uPDFs

    if(num_of_uPDFs>3) then
      CLOSE (51, STATUS='KEEP')
      ERROR STOP ErrorString('Maximum allowed number of uPDFs is 3. Requested '//trim(inttostr(num_of_uPDFs)),moduleName)
    else if(num_of_uPDFs>0) then

      allocate(names_uPDF(1:num_of_uPDFs))
      call MoveTO(51,'*p2  ')
      do i=1,num_of_uPDFs
        read(51,"(A)") lineToRead
        names_uPDF(i)=trim(adjustl(lineToRead))
      end do

    else
      !!! initialization is not needed
      if(outputLevel>2) write(*,*)'    no uPDFs to initialize...'
    end if
    
    !!! -------------------------------------Search for uFF initialization options
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) num_of_uFFs

    if(num_of_uFFs>6) then
      CLOSE (51, STATUS='KEEP')
      ERROR STOP ErrorString('Maximum allowed number of uFFs is 6. Requested '//trim(inttostr(num_of_uFFs)),moduleName)
    else if(num_of_uFFs>0) then

      allocate(names_uFF(1:num_of_uFFs))
      call MoveTO(51,'*p2  ')
      do i=1,num_of_uFFs
        read(51,"(A)") lineToRead
        names_uFF(i)=trim(adjustl(lineToRead))
      end do
      
    else
      !!! initialization is not needed
      if(outputLevel>2) write(*,*)'    no uFFs to initialize...'
    end if
    
     !!!---------------------------------------- Search for lpPDF initialization options
    call MoveTO(51,'*D   ')
    call MoveTO(51,'*p1  ')
    read(51,*) num_of_lpPDFs
    if(num_of_lpPDFs>1) then
      CLOSE (51, STATUS='KEEP')
      ERROR STOP ErrorString('Maximum allowed number of lpPDFs is 1. Requested '//trim(inttostr(num_of_uPDFs)),moduleName)
    else if(num_of_lpPDFs>0) then

      allocate(names_lpPDF(1:num_of_lpPDFs))
      call MoveTO(51,'*p2  ')
      do i=1,num_of_lpPDFs
        read(51,"(A)") lineToRead
        names_lpPDF(i)=trim(adjustl(lineToRead))
      end do

    else
      !!! initialization is not needed
      if(outputLevel>2) write(*,*)'    no lpPDFs to initialize...'
    end if

    
    !!!---------------------------------------- Search for gPDF initialization options
    call MoveTO(51,'*E   ')
    call MoveTO(51,'*p1  ')
    read(51,*) num_of_gPDFs
    if(num_of_gPDFs>2) then
      CLOSE (51, STATUS='KEEP')
      ERROR STOP ErrorString('Maximum allowed number of gPDFs is 2. Requested '//trim(inttostr(num_of_uPDFs)),moduleName)
    else if(num_of_gPDFs>0) then
      allocate(names_gPDF(1:num_of_gPDFs))
      call MoveTO(51,'*p2  ')
      do i=1,num_of_gPDFs
        read(51,"(A)") lineToRead
        names_gPDF(i)=trim(adjustl(lineToRead))
      end do

    else
      !!! initialization is not needed
      if(outputLevel>2) write(*,*)'    no gPDFs to initialize...'
    end if

    !!!---------------------------------------- Search for hPDF initialization options
    call MoveTO(51,'*F   ')
    call MoveTO(51,'*p1  ')
    read(51,*) num_of_hPDFs
    if(num_of_hPDFs>2) then
      CLOSE (51, STATUS='KEEP')
      ERROR STOP ErrorString('Maximum allowed number of hPDFs is 2. Requested '//trim(inttostr(num_of_uPDFs)),moduleName)
    else if(num_of_hPDFs>0) then
      allocate(names_hPDF(1:num_of_hPDFs))
      call MoveTO(51,'*p2  ')
      do i=1,num_of_hPDFs
        read(51,"(A)") lineToRead
        names_hPDF(i)=trim(adjustl(lineToRead))
      end do

    else
      !!! initialization is not needed
      if(outputLevel>2) write(*,*)'    no hPDFs to initialize...'
    end if
 
  CLOSE (51, STATUS='KEEP') 

!   !-----------------------------------------------------
!   !=====LHAPDF======
!   call InitPDFsetByNameM(1,alphaNAME)
!   call InitPDFM(1,0)
!   !=================
!   !-----------------------------------------------------

  !!!! initialize the alphaS-LHA table
  call ReadInfo_alpha(alphaNAME,pathToLHA,outputLevel)

  !!!! initialization of uPDFs (only 5 maximum PDFs allowed)
  allocate(current_replica_uPDFs(1:num_of_uPDFs))
  if(num_of_uPDFs>0) then
    current_replica_uPDFs(1)=0
    call ReadInfo_uPDF1(trim(names_uPDF(1)),pathToLHA,outputLevel)
  end if
  if(num_of_uPDFs>1) then
    current_replica_uPDFs(2)=0
    call ReadInfo_uPDF2(trim(names_uPDF(2)),pathToLHA,outputLevel)
  end if
  if(num_of_uPDFs>2) then
    current_replica_uPDFs(3)=0
    call ReadInfo_uPDF3(trim(names_uPDF(3)),pathToLHA,outputLevel)
  end if

  !!!! initialization of uFFs (only 5 maximum PDFs allowed)
  allocate(current_replica_uFFs(1:num_of_uFFs))
  if(num_of_uFFs>0) then
    current_replica_uFFs(1)=0
    call ReadInfo_uFF1(trim(names_uFF(1)),pathToLHA,outputLevel)
  end if
  if(num_of_uFFs>1) then
    current_replica_uFFs(2)=0
    call ReadInfo_uFF2(trim(names_uFF(2)),pathToLHA,outputLevel)
  end if
  if(num_of_uFFs>2) then
    current_replica_uFFs(3)=0
    call ReadInfo_uFF3(trim(names_uFF(3)),pathToLHA,outputLevel)
  end if
  if(num_of_uFFs>3) then
    current_replica_uFFs(4)=0
    call ReadInfo_uFF4(trim(names_uFF(4)),pathToLHA,outputLevel)
  end if
  if(num_of_uFFs>4) then
    current_replica_uFFs(5)=0
    call ReadInfo_uFF5(trim(names_uFF(5)),pathToLHA,outputLevel)
  end if
  if(num_of_uFFs>5) then
    current_replica_uFFs(6)=0
    call ReadInfo_uFF6(trim(names_uFF(6)),pathToLHA,outputLevel)
  end if

  !!!! initialization of lpPDFs (only 1 maximum PDFs allowed)
  allocate(current_replica_lpPDFs(1:num_of_lpPDFs))
  if(num_of_lpPDFs>0) then
    current_replica_lpPDFs(1)=0
    call ReadInfo_lpPDF1(trim(names_lpPDF(1)),pathToLHA,outputLevel)
  end if

  !!!! initialization of gPDFs (only 2 maximum PDFs allowed)
  allocate(current_replica_gPDFs(1:num_of_gPDFs))
  if(num_of_gPDFs>0) then
    current_replica_gPDFs(1)=0
    call ReadInfo_gPDF1(trim(names_gPDF(1)),pathToLHA,outputLevel)
  end if
  if(num_of_gPDFs>1) then
    current_replica_gPDFs(2)=0
    call ReadInfo_gPDF2(trim(names_gPDF(2)),pathToLHA,outputLevel)
  end if

  !!!! initialization of hPDFs (only 2 maximum PDFs allowed)
  allocate(current_replica_hPDFs(1:num_of_hPDFs))
  if(num_of_hPDFs>0) then
    current_replica_hPDFs(1)=0
    call ReadInfo_hPDF1(trim(names_hPDF(1)),pathToLHA,outputLevel)
  end if
  if(num_of_hPDFs>1) then
    current_replica_hPDFs(2)=0
    call ReadInfo_hPDF2(trim(names_hPDF(2)),pathToLHA,outputLevel)
  end if
     
   if(outputLevel>0) write(*,*) color('----- arTeMiDe.QCDinput '//trim(version)//': .... initialized',c_green)
   if(outputLevel>1) write(*,*) ' '
  started=.true.
 end subroutine QCDinput_Initialize

!!! set a different replica of uPDF (pointing to hadron)
!!! check whatever its the same or not.
!!! in the case the replica is same as stored -- does not change PDF (newPDF=false)
!!! in the case the replica is different -- does change PDF (newPDF=true)
subroutine QCDinput_SetPDFreplica(rep,hadron,newPDF)
    integer,intent(in):: rep,hadron
    logical,intent(out)::newPDF

    if(hadron<1 .or. hadron>num_of_uPDFs) &
      ERROR STOP ErrorString('SetPDFreplica. Called unexisting hadron. h= '//trim(inttostr(hadron)),moduleName)
    !!! if the number of replica to change coincides with the already used. Do not change it
    if(current_replica_uPDFs(hadron)==rep) then
      newPDF=.false.
    else
      current_replica_uPDFs(hadron)=rep

      SELECT CASE(hadron)
        CASE (1)
          call SetReplica_uPDF1(rep)
        CASE (2)
          call SetReplica_uPDF2(rep)
        CASE (3)
          call SetReplica_uPDF3(rep)
        END SELECT

      newPDF=.true.
    end if
end subroutine QCDinput_SetPDFreplica

!!! set a different replica of uFF (pointing to hadron)
!!! check whatever its the same or not.
!!! in the case the replica is same as stored -- does not change PDF (newPDF=false)
!!! in the case the replica is different -- does change PDF (newPDF=true)
subroutine QCDinput_SetFFreplica(rep,hadron,newPDF)
    integer,intent(in):: rep,hadron
    logical,intent(out)::newPDF

    if(hadron<1 .or. hadron>num_of_uFFs) &
      ERROR STOP ErrorString('SetFFreplica. Called unexisting hadron. h= '//trim(inttostr(hadron)),moduleName)
    !!! if the number of replica to change coincides with the already used. Do not change it
    if(current_replica_uFFs(hadron)==rep) then
      newPDF=.false.
    else
      current_replica_uFFs(hadron)=rep

      SELECT CASE(hadron)
        CASE (1)
          call SetReplica_uFF1(rep)
        CASE (2)
          call SetReplica_uFF2(rep)
        CASE (3)
          call SetReplica_uFF3(rep)
        CASE (4)
          call SetReplica_uFF4(rep)
        CASE (5)
          call SetReplica_uFF5(rep)
        CASE (6)
          call SetReplica_uFF6(rep)
        END SELECT

      newPDF=.true.
    end if
end subroutine QCDinput_SetFFreplica

!!! set a different replica of lp_PDF (pointing to hadron)
!!! check whatever its the same or not.
!!! in the case the replica is same as stored -- does not change PDF (newPDF=false)
!!! in the case the replica is different -- does change PDF (newPDF=true)
subroutine QCDinput_SetlpPDFreplica(rep,hadron,newPDF)
    integer,intent(in):: rep,hadron
    logical,intent(out)::newPDF

    if(hadron<1 .or. hadron>num_of_lpPDFs) &
      ERROR STOP ErrorString('SetlpPDFreplica. Called unexisting hadron. h= '//trim(inttostr(hadron)),moduleName)
    !!! if the number of replica to change coincides with the already used. Do not change it
    if(current_replica_lpPDFs(hadron)==rep) then
      newPDF=.false.
    else
      current_replica_lpPDFs(hadron)=rep

      SELECT CASE(hadron)
        CASE (1)
          call SetReplica_lpPDF1(rep)
        END SELECT

      newPDF=.true.
    end if
end subroutine QCDinput_SetlpPDFreplica

!!! set a different replica of gPDF (pointing to hadron)
!!! check whatever its the same or not.
!!! in the case the replica is same as stored -- does not change PDF (newPDF=false)
!!! in the case the replica is different -- does change PDF (newPDF=true)
subroutine QCDinput_SetgPDFreplica(rep,hadron,newPDF)
    integer,intent(in):: rep,hadron
    logical,intent(out)::newPDF

    if(hadron<1 .or. hadron>num_of_gPDFs) &
      ERROR STOP ErrorString('SetgPDFreplica. Called unexisting hadron. h= '//trim(inttostr(hadron)),moduleName)
    !!! if the number of replica to change coincides with the already used. Do not change it
    if(current_replica_gPDFs(hadron)==rep) then
      newPDF=.false.
    else
      current_replica_gPDFs(hadron)=rep

      SELECT CASE(hadron)
        CASE (1)
          call SetReplica_gPDF1(rep)
        CASE (2)
          call SetReplica_gPDF2(rep)
        END SELECT

      newPDF=.true.
    end if
end subroutine QCDinput_SetgPDFreplica

!!! set a different replica of hPDF (pointing to hadron)
!!! check whatever its the same or not.
!!! in the case the replica is same as stored -- does not change PDF (newPDF=false)
!!! in the case the replica is different -- does change PDF (newPDF=true)
subroutine QCDinput_SethPDFreplica(rep,hadron,newPDF)
    integer,intent(in):: rep,hadron
    logical,intent(out)::newPDF

    if(hadron<1 .or. hadron>num_of_hPDFs) &
      ERROR STOP ErrorString('SethPDFreplica. Called unexisting hadron. h= '//trim(inttostr(hadron)),moduleName)
    !!! if the number of replica to change coincides with the already used. Do not change it
    if(current_replica_hPDFs(hadron)==rep) then
      newPDF=.false.
    else
      current_replica_hPDFs(hadron)=rep

      SELECT CASE(hadron)
        CASE (1)
          call SetReplica_hPDF1(rep)
        CASE (2)
          call SetReplica_hPDF2(rep)
        END SELECT

      newPDF=.true.
    end if
end subroutine QCDinput_SethPDFreplica
 
!!number of active flavor at given mu
function activeNf(mu)
  real(dp),intent(in)::mu
  integer::activeNf

  if(mu>mBOTTOM) then
     activeNf=5
    else if(mu>mCHARM) then
     activeNf=4
    else
     activeNf=3
  end if
end function activeNf
 
!!!!alphas(Q)/4pi
!!! NOT FORGET 4 PI !!!
function As(Q)
real(dp),intent(in)::Q
real(dp)::As

As=AlphaS_fromLHA(Q)/pix4
!  as=1d0/(2d0*23d0/3d0*Log(Q/0.08782708014552364d0))  !!! Nf=5 LO solution


! !-------------------------
! !======LHAPDF========
! real(dp)::alphasPDF
! As=alphasPDF(Q)/pix4
! !======LHAPDF========
! !-------------------------


end function As

!!!!array of x times PDF(x,Q) for hadron 'hadron'
!!! unpolarized PDF used in uTMDPDF
!!!! array is (-5:5) (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
function xPDF(x,Q,hadron)
    real(dp),intent(in) :: x,Q
    integer,intent(in):: hadron
    real(dp), dimension(-5:5):: xPDF

!     !-------------------------
!     !======LHAPDF========
!     real(dp), dimension(-6:6):: inputPDF
!     call evolvePDFM(1,x,Q,inputPDF)
!     xPDF=inputPDF(-5:5)
!     !======LHAPDF========
!     !-------------------------

    SELECT CASE(hadron)
      CASE(1)
        xPDF=xPDF_uPDF1(x,Q)
      CASE(2)
        xPDF=xPDF_uPDF2(x,Q)
      CASE(3)
        xPDF=xPDF_uPDF3(x,Q)
      CASE DEFAULT
        ERROR STOP ErrorString('xPDF. Called unexisting hadron. h= '//trim(inttostr(hadron)),moduleName)
    END SELECT

end function xPDF
  
!!!! return x*F(x,mu)
!!!! enumeration of flavors
!!!!  f = -5,-4, -3,  -2,  -1,0,1,2,3,5,4
!!!!    = bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b
!!! unpolarized FF used in uTMDFF
!!!! enumeration of hadrons
function xFF(x,Q,hadron)
    integer,intent(in) :: hadron
    real(dp),intent(in) :: x,Q
    real(dp),dimension(-5:5):: xFF

    SELECT CASE(hadron)
      CASE(1)
        xFF=xPDF_uFF1(x,Q)
      CASE(2)
        xFF=xPDF_uFF2(x,Q)
      CASE(3)
        xFF=xPDF_uFF3(x,Q)
      CASE(4)
        xFF=xPDF_uFF4(x,Q)
      CASE(5)
        xFF=xPDF_uFF5(x,Q)
      CASE(6)
        xFF=xPDF_uFF6(x,Q)
      CASE DEFAULT
        ERROR STOP ErrorString('xFF. Called unexisting hadron. h= '//trim(inttostr(hadron)),moduleName)
    END SELECT

end function xFF
 
!!!!array of x times PDF(x,Q) for hadron 'hadron'
!!! unpolarized PDF used in lpTMDPDF
!!!! array is (-5:5) (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
function x_lp_PDF(x,Q,hadron)
    real(dp),intent(in) :: x,Q
    integer,intent(in):: hadron
    real(dp), dimension(-5:5):: x_lp_PDF

    SELECT CASE(hadron)
    CASE(1)
      x_lp_PDF=xPDF_lpPDF1(x,Q)
    CASE DEFAULT
      ERROR STOP ErrorString('x_lp_PDF. Called unexisting hadron. h= '//trim(inttostr(hadron)),moduleName)
  END SELECT

end function x_lp_PDF

!!!!array of x times PDF(x,Q) for hadron 'hadron'
!!! helicity PDF
!!!! array is (-5:5) (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
function x_gPDF(x,Q,hadron)
    real(dp),intent(in) :: x,Q
    integer,intent(in):: hadron
    real(dp), dimension(-5:5):: x_gPDF

    SELECT CASE(hadron)
    CASE(1)
      x_gPDF=xPDF_gPDF1(x,Q)
    CASE(2)
      x_gPDF=xPDF_gPDF2(x,Q)
    CASE DEFAULT
      ERROR STOP ErrorString('xgPDF. Called unexisting hadron. h= '//trim(inttostr(hadron)),moduleName)
  END SELECT

end function x_gPDF

!!!!array of x times PDF(x,Q) for hadron 'hadron'
!!! transversity PDF
!!!! array is (-5:5) (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
function x_hPDF(x,Q,hadron)
    real(dp),intent(in) :: x,Q
    integer,intent(in):: hadron
    real(dp), dimension(-5:5):: x_hPDF

    SELECT CASE(hadron)
    CASE(1)
      x_hPDF=xPDF_hPDF1(x,Q)
    CASE(2)
      x_hPDF=xPDF_hPDF2(x,Q)
    CASE DEFAULT
      ERROR STOP ErrorString('xhPDF. Called unexisting hadron. h= '//trim(inttostr(hadron)),moduleName)
  END SELECT

end function x_hPDF
 
end module QCDinput
