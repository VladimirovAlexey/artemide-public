!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 2.01
!
! Interface module to the user defined alpha-s, PDF, FF, etc.
! Could be interfaced to LHAPDF
!
!	ver.2.0: 28.03.2019 AV
!	ver.2.01: 14.06.2019 AV (+lpPDF)
!				A.Vladimirov
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module QCDinput
use aTMDe_Numerics
use IO_functions
implicit none

private



public::QCDinput_Initialize,QCDinput_IsInitialized
public::As,activeNf
public::xPDF,xFF,x_lp_PDF
public:: QCDinput_SetPDFreplica

character (len=8),parameter :: moduleName="QCDinput"
!Current version of module
character (len=5),parameter :: version="v2.03"
!Last appropriate verion of constants-file
integer,parameter::inputver=4
!--- general
logical:: started=.false.
integer::outputLevel

real(dp),public::mCHARM,mBOTTOM,mTOP

!---uPDFs
integer::num_of_uPDFs,startPDFindex
integer,allocatable::enumeration_of_uPDFs(:)

!---uFFs
integer::num_of_uFFs,startFFindex
integer,allocatable::enumeration_of_uFFs(:)

!---lpPDFs
integer::num_of_lpPDFs,startlpPDFindex
integer,allocatable::enumeration_of_lpPDFs(:)

 contains 
 
 
 function QCDinput_IsInitialized()
  logical::QCDinput_IsInitialized
  QCDinput_IsInitialized=started 
 end function QCDinput_IsInitialized


 !------------------------------Functions below are to be changed by user (if needed)
 subroutine QCDinput_Initialize(file,prefix)
  character(len=*)::file
  character(len=*),optional::prefix
  character(len=300)::path
  character(len=64),allocatable::names(:)
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
      write(*,*) 'artemide.'//trim(moduleName)//': const-file version is too old.'
      write(*,*) '		     Update the const-file with artemide.setup'
      write(*,*) '  '
      stop
    end if
    call MoveTO(51,'*p2  ')
    read(51,*) outputLevel    
    if(outputLevel>1) write(*,*) '--------------------------------------------- '
    if(outputLevel>1) write(*,*) 'artemide.QCDinput: initialization started ... '
    
    !!! Search for QCDinput initialization options
    call MoveTO(51,'*1   ')
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
    if(num_of_uPDFs>0) then
      !! initialization of LHAPDF grids
      allocate(enumeration_of_uPDFs(1:num_of_uPDFs))
      allocate(names(1:num_of_uPDFs))
      allocate(replicas(1:num_of_uPDFs))
      call MoveTO(51,'*p2  ')
      read(51,*) enumeration_of_uPDFs
      call MoveTO(51,'*p3  ')
      do i=1,num_of_uPDFs
       read(51,*) names(i)
      end do
      call MoveTO(51,'*p4  ')
      read(51,*) replicas
      
      !!! actually initialization
      startPDFindex=0
      do i=1,num_of_uPDFs
        call InitPDFsetByNameM(i+startPDFindex,names(i))
        call InitPDFM(i+startPDFindex,replicas(i))
        if(outputLevel>2) write(*,"('     uPDF(hadron=',I3,') initialized by : ',A,' (replica= ',I5,')')") &
		    enumeration_of_uPDFs(i),trim(names(i)),replicas(i)
      end do
      
      deallocate(names,replicas)
    else
      !!! initialization is not needed
      if(outputLevel>2)	write(*,*)'    no uPDFs to initialize...'
    end if
    
    !!! -------------------------------------Search for uFF initialization options
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) num_of_uFFs
    if(num_of_uFFs>0) then
      allocate(enumeration_of_uFFs(1:num_of_uFFs))
      allocate(names(1:num_of_uFFs))
      allocate(replicas(1:num_of_uFFs))
      call MoveTO(51,'*p2  ')
      read(51,*) enumeration_of_uFFs
      call MoveTO(51,'*p3  ')
      do i=1,num_of_uFFs
       read(51,*) names(i)
      end do
      call MoveTO(51,'*p4  ')
      read(51,*) replicas
      
      !!! actually initialization
      startFFindex=num_of_uPDFs+startPDFindex
      do i=1,num_of_uFFs
	call InitPDFsetByNameM(startFFindex+i,names(i))
	call InitPDFM(startFFindex+i,replicas(i))
	if(outputLevel>2) write(*,"('      uFF(hadron=',I3,') initialized by : ',A,' (replica= ',I5,')')") &
		    enumeration_of_uFFs(i),trim(names(i)),replicas(i)
      end do
      
      deallocate(names,replicas)
      
    else
      !!! initialization is not needed
      if(outputLevel>2)	write(*,*)'    no uFFs to initialize...'
    end if
    
     !!!---------------------------------------- Search for lpPDF initialization options
    call MoveTO(51,'*D   ')
    call MoveTO(51,'*p1  ')
    read(51,*) num_of_lpPDFs    
    if(num_of_lpPDFs>0) then
      !! initialization of LHAPDF grids
      allocate(enumeration_of_lpPDFs(1:num_of_lpPDFs))
      allocate(names(1:num_of_lpPDFs))
      allocate(replicas(1:num_of_lpPDFs))
      call MoveTO(51,'*p2  ')
      read(51,*) enumeration_of_lpPDFs
      call MoveTO(51,'*p3  ')
      do i=1,num_of_lpPDFs
       read(51,*) names(i)
      end do
      call MoveTO(51,'*p4  ')
      read(51,*) replicas
      
      !!! actually initialization
      startlpPDFindex=0
      do i=1,num_of_lpPDFs
	call InitPDFsetByNameM(i+startlpPDFindex,names(i))
	call InitPDFM(i+startlpPDFindex,replicas(i))
	if(outputLevel>2) write(*,"('     lpPDF(hadron=',I3,') initialized by : ',A,' (replica= ',I5,')')") &
		    enumeration_of_lpPDFs(i),trim(names(i)),replicas(i)
      end do
      
      deallocate(names,replicas)
    else
      !!! initialization is not needed
      if(outputLevel>2)	write(*,*)'    no lpPDFs to initialize...'
    end if
 
  CLOSE (51, STATUS='KEEP') 
     
   if(outputLevel>0) write(*,*) color('----- arTeMiDe.QCDinput '//trim(version)//': .... initialized',c_green)
   if(outputLevel>1) write(*,*) ' '
  started=.true.
 end subroutine QCDinput_Initialize
 
 !!! provide the index of grid associated with uPDF(hadron)
 function index_of_uPDF(hadron)
  integer,intent(in)::hadron
  integer::index_of_uPDF
  integer::i
  
  if(num_of_uPDFs==0) then
    write(*,*) ErrorString('no uPDFs are initialized',moduleName)
    stop
  else
    do i=1,num_of_uPDFs
      if(enumeration_of_uPDFs(i)==hadron) then
	index_of_uPDF=i+startPDFindex
	return
      end if
    end do
    !!! if we exit from the loop it means index is not found
      write(*,*) ErrorString('uPDF is not found',moduleName)
      write(*,"('no uPDF for hadron ',I3,' is initialized. Set PDF for hadron (',I3,')')") &
		  hadron,enumeration_of_uPDFs(1)
      index_of_uPDF=1+startPDFindex
  end if
 end function index_of_uPDF
 
  !!! provide the index of grid associated with uFF(hadron)
 function index_of_uFF(hadron)
  integer,intent(in)::hadron
  integer::index_of_uFF
  integer::i
  
  if(num_of_uFFs==0) then
    write(*,*) ErrorString('no uFFs are initialized',moduleName)
    stop
   else
    do i=1,num_of_uFFs
      if(enumeration_of_uFFs(i)==hadron) then
	index_of_uFF=i+startFFindex
	return
      end if
    end do
    !!! if we exit from the loop it means index is not found
     write(*,*) ErrorString('uFF is not found',moduleName)
      write(*,"('no uFF for hadron ',I3,' is initialized. Set FF for hadron (',I3,')')") &
		  hadron,enumeration_of_uPDFs(1)
      index_of_uFF=1+startFFindex
  end if
 end function index_of_uFF
 
  !!! provide the index of grid associated with lpPDF(hadron)
 function index_of_lpPDF(hadron)
  integer,intent(in)::hadron
  integer::index_of_lpPDF
  integer::i
  
  if(num_of_lpPDFs==0) then
    write(*,*) ErrorString('no lpPDFs are initialized',moduleName)
    stop
  else
    do i=1,num_of_lpPDFs
      if(enumeration_of_lpPDFs(i)==hadron) then
	index_of_lpPDF=i+startlpPDFindex
	return
      end if
    end do
    !!! if we exit from the loop it means index is not found
    write(*,*) ErrorString('lpPDF is not found',moduleName)
    write(*,"('no lpPDF for hadron ',I3,' is initialized. Set PDF for hadron (',I3,')')") &
		  hadron,enumeration_of_lpPDFs(1)
      index_of_lpPDF=1+startlpPDFindex
  end if
 end function index_of_lpPDF
 
 !!! set a different replica number for PDF.
 subroutine QCDinput_SetPDFreplica(rep)
 integer:: rep
  call InitPDFM(1,rep)
 end subroutine QCDinput_SetPDFreplica
 
!!number of active flavor at given mu
function activeNf(mu)
  integer::activeNf
  real(dp)::mu
  
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
 real(dp)::as,Q,alphasPDF
 
 As=alphasPDF(Q)/pix4
!  as=1d0/(2d0*23d0/3d0*Log(Q/0.08782708014552364d0))  !!! Nf=5 LO solution
 
 end function As

 !!!!array of x times PDF(x,Q) for hadron 'hadron'
 !!! unpolarized PDF used in uTMDPDF
 !!!! array is (-5:5) (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
 function xPDF(x,Q,hadron)
      real(dp),intent(in) :: x,Q
      integer,intent(in):: hadron
      real(dp), dimension(-5:5):: xPDF
      real(dp), dimension (-6:6)::inputPDF
      
      call evolvePDFM(index_of_uPDF(hadron),x,Q,inputPDF)
      
      xPDF=inputPDF(-5:5)
      
      
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
      real(dp), dimension (-6:6)::inputFF
      
      call evolvePDFM(index_of_uFF(hadron),x,Q,inputFF)
      
      xFF=inputFF(-5:5)
      
  end function xFF
 
  !!!!array of x times PDF(x,Q) for hadron 'hadron'
  !!! unpolarized PDF used in lpTMDPDF
 !!!! array is (-5:5) (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
 function x_lp_PDF(x,Q,hadron)
      real(dp),intent(in) :: x,Q
      integer,intent(in):: hadron
      real(dp), dimension(-5:5):: x_lp_PDF
      real(dp), dimension (-6:6)::inputPDF
      
      call evolvePDFM(index_of_lpPDF(hadron),x,Q,inputPDF)
      
      x_lp_PDF=inputPDF(-5:5)
      
  end function x_lp_PDF
 
end module QCDinput
