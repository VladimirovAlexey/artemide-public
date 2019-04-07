!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 2.0
!
! Interface module to the user defined alpha-s, PDF, FF, etc.
! Could be interfaced to LHAPDF
!
!	ver.2.0: 28.03.2019 AV
!				A.Vladimirov
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module QCDinput

implicit none

private



public::QCDinput_Initialize,As,QCDinput_IsInitialized
public::xPDF,xFF
public:: QCDinput_SetPDFreplica

!Current version of module
 character (len=5),parameter :: version="v2.00"
!--- general
logical:: started=.false.
integer::outputLevel

real*8,public::mCHARM,mBOTTOM

!---uPDFs
integer::num_of_uPDFs,startPDFindex
integer,allocatable::enumeration_of_uPDFs(:)

!---uFFs
integer::num_of_uFFs,startFFindex
integer,allocatable::enumeration_of_uFFs(:)

 contains 
 
 
 function QCDinput_IsInitialized()
  logical::QCDinput_IsInitialized
  QCDinput_IsInitialized=started 
 end function QCDinput_IsInitialized
 
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
 
 !------------------------------Functions below are to be changed by user (if needed)
 subroutine QCDinput_Initialize(file,prefix)
  character(len=*)::file
  character(len=*),optional::prefix
  character(len=300)::path,line
  character(len=64),allocatable::names(:)
  integer::i
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
    call MoveTO(51,'*p2  ')
    read(51,*) outputLevel    
    if(outputLevel>2) write(*,*) '--------------------------------------------- '
    if(outputLevel>2) write(*,*) 'artemide.QCDinput: initialization started ... '
    
    !!! Search for QCDinput initialization options
    call MoveTO(51,'*1   ')
    !!! Search for parameters
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) mCHARM
    call MoveTO(51,'*p2  ')
    read(51,*) mBOTTOM
    
    if(outputLevel>2) write(*,*) '    mass of charm:  ',mCHARM
    if(outputLevel>2) write(*,*) '    mass of bottom: ',mBOTTOM
    
    
    !!! Search for uPDF initialization options
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
    
    !!! Search for uFF initialization options
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
 
  CLOSE (51, STATUS='KEEP') 
  
  if(outputLevel>2)	write(*,*)'QCDinput succesfully initialized.'
  started=.true.
 end subroutine QCDinput_Initialize
 
 !!! provide the index of grid associated with uPDF(hadron)
 function index_of_uPDF(hadron)
  integer,intent(in)::hadron
  integer::index_of_uPDF
  integer::i
  
  if(num_of_uPDFs==0) then
    write(*,*) 'artemide.QCDinput: CRITICAL ERROR: no uPDFs are initialized'
    stop
  else
    do i=1,num_of_uPDFs
      if(enumeration_of_uPDFs(i)==hadron) then
	index_of_uPDF=i+startPDFindex
	return
      end if
    end do
    !!! if we exit from the loop it means index is not found
      write(*,"('artemide.QCDinput: ERROR: no uPDF for hadron ',I3,'is initialized. Set PDF for hadron (',I3,')')") &
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
    write(*,*) 'artemide.QCDinput: CRITICAL ERROR: no uFFs are initialized'
    stop
   else
    do i=1,num_of_uFFs
      if(enumeration_of_uFFs(i)==hadron) then
	index_of_uFF=i+startFFindex
	return
      end if
    end do
    !!! if we exit from the loop it means index is not found
      write(*,"('artemide.QCDinput: ERROR: no uFF for hadron ',I3,'is initialized. Set FF for hadron (',I3,')')") &
		  hadron,enumeration_of_uPDFs(1)
      index_of_uFF=1+startFFindex
  end if
 end function index_of_uFF
 
 !!! set a different replica number for PDF.
 subroutine QCDinput_SetPDFreplica(rep)
 integer:: rep
  call InitPDF(rep)
 end subroutine QCDinput_SetPDFreplica
 
 !!!!alphas(Q)/4pi
 !!! NOT FORGET 4 PI !!!
 function As(Q)
 real*8::as,Q,alphasPDF
 real*8,parameter::Pi4=12.566370614359172d0
 
 As=alphasPDF(Q)/Pi4
 
 end function As

 !!!!array of x times PDF(x,Q) for hadron 'hadron'
 !!!! array is (-5:5) (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
 function xPDF(x,Q,hadron)
      real*8,intent(in) :: x,Q
      integer,intent(in):: hadron
      real*8, dimension(-5:5):: xPDF
      real*8, dimension (-6:6)::inputPDF
      
      call evolvePDFM(index_of_uPDF(hadron),x,Q,inputPDF)
      
      xPDF=inputPDF(-5:5)
      
  end function xPDF
  
  
    !!!! return x*F(x,mu)
  !!!! enumeration of flavors
  !!!!  f = -5,-4, -3,  -2,  -1,0,1,2,3,5,4
  !!!!    = bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b
  !!!! enumeration of hadrons 
  function xFF(x,Q,hadron)
      integer,intent(in) :: hadron
      real*8,intent(in) :: x,Q
      real*8,dimension(-5:5):: xFF
      real*8, dimension (-6:6)::inputFF
      
      call evolvePDFM(index_of_uFF(hadron),x,Q,inputFF)
      
      xFF=inputFF(-5:5)
  end function xFF
 
end module QCDinput