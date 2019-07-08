!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Model interface for linearly-polarized gluon TMD PDF
!
!			the files of model it-self are defined by user in /model
!
!				A.Vladimirov (13.06.2018)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! USER DEFINED FUNCTIONS   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
module lpTMDPDF_model
implicit none

private
  character(50)::name='NONAME'
  character (len=5),parameter :: version="v2.01"
  
  integer::lengthNP

public::ModelInit,FNP,mu_OPE,TestFNP,TestMU,GiveReplicaParameters,TestbSTAR,bSTAR

  INCLUDE 'Tables/NumConst.f90'

contains

  INCLUDE 'Model/lpTMDPDF_model.f90'
 
  !!!!!! Write nessecery model intitialization.
  subroutine ModelInit(outputlevel,length)
    integer:: outputlevel,length
    
    if(outputLevel>1) write(*,*) 'linear-pol. gluon TMDPDF Model initialization (ver.',version,') ...'
    
    lengthNP=length
    
    call ModelInitialization()
    if(outputlevel>1) write(*,*) 'linear-pol. gluon PDF Model name: ',name,'  ... initialized'
  end subroutine ModelInit
  
  !!! test FNP for z-dependance
  function TestFNP(hadronsInGRID,lambdaNPlength)
  logical::TestFNP
  real*8::xR,bR
  real*8,dimension(-5:5)::test1,test2
  integer::h,numberOfHadrons,i,j,lambdaNPlength
  integer,intent(in)::hadronsInGRID(:)
  real*8,allocatable::lambdaNP(:)
  
  allocate(lambdaNP(1:lambdaNPlength))
  
  numberOfHadrons=size(hadronsInGRID)
  
  call RANDOM_NUMBER(lambdaNP)
  lambdaNP=3d0*lambdaNP
    
    TestFNP=.false.
    do h=1,numberOfHadrons
    do i=1,10
        !!! generate some random input
	call RANDOM_NUMBER(xR)
	if(xR>0.99d0) xR=xR/2d0
	if(xR<0.00001d0) xR=0.0001d0+xR
  
	call RANDOM_NUMBER(bR)  
	bR=3d0*bR
    
      test1=FNP(xR,0.9d0,bR,hadronsInGRID(h),lambdaNP)
      do j=1,10
	test2=FNP(xR,0.8d0/REAL(j),bR,hadronsInGRID(h),lambdaNP)
	if(SUM(ABS(test1-test2))>1d-10) then
	  TestFNP=.true.
	  exit
	end if	
      end do
    end do
    end do
    
    deallocate(lambdaNP)
    
  end function TestFNP
  
  
  !!! test MU for x-dependance
  function TestMU()
  logical::TestMU
  real*8::xR,bR
  real*8::test1,test2
  integer::i,j
  TestMU=.false.
    do i=1,10
    call RANDOM_NUMBER(bR)  
    bR=3d0*bR
        !!! generate some random input
    call RANDOM_NUMBER(xR)
    if(xR>0.99d0) xR=xR/2d0
    if(xR<0.00001d0) xR=0.0001d0+xR
    test1=mu_OPE(xR,bR)
      
    !!! generate some random input
    call RANDOM_NUMBER(xR)
    if(xR>0.99d0) xR=xR/2d0
    if(xR<0.00001d0) xR=0.0001d0+xR
    test2=mu_OPE(xR,bR)
    
    if(ABS(test1-test2)>1d-10) then
	  TestMU=.true.
	  exit
    end if	
    end do
  end function TestMU
  
  
    !!! test bSTAR for lambdaNP-dependance
  !!! .true. = lambda-dependent
  function TestbSTAR(lambdaNPlength)
    integer::lambdaNPlength,i,j,k
    logical::TestbSTAR
    real*8::bT,dummy1,dummy2
    real*8,allocatable::lambdaNP(:)
    
    allocate(lambdaNP(1:lambdaNPlength))
    TestbSTAR=.false.
    do i=1,lambdaNPlength
      lambdaNP=0.5d0+0d0*lambdaNP
      do k=0,3
	bT=0.1d0+k
	dummy1=bSTAR(bT,lambdaNP)
	
	do j=1,3
	  lambdaNP(i)=0.9d0*j
	  dummy2=bSTAR(bT,lambdaNP)
	  if(abs(dummy2-dummy1)>1d-10) then
	    TestbSTAR=.true.
	    exit
	  end if
	end do
      end do
    end do
    
  end function TestbSTAR
  
  
  !!! transfer the replica paramters to main function
  function GiveReplicaParameters(num)
  real*8::GiveReplicaParameters(1:lengthNP)
  integer::num,i,ss
  ss=size(ReplicaParameters(num))
  
  if(lengthNP<ss) then
   write(*,*) 'ERROR: arTeMiDe.lpTMDPDF_MODEL: number of non-pertrubative parameters for model ',name,&
		    ' should be >= ',ss
   write(*,*) 'EVALUATION STOP'
   stop
  else if(lengthNP==ss) then
   GiveReplicaParameters=ReplicaParameters(num)
  else
   GiveReplicaParameters(1:ss)=ReplicaParameters(num)
   do i=ss+1,lengthNP
    GiveReplicaParameters(i)=0d0
   end do
   
  end if
  
  end function GiveReplicaParameters

end module lpTMDPDF_model
  