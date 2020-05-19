!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Model for unpolarized TMD PDF  BSV19  [1902.08474]
!
!				A.Vladimirov (11.07.2019)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! USER DEFINED FUNCTIONS   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !!!!!! Write nessecery model intitialization.
  subroutine ModelInitialization()
    logical::file_exists
    integer::num_NP,i,j  

     name="SV19"
     write(*,*) 'Model SV19 is used. Please, cite 1902.08474'
     
  end subroutine ModelInitialization
  
  
  !!! This is  non-pertrubative function
  !!! non=pertrubative parameters are lambdaNP()
  !!! x-- is the bjorken variable of TMD
  !!! z-- is convolution variable
  function FNP(x,z,bT,hadron,lambdaNP)
  real*8::x,z,bT
  real*8,dimension(-5:5)::FNP
  real*8::FNP0
  integer::hadron
  real*8,intent(in)::lambdaNP(:)

   real*8::bb,w1,w2,w3

   
!!!    SV19 MODEL
    bb=bT**2
    w1=lambdaNP(1)*(1-x)+x*lambdaNP(2)+x*(1-x)*lambdaNP(5)
    w2=lambdaNP(3)*x**lambdaNP(4)+lambdaNP(6)
   
    if(w2<0d0 .or. w1<0d0) then !!! case of negative power, we return absolutely incorrect expression.
      if(bT<1d0) then
	FNP0=-1d0
      else
	FNP0=0d0
      end if
    else
    FNP0=Exp(-w1*bb/sqrt(1+w2*bb))
    end if

!   bb=bT**2
!   w1=x*lambdaNP(2)+(1-x)*lambdaNP(3)
!   w2=lambdaNP(4)*x**lambdaNP(6)+lambdaNP(5)*(1-x)**lambdaNP(7)
!   FNP0=Exp(-w1*bb/Sqrt(1+w1**2*bb/lambdaNP(1)**2))*(1d0+w2*bb/(1d0+lambdaNP(1)**2*bb))  

!   FNP0=Exp(-lambdaNP(1)*bb)
   
  FNP=FNP0*(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)

  end function FNP
  
  !!!! This is the function b* that enter the logarithms of coefficient function
  !!!! at small-b it should be ~b to match the collinear regime
  !!!! at large-b it is a part of model
  !!!! NOTE: if it is lambda-dependent, the grid will be recalculate each reset of lambdaNP
  function bSTAR(bT,lambdaNP)
    real*8,intent(in)::bT
    real*8,intent(in)::lambdaNP(:)
    real*8::bSTAR
    
    bSTAR=bT/sqrt(1d0+(bT/C0_const)**2)
    
  end function bSTAR
  
  !!!!This function is the mu(x,b), which is used inside the OPE
  function mu_OPE(x,bt)
  real*8::bt,mu_OPE,x
  
  !mu_OPE=C0_const*1d0/bT+2d0
  mu_OPE=C0_const/bt*Sqrt(1+bT**2/(C0_const)**2)
  
  if(mu_OPE>1000d0) then
    mu_OPE=1000d0
  end if
  end function mu_OPE
  
  
  !!!! if the option UseComposite TMD is OFF, this function is ignored
  !!!! If the option UseComposite TMD is ON,
  !!!! than the TMD for hardon is build as TMD(hadron)=Sum_c CA(h,c) TMD(c)
  !!!! where h=hadron, CA=coefficientArray
  !!!! coefficientArray real*8 list of coefficeints
  !!!! includeArray is logical array list (true=TMD(c) is computed, false TMD(c) ignored)
  subroutine GetCompositionArray(hadron,lambdaNP,includeArray,coefficientArray)  
  real*8,intent(in)::lambdaNP(:)
  integer::hadron
  logical,allocatable,intent(out)::includeArray(:)
  real*8,allocatable,intent(out)::coefficientArray(:)
  
  integer::k
  integer,parameter::kStart=7!!! 
  real*8::alpha
   
!    allocate(includeArray(1:1))
!    allocate(coefficientArray(1:1))
!
!    includeArray=(/.true./)
!    coefficientArray=(/1d0/)

!   HERA CASE
   allocate(includeArray(1:29))
   allocate(coefficientArray(1:29))
   
   includeArray(1)=.true.
   coefficientArray(1)=1d0
   do k=1,14
    alpha=lambdaNP(kStart+k)
    if(lambdaNP(kStart+k)>=0d0) then
      includeArray(2*k)=.true.
      includeArray(2*k+1)=.false.
      coefficientArray(1)=coefficientArray(1)-alpha
      coefficientArray(2*k)=alpha
      coefficientArray(2*k+1)=0d0
    else
      includeArray(2*k)=.false.
      includeArray(2*k+1)=.true.
      coefficientArray(1)=coefficientArray(1)+alpha
      coefficientArray(2*k)=0d0
      coefficientArray(2*k+1)=-alpha
    end if
   
   end do
   
  end subroutine GetCompositionArray
  
 !!! reads the replica-file 
 !!! -1 is suggested for initialization replica
 !!! 0 is the mean reaplics
 !!! 1 -- 100 replicas
 function ReplicaParameters(rep)
 integer::rep
 real*8::ReplicaParameters(1:lengthNP)
 integer::i
 
 write(*,*) warningstring("set model replica via artemide-control module",moduleName)
 write(*,*) warningstring("some generic NP values returned",moduleName)
 ReplicaParameters(1)=1d0
 do i=2,lengthNP
  ReplicaParameters(1)=0.001d0
 end do
 
 end function ReplicaParameters
  
