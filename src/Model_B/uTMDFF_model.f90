!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Model for unpolarized TMD FF
!
!			corresponds to model 1
!			FNP=Cosh((l1/l2-l1/2)b)/Cosh((l1/l2+l1/2)b)
!			muOPE=C0/b+2
!
!			Requres two NP parameters (initated by PDF values)
!			Initialized by best values with bb* model for evolution
!
!				A.Vladimirov (25.04.2018)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! USER DEFINED FUNCTIONS   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
 
  !!!!!! Write nessecery model intitialization.
  subroutine ModelInitialization()  
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
  
  real*8::bb,w1,w2
  
  bb=bT**2/x**2
  
!   if(hadron==1) then
    w1=lambdaNP(1)*x+lambdaNP(2)*(1d0-x)
    w2=lambdaNP(3)
    FNP0=Exp(-bb*w1/sqrt(1d0+w2*bb))*(1+lambdaNP(4)*bb)
!   else
!     w1=lambdaNP(5)*x+lambdaNP(6)*(1d0-x)
!     w2=lambdaNP(7)
!     FNP0=Exp(-bb*w1/sqrt(1d0+w2*bb))*(1+lambdaNP(8)*bb)
!   end if
  
  FNP=FNP0*(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
!   FNP=(/FNP1,FNP1,FNP1,FNP1,FNP1,FNP0,FNP0,FNP0,FNP1,FNP1,FNP1/)
!   FNP=(/FNP2,FNP2,FNP2,FNP2,FNP2,FNP2,FNP0,FNP1,FNP2,FNP2,FNP2/)
  end function FNP
  
  !!!! This is the function b* that enter the logarithms of coefficient function
  !!!! at small-b it should be ~b to match the collinear regime
  !!!! at large-b it is a part of model
  !!!! NOTE: if it is lambda-dependent, the grid will be recalculate each reset of lambdaNP
  function bSTAR(bT,lambdaNP)
    real*8,intent(in)::bT
    real*8,intent(in)::lambdaNP(:)
    real*8::bSTAR
    
    bSTAR=bT/sqrt(1d0+(bT/500d0)**2)
    
  end function bSTAR
  
  !!!!This function is the mu(x,b), which is used inside the OPE
  function mu_OPE(x,bT)
    real*8,intent(in)::bT,x
    real*8::mu_OPE
    
    mu_OPE=C0_const*x/bT+2d0
    
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
   
   allocate(includeArray(1:1))
   allocate(coefficientArray(1:1))
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