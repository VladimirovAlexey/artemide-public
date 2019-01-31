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
    name="BASE MODEL 1"
    
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

   real*8:: bb,c1,c0
   real*8::ZZ,AA
   bb=bT/lambdaNP(1)
   AA=(1-x)*lambdaNP(2)+x*lambdaNP(3)
   c0=(AA-0.5d0)*bb
   c1=(AA+0.5d0)*bb
   ZZ=1d0
   
   if(AA<0.0001d0) then
   FNP0=1d6
   else if(c1>50d0) then
   FNP0=(1-lambdaNP(4))*Exp(c0-c1)*Abs(ZZ)+lambdaNP(4)*Exp(-bT**2/lambdaNP(5)**2)
   else
   FNP0=(1-lambdaNP(4))*Cosh(c0)/Cosh(c1)*Abs(ZZ)+lambdaNP(4)*Exp(-bT**2/lambdaNP(5)**2)
   end if


  
  FNP=FNP0*(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
!   FNP=(/FNP1,FNP1,FNP1,FNP1,FNP1,FNP0,FNP0,FNP0,FNP1,FNP1,FNP1/)
!   FNP=(/FNP2,FNP2,FNP2,FNP2,FNP2,FNP2,FNP0,FNP1,FNP2,FNP2,FNP2/)
  end function FNP
  
    !!!!This function is the mu(x,b), which is used inside the OPE
  function mu_OPE(x,bt)
  real*8::bt,mu_OPE,x
  !mu_OPE=C0_const*SQRT(1+bT**2)/bT+1d0
  mu_OPE=C0_const*1d0/bT+2d0
  
  if(mu_OPE>1000d0) then
    mu_OPE=1000d0
  end if
  end function mu_OPE
  
 function ReplicaParameters(rep)
 integer::rep
 real*8::ReplicaParameters(1:1)
 ReplicaParameters=(/0d0/) 
 end function ReplicaParameters