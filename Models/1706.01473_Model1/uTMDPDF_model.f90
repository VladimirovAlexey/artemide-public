!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Model for unpolarized TMD PDF
!
!			corresponds to model 1
!			FNP=Cosh((l1/l2-l1/2)b)/Cosh((l1/l2+l1/2)b)
!			muOPE=C0/b+2
!
!			Requres two NP parameters
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
    if(lambdaNPlength<2) then
      write(*,*) 'arTeMiDe.uTMDPDF-model: Number NP parameters for TMDPDF has less then 2'
      write(*,*) 'Evaluation STOP'
      stop
    end if
    
    SELECT CASE(order_global)
      CASE(0)
	lambdaNP(1)=0.112d0
	lambdaNP(2)=0.828d0
      CASE(1)
	lambdaNP(1)=0.179d0
	lambdaNP(2)=0.354d0
      CASE(2)
	lambdaNP(1)=0.246d0
	lambdaNP(2)=0.311d0
    END SELECT
  end subroutine ModelInitialization
  
  
  !!! This is (flavor-independent) non-pertrubative function
  !!! non=pertrubative parameters are lambdaNP()
  !!! if renormalon correction is on, the parameter lambdaNP(2) stays infront of it.
  function FNP(z,bT,hadron)
  real*8::z,bT
  real*8,dimension(-5:5)::FNP
  real*8::y2,FNP0,FNP1
  real*8::a,b
  integer::hadron
  
  !!!!MODEL 0
!   FNP=EXP(-lambdaNP(1)*bT**2)
  
  !!!!MODEL 1
  a=lambdaNP(2)/lambdaNP(1)
  b=lambdaNP(1)/2d0
  
  If((a==0).or.(b==0)) then
  FNP0=0
  else
  If((a+b)*bT>200d0) then
   If(a>b) then
    FNP0=EXP(-lambdaNP(1)*bT)
   else
    FNP0=EXP(-2d0*a*bT)
   end if
  else
   FNP0=COSH((a-b)*bT)/COSH((a+b)*bT)
  end if
  end if
  
!   if(FNP0==0) write(*,*) bT, lambdaNP
   !!!! MODEL 2
!    FNP0=EXP(-lambdaNP(3)*z*bT**2/SQRT(1+(lambdaNP(3)*z*bT/lambdaNP(1))**2))
  
  FNP=FNP0*(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
!   FNP=(/FNP1,FNP1,FNP1,FNP1,FNP1,FNP0,FNP0,FNP0,FNP1,FNP1,FNP1/)
  end function FNP
  
    !!!!This function is the mu(x,b), which is used inside the OPE
  function mu_OPE(x,bt)
  real*8::bt,mu_OPE,x
  !mu_OPE=C0_const*SQRT(1+bT**2)/bT+1d0
  mu_OPE=C0_const*1d0/bT+2d0
  !write(*,*) as(1000d0),as(5000d0),as(10000d0),as(15000d0),as(20000d0),as(25000d0)
  if(mu_OPE>1000d0) then
    mu_OPE=1000d0
  end if
  end function mu_OPE
