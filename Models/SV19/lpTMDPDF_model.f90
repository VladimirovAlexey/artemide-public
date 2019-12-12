!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Model for linearly polarized gluons
!
!			corresponds to 
!				A.Vladimirov (13.06.2019)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! USER DEFINED FUNCTIONS   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !!!!!! Write nessecery model intitialization.
  subroutine ModelInitialization()  
    name="NONAME"
    
    !!write(*,*) 'Model BSV19.bFIT is used. Please, cite 1902.08474'
    
  end subroutine ModelInitialization
  
  
  !!! This is  non-pertrubative function
  !!! non=pertrubative parameters are lambdaNP()
  !!! x-- is the bjorken variable of TMD
  !!! z-- is convolution variable
  function FNP(x,z,bT,hadron,lambdaNP)
    real*8,intent(in)::x,z,bT
    integer,intent(in)::hadron
    real*8,intent(in)::lambdaNP(:)
    real*8,dimension(-5:5)::FNP
    
    real*8::FNP0,p1,p2
    
    !!! NOT USED
    FNP0=Exp(-0.5d0*bT**2)


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
    
    bSTAR=bT/sqrt(1d0+(bT/500d0)**2)
    
  end function bSTAR
  
  !!!!This function is the mu(x,b), which is used inside the OPE
  function mu_OPE(x,bT)
    real*8,intent(in)::bT,x
    real*8::mu_OPE
    
    mu_OPE=C0_const*1d0/bT+2d0
    
    if(mu_OPE>1000d0) then
      mu_OPE=1000d0
    end if
  end function mu_OPE
  
   !!! this is the table of replica prameters extracted in fit BSV19.
 !!! -2 is suggested for initialization replica
 !!! -1 is the best fit
 !!! 0 is the mean reaplics
 !!! 1 -- 100 replicas
 function ReplicaParameters(rep)
 integer::rep
 real*8::ReplicaParameters(1:6)    
 
 Write(*,*) 'THERE IS NOT YET FIT OF LINEARLY POLARIZED GLUONS'
 stop
 
 end function ReplicaParameters
  