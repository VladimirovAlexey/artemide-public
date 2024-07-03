!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Model for unpolarized TMD FF  SV19  [1912.06532]
!
!				A.Vladimirov (11.07.2019)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module uTMDFF_model
use aTMDe_Numerics
use IO_functions
implicit none

private

!!!!!------------------------------------------------------------------------------------
!!!!! These functions MUST defined in module  !!
!!!!!
!!!!! 1) The subroutine is called during the initialization of TMD-module
!!!!!    arg=array of initial NP-parameters
public:: ModelInitialization
!!!!! 2) The subroutine that is called on reset of NP-parameters in TMD-module
!!!!!    arg=array of new NP-parameters
public:: ModelUpdate
!!!!! 3) Function which returns FNP function
!!!!!    arg=(x,z,b,hadron,lambdaNP) with x=x_Bj for TMD (real_dp), z=convolution variable(real_dp), 
!!!!!    b=transverse distance(real_dp), hadron=number of the hadron in grid(integer)
!!!!!    lambdaNP = array of NP parameters (real_dp(:))
real(dp),public,dimension(-5:5):: FNP
!!!!! 4) Function which returns the value of b used as argument of convolution integrals
!!!!!    arg=(b,lambdaNP) with b=transverse distance(real_dp), lambdaNP = array of NP parameters (real_dp(:))
real(dp),public:: bSTAR
!!!!! 5) Function which returns the scale of matching (OPE scale)
!!!!!    arg=(z,bt) with z=convolution variable(real_dp), b=transverse distance(real_dp)
real(dp),public:: muOPE
!!!!!------------------------------------------------------------------------------------

real(dp),allocatable::NPparam(:)

contains  

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! USER DEFINED FUNCTIONS   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
!!!!!! Write nessecery model intitialization.
subroutine ModelInitialization(lengthNP)
    integer,intent(in)::lengthNP
    !!!!!! here are the initial parameters!!
    allocate(NPparam(1:lengthNP))
    NPparam=0._dp
    
    write(*,*) color(">>>  The model for uTMDFF is SV19. Please, cite [1912.06532]   <<<",c_cyan)
    
end subroutine ModelInitialization

!!!!!! Write nessecery model update (e.g. save current NP-parameters)
!!!!!! newNPParams is the new NP-array
subroutine ModelUpdate(newNPParams)  
    real(dp),intent(in):: newNPParams(:)
    
    NPparam=newNPParams !! save new vector of NP-parameters

end subroutine ModelUpdate
  
!!! This is  non-pertrubative function
!!! non=pertrubative parameters are lambdaNP()
!!! x-- is the bjorken variable of TMD
function FNP(x,bT,hadron,lambdaNP)
    real(dp),intent(in)::x,bT
    integer,intent(in)::hadron
    real(dp),intent(in)::lambdaNP(:)

    real(dp)::FNP0  
    real(dp)::bb,w1,w2,M

    bb=bT**2/x**2
    
!     if(hadron==1) then
        w1=lambdaNP(1)*x+lambdaNP(2)*(1d0-x)
        w2=lambdaNP(3)
        FNP0=Exp(-bb*w1/sqrt(1d0+w2*bb))*(1+lambdaNP(4)*bb)
        
!     else
!         w1=lambdaNP(5)*x+lambdaNP(6)*(1d0-x)
!         w2=lambdaNP(7)
!         FNP0=Exp(-bb*w1/sqrt(1d0+w2*bb))*(1+lambdaNP(8)*bb)
!     end if

    FNP=FNP0*(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/) 

end function FNP
  
!!!! This is the function b* that enters the logarithms of coefficient function
!!!! at small-b it should be ~b to match the collinear regime
!!!! at large-b it is a part of model
!!!! x -- is the global x for TMDPDF,
!!!! y -- is the convolution variable in the definition \int dy/y C(y) PDF(x/y)
pure function bSTAR(bT,x,y)
    real(dp),intent(in)::bT,x,y

    bSTAR=bT/sqrt(1d0+(bT/500d0)**2)
end function bSTAR

!!!!This function is the mu(x,b), which is used inside the OPE
!!!! x -- is the global x for TMDPDF,
!!!! y -- is the convolution variable in the definition \int dy/y C(y) PDF(x/y)
!!!! c4-- is the scale variation variable
pure function muOPE(bt,x,y,c4)
    real(dp),intent(in)::bt,x,y,c4

    muOPE=C0_const*x*c4/bT+2d0

    if(muOPE>1000d0) then
        muOPE=1000d0
    end if
end function muOPE

end module uTMDFF_model
