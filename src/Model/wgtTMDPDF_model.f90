!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Model function for worm-gear-T TMDPDF
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module wgtTMDPDF_model
use aTMDe_Numerics
use IO_functions
implicit none

private

!!!!!------------------------------------------------------------------------------------
!!!!! These functions MUST defined in module  !!
!!!!!
!!!!! 1) The subroutine is called during the initialization of TMD-module
public:: ModelInitialization
!!!!! 2) The subroutine that is called on reset of NP-parameters in TMD-module
!!!!!    arg=array of new NP-parameters
public:: ModelUpdate
!!!!! 3) Function which returns FNP function
!!!!!    arg=(x,b,hadron,lambdaNP) with x=x_Bj for TMD (real_dp), 
!!!!!    b=transverse distance(real_dp), hadron=number of the hadron in grid(integer)
!!!!!    lambdaNP = array of NP parameters (real_dp(:))
real(dp),public,dimension(-5:5):: FNP
!!!!! 4) Function which returns the value of b used as argument of convolution integrals
!!!!!    arg=(bT,x,y) with br=transverse distance(real_dp), x being Bjorken x (real_dp), and y being convolution variable (real_dp)
real(dp),public:: bSTAR
!!!!! 5) Function which returns the scale of matching (OPE scale)
!!!!!    arg=(bT,x,y,c4) with bT, x, y same as in bSTAR, and c4(real_dp) is scale variation parameter
real(dp),public:: muOPE
!!!!! 6) Function which returns FNP function that accompanies the tw3-part of OPE
!!!!!    arg=(x,b,hadron,lambdaNP) with x=x_Bj for TMD (real_dp),
!!!!!    b=transverse distance(real_dp), hadron=number of the hadron in grid(integer)
!!!!!    lambdaNP = array of NP parameters (real_dp(:))
real(dp),public,dimension(-5:5):: FNP_tw3
!!!!! 7) Function which returns the value of b used as argument of convolution integrals of tw3-part
!!!!!    arg=(bT,x,y) with br=transverse distance(real_dp), x being Bjorken x (real_dp), and y being convolution variable (real_dp)
real(dp),public:: bSTAR_tw3
!!!!! 8) Function which returns the scale of matching (OPE scale) in the tw3-part
!!!!!    arg=(bT,x,y,c4) with bT, x, y same as in bSTAR, and c4(real_dp) is scale variation parameter
real(dp),public:: muOPE_tw3
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
    
    !!write(*,*) color(">>>  The model for wgtTMDPDF for ART23   <<<",c_cyan)
    
end subroutine ModelInitialization

!!!!!! Write nessecery model update (e.g. save current NP-parameters)
!!!!!! newNPParams is the new NP-array
subroutine ModelUpdate(newNPParams)  
    real(dp),intent(in):: newNPParams(:)
    
    NPparam=newNPParams !! save new vector of NP-parameters

end subroutine ModelUpdate

!!! This is  non-perturbative function that multiplies tw2-part
!!! non=perturbative parameters are lambdaNP()
!!! x-- is the bjorken variable of TMD
function FNP(x,bT,hadron,lambdaNP)
  real(dp),intent(in)::x,bT    
  integer,intent(in)::hadron
  real(dp),intent(in)::lambdaNP(:)
    real(dp)::FNP0

    FNP0=Exp(-0.5d0*bT**2)

    FNP=FNP0*(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
end function FNP
  
!!!! This is the function b* that enters the logarithms of coefficient function in tw2-part
!!!! at small-b it should be ~b to match the collinear regime
!!!! at large-b it is a part of model
!!!! x -- is the global x for TMDPDF,
!!!! y -- is the convolution variable in the definition \int dy/y C(y) PDF(x/y)
pure function bSTAR(bT,x,y)
    real(dp),intent(in)::bT,x,y

    bSTAR=bT/sqrt(1d0+(bT/500d0)**2)
end function bSTAR
  
!!!!This function is the mu(x,b), which is used inside the OPE in tw2-part
!!!! x -- is the global x for TMDPDF,
!!!! y -- is the convolution variable in the definition \int dy/y C(y) PDF(x/y)
!!!! c4-- is the scale variation variable
pure function muOPE(bt,x,y,c4)
    real(dp),intent(in)::bt,x,y,c4

    muOPE=C0_const*c4/bT+2d0
    
    if(muOPE>1000d0) then
        muOPE=1000d0
    end if
end function muOPE

!!! This is  non-perturbative function that multiplies tw3-part
!!! non=perturbative parameters are lambdaNP()
!!! x-- is the bjorken variable of TMD
function FNP_tw3(x,bT,hadron,lambdaNP)
  real(dp),intent(in)::x,bT
  integer,intent(in)::hadron
  real(dp),intent(in)::lambdaNP(:)
  real(dp)::FNP0

    FNP0=Exp(-0.5d0*bT**2)

    FNP_tw3=FNP0*(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
end function FNP_tw3

!!!! This is the function b* that enters the logarithms of coefficient function in tw3-part
!!!! at small-b it should be ~b to match the collinear regime
!!!! at large-b it is a part of model
!!!! x -- is the global x for TMDPDF,
!!!! y -- is the convolution variable in the definition \int dy/y C(y) PDF(x/y)
pure function bSTAR_tw3(bT,x,y)
    real(dp),intent(in)::bT,x,y

    bSTAR_tw3=bT/sqrt(1d0+(bT/500d0)**2)
end function bSTAR_tw3

!!!!This function is the mu(x,b), which is used inside the OPE in tw3-part
!!!! x -- is the global x for TMDPDF,
!!!! y -- is the convolution variable in the definition \int dy/y C(y) PDF(x/y)
!!!! c4-- is the scale variation variable
pure function muOPE_tw3(bt,x,y,c4)
    real(dp),intent(in)::bt,x,y,c4

    muOPE_tw3=C0_const*c4/bT+2d0

    if(muOPE_tw3>1000d0) then
        muOPE_tw3=1000d0
    end if
end function muOPE_tw3

end module wgtTMDPDF_model
