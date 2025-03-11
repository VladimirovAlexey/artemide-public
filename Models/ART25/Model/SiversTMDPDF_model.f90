!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Model for Sivers TMD PDF [20??.????]
!
!				A.Vladimirov (21.05.2020)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module SiversTMDPDF_model
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
    
    write(*,*) color(">>>  The model for Sivers function is BPV20. Please, cite [20??.????]   <<<",c_cyan)
    
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

    real(dp)::bProfile
    real(dp)::FNPu,FNPd,FNPs,FNPsea,Normu,Normd,Normsea,YY

    !!! profile in b is common for all (5 parameters)    
    YY=(lambdaNP(1)+x*lambdaNP(2))*(bT**2)/sqrt(1d0+Abs(lambdaNP(3))*x**2*bT**2)
    bProfile=exp(-YY)
    !bProfile=1d0/cosh((lambdaNP(1)+x**2*lambdaNP(2))*bT)
    
    !!! u-quark(3 parameters)
    Normu=(3d0+lambdaNP(7)+lambdaNP(8)*(1+lambdaNP(7)))/((lambdaNP(7)+1d0)*(lambdaNP(7)+2d0)*(lambdaNP(7)+3d0))    
    FNPu=lambdaNP(6)*(1-x)*x**lambdaNP(7)*(1+lambdaNP(8)*x)/Normu
    !!! d-quark(3 parameters)
    Normd=(3d0+lambdaNP(10)+lambdaNP(11)*(1+lambdaNP(10)))/((lambdaNP(10)+1d0)*(lambdaNP(10)+2d0)*(lambdaNP(10)+3d0))    
    FNPd=lambdaNP(9)*(1-x)*x**lambdaNP(10)*(1+lambdaNP(11)*x)/Normd
    !!! sea-quark(3 parameters)
    Normsea=1d0/((lambdaNP(13)+1d0)*(lambdaNP(13)+2d0))    
    FNPs=lambdaNP(12)*(1-x)*x**lambdaNP(13)/Normsea
    FNPsea=lambdaNP(14)*(1-x)*x**lambdaNP(13)/Normsea
    
    FNP=bProfile*(/0d0,0d0,FNPsea,FNPsea,FNPsea,0d0,FNPd,FNPu,FNPs,0d0,0d0/)

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

    muOPE=C0_const*c4/bT+2d0

    if(muOPE>1000d0) then
        muOPE=1000d0
    end if
end function muOPE
  
end module SiversTMDPDF_model
