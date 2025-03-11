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

    real(dp)::FNP0,FNP1,FNP2,FNP3
    real(dp)::bb,w1,w2,Fmain

    bb=bT**2/x**2
    !!ART25
    if(hadron==1 .or. hadron==3) then
        !!! pion
        Fmain=1._dp/cosh(lambdaNP(1)*bT/x)
        !!! u
        FNP1=(1+lambdaNP(2)*bb)*Fmain
        !!! dbar
        FNP2=(1+lambdaNP(3)*bb)*Fmain
        !!! sea
        FNP0=(1+lambdaNP(4)*bb)*Fmain
        !!! uBar
        FNP3=(1+lambdaNP(9)*bb)*Fmain

        if(hadron==1) FNP=(/FNP0,FNP0,FNP0,FNP3,FNP2, Fmain  ,FNP0,FNP1,FNP0,FNP0,FNP0/) !!!! pion+ [u dBar]
        if(hadron==3) FNP=(/FNP0,FNP0,FNP0,FNP1,FNP0, Fmain  ,FNP2,FNP3,FNP0,FNP0,FNP0/) !!!! pion- [d uBar]
    else if(hadron==2 .or. hadron==4) then
        !!! kaon
        Fmain=1._dp/cosh(lambdaNP(5)*bT/x)
        !!! u
        FNP1=(1+lambdaNP(6)*bb)*Fmain
        !!! sbar
        FNP2=(1+lambdaNP(7)*bb)*Fmain
        !!! sea
        FNP0=(1+lambdaNP(8)*bb)*Fmain
        !!! uBar
        FNP3=(1+lambdaNP(10)*bb)*Fmain

        if(hadron==2) FNP=(/FNP0,FNP0,FNP2,FNP3,FNP0, Fmain ,FNP0,FNP1,FNP0,FNP0,FNP0/) !!!! kaon+ [u sBar]
        if(hadron==4) FNP=(/FNP0,FNP0,FNP0,FNP1,FNP0, Fmain ,FNP0,FNP3,FNP2,FNP0,FNP0/) !!!! kaon- [s uBar]
    else
        FNP=1._dp/cosh(lambdaNP(1)*bT/x)*(/1.d0,1.d0,1.d0,1.d0,1.d0,0.d0,1.d0,1.d0,1.d0,1.d0,1.d0/)
    end if


end function FNP
  
!!!! This is the function b* that enters the logarithms of coefficient function
!!!! at small-b it should be ~b to match the collinear regime
!!!! at large-b it is a part of model
!!!! x -- is the global x for TMDPDF,
!!!! y -- is the convolution variable in the definition \int dy/y C(y) PDF(x/y)
pure function bSTAR(bT,x,y)
    real(dp),intent(in)::bT,x,y
    real(dp)::ee

    ee=exp(-0.04d0*bT**2)

    !bSTAR=bT/sqrt(1d0+(bT/500d0)**2)
    !bSTAR=bT/sqrt(1d0+(bT/1.d0)**2)
    bSTAR=bT*ee+(1-ee)*C0_const/muOPE(bT,x,y,1.d0)

end function bSTAR

!!!!This function is the mu(x,b), which is used inside the OPE
!!!! x -- is the global x for TMDPDF,
!!!! y -- is the convolution variable in the definition \int dy/y C(y) PDF(x/y)
!!!! c4-- is the scale variation variable
pure function muOPE(bt,x,y,c4)
    real(dp),intent(in)::bt,x,y,c4

    muOPE=C0_const*c4*x/bT+5d0

    if(muOPE>100d0) then
        muOPE=99d0
    end if
end function muOPE

end module uTMDFF_model
