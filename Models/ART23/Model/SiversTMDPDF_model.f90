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
!!!!!    arg=array of initial NP-parameters
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
!!!!!    arg=(b,lambdaNP) with b=transverse distance(real_dp), lambdaNP = array of NP parameters (real_dp(:))
real(dp),public:: bSTAR
!!!!! 5) Function which returns the scale of matching (OPE scale)
!!!!!    arg=(bt) with b=transverse distance(real_dp)
real(dp),public:: mu_OPE
!!!!! 6) Subroutine which returns the array of parameters CA which compose the TMDs into a single one
!!!!!    i.e. the TMD for hardon=h is build as TMD(h)=Sum_c CA(h,c) TMD(c)
!!!!!    it is used only if the option UseComposite TMD is ON,
!!!!!    arg=(h,lambdaNP,includeArray,CA) with h=hadron(integer),lambdaNP = array of NP parameters (real_dp(:))
!!!!!    includeArray=logical array with .true. for terms included in the sum (logical(:),allocatable,intent(out))
!!!!!    CA=coefficient CA (real_dp(:),allocatable,intent(out))
public:: GetCompositionArray
!!!!! 7) Subroutine which returns the array of NP-parameters corresponding to certain integer (replica)
!!!!!    arg=rep input integer,  NParray (real_dp(:), allocatable, intent(out))  returned array
public:: GetReplicaParameters
!!!!!------------------------------------------------------------------------------------

real(dp),allocatable::NPparam(:)

contains  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! USER DEFINED FUNCTIONS   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!! Write nessecery model intitialization.
subroutine ModelInitialization(NPstart)
    real(dp),intent(in)::NPstart(:)
    allocate(NPparam(1:size(NPstart)))
    NPparam=NPstart
    
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
  
!!!! This is the function b* that enter the logarithms of coefficient function
!!!! at small-b it should be ~b to match the collinear regime
!!!! at large-b it is a part of model
!!!! NOTE: if it is lambda-dependent, the grid will be recalculate each reset of lambdaNP
pure function bSTAR(bT,lambdaNP)
    real(dp),intent(in)::bT
    real(dp),intent(in)::lambdaNP(:)

    bSTAR=bT/sqrt(1d0+(bT/500d0)**2)

end function bSTAR
  
!!!!This function is the mu(x,b), which is used inside the OPE
pure function mu_OPE(bt)
    real(dp),intent(in)::bt

    mu_OPE=C0_const*1d0/bT+2d0

    if(mu_OPE>1000d0) then
        mu_OPE=1000d0
    end if
end function mu_OPE
  
!!!! if the option UseComposite TMD is OFF, this function is ignored
!!!! If the option UseComposite TMD is ON,
!!!! than the TMD for hardon is build as TMD(hadron)=Sum_c CA(h,c) TMD(c)
!!!! where h=hadron, CA=coefficientArray
!!!! coefficientArray real(dp) list of coefficeints
!!!! includeArray is logical array list (true=TMD(c) is computed, false TMD(c) ignored)
subroutine GetCompositionArray(hadron,lambdaNP,includeArray,coefficientArray)  
    real(dp),intent(in)::lambdaNP(:)
    integer::hadron
    logical,allocatable,intent(out)::includeArray(:)
    real(dp),allocatable,intent(out)::coefficientArray(:)

    allocate(includeArray(1:1))
    allocate(coefficientArray(1:1))
end subroutine GetCompositionArray
  
  
!!! In SV19 model the replica parameters are stored in separate file.
subroutine GetReplicaParameters(rep,NParray)
    integer,intent(in)::rep
    real(dp),allocatable,intent(out)::NParray(:)
    integer::i
    
    allocate(NParray(1:size(NPparam)))

    write(*,*) warningstring("set model replica via artemide-control module","Sivers")
    write(*,*) warningstring("some generic NP values returned","Sivers")
    NParray(1)=1d0
    do i=2,size(NPparam)
        NParray(1)=0.001d0
    end do

end subroutine GetReplicaParameters
  
end module SiversTMDPDF_model
