!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Model for Sivers TMD PDF [20??.????]
!
!				A.Vladimirov (21.05.2020)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module SiversTMDPDF_model
use aTMDe_Numerics
use IO_functions
use SnowFlake
!use SnowFlake_Model
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
real(dp)::FNPu,FNPd,FNPs,FNPc,FNPb,FNPuB,FNPdB,FNPsB,FNPcB,FNPbB

real(dp)::mu

!!! profile in b is common for all
bProfile=1._dp/cosh(lambdaNP(1)*bT)

mu=muOPE(bT,x,1.d0,1.d0)

FNPd=-pi*GetPDF(-x,0.0000000001d0,mu,1,outputT="T")
FNPu=-pi*GetPDF(-x,0.0000000001d0,mu,2,outputT="T")
FNPs=-pi*GetPDF(-x,0.0000000001d0,mu,3,outputT="T")
FNPc=-pi*GetPDF(-x,0.0000000001d0,mu,4,outputT="T")
FNPb=-pi*GetPDF(-x,0.0000000001d0,mu,5,outputT="T")

FNPdB=-pi*GetPDF(x,0.0000000001d0,mu,1,outputT="T")
FNPuB=-pi*GetPDF(x,0.0000000001d0,mu,2,outputT="T")
FNPsB=-pi*GetPDF(x,0.0000000001d0,mu,3,outputT="T")
FNPcB=-pi*GetPDF(x,0.0000000001d0,mu,4,outputT="T")
FNPbB=-pi*GetPDF(x,0.0000000001d0,mu,5,outputT="T")

FNP=bProfile*(/FNPbB,FNPcB,FNPsB,FNPuB,FNPdB,0d0,FNPd,FNPu,FNPs,FNPc,FNPb/)


if(ISNAN(FNP(5)) .or.ISNAN(FNP(4)) .or.ISNAN(FNP(3)) .or.ISNAN(FNP(2)) .or.ISNAN(FNP(1))) then
    write(*,*) "------ NAN INSIDE THE SIVERTMDPDF-model -----"
    write(*,*) "x,bT,hadron,lambdaNP--->",x,bT,hadron,lambdaNP
    write(*,*) "mu--->",mu
    write(*,*) "--->",FNP
    write(*,*) "bProfile--->",bProfile
    !call  PrintParameters()
end if

if(abs(FNP(1))>10d10 .or.abs(FNP(2))>10d10 .or.abs(FNP(3))>10d10) then
    write(*,*) "------ HUGE NUMBER INSIDE THE SIVERTMDPDF-model -----"
    write(*,*) "x,bT,hadron,lambdaNP--->",x,bT,hadron,lambdaNP
    write(*,*) "mu--->",mu
    write(*,*) "--->",FNP
    write(*,*) "bProfile--->",bProfile
    !call  PrintParameters()
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

!!!! ART25
bSTAR=bT*ee+(1-ee)*C0_const/muOPE(bT,x,y,1.d0)

!!!! ART23
!bSTAR=bT/sqrt(1d0+(bT/500d0)**2)
!bSTAR=bT/sqrt(1d0+(bT/1.d0)**2)

end function bSTAR

!!!!This function is the mu(x,b), which is used inside the OPE
!!!! x -- is the global x for TMDPDF,
!!!! y -- is the convolution variable in the definition \int dy/y C(y) PDF(x/y)
!!!! c4-- is the scale variation variable
pure function muOPE(bt,x,y,c4)
real(dp),intent(in)::bt,x,y,c4

muOPE=C0_const*c4/bT+5d0
if(muOPE>100d0) then
    muOPE=100d0
end if
end function muOPE
  
end module SiversTMDPDF_model
