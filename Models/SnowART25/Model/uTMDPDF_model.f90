!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Model function for unpolarized TMDPDF 
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module uTMDPDF_model
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
    
    NPparam=(/0.23d0, 0.01d0, 0.32d0, 9.3d0, 5.5d0, 1.8d0, 6.8d0, 0.1d0, 1.1d0, 3.8d0, 0.d0, 0.1d0/)

    write(*,*) color(">>>  The model for uTMDPDF for ART23   <<<",c_cyan)
    
end subroutine ModelInitialization

!!!!!! Write nessecery model update (e.g. save current NP-parameters)
!!!!!! newNPParams is the new NP-array
subroutine ModelUpdate(newNPParams)  
    real(dp),intent(in):: newNPParams(:)
    
    NPparam=newNPParams !! save new vector of NP-parameters

end subroutine ModelUpdate

!!! This is  non-perturbative function
!!! non=perturbative parameters are lambdaNP()
!!! x-- is the bjorken variable of TMD
function FNP(x,bT,hadron,lambdaNP)
  real(dp),intent(in)::x,bT    
  integer,intent(in)::hadron
  real(dp),intent(in)::lambdaNP(:)
  real*8::FNP0,FNPu,FNPd,FNPubar,FNPdbar,FNPr

   real*8::bb,w1,w2,wu,wd,wubar,wdbar,wr
   
   if(hadron==1) then
   
    !bb=bT**2
! ! !    ART23
!     wu=lambdaNP(1)*(1-x)+x*lambdaNP(2)
!     wd=lambdaNP(3)*(1-x)+x*lambdaNP(4)
!     wubar=lambdaNP(5)*(1-x)+x*lambdaNP(6)
!     wdbar=lambdaNP(7)*(1-x)+x*lambdaNP(8)
!     wr=lambdaNP(9)*(1-x)+x*lambdaNP(10)
!

!    ART25
   bb=bT**2
     wu=lambdaNP(1)*(1-x)**lambdaNP(5)+x*lambdaNP(2)
     wd=lambdaNP(3)*(1-x)**lambdaNP(7)+x*lambdaNP(4)
     wubar=lambdaNP(1)*(1-x)+x*lambdaNP(6)
     wdbar=lambdaNP(3)*(1-x)+x*lambdaNP(8)
     wr=lambdaNP(9)*(1-x)+x*lambdaNP(10)
    if(wu<0d0 .or. wd<0d0 .or. wubar<0d0 .or. wdbar<0d0 .or. wr<0d0) then
        FNPu=Exp(-10d0*bb)
        FNPd=Exp(-10d0*bb)
        FNPubar=Exp(-10d0*bb)
        FNPdbar=Exp(-10d0*bb)
        FNPr=Exp(-10d0*bb)
    else
        FNPu=1d0/cosh(wu*bT)
        FNPd=1d0/cosh(wd*bT)
        FNPubar=1d0/cosh(wubar*bT)
        FNPdbar=1d0/cosh(wdbar*bT)
        FNPr=1d0/cosh(wr*bT)
    end if
!
!     FNP0=1/cosh((lambdaNP(1)*x+(1-x)*lambdaNP(2))*bT)
!     bb=bT**2/(lambdaNP(3)**2+bT**2)
!
!     FNPu=   FNP0*(1d0+(lambdaNP(4)*x+lambdaNP(5)*(1-x)-lambdaNP(6)*log(x))*bb)
!     FNPd=   FNP0*(1d0+(lambdaNP(7)*x+lambdaNP(8)*(1-x)-lambdaNP(9)*log(x))*bb)
!     FNPubar=FNP0*(1d0+(lambdaNP(10)*x+lambdaNP(5)*(1-x)-lambdaNP(6)*log(x))*bb)
!     FNPdbar=FNP0*(1d0+(lambdaNP(11)*x+lambdaNP(8)*(1-x)-lambdaNP(9)*log(x))*bb)
!     FNPr=   FNP0*(1d0+(lambdaNP(12)*x+lambdaNP(13)*(1-x)-lambdaNP(14)*log(x))*bb)


    FNP=(/&
    FNPr,FNPr,FNPr,FNPubar,FNPdbar,&
    exp(-0.25d0*bb),&
    FNPd,FNPu,FNPr,FNPr,FNPr/)

  else 
      bb=bT**2
      w1=(lambdaNP(7)+(1-x)**2*lambdaNP(8))
      w2=lambdaNP(9)
      if(w2<0d0 .or. w1<0d0) then
      FNP0=-1d0
      else
      FNP0=Exp(-w1*bb/sqrt(1+w2*bb))
      end if

      FNP=FNP0*(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
      
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
    !muOPE=C0_const*c4/bT+2d0

    !!!! like traditional b*
    !muOPE=C0_const*sqrt(1+bT**2/(c4*C0_const)**2)*c4/bT
    !!!! like MAP24
    !muOPE=(1-Exp(-(c4*bT/C0_const)**4))**(-0.25)
    
    if(muOPE>1000d0) then
        muOPE=1000d0
    end if
end function muOPE

end module uTMDPDF_model
