!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!
!!!!!!!!!       part of TMDR module for artemide
!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!-------------------------------------------------------------------------------------------------------
!!!-----------------------------Routines for the evolution of type 1--------------------------------------
!!!--------------------------(ordinary CS evolution=improved-D evolution)---------------------------------
!!!-------------------------------------------------------------------------------------------------------

!!! Evolution exponent in the improved D-picture
!!! from point to point, via intemidiate scale mu0
function TMDR_R_type1(b,muf,zetaf,mui,zetai,mu0,f)
  real(dp)::TMDR_R_type1,b,muf,zetaf,mui,zetai,mu0
  integer::f
  
  counter=0
  
  if(b<1d-6) b=1d-6
  
  TMDR_R_type1=EXP(IntegralG1(muf,mui,zetaf,f)-(IntegralG2(mui,mu0,f)+DNP(mu0,b,f))*Log(zetaf/zetai))
  
  !write(*,*) 'TMDR_R_type1: number of AD calls ',counter

end function TMDR_R_type1


!!! Evolution exponent in the improved D-picture to zeta-line
!!! The evolution from point (muf,zetaf) -> (zeta_mui) via scale mu0
 function TMDR_Rzeta_type1(b,muf,zetaf,mui,mu0,f)
  real(dp)::TMDR_Rzeta_type1,b,muf,zetaf,mui,mu0,zetaP
  integer::f
  
  counter=0
  
  if(b<1d-6) b=1d-6
  
  zetaP=zetaNP(mui,b,f)
  
  TMDR_Rzeta_type1=EXP(IntegralG1(muf,mui,zetaf,f)-(IntegralG2(mui,mu0,f)+DNP(mu0,b,f))*Log(zetaf/zetaP))
  
  !write(*,*) 'TMDR_Rzeta_type1: number of AD calls ',counter
  
 end function TMDR_Rzeta_type1

 !!!This is the integral   \int_{mu_i}^{muf} dmu/mu (Gamma(mu)Log[mu^2/zeta]-gammaV)
!!! evaluation by adaptive simpson
 function IntegralG1(muf,mui,zeta,f)
  real(dp)::IntegralG1,muf,mui,zeta,lnz
  integer::f
  
  real(dp)::y1,y3,y5,X1,X3,X5,maxV
  
  lnz=Log(zeta)
  
  y1=Log(muf)
  y5=Log(mui)
  if(y1<y5) then
   y3=y1
   y1=y5
   y5=y3
  end if
  y3=(y1+y5)/2d0
  
  X1=gammaCUSP(Exp(y1),f)*(2d0*y1-lnz)-gammaV(Exp(y1),f)
  X3=gammaCUSP(Exp(y3),f)*(2d0*y3-lnz)-gammaV(Exp(y3),f)
  X5=gammaCUSP(Exp(y5),f)*(2d0*y5-lnz)-gammaV(Exp(y5),f)
  
  counter=counter+3
  
  !maxV is approximate value of the integral, just to normalize the errors
  maxV=ABS((y1-y5)*(X1+4d0*X3+X5)/6d0)

  !!!! swithc sign if the integral is opposite
  if(muf>mui) then
  IntegralG1=integral1_S(y1,y5,lnz,f,X1,X3,X5,maxV)
  else
  IntegralG1=-integral1_S(y1,y5,lnz,f,X1,X3,X5,maxV)
  end if
  
 end function IntegralG1

!!! expected that y1>y2 (logarithm scale)
 recursive function integral1_S(y1,y5,lnz,f,X1,X3,X5,maxV) result(interX)
   real(dp) ::y1,y5,lnz
   integer::f
   real(dp) :: interX,X1,X2,X3,X4,X5
   real(dp) :: value,valueAB,valueACB
   real(dp) :: y2,y3,y4,deltay,maxV
   
   deltay=y1-y5
   y4=y5+deltay/4d0
   y3=y5+deltay/2d0
   y2=y1-deltay/4d0
   
   
   X2=gammaCUSP(Exp(y2),f)*(2d0*y2-lnz)-gammaV(Exp(y2),f)
   X4=gammaCUSP(Exp(y4),f)*(2d0*y4-lnz)-gammaV(Exp(y4),f)
   
   counter=counter+2
   
   valueAB=deltay*(X1+4d0*X3+X5)/6d0
   valueACB=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
!    write(*,*) y1,y3,y5,valueACB-valueAB
   
   If(ABS(valueACB-valueAB)/15d0>tolerance*maxV) then
    interX=integral1_S(y1,y3,lnz,f,X1,X2,X3,maxV)&
	  +integral1_S(y3,y5,lnz,f,X3,X4,X5,maxV)
   else
    interX=valueACB
   end if
   
 end function integral1_S


 !!!This is the integral   \int_{mu_i}^{muf} dmu/mu Gamma(mu)
!!! evaluation by adaptive simpson
function IntegralG2(muf,mui,f)
  real(dp)::IntegralG2,muf,mui
  integer::f
  
  real(dp)::y1,y3,y5,X1,X3,X5,maxV
  
  if(Abs(muf-mui)<tolerance) then
    IntegralG2=0d0
    return
  end if
  
  y1=Log(muf)
  y5=Log(mui)
  if(y1<y5) then
   y3=y1
   y1=y5
   y5=y3
  end if
  y3=(y1+y5)/2d0
  
  X1=gammaCUSP(Exp(y1),f)
  X3=gammaCUSP(Exp(y3),f)
  X5=gammaCUSP(Exp(y5),f)
  
  counter=counter+3
  
  !maxV is approximate value of the integral, just to normalize the errors
  maxV=ABS((y1-y5)*(X1+4d0*X3+X5)/6d0)

  !!!! swithc sign if the integral is opposite
  if(muf>mui) then
  IntegralG2=integral2_S(y1,y5,f,X1,X3,X5,maxV)
  else
  IntegralG2=-integral2_S(y1,y5,f,X1,X3,X5,maxV)
  end if
  
 end function IntegralG2

!!! expected that y1>y2 (logarithm scale)
 recursive function integral2_S(y1,y5,f,X1,X3,X5,maxV) result(interX)
   real(dp) ::y1,y5
   integer::f
   real(dp) :: interX,X1,X2,X3,X4,X5
   real(dp) :: value,valueAB,valueACB
   real(dp) :: y2,y3,y4,deltay,maxV
   
   deltay=y1-y5
   y4=y5+deltay/4d0
   y3=y5+deltay/2d0
   y2=y1-deltay/4d0
   
   
   X2=gammaCUSP(Exp(y2),f)
   X4=gammaCUSP(Exp(y4),f)
   
   counter=counter+2
   
   valueAB=deltay*(X1+4d0*X3+X5)/6d0
   valueACB=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
!    write(*,*) y1,y3,y5,valueACB-valueAB
   
   If(ABS(valueACB-valueAB)/15d0>tolerance*maxV) then
    interX=integral2_S(y1,y3,f,X1,X2,X3,maxV)&
	  +integral2_S(y3,y5,f,X3,X4,X5,maxV)
   else
    interX=valueACB
   end if
   
 end function integral2_S
