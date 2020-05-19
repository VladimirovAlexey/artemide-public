!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!
!!!!!!!!!       part of TMDR module for artemide
!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!-------------------------------------------------------------------------------------------------------
!!!-----------------------------Routines for the evolution of type 2--------------------------------------
!!!----------------------------------(improved-gamma evolution)-------------------------------------------
!!!-------------------------------------------------------------------------------------------------------

!!! Evolution exponent in the improved gamma-picture
!!! Evolution from (muf,zetaf) -> (mui,zetai)
 function TMDR_R_type2(b,muf,zetaf,mui,zetai,f)
  real(dp)::TMDR_R_type2,b,muf,zetaf,mui,zetai
  integer::f
  
  counter=0
  
  if(b<1d-6) b=1d-6
  
  TMDR_R_type2=EXP(-IntegralG3(muf,mui,b,f)+DNP(muf,b,f)*Log(muf**2/zetaf)-DNP(mui,b,f)*Log(mui**2/zetai))
  
  !write(*,*) 'TMDR_R_type2: number of AD calls ',counter
  
  if(TMDR_R_type2>1d6) then
    write(*,*) ErrorString('Evolution factor (type2-1) is TOO HUGE check the formula',moduleName)
    write(*,*) 'NP parameters =', NPparam
    write(*,*) 'b=',b,'zetaf=',zetaf,'muf=',muf,'zetai=',zetai,'mui=',mui
    write(*,*) 'int=',IntegralG3(muf,mui,b,f),'t1=',DNP(muf,b,f)*Log(muf**2/zetaf),'t2=',DNP(mui,b,f)*Log(mui**2/zetai)
    write(*,*) 'Evaluation continue with R=10^6'
    TMDR_R_type2=1d6
  end if

 end function TMDR_R_type2

 !!! Evolution exponent in the improved gamma-picture to zeta-line
 !!! Evolution from (muf,zetaf) -> (zeta_mui)
 function TMDR_Rzeta_type2(b,muf,zetaf,mui,f)
  real(dp)::TMDR_Rzeta_type2,b,muf,zetaf,mui,zetaP
  integer::f
  
  counter=0
  
  if(b<1d-6) b=1d-6
  
  zetaP=zetaNP(mui,b,f)
  
  TMDR_Rzeta_type2=EXP(-IntegralG3(muf,mui,b,f)+DNP(muf,b,f)*Log(muf**2/zetaf)-DNP(mui,b,f)*Log(mui**2/zetaP))
  
  !write(*,*) 'TMDR_Rzeta_type2: number of AD calls ',counter
  
  if(TMDR_Rzeta_type2>1d6) then
    write(*,*) ErrorString('Evolution factor (type2-2) is TOO HUGE check the formula',moduleName)
    write(*,*) 'b=',b,'zetaf=',zetaf,'muf=',muf,'zetaP=',zetaP,'mui=',mui
    write(*,*) 'int=',IntegralG3(muf,mui,b,f),'t1=',DNP(muf,b,f)*Log(muf**2/zetaf),'t2=',DNP(mui,b,f)*Log(mui**2/zetaP)
    write(*,*) 'Evaluation continue with R=10^6'
    TMDR_Rzeta_type2=1d6
  end if

 end function TMDR_Rzeta_type2
 
 !!!This is the integral   \int_{mu_i}^{muf} dmu/mu (2 D(mu,b)+gammaV)
!!! evaluation by adaptive simpson
 function IntegralG3(muf,mui,b,f)
  real(dp)::IntegralG3,muf,mui,b
  integer::f
  
  real(dp)::y1,y3,y5,X1,X3,X5,maxV
  
  y1=Log(muf)
  y5=Log(mui)
  if(y1<y5) then
   y3=y1
   y1=y5
   y5=y3
  end if
  y3=(y1+y5)/2d0
  
  X1=2d0*DNP(Exp(y1),b,f)+gammaV(Exp(y1),f)
  X3=2d0*DNP(Exp(y3),b,f)+gammaV(Exp(y3),f)
  X5=2d0*DNP(Exp(y5),b,f)+gammaV(Exp(y5),f)
  
  counter=counter+3
  
  !maxV is approximate value of the integral, just to normalize the errors
  maxV=ABS((y1-y5)*(X1+4d0*X3+X5)/6d0)

  !!!! swithc sign if the integral is opposite
  if(muf>mui) then
  IntegralG3=integral3_S(y1,y5,b,f,X1,X3,X5,maxV)
  else
  IntegralG3=-integral3_S(y1,y5,b,f,X1,X3,X5,maxV)
  end if
  
 end function IntegralG3

!!! expected that y1>y2 (logarithm scale)
 recursive function integral3_S(y1,y5,b,f,X1,X3,X5,maxV) result(interX)
   real(dp) ::y1,y5,b
   integer::f
   real(dp) :: interX,X1,X2,X3,X4,X5
   real(dp) :: value,valueAB,valueACB
   real(dp) :: y2,y3,y4,deltay,maxV
   
   deltay=y1-y5
   y4=y5+deltay/4d0
   y3=y5+deltay/2d0
   y2=y1-deltay/4d0
   
   
   X2=2d0*DNP(Exp(y2),b,f)+gammaV(Exp(y2),f)
   X4=2d0*DNP(Exp(y4),b,f)+gammaV(Exp(y4),f)
   
   counter=counter+2
   
   valueAB=deltay*(X1+4d0*X3+X5)/6d0
   valueACB=deltay*(X1+4d0*X2+2d0*X3+4d0*X4+X5)/12d0
   
!    write(*,*) y1,y3,y5,valueACB-valueAB
   
   If(ABS(valueACB-valueAB)/15>tolerance*maxV) then
    interX=integral3_S(y1,y3,b,f,X1,X2,X3,maxV)&
	  +integral3_S(y3,y5,b,f,X3,X4,X5,maxV)
   else
    interX=valueACB
   end if
   
 end function integral3_S
  
  
