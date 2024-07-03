!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 2.04
!
!	The module that contains standard integration routines (mainly used for bin integration)
!       the main difference from library routines is the absence of checks (boundaries, etc.), 
!           the checks are done in the artemide main code.
!
!   Currently it incorporates:
!       S5: 5-point Simpsons
!       SN: N-point Simpsons (N=even)
!       SA: Adaptive Simpsons
!       G7: 7-point Gauss
!       K15: 15-point Kronrod
!       GK: Adaptive Gauss Kronrod 7/15
!
!				A.Vladimirov (17.04.2020)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module IntegrationRoutines
use aTMDe_numerics
implicit none

private
INCLUDE '../Tables/G7K15.f90'
INCLUDE '../Tables/G20K41.f90'

public::Integrate_S5,Integrate_SN,Integrate_SA,Integrate_SA_2D
public::Integrate_G3,Integrate_G7,Integrate_K15,Integrate_K41,Integrate_GK,Integrate_GK2041
public::Integrate_GK_array5


real(dp),dimension(1:3,1:3),parameter::SA3_2D_w=reshape((/1._dp,4._dp,1._dp,&
                                                  4._dp,16._dp,4._dp,&
                                                  1._dp,4._dp,1._dp/)/36._dp,shape(SA3_2D_w))

real(dp),dimension(1:5,1:5),parameter::SA5_2D_w=reshape((/1._dp,4._dp,2._dp,4._dp,1._dp,&
                                                  4._dp,16._dp,8._dp,16._dp,4._dp,&
                                                  2._dp,8._dp,4._dp,8._dp,2._dp,&
                                                  4._dp,16._dp,8._dp,16._dp,4._dp,&
                                                  1._dp,4._dp,2._dp,4._dp,1._dp/)/144._dp,shape(SA5_2D_w))

!!! parameters for Gauss quadrature with 3 points
real(dp), parameter, dimension(1:3) :: Xi_g3 = (/-sqrt(0.6_dp),0._dp,sqrt(0.6_dp)/)
real(dp), parameter, dimension(1:3) :: Wi_g3 = (/5._dp/9._dp,8._dp/9._dp,5._dp/9._dp/)

!!! this is interface for function (-5:5) in the integration
abstract interface 
    function func_array5(x)
        import::dp
        real(dp),dimension(-5:5) :: func_array5
        real(dp), intent(in) ::x
    end function func_array5
end interface
  
!!! This is parameter to which the integral is compared.
!!! if absolute value of integral <zero. Its returned (otherwise, there could be infinite loop of precision)
real(dp),parameter:: zero=10.d-8

contains


!------------------------------------------- S5 --------------------------------------------

!!! Simpson by 5 points
!!! Use it for estimations only!
!!! f::  function of 1 variable
!!! xMin, and xMax boundaries of the integral. xMax>xMin !!
function Integrate_S5(f,xMin,xMax)
    real(dp)::f,Integrate_S5
    real(dp),intent(in)::xMin,xMax
    real(dp)::delta,x2,x3,x4,f1,f2,f3,f4,f5

    delta=xMax-xMin
    x2=xMin+delta/4._dp
    x3=xMin+delta/2._dp
    x4=xMax-delta/4._dp

    f1=f(xMin)
    f2=f(x2)
    f3=f(x3)
    f4=f(x4)
    f5=f(xMax)
        
    Integrate_S5=delta*(f1+4._dp*f2+2._dp*f3+4._dp*f4+f5)/12._dp
    
end function Integrate_S5

!------------------------------------------- SN --------------------------------------------

!!! Simpson by N points
!!! f::  function of 1 variable
!!! xMin, and xMax boundaries of the integral. xMax>xMin !!
!!! N is an even integer N>2
function Integrate_SN(f,xMin,xMax,N)
    real(dp)::f,Integrate_SN
    real(dp),intent(in)::xMin,xMax
    real(dp)::delta,inter,xCur
    integer::N,i

    !!!!!!!!!!!!!!!!!!!fixed number Simpsons
    delta=(xMax-xMin)/N
    inter=f(xMin)   !!!! first term
    
    !!!! even terms
    do i=1,N-1,2
        xCur=xMin+i*delta
        inter=inter+4._dp*f(xCur)    
    end do
    !!!! odd term
    do i=2,N-2,2
        xCur=xMin+i*delta
        inter=inter+2._dp*f(xCur)    
    end do
    
    inter=inter+f(xMax)!!!! last term
        
    Integrate_SN=delta/3._dp*inter
    
end function Integrate_SN


!------------------------------------------- SA --------------------------------------------

!!! Simpson adaptive
!!! f::  function of 1 variable
!!! xMin, and xMax boundaries of the integral. xMax>xMin !!
!!! tolerance is relative (it is weighted by approximate value of integral)
function Integrate_SA(f,xMin,xMax,tolerance)
    real(dp)::f,Integrate_SA,tolerance,eps
    real(dp),intent(in)::xMin,xMax
    real(dp)::delta,x2,x3,x4,f1,f2,f3,f4,f5

    
    delta=xMax-xMin
    x2=xMin+delta/4._dp
    x3=xMin+delta/2._dp
    x4=xMax-delta/4._dp

    f1=f(xMin)
    f2=f(x2)
    f3=f(x3)
    f4=f(x4)
    f5=f(xMax)
    
    !!! the error parameter is weighted with the approximate integral size
    eps=tolerance*abs(delta*(f1+4._dp*f2+2._dp*f3+4._dp*f4+f5)/12._dp)
        
    Integrate_SA=SA_Rec(f,xMin,x2,x3,f1,f2,f3,eps)+SA_Rec(f,x3,x4,xMax,f3,f4,f5,eps)
    
end function Integrate_SA

recursive function SA_Rec(f,x1,x3,x5,f1,f3,f5,eps) result(res)
    real(dp)::f,x1,x2,x3,x4,x5
    real(dp)::f1,f2,f3,f4,f5,eps,res
    real(dp)::value15,value135
    
    x2=(x1+x3)/2._dp
    f2=f(x2)
    x4=(x3+x5)/2._dp
    f4=f(x4)
    
    value15=(x5-x1)*(f1+4._dp*f3+f5)/6._dp
    value135=(x5-x1)*(f1+4._dp*f2+2._dp*f3+4._dp*f4+f5)/12._dp
        
    If(ABS(value135-value15)>eps) then
        res=SA_Rec(f,x1,x2,x3,f1,f2,f3,eps)+SA_Rec(f,x3,x4,x5,f3,f4,f5,eps)
    else
        res=value135
    end if
    
end function SA_Rec

!------------------------------------------- SA 2D--------------------------------------------

!!! Simpson adaptive
!!! f::  function of 1 variable
!!! xMin, and xMax boundaries of the integral. xMax>xMin !!
!!! tolerance is relative (it is weighted by approximate value of integral)
function Integrate_SA_2D(f,xMin,xMax,yMin,yMax,tolerance)
    real(dp)::f,Integrate_SA_2D,tolerance,eps
    real(dp),intent(in)::xMin,xMax,yMin,yMax
    real(dp)::deltaX,deltaY,x2,x3,x4,y2,y3,y4,rInt
    real(dp),dimension(1:5,1:5)::ff

    !write(*,*) "INTEGRATE_SA_2D: this routine is not tested"

    deltaX=(xMax-xMin)/4._dp
    x2=xMin+deltaX
    x3=xMin+deltaX*2._dp
    x4=xMax-deltaX

    deltaY=(yMax-yMin)/4._dp
    y2=yMin+deltaY
    y3=yMin+deltaY*2._dp
    y4=yMax-deltaY

    ff=reshape((/&
    f(xMin,yMin),f(xMin,y2),f(xMin,y3),f(xMin,y4),f(xMin,yMax),&
    f(x2,yMin),f(x2,y2),f(x2,y3),f(x2,y4),f(x2,yMax),&
    f(x3,yMin),f(x3,y2),f(x3,y3),f(x3,y4),f(x3,yMax),&
    f(x4,yMin),f(x4,y2),f(x4,y3),f(x4,y4),f(x4,yMax),&
    f(xMax,yMin),f(xMax,y2),f(xMax,y3),f(xMax,y4),f(xMax,yMax)/),shape(SA5_2D_w))

    rInt=sum(ff*SA5_2D_w)

    !!! the error parameter is weighted with the approximate integral size
    eps=tolerance*abs(deltaX*deltaY*rInt)

    Integrate_SA_2D=SA2D_Rec(f,xMin,x2,x3,yMin,y2,y3,ff(1:3,1:3),eps)&
                    +SA2D_Rec(f,x3,x4,xMax,yMin,y2,y3,ff(3:5,1:3),eps)&
                    +SA2D_Rec(f,xMin,x2,x3,y3,y4,yMax,ff(1:3,3:5),eps)&
                    +SA2D_Rec(f,x3,x4,xMax,y3,y4,yMax,ff(3:5,3:5),eps)

end function Integrate_SA_2D

recursive function SA2D_Rec(f,x1,x3,x5,y1,y3,y5,fIn,eps) result(res)
    real(dp),dimension(1:3,1:3)::fIn
    real(dp)::f,x1,x2,x3,x4,x5,y1,y2,y3,y4,y5
    real(dp)::eps,res
    real(dp)::value3,value5
    real(dp),dimension(1:5,1:5)::ff

    x2=(x1+x3)/2._dp
    x4=(x3+x5)/2._dp
    y2=(y1+y3)/2._dp
    y4=(y3+y5)/2._dp

    ff=reshape((/&
    fIn(1,1),f(x1,y2),fIn(1,2),f(x1,y4),fIn(1,3),&
    f(x2,y1),f(x2,y2),f(x2,y3),f(x2,y4),f(x2,y5),&
    fIn(2,1),f(x3,y2),fIn(2,2),f(x3,y4),fIn(2,3),&
    f(x4,y1),f(x4,y2),f(x4,y3),f(x4,y4),f(x4,y5),&
    fIn(3,1),f(x5,y2),fIn(3,2),f(x5,y4),fIn(3,3)/),shape(SA5_2D_w))

    value3=(x5-x1)*(y5-y1)*sum(fIn*SA3_2D_w)
    value5=(x5-x1)*(y5-y1)*sum(ff*SA5_2D_w)

    If(ABS(value5-value3)>eps) then
        res=SA2D_Rec(f,x1,x2,x3,y1,y2,y3,ff(1:3,1:3),eps)&
                    +SA2D_Rec(f,x3,x4,x5,y1,y2,y3,ff(3:5,1:3),eps)&
                    +SA2D_Rec(f,x1,x2,x3,y3,y4,y5,ff(1:3,3:5),eps)&
                    +SA2D_Rec(f,x3,x4,x5,y3,y4,y5,ff(3:5,3:5),eps)
    else
        res=value5
    end if

end function SA2D_Rec

!------------------------------------------- G3 --------------------------------------------

!!! Gauss 3-points (very inaccurate only for estimations)
!!! f::  function of 1 variable
!!! xMin, and xMax boundaries of the integral. xMax>xMin !!
function Integrate_G3(f,xMin,xMax)
    real(dp)::f,Integrate_G3
    real(dp),intent(in)::xMin,xMax
    real(dp)::delta,av

    delta=(xMax-xMin)/2._dp
    av=(xMax+xMin)/2._dp

    Integrate_G3=delta*(Wi_g3(1)*f(Xi_g3(1)*delta+av)&
                        +Wi_g3(2)*f(Xi_g3(2)*delta+av)&
                        +Wi_g3(3)*f(Xi_g3(3)*delta+av))

end function Integrate_G3

!------------------------------------------- G7 --------------------------------------------

!!! Gauss 7-points
!!! f::  function of 1 variable
!!! xMin, and xMax boundaries of the integral. xMax>xMin !!
function Integrate_G7(f,xMin,xMax)
    real(dp)::f,Integrate_G7
    real(dp),intent(in)::xMin,xMax
    real(dp)::delta,av,inter
    integer::i

    delta=(xMax-xMin)/2._dp
    av=(xMax+xMin)/2._dp
    
    inter=0._dp
    do i=1,7
        inter=inter+Wi_g77(i)*f(Xi_g7(i)*delta+av)
    end do
        
    Integrate_G7=delta*inter
    
end function Integrate_G7

!------------------------------------------- K15 --------------------------------------------

!!! Kronrod 15-points
!!! f::  function of 1 variable
!!! xMin, and xMax boundaries of the integral. xMax>xMin !!
function Integrate_K15(f,xMin,xMax)
    real(dp)::f,Integrate_K15
    real(dp),intent(in)::xMin,xMax
    real(dp)::delta,av,inter
    integer::i

    delta=(xMax-xMin)/2._dp
    av=(xMax+xMin)/2._dp
    
    inter=0._dp
    do i=1,15
        inter=inter+Wi_k15(i)*f(Xi_k15(i)*delta+av)
    end do
        
    Integrate_K15=delta*inter
    
end function Integrate_K15

!------------------------------------------- K41 --------------------------------------------

!!! Kronrod 41-points
!!! f::  function of 1 variable
!!! xMin, and xMax boundaries of the integral. xMax>xMin !!
function Integrate_K41(f,xMin,xMax)
    real(dp)::f,Integrate_K41
    real(dp),intent(in)::xMin,xMax
    real(dp)::delta,av,inter
    integer::i

    delta=(xMax-xMin)/2._dp
    av=(xMax+xMin)/2._dp

    inter=0._dp
    do i=1,41
        inter=inter+Wi_k41(i)*f(Xi_k41(i)*delta+av)
    end do

    Integrate_K41=delta*inter

end function Integrate_K41

!------------------------------------------- GK7/15 --------------------------------------------

!!! Gauss-Kronrod 7/15 adaptive
!!! f::  function of 1 variable
!!! xMin, and xMax boundaries of the integral. xMax>xMin !!
!!! tolerance is relative (it is wieghted by approximate value of integral)
function Integrate_GK(f,xMin,xMax,tolerance)
    real(dp)::f,Integrate_GK
    real(dp),intent(in)::xMin,xMax
    real(dp)::delta,av,g7,k15,eps,tolerance,fI
    integer::i

    delta=(xMax-xMin)/2._dp
    av=(xMax+xMin)/2._dp
    
    g7=0._dp
    k15=0._dp
    do i=1,15
        fI=f(Xi_k15(i)*delta+av)
        g7=g7+Wi_g7(i)*fI
        k15=k15+Wi_k15(i)*fI
    end do
    
    eps=delta*abs(k15)*tolerance

    if(abs(eps)<zero) eps=zero
    
    if(delta*abs(k15-g7)>eps) then
        Integrate_GK=GK_Rec(f,xMin,av,eps)+GK_Rec(f,av,xMax,eps)
    else
        Integrate_GK=delta*k15
    end if
    
end function Integrate_GK

recursive function GK_Rec(f,xMin,xMax,eps) result(res)
    real(dp)::f,res
    real(dp),intent(in)::xMin,xMax
    real(dp)::delta,av,g7,k15,eps,fI
    integer::i

    delta=(xMax-xMin)/2._dp
    av=(xMax+xMin)/2._dp
    
    g7=0._dp
    k15=0._dp
    do i=1,15
        fI=f(Xi_k15(i)*delta+av)
        g7=g7+Wi_g7(i)*fI
        k15=k15+Wi_k15(i)*fI
    end do

    if(delta*abs(k15-g7)>eps) then
        res=GK_Rec(f,xMin,av,eps)+GK_Rec(f,av,xMax,eps)
    else
        res=delta*k15
    end if
    
end function GK_Rec

!------------------------------------------- GK7/15 (-5:5) --------------------------------------------

!!! Gauss-Kronrod 7/15 adaptive for array
!!! f::  array(-5:5) function of 1 variable
!!! xMin, and xMax boundaries of the integral. xMax>xMin !!
!!! tolerance is relative (it is wieghted by approximate value of integral)
!!! the integral is considered convergent once all components are convergent
function Integrate_GK_array5(f,xMin,xMax,tolerance)
    real(dp),dimension(-5:5)::Integrate_GK_array5
    procedure(func_array5)::f
    real(dp),intent(in)::xMin,xMax,tolerance
    real(dp)::delta,av
    real(dp), dimension(-5:5)::fI,g7,k15,eps
    integer::i
    logical::ISconvergent
    
    
    
    delta=(xMax-xMin)/2._dp
    av=(xMax+xMin)/2._dp
    g7=0._dp
    k15=0._dp
    do i=1,15
        fI=f(Xi_k15(i)*delta+av)
        g7=g7+Wi_g7(i)*fI
        k15=k15+Wi_k15(i)*fI
    end do
    
    !!! check convergence
    eps=delta*abs(k15)*tolerance
    !!! if integral is almost zero, I compare to "zero"
    do i=-5,5
        if(abs(eps(i))<zero) eps(i)=zero
    end do
    
    ISconvergent=.true.
    do i=-5,5
       if(delta*abs(k15(i)-g7(i))>eps(i)) then
        ISconvergent=.false.
        exit
       end if
    end do
    
    if(ISconvergent) then
        Integrate_GK_array5=delta*k15                
    else
        Integrate_GK_array5=GK_array5_Rec(f,xMin,av,eps)+GK_array5_Rec(f,av,xMax,eps)
    end if
    
end function Integrate_GK_array5

recursive function GK_array5_Rec(f,xMin,xMax,eps) result(res)
    real(dp),dimension(-5:5)::res
    procedure(func_array5)::f
    real(dp),intent(in)::xMin,xMax
    real(dp),dimension(-5:5),intent(in)::eps
    real(dp)::delta,av
    real(dp),dimension(-5:5)::g7,k15,fI
    integer::i
    logical::ISconvergent
    
    delta=(xMax-xMin)/2._dp
    av=(xMax+xMin)/2._dp

    g7=0._dp
    k15=0._dp
    do i=1,15
        fI=f(Xi_k15(i)*delta+av)
        g7=g7+Wi_g7(i)*fI
        k15=k15+Wi_k15(i)*fI
    end do
    
    !!! check convergence
    ISconvergent=.true.
    do i=-5,5
       if(delta*abs(k15(i)-g7(i))>eps(i)) then
        ISconvergent=.false.
        exit
       end if
    end do
    
    if(ISconvergent) then
        res=delta*k15                
    else
        res=GK_array5_Rec(f,xMin,av,eps)+GK_array5_Rec(f,av,xMax,eps)
    end if
    
end function GK_array5_Rec



!------------------------------------------- GK20/41 adaptive --------------------------------------------
!!! Gauss-Kronrod 20/41 adaptive
!!! f::  function of 1 variable
!!! xMin, and xMax boundaries of the integral. xMax>xMin !!
!!! tolerance is relative (it is wieghted by approximate value of integral)
function Integrate_GK2041(f,xMin,xMax,tolerance)
    real(dp)::f,Integrate_GK2041
    real(dp),intent(in)::xMin,xMax,tolerance
    real(dp)::delta,av,g7,k15,eps,fI
    integer::i

    delta=(xMax-xMin)/2._dp
    av=(xMax+xMin)/2._dp

    if(delta<zero) then
        Integrate_GK2041=0._dp
        return
    end if

    g7=0._dp
    k15=0._dp
    do i=1,41
        fI=f(Xi_k41(i)*delta+av)
        g7=g7+Wi_g20(i)*fI
        k15=k15+Wi_k41(i)*fI
    end do

    eps=delta*abs(k15)*tolerance

!    write(*,*) "-->",delta*abs(k15-g7),eps,delta*k15
    if(abs(delta*k15)<zero) then
        Integrate_GK2041=delta*k15
    else if(delta*abs(k15-g7)>eps) then
        Integrate_GK2041=GK2041_Rec(f,xMin,av,eps)+GK2041_Rec(f,av,xMax,eps)
    else
        Integrate_GK2041=delta*k15
    end if

end function Integrate_GK2041

recursive function GK2041_Rec(f,xMin,xMax,eps) result(res)
    real(dp)::f,res
    real(dp),intent(in)::xMin,xMax,eps
    real(dp)::delta,av,g7,k15,fI
    integer::i

    delta=(xMax-xMin)/2._dp
    av=(xMax+xMin)/2._dp

    if(delta<zero) then
        res=0._dp
        return
    end if

    g7=0._dp
    k15=0._dp
    do i=1,41
        fI=f(Xi_k41(i)*delta+av)
        g7=g7+Wi_g20(i)*fI
        k15=k15+Wi_k41(i)*fI
    end do

    if(abs(delta*k15)<zero) then
        res=delta*k15
    else if(delta*abs(k15-g7)>eps) then
        res=GK_Rec(f,xMin,av,eps)+GK_Rec(f,av,xMax,eps)
    else
        res=delta*k15
    end if

end function GK2041_Rec
!!!!-------------------------------------------------------------------------------------------------------------

end module IntegrationRoutines
