!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use IntegrationRoutines
implicit none

real*8:: re, r1,r2,r3,r4,r5,r6,r7,r8

write(*,*) '-------------------------------------------------------------------'
write(*,*) '---   Test of the IntegrationRoutines submodule                 ---'
write(*,*) '-------------------------------------------------------------------'
write(*,*) ' '
write(*,*) '| Exact          | S5             | SN(n=10)       | SN(n=100)      |',&
            ' SA(eps=0.1)    | SA(eps=0.001)  | G7             | K15            | GK (eps=0.001) |'
write(*,*) '|                |                |                |                |',&
'                |                |'
write(*,*) '1) (0.5+x)^3   [-1,1]  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -',&
'  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -|'

re=5d0/4d0
r1=Integrate_S5(f1,-1d0,1d0)
r2=Integrate_SN(f1,-1d0,1d0,10)
r3=Integrate_SN(f1,-1d0,1d0,100)
r4=Integrate_SA(f1,-1d0,1d0,0.1d0)
r5=Integrate_SA(f1,-1d0,1d0,0.001d0)
r6=Integrate_G7(f1,-1d0,1d0)
r7=Integrate_K15(f1,-1d0,1d0)
r8=Integrate_GK(f1,-1d0,1d0,0.001d0)
write(*,'(" | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | "&
,F14.10," | ")')&
re,r1,r2,r3,r4,r5,r6,r7,r8
write(*,'(" | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | "&
,F14.10," | ")')&
1d0,r1/re,r2/re,r3/re,r4/re,r5/re,r6/re,r7/re,r8/re

write(*,*) '|                |                |                |                |',&
'                |                |                |                |'
write(*,*) '2) f= (0.5+x)^6   [-1,1]  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -',&
'  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -|'

re=547d0/224d0
r1=Integrate_S5(f2,-1d0,1d0)
r2=Integrate_SN(f2,-1d0,1d0,10)
r3=Integrate_SN(f2,-1d0,1d0,100)
r4=Integrate_SA(f2,-1d0,1d0,0.1d0)
r5=Integrate_SA(f2,-1d0,1d0,0.001d0)
r6=Integrate_G7(f2,-1d0,1d0)
r7=Integrate_K15(f2,-1d0,1d0)
r8=Integrate_GK(f2,-1d0,1d0,0.001d0)
write(*,'(" | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | "&
,F14.10," | ")')&
re,r1,r2,r3,r4,r5,r6,r7,r8
write(*,'(" | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | "&
,F14.10," | ")')&
1d0,r1/re,r2/re,r3/re,r4/re,r5/re,r6/re,r7/re,r8/re

write(*,*) '|                |                |                |                |',&
'                |                |                |                |'
write(*,*) '3) f= sqrt(1+x^2)   [-1,1]   -  -  -  -  -  -  -  -  -  -  -  -  -  -',&
'  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -|'

re=1.5707963267948966d0
r1=Integrate_S5(f3,-1d0,1d0)
r2=Integrate_SN(f3,-1d0,1d0,10)
r3=Integrate_SN(f3,-1d0,1d0,100)
r4=Integrate_SA(f3,-1d0,1d0,0.1d0)
r5=Integrate_SA(f3,-1d0,1d0,0.001d0)
r6=Integrate_G7(f3,-1d0,1d0)
r7=Integrate_K15(f3,-1d0,1d0)
r8=Integrate_GK(f3,-1d0,1d0,0.001d0)
write(*,'(" | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | "&
,F14.10," | ")')&
re,r1,r2,r3,r4,r5,r6,r7,r8
write(*,'(" | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | "&
,F14.10," | ")')&
1d0,r1/re,r2/re,r3/re,r4/re,r5/re,r6/re,r7/re,r8/re

write(*,*) '|                |                |                |                |',&
'                |                |                |                |'
write(*,*) '4) f= Exp(-5 x^2)   [-1,1]   -  -  -  -  -  -  -  -  -  -  -  -  -  -',&
'  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -|'

re=0.791424619221027d0
r1=Integrate_S5(f4,-1d0,1d0)
r2=Integrate_SN(f4,-1d0,1d0,10)
r3=Integrate_SN(f4,-1d0,1d0,100)
r4=Integrate_SA(f4,-1d0,1d0,0.1d0)
r5=Integrate_SA(f4,-1d0,1d0,0.001d0)
r6=Integrate_G7(f4,-1d0,1d0)
r7=Integrate_K15(f4,-1d0,1d0)
r8=Integrate_GK(f4,-1d0,1d0,0.001d0)
write(*,'(" | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | "&
,F14.10," | ")')&
re,r1,r2,r3,r4,r5,r6,r7,r8
write(*,'(" | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | "&
,F14.10," | ")')&
1d0,r1/re,r2/re,r3/re,r4/re,r5/re,r6/re,r7/re,r8/re

write(*,*) '|                |                |                |                |',&
'                |                |                |                |'
write(*,*) '5) f= (0.5+x)^2 Exp(-5 x^2)   [-1,1]  -  -  -  -  -  -  -  -  -  -  -',&
'  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -|'

re=0.2756510273275424d0
r1=Integrate_S5(f5,-1d0,1d0)
r2=Integrate_SN(f5,-1d0,1d0,10)
r3=Integrate_SN(f5,-1d0,1d0,100)
r4=Integrate_SA(f5,-1d0,1d0,0.1d0)
r5=Integrate_SA(f5,-1d0,1d0,0.001d0)
r6=Integrate_G7(f5,-1d0,1d0)
r7=Integrate_K15(f5,-1d0,1d0)
r8=Integrate_GK(f5,-1d0,1d0,0.001d0)
write(*,'(" | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | "&
,F14.10," | ")')&
re,r1,r2,r3,r4,r5,r6,r7,r8
write(*,'(" | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | "&
,F14.10," | ")')&
1d0,r1/re,r2/re,r3/re,r4/re,r5/re,r6/re,r7/re,r8/re

write(*,*) '|                |                |                |                |',&
'                |                |                |                |'
write(*,*) '6) f= sin(x)+0.5 Cos(4 x)+0.1 Cos(8x)   [-3.5,3.5]   -  -  -  -  -  -',&
'  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -  -|'

re=0.2544244836314143d0
r1=Integrate_S5(f6,-3.5d0,3.5d0)
r2=Integrate_SN(f6,-3.5d0,3.5d0,10)
r3=Integrate_SN(f6,-3.5d0,3.5d0,100)
r4=Integrate_SA(f6,-3.5d0,3.5d0,0.1d0)
r5=Integrate_SA(f6,-3.5d0,3.5d0,0.001d0)
r6=Integrate_G7(f6,-3.5d0,3.5d0)
r7=Integrate_K15(f6,-3.5d0,3.5d0)
r8=Integrate_GK(f6,-3.5d0,3.5d0,0.001d0)
write(*,'(" | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | "&
,F14.10," | ")')&
re,r1,r2,r3,r4,r5,r6,r7,r8
write(*,'(" | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | ",F14.10," | "&
,F14.10," | ")')&
1d0,r1/re,r2/re,r3/re,r4/re,r5/re,r6/re,r7/re,r8/re

write(*,*) '-----------------------------------------------------------------------',&
'-----------------------------------------------------------------------------------'
contains

!!! x^3 Polynomial function
function f1(x)
    real*8::f1,x
    
    f1=(0.5d0+x)**3

end function f1

!!! x^6 Polynomial function
function f2(x)
    real*8::f2,x
    
    f2=(0.5d0+x)**6

end function f2

!!! Squre root
function f3(x)
    real*8::f3,x
    
    f3=Sqrt(1d0-x**2)

end function f3

!!! Gauss
function f4(x)
    real*8::f4,x
    
    f4=Exp(-5d0*x**2)

end function f4

!!! Gauss
function f5(x)
    real*8::f5,x
    
    f5=(0.5d0+x)**2*Exp(-5d0*x**2)

end function f5

!!! Ugly
function f6(x)
    real*8::f6,x
    
    f6=Sin(x)+0.5d0*Cos(4d0*x)+0.1d0*Cos(8d0*x)

end function f6

end program example
