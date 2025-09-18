!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                          !!
!! Example code for computation of evolution with Snowflake library                         !!
!!                                          A.Vladimirov 01.04.2024                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program example
use SnowFlake
implicit none

real*8::x1,x2,x3
real*8::mu0,mu1,Q
real*8::r1,r2,r3

!!! Prior to run Snowslace initialize it with INI-file.
!!! This procedure call for computation of integration kernels. It can take up to several minutes (depending on grid size)
!!! Default-configuration of INI-file can be found in the root directory of SnowFlake
!!! Here I initialize with TEST.ini file, which uses unreliably small grid
!!!
!!! NOTE: initialization procedure prints list of all parameters in terminal. It can be switched-off in the INI-file (option 0).
call  SnowFlake_Initialize("TEST.ini","prog/")
!call  SnowFlake_Initialize("Snowflake.ini")

!!! Request to compute the evolution of a configuration from mu0 to mu1. Both mu's are in GeV.
!!!
!!! Argument alpha is alpha_s=g^2/4pi for QCD. It must be a function with following inteface (see example in the end of the code)
!!!function alpha(mu)
!!!        real*8,intent(in)::mu
!!!        real*8::alpha
!!!
!!! Rest arguments represent the boundary conditions each of them must be a function of two real*8 variables f(x,y)
!!! These arguments are optional. If an argument is not presented it is interpreted as f(x,y)=0
!!! There are following arguments:
!!! G1,U1,D1,S1,C1,B1,G2,U2,D2,S2,C2,B2
!!! they correspond to different flavors and symmetry.
!!! G=gluon, U=u-quark, D=d-quark, S=s-quark, C=c-quark, B=b-quark
!!! The type 1 or 2 is determied by inputQ and inputG:
!!! If inputQ="T" 1=T(x1,x2,x3), 2=\DeltaT(x1,x2,x3)
!!! If inputQ="S" 1=S^+(x1,x2,x3), 2=S^-(x1,x2,x3)
!!! If inputQ="C" 1=\frak{S}^+(x1,x2,x3), 2=\frak{S}^-(x1,x2,x3) (default)
!!! If inputG="T" 1=T^+_{3F}(x1,x2,x3), 2=T^{-}_{3F}(x1,x2,x3)
!!! If inputG="C" 1=\frak{T}^+(x1,x2,x3), 2=\frak{T}^-(x1,x2,x3) (default)
!!!
!!! IMPORTANT: boundary conditions MUST satisfy physical symmetries. Otherwise, result is not correct.
!!! There is no check for the symmetry.
!!!
!!! The evolution results are stored in memory (from mu0 to mu1 with a step sepcified in INI)
!!! such that aftewards the value at some Q, can be extracted by interpolation
mu0=1.d0
mu1=25.d0
call ComputeEvolution(mu0,mu1,alpha,U1=initialF,U2=initialA,G1=initialG,inputQ="T",inputG="T")

!!! The grid resulting after evolution is stored in memory.
!!! To access it call GetPDF(x1,x2,f)
!!! where arguments x1,x2 are x1,x2 (real*8*)
!!! f is integer indicating the type and flavor of PDF.
!!! It follows the LHAPDF pattern (1=d,2=u,3=s,4=c,5=b, 0 or 10=gluon), and depends on outputT parameter (optional)
!!! for outputT="T"
!!! f = -10 [T_{3F}^-(gluon)], -5...-1 [\Delta T(quark)], 0=10 [T_{3F}^+(gluon)], 1...5 [T(quark)]
!!! for outputT="S"
!!! f = -10 [T_{3F}^-(gluon)], -5...-1 [S^-(quark)], 0=10 [T_{3F}^+(gluon)], 1...5 [S^+(quark)]
!!! for outputT="C" (default)
!!! f = -10 [\frak{T}^-(gluon)], -5...-1 [\frak{S}^-(quark)], 0=10 [\frak{T}^+(gluon)], 1...5 [\frak{T}^+(quark)]

!!! Here I call for result of Evolution for T(0.2,-0.3,0.1) for u quark
x1=0.2d0
x2=-0.3d0
Q=10.d0
r1=GetPDF(x1,x2,Q,2,outputT='T')
!!! Here I call for result of Evolution for S^+(0.2,-0.3,0.1) for u quark
r2=GetPDF(x1,x2,Q,2,outputT='S')
!!! Here I call for result of Evolution for frak{S}^+(0.2,-0.3,0.1) for u quark
r3=GetPDF(x1,x2,Q,2,outputT='C')
!!! The values are
write(*,*) "1 GeV->10 GeV. U-quark >>>",r1,r2,r3

Q=20
r1=GetPDF(x1,x2,Q,2,outputT='T')
!!! Here I call for result of Evolution for S^+(0.2,-0.3,0.1) for u quark
r2=GetPDF(x1,x2,Q,2,outputT='S')
!!! Here I call for result of Evolution for frak{S}^+(0.2,-0.3,0.1) for u quark
r3=GetPDF(x1,x2,Q,2,outputT='C')
!!! The values are
write(*,*) "1 GeV->20 GeV. U-quark >>>",r1,r2,r3

!!! Because of the symmetry T(x1,x2,x3)=T(-x3,-x2,-x1) and DeltaT(x1,x2,x3)=DeltaT(-x3,-x2,-x1)
!!! let me check it
r1=GetPDF(0.2d0,-0.3d0,Q,2,outputT='T')
r2=GetPDF(-0.1d0,0.3d0,Q,2,outputT='T')
write(*,*) "These numbers should be same",r1,r2
r1=GetPDF(0.2d0,-0.3d0,Q,-2,outputT='T')
r2=GetPDF(-0.1d0,0.3d0,Q,-2,outputT='T')
write(*,*) "These numbers should be of opposite sign",r1,r2

!!!!Same procedure applies to chiral-odd quark distributions
!!!! In this case one does not need to specify the type of input, since it can be only of single type.
!!!! Also the boundary conditions are different only by flavors. Because evolution equation for H and E is the same.
!!!! thus boundary conditions arguments are only U1,D1,S1,C1,B1

call ComputeEvolutionChiralOdd(mu0,mu1,alpha,D1=initialF)

!!! To get the result call GetPDFChiralOdd(x1,x2,f)
!!! with f=1,2,3,4,5 only (1=d,2=u,3=s,4=c,5=b)
x1=0.2d0
x2=-0.3d0
Q=10.d0
r1=GetPDFChiralOdd(x1,x2,Q,1)
!!!! Since u-quark was not introduced this number should be zero
r2=GetPDFChiralOdd(x1,x2,Q,2)
write(*,*) "1 GeV->10 GeV. (chiral odd) D,U-quarks >>>",r1,r2

Q=20.d0
r1=GetPDFChiralOdd(x1,x2,Q,1)
!!!! Since u-quark was not introduced this number should be zero
r2=GetPDFChiralOdd(x1,x2,Q,2)
write(*,*) "1 GeV->20 GeV. (chiral odd) D,U-quarks >>>",r1,r2

contains

    !!!! some symmetric function
    function initialF(x,y)
        real*8::x,y,initialF

        if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
            initialF=(1-x**2)*(1-y**2)*(1-(x+y)**2)
        else
            initialF=0.d0
        end if
    end function initialF

    !!!! some asymmetric function
    function initialA(x,y)
        real*8::x,y,initialA
        real*8,parameter::pi=3.141592653589793d0

        if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
            initialA=sin(pi*y)*cos(pi/2*x)*cos(pi/2*(x+y))
        else
            initialA=0.d0
        end if
    end function initialA

    !!!! some symmetric function
    function initialG(x,y)
        real*8::x,y,initialG

        if(abs(x)<1 .and. abs(y)<1 .and. abs(x+y)<1) then
            initialG=sin(2*x+y)*(1-x**2)*(1-y**2)*(1-(x+y)**2)/sqrt(max(abs(x),abs(y),abs(x+y)))
        else
            initialG=0.d0
        end if
    end function initialG

    !!!!!! Function for alpha_s of QCD.
    !!!!!! Here is a simple model for alpha-s. You can use any other model, or interface it with LHAPDF, or other code
    !!!!!! Preserve the interface.
    pure function alpha(mu)
        real*8,intent(in)::mu
        real*8::alpha
        alpha=12.566370614359172d0/11/(2*log(mu)+3.d0)
    end function alpha
end program example
