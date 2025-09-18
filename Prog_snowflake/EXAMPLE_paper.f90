!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                          !!
!! Example code for computation of example figure in the paper using  Snowflake library     !!
!!                                          A.Vladimirov 01.04.2024                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program example
use SnowFlake
implicit none

real*8::x1,x2,x3,r1
real*8::mu0,mu1
integer::i

!!! Prior to run Snowslace initialize it with INI-file.
!!! This procedure call for computation of integration kernels. It can take up to several minutes (depending on grid size)
!!! Default-configuration of INI-file can be found in the root directory of SnowFlake
!!! Here I initialize with TEST.ini file, which uses unreliably small grid
!!!
!!! NOTE: initialization procedure prints list of all parameters in terminal. It can be switched-off in the INI-file (option 0).
!call  SnowFlake_Initialize("TEST.ini","prog/")
call  SnowFlake_Initialize("Snowflake.ini")

!!! Request evolution of a configuration from mu0 to mu1. Both mu's are in GeV.
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
mu0=1.d0
mu1=100.d0
call ComputeEvolution(mu0,mu1,alpha,&
    U1=inTU,U2=inDeltaTU,D1=inTD,D2=inDeltaTD,S1=inTS,S2=inDeltaTS,G1=inGp,G2=inGm,inputQ="T",inputG="T")

do i=-50,50
    if(i==0) cycle
    x1=i/100.d0
    x2=-2*x1
    x3=x1

    r1=GetPDF(x1,x2,mu1,2,outputT='T')
    write(*,'("{",F6.3,", ",F12.8,"},")',advance='no') x1,r1
end do
write(*,*) " "
contains

    !!!! some symmetric function
    function phi(x1,x2,x3)
        real*8::x1,x2,x3,phi

        if(abs(x1)<1 .and. abs(x2)<1 .and. abs(x3)<1) then
            phi=(1-x1**2)*(1-x2**2)*(1-x3**2)
        else
            phi=0.d0
        end if
    end function phi

    !!!! radious |x|_inf
    function r(x1,x2,x3)
        real*8::x1,x2,x3,r

        r=max(abs(x1),abs(x2),abs(x3))

    end function r

    !!!! Tu initial
    function inTU(x,y)
        real*8::x,y,inTU

        inTU=cos(4*y)*phi(x,y,-x-y)
    end function inTU

    !!!! Td initial
    function inTD(x,y)
        real*8::x,y,inTD
        real*8,parameter::pi=3.141592653589793d0

        inTD=(2-cos(3*pi*phi(x,y,-x-y)))*phi(x,y,-x-y)
    end function inTD

    !!!! DeltaTu initial
    function inDeltaTU(x,y)
        real*8::x,y,inDeltaTU
        real*8,parameter::pi=3.141592653589793d0

        inDeltaTU=(sin(pi*y)+4*(x**2-(x+y)**2))*phi(x,y,-x-y)
    end function inDeltaTU

    !!!! DeltaTd initial
    function inDeltaTD(x,y)
        real*8::x,y,inDeltaTD
        real*8,parameter::pi=3.141592653589793d0

        inDeltaTD=2*sin(pi*y)/r(x,y,-x-y)*(1-cos(phi(x,y,-x-y)))
    end function inDeltaTD

    !!!! Ts initial
    function inTS(x,y)
        real*8::x,y,inTS

        inTS=-0.3d0*inTD(x,y)
    end function inTS

    !!!! DeltaTs initial
    function inDeltaTS(x,y)
        real*8::x,y,inDeltaTS

        inDeltaTS=-0.3d0*inDeltaTD(x,y)
    end function inDeltaTS

    !!!! plus function for gluon
    function inGp(x,y)
        real*8::x,y,inGp
        inGp=phi(x,y,-x-y)*r(x,y,-x-y)*sin(x-(-x-y))
    end function inGp

    !!!! plus function for gluon
    function inGm(x,y)
        real*8::x,y,inGm
        inGm=phi(x,y,-x-y)*r(x,y,-x-y)*cos(x-(-x-y))
    end function inGm

    !!!!!! Function for alpha_s of QCD.
    !!!!!! Here is a simple model for alpha-s. You can use any other model, or interface it with LHAPDF, or other code
    !!!!!! Preserve the interface.
    pure function alpha(mu)
        real*8,intent(in)::mu
        real*8::alpha
        alpha=12.566370614359172d0/11/(2*log(mu)+3.d0)
        !alpha=1.d0/11/(2*log(mu)+3.d0)
    end function alpha
end program example
