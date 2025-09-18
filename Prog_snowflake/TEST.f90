!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                          !!
!! Simple test-code for Snowlake. Checks compilation integrity and simple evaluation        !!
!!                                          A.Vladimirov 01.04.2024                         !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program SnowflakeTEST
use SnowFlake

implicit none

real*8::t1,t2
!$ real*8::omp_get_wtime
real*8::x1,x2,x3,Q

call  SnowFlake_Initialize("TEST.ini","Prog_snowflake/")

Q=25.d0

call cpu_time(t1)
!$ t1=omp_get_wtime()
call ComputeEvolution(1.d0,Q,alpha,U1=initialF,U2=initialA,G1=initialG,inputQ="T",inputG="T")

call cpu_time(t2)
!$ t2=omp_get_wtime()
write(*,*) "Time for computation of evolution",t2-t1

write(*,*) "  "

write(*,*) "                |    U-QUARK  T(x1,x2,x3)    ||    U-QUARK  deltaT(x1,x2,x3)"
write(*,*) "   (x1,x2,x3)   |  1 GeV       | 25 GeV      ||   1 GeV      | 25 GeV     "
x1=0.3d0
x2=-0.4d0
x3=-x1-x2
write(*,'("(",F4.2,",",F4.2,",",F4.2,") | ",F12.5," | ",F12.5,"|| ",F12.5," | ",F12.5)') x1,x2,x3, &
initialF(x1,x2),GetPDF(x1,x2,Q,2,outputT='T'),initialA(x1,x2),GetPDF(x1,x2,Q,-2,outputT='T')
x1=-0.1d0
x2=.4d0
x3=-x1-x2
write(*,'("(",F4.2,",",F4.2,",",F4.2,") | ",F12.5," | ",F12.5,"|| ",F12.5," | ",F12.5)') x1,x2,x3, &
initialF(x1,x2),GetPDF(x1,x2,Q,2,outputT='T'),initialA(x1,x2),GetPDF(x1,x2,Q,-2,outputT='T')
x1=0.45d0
x2=0.00d0
x3=-x1-x2
write(*,'("(",F4.2,",",F4.2,",",F4.2,") | ",F12.5," | ",F12.5,"|| ",F12.5," | ",F12.5)') x1,x2,x3, &
initialF(x1,x2),GetPDF(x1,x2,Q,2,outputT='T'),initialA(x1,x2),GetPDF(x1,x2,Q,-2,outputT='T')
x1=0.012d0
x2=0.00d0
x3=-x1-x2
write(*,'("(",F4.2,",",F4.2,",",F4.2,") | ",F12.5," | ",F12.5,"|| ",F12.5," | ",F12.5)') x1,x2,x3, &
initialF(x1,x2),GetPDF(x1,x2,Q,2,outputT='T'),initialA(x1,x2),GetPDF(x1,x2,Q,-2,outputT='T')
write(*,*) "----------------------------------------------------------------------------"
write(*,*) "                |    D-QUARK  T(x1,x2,x3)    ||    D-QUARK  deltaT(x1,x2,x3)"
write(*,*) "   (x1,x2,x3)   |  1 GeV       | 25 GeV      ||   1 GeV      | 25 GeV     "
x1=0.3d0
x2=-0.4d0
x3=-x1-x2
write(*,'("(",F4.2,",",F4.2,",",F4.2,") | ",F12.5," | ",F12.5,"|| ",F12.5," | ",F12.5)') x1,x2,x3, &
0.d0,GetPDF(x1,x2,Q,1,outputT='T'),0.d0,GetPDF(x1,x2,Q,-1,outputT='T')
x1=-0.1d0
x2=.4d0
x3=-x1-x2
write(*,'("(",F4.2,",",F4.2,",",F4.2,") | ",F12.5," | ",F12.5,"|| ",F12.5," | ",F12.5)') x1,x2,x3, &
0.d0,GetPDF(x1,x2,Q,1,outputT='T'),0.d0,GetPDF(x1,x2,Q,-1,outputT='T')
x1=0.45d0
x2=0.00d0
x3=-x1-x2
write(*,'("(",F4.2,",",F4.2,",",F4.2,") | ",F12.5," | ",F12.5,"|| ",F12.5," | ",F12.5)') x1,x2,x3, &
0.d0,GetPDF(x1,x2,Q,1,outputT='T'),0.d0,GetPDF(x1,x2,Q,-1,outputT='T')
x1=0.012d0
x2=0.00d0
x3=-x1-x2
write(*,'("(",F4.2,",",F4.2,",",F4.2,") | ",F12.5," | ",F12.5,"|| ",F12.5," | ",F12.5)') x1,x2,x3, &
0.d0,GetPDF(x1,x2,Q,1,outputT='T'),0.d0,GetPDF(x1,x2,Q,-1,outputT='T')
write(*,*) "----------------------------------------------------------------------------"
write(*,*) "                |    Gluon T_{F+}(x1,x2,x3)  ||   Gluon T_{F-}(x1,x2,x3)"
write(*,*) "   (x1,x2,x3)   |  1 GeV       | 25 GeV      ||   1 GeV      | 25 GeV     "
x1=0.3d0
x2=-0.4d0
x3=-x1-x2
write(*,'("(",F4.2,",",F4.2,",",F4.2,") | ",F12.5," | ",F12.5,"|| ",F12.5," | ",F12.5)') x1,x2,x3, &
initialG(x1,x2),GetPDF(x1,x2,Q,0,outputT='T'),0.d0,GetPDF(x1,x2,Q,-10,outputT='T')
x1=-0.1d0
x2=.4d0
x3=-x1-x2
write(*,'("(",F4.2,",",F4.2,",",F4.2,") | ",F12.5," | ",F12.5,"|| ",F12.5," | ",F12.5)') x1,x2,x3, &
initialG(x1,x2),GetPDF(x1,x2,Q,0,outputT='T'),0.d0,GetPDF(x1,x2,Q,-10,outputT='T')
x1=0.45d0
x2=0.00d0
x3=-x1-x2
write(*,'("(",F4.2,",",F4.2,",",F4.2,") | ",F12.5," | ",F12.5,"|| ",F12.5," | ",F12.5)') x1,x2,x3, &
initialG(x1,x2),GetPDF(x1,x2,Q,0,outputT='T'),0.d0,GetPDF(x1,x2,Q,-10,outputT='T')
x1=0.012d0
x2=0.00d0
x3=-x1-x2
write(*,'("(",F4.2,",",F4.2,",",F4.2,") | ",F12.5," | ",F12.5,"|| ",F12.5," | ",F12.5)') x1,x2,x3, &
initialG(x1,x2),GetPDF(x1,x2,Q,0,outputT='T'),0.d0,GetPDF(x1,x2,Q,-10,outputT='T')

write(*,*) "  "
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

    pure function alpha(mu)
        real*8,intent(in)::mu
        real*8::alpha
        alpha=12.566370614359172d0/11/(2*log(mu)+3.d0)
    end function alpha
end program SnowflakeTEST
