!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a MSHT plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module Grid_Test
INCLUDE '/data/arTeMiDe_Repository/artemide/src/Code/Twist2/Twist2_ChGrid.f90'
end module Grid_Test



program example
use Grid_Test
use aTMDe_control
use uTMDPDF_OPE
implicit none

real*8::qq(-5:5),ff(-5:5)

integer::n,k
real*8::x,b,time1,time2
real*8::bList(1:6),xlist(1:5)

call artemide_Initialize('ART23_MSHT_N4LL.atmde',prefix='Prog/Tests/const-files/')
call artemide_SetNPparameters_uTMDPDF((/&
0.874245d0,0.913883d0,0.991563d0,6.05412d0,&
0.353908d0,46.6064d0,0.115161d0,1.53235d0,&
1.31966d0,0.434833d0, 0.0d0, 0.0d0/))

call Twist2_ChGrid_Initialize((/0.00001d0,0.001d0,0.1d0,0.5d0,1.d0/),&
                              (/0.00001d0,0.01d0,0.2d0,2.d0,25.d0/),32,32,1,.true.,"JAJA",2)
call Twist2_ChGrid_MakeGrid(F)

bList=(/0.01d0,0.1d0,0.2d0,0.7d0,1.d0,10d0/)
xList=(/0.0001d0,0.001d0,0.1d0,0.2d0,0.5d0/)

call cpu_time(time1)
write(*,*) "---------------- vsX --------------------------"
!open (22, FILE='/data/WorkingFiles/TMD/Fit_Notes/MathFiles/Tests/data/ART23_grid(MSHT_N4LL).dat', STATUS='REPLACE')
!open (22, FILE='/data/WorkingFiles/TMD/Fit_Notes/MathFiles/Tests/data/exact(MSHT_N4LL)_vsX.dat', STATUS='REPLACE')
open (22, FILE='/data/WorkingFiles/TMD/Fit_Notes/MathFiles/Tests/data/CH1_grid(MSHT_N4LL)_vsX.dat', STATUS='REPLACE')
do k=1,6
do n=0,800
    !b=exp(0.6*k)-1.
    b=bList(k)
    x=10**(-0.005*n)

    !FF=uTMDPDF_OPE_convolution(x,b,1,.true.)
    FF=ExtractFromGrid(x,b,1)/x
    write(22,*) x,b,FF(-5),FF(-4),FF(-3),FF(-2),FF(-1),FF(0),FF(1),FF(2),FF(3),FF(4),FF(5)
end do
end do
close(22)
call cpu_time(time2)
write(*,*) "--->TIME<---",time2-time1

call cpu_time(time1)
write(*,*) "---------------- vsB --------------------------"
!open (22, FILE='/data/WorkingFiles/TMD/Fit_Notes/MathFiles/Tests/data/ART23_grid(MSHT_N4LL).dat', STATUS='REPLACE')
!open (22, FILE='/data/WorkingFiles/TMD/Fit_Notes/MathFiles/Tests/data/exact(MSHT_N4LL)_vsB.dat', STATUS='REPLACE')
open (22, FILE='/data/WorkingFiles/TMD/Fit_Notes/MathFiles/Tests/data/CH1_grid(MSHT_N4LL)_vsB.dat', STATUS='REPLACE')
do k=1,5
do n=0,800
    !b=exp(0.6*k)-1.
    x=xList(k)
    b=0.01d0*n

    !FF=uTMDPDF_OPE_convolution(x,b,1,.true.)
    FF=ExtractFromGrid(x,b,1)/x
    write(22,*) x,b,FF(-5),FF(-4),FF(-3),FF(-2),FF(-1),FF(0),FF(1),FF(2),FF(3),FF(4),FF(5)
end do
end do
close(22)
call cpu_time(time2)
write(*,*) "--->TIME<---",time2-time1

call cpu_time(time1)
write(*,*) "---------------- vsB --------------------------"
!open (22, FILE='/data/WorkingFiles/TMD/Fit_Notes/MathFiles/Tests/data/ART23_grid(MSHT_N4LL).dat', STATUS='REPLACE')
!open (22, FILE='/data/WorkingFiles/TMD/Fit_Notes/MathFiles/Tests/data/exact(MSHT_N4LL).dat', STATUS='REPLACE')
open (22, FILE='/data/WorkingFiles/TMD/Fit_Notes/MathFiles/Tests/data/CH1_grid(MSHT_N4LL).dat', STATUS='REPLACE')
do k=0,200
do n=0,200
    !b=exp(0.6*k)-1.
    x=10**(-k*0.02)
    b=0.01d0+0.05d0*n

    !FF=uTMDPDF_OPE_convolution(x,b,1,.true.)
    FF=ExtractFromGrid(x,b,1)/x
    write(22,*) x,b,FF(-5),FF(-4),FF(-3),FF(-2),FF(-1),FF(0),FF(1),FF(2),FF(3),FF(4),FF(5)
end do
end do
close(22)
call cpu_time(time2)
write(*,*) "--->TIME<---",time2-time1


contains

function F(x,b,h,wG)
  real*8,dimension(-5:5)::F
  real*8,intent(in)::x,b
  integer,intent(in)::h
  logical,intent(in)::wG

  !F=(4.34*x**(-0.15)*(1-x)**9.11-1.048*x**(-0.167)*(1-x)**25.)*(1+log(b)+3.4*log(b)**2)/cosh(0.3*b)
  F=x*uTMDPDF_OPE_convolution(x,b,h,.true.)
end function F

end program example
