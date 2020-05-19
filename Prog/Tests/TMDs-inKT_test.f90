!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program runs the test functions of TMDs_inKT and compare to exact values
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_Numerics
use aTMDe_control
use TMDs_inKT
implicit none

integer::i,iMax
real*8::qTMax,qTStep
real*8::qT(1:45)
real*8::TT(-5:5)
real*8::ff(-5:5),x

call artemide_Initialize('const-TMD-inKT',prefix='Prog/Tests/const-files/')
! call artemide_Initialize('const-DY_LO',prefix='/home/alexey/artemide_Repository/Constants-files/')!
call artemide_SetNPparameters_TMDR((/1.86041d0, 0.029551d0/))
call artemide_SetNPparameters_uTMDPDF((/0.253434d0, 9.04351d0, 346.999d0, 2.47992d0, -5.69988d0, 0.d0, 0.d0/))

x=0.1d0

iMax=45

qt=(/0.001d0,0.005d0,0.01d0,0.025d0,0.05d0,0.1d0,0.25d0,0.5d0,0.75d0,1d0,&
1.5d0,2d0,3d0,4d0,5d0,6d0,7d0,8d0,9d0,10d0,&
12d0,14d0,16d0,18d0,20d0,23d0,26d0,29d0,32d0,35d0,38d0,41d0,45d0,49d0,&
53d0,57d0,61d0,65d0,70d0,75d0,80d0,85d0,90d0,95d0,100d0/)

write(*,*) '-----------------------------------------------------------------------------------------'
write(*,*) '--------------------Test program for TMDs_inKT-------------------------------------------'
write(*,*) '--------------------Evaluates Fourier integrals and compare to exact values--------------'
write(*,*) '--------------------Shows ratios of test/exact ------------------------------------------'
write(*,*) '-----------------------------------------------------------------------------------------'


write(*,*) '------------------ Exponential functions --(ratio)-----------------'
do i=1,iMax
    TT=testTMD_kT(x,qT(i))
    ff(-5)=3d0/pix2*(3d0*qt(i)**4-24d0*qt(i)**2*x**2+8d0*x**4)/(qT(i)**2d0+x**2)**(4.5d0)
    ff(-4)=3d0*x/pix2*(2d0*x**2-3d0*qt(i)**2)/(qT(i)**2d0+x**2)**(3.5d0)
    ff(-3)=1d0/pix2*(2d0*x**2-qt(i)**2)/(qT(i)**2d0+x**2)**(2.5d0)
    ff(-2)=x/pix2/(qT(i)**2d0+x**2)**(1.5d0)
    write(*,*) qT(i),TT(-5)/ff(-5),TT(-4)/ff(-4),TT(-3)/ff(-3),TT(-2)/ff(-2)
    
    call TMDs_inKT_ResetCounters()
end do


write(*,*) '------------------ Gaussian functions --(difference is show, since gauss decay too fast)----'
do i=1,15
    TT=testTMD_kT(x,qT(i))
    ff(2)=(4d0*x-qt(i)**2)/(16d0*pi*x**3)*exp(-qt(i)**2/(4d0*x))
    ff(3)=(1d0)/(4d0*pi*x)*exp(-qt(i)**2/(4d0*x))
    write(*,*) qT(i),TT(2)-ff(2),TT(3)-ff(3)
    
    call TMDs_inKT_ResetCounters()
end do

write(*,*) '------------------ (oscilating+Weakly decaying) and (Weakly decaying) --(ratio)-------------'
do i=1,iMax
    TT=testTMD_kT(x,qT(i))
    ff(4)=Bessel_J0(2*sqrt(qt(i)*x))/(pix2*qt(i))
    ff(5)=1/(pix2*qt(i))
    write(*,*) qT(i),TT(4)/ff(4),TT(5)/ff(5)
    
    call TMDs_inKT_ResetCounters()
end do

! x=0.01d0
! do i=1,iMax
! 
!     TT=uTMDPDF_kT_5(x,qT(i),1)
!     write(*,*) "{",qT(i),",", TT(1),"},"
! end do


end program example
