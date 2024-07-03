!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot for AS-term
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


! module KTtest
! INCLUDE '/data/arTeMiDe_Repository/artemide/src/Code/KTspace/Fourier_Levin.f90'
! end module KTtest


! module KTtest2
! INCLUDE '/data/arTeMiDe_Repository/artemide/src/Code/KTspace/grid_inKT(new).f90'
! end module KTtest2

program example
use aTMDe_control
use uTMDPDF
use Fourier_Levin_uTMDPDF

implicit none

integer::i,j
real*8,dimension(-5:5)::FinQ
real*8,dimension(1:5,0:16,-5:5)::FinQ2
real*8::kT,x,mu
integer::f

call artemide_Initialize('ART23_MSHT_N4LL.atmde',prefix='Prog/Tests/const-files/')
call artemide_SetNPparameters_uTMDPDF((/&
0.874245d0,0.913883d0,0.991563d0,6.05412d0,&
0.353908d0,46.6064d0,0.115161d0,1.53235d0,&
1.31966d0,0.434833d0, 0.0d0, 0.0d0/))


! call testGrid_inKT()
! stop

! ! FinQ2=Fourier_Levin_array(toFourier)
! ! !write(*,*) Fourier_Levin_array(toFourier)
! ! !write(*,*) FinQ2
! !
! ! write(*,*) "-------- kT values"
! ! do i=1,5
! ! do j=0,16
! ! write(*,'(F8.4,",")',advance="no") FinQ2(i,j,2)
! ! end do
! ! end do
! ! write(*,*)
! !
! ! contains
! !
! ! function toFourier(b)
! ! real*8::toFourier(-5:5)
! ! real*8,intent(in)::b
! !
! ! toFourier=Exp(-b)
! ! end function toFourier
x=0.003d0
f=-2
mu=38.50d0


write(*,*) "--------------------in KT -------- mu=",mu
do i=1,20
kT=0.1*i
FinQ=uTMDPDF_inKT(x,kT,mu,mu**2,1)
write(*,'("{",F5.2,",",F12.8,"},")',advance="no") kT,kT**2*FinQ(f)
end do

do i=2,99
kT=1.d0+i
FinQ=uTMDPDF_inKT(x,kT,mu,mu**2,1)
write(*,'("{",F5.1,",",F12.8,"},")',advance="no") kT,kT**2*FinQ(f)
end do
write(*,*) " "

write(*,*) "--------------------in KT from grid -- mu=",mu
do i=1,20
kT=0.1*i
FinQ=uTMDPDF_inKT(x,kT,mu,1)
write(*,'("{",F5.2,",",F12.8,"},")',advance="no") kT,kT**2*FinQ(f)
end do

do i=2,99
kT=1.d0+i
FinQ=uTMDPDF_inKT(x,kT,mu,1)
write(*,'("{",F5.1,",",F12.8,"},")',advance="no") kT,kT**2*FinQ(f)
end do
write(*,*) " "

end program example
