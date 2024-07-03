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

implicit none

integer::i,j
real*8,dimension(-5:5)::FinQ
real*8::kT,x,mu
integer::f

call artemide_Initialize('ART23_MSHT_N4LL.atmde',prefix='Prog/Tests/const-files/')
call artemide_SetNPparameters_uTMDPDF((/&
0.874245d0,0.913883d0,0.991563d0,6.05412d0,&
0.353908d0,46.6064d0,0.115161d0,1.53235d0,&
1.31966d0,0.434833d0, 0.0d0, 0.0d0/))

x=0.01d0
f=-2
mu=10.d0



write(*,*) "--------------------in B -------- OPT"
kT=0.00001
FinQ=uTMDPDF_inB(x,kT,1)
write(*,'("{",F5.2,",",F12.8,"},")',advance="no") kT,FinQ(f)

do i=1,20
kT=0.01*i
FinQ=uTMDPDF_inB(x,kT,1)
write(*,'("{",F5.2,",",F12.8,"},")',advance="no") kT,FinQ(f)
end do

do i=1,98
kT=0.2d0+i*0.1
FinQ=uTMDPDF_inB(x,kT,1)
write(*,'("{",F5.1,",",F12.8,"},")',advance="no") kT,FinQ(f)
end do
do i=1,20
kT=10.d0+i
FinQ=uTMDPDF_inB(x,kT,1)
write(*,'("{",F5.1,",",F12.8,"},")',advance="no") kT,FinQ(f)
end do
write(*,*) " "

write(*,*) "--------------------in KT -------- OPT"
do i=1,20
kT=0.1*i
FinQ=uTMDPDF_inKT(x,kT,1)
write(*,'("{",F5.2,",",F12.8,"},")',advance="no") kT,kT**2*FinQ(f)
end do

do i=2,99
kT=1.d0+i
FinQ=uTMDPDF_inKT(x,kT,1)
write(*,'("{",F5.1,",",F12.8,"},")',advance="no") kT,kT**2*FinQ(f)
end do
write(*,*) " "

write(*,*) "--------------------in B -------- mu=",mu
kT=0.00001
FinQ=uTMDPDF_inB(x,kT,mu,mu**2,1)
write(*,'("{",F5.2,",",F12.8,"},")',advance="no") kT,FinQ(f)

do i=1,20
kT=0.01*i
FinQ=uTMDPDF_inB(x,kT,mu,mu**2,1)
write(*,'("{",F5.2,",",F12.8,"},")',advance="no") kT,FinQ(f)
end do

do i=1,98
kT=0.2d0+i*0.1
FinQ=uTMDPDF_inB(x,kT,mu,mu**2,1)
write(*,'("{",F5.1,",",F12.8,"},")',advance="no") kT,FinQ(f)
end do
do i=1,20
kT=10.d0+i
FinQ=uTMDPDF_inB(x,kT,mu,mu**2,1)
write(*,'("{",F5.1,",",F12.8,"},")',advance="no") kT,FinQ(f)
end do
write(*,*) " "

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

write(*,*) "- AS function of Q =-"
x=0.001d0
f=2
kT=1.d0
do i=0,40
mu=exp(0.1324579341637009*i)
FinQ=uTMDPDF_inKT(x,kT,mu,mu**2,1)
write(*,'("{",F5.1,",",F12.8,"},")',advance="no") mu,FinQ(f)
end do
write(*,*) " "

end program example
