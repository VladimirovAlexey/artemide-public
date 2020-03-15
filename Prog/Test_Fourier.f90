!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program computes Fouries of test cases and compare to exact values
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_control
use TMDF
implicit none

integer,parameter::Nt=18
real*8::testString(1:Nt),qT(1:Nt)
real*8::par1,par2,par3,tt,tt1
integer::j

call artemide_Initialize('const-TEST',prefix='/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Prog/')

!----------------------------------------------------------------------------------------------------------------------------------
qT(1)=0.001d0
qT(2)=0.01d0
qT(3)=0.05d0
qT(4)=0.1d0
qT(5)=0.5d0
qT(6)=1d0
qT(7)=2d0
qT(8)=4d0
qT(9)=6d0
qT(10)=10d0
qT(11)=15d0
qT(12)=20d0
qT(13)=25d0
qT(14)=30d0
qT(15)=40d0
qT(16)=50d0
qT(17)=75d0
qT(18)=90d0

!----------------------------------------------------------------------------------------------------------------------------------
write(*,*) '------------------------Exponential test-------------------------------------------------------------------------'
write(*,*) '>>>> mu=0.2, x1=0, x2=0 '
par1=0.2d0
par2=0d0
par3=0d0

! tt=TMDF_F(1d0,2d0,par2,par3,par1,1d0,1d0,9999)
! 
! stop

do j=1,Nt
  tt=TMDF_F(1d0,qT(j),par2,par3,par1,1d0,1d0,9999)
  tt1=Exact9999(qT(j),par1,par2,par3)
  write(*,"(F14.6,'   ')",advance='no') tt/tt1
  
  call TMDF_ResetCounters()
end do
write(*,*)

write(*,*) '>>>> mu=0.5, x1=0, x2=0 '
par1=0.5d0
par2=0d0
par3=0d0

do j=1,Nt
  tt=TMDF_F(1d0,qT(j),par2,par3,par1,1d0,1d0,9999)
  tt1=Exact9999(qT(j),par1,par2,par3)
  write(*,"(F14.6,'   ')",advance='no') tt/tt1
  
  call TMDF_ResetCounters()
end do
write(*,*)

write(*,*) '------------------------Gaussian test-------------------------------------------------------------------------'
write(*,*) '>>>> mu=0.2, x1=0, x2=0 '
par1=0.2d0
par2=0d0
par3=0d0
do j=1,Nt
  tt=TMDF_F(1d0,qT(j),par2,par3,par1,1d0,1d0,9998)
  tt1=Exact9998(qT(j),par1,par2,par3)
  if(abs(tt1)>0d0) then
    write(*,"(F14.6,'   ')",advance='no') tt/tt1
  else 
    write(*,"(F14.6,'/0   ')",advance='no') tt
  end if
  
  call TMDF_ResetCounters()
end do
write(*,*)

write(*,*) '>>>> mu=0.5, x1=0, x2=0 '
par1=0.5d0
par2=0d0
par3=0d0
do j=1,Nt
  tt=TMDF_F(1d0,qT(j),par2,par3,par1,1d0,1d0,9998)
  tt1=Exact9998(qT(j),par1,par2,par3)
  if(abs(tt1)>0d0) then
    write(*,"(F14.6,'   ')",advance='no') tt/tt1
  else 
    write(*,"(F14.6,'/0   ')",advance='no') tt
  end if
  
  call TMDF_ResetCounters()
end do
write(*,*)

contains

!!! the exact value of test function 9999 at qt,mu,x1,x2
function Exact9999(qT_in,mu,x1,x2)
  real*8::Exact9999,qT_in,mu,x1,x2
  real*8::X
  X=(qT_in/mu)**2
  
  Exact9999=1d0/(2d0*(mu**2)*(1+X)**(3d0/2d0))*(1d0+x1/mu**2*(6d0-9d0*X)/(1+X)**2+15d0*x2/mu**4*(8d0-40d0*X+15d0*X**2)/(1+X)**4)

end function Exact9999

!!! the exact value of test function 9998 at qt,mu,x1,x2
function Exact9998(qT_in,mu,x1,x2)
  real*8::Exact9998,qT_in,mu,x1,x2
  real*8::Y
  Y=qT_in**2/4d0/mu
  Exact9998=Exp(-Y)/4d0/mu*(1d0+x1/mu*(1d0-Y)+x2/mu**2*(2d0-4d0*Y+Y**2))

end function Exact9998

end program example
  
