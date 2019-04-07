program example
use aTMDe_control
use TMDX_DY
implicit none

!--------------------------------------------------------
!! CALCULATING THE unpolarized DY-cross-section for RHIC
!!
!!
!!leptons-PT > 25 GeV 
!!|leptons-Eta| < 1.1 
!!73 < Minv < 126 GeV 
!!sqrt(s) = 510 GeV
!--------------------------------------------------------
! pt-bins
integer,parameter::numPoints=11
integer,parameter::numReplica=300
!
!-------Kinematics definition
real*8,parameter::ptABSOLUTEMAX=25d0
real*8,parameter::s=510d0**2
integer,dimension(1:3),parameter::proc=(/1,1,5/) !p+p->gamma*+Z
real*8,parameter::etamax=1.1d0
real*8,parameter::ptCut=25d0
real*8,parameter::Qmin=73d0
real*8,parameter::Qmax=126d0
!------
real*8::pmin(1:numPoints),pmax(1:numPoints),X(1:numPoints),Xrep(1:numReplica,1:numPoints)
real*8::Xmean(1:numPoints), Xdev(1:numPoints)
integer::i,j
!$  real*8::OMP_get_wtime,t1,t2,t3,t4

call artemide_Initialize('const-DYfit18_NNLO',prefix='/misc/data2/braun/vla18041/arTeMiDe_Repository/')

!$  t3=OMP_get_wtime()
call artemide_SetReplica_uTMDPDF(0)
call artemide_SetReplica_TMDR(0)

!! make pt-bins
! do i=1,numPoints
!   pmin(i)=ptABSOLUTEMAX*real(i-1)/real(numPoints)
!   pmax(i)=ptABSOLUTEMAX*real(i)/real(numPoints)
! end do
pmin=(/0d0, 1.25d0, 2.5d0, 3.75d0, 5d0, 7.5d0, 10d0, 12.5d0, 15d0, 17.5d0, 20d0/)
pmax=(/1.25d0, 2.5d0, 3.75d0, 5d0, 7.5d0, 10d0, 12.5d0, 15d0, 17.5d0, 20d0,25d0/)

! ! just for check
! do i=1,numPoints
!   call xSec_DY(X1(i),proc,s,(/pmin(i),pmax(i)/),(/Qmin,Qmax/),(/-2d0,2d0/),.true.,(/ptCut,ptCut,-etamax,etamax/))
!   write(*,*) pmin(i),'--',pmax(i),': ',X1(i)
! end do

!$  t1=OMP_get_wtime()
!! setting parameters
call TMDX_DY_SetProcess(proc)
call TMDX_DY_XSetup(s,91d0,0d0)
call TMDX_DY_SetCuts(.true.,ptCut,-etamax,etamax)

!! main replica
call CalcXsec_DY_PTint_Qint_Yint(X,pmin,pmax,Qmin,Qmax)
!$  t2=OMP_get_wtime()
!$  write(*,*) 'Evaluation took', t2-t1, 'sec. (wallclock time)'

do j=1,numPoints
  X(j)=X(j)/(pmax(j)-pmin(j))
end do

!!! fit replicas
do i=1,numReplica
  write(*,*) '----------------------REPLICA',i,'------------------------'
  !$  t1=OMP_get_wtime()
  call artemide_SetReplica_uTMDPDF(i)
  call artemide_SetReplica_TMDR(i)
  !! setting parameters
  call TMDX_DY_SetProcess(proc)
  call TMDX_DY_XSetup(s,91d0,0d0)
  call TMDX_DY_SetCuts(.true.,ptCut,-etamax,etamax)

  call CalcXsec_DY_PTint_Qint_Yint(Xrep(i,:),pmin,pmax,Qmin,Qmax)
  do j=1,numPoints
  Xrep(i,j)=Xrep(i,j)/(pmax(j)-pmin(j))
end do
  
  !$  t2=OMP_get_wtime()
  !$  write(*,*) 'Evaluation took', t2-t1, 'sec. (wallclock time)'
end do

call artemide_ShowStatistics()
!$  t4=OMP_get_wtime()
!$  write(*,*) 'Total Evaluation took', t4-t3, 'sec. (wallclock time)'

!! caluculating mean value
do i=1,numPoints
  Xmean(i)=sum(Xrep(:,i))/numReplica
end do
write(*,*) 'check:',Xmean-X

!! caluculating deviation
do i=1,numPoints
  Xdev(i)=sqrt(sum(Xrep(:,i)**2)/numReplica-Xmean(i)**2)
end do


open(7, file='RHIC-01-bins', status="replace", action="write")
do i=1,numPoints
  write(7,"(F8.2,F8.2,F12.6,F12.6,F12.6)") pmin(i),pmax(i), X(i), Xmean(i)-Xdev(i)-X(i),Xmean(i)+Xdev(i)-X(i)
end do
CLOSE (7, STATUS='KEEP')

end program example