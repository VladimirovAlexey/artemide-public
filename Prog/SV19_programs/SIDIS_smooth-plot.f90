!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at given s,Q,x,z, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_control
use TMDX_SIDIS
implicit none

integer,dimension(1:3)::proc
real*8::x,z,Q,s
real*8::stepqT,maxqT
integer::i,maxI
real*8,allocatable::Xsec(:),qT(:)

call artemide_Initialize('const-SIDIS_LO',prefix='/home/vla18041/LinkData2/arTeMiDe_Repository/Constants-files/')

call artemide_SetNPparameters_TMDR((/1.86041d0, 0.029551d0/))
call artemide_SetNPparameters_uTMDPDF((/0.253434d0, 9.04351d0, 346.999d0, 2.47992d0, -5.69988d0, 0.d0, 0.d0/))
call artemide_SetNPparameters_uTMDFF((/0.264d0,0.479d0,0.459d0,0.539d0/))

proc=(/1,1,2001/) !!!! p->pi+ SIDIS
s=300d0
Q=4d0
x=0.1d0
z=0.4d0

call TMDX_SIDIS_setProcess(proc)
call TMDX_SIDIS_XSetup(s,z,x,Q,0.938d0,0.d0)
call TMDX_SIDIS_SetCuts(.false.,0d0,-10d0,10d0,0d0)!!! no cuts

maxqT=0.4d0*Q*z
stepqT=maxqT/25d0
maxI=Int(maxqT/stepqT)

allocate(qT(1:maxI))
allocate(Xsec(1:maxI))

do i=1,maxI
  qT(i)=stepqT*i
end do

call CalcXsec_SIDIS(Xsec,qT)

do i=1,maxI
  write(*,'("{",F6.3,",",F12.8,"},")') qT(i),Xsec(i)/1000d0
  !write(*,*) qT(i),X(i)
end do

end program example