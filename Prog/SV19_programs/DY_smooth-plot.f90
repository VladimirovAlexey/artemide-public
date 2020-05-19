!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_control
use TMDX_DY
implicit none

integer,dimension(1:3)::proc
real*8::y,Q,s
real*8::stepqT,maxqT,r
integer::i,maxI
real*8,allocatable::X(:),qT(:),X1(:),X2(:),X3(:),X4(:)

call artemide_Initialize('const-DY_LO',prefix='/home/vla18041/LinkData2/arTeMiDe_Repository/Constants-files/')
! call artemide_Initialize('const-DY_LO',prefix='/home/alexey/artemide_Repository/Constants-files/')!
call artemide_SetNPparameters_TMDR((/1.86041d0, 0.029551d0/))
call artemide_SetNPparameters_uTMDPDF((/0.253434d0, 9.04351d0, 346.999d0, 2.47992d0, -5.69988d0, 0.d0, 0.d0/))

! call artemide_SetNPparameters_TMDR((/2d0,0.001d0/))
! call artemide_SetNPparameters_uTMDPDF((/0.01d0, 0.01d0, 0d0, 0d0, 0d0, 0.d0, 0.d0/))


proc=(/1,1,5/) !!!! pp DY
s=8000d0**2
Q=40.d0
y=2d0

s=200d0**2
Q=6d0
y=1.5d0

call TMDX_DY_setProcess(proc)
call TMDX_DY_XSetup(s,Q,y)
call TMDX_DY_SetCuts(.false.,0d0,-10d0,10d0)!!! no cuts

maxqT=0.3d0*Q
stepqT=0.025d0
maxI=Int(maxqT/stepqT)

allocate(qT(1:maxI))
allocate(X(1:maxI))

do i=1,maxI
  qT(i)=stepqT*i
end do

call CalcXsec_DY(X,qT)

!!!!!!!!!!!!!!!!!!!!!!!!!!! to compute variation band do this
!!! else
! goto 100
allocate(X1(1:maxI))
allocate(X2(1:maxI))
allocate(X3(1:maxI))
allocate(X4(1:maxI))

r=2d0

call artemide_SetScaleVariations(1d0,1d0/r,1d0,1d0)
call CalcXsec_DY(X1,qT)
call artemide_SetScaleVariations(1d0,1d0*r,1d0,1d0)
call CalcXsec_DY(X2,qT)
call artemide_SetScaleVariations(1d0,1d0,1d0,1d0/r)
call CalcXsec_DY(X3,qT)
call artemide_SetScaleVariations(1d0,1d0,1d0,1d0*r)
call CalcXsec_DY(X4,qT)


100 do i=1,maxI
    if(allocated(X1)) then
        write(*,'("{",F6.3,",",F14.12,",",F14.12,",",F14.12,"},")') &
            qT(i),X(i),max(X(i),X1(i),X2(i),X3(i),X4(i)),min(X(i),X1(i),X2(i),X3(i),X4(i))
    else
        write(*,'("{",F6.3,",",F14.12,"},")') qT(i),X(i)
    end if
  !write(*,*) qT(i),X(i)
end do

end program example