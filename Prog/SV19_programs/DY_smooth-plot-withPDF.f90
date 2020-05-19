!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_control
use uTMDPDF
use TMDX_DY
implicit none

integer::numR=100
integer,dimension(1:3)::proc
real*8::y,Q,s
real*8::stepqT,maxqT
integer::i,j,maxI,k
real*8,allocatable::X(:),qT(:),mean(:),deviation(:)
real*8 :: a1,a2,a3,a4,a5,a6,a7,a8,a9

call artemide_Initialize('const-DY_TEST',prefix='/home/vla18041/LinkData2/arTeMiDe_Repository/Constants-files/')

call artemide_SetNPparameters_TMDR((/1.86041d0, 0.029551d0/))
call artemide_SetNPparameters_uTMDPDF((/0.253434d0, 9.04351d0, 346.999d0, 2.47992d0, -5.69988d0, 0.d0, 0.d0/))

proc=(/1,1,5/) !!!! pp DY


! s=25d0**2
! Q=8d0
! y=0.2d0
! stepqT=0.1d0

s=8000d0**2
Q=91.d0
y=0d0
stepqT=0.25d0

maxqT=0.25d0*Q

maxI=Int(maxqT/stepqT)


call TMDX_DY_setProcess(proc)
call TMDX_DY_XSetup(s,Q,y)
call TMDX_DY_SetCuts(.false.,0d0,-10d0,10d0)!!! no cuts



allocate(qT(1:maxI))
allocate(X(1:maxI))
allocate(deviation(1:maxI))
allocate(mean(1:maxI))

do i=1,maxI
  qT(i)=stepqT*i
end do

call CalcXsec_DY(X,qT)

write(*,*) '-------------------------Central------------------------'

do i=1,maxI
  write(*,'("{",F6.3,",",F12.8,"},")') qT(i),X(i)
  !write(*,*) qT(i),X(i)
end do
!goto 10
write(*,*) '-------------------------PDF distr------------------------'

do i=1,maxI        
    mean(i)=0d0
    deviation(i)=0d0
  end do

OPEN(UNIT=51, FILE="/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/FIGURES_DY+SIDIS_2019/PDF-replicas/NNPDF_dist.txt",&
      ACTION="read", STATUS="old")
 !! replicas
 do j=1,numR
  read(51,*) k,a1,a2,a3,a4,a5,a6,a7,a8,a9
  call uTMDPDF_SetPDFReplica(k)
  call uTMDPDF_SetLambdaNP((/a3,a4,a5,a6,a7,a8,a9/),.false.,.false.)
  
  call CalcXsec_DY(X,qT)
  
  do i=1,maxI        
    mean(i)=mean(i)+X(i)
    deviation(i)=deviation(i)+X(i)**2
  end do
 end do
 CLOSE (51, STATUS='KEEP')
 
 do i=1,maxI
  mean(i)=mean(i)/numR
  deviation(i)=Sqrt(deviation(i)/numR-mean(i)**2)
 end do
 
 do i=1,maxI
  write(*,"('{',F10.7,',',F10.7,',',F10.7,',',F10.7,'},')") qT(i),mean(i),mean(i)-deviation(i),mean(i)+deviation(i)
 end do
 
10 write(*,*) '-------------------------PDF distr+------------------------'
 do i=1,maxI        
    mean(i)=0d0
    deviation(i)=0d0
  end do

OPEN(UNIT=51, FILE="/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/FIGURES_DY+SIDIS_2019/PDF-replicas/NNPDF_dist+.txt",&
      ACTION="read", STATUS="old")
 !! replicas
 numR=46
 do j=1,numR
  read(51,*) k,a1,a2,a3,a4,a5,a6,a7,a8,a9
  call uTMDPDF_SetPDFReplica(k)
  call uTMDPDF_SetLambdaNP((/a3,a4,a5,a6,a7,a8,a9/),.false.,.false.)
  
  call CalcXsec_DY(X,qT)
  
  do i=1,maxI        
    mean(i)=mean(i)+X(i)
    deviation(i)=deviation(i)+X(i)**2
  end do
 end do
 CLOSE (51, STATUS='KEEP')
 
 do i=1,maxI
  mean(i)=mean(i)/numR
  deviation(i)=Sqrt(deviation(i)/numR-mean(i)**2)
 end do
 
 do i=1,maxI
  write(*,"('{',F10.7,',',F10.7,',',F10.7,',',F10.7,'},')") qT(i),mean(i),mean(i)-deviation(i),mean(i)+deviation(i)
 end do
 
 write(*,*) '-------------------------fNP distir------------------------'
  do i=1,maxI        
    mean(i)=0d0
    deviation(i)=0d0
  end do
 
OPEN(UNIT=51, FILE="/home/vla18041/LinkData2/WorkingFiles/TMD/Fit_Notes/FIGURES_DY+SIDIS_2019/PDF-replicas/NNPDF_stat.txt",&
      ACTION="read", STATUS="old")
 !! replicas
 numR=100
 do j=1,numR
  read(51,*) k,a1,a2,a3,a4,a5,a6,a7,a8,a9
  call uTMDPDF_SetPDFReplica(k)
  call uTMDPDF_SetLambdaNP((/a3,a4,a5,a6,a7,a8,a9/),.false.,.false.)
  
  call CalcXsec_DY(X,qT)
  
  do i=1,maxI        
    mean(i)=mean(i)+X(i)
    deviation(i)=deviation(i)+X(i)**2
  end do
 end do
 CLOSE (51, STATUS='KEEP')
 
 do i=1,maxI
  mean(i)=mean(i)/numR
  deviation(i)=Sqrt(deviation(i)/numR-mean(i)**2)
 end do
 
 do i=1,maxI
  write(*,"('{',F10.7,',',F10.7,',',F10.7,',',F10.7,'},')") qT(i),mean(i),mean(i)-deviation(i),mean(i)+deviation(i)
 end do

end program example