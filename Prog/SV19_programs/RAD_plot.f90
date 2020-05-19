program example
use aTMDe_control
use TMDR
implicit none

character(*),parameter::prefix='/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/'
character(*),parameter::repFILE=prefix//'Replicas/DY+SIDIS/DY+SIDIS_NNPDF31+DSS_nnlo.rep'
character(*),parameter::constFILE=prefix//'constants-files/DY+SIDIS_nnlo/const-DY+SIDIS_NNPDF31+DSS_nnlo'

real*8,allocatable::NParray(:)
integer::numR,numB,i,j
real*8,allocatable::central(:),mean(:),deviation(:),b(:)
real*8::bMax,step,mu,dd

call artemide_Initialize(constFILE)

numR=artemide_NumOfReplicasInFile(repFILE)

call artemide_GetReplicaFromFile(repFILE,0,NParray)
call artemide_SetNPparameters_TMDR(NParray(1:2))

numB=40
bMax=5d0
!mu=91d0
mu=2.39d0

allocate(b(1:numB))
do i=1,20
  b(i)=0.005d0+i/20d0
end do
do i=21,numB
  b(i)=bMax+(bMax-1d0)*(i-numB)/(numB-20)
end do

allocate(central(1:numB))
allocate(mean(1:numB))
allocate(deviation(1:numB))

do i=1,numB
  central(i)=DNP(mu, b(i),1)
  mean(i)=0d0
  deviation=0d0
end do

do j=1,numR
  call artemide_GetReplicaFromFile(repFILE,j,NParray)
  call artemide_SetNPparameters_TMDR(NParray(1:2))
  do i=1,numB
    dd=DNP(mu, b(i),1)
    mean(i)=mean(i)+dd
    deviation(i)=deviation(i)+dd**2
  end do  
  
end do

 do i=1,numB
  mean(i)=mean(i)/numR
  deviation(i)=Sqrt(deviation(i)/numR-mean(i)**2)
 end do


do i=1,numB
  write(*,"('{',F10.7,',',F10.7,',',F10.7,',',F10.7,'},')") b(i),central(i),mean(i)-deviation(i),mean(i)+deviation(i)
end do


!  


end program example
