!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot at givem s,Q,y, in qT (no bin-integrations)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program example
use aTMDe_control
use TMDX_DY
implicit none

character(*),parameter::prefix='/home/vla18041/LinkData2/arTeMiDe_Repository/artemide/Models/SV19/'
character(*),parameter::repFILE=prefix//'Replicas/DY+SIDIS/DY+SIDIS_NNPDF31+DSS_nnlo.rep'
character(*),parameter::constFILE=prefix//'constants-files/DY+SIDIS_nnlo/const-DY+SIDIS_NNPDF31+DSS_nnlo'
real*8,allocatable::NParray(:)
integer::numR

integer,dimension(1:3)::proc
real*8::s,yMin,yMax,Qmin,Qmax
real*8::stepqT,maxqT,r
integer::i,maxI
real*8,allocatable::X0(:),Xplus(:),Xminus(:),qTmin(:),qTmax(:)


call artemide_Initialize(constFILE)

numR=artemide_NumOfReplicasInFile(repFILE)

call artemide_GetReplicaFromFile(repFILE,0,NParray)
!call artemide_SetNPparameters_TMDR(NParray(1:2))
call artemide_SetNPparameters(NParray)



proc=(/1,1,5/) !!!! pp DY
s=7000d0**2
Qmin=20.d0
Qmax=40.d0
ymin=2.d0
ymax=4.5d0

call TMDX_DY_setProcess(proc)
call TMDX_DY_XSetup(s,(Qmax+Qmin)/2d0,(ymax+ymin)/2d0)
call TMDX_DY_SetCuts(.true.,15d0,2.d0,4.5d0)!!! no cuts

maxqT=0.3d0*(Qmax+Qmin)/2d0
stepqT=0.25d0
maxI=Int(maxqT/stepqT)

allocate(qTmin(1:maxI))
allocate(qTmax(1:maxI))

do i=1,maxI
  qTmin(i)=stepqT*(i-1)
  qTmax(i)=stepqT*i
end do

!call ComputeXsecWithStatErr(qTmin,qTmax,X0,Xminus,Xplus)
call ComputeXsecWithScaleErr(qTmin,qTmax,X0,Xminus,Xplus)


write(*,*) '---------------------------S=7TeV---------------------------------------'
do i=1,maxI
        write(*,'("{",F6.3,",",F6.3,",",F14.12,",",F14.12,",",F14.12,"},")') qTmin(i),qTmax(i),X0(i),Xminus(i),Xplus(i)
end do


s=8000d0**2
call TMDX_DY_XSetup(s,(Qmax+Qmin)/2d0,(ymax+ymin)/2d0)

!call ComputeXsecWithStatErr(qTmin,qTmax,X0,Xminus,Xplus)
call ComputeXsecWithScaleErr(qTmin,qTmax,X0,Xminus,Xplus)
write(*,*) '---------------------------S=8TeV---------------------------------------'
do i=1,maxI
        write(*,'("{",F6.3,",",F6.3,",",F14.12,",",F14.12,",",F14.12,"},")') qTmin(i),qTmax(i),X0(i),Xminus(i),Xplus(i)
end do

contains 

!!! compute statistical error due to replica distribution
subroutine ComputeXsecWithStatErr(qTmin_in,qTmax_in,Xout,XminusOut,XplusOut)
  real*8,intent(in)::qTmin_in(:),qTmax_in(:)
  real*8,allocatable,intent(out)::XOut(:),XminusOut(:),XplusOut(:)
  real*8,allocatable::mean(:),deviation(:),X(:),Err(:)
  
  integer::l,i,j
  
  l=size(qTmin_in)
  allocate(XOut(1:l))
  allocate(XminusOut(1:l))
  allocate(XplusOut(1:l))
  allocate(mean(1:l))
  allocate(deviation(1:l))
  allocate(X(1:l))
  allocate(Err(1:l))
  
  do i=1,l
    mean(i)=0d0
    deviation(i)=0d0
  end do
  
  do j=1,numR
    call artemide_GetReplicaFromFile(repFILE,j,NParray)
    call artemide_SetNPparameters(NParray)
    
    call CalcXsec_DY_PTint_Qint_Yint(X,qTmin_in,qTmax_in,Qmin,Qmax,yMin,yMax)
    
    mean=mean+X
    deviation=deviation+X**2
    
  end do

  do i=1,l
    Xout(i)=mean(i)/numR
    Err(i)=Sqrt(deviation(i)/numR - Xout(i)**2)
    XminusOut(i)=Xout(i)-Err(i)
    XplusOut(i)=Xout(i)+Err(i)
  end do
  
end subroutine ComputeXsecWithStatErr

!!! compute scle-variation error
subroutine ComputeXsecWithScaleErr(qTmin_in,qTmax_in,XOut,XminusOut,XplusOut)
  real*8,intent(in)::qTmin_in(:),qTmax_in(:)
  real*8,allocatable,intent(out)::XOut(:),XminusOut(:),XplusOut(:)
  real*8,allocatable::X(:),X1(:),X2(:),X3(:),X4(:)
  
  integer::l,i,j
  
  l=size(qTmin_in)
  allocate(XOut(1:l))
  allocate(XminusOut(1:l))
  allocate(XplusOut(1:l))
  
  allocate(X(1:l))
  allocate(X1(1:l))
  allocate(X2(1:l))
  allocate(X3(1:l))
  allocate(X4(1:l))
  
  call artemide_SetScaleVariations(1d0,1d0,1d0,1d0)
  call CalcXsec_DY_PTint_Qint_Yint(X,qTmin_in,qTmax_in,Qmin,Qmax,yMin,yMax)

  r=2d0
  
  call artemide_SetScaleVariations(1d0,1d0/r,1d0,1d0)
  call CalcXsec_DY_PTint_Qint_Yint(X1,qTmin_in,qTmax_in,Qmin,Qmax,yMin,yMax)
  call artemide_SetScaleVariations(1d0,1d0*r,1d0,1d0)
  call CalcXsec_DY_PTint_Qint_Yint(X2,qTmin_in,qTmax_in,Qmin,Qmax,yMin,yMax)
  call artemide_SetScaleVariations(1d0,1d0,1d0,1d0/r)
  call CalcXsec_DY_PTint_Qint_Yint(X3,qTmin_in,qTmax_in,Qmin,Qmax,yMin,yMax)
  call artemide_SetScaleVariations(1d0,1d0,1d0,1d0*r)
  call CalcXsec_DY_PTint_Qint_Yint(X4,qTmin_in,qTmax_in,Qmin,Qmax,yMin,yMax)
  
  do i=1,l
    Xout(i)=X(i)
    XminusOut(i)=min(X(i),X1(i),X2(i),X3(i),X4(i))
    XplusOut(i)=max(X(i),X1(i),X2(i),X3(i),X4(i))
  end do

  
end subroutine ComputeXsecWithScaleErr

end program example