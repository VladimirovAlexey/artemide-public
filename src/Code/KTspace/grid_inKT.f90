!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 3.0
!
!    This file contains the module, which is common for all TMD-evaluation modules
!
!                A.Vladimirov (16.03.2024)
!---------------------------------------------------------------------------------------

!!!!
!!!! Store the grid in kT-space.
!!!! the grid is 4D (x,kT,Q,f)
!!!! along x and kT, it is Chebyshev, along Q it is Lagrange

!     In the file that uses it add
!module NAME
! INCLUDE this_file
!end module NAME
!

!module grid_inKT

use aTMDe_Numerics
use IO_functions

implicit none

private

character(len=9),parameter :: moduleName="grid-inKT"
character(len=11):: parentModuleName
integer::outputLevel

!!!!---------- Variables about the grid in k-space
!!!! kRanges are the list of values of subgrids for kT
!!!! Nodes are the list of values of nodes in the terms of T=[-1,1]
!!!! NodeFactors are the list of factors (-1)^i/2^%., there %=1 for i=0,NUM, and 0 otherwice
real(dp),allocatable::kRanges(:),kNodes(:),kNodeFactors(:)
!!!! number of Subgrids
integer::numKsubgrids
!!!! number of nodes, it is the same for all subgrids of given class
integer::kGridSize

!!!!---------- Variables about the grid in k-space
!!!! kRanges are the list of values of subgrids for kT
!!!! Nodes are the list of values of nodes in the terms of T=[-1,1]
!!!! NodeFactors are the list of factors (-1)^i/2^%., there %=1 for i=0,NUM, and 0 otherwice
real(dp),allocatable::xRanges(:),xNodes(:),xNodeFactors(:)
!!!! number of Subgrids
integer::numXsubgrids
!!!! number of nodes, it is the same for all subgrids of given class
integer::xGridSize

!!!!---------- Variables about the grid in Q-space
!!!! ordinary qubic interpolation over logarithm scale
!!!! number of nodes in the grid
integer::QGridSize
!!!! step size in the logarithmic space
real(dp)::Qstep
!!!! nodes in Q, to save the time of computation
real(dp),allocatable::QNodes(:)

!!!!---------- General parameters
!!!! number of hadrons
integer::numH
!!!! include Gluon
logical::withGluon
!!!! utmost values of the grids
real(dp)::QMIN,QMAX,kMIN,kMAX,xMIN,lnQMIN,lnQMAX

real(dp)::zero=10.d-12

!!!! Intervals are the lists of (u_{k+1}-u_k)/2 for subgrids
!!!! Means are the lists of (u_{k+1}+u_k)/2 for subgrids
!!!! these are used to speed-up transformation from and to grids (the list are in the trnasformed variables)
real(dp),allocatable::xIntervals(:),xMeans(:),kIntervals(:),kMeans(:)

!!!!--------- MAIN GRID VARIABLE
!!!! (subX,nodeX,subK,nodeK,nodeQ,f,h)
real(dp),allocatable::mainGRID(:,:,:,:,:,:,:)

!!! this is interface for input in grid-function (subgridsinkT,-5:5)
abstract interface
    function function_inKT(x,Q,h,arraySize1,arraySize2)
        import::dp
        integer,intent(in)::arraySize1,arraySize2
        real(dp),dimension(1:arraySize1,0:arraySize2,-5:5) :: function_inKT
        real(dp), intent(in) ::x,Q
        integer,intent(in)::h

    end function function_inKT
end interface

!!! this is interface for input into grid-testing function (-5:5)
abstract interface
    function TMD_atKT(x,kT,Q,h)
        import::dp
        real(dp),dimension(-5:5) :: TMD_atKT
        real(dp), intent(in) ::x,Q,kT
        integer,intent(in)::h

    end function TMD_atKT
end interface

public:: Initialize_GridInKT, PrepareGrid_inKT, ExtractFromGrid_inKT, TestGrid_inKT_internal

contains

!!! initialization
subroutine Initialize_GridInKT(path,moduleLine,gridLine,numH_in,withGluon_in,name,outLevel)
character(len=*)::path
character(len=5),intent(in)::moduleLine,gridLine
integer,intent(in)::outLevel
character(*),intent(in)::name
integer,intent(in)::numH_in
logical,intent(in)::withGluon_in

integer::i

parentModuleName=name
outputLevel=outLevel

!!!! read input about b and kT-spaces
OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")
    !-------------Parameters of grid in X
    call MoveTO(51,moduleLine)
    call MoveTO(51,gridLine)
    call MoveTO(51,'*p3  ')
    read(51,*) numXsubgrids
    allocate(xRanges(0:numXsubgrids))
    call MoveTO(51,'*p4  ')
    read(51,*) xRanges
    call MoveTO(51,'*p5  ')
    read(51,*) xGridSize

    !-------------Parameters of grid in KT
    call MoveTO(51,'*p6  ')
    read(51,*) numKsubgrids
    allocate(kRanges(0:numKsubgrids))
    call MoveTO(51,'*p7  ')
    read(51,*) kRanges
    call MoveTO(51,'*p8  ')
    read(51,*) kGridSize

    !-------------Parameters of grid in Q
    call MoveTO(51,'*p9  ')
    read(51,*) QMIN
    call MoveTO(51,'*p10 ')
    read(51,*) QMAX
    call MoveTO(51,'*p11 ')
    read(51,*) QGridSize
CLOSE (51, STATUS='KEEP')

xMIN=xRanges(0)
kMIN=kRanges(0)
kMAX=kRanges(numKsubgrids)
lnQMIN=log(QMIN)
lnQMAX=log(QMAX)

numH=numH_in

withGluon=withGluon_in

!!!!allocation of lists

allocate(xNodes(0:xGridSize))
allocate(kNodes(0:kGridSize))
allocate(QNodes(0:QGridSize))
allocate(xNodeFactors(0:xGridSize))
allocate(kNodeFactors(0:kGridSize))

allocate(xIntervals(1:numXsubgrids),xMeans(1:numXsubgrids))
allocate(kIntervals(1:numKsubgrids),kMeans(1:numKsubgrids))

!!!!filing the working variables

xIntervals=(log(xRanges(1:numXsubgrids))-log(xRanges(0:numXsubgrids-1)))/2._dp
kIntervals=(log(kRanges(1:numKsubgrids))-log(kRanges(0:numKsubgrids-1)))/2._dp
xMeans=(log(xRanges(1:numXsubgrids))+log(xRanges(0:numXsubgrids-1)))/2._dp
kMeans=(log(kRanges(1:numKsubgrids))+log(kRanges(0:numKsubgrids-1)))/2._dp

xNodeFactors=1._dp

do i=0,xGridSize
  xNodes(i)=cos(i*pi/xGridSize)

  if(i==0 .or. i==xGridSize) then
    xNodeFactors(i)=xNodeFactors(i)/2
  end if
  if(mod(i,2)==1) xNodeFactors(i)=-xNodeFactors(i)

end do

kNodeFactors=1._dp

do i=0,kGridSize
  kNodes(i)=cos(i*pi/kGridSize)

  if(i==0 .or. i==kGridSize) then
    kNodeFactors(i)=kNodeFactors(i)/2
  end if
  if(mod(i,2)==1) kNodeFactors(i)=-kNodeFactors(i)
end do


Qstep=(lnQMAX-lnQMIN)/QGridSize
!!!!! store the values of Q
do i=0,QGridSize
    QNodes(i)=exp(lnQMIN+i*Qstep)
end do

!!!!! allocating our huge grid.
allocate(mainGRID(1:numXsubgrids,0:xGridSize,1:numKsubgrids,0:kGridSize,0:QGridSize,-5:5,1:numH))

end subroutine Initialize_GridInKT

!!!!! values of x computed from the nodes
!!!!! n=subgrid, k=node
pure function XfromNode(n,k)
integer,intent(in)::n,k
real(dp)::XfromNode
XfromNode=exp(xIntervals(n)*xNodes(k)+xMeans(n))
end function XfromNode

!!!!! values of K computed from the nodes
!!!!! n=subgrid, k=node
pure function KfromNode(n,k)
integer,intent(in)::n,k
real(dp)::KfromNode
KfromNode=exp(KIntervals(n)*KNodes(k)+KMeans(n))
end function KfromNode

!!!!! compute and store grid.
!!!!! this function is defined together with Fourier_Levin_array in Fourier_Levin module
!!!!! they exchange the values via the TMD-module.
!!!!! the Gird is TMD*x*kT^2 (factor kT^2 is in the Fourier_Levin_array)
subroutine PrepareGrid_inKT(F)
procedure(function_inKT)::F
real(dp),dimension(1:numKsubgrids,0:kGridSize,-5:5):: receivedValues
integer::h,iQ,iX,jX,ff,iK,jK
real(dp)::time1,time2
!$ real*8::omp_get_wtime

call cpu_time(time1)
!$ time1=omp_get_wtime()
if(outputlevel>1) write(*,*) 'arTeMiDe.',parentModuleName,' starts to compute grid in KT.'
do h=1,numH
do iQ=0,QGridSize
do iX=1,numXsubgrids
do jX=0,xGridSize
    !!!! requesting the values of Fourier
    receivedValues=F(XfromNode(iX,jX),QNodes(iQ),h,numKsubgrids,kGridSize)
    mainGRID(iX,jX,1:numKsubgrids,0:kGridSize,iQ,-5:5,h)=receivedValues*XfromNode(iX,jX)
end do
end do
end do
end do


call cpu_time(time2)
!$ time2=omp_get_wtime()

do h=1,numH
do iQ=0,QGridSize
do iX=1,numXsubgrids
do jX=0,xGridSize
do iK=1,numKsubgrids
do jK=0,kGridSize
do ff=-5,5
    if(ISNAN(mainGRID(iX,jX,iK,jK,iQ,ff,h))) then
        write(*,*) ErrorString("Element of the grid-inKT is NaN. Evaluation STOP",parentModuleName)
        stop
    end if
    if(abs(mainGRID(iX,jX,iK,jK,iQ,ff,h))>1.d8) then
        write(*,*) ErrorString("Element of the grid-inKT is greater than 10^8. Evaluation STOP",parentModuleName)
        write(*,*) "At X=",XfromNode(iX,jX)," kT=",KfromNode(iK,jK)," Q=",QNodes(iQ), "h=",h
        write(*,*) "Value=",mainGRID(iX,jX,iK,jK,iQ,:,h)
        stop
    end if
end do
end do
end do
end do
end do
end do
end do

if(outputlevel>1) then
if(numH>1) then
    write(*,'(" ",A,": Grids in KT are built  (",I3," x",I3,")x(",I3," x",I3,")x",I3,"  calc.time=",F6.2,"s. ")')&
    parentModuleName, numXsubgrids,xGridSize,numKsubgrids,kGridSize,QGridSize, time2-time1
else
    write(*,'(" ",A,": Grid in KT is built  (",I3," x",I3,")x(",I3," x",I3,")x",I3,"  calc.time=",F6.2,"s. ")')&
    parentModuleName, numXsubgrids,xGridSize,numKsubgrids,kGridSize,QGridSize, time2-time1
end if
end if

end subroutine PrepareGrid_inKT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!  INTERPOLATION PART !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! return the value of t=[-1,1] from value of x, and number of subgrid
pure function TfromX(x,n)
integer,intent(in)::n
real(dp),intent(in)::x
real(dp)::TfromX

TfromX=(log(x)-xMeans(n))/xIntervals(n)

end function TfromX

!!!! return the value of t=[-1,1] from value of b, and number of subgrid
pure function TfromK(kT,n)
integer,intent(in)::n
real(dp),intent(in)::kT
real(dp)::TfromK

TfromK=(log(kT)-kMeans(n))/kIntervals(n)
end function TfromK

!!!!! this function interpolates the the grid to the element t
!!!!! it is the baricentric formula
!!!!! f(x)=sum b(i)f(i))/(t-t(i))/sum b(i)/(t-t(i))
function interpolateInX(t,grid)
real(dp),intent(in),dimension(0:xGridSize,-5:5)::grid
real(dp),intent(in)::t
real(dp),dimension(-5:5)::interpolateInX
real(dp),dimension(0:xGridSize)::deltaT
integer::i

deltaT=t-xNodes
!!! pass though each term, if it is close to 0, then the value of fucntion is given by this node
do i=0,xGridSize
    if(abs(deltaT(i))<zero) then
      interpolateInX(-5:5)=grid(i,-5:5)
    return
  end if
end do

deltaT=xNodeFactors/deltaT

interpolateInX=matmul(deltaT,grid)/sum(deltaT)

end function interpolateInX

!!!!! this function interpolates the the grid to the element t
!!!!! just as the previous formula but it operates over list (size of x-nodes) of grids (speed up by factor n)
!!!!! returns the list of interpolation values
function interpolateInK_array(t,grid)
real(dp),intent(in),dimension(0:xGridSize,0:kGridSize,-5:5)::grid
real(dp),intent(in)::t
real(dp),dimension(0:xGridSize,-5:5)::interpolateInK_array
real(dp),dimension(0:kGridSize)::deltaT
integer::i,j

deltaT=t-kNodes
!!! pass though each term, if it is close to 0, then the value of fucntion is given by this node
do i=0,kGridSize
    if(abs(deltaT(i))<zero) then
      interpolateInK_array(0:xGridSize,-5:5)=grid(0:xGridSize,i,-5:5)
    return
  end if
end do

deltaT=kNodeFactors/deltaT
do i=0,xGridSize
do j=-5,5
 interpolateInK_array(i,j)=sum(deltaT(:)*grid(i,:,j))
end do
end do
interpolateInK_array=interpolateInK_array/sum(deltaT)

end function interpolateInK_array

!!!! Interpolation from the grid.
!!!! for values below xMin, above x=1 terminates
!!!! for values below kMin freeze value
!!!! for values above kMax extrapolate with 1/kT^2
function ExtractFromGrid_inKT(x,kT,Q,h)
real(dp),intent(in)::x,kT,Q
integer,intent(in)::h
real(dp),dimension(-5:5)::ExtractFromGrid_inKT

integer::i,nX,nK,nQ
real(dp)::tX,tK

real(dp),dimension(0:xGridSize,-5:5)::interGrid
real(dp),dimension(1:4,-5:5)::interQ
real(dp),dimension(1:4)::deltaTQ

!!! checking exeptions
if(h==0 .or. h>numH) then
write(*,*) ErrorString('the hadron '//numToStr(h)//' is not found in the grid',parentModuleName//"."//moduleName)
write(*,*) 'arTeMiDe: evaluation STOP'
stop
end if

if(x+zero<XMin) then
write(*,*) ErrorString('The TMD with x ='//numToStr(x)//' is called. Current grid size is up to '//&
numToStr(XMin)//'. Enlarge boundaries.',parentModuleName//"."//moduleName)
write(*,*) 'arTeMiDe: evaluation STOP'
stop
end if
if(x>1.d0) then
write(*,*) ErrorString('The TMD with x >1 ('//numToStr(x)//') is called.',parentModuleName//"."//moduleName)
write(*,*) 'arTeMiDe: evaluation STOP'
stop
end if
if(kT<0d0) then
write(*,*) ErrorString('The TMD with kT <0 ('//numToStr(kT)//') is called.',parentModuleName//"."//moduleName)
write(*,*) 'arTeMiDe: evaluation STOP'
stop
end if
if(Q<QMIN) then
write(*,*) ErrorString('The TMD with Q <QMIN ('//numToStr(Q)//'<'//numToStr(QMIN)//') is called.',&
            parentModuleName//"."//moduleName)
write(*,*) 'arTeMiDe: evaluation STOP'
stop
end if
if(Q>QMAX) then
write(*,*) ErrorString('The TMD with Q <QMAX ('//numToStr(Q)//">"//numToStr(QMAX)//') is called.',&
            parentModuleName//"."//moduleName)
write(*,*) 'arTeMiDe: evaluation STOP'
stop
end if

if(x==1.d0) then
ExtractFromGrid_inKT=0._dp
return
end if

!!!! determining Q-node
nQ=int((log(Q)-lnQMIN)/Qstep)
if(nQ==0) then
    nQ=1
else if(nQ==QGridSize) then
    nQ=QGridSize-1
end if
!tQ=((log(Q)-log(QMIN))/Qstep)

!!!! searching for the subgrid in X
nX=0
do i=1,numXsubgrids
    if(x<xRanges(i)) then
        nX=i
        exit
    end if
end do

!!!!! Reminder::  mainGRID(1:numXsubgrids,0:xGridSize,1:numKsubgrids,0:kGridSize,0:QGridSize,-5:5,1:numH)

!!!!! for kT>kTMax extract the frozen value (it will be devided by kT^2 in the end)
if(kT>=kMAX) then
    tX=TfromX(x,nX)
    interQ(1,:)=interpolateInX(tX,mainGRID(nX,0:xGridSize,numKsubgrids,0,nQ-1,-5:5,h))
    interQ(2,:)=interpolateInX(tX,mainGRID(nX,0:xGridSize,numKsubgrids,0,nQ-0,-5:5,h))
    interQ(3,:)=interpolateInX(tX,mainGRID(nX,0:xGridSize,numKsubgrids,0,nQ+1,-5:5,h))
    interQ(4,:)=interpolateInX(tX,mainGRID(nX,0:xGridSize,numKsubgrids,0,nQ+2,-5:5,h))

else
!!!! searching for the subgrid in B
!!!! note that nB=0 implies that 0<b<bMIN
nK=0
do i=0,numKsubgrids
    if(kT<kRanges(i)) then
    nK=i
    exit
    end if
end do

!!!! case if 0<k<kMIN, return value at kMIN
if(nK==0) then
    tX=TfromX(x,nX)

    interQ(1,:)=interpolateInX(tX,mainGRID(nX,0:xGridSize,1,kGridSize,nQ-1,-5:5,h))
    interQ(2,:)=interpolateInX(tX,mainGRID(nX,0:xGridSize,1,kGridSize,nQ-0,-5:5,h))
    interQ(3,:)=interpolateInX(tX,mainGRID(nX,0:xGridSize,1,kGridSize,nQ+1,-5:5,h))
    interQ(4,:)=interpolateInX(tX,mainGRID(nX,0:xGridSize,1,kGridSize,nQ+2,-5:5,h))

else

    tX=TfromX(x,nX)
    tK=TfromK(kT,nK)

    !!! first in b then in x
    interGrid=interpolateInK_array(tK,mainGRID(nX,0:xGridSize,nK,0:kGridSize,nQ-1,-5:5,h))
    interQ(1,:)=interpolateInX(tX,interGrid)

    interGrid=interpolateInK_array(tK,mainGRID(nX,0:xGridSize,nK,0:kGridSize,nQ-0,-5:5,h))
    interQ(2,:)=interpolateInX(tX,interGrid)

    interGrid=interpolateInK_array(tK,mainGRID(nX,0:xGridSize,nK,0:kGridSize,nQ+1,-5:5,h))
    interQ(3,:)=interpolateInX(tX,interGrid)

    interGrid=interpolateInK_array(tK,mainGRID(nX,0:xGridSize,nK,0:kGridSize,nQ+2,-5:5,h))
    interQ(4,:)=interpolateInX(tX,interGrid)
end if
end if

!!!!! barycentric interpolation with Lagrange polynomials

deltaTQ=Q-QNodes(nQ-1:nQ+2)

!!! pass though each term, if it is close to 0, then the value of fucntion is given by this node
if(abs(deltaTQ(1))<zero) then
      ExtractFromGrid_inKT=interQ(1,:)
else if(abs(deltaTQ(2))<zero) then
      ExtractFromGrid_inKT=interQ(2,:)
else if(abs(deltaTQ(3))<zero) then
      ExtractFromGrid_inKT=interQ(3,:)
else if(abs(deltaTQ(4))<zero) then
      ExtractFromGrid_inKT=interQ(4,:)
else
!!!!! interpolation
!!!!! these are barycentric weights of Lagrange interpolation
deltaTQ=(/-1,3,-3,1/)/deltaTQ/6._dp
ExtractFromGrid_inKT=matmul(deltaTQ,interQ)/sum(deltaTQ)
end if

ExtractFromGrid_inKT=ExtractFromGrid_inKT/x/kT**2


!
!
!   do i=-5,5
!     if(ISNAN(ExtractFromGrid_inKT(i))) then
!      write(*,'(" Extracted value (x,kT,Q,f,h) =(",F6.5,", ",F6.2,", ",F6.2,", ",I2,", ",I2,") is computed to NaN")') &
!             x,kT,Q,h,i
!       write(*,*) "==>",ExtractFromGrid_inKT
!       stop
!     end if
!     if(abs(ExtractFromGrid_inKT(i))>1.d12) then
!      write(*,'(" Extracted value (x,kT,Q,f,h) =(",F6.5,", ",F12.8,", ",F6.2,", ",I2,", ",I2,") is computed >10^12")') &
!             x,kT,Q,h,i
!       write(*,*) "==>",ExtractFromGrid_inKT
!       stop
!     end if
!   end do

end function ExtractFromGrid_inKT


subroutine TestGrid_inKT_internal(F)
procedure(TMD_atKT)::F
real(dp),dimension(1:xGridSize)::xTestNodes,xTestValues
real(dp),dimension(1:kGridSize)::kTestNodes,kTestValues
real(dp),dimension(1:QGridSize)::QTestValues
real(dp),dimension(1:xGridSize,1:kGridSize,1:QGridSize,-5:5)::fromGrid,fromExact
integer::i,j,k,iX,jK,kQ,h,fv,totN
real(dp),dimension(1:xGridSize)::res
logical::showLine
real(dp)::avBadness

write(*,*) "---------------------------------INITIATE THE GRID-inKT TEST ------------------------------------"
write(*,*) "                                 ",parentModuleName
write(*,*) "-------------------------------------------------------------------------------------------------"
write(*,*) " The comparison is done by formula R=|X-Y|/(|Y|+0.0001),"
write(*,*) " where X is values from the grid, Y is computed value."
write(*,*) " The computation is done for central points in-between nodes."
write(*,*) " For values R<0.01 the log10(R) is shown in []."
write(*,*) " To save the monito space, only the possibly problematic values are shown."
write(*,*) " For the rest the average value is demonstrated."
write(*,*) "-------------------------------------------------------------------------------------------------"

!!!! locate the test nodes in the middle of the grid.
do i=1,xGridSize
xTestNodes(i)=cos((i-0.5d0)*pi/xGridSize)
end do
do i=1,kGridSize
kTestNodes(i)=cos((i-0.5d0)*pi/kGridSize)
end do
do i=1,QGridSize
QTestValues(i)=exp(lnQMIN+(i-0.5d0)*Qstep)
end do

do h=1,numH


write(*,*) "-------------------------------------------------------------------------------------------------"
write(*,*) color("                         HADRON =  ",c_yellow),h
write(*,*) "-------------------------------------------------------------------------------------------------"
do i=1,numXsubgrids
do j=1,numKsubgrids
do k=1,QGridSize
  write(*,'(A,I3,", ",I3)') "----- Subgrid :",i,j
  write(*,'(A,F10.6," <x< ",F10.6)') "----- ",xRanges(i-1),xRanges(i)
  write(*,'(A,F10.6,"<kT< ",F10.6)') "----- ",kRanges(j-1),kRanges(j)

  xTestValues=exp(xIntervals(i)*xTestNodes+xMeans(i))
  kTestValues=exp(kIntervals(j)*kTestNodes+kMeans(j))

  do iX=1,xGridSize
  do jK=1,kGridSize
  do kQ=1,QGridSize
    fromGrid(iX,jK,kQ,-5:5)=ExtractFromGrid_inKT(xTestValues(iX),kTestValues(jK),QTestValues(kQ),h)
    fromExact(iX,jK,kQ,-5:5)=F(xTestValues(iX),kTestValues(jK),QTestValues(kQ),h)
  end do
  end do
  end do

  do fv=-3,3
  if(fv==0 .and. (.not.withGluon)) cycle !!!! skip gluons if they are not included

  write(*,*) color("                         FLAVOR =  ",c_yellow),fv

  write(*,'(A)',advance='no') "    kT   |   Q \ x ||"
  do iX=xGridSize,2,-1
    write(*,'(F8.4,"  ")',advance='no') xTestValues(iX)
  end do
  write(*,'(F8.4,"  ")') xTestValues(1)
  write(*,'(A)',advance='no') "  -------|"
  do iX=1,xGridSize-1
    write(*,'(A)',advance='no') "----------"
  end do
  write(*,'(A)') "----------"

    totN=1
    avBadness=0._dp
  do jK=kGridSize,2,-1
  do kQ=1,QGridSize

    showLine=.false.
    do iX=xGridSize,1,-1
        res(iX)=abs(fromGrid(iX,jK,kQ,fv)-fromExact(iX,jK,kQ,fv))/(abs(fromExact(iX,jK,kQ,fv))+0.0001d0)
        avBadness=avBadness+res(iX)
        totN=totN+1
        if(res(iX)>0.01) showLine=.true.
    end do

    if(showLine) then
        write(*,'(F8.4," |",F8.2," ||")',advance='no')  kTestValues(jK),QTestValues(kQ)

        do iX=xGridSize,1,-1

        if(res(iX)>0.01) then
            write(*,'(F8.4)',advance='no') res(iX)
        else
            write(*,'("[",F6.2,"]")',advance='no') log10(res(iX))
        end if

        if(res(iX)>1.) then
            write(*,'(A)',advance='no') color("!!",c_red_bold)
        else if(res(iX)>0.1) then
            write(*,'(A)',advance='no') color("!!",c_red)
        else if(res(iX)>0.01) then
            write(*,'(A)',advance='no') color("! ",c_yellow)
        else
            write(*,'(A)',advance='no') "  "
        end if
        end do

        write(*,'(A)') " "
    end if

  end do
  end do

  write(*,'("  AVARAGE BADNESS FOR THIS SECTOR ",F16.12)')avBadness/totN

  end do

end do
end do
end do

end do

end subroutine TestGrid_inKT_internal

!end module grid_inKT
