!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 3.0
!
!    This file contains the module, which is common for all TMD-evaluation modules
!
!                A.Vladimirov (11.03.2024)
!---------------------------------------------------------------------------------------

!     In the file that uses it add
!module NAME
! INCLUDE this_file
!end module NAME
!

!module Twist2_ChGrid

use aTMDe_Numerics
use IO_functions

implicit none

private

character(len=8),parameter :: moduleName="tw2-grid"
character(len=11):: parentModuleName
integer::outputLevel

!!!! xRanges & bRanges are the list of values of subgrids for x and b
real(dp),allocatable::xRanges(:),bRanges(:)
!!!! number of Subgrids
integer::numXsubgrids,numBsubgrids
!!!! xNodes & bNodes are the list of values of nodes in the terms of T=[-1,1]
real(dp),allocatable::xNodes(:),bNodes(:)
!!!! xNodes & bNodes are the list of factors (-1)^i/2^%., there %=1 for i=0,NUM, and 0 otherwice
real(dp),allocatable::xNodeFactors(:),bNodeFactors(:)
!!!! xGridSize & bGridSize are the number of nodes in the subgrids
!!!! number of nodes is made same for all subgrids, in order to simplify memory operation (store all nodes in single multi-array)
integer::xGridSize,bGridSize
!!!! utmost values of the grids
real(dp)::xMIN,bMIN,bMAX
!!!! number of hadrons
integer::numH
!!!! include Gluon
logical::withGluon

!!! parameter of tolerance
real(dp)::zero=10.d-8

!!!! xIntervals & bIntervals are the list of (u_{k+1}-u_k)/2 for subgrids
!!!! xMeans & bMeanss are the list of (u_{k+1}+u_k)/2 for subgrids
!!!! these are used to speed-up transformation from and to grids (the list are in the trnasformed variables)
real(dp),allocatable::xIntervals(:),bIntervals(:),xMeans(:),bMeans(:)

!!!! the main grid variable
!!!! (subgrid in X,subgrid in B, grid in X, grid in B, flavor(-5:5), hadron)
real(dp),allocatable::gridMain(:,:,:,:,:,:)

public::Twist2_ChGrid_Initialize,Twist2_ChGrid_MakeGrid,ExtractFromGrid,TestGrid

!!! this is interface for twist2-convolution function (-5:5)
abstract interface
    function tw2_convolution(x,b,h,wG)
        import::dp
        real(dp),dimension(-5:5) :: tw2_convolution
        real(dp), intent(in) ::x,b
        integer,intent(in)::h
        logical,intent(in)::wG
    end function tw2_convolution
end interface

contains

!!!! prepare variables


!!!! N
!subroutine Twist2_ChGrid_Initialize(xRanges_in,bRanges_in,xGridSize_in,bGridSize_in,numH_in,withGluon_in,name,outLevel)
subroutine Twist2_ChGrid_Initialize(path,moduleLine,gridLine,numH_in,withGluon_in,name,outLevel)
!real(dp),intent(in)::xRanges_in(:),bRanges_in(:)
!integer,intent(in)::xGridSize_in,bGridSize_in,outLevel,numH_in
character(len=300)::path
character(len=5),intent(in)::moduleLine,gridLine
integer,intent(in)::outLevel,numH_in
logical,intent(in)::withGluon_in
character(*),intent(in)::name

integer::i


parentModuleName=name
outputLevel=outLevel

OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")
    !-------------Parameters of grid
    call MoveTO(51,moduleLine)
    call MoveTO(51,gridLine)
    call MoveTO(51,'*p1  ')
    read(51,*) numXsubgrids
    allocate(xRanges(0:numXsubgrids))
    call MoveTO(51,'*p2  ')
    read(51,*) xRanges
    call MoveTO(51,'*p3  ')
    read(51,*) xGridSize
    call MoveTO(51,'*p4  ')
    read(51,*) numBsubgrids
    allocate(bRanges(0:numBsubgrids))
    call MoveTO(51,'*p5  ')
    read(51,*) bRanges
    call MoveTO(51,'*p6  ')
    read(51,*) bGridSize
CLOSE (51, STATUS='KEEP')


! !!!processing input parameters
xMIN=xRanges(0)
bMIN=bRanges(0)
bMAX=bRanges(numBsubgrids)
numH=numH_in

withGluon=withGluon_in

!!!!allocation of lists

allocate(xNodes(0:xGridSize))
allocate(bNodes(0:bGridSize))
allocate(xNodeFactors(0:xGridSize))
allocate(bNodeFactors(0:bGridSize))

allocate(xIntervals(1:numXsubgrids),xMeans(1:numXsubgrids))
allocate(bIntervals(1:numBsubgrids),bMeans(1:numBsubgrids))

allocate(gridMain(1:numXsubgrids,1:numBsubgrids,0:xGridSize,0:bGridSize,-5:5,1:numH))

!!!!filing the working variables

xIntervals=(log(xRanges(1:numXsubgrids))-log(xRanges(0:numXsubgrids-1)))/2._dp
bIntervals=(log(bRanges(1:numBsubgrids))-log(bRanges(0:numBsubgrids-1)))/2._dp
xMeans=(log(xRanges(1:numXsubgrids))+log(xRanges(0:numXsubgrids-1)))/2._dp
bMeans=(log(bRanges(1:numBsubgrids))+log(bRanges(0:numBsubgrids-1)))/2._dp

xNodeFactors=1._dp

do i=0,xGridSize
  xNodes(i)=cos(i*pi/xGridSize)

  if(i==0 .or. i==xGridSize) then
    xNodeFactors(i)=xNodeFactors(i)/2
  end if
  if(mod(i,2)==1) xNodeFactors(i)=-xNodeFactors(i)

end do

bNodeFactors=1._dp

do i=0,bGridSize
  bNodes(i)=cos(i*pi/bGridSize)

  if(i==0 .or. i==bGridSize) then
    bNodeFactors(i)=bNodeFactors(i)/2
  end if
  if(mod(i,2)==1) bNodeFactors(i)=-bNodeFactors(i)
end do


if(outputLevel>2) then
    write(*,*) 'Grid options:'
    write(*,'(A)',advance="no")  ' |  xSubGrid                 ='
    write(*,'(ES10.3)') xRanges
    write(*,'(A)',advance='no')  ' |  bSubGrid                 ='
    write(*,'(ES10.3)') bRanges
    write(*,'(A,I6,A,I6,A)') ' |  (GridSizeX,GridSizeB)    =(',xGridSize,',',bGridSize,')'
    write(*,'(A,I3)')   ' |  hadrons to grid           =',numH
end if

end subroutine Twist2_ChGrid_Initialize

!!!!! values of x computed from the nodes
!!!!! n=subgrid, k=node
pure function XfromNode(n,k)
integer,intent(in)::n,k
real(dp)::XfromNode
XfromNode=exp(xIntervals(n)*xNodes(k)+xMeans(n))
end function XfromNode

!!!!! values of b computed from the nodes
!!!!! n=subgrid, k=node
pure function BfromNode(n,k)
integer,intent(in)::n,k
real(dp)::BfromNode
BfromNode=exp(bIntervals(n)*bNodes(k)+bMeans(n))
end function BfromNode

!!! this subroutine creates the grid.
!!! F is the function of CxF_compute(x_local,b_local,h,withGluon)
subroutine Twist2_ChGrid_MakeGrid(F)
  procedure(tw2_convolution)::F
  real(dp):: x_local,b_local
  integer:: iX,iB,h,jX,jB,n
  real(dp)::time1,time2
  !$ real*8::omp_get_wtime

  call cpu_time(time1)
  !$ time1=omp_get_wtime()
  if(outputlevel>2) write(*,*) 'arTeMiDe.',parentModuleName,' starts to compute grid.'

!   write(*,*) "1-->",bIntervals
!   write(*,*) "2-->",bMeans
!
!   do iB=1,numBsubgrids
!   do jB=0,bGridSize
!   b_local=BfromNode(iB,jB)
!   write(*,*) iB,jB,b_local
!   end do
!   end do
!   stop

  do h=1,numH
   !$OMP DO
   do iX=1,numXsubgrids
    do iB=1,numBsubgrids
      do jX=0,xGridSize
        do jB=0,bGridSize
          x_local=XfromNode(iX,jX)
          b_local=BfromNode(iB,jB)

          gridMain(iX,iB,jX,jB,-5:5,h)=F(x_local,b_local,h,withGluon)

          do n=-5,5
           if(ISNAN(gridMain(iX,iB,jX,jB,n,h))) then
            write(*,'(" Function value at (x,b,f,h) =(",F6.5,", ",F6.2,", ",I2,", ",I2,") is computed to NaN")') &
            x_local,b_local,h,n
            end if

            if(abs(gridMain(iX,iB,jX,jB,n,h))>10.d6) then
            write(*,'(" Function value at (x,b,f,h) =(",F6.5,", ",F6.2,", ",I2,", ",I2,") is > 10^6")') &
            x_local,b_local,h,n
            write(*,*) "-->",gridMain(iX,iB,jX,jB,-5:5,h)
            end if
          end do

          end do
        end do
      end do
    end do
    !$OMP END DO
    if(outputLevel>1 .and. numH>1) write(*,'(" ",A,": Grid for hadron ",I3," is done")') parentModuleName,h
   end do

  call cpu_time(time2)
  !$ time2=omp_get_wtime()

  if(outputlevel>1) then
    if(numH>1) then
      write(*,'(" ",A,": Grids are built  (",I3," x",I3,")x(",I3," x",I3,")  calc.time=",F6.2,"s. ")')&
        moduleName, numXsubgrids,xGridSize,numBsubgrids,bGridSize, time2-time1
    else
      write(*,'(" ",A,": Grid is built  (",I3," x",I3,")x(",I3," x",I3,")  calc.time=",F6.2,"s. ")')&
        moduleName, numXsubgrids,xGridSize,numBsubgrids,bGridSize, time2-time1
    end if
  end if
    
end subroutine Twist2_ChGrid_MakeGrid


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
!!!!! just as the previous formula but it operates over list (size of b-nodes) of grids (speed up by factor n)
!!!!! returns the list of interpolation values
function interpolateInX_array(t,grid)
real(dp),intent(in),dimension(0:xGridSize,0:bGridSize,-5:5)::grid
real(dp),intent(in)::t
real(dp),dimension(0:bGridSize,-5:5)::interpolateInX_array
real(dp),dimension(0:xGridSize)::deltaT
integer::i,j

deltaT=t-xNodes
!!! pass though each term, if it is close to 0, then the value of fucntion is given by this node
do i=0,xGridSize
    if(abs(deltaT(i))<zero) then
      interpolateInX_array(0:bGridSize,-5:5)=grid(i,0:bGridSize,-5:5)
    return
  end if
end do

deltaT=xNodeFactors/deltaT
do i=0,bGridSize
do j=-5,5
 interpolateInX_array(i,j)=sum(deltaT(:)*grid(:,i,j))
end do
end do
interpolateInX_array=interpolateInX_array/sum(deltaT)

end function interpolateInX_array

!!!!! this function interpolates the the grid to the element t
!!!!! it is the baricentric formula
!!!!! f(x)=sum b(i)f(i))/(t-t(i))/sum b(i)/(t-t(i))
function interpolateInB(t,grid)
real(dp),intent(in),dimension(0:bGridSize,-5:5)::grid
real(dp),intent(in)::t
real(dp),dimension(-5:5)::interpolateInB
real(dp),dimension(0:bGridSize)::deltaT
integer::i

deltaT=t-bNodes
!!! pass though each term, if it is close to 0, then the value of fucntion is given by this node
do i=0,bGridSize
    if(abs(deltaT(i))<zero) then
      interpolateInB(-5:5)=grid(i,-5:5)
    return
  end if
end do

deltaT=bNodeFactors/deltaT
interpolateInB=matmul(deltaT,grid)/sum(deltaT)

end function interpolateInB

!!!!! this function interpolates the the grid to the element t
!!!!! just as the previous formula but it operates over list (size of x-nodes) of grids (speed up by factor n)
!!!!! returns the list of interpolation values
function interpolateInB_array(t,grid)
real(dp),intent(in),dimension(0:xGridSize,0:bGridSize,-5:5)::grid
real(dp),intent(in)::t
real(dp),dimension(0:xGridSize,-5:5)::interpolateInB_array
real(dp),dimension(0:bGridSize)::deltaT
integer::i,j

deltaT=t-bNodes
!!! pass though each term, if it is close to 0, then the value of fucntion is given by this node
do i=0,bGridSize
    if(abs(deltaT(i))<zero) then
      interpolateInB_array(0:xGridSize,-5:5)=grid(0:xGridSize,i,-5:5)
    return
  end if
end do

deltaT=bNodeFactors/deltaT
do i=0,xGridSize
do j=-5,5
 interpolateInB_array(i,j)=sum(deltaT(:)*grid(i,:,j))
end do
end do
interpolateInB_array=interpolateInB_array/sum(deltaT)

end function interpolateInB_array

!!!! return the value of t=[-1,1] from value of x, and number of subgrid
pure function TfromX(x,n)
integer,intent(in)::n
real(dp),intent(in)::x
real(dp)::TfromX

TfromX=(log(x)-xMeans(n))/xIntervals(n)

end function TfromX

!!!! return the value of t=[-1,1] from value of b, and number of subgrid
pure function TfromB(b,n)
integer,intent(in)::n
real(dp),intent(in)::b
real(dp)::TfromB

TfromB=(log(b)-bMeans(n))/bIntervals(n)
end function TfromB

!!!! Interpolation from the grid.
!!!! for values below xMin, above x=1 terminates
!!!! for values below bMin freeze value
!!!! for values above bMax =0
function ExtractFromGrid(x,bT,h)
  real(dp),intent(in)::x,bT
  integer,intent(in)::h
  real(dp),dimension(-5:5)::ExtractFromGrid

  integer::i,nX,nB
  real(dp)::tX,tB
  !real(dp),dimension(0:bGridSize,-5:5)::interGrid
  real(dp),dimension(0:xGridSize,-5:5)::interGrid

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
  if(bT<0d0) then
   write(*,*) ErrorString('The TMD with bT <0 ('//numToStr(bT)//') is called.',parentModuleName//"."//moduleName)
   write(*,*) 'arTeMiDe: evaluation STOP'
   stop
  end if

  if(x==1.d0) then
    ExtractFromGrid=0._dp
    return
  end if

  !!!! searching for the subgrid in X
  nX=0
  do i=1,numXsubgrids
    if(x<xRanges(i)) then
      nX=i
      exit
    end if
  end do

  !!!!! for b>bMAX freeze the value
  if(bT>=BMAX) then
   tX=TfromX(x,nX)
   ExtractFromGrid=interpolateInX(tX,gridMain(nX,numBsubgrids,0:xGridSize,0,-5:5,h))

  else
    !!!! searching for the subgrid in B
    !!!! note that nB=0 implies that 0<b<bMIN
    nB=0
    do i=0,numBsubgrids
      if(bT<bRanges(i)) then
        nB=i
        exit
      end if
    end do

    !!!! case if 0<b<bMIN, return value at bMIN
    if(nB==0) then
      tX=TfromX(x,nX)
      ExtractFromGrid=interpolateInX(tX,gridMain(nX,1,0:xGridSize,bGridSize,-5:5,h))

    else

      tX=TfromX(x,nX)
      tB=TfromB(bT,nB)

      !!! first in b then in x
      interGrid=interpolateInB_array(tB,gridMain(nX,nB,0:xGridSize,0:bGridSize,-5:5,h))
      ExtractFromGrid=interpolateInX(tX,interGrid)
    end if
  end if

!   do i=-5,5
!     if(ISNAN(ExtractFromGrid(i))) then
!      write(*,'(" Extracted value (x,b,f,h) =(",F6.5,", ",F6.2,", ",I2,", ",I2,") is computed to NaN")') &
!             x,bT,h,i
!       write(*,*) "==>",ExtractFromGrid
!       write(*,*) " "
!       write(*,*) "-->",interGrid(:,i)
!       write(*,*) " "
!       write(*,*) "T->",tX,tB
!       write(*,*) " "
!       write(*,*) "X->",tX-xNodes
!       write(*,*) " "
!       write(*,*) "B->",tB-bNodes
!       stop
!     end if
!     if(abs(ExtractFromGrid(i))>1000000.d0) then
!      write(*,'(" Extracted value (x,b,f,h) =(",F6.5,", ",F6.2,", ",I2,", ",I2,") is computed >10^6")') &
!             x,bT,h,i
!       write(*,*) "==>",ExtractFromGrid
!       write(*,*) " "
!       write(*,*) "-->",interGrid(:,i)
!       write(*,*) " "
!       write(*,*) "T->",tX,tB
!       write(*,*) " "
!       write(*,*) "X->",tX-xNodes
!       write(*,*) " "
!       write(*,*) "B->",tB-bNodes
!       stop
!     end if
!   end do

end function ExtractFromGrid

subroutine TestGrid(F)
procedure(tw2_convolution)::F
real(dp),dimension(1:xGridSize)::xTestNodes,xTestValues
real(dp),dimension(1:bGridSize)::bTestNodes,bTestValues
real(dp),dimension(1:xGridSize,1:bGridSize,-5:5)::fromGrid,fromExact
integer::i,j,iX,jB,h,fv
real(dp)::res

write(*,*) "---------------------------------INITIATE THE GRID TEST------------------------------------------"
write(*,*) "                                 ",parentModuleName
write(*,*) "-------------------------------------------------------------------------------------------------"
write(*,*) " The comparison is done by formula R=|X-Y|/(|Y|+0.0001),"
write(*,*) " where X is values from the grid, Y is computed value."
write(*,*) " The computation is done for central points in-between nodes."
write(*,*) " For values R<0.01 the log10(R) is shown in []."
write(*,*) "-------------------------------------------------------------------------------------------------"

!!!! locate the test nodes in the middle of the grid.
do i=1,xGridSize
xTestNodes(i)=cos((i-0.5d0)*pi/xGridSize)
end do
do i=1,bGridSize
bTestNodes(i)=cos((i-0.5d0)*pi/bGridSize)
end do

do h=1,numH


write(*,*) "-------------------------------------------------------------------------------------------------"
write(*,*) color("                         HADRON =  ",c_yellow),h
write(*,*) "-------------------------------------------------------------------------------------------------"
do i=1,numXsubgrids
do j=1,numBsubgrids
  write(*,'(A,I3,", ",I3)') "----- Subgrid :",i,j
  write(*,'(A,F10.6," <x< ",F10.6)') "----- ",xRanges(i-1),xRanges(i)
  write(*,'(A,F10.6," <b< ",F10.6)') "----- ",bRanges(j-1),bRanges(j)

  xTestValues=exp(xIntervals(i)*xTestNodes+xMeans(i))
  bTestValues=exp(bIntervals(j)*bTestNodes+bMeans(j))

  do iX=1,xGridSize
  do jB=1,bGridSize
    fromGrid(iX,jB,-5:5)=ExtractFromGrid(xTestValues(iX),bTestValues(jB),h)
    fromExact(iX,jB,-5:5)=F(xTestValues(iX),bTestValues(jB),h,withGluon)
  end do
  end do

  do fv=-3,3
  if(fv==0 .and. (.not.withGluon)) cycle !!!! skip gluons if they are not included

  write(*,*) color("                         FLAVOR =  ",c_yellow),fv

  write(*,'(A)',advance='no') "  b \ x  |"
  do iX=xGridSize,2,-1
    write(*,'(F8.4,"  ")',advance='no') xTestValues(iX)
  end do
  write(*,'(F8.4,"  ")') xTestValues(1)
  write(*,'(A)',advance='no') "  -------|"
  do iX=1,xGridSize-1
    write(*,'(A)',advance='no') "----------"
  end do
  write(*,'(A)') "----------"

  do jB=bGridSize,2,-1
    write(*,'(F8.4," |")',advance='no')  bTestValues(jB)
    do iX=xGridSize,2,-1
      res=abs(fromGrid(iX,jB,fv)-fromExact(iX,jB,fv))/(abs(fromExact(iX,jB,fv))+0.0001d0)
      if(res>0.01) then
        write(*,'(F8.4)',advance='no') res
      else
        write(*,'("[",F6.2,"]")',advance='no') log10(res)
      end if

      if(res>1.) then
        write(*,'(A)',advance='no') color("!!",c_red_bold)
      else if(res>0.1) then
        write(*,'(A)',advance='no') color("!!",c_red)
      else if(res>0.01) then
        write(*,'(A)',advance='no') color("! ",c_yellow)
      else
        write(*,'(A)',advance='no') "  "
      end if
    end do

    res=abs(fromGrid(1,jB,fv)-fromExact(1,jB,fv))/(abs(fromExact(1,jB,fv))+0.0001d0)
    if(res>0.01) then
      write(*,'(F8.4)') res
    else
      write(*,'("[",F6.2,"]")') log10(res)
    end if

    if(res>1.) then
      write(*,'(A)') color("!!",c_red_bold)
    else if(res>0.1) then
      write(*,'(A)') color("!!",c_red)
    else if(res>0.01) then
      write(*,'(A)') color("! ",c_yellow)
    else
      write(*,'(A)') "  "
    end if

  end do

  end do

end do
end do

end do

end subroutine TestGrid

!end module Twist2_ChGrid
