!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This file contains the definition of 1D Lagrange grid, which is filled upon request !!!
!!!                             8.10.2025                                               !!!
!!!                                 A.Vladimirov                                        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module grid_1D
!use aTMDe_Numerics
implicit none
!!!!! The definition of a grid includes
!!!!! ------- These are defined upon inicialization of the class
!!!!! 1) Size of grid: N
!!!!! 2) Limits of grid: xMIN, xMAX
!!!!! 3) Function FtoGrid(x) : the function to be stored in grid
!!!!! ------- These are defined for each implementation
!!!!! 4) Function factorR(x), that multiplies the input function to make it smoother in grid (and divide on restoration)
!!!!! 5) Function XtoN(x) which transforms input x (with xMIN and xMAX) into [0,N]
!!!!! 6) Function NtoX(n) which transforms input n (with [0,N]) into [xMIN,xMAX]
!!!!! ------- These are static
!!!!! 7) subroutine initialize(xMin,xMax,N,FtoGrid): which creates an empty grid
!!!!! 8) subroutine storePoint(n): it calles FtoGrid for requested point and stores it
!!!!! 9) function getValue(x): it attempts to reconstract the point from the grid, if some node of grid a un-known it computes them

type, abstract::grid1D !!!!! name of abstract class
    integer:: N !!!! size of the grid
    real*8::xMIN,xMAX !!!! limits of the grid
    real*8,allocatable,dimension(:)::mainGrid !!!!! main grid
contains
    !!!!! abstract procedures
    procedure(iF1D), deferred::FtoGrid !!!! function to be stored in grid
    procedure(iF1D), deferred::factorR !!!! multiplies the input function to make it smoother in grid
    procedure(iF1D), deferred::XtoN   !!!! transforms input x (with xMIN and xMAX) into [0,N]
    procedure(iF1D), deferred::NtoX   !!!! transforms input n (with [0,N]) into [xMIN,xMAX]
    !!!!! concrete procedures
    procedure::Initialize
    procedure::storePoint
    procedure::getValue
    procedure::getValueAtNode
end type

abstract interface
    function iF1D(x)
        real*8,intent(out)::iF1D
        real*8,intent(in)::x
    end function iF1D
end interface

contains

!!!!!!! initialize the grid, by allocating it, and settin main constants of the object
subroutine Initialize(self,xMIN_in,xMAX_in,N_in)
class(grid1D), intent(inout)::self
real*8,intent(in)::xMIN_in,xMAX_in
integer,intent(in)::N_in

self % xMIN=xMIN_in
self % xMAX=xMAX_in
self % N=N_in

if(allocated(self % mainGrid)) then deallocate(self % mainGrid)

allocate(self % mainGrid(0:N_in), source=-50.d0)
end subroutine Initialize


!!!!!!! computes the value of function at node, n and store it in the MainGrid
!!!!!!! the stored function is multipliedby factorR(x)
subroutine storePoint(self,n_in)
class(grid1D), intent(inout)::self
integer,intent(in)::n_in

real*8::x

x=NtoX(n_in)
self % mainGrid(n_in)=FtoGrid(x)*factorR(x)

end subroutine storePoint

!!!!!!! Intents to extract the value of function at the node. If it is not setup compute it and stores
function getValueAtNode(self,n_in)
class(grid1D), intent(inout)::self
real*8,intent(out)::getValueAtNode
integer,intent(in)::n_in

!!!! check that function is stored comparing to the default value
if(abs(self % MainGrid(n_in)+50.d0)<0.0000000001d0) then
    call self % storePoint(n_in)
end if
!!!! return value from the grid.
getValueAtNode=self % MainGrid(n_in)
end function getValueAtNode

!!!!!!! Intents to extract the value of function from the grid.
function getValue(self,x_in)
class(grid1D), intent(inout)::self
real*8,intent(out)::getValue
real*8,intent(in)::x_in

real*8::t
integer::t0,t1,t2,t3
real*8::P0,P1,P2,P3, P01,P12,P23, P012, P123, Px

t=self % XtoN(x_in)
!!!!the case of utmost left point
if(t<1) then
    t0=0
    t1=1
    t2=2
    t3=3
else if(t>((self % N)-1)) then !!!! utmost right point
    t3=self % N
    t2=self % N-1
    t1=self % N-2
    t0=self % N-3
else
    t1=int(t)
    t0=t1-1
    t2=t1+1
    t3=t1+2
end if

!!!! at nodes
P0=self % getValueAtNode(t0)
P1=self % getValueAtNode(t1)
P2=self % getValueAtNode(t2)
P3=self % getValueAtNode(t3)

!!! liniar interpolation
P01=(t1-t)*P0+(t-t0)*P1
P12=(t2-t)*P1+(t-t1)*P2
P23=(t3-t)*P2+(t-t2)*P3

!!! cuadratic interpolation
P012=((t2-t)*P01+(t-t0)*P12)/2
P123=((t3-t)*P12+(t-t1)*P23)/2

!!! cubic
Px=((t3-t)*P012+(t-t0)*P123)/3

getValue=Px/factorR(x_in)

end function getValue

end module grid_1D
