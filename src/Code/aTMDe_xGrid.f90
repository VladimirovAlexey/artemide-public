!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! This file contains the definition of 4D grid which is filled upon request           !!!
!!!   the grid is used to compute the structure functions is a function of (Q2,qT,x1,x2)!!!
!!! As matter of fact it is a grid of grids. Each 4-cude is a Chebyschev grid           !!!
!!!                             8.10.2025                                               !!!
!!!                                 A.Vladimirov                                        !!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module aTMDe_xGrid
use aTMDe_Numerics
use aTMDe_interfaces
use aTMDe_IO
use EWinput
implicit none

private

!!!!!!!! parameter to check the absence of the value
real(dp),parameter::zero=0.00000000001_dp

!!!!! the number of points inside of Ch-grid
logical::isInitiated=.false.
integer,parameter::ChN=8
real(dp)::ChNodes(0:ChN)  !!!nodes of Ch-grid
real(dp)::ChNodeFactors(0:ChN) !!! node factors for baricentric formula

type:: ChebyschevCube
    logical::isSet=.false.
    integer::N1,N2,N3,N4
    real(dp)::cube(0:ChN,0:ChN,0:ChN,0:ChN)
contains
    procedure::setCube
    procedure::getValue_fromCube
end type

type, public::Xgrid !!!!! name of abstract class
    private
     !!!!!! this parameters needed for generation of messages
    character(:), allocatable::parentName
    integer::outputLevel
    type(Warning_OBJ)::Warning_Handler
    !!!!!! this is function to grid
    procedure(strFUNC), pointer, nopass :: FtoGrid
    !!!!!! the process number
    integer,dimension(1:3)::process0
    !!!!!! parameters of the grid
    !!!!!! size of grids
    integer::N1,N2,N3,N4
    !!!!!! minValues and maxValues of grids
    real(dp)::v1MIN,v2MIN,v3MIN,v4MIN
    real(dp)::v1MAX,v2MAX,v3MAX,v4MAX
    !!!!!! extra parameters pf grids
    real(dp)::p11,p21,p31,p41

    type(ChebyschevCube),allocatable::mainGrid(:,:,:,:)

    !!!! counters
    integer::pointsSet=0,pointsCalled=0
contains
    procedure,public::getValue =>getValue_this
    procedure,public::reset =>reset_this
    procedure,private::storePoint
end type

interface Xgrid
    procedure :: constructor
end interface Xgrid

contains

!!!!!!!!!----------------------------------------------------------------------
!!!!!     transformations for arguments to the space of grid [0,N]
!!!!! the trnaformation must satify v=vMin <=>t=0; v=vMax <=>t=N
!!!!! For the moment I use the same transformation for all directions. Namely logarithmic grid
!!!!! t= N* log(v/vMin)/log(vMax/vMin)
!!!!! So, i defined par1=vMin and par2= N/log(vMax/vMin)
pure function VtoT(v,par1,par2)
real(dp),intent(in)::v,par1,par2
real(dp)::VtoT
VtoT=log(v/par1)*par2
end function VtoT

pure function TtoV(t,par1,par2)
real(dp),intent(in)::t,par1,par2
real(dp)::TtoV
TtoV=par1*exp(t/par2)
end function TtoV

!!!!! transformation for Q direction (it is different because qT can be 0) It is quadratic
!!!!! t= N*sqrt((t-vMin)/(vMin-vMax))
!!!!! So, i defined par1=vMin and par2= N/sqrt(vMax-vMin)
pure function VtoT1(v,par1,par2)
real(dp),intent(in)::v,par1,par2
real(dp)::VtoT1
VtoT1=par2*sqrt(v-par1)
end function VtoT1

pure function TtoV1(t,par1,par2)
real(dp),intent(in)::t,par1,par2
real(dp)::TtoV1
TtoV1=(t/par2)**2+par1
end function TtoV1

!!!!! transformation for qT direction (it is different because qT can be 0)
pure function VtoT2(v,par1,par2)
real(dp),intent(in)::v,par1,par2
real(dp)::VtoT2
!VtoT2=par2*sqrt(v-par1)
VtoT2=par2*(v-par1)
end function VtoT2

pure function TtoV2(t,par1,par2)
real(dp),intent(in)::t,par1,par2
real(dp)::TtoV2
!TtoV2=(t/par2)**2+par1
TtoV2=(t/par2)+par1
end function TtoV2

!!!!!! transformation of t \in [n,n+1] to u \in [-1,1]
pure function TtoU(t)
real(dp),intent(in)::t
real(dp)::TtoU
TtoU=2*(int(t)-t)+1
end function TtoU
!!!!!! transfomration back
pure function UtoT(u,n)
real(dp),intent(in)::u
integer,intent(in)::n
real(dp)::UtoT
UtoT=(1+2*n-u)*0.5_dp
end function UtoT

!!!!!! transformation of t \in [n,n+1] to u \in [-1,1] and directly to the beta/(u-u_j)
pure function TtoW(t)
real(dp),intent(in)::t
real(dp),dimension(0:ChN)::TtoW

integer::i

TtoW=TtoU(t)-ChNodes
do i=0,ChN
    if(Abs(TtoW(i))<zero) then
        TtoW=0._dp
        TtoW(i)=1._dp
        return
    end if
end do
TtoW=ChNodeFactors/TtoW
end function TtoW

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Functions for the 4D cube !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
subroutine setCube(this,origin,n1,n2,n3,n4)
class(ChebyschevCube),intent(inout)::this
type(Xgrid),intent(in)::origin
integer,intent(in)::n1,n2,n3,n4
integer::i1,i2,i3,i4
real(dp)::v1,v2,v3,v4

this%N1=n1
this%N2=n2
this%N3=n3
this%N4=n4

do i1=0,ChN
do i2=0,ChN
do i3=0,ChN
do i4=0,ChN
    v1=TtoV1(UtoT(ChNodes(i1),n1),origin%v1MIN,origin%p11)
    v2=TtoV2(UtoT(ChNodes(i2),n2),origin%v2MIN,origin%p21)
    v3=TtoV(UtoT(ChNodes(i3),n3),origin%v3MIN,origin%p31)
    v4=TtoV(UtoT(ChNodes(i4),n4),origin%v4MIN,origin%p41)
    this%cube(i1,i2,i3,i4)=origin%FtoGrid(v1,v2,v3,v4,origin%process0)*weigthFunction(v1,v2,v3,v4,origin%process0)
end do
end do
end do
end do
this%isSet=.true.

end subroutine setCube

!!!!!!! interpolation by baricentric formula
function getValue_fromCube(this,t1,t2,t3,t4)
class(ChebyschevCube),intent(in)::this
real(dp), intent(in)::t1,t2,t3,t4
real(dp)::getValue_fromCube

real(dp),dimension(0:ChN)::delta1,delta2,delta3,delta4
integer::i1,i2,i3,i4


delta1=TtoW(t1)
delta2=TtoW(t2)
delta3=TtoW(t3)
delta4=TtoW(t4)

getValue_fromCube=0._dp
do i1=0,ChN
do i2=0,ChN
do i3=0,ChN
do i4=0,ChN
    getValue_fromCube=getValue_fromCube+this%cube(i1,i2,i3,i4)*delta1(i1)*delta2(i2)*delta3(i3)*delta4(i4)
end do
end do
end do
end do

getValue_fromCube=getValue_fromCube/sum(delta1)/sum(delta2)/sum(delta3)/sum(delta4)

end function getValue_fromCube

!!!!!!!!!----------------------------------------------------------------------
!!!!!! Weight function
!!!!!! this function specifies the weight depending on the process
!!!!!! such that grid saves F*W, and this should be a smooth function
function weigthFunction(Q2,qT,x1,x2,process0)
real(dp)::weigthFunction
real(dp),intent(in)::Q2,qT,x1,x2
integer,dimension(1:3),intent(in)::process0

SELECT CASE(process0(3))
CASE(2,3)!!!!!!!! Z-boson peak
    weigthFunction=x1*x2*(qT+1)**2*((Q2-MZ2)**2+GammaZ2*MZ2)/(Q2*Q2)
CASE(4,5,6)!!!!!!!! W-boson peak
    weigthFunction=x1*x2*(qT+1)**2*((Q2-MW2)**2+GammaW2*MW2)/(Q2*Q2)
CASE DEFAULT
    weigthFunction=x1*x2*(qT+1)**2
end select

end function weigthFunction

!!!!! constructor for the Xgrid.
!!!!! one should specify process0, min,max,N for each variable
function constructor(F,proc_in,v1min,v1max,v1N,v2min,v2max,v2N,v3min,v3max,v3N,v4min,v4max,v4N,name,outLevel) result(this)
type(Xgrid)::this
procedure(strFUNC)::F
integer,dimension(1:3),intent(in)::proc_in
integer,intent(in)::v1N,v2N,v3N,v4N
real(dp),intent(in)::v1Min,v1Max,v2Min,v2Max,v3Min,v3Max,v4Min,v4Max
character(*),intent(in)::name
integer::outLevel
integer::i

this%parentName=name
this%outputLevel=outLevel
if(this%outputLevel>1) write(*,*) trim(this%parentName)//": Initializing the structure-function grid for process ("// &
numToStr(proc_in(1))//","//int4ToStr(proc_in(2))//","//int4ToStr(proc_in(3))//")"

this%Warning_Handler=Warning_OBJ(moduleName=this%parentName,messageCounter=0,messageTrigger=3)

!!!! the common parameter should be computed
if(.not.isInitiated) then
ChNodeFactors=1._dp
do i=0,ChN
    ChNodes(i)=cos(i*pi/ChN)

  if(i==0 .or. i==ChN) then
    ChNodeFactors(i)=ChNodeFactors(i)/2
  end if
  if(mod(i,2)==1) ChNodeFactors(i)=-ChNodeFactors(i)
end do
end if

this%pointsSet=0
this%pointsCalled=0

this%FtoGrid=>F
this%process0=proc_in

this%v1MIN=v1MIN
this%v1MAX=v1MAX
this%N1=v1N

this%v2MIN=v2MIN
this%v2MAX=v2MAX
this%N2=v2N

this%v3MIN=v3MIN
this%v3MAX=v3MAX
this%N3=v3N

this%v4MIN=v4MIN
this%v4MAX=v4MAX
this%N4=v4N

!this%p11=v1N/Log(v1Max/v1Min)
this%p11=v1N/sqrt(v1Max-v1Min)
!this%p21=v2N/Log(v2Max/v2Min)
!this%p21=v2N/sqrt(v2Max-v2Min)
this%p21=v2N/(v2Max-v2Min)
this%p31=v3N/Log(v3Max/v3Min)
this%p41=v4N/Log(v4Max/v4Min)

allocate(this%mainGrid(0:this%N1,0:this%N2,0:this%N3,0:this%N3))

end function constructor

!!!!!!! reset the grid such that it will be computed for anew
subroutine reset_this(this)
class(Xgrid), intent(inout)::this

call this%Warning_Handler%reset()

if(allocated(this%mainGrid)) deallocate(this%mainGrid)
allocate(this%mainGrid(0:this%N1,0:this%N2,0:this%N3,0:this%N3))

if(this%outputLevel>2) write(*,*) trim(this%parentName)//": grid reset for process ("// &
numToStr(this%process0(1))//","//int4ToStr(this%process0(2))//","//int4ToStr(this%process0(3))//")"
if(this%outputLevel>2) write(*,*) trim(this%parentName)//": nodes used   : "//numToStr(this%pointsSet)
if(this%outputLevel>2) write(*,*) trim(this%parentName)//": nodes called : "//numToStr(this%pointsCalled)

this%pointsSet=0
this%pointsCalled=0

end subroutine reset_this

!!!!!!! computes the value of function at node, n and store it in the MainGrid
!!!!!!! the stored function is multipliedby weigthFunction(x)
subroutine storePoint(this,n1,n2,n3,n4)
class(Xgrid), intent(inout)::this
integer,intent(in)::n1,n2,n3,n4

this%pointsSet=this%pointsSet+4*(ChN+1)
call this%mainGrid(n1,n2,n3,n4)%setCube(this,n1,n2,n3,n4)

end subroutine storePoint

!!!!!!! Computes the interpolation from the grid
!!!!!!! using the Neuvile algorithm for qubic Legandre interpolation
!!!!!!! it is assumed that distance between node is =1
function getValue_this(this,v1,v2,v3,v4)
class(Xgrid), intent(inout)::this
real(dp)::getValue_this
real(dp),intent(in)::v1,v2,v3,v4

real(dp)::t1,t2,t3,t4
integer::n1,n2,n3,n4
real(dp)::subGrid4D(0:3,0:3,0:3,0:3),subGrid3D(0:3,0:3,0:3),subGrid2D(0:3,0:3),subGrid1D(0:3)

this%pointsCalled=this%pointsCalled+1

!!!!! check tha point inside of the grid
if(v1<this%v1MIN .or. v2<this%v2MIN .or. v3<this%v3MIN .or. v4<this%v4MIN &
    .or. v1>this%v1MAX .or. v2>this%v2MAX .or. v3>this%v3MAX .or. v4>this%v4MAX) then
    if(v1<this%v1MIN .or. v1>this%v1MAX) call this%Warning_Handler%WarningRaise("Point Q2 outside of specified grid !")
    if(v2<this%v2MIN .or. v2>this%v2MAX) call this%Warning_Handler%WarningRaise("Point qT outside of specified grid !")
    if(v3<this%v3MIN .or. v3>this%v3MAX) call this%Warning_Handler%WarningRaise("Point x1 outside of specified grid !")
    if(v4<this%v4MIN .or. v4>this%v4MAX) call this%Warning_Handler%WarningRaise("Point x2 outside of specified grid !")


    if(this%outputLevel>2) then
        write(*,*) "  Requested point:", v1,v2,v3,v4
        write(*,*) "  Min value      :",this%v1MIN,this%v2MIN,this%v3MIN,this%v4MIN
        write(*,*) "  Max value      :",this%v1MAX,this%v2MAX,this%v3MAX,this%v4MAX
    end if
    getValue_this=this%FtoGrid(v1,v2,v3,v4,this%process0)
end if

!write(*,*) "ask for:",v1,v2,v3,v4

!!!!! getting the values of interpolator variable
t1=VtoT1(v1,this%v1Min,this%p11)
t2=VtoT2(v2,this%v2Min,this%p21)
t3=VtoT(v3,this%v3Min,this%p31)
t4=VtoT(v4,this%v4Min,this%p41)

n1=int(t1)
n2=int(t2)
n3=int(t3)
n4=int(t4)
!!!! check that function is stored comparing to the default value
if(.not.(this%MainGrid(n1,n2,n3,n4)%isSet)) call this%storePoint(n1,n2,n3,n4)


!!!! return value from the grid.
getValue_this=this%MainGrid(n1,n2,n3,n4)%getValue_fromCube(t1,t2,t3,t4)

getValue_this=getValue_this/weigthFunction(v1,v2,v3,v4,this%process0)
end function getValue_this


end module aTMDe_xGrid
