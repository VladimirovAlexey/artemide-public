!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 3.03
!
!	Contains the module that realizes the Fourier integration with Ogata quadrature
!
!				A.Vladimirov (14.10.2025)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module aTMDe_Ogata
use aTMDe_Numerics
use aTMDe_interfaces
use aTMDe_IO
implicit none

private

!------------------------------------------Tables of zeros of Bessel------------------------------------------------------------
integer,parameter::Nmax=1000
INCLUDE 'Tables/BesselZero1000.f90'

!!!!! I split the qT over runs qT<qTSegmentationBoundary
!!!!! In each segment I have the ogata quadrature with h=hOGATA*hSegmentationWeight
!!!!! It helps to convergen integrals, since h(optimal) ~ qT
integer,parameter::hSegmentationNumber=7
real(dp),dimension(1:hSegmentationNumber),parameter::hSegmentationWeight=(/0.001d0,0.01d0,0.1d0,1d0,2d0,5d0,10d0/)
real(dp),dimension(1:hSegmentationNumber),parameter::qTSegmentationBoundary=(/0.001d0,0.01d0,0.1d0,1d0,10d0,50d0,200d0/)

!!!!!! The object of Ogata integration --------
type, public :: OgataIntegrator
  !!!------ specified upon construction
  !!!!!! this parameters needed for generation of messages
  character(:), allocatable::parentName
  integer::outputLevel

  !!!------ specified by Prepare tables
  !!!!!! the order of the transform assigned to this object (it is used to defined trnsformation of TMDs only)
  integer::N
  !!!!!! general mass scale
  real(dp)::TMDmass,transformationFactor
  !!!!!! minimal kT == the value of kT below which the transfomation is frozen
  real(dp):: kTmin
  !!!!!! numerical parameters of the evaluation
  real(dp)::hOGATA,tolerance
  !!!weights of ogata quadrature
  real(dp),dimension(1:hSegmentationNumber,0:3,1:Nmax)::ww
  !!!nodes of ogata quadrature
  real(dp),dimension(1:hSegmentationNumber,0:3,1:Nmax)::bb

contains
  !!!!! this is ordinary hankel transformation 1D->1D (with 1/2 b^{n+1}J_n)
  procedure,public:: Transform => Transform_def
  !!!!! this is general hankel transformation (-5:5)->(-5:5) (with 1/2 b^n J_k)
  !procedure,public:: Moment_G => G_def
  !!!!! The transformation of TMD to kT-space
  !!!!! this is G[N+1,N] with extra factor M^{2num}/num!/qT^num.
  procedure,public:: TransformTMD => TransformTMD_def

  !!!! The moment G=G_{n,n} in [2402.01836]
  procedure,public:: Moment_G => Moment_G_def
  !!!! The moment X=G_{n+1,n} in [2402.01836]
  procedure,public:: Moment_X => Moment_X_def
end type

interface OgataIntegrator
    procedure :: constructor
end interface OgataIntegrator

contains

!!!!! constructor for the OgataIntegrator. It automatically prepare the tables with given h
!!!!! IMPORTANT NOTE: the tables are generated for trnaformation with J_n/2 (this 1/2 is for historical reasons)
function constructor(parentName,outputLevel,order_in, tolerance_in,hOGATA_in,TMDmass_in,kTmin) result(this)
type(OgataIntegrator)::this
character(len=*),intent(in)::parentName
integer,intent(in)::outputLevel,order_in
real(dp),intent(in)::tolerance_in,hOGATA_in,TMDmass_in,kTmin

integer::i,j,k
real(dp)::xi,hS!=h*hSegmentationWeight

this%parentName=parentName
this%outputLevel=outputLevel
if(this%outputLevel>2) write(*,*) this%parentName,': preparing Ogata tables'

this%TMDmass=TMDmass_in
this%N=order_in
this%tolerance=tolerance_in
this%hOGATA=hOGATA_in
this%kTmin=kTmin

SELECT CASE(this%N)
CASE(0)
    this%transformationFactor=1._dp
CASE(1)
    this%transformationFactor=this%TMDmass**2
CASE(2)
    this%transformationFactor=this%TMDmass**4/2
CASE(3)
    this%transformationFactor=this%TMDmass**6/6
CASE DEFAULT
    write(*,*) ErrorString('aTMDe_Ogata: Ogata transformation with N>3 is not defined. ',this%parentName),this%N
    error stop
END SELECT

do j=1,hSegmentationNumber
do k=0,3
do i=1,Nmax

  hS=this%hOGATA*hSegmentationWeight(j)
  xi=JZero(k,i)

  !!! if we too far away in xI*hS, the double exponential grow rapidly.
  !!! and for >6, it generates term 10^{300} and exceed the presision

  if(xi*hS>6.d0) then
      this%bb(j,k,i)=xi*Tanh(piHalf*Sinh(xi*hS))
      this%ww(j,k,i)=BESSEL_JN(k,this%bb(j,k,i))/xi/(BESSEL_JN(k+1,xi)**2)

  else
      this%bb(j,k,i)=xi*Tanh(piHalf*Sinh(xi*hS))
      this%ww(j,k,i)=BESSEL_JN(k,this%bb(j,k,i))/xi/(BESSEL_JN(k+1,xi)**2)&
      *(pi*xi*hS*Cosh(xi*hS)/(2d0*Cosh(piHalf*Sinh(xi*hS))**2)+Tanh(piHalf*Sinh(xi*hS)))
  end if

end do
end do
end do

if(this%outputLevel>2) write(*,'(A,I4)') ' | Maximum number of nodes    :',Nmax
if(this%outputLevel>1) write(*,*) this%parentName,': Ogata tables (transformation order ',int4ToStr(this%N),') are prepared'
if(outputLevel>2) write(*,*) color("If you use TMMs, please, cite [2402.01836]",c_cyan)
end function constructor

!!!This is the defining module function
!!! It evaluates the integral
!!!  int_0^infty   b^(n+1) db/2  Jn(b qT) F(b)
!!!
!!! or order of hankel transform is defined by this%N,
!!! the function to integrate is "func_1D"
subroutine Transform_def(this,F,qT_in,n,res,ISconvergent)
class(OgataIntegrator), intent(inout)::this
procedure(func_1D)::F
real(dp),intent(in)::qT_in
integer,intent(in)::n
real(dp),intent(out)::res
logical,intent(out)::ISconvergent

real(dp)::integral,eps,delta,qT
real(dp)::v1,v2,v3,v4
integer::k,j,Nsegment

integral=0.d0
ISconvergent=.true.

if(qT_in<this%kTmin) then
    qT=this%kTmin
else
    qT=qT_in
end if

v1=1d0
v2=1d0
v3=1d0
v4=1d0

!!! define segment of qT
do j=1,hSegmentationNumber
  if(qT<qTSegmentationBoundary(j)) exit
end do
if(j>hSegmentationNumber) then
  Nsegment=hSegmentationNumber
else
  Nsegment=j
end if

!!! sum over OGATA nodes
do k=1,Nmax!!! maximum of number of bessel roots preevaluated in the head
  eps=this%ww(Nsegment,n,k)*(this%bb(Nsegment,n,k)**(n+1))*F(this%bb(Nsegment,n,k)/qT)

  v4=v3
  v3=v2
  v2=v1
  v1=ABS(eps)

  delta=(v1+v2+v3+v4)
  integral=integral+eps

  !!! here we check that residual term is smaller than already collected integral
  !!! also checking the zerothness of the integral. If already collected integral is null it is null
  !!! Here is potential bug. If the first 10 points give zero (whereas some later points do not), the integral will be zero
  if((delta<this%tolerance*abs(integral) .or. abs(integral)<1d-32) .and. k>=10) exit
end do

!!!!! if we run out of maximum number of nodes.
if(k>=Nmax) then
  if(this%outputlevel>0) WRITE(*,*) WarningString('OGATA quadrature diverge. W decaing too slow? ',this%parentName)
  if(this%outputlevel>1) then
    write(*,*) 'Information over the last call ----------'
    write(*,*) 'bt/qT= ',this%bb(Nsegment,n,Nmax)/qT, 'qT=',qT, '| segmentation zone=',Nsegment,&
        ' ogata h=',this%hOGATA*hSegmentationWeight(Nsegment)
    write(*,*) 'W=',F(this%bb(Nsegment,n,Nmax)/qT), 'eps/integral =', eps/integral
    write(*,*) 'residual term=',delta, '>',this%tolerance
    write(*,*) '------------------------------------------'
  end if
  ISconvergent=.false.
end if

!!! result is scaled by qT [because the argument of Bessel was scaled bqT-> B]
res=integral/(qT**(n+2))
end subroutine Transform_def

!!!This is the defining module function
!!! It evaluates the integral
!!! int_0^infty   db/(2pi) b^n  J_k(b qT) F1
!!! note the factor 1/(2pi). 1/2 of it is included into ww and 1/pi is included here
!!!
!!! the function to integrate is "func_1D_array5" i.e. (-5:5) function
function G_def(this,F,qT,n,k)
class(OgataIntegrator), intent(inout)::this
procedure(func_1D_array5)::F
real(dp),intent(in)::qT
integer,intent(in)::n,k
real(dp),dimension(-5:5)::G_def

real(dp)::integral(-5:5),eps(-5:5)
real(dp)::v1(-5:5),v2(-5:5),v3(-5:5),v4(-5:5),delta(-5:5)
logical:: partDone(-5:5)
integer::r,j,Nsegment

integral=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)

v1=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
v2=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
v3=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
v4=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
partDone=(/.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false./)

!!! define segment of qT
do j=1,hSegmentationNumber
    if(qT<qTSegmentationBoundary(j)) exit
end do
if(j>hSegmentationNumber) then
    Nsegment=hSegmentationNumber
else
    Nsegment=j
end if

do r=1,Nmax!!! maximum of number of bessel roots preevaluated in the head
    eps=this%ww(Nsegment,k,r)*(this%bb(Nsegment,k,r)**n)*F(this%bb(Nsegment,k,r)/qT)

    v4=v3
    v3=v2
    v2=v1
    v1=abs(eps)

    delta=(v1+v2+v3+v4)
    integral=integral+eps

    !!! here we check that residual term is smaller than already collected integral
    !!! also checking the zerothness of the integral. If already collected integral is null it is null
    !!! Here is potential bug. If the first 10 points give zero (whereas some later points do not), the integral will be zero
    !!! I check for each separate flavor
    do j=-5,5
        if((delta(j)<this%tolerance*ABS(integral(j)) .or. ABS(integral(j))<1d-32) .and. r>=10) partDone(j)=.true.
    end do
    if(partDone(-5).and.partDone(-4).and.partDone(-3).and.partDone(-2).and.partDone(-1)&
        .and.partDone(0).and.partDone(1).and.partDone(2).and.partDone(3).and.partDone(4).and.partDone(5)) exit

end do

if(r>=Nmax) then
    if(this%outputlevel>0) write(*,*) WarningString('OGATA quadrature diverge. TMD decaing too slow? ',this%parentName)
        if(this%outputlevel>2) then
        write(*,*) 'Information over the last call ----------'
        write(*,*) partDone
        write(*,*) 'bt/qT= ',this%bb(Nsegment,k,Nmax)/qT, 'qT=',qT, '| segmentation zone=',Nsegment,&
            ' ogata h=',this%hOGATA*hSegmentationWeight(Nsegment)
        write(*,*) 'W=',F(this%bb(Nsegment,k,r)/qT), 'eps/integral =', eps/integral
        write(*,*) 'v1+v2+v3+v4=',v1+v2+v3+v4, '>',this%tolerance*(ABS(integral(1))+ABS(integral(2)))
        write(*,*) '------------------------------------------'
        end if
end if
!!! result is scaled by qT [because the argument of Bessel was scaled bqT-> B]
G_def=integral/(qT**(n+1)*pi)

end function G_def

!!! It evaluates the integral
!!!  int_0^infty   b db/2pi  J_num(b qT) F1  (b/qT)^num M^{2num}/num! = G[n+1,n] (M^{2num}/num!/pi)/qT^num
!!!
!!! or order of hankel transform is defined by this%N,
!!! the function to integrate is "func_1D_array5" i.e. (-5:5) function
function TransformTMD_def(this,F,qT)
class(OgataIntegrator), intent(inout)::this
procedure(func_1D_array5)::F
real(dp),intent(in)::qT
real(dp),dimension(-5:5)::TransformTMD_def

TransformTMD_def=G_def(this,F,qT,this%N+1,this%N)*this%transformationFactor/(qT)**this%N

end function TransformTMD_def

!!!! The moment G=G_{n,n} in [2402.01836], defined as
!!!! mu^{n+1}/2^n/n! int b^n J_{n+1} (b mu) F
!!!! note: there is no factor (2pi)
function Moment_G_def(this,F,mu)
class(OgataIntegrator), intent(inout)::this
procedure(func_1D_array5)::F
real(dp),intent(in)::mu
real(dp),dimension(-5:5)::Moment_G_def(-5:5)

SELECT CASE(this%N)
    CASE(0)
        Moment_G_def=pix2*mu*G_def(this,F,mu,0,1)
    CASE(1)
        Moment_G_def=pix2*mu**2/2*G_def(this,F,mu,1,2)
    CASE(2)
        Moment_G_def=pix2*mu**3/8*G_def(this,F,mu,2,3)
    CASE DEFAULT
        error stop ErrorString("MOMENT G is defined only for n=0,1,2",this%parentName)
END SELECT

end function Moment_G_def

!!!! The moment X=G_{n+1,n} in [2402.01836], defined as
!!!! mu^{n+3}/2^n/n!/(2M^2) int b^n ((n+1)J_{n+1} (b mu)-J_{n+3} (b mu))/(n+2)
!!!! =
!!!! mu^{n+3}/2^n/n!/(2M^2) int b^n (J_{n+1} (b mu)-J_{n+2}(b mu)/(b mu))
!!!! I use the last formula because it utilizes lower order of Bessel-function
!!!! note: there is no factor (2pi)
function Moment_X_def(this,F,mu)
class(OgataIntegrator), intent(inout)::this
procedure(func_1D_array5)::F
real(dp),intent(in)::mu
real(dp),dimension(-5:5)::Moment_X_def(-5:5)

SELECT CASE(this%N)
    CASE(0)
        Moment_X_def=pix2*mu**2/(2*this%TMDmass**2)*(mu*G_def(this,F,mu,0,1)-2*G_def(this,F,mu,-1,2))
        !Moment_X=pix2*mu**3/(2*2*this%TMDmass**2)*(G_def(this,F,mu,0,1)-G_def(this,F,mu,0,3))
    CASE(1)
        Moment_X_def=pix2*mu**3/(4*this%TMDmass**2)*(mu*G_def(this,F,mu,1,2)-2*G_def(this,F,mu,0,3))
    CASE DEFAULT
        error stop ErrorString("MOMENT X is defined only for n=0,1",this%parentName)
END SELECT

end function Moment_X_def

end module aTMDe_Ogata
