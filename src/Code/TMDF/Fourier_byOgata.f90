!            arTeMiDe 3.0
!
!    This file contains the module, which realises the Fourier integration by Ogata quadrature
!    The core of the module is inhereted from artemide v.2.06
!
!                A.Vladimirov (15.03.2024)
!---------------------------------------------------------------------------------------

!     In the file that uses it add
!module NAME
! INCLUDE this_file
!end module NAME
!

!module Fourier_byOgata

use aTMDe_Numerics
use IO_functions

implicit none

private


!------------------------------------------Tables-----------------------------------------------------------------------
integer,parameter::Nmax=1000
INCLUDE 'Tables/BesselZero1000.f90'

character(len=8),parameter :: moduleName="Ogata"
character(len=11):: parentModuleName
integer::outputLevel

!!!!! I split the qT over runs qT<qTSegmentationBoundary
!!!!! In each segment I have the ogata quadrature with h=hOGATA*hSegmentationWeight
!!!!! It helps to convergen integrals, since h(optimal) ~ qT
integer,parameter::hSegmentationNumber=6
real(dp),dimension(1:hSegmentationNumber),parameter::hSegmentationWeight=(/0.001d0,0.01d0,0.1d0,1d0,2.d0,5d0/)
real(dp),dimension(1:hSegmentationNumber),parameter::qTSegmentationBoundary=(/0.001d0,0.01d0,0.1d0,1.d0,10.d0,50.d0/)

real(dp)::hOGATA,tolerance
!!!weights of ogata quadrature
real(dp),dimension(1:hSegmentationNumber,0:3,1:Nmax)::ww
!!!nodes of ogata quadrature
real(dp),dimension(1:hSegmentationNumber,0:3,1:Nmax)::bb

public::PrepareTables,Fourier_byOgata

!!! this is interface for a function of b
abstract interface
    function W_inB(b)
        import::dp
        real(dp):: W_inB
        real(dp), intent(in) ::b
    end function W_inB
end interface


contains


!!!Prepare tables for Ogata quadrature with given h
subroutine PrepareTables(tolerance_in,hOGATA_in)
real(dp),intent(in)::tolerance_in,hOGATA_in
integer::i,k,j
real(dp)::hS!=h*hSegmentationWeight
real(dp)::xi

if(outputLevel>1) write(*,*) parentModuleName,': preparing Ogata tables'

tolerance=tolerance_in
hOGATA=hOGATA_in

do j=1,hSegmentationNumber
do k=0,3
do i=1,Nmax

  hS=hOGATA*hSegmentationWeight(j)
  xi=JZero(k,i)

  !     ww(j,k,i)=BESSEL_JN(k,bb(j,k,i))/xi/(BESSEL_JN(k+1,xi)**2)&
  !         *(pi*xi*hS*Cosh(xi*hS)+Sinh(pi*Sinh(xi*hS)))/(1d0+Cosh(pi*Sinh(xi*hS)))

  !!! if we too far away in xI*hS, the double exponential grow rapidly.
  !!! and for >6, it generates term 10^{300} and exceed the presision

  if(xi*hS>6.d0) then
      bb(j,k,i)=xi*Tanh(piHalf*Sinh(xi*hS))
      ww(j,k,i)=BESSEL_JN(k,bb(j,k,i))/xi/(BESSEL_JN(k+1,xi)**2)

  else
      bb(j,k,i)=xi*Tanh(piHalf*Sinh(xi*hS))
      ww(j,k,i)=BESSEL_JN(k,bb(j,k,i))/xi/(BESSEL_JN(k+1,xi)**2)&
      *(pi*xi*hS*Cosh(xi*hS)/(2d0*Cosh(piHalf*Sinh(xi*hS))**2)+Tanh(piHalf*Sinh(xi*hS)))
  end if

end do
end do
end do


if(outputLevel>2) write(*,'(A,I4)') ' | Maximum number of nodes    :',Nmax
if(outputLevel>1) write(*,*) parentModuleName,': Ogata tables prepared'

end subroutine PrepareTables
 
!!!This is the defining module function
!!! It evaluates the integral
!!!  int_0^infty   b^(n+1) db/2  Jn(b qT) F(b)
!!!
!!! the parameter n is passed in n, function F is F(b)
subroutine Fourier_byOgata(n,F,qT,res,ISconvergent)
procedure(W_inB)::F
integer,intent(in)::n
real(dp),intent(in)::qT
real(dp),intent(out)::res
logical,intent(out)::ISconvergent

real(dp)::integral,eps,delta
real(dp)::v1,v2,v3,v4
integer::k,j,Nsegment

integral=0.d0
ISconvergent=.true.

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
  eps=ww(Nsegment,n,k)*(bb(Nsegment,n,k)**(n+1))*F(bb(Nsegment,n,k)/qT)

  v4=v3
  v3=v2
  v2=v1
  v1=ABS(eps)

  delta=(v1+v2+v3+v4)
  integral=integral+eps

  !!! here we check that residual term is smaller than already collected integral
  !!! also checking the zerothness of the integral. If already collected integral is null it is null
  !!! Here is potential bug. If the first 10 points give zero (whereas some later points do not), the integral will be zero
  if((delta<tolerance*abs(integral) .or. abs(integral)<1d-32) .and. k>=10) exit
end do

!!!!! if we run out of maximum number of nodes.
if(k>=Nmax) then
  if(outputlevel>0) WRITE(*,*) WarningString('OGATA quadrature diverge. TMD decaing too slow? ',parentModuleName)
  if(outputlevel>1) then
    write(*,*) 'Information over the last call ----------'
    write(*,*) 'bt/qT= ',bb(Nsegment,n,Nmax)/qT, 'qT=',qT, '| segmentation zone=',Nsegment,&
        ' ogata h=',hOGATA*hSegmentationWeight(Nsegment)
    write(*,*) 'W=',F(bb(Nsegment,n,Nmax)/qT), 'eps/integral =', eps/integral
    write(*,*) 'residual term=',delta, '>',tolerance
    write(*,*) '------------------------------------------'
  end if
  ISconvergent=.false.
end if

!!! result is scaled by qT [because the argument of Bessel was scaled bqT-> B]
res=integral/(qT**(n+2))

end subroutine Fourier_byOgata

!end module Fourier_byOgata
