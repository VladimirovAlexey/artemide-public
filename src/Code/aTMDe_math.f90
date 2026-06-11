!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 3.03
!
!	The module that contains standard routines for mathematics (polynomials, spec.functions, etc)
!
!				A.Vladimirov (20.01.2026)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module aTMDe_math
use aTMDe_numerics
use aTMDe_interfaces
implicit none

private

public::ChebyshevT,ChebyshevT_array,ChebyshevT_int_array

contains

!!!!!! Return Chebyshev polynomial of order |n| at value x
function ChebyshevT(n,x)
integer,intent(in)::n
real(dp),intent(in)::x
real(dp)::ChebyshevT

real(dp),dimension(0:abs(n))::tt

tt=ChebyshevT_array(n,x)

ChebyshevT=tt(abs(n))

end function ChebyshevT

!!!!!! Return array of Chebyshev polynomial at value x, [T0,T1,T2,...,Tn]
pure function ChebyshevT_array(n,x)
integer,intent(in)::n
real(dp),intent(in)::x
real(dp),dimension(0:abs(n))::ChebyshevT_array

integer::i

ChebyshevT_array(0)=1._dp
if(n>0) ChebyshevT_array(1)=x

do i=2,abs(n)
    ChebyshevT_array(i)=2*x*ChebyshevT_array(i-1)-ChebyshevT_array(i-2)
end do

end function ChebyshevT_array

!!!!!! Return array of the indefinite integral of Chebyshev polynomial at value x, T_k -> (T_{k+1}(x)/(k+1)-T_{k-1}(x)/(k-1))/2
pure function ChebyshevT_int_array(n,x)
integer,intent(in)::n
real(dp),intent(in)::x
real(dp),dimension(0:abs(n))::ChebyshevT_int_array
real(dp),dimension(0:abs(n)+1)::T

integer::i

T(0)=1._dp
ChebyshevT_int_array(0)=x
if(n>0) then
    T(1)=x
    ChebyshevT_int_array(1)=x**2/2
end if
if(n>1) T(2)=2*x**2-1

do i=2,abs(n)
    T(i+1)=2*x*T(i)-T(i-1)
    ChebyshevT_int_array(i)=(T(i+1)/(i+1)-T(i-1)/(i-1))/2
end do

end function ChebyshevT_int_array

end module aTMDe_math
