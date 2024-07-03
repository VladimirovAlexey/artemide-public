!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 3.0
!
!    This file contains the module, which is common for all TMD-evaluation modules
!
!                A.Vladimirov (16.03.2024)
!---------------------------------------------------------------------------------------

!!!!
!!!! This module performes the Fourier tranfomation to kT-space.
!!!! It does it by the Levin method, and transform the Chebyshev grid to the Chebyshev grid.
!!!! It precomputes the variable bValues (which is a list of b's)
!!!! And the TranformationArray. Applying trnaformation array to the f(bValues), one get a grid in kT-space

!     In the file that uses it add
!module NAME
! INCLUDE this_file
!end module NAME
!

!module Fourier_Levin

use aTMDe_Numerics
use IO_functions
use InverseMatrix, only: Inverse

implicit none

private

character(len=7),parameter :: moduleName="Fourier"
character(len=11):: parentModuleName
integer::outputLevel

!!!! General parameters
real(dp)::TMDmass=1._dp         !! mass parameter used as mass-scale
integer::TMDtypeN !!!!! this is the order of Bessel-transform (IT IS STRICT FOR TMD)

!!!!---------- Variables about the grid in b-space
!!!! bValues is the list of b-values from all subgrids, collected in the single list
!!!! It is used for the transform, and has the shape (1:lengthOfbValues)
real(dp),allocatable::bValues(:)
!!!! the length of bValues is =(numBsubgrids*(bGridSize+1))
integer::lengthOfbValues
!!!! bRanges are the list of values of subgrids for b
!!!! Nodes are the list of values of nodes in the terms of T=[-1,1]
!!!! NodeFactors are the list of factors (-1)^i/2^%., there %=1 for i=0,NUM, and 0 otherwice
real(dp),allocatable::bRanges(:),bNodes(:),bNodeFactors(:)
!!!! number of Subgrids
integer::numBsubgrids
!!!! number of nodes, it is the same for all subgrids of given class
integer::bGridSize
!!!! Intervals are the lists of (u_{k+1}-u_k)/2 for subgrids
!!!! Means are the lists of (u_{k+1}+u_k)/2 for subgrids
!!!! these are used to speed-up transformation from and to grids (the list are in the trnasformed variables)
real(dp),allocatable::bIntervals(:),bMeans(:)


!!!!---------- Variables about the grid in k-space
!!!! kRanges are the list of values of subgrids for kT
!!!! Nodes are the list of values of nodes in the terms of T=[-1,1]
!!!! NodeFactors are the list of factors (-1)^i/2^%., there %=1 for i=0,NUM, and 0 otherwice
real(dp),allocatable::kRanges(:),kNodes(:),kNodeFactors(:)
!!!! number of Subgrids
integer::numKsubgrids
!!!! number of nodes, it is the same for all subgrids of given class
integer::kGridSize
!!!! Intervals are the lists of (u_{k+1}-u_k)/2 for subgrids
!!!! Means are the lists of (u_{k+1}+u_k)/2 for subgrids
!!!! these are used to speed-up transformation from and to grids (the list are in the trnasformed variables)
real(dp),allocatable::kIntervals(:),kMeans(:)

!!!!--------- Variables concerning trnaformation to the kT-space
!!!! The matrix or arrays (for each node in kT) to transform f(b)
real(dp),allocatable::TranformationArray(:,:,:)


real(dp),allocatable::CCweight(:)!!!!! weights of CC-quadrature
real(dp)::zero=10.d-12
!!!! first zero of bessel function J0,J1,...,J5
real(dp),dimension(0:5),parameter::besselZERO=(/&
2.4048255576957727686216318793265_dp, 3.8317059702075123156144358863082_dp, &
5.1356223018406825563014016901378_dp, 6.3801618959239835062366146419427_dp, &
7.5883424345038043850696300079857_dp, 8.7714838159599540191228671334096_dp/)

public:: Initialize_Fourier_Levin,Fourier_Levin,Fourier_Levin_array,TestFourier

!!! this is interface for abstract TMD function (-5:5)
abstract interface
    function TMD_distribution(b)
        import::dp
        real(dp),dimension(-5:5) :: TMD_distribution
        real(dp), intent(in) ::b
    end function TMD_distribution
end interface

contains


!!!!!! Initialization of parameters from the const-file
!!! path: path to INI-file
!!! moduleLine: "*12  " the line for the module input
!!! gridLine: "*F  " the line for the KT-grid input within the module section
!!! name: name of the parent module (for messages)
!!! outLevel: output Level (for messages)
!!! TMDt: is the type of KT-trnasform for TMD (n=0, b J_0;  n=1 b^2 J_1)
subroutine Initialize_Fourier_Levin(path,moduleLine,gridLine,name,outLevel,TMDt)
character(len=*)::path
character(len=5),intent(in)::moduleLine,gridLine
integer,intent(in)::outLevel,TMDt
character(*),intent(in)::name


logical::prepareGrid
integer::i,j,k

parentModuleName=name
outputLevel=outLevel

TMDtypeN=TMDt

!!!! read input about b and kT-spaces
OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")
    !------------- General parameters
    call MoveTO(51,'*0   ')
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p2  ')
    read(51,*) TMDmass
    !-------------Parameters of grid in kT
    call MoveTO(51,moduleLine)
    call MoveTO(51,gridLine)
    call MoveTO(51,'*p1  ')
    read(51,*) prepareGrid
    call MoveTO(51,'*p6  ')
    read(51,*) numKsubgrids
    allocate(kRanges(0:numKsubgrids))
    call MoveTO(51,'*p7  ')
    read(51,*) kRanges
    call MoveTO(51,'*p8  ')
    read(51,*) kGridSize

    !-------------Parameters of grid in b
    call MoveTO(51,'*p12 ')
    read(51,*) numBsubgrids
    allocate(bRanges(0:numBsubgrids))
    call MoveTO(51,'*p13 ')
    read(51,*) bRanges
    call MoveTO(51,'*p14 ')
    read(51,*) bGridSize
CLOSE (51, STATUS='KEEP')


!!!!allocation of lists
allocate(kNodes(0:kGridSize))
allocate(kNodeFactors(0:kGridSize))
allocate(kIntervals(1:numKsubgrids),kMeans(1:numKsubgrids))

allocate(bNodes(0:bGridSize))
allocate(bNodeFactors(0:bGridSize))
allocate(bIntervals(1:numBsubgrids),bMeans(1:numBsubgrids))

allocate(CCweight(0:bGridSize))


!!!!filing the working variables for K-grid
kIntervals=(log(kRanges(1:numKsubgrids))-log(kRanges(0:numKsubgrids-1)))/2._dp
kMeans=(log(kRanges(1:numKsubgrids))+log(kRanges(0:numKsubgrids-1)))/2._dp

kNodeFactors=1._dp

do i=0,kGridSize
  kNodes(i)=cos(i*pi/kGridSize)

  if(i==0 .or. i==kGridSize) then
    kNodeFactors(i)=kNodeFactors(i)/2
  end if
  if(mod(i,2)==1) kNodeFactors(i)=-kNodeFactors(i)
end do


!!!!!!! now for the B-space
numBsubgrids=size(bRanges)-1

!!!!! filling working variables for b-grid
bIntervals=(log(bRanges(1:numBsubgrids))-log(bRanges(0:numBsubgrids-1)))/2._dp
bMeans=(log(bRanges(1:numBsubgrids))+log(bRanges(0:numBsubgrids-1)))/2._dp

bNodeFactors=1._dp

do i=0,bGridSize
  bNodes(i)=cos(i*pi/bGridSize)

  if(i==0 .or. i==bGridSize) then
    bNodeFactors(i)=bNodeFactors(i)/2
  end if
  if(mod(i,2)==1) bNodeFactors(i)=-bNodeFactors(i)
end do

!!!!!! make the list of required values of b
lengthOfbValues=numBsubgrids*(bGridSize+1)
allocate(bValues(1:lengthOfbValues))
k=1
do i=1,numBsubgrids
do j=0,bGridSize
    bValues(k)=BfromNode(i,j)
    k=k+1
end do
end do

!!!! precompute CCweights
call CCweight_compute()

if(prepareGrid) then
    allocate(TranformationArray(1:numKsubgrids,0:kGridSize,1:lengthOfbValues))
    call PrepareTransformationMatrix()
    if(outLevel>=2) write(*,*) "Initialization of tables for Hankel transform over grid is done for "//trim(parentModuleName)
end if

end subroutine Initialize_Fourier_Levin

!!!!! values of K computed from the nodes
!!!!! n=subgrid, k=node
pure function KfromNode(n,k)
integer,intent(in)::n,k
real(dp)::KfromNode
KfromNode=exp(KIntervals(n)*KNodes(k)+KMeans(n))
end function KfromNode

!!!!! values of b computed from the nodes
!!!!! n=subgrid, k=node
pure function BfromNode(n,k)
integer,intent(in)::n,k
real(dp)::BfromNode
BfromNode=exp(bIntervals(n)*bNodes(k)+bMeans(n))
end function BfromNode

!!!!! This subroutine creates the transformation matrix from b to kT space
!!!!! the transformation is
!!!!! for TMDtypeN=0:  kT^2\int_0^\infty db/(2pi) J_0(b kT) [b F(b)]
!!!!! for TMDtypeN=1:  kT^2 M^2/kT\int_0^\infty db/(2pi) J_1(b kT) [b^2 F(b)]
!!!!! It works over the grid
subroutine PrepareTransformationMatrix()
integer::i,j

    !!!!! actually computation of the arrays for Fourier transform
SELECT CASE(TMDtypeN)
    CASE(0) !!!! uTMDPDF, and similar
        do i=1,numKsubgrids
        do j=0,kGridSize
            TranformationArray(i,j,:)=FourierArrayAtQ(KfromNode(i,j))/pix2*KfromNode(i,j)**2
        end do
        end do
    CASE(1) !!!! SiversTMDPDF, and similar
        do i=1,numKsubgrids
        do j=0,kGridSize
            TranformationArray(i,j,:)=FourierArrayAtQ(KfromNode(i,j))/pix2*(TMDMass**2)*KfromNode(i,j)
        end do
        end do
    CASE DEFAULT
        write(*,*) ErrorString("Fourier_Levin: Unknown TMDtype. Presently impleneted only types 0,1",parentModuleName)
        stop

END SELECT

end subroutine PrepareTransformationMatrix

!!!!!! make the Fourier of the function F with the interface (-5:5)
!!!!!! at the point kT (in GeV)
!!!!!! the transformation is
!!!!!! for TMDtypeN=0:  \int_0^\infty db/(2pi) b J_0(b q) F(b)
!!!!!! for TMDtypeN=1:  M^2/kT\int_0^\infty db/(2pi) b^2 J_1(b q) F(b)
function Fourier_Levin(F,kT)
real(dp),intent(in)::kT
real(dp),dimension(-5:5)::Fourier_Levin
procedure(TMD_distribution)::F
integer::i
real(dp),dimension(1:lengthOfbValues)::Levin_Array

Levin_Array=FourierArrayAtQ(kT)

SELECT CASE(TMDtypeN)
    CASE(0) !!! uTMDPDF, and similar
        Fourier_Levin=0._dp
        do i=1,lengthOfbValues
        Fourier_Levin=Fourier_Levin+Levin_Array(i)*F(bValues(i))*bValues(i)
        end do
        Fourier_Levin=Fourier_Levin/pix2

    CASE(1) !!! SiversTMDPDF, and similar
        Fourier_Levin=0._dp
        do i=1,lengthOfbValues
        Fourier_Levin=Fourier_Levin+Levin_Array(i)*F(bValues(i))*bValues(i)*bValues(i)
        end do
        Fourier_Levin=Fourier_Levin/pix2*TMDMass*TMDMass/kT
    CASE DEFAULT
        write(*,*) ErrorString("Fourier_Levin: Unknown TMDtype. Presently impleneted only types 0,1",parentModuleName)
        stop

END SELECT

end function Fourier_Levin

!!!!!! make the Fourier of the function F with the interface (-5:5)
!!!!!! at the modes of kT-grid. Returns the list (1:numKsubgrids, 0:kGridSize,-5:5)
!!!!!! the transformation is
!!!!!! TMDtypeN=0:   kT^2\int_0^\infty db/(2pi) b J_0(b kT) F(b)
!!!!!! TMDtypeN=1:   kT^2 M**2/kT\int_0^\infty db/(2pi) b^2 J_1(b kT) F(b)
!!!!!! IMPORTANT THE RESULT IS MULTIPLIED by kT^2
function Fourier_Levin_array(F2)
real(dp),dimension(1:numKsubgrids, 0:kGridSize,-5:5)::Fourier_Levin_array
procedure(TMD_distribution)::F2
integer::i,j,ff
real(dp),dimension(1:lengthOfbValues,-5:5)::Function_Array

!!!! request array
if(TMDtypeN==0) then
    do i=1,lengthOfbValues
    Function_Array(i,:)=F2(bValues(i))*bValues(i)
    end do
else if(TMDtypeN==1) then
    do i=1,lengthOfbValues
    Function_Array(i,:)=F2(bValues(i))*bValues(i)*bValues(i)
    end do
!!!!!
end if

!!!!! TMDtypeN=0: factor kT^2/2pi is taken into account in the TranformationArray
!!!!! TMDtypeN=1: factor kT^2/2pi M^2/kT is taken into account in the TranformationArray
do i=1,numKsubgrids
do j=0,kGridSize
do ff=-5,5
    Fourier_Levin_array(i,j,ff)=dot_product(TranformationArray(i,j,:),Function_Array(:,ff))
end do
end do
end do

end function Fourier_Levin_array

!!!!! returns the array of weights, which should be multiplied by array of function at bNodes to get the Integral
!!!!! for TMDtypeN=0: \int_bMIN^bMAX J_0(b q) f(b) db
!!!!! for TMDtypeN=1: \int_bMIN^bMAX J_1(b q) f(b) db
!!!!! for the grids with bMAX*q<besselZERO*a, CC-quadrature is used, otherwise Levin
!!!!! the parameter a is any reasonable number
!!!!! --> if q is too small, the Levin matrix became unstable
!!!!! --> if q is too large, the CC quadrature does not work
!!!!! empirically I found that 1/2 workds good.
function FourierArrayAtQ(q)
real(dp),intent(in)::q
real(dp),dimension(1:lengthOfbValues)::FourierArrayAtQ

real(dp)::inv(0:bGridSize)
integer::i
real(dp)::valueTOswitch

!!!! this is value at which I change CC to Levin.
!!!! if (bq) is very small, the integral is not oscilating, and is better estimated by CC
!!!! It should be compared with BesselZero.
!!!! Empirically I found that it is better to put  as
valueTOswitch=besselZERO(TMDtypeN)*bGridSize**2/128.d0

!!!! if the integral is not oscilating, one better use CC-quadrature.
do i=1,numBsubgrids
    if(q*bRanges(i)<valueTOswitch) then
        inv=InverseCCOperator(i,q)
    else
        inv=InverseLevinOperator(i,q)
    end if
    FourierArrayAtQ(1+(i-1)*(bGridSize+1):i*(bGridSize+1))=inv(0:bGridSize)
end do

!!!! This is the correction for the extremely small-b
!!!! It is needed only for type=0, because in other cases, the integral drops fast at b->0.
if(TMDtypeN==0) then
    FourierArrayAtQ(bGridSize+1)=FourierArrayAtQ(bGridSize+1)+bessel_j1(bRanges(0)*q)/q
end if

contains

!!!! returns the derivative matrix for Chebyshev grid (see 2.19 in [2112.09703])
!!!! it does not contains the prefactor 1/b/([B-A]/2) recuired for proper derivative
pure function DerivativeMatrix(j,k)
integer,intent(in):: j,k
real(dp)::DerivativeMatrix
real(dp)::dummy
if(j==k) then
    if(j==0) then
        DerivativeMatrix=real(2*bGridSize**2+1,dp)/6
    else if (j==bGridSize) then
        DerivativeMatrix=-real(2*bGridSize**2+1,dp)/6
    else
        dummy=Cos(j*pi/bGridSize)
        DerivativeMatrix=-dummy/(2*(1-dummy**2))
    end if
else
    if(j==0 .or. j==bGridSize) then
        dummy=2._dp
    else
        dummy=1._dp
    end if
    if(k==0 .or. k==bGridSize) then
        dummy=dummy/2._dp
    end if
    if(mod(k+j,2)==1) then
        dummy=-dummy
    end if
    DerivativeMatrix=dummy/(Cos(j*pi/bGridSize)-Cos(k*pi/bGridSize))
end if
end function DerivativeMatrix

!!!! computes the element of the Levin matrix for the b-subgrid n, and q
!!!! matrix is (0:2*bGridSize+1)
!!!! the matrix has block form
!!!! |  -qI  D-I/b |
!!!! |   D    qI   |
!!!! where D is derivative matrix (divided by interval*b[i] because of derivative) and I is identiy matrix
!!!! each block is (bGridSize+1) x (bGridSize+1)
!!!! This matrix is for {h2,h1} and additionally trnasposed (because the matrix with {h1,h2} cannot be Gaussian Eliminated simply)
function ElementOfLevinMatrix(i,j,n,q)
real(dp)::ElementOfLevinMatrix
real(dp),intent(in)::q
integer,intent(in)::n
integer::i,j


! !!! |  -qI  D-I/b |
! !!! |   D    qI   |
! if(i<bGridSize+1 .and. j<bGridSize+1) then                      !!!!! upper-left block
!     if(i==j) then
!         ElementOfLevinMatrix=-q
!     else
!         ElementOfLevinMatrix=0._dp
!     end if
! else if (i<bGridSize+1 .and. j>bGridSize) then                !!!!! upper-right block
!     ElementOfLevinMatrix=DerivativeMatrix(i,j-bGridSize-1)/bIntervals(n)/BfromNode(n,i)
!     if(i==j-bGridSize-1) then
!         ElementOfLevinMatrix=ElementOfLevinMatrix-1/BfromNode(n,i)
!     end if
!
! else if (i>bGridSize .and. j<bGridSize+1) then                !!!!! lower-left block
!     ElementOfLevinMatrix=DerivativeMatrix(i-bGridSize-1,j)/bIntervals(n)/BfromNode(n,i-bGridSize-1)
! else
!     if(i==j) then
!         ElementOfLevinMatrix=q
!     else
!         ElementOfLevinMatrix=0._dp
!     end if
! end if

!!!! |   D         q I   |
!!!! |   -q I    D-I/b   |
if(i<bGridSize+1 .and. j<bGridSize+1) then                      !!!!! upper-left block
    ElementOfLevinMatrix=DerivativeMatrix(i,j)/bIntervals(n)/BfromNode(n,i)
else if (i<bGridSize+1 .and. j>bGridSize) then                !!!!! upper-right block
    if(i==j-bGridSize-1) then
        ElementOfLevinMatrix=q
    else
        ElementOfLevinMatrix=0._dp
    end if

else if (i>bGridSize .and. j<bGridSize+1) then                !!!!! lower-left block
    if(i-bGridSize-1==j) then
        ElementOfLevinMatrix=-q
    else
        ElementOfLevinMatrix=0._dp
    end if
else
    ElementOfLevinMatrix=DerivativeMatrix(i-bGridSize-1,j-bGridSize-1)/bIntervals(n)/BfromNode(n,i-bGridSize-1)
    if(i==j) then
        ElementOfLevinMatrix=ElementOfLevinMatrix-1/BfromNode(n,i-bGridSize-1)
    end if
end if

end function ElementOfLevinMatrix

!!!! finds the inverse of the Levin Matrix by Gaussian elimination
!!!! then builds the array A_{i}J_1(B q)-A_{i}J_1(A q)+A_{i}J_0(B q)-A_{i}J_0(A q)
!!!! which corresponds to the integral \int_A^B J_0(qb)
function InverseLevinOperator(n,q)
real(dp),intent(in)::q
integer,intent(in)::n
integer::NUM
real(dp),dimension(0:bGridSize)::InverseLevinOperator

real(dp),dimension(0:2*bGridSize+1,0:2*bGridSize+1)::MM!!!! matrix to inverse
real(dp),dimension(0:2*bGridSize+1,0:2*bGridSize+1)::invMM!!!! result of inverstion
integer::i,j

NUM=2*bGridSize+1


!!!! make a matrix to invert
do i=0,NUM
do j=0,NUM
    MM(i,j)=ElementOfLevinMatrix(i,j,n,q)
end do
end do

!!!!! inversion is done by Crout (realisation by https://github.com/Beliavsky/Matrix_Inversion)
invMM=inverse(MM)

!!!!! constract the transformation array
!!!!! NOTE: that t->1 corresponds to the lower limit, and t->-1 corresponds to the upper limit

!!!! for type=0, I take the upper part of the matrix
!!!! for type=1, I take the lower part of the matrix
if(TMDtypeN==0) then
    InverseLevinOperator=invMM(0,0:bGridSize)*bessel_j0(BfromNode(n,0)*q) &
                        -invMM(bGridSize,0:bGridSize)*bessel_j0(BfromNode(n,bGridSize)*q) &
                        +invMM(bGridSize+1,0:bGridSize)*bessel_j1(BfromNode(n,0)*q) &
                        -invMM(NUM,0:bGridSize)*bessel_j1(BfromNode(n,bGridSize)*q)
else if(TMDtypeN==1) then
    InverseLevinOperator=invMM(0,bGridSize+1:NUM)*bessel_j0(BfromNode(n,0)*q) &
                        -invMM(bGridSize,bGridSize+1:NUM)*bessel_j0(BfromNode(n,bGridSize)*q) &
                        +invMM(bGridSize+1,bGridSize+1:NUM)*bessel_j1(BfromNode(n,0)*q) &
                        -invMM(NUM,bGridSize+1:NUM)*bessel_j1(BfromNode(n,bGridSize)*q)
!!!
end if
end function InverseLevinOperator


!!!!! in some cases (very small q), the Levin-integration is not good.
!!!!! in these caes it is convenient to use Clenshaw Curtis since it uses the same nodes
!!!!! this function computes the cuadrature for the integral int J_0,1(bq)f[b] in the range n.
function InverseCCOperator(n,q)
real(dp),intent(in)::q
integer,intent(in)::n
real(dp),dimension(0:bGridSize)::InverseCCOperator

integer::i

if(TMDtypeN==0) then
    do i=0,bGridSize
        InverseCCOperator(i)=bIntervals(n)*BfromNode(n,i)*bessel_j0(q*BfromNode(n,i))*CCweight(i)
    end do
else if(TMDtypeN==1) then
    do i=0,bGridSize
        InverseCCOperator(i)=bIntervals(n)*BfromNode(n,i)*bessel_j1(q*BfromNode(n,i))*CCweight(i)
    end do
!!!!!
end if

end function InverseCCOperator

end function FourierArrayAtQ


!!!! precomputes the CCweights
!!!! the formula is
!!!! w_j=beta(j)4/N sum_{k=0,2,..}^n beta[k] Cos[j k pi/n]/(1-k^2)
!!!! where n is number of nodes, beta=1/2 for 0 or n.
subroutine CCweight_compute()
integer::j,k
do j=0,bGridSize
    CCweight(j)=0._dp
    do k=0,bGridSize,2
        if(k==0 .or. k==bGridSize) then
            CCweight(j)=CCweight(j)+Cos(j*k*pi/bGridSize)/2/(1._dp-k**2)
        else
            CCweight(j)=CCweight(j)+Cos(j*k*pi/bGridSize)/(1._dp-k**2)
        end if
    end do
    if(j==0 .or. j==bGridSize) CCweight(j)=CCweight(j)/2
    CCweight(j)=CCweight(j)*4._dp/bGridSize
end do
end subroutine CCweight_compute


!!!!! Some elementery check of Fourier
subroutine TestFourier()
real*8,dimension(-5:5)::FinQ
real*8,dimension(-5:5)::fromMath

write(*,*) " Checking some functions agains numerical integration in Mathematica "

write(*,*) " ----- at kt=0.1"
FinQ=Fourier_Levin(ToCheck,0.1d0)

if(TMDtypeN==0) then
    fromMath=(/0.3929447215912037d0, 0.02540379739979786d0, 0.908411878110134d0, &
    0.0005968092778014869d0, 0.6460648742389542d0, 0.06437676847530796d0, &
    0.7524510443091221d0, 0.49655034487385374d0, 0.34170225292532014d0, &
    -0.20752205491870596d0, -0.8852238038904702d0/)
else if(TMDtypeN==1) then
    fromMath=(/0.9823618039802421d0, 0.012174343801820444d0, 9.06266059229002d0, &
    -0.11935315213851815d0, 4.541014304087042d0, 0.15178780213084256d0, &
    5.168792156074321d0, 0.09872312114251779d0, 2.2245572492432713d0, &
    -2.6574335853522753d0, -15.058420056309387d0/)
end if

write(*,*) "Relative precision --->", FinQ/fromMath-1
write(*,*) "Absolute precision --->", FinQ-fromMath

write(*,*) " ----- at kt=1."

FinQ=Fourier_Levin(ToCheck,1.d0)

if(TMDtypeN==0) then
    fromMath=(/0.11399663659959798d0, 0.0203822972275855d0, -0.042202327319864175d0, &
    0.04313149087386549d0, 0.07117361597522101d0, 0.02972110683414643d0, &
    0.06294594297527262d0, 0.449934673142791d0, 0.042891077569203635d0, &
    0.02520683225073642d0, 0.052674083664594376d0/)
else if(TMDtypeN==1) then
    fromMath=(/0.28499159149899156d0, 0.008434054025207796d0, 0.10550581829966166d0, &
    -0.036137949441505786d0, 0.20024316087784436d0, 0.0341951163735701d0, &
    0.28005386757983064d0, 0.0898582807695214d0, 0.12755305478600654d0, &
    -0.00580021860141175d0, 0.166640906427954d0/)
end if

write(*,*) "Relative precision --->", FinQ/fromMath-1
write(*,*) "Absolute precision --->", FinQ-fromMath


write(*,*) " ----- at kt=10."

FinQ=Fourier_Levin(ToCheck,10.d0)

if(TMDtypeN==0) then
    fromMath=(/1.8747707813265846d-11, 0.0003633011162405713d0, &
    -0.000013741292555857953d0, 0.00016380750398139264d0, &
    1.2046139297866576d-8, 0.00005000066812517099d0, &
    -0.001859693917933602d0, 0.004722227621666028d0, &
    -0.0002109411666441108d0, -0.000042730758216799126d0, &
    -0.0003707316746467484d0/)
else if(TMDtypeN==1) then
    fromMath=(/8.099077398321386d-12, 0.000010257913870971234d0, &
    -6.756609413191178d-7, 5.004781471890181d-6, &
    1.9618241621743556d-9, 1.5567616589336672d-6, &
    -0.000036821940886544876d0, -0.0005685540405296604d0, &
    -4.5882687402352425d-6, 2.0882166810469566d-7, -0.000014544430739799435d0/)
end if

write(*,*) "Relative precision --->", FinQ/fromMath-1
write(*,*) "Absolute precision --->", FinQ-fromMath

write(*,*) " ----- at kt=100."

FinQ=Fourier_Levin(ToCheck,100.d0)

if(TMDtypeN==0) then
    fromMath=(/6.425713317692464d-12, 3.9751415835768074d-7, &
    -1.5915162873903991d-10, 1.592026776123332d-7, &
    2.7850115402326854d-12, 4.776807258675328d-8, &
    -0.000017098111555254256d0, 4.710252392769662d-9, &
    -9.231843125300607d-7, -7.159579504948952d-7, &
    1.5037898670724917d-6/)
else if(TMDtypeN==1) then
    fromMath=(/-1.881595114202385d-13, 1.19195798619532d-10, &
    -5.81190425910394d-13, 4.777034783574835d-11, -2.08572815045325d-13, 1.4363876501514724d-11, &
    -3.4829308907710516d-9, 2.6416349924401557d-11, -2.2449589346930395d-10, -1.6357764325066223d-10, &
    3.277343692223922d-10/)
end if

write(*,*) "Relative precision --->", FinQ/fromMath-1
write(*,*) "Absolute precision --->", FinQ-fromMath

contains
function ToCheck(b)
real*8,dimension(-5:5)::ToCheck
real*8,intent(in)::b

ToCheck=(/exp(-0.2*b**2),exp(-2.5*b),exp(-b)*b**2,exp(-b)*Cos(b),&
            (1.d0+0.2*b**2)/Cosh(b),1/(b**2+2.d0)*exp(-0.6*b), &
            exp(-0.8*b)*(1+sqrt(b)+log(b)),1.d0/(b**8+0.1d0),&
            Exp(-0.8*b)*b**(0.45), &
            Exp(-0.8*b)*(b**(0.45)-b**(0.93)),&
            Exp(-0.8*b)*(b**(0.45)*cos(3*b)-b**(0.93)*log(b)**2)            /)
end function ToCheck

end subroutine TestFourier


!end module Fourier_Levin
