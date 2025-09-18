!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                          !!
!! This module operates with the grid and grid variables on the hexagon                     !!
!!                                                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module HexGrid
use IO_snowflake
implicit none

!!! 0=linear
!!! 1=logarithmic (!!!! use only this!!!! rest only for tests)
!!! 2=1/r-grid
#define GRIDTYPE 1

private

INCLUDE 'commonVariables.f90'

real(dp),parameter::pi=3.141592653589793238462643383279502884197169399375105820974944592_dp        !!pi

!!!! this parameter is used to scale grid [xmin,1]->[0,1] logarithmic
!!!! = log(xmin)
real(dp):: grid_param0

!!!! the list of chebyshev nodes, and factors for RHO-subgrid
real(dp),allocatable,dimension(:)::nodes_inRHO,nodeFactors_inRHO
!!!! the list of chebyshev nodes, and factors for PHI-subgrid
real(dp),allocatable,dimension(:)::nodes_inPHI,nodeFactors_inPHI

!!!! the matrix of differentiation of Chebyschev grid over phi, and rho
real(dp),allocatable,dimension(:,:)::D_der_RHO,D_der_PHI

public::Initialize_HexGrid
public::RP_fromX12,PHI_fromX
public::get_RhoP_from_1Dindex,get_X123_from_1Dindex,getSubgrid_PHI,RhoP_fromX12,X123_fromRP
public::get_2Dindex_from_1Dindex,subgridNUM_from_NUM_PHI
public::GETgrid,GETinterpolation,GETinterpolatorB,GETinterpolation_RHO,GETinterpolatorB_dX2,GETinterpolatorRHO_X
public::interpolatorB_in_grid

public::LimitsX2,LimitsX1,LimitsX3

public::Dgrid_dX2

contains

subroutine Initialize_HexGrid(path)
character(len=*)::path
integer::i,j

!----------------- reading ini-file --------------------------------------
OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")

call MoveTO(51,'0   :')
read(51,*) showINI
call MoveTO(51,'1   :')
read(51,*) showPROCESS

!!! Minimal value of X
call MoveTO(51,'A.1 :')
read(51,*) xMIN

!!! Number of nodes in each subgrid in R
call MoveTO(51,'A.2 :')
read(51,*) NUM_RHO

!!! Number of nodes in each subgrid in PHI
call MoveTO(51,'A.3 :')
read(51,*) NUM_PHI

!!!!------------------ numerical parameters
call MoveTO(51,'B.1 :')
read(51,*) zero

CLOSE (51, STATUS='KEEP')

grid_param0=log(xMIN)

!!!! the last integer for grid in RHO and PHI, 0... N (i.e. total size is +1)
NUM_TOT_RHO=NUM_RHO
NUM_TOT_PHI=6*NUM_PHI-1

!!!! the last integer for 1D grid in RHO and PHI, 0... N (i.e. total size is +1)
NUM_TOT=(NUM_TOT_RHO+1)*(NUM_TOT_PHI+1)-1


if(showINI) write(*,*) "Initialization of hex-grid ..."

if(showINI) write(*,'("  -- size of a subgrid      : (0,..., ", I3 ,")x(0,..., ", I3, ")")') NUM_RHO,NUM_PHI
if(showINI) write(*,'("  -- size of a 1D index     : 0, ...,", I5)') NUM_TOT

!!!! specifying nodes and factors in the subgrids
allocate(nodes_inRHO(0:NUM_RHO),nodeFactors_inRHO(0:NUM_RHO))
allocate(nodes_inPHI(0:NUM_PHI),nodeFactors_inPHI(0:NUM_PHI))

!!!! chebyshev node =cos(i pi/N)
!!!! chebyshev factor = (-1)^i beta_i
nodeFactors_inRHO=1._dp
do i=0,NUM_RHO
  nodes_inRHO(i)=cos(i*pi/NUM_RHO)

  if(i==0 .or. i==NUM_RHO) then
    nodeFactors_inRHO(i)=nodeFactors_inRHO(i)/2
  end if
  if(mod(i,2)==1) nodeFactors_inRHO(i)=-nodeFactors_inRHO(i)
end do

nodeFactors_inPHI=1._dp
do i=0,NUM_PHI
  nodes_inPHI(i)=cos(i*pi/NUM_PHI)

  if(i==0 .or. i==NUM_PHI) then
    nodeFactors_inPHI(i)=nodeFactors_inPHI(i)/2
  end if
  if(mod(i,2)==1) nodeFactors_inPHI(i)=-nodeFactors_inPHI(i)
end do

!!!! setting up derivative matrices
allocate(D_der_RHO(0:NUM_RHO,0:NUM_RHO),D_der_PHI(0:NUM_PHI,0:NUM_PHI))

do i=0,NUM_RHO
do j=0,NUM_RHO
    if(i==0 .and. j==0) then
        D_der_RHO(i,j)=(1._dp+2*NUM_RHO**2)/6._dp
    else if(i==NUM_RHO .and. j==NUM_RHO) then
        D_der_RHO(i,j)=-(1._dp+2*NUM_RHO**2)/6._dp
    else if(i==j) then
        D_der_RHO(i,j)=-nodes_inRHO(i)/2/(1-nodes_inRHO(i)**2)
    else
        D_der_RHO(i,j)=nodeFactors_inRHO(j)/nodeFactors_inRHO(i)/(nodes_inRHO(i)-nodes_inRHO(j))
    end if
end do
end do

do i=0,NUM_PHI
do j=0,NUM_PHI
    if(i==0 .and. j==0) then
        D_der_PHI(i,j)=(1._dp+2*NUM_PHI**2)/6._dp
    else if(i==NUM_PHI .and. j==NUM_PHI) then
        D_der_PHI(i,j)=-(1._dp+2*NUM_PHI**2)/6._dp
    else if(i==j) then
        D_der_PHI(i,j)=-nodes_inPHI(i)/2/(1-nodes_inPHI(i)**2)
    else
        D_der_PHI(i,j)=nodeFactors_inPHI(j)/nodeFactors_inPHI(i)/(nodes_inPHI(i)-nodes_inPHI(j))
    end if
end do
end do
!!! the factor -2 for phi-derivative is because of translation of phi->t
D_der_PHI=-2*D_der_PHI
!!! the factor -2 for phi-derivative is because of translation of rho->t
D_der_RHO=-2*D_der_RHO

if(showINI) write(*,*) "... initialization of hex-grid completed."
end subroutine Initialize_HexGrid

!!!!!!-----------------------------------------------------------------------------------------------------
!!!! --------------------------- OPERATIONS WITH INDICES AND COORDINATES-----------------------------------
!!!!!!-----------------------------------------------------------------------------------------------------

!!!! transformation from r to logarithm scale rho
!!!! r MUST be in [xmin,1]
!!!! r= xmin ->   rho=0
!!!! r= 1    ->   rho=1
pure function rho_fromR(r)
real(dp),intent(in)::r
real(dp)::rho_fromR
#if GRIDTYPE==1
rho_fromR=abs(1._dp-log(r)/grid_param0)
#elif GRIDTYPE==2
rho_fromR=(1-xMin/r)/(1-xMin)
#else
!!=0
rho_fromR=(r-xMin)/(1-xMin)
#endif

end function rho_fromR

!!!! transformation from the logarithm scale rho to R
!!!! rho MUST be in [0,1]
!!!! rho= 0  ->   R=xMIN
!!!! rho= 1  ->   R=1
pure function r_fromRHO(rho)
real(dp),intent(in)::rho
real(dp)::r_fromRHO

#if GRIDTYPE==1
r_fromRHO=exp(-grid_param0*(rho-1))
#elif GRIDTYPE==2
r_fromRHO=xMin/(1-rho+xMin*rho)
#else
!!=0
r_fromRHO=xMin+(1-xMin)*rho
#endif

end function r_fromRHO

!!!! Radius of the X-point on the hexagon
pure function R_fromX(x1,x2)
real(dp),intent(in)::x1,x2
real(dp)::R_fromX
R_fromX=max(abs(x1),abs(x2),abs(x1+x2))
end function R_fromX

!!!! Angle of the point on the hexagon
pure function PHI_fromX(x1,x2)
real(dp),intent(in)::x1,x2
real(dp)::PHI_fromX
real(dp)::r,x3
x3=-x1-x2
r=max(abs(x1),abs(x2),abs(x3))

if(x1>0 .and. x2>=0 .and. x3<0) then
    PHI_fromX=x2/r
else if(x1<=0 .and. x2>0 .and. x3<0) then
    PHI_fromX=(1-x1/r)
else if(x1<0 .and. x2>0 .and. x3>=0) then
    PHI_fromX=(3-x2/r)
else if(x1<0 .and. x2<=0 .and. x3>0) then
    PHI_fromX=(3-x2/r)
else if(x1>=0 .and. x2<0 .and. x3>0) then
    PHI_fromX=(4+x1/r)
else if(x1>0 .and. x2<0 .and. x3<=0) then
    PHI_fromX=(6+x2/r)
else
    PHI_fromX=0._dp
end if
end function PHI_fromX

!!!! returns (r,phi) from X12
pure subroutine RP_fromX12(x1,x2,r,phi)
real(dp),intent(in)::x1,x2
real(dp),intent(out)::r,phi
real(dp)::x3
x3=-x1-x2
r=max(abs(x1),abs(x2),abs(x3))

if(x1>0 .and. x2>=0 .and. x3<0) then
    phi=x2/r
else if(x1<=0 .and. x2>0 .and. x3<0) then
    phi=(1-x1/r)
else if(x1<0 .and. x2>0 .and. x3>=0) then
    phi=(3-x2/r)
else if(x1<0 .and. x2<=0 .and. x3>0) then
    phi=(3-x2/r)
else if(x1>=0 .and. x2<0 .and. x3>0) then
    phi=(4+x1/r)
else if(x1>0 .and. x2<0 .and. x3<=0) then
    phi=(6+x2/r)
else
    phi=0._dp
end if
end subroutine RP_fromX12

!!!! returns (rho,phi) from X12
pure subroutine RhoP_fromX12(x1,x2,rho,phi)
real(dp),intent(in)::x1,x2
real(dp),intent(out)::rho,phi
real(dp)::r

call RP_fromX12(x1,x2,r,phi)
rho=rho_fromR(r)

end subroutine RhoP_fromX12


!!! computes x123 from gives values of rho and phi
subroutine X123_fromRhoP(rho,phi,x1,x2,x3)
real(dp),intent(in)::rho,phi
real(dp),intent(out)::x1,x2,x3
real(dp)::r
r=r_fromRHO(rho)

call X123_fromRP(r,phi,x1,x2,x3)

end subroutine X123_fromRhoP

!!! computes x123 from gives values of R and phi
subroutine X123_fromRP(r,phi,x1,x2,x3)
real(dp),intent(in)::r,phi
real(dp),intent(out)::x1,x2,x3
real(dp)::pp

!!!! turn phi to [0,6]
if(phi<0) then
    pp=phi+6
else if (phi>6) then
    pp=phi-6
else
    pp=phi
end if

if(pp<1) then
    x1=r*(1-pp)
    x2=r*pp
    x3=-r
else if(pp<2) then
    x1=r*(1-pp)
    x2=r
    x3=r*(pp-2)
else if(pp<3) then
    x1=-r
    x2=r*(3-pp)
    x3=r*(pp-2)
else if(pp<4) then
    x1=r*(pp-4)
    x2=r*(3-pp)
    x3=r
else if(pp<5) then
    x1=r*(pp-4)
    x2=-r
    x3=r*(5-pp)
else
    x1=r
    x2=r*(pp-6)
    x3=r*(5-pp)
end if

end subroutine X123_fromRP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Ooperations with grid numbering !!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! returns the number of rho-subgrid, and the number of node within it from the single-rho number
!!!! Note, that the subgrid-number is computed as
!!!! mG=0 N=0,...,Nf-1
!!!! mG=2 N=Nf,...,2Nf-1
!!!! ....
!!!! mG=5 N=5Nf,...,6Nf-1 (cyclic)
pure subroutine subgridNUM_from_NUM_PHI(m,subgrid,node)
integer,intent(in)::m
integer,intent(out)::subgrid,node
subgrid=int(m/NUM_PHI)
node=mod(m,NUM_PHI)
end subroutine subgridNUM_from_NUM_PHI

!!!! returns the phi number from number of subgrid, and the number of node within it
pure subroutine NUM_from_subgridNUM_PHI(subgrid,node,m)
integer,intent(in)::subgrid,node
integer,intent(out)::m
m=subgrid*NUM_PHI+node
end subroutine NUM_from_subgridNUM_PHI

!!!! returns the coordinates 2D index from 1D
pure subroutine get_2Dindex_from_1Dindex(s,n,m)
integer,intent(in)::s
integer,intent(out)::n,m

n=int(s/(NUM_TOT_PHI+1))
m=mod(s,NUM_TOT_PHI+1)

end subroutine get_2Dindex_from_1Dindex

!!!! returns the coordinates 1D index from 2D
pure subroutine get_1Dindex_from_2Dindex(n,m,s)
integer,intent(in)::n,m
integer,intent(out)::s

s=(NUM_TOT_PHI+1)*n+m

end subroutine get_1Dindex_from_2Dindex

!!!! returns the coordinates (rho,phi) from 2D index
subroutine get_RhoP_from_1Dindex(n_in,rho,phi)
integer,intent(in)::n_in
real(dp),intent(out)::rho,phi

integer::n_out,m_out,mG,mm

call get_2Dindex_from_1Dindex(n_in,n_out,m_out)

!!!!! get numbers of subgrid and node in phi
call subgridNUM_from_NUM_PHI(m_out,mG,mm)

rho=(-nodes_inRHO(n_out)+1)/2
phi=(-nodes_inPHI(mm)+1)/2+mG

end subroutine get_RhoP_from_1Dindex

!!!! returns the number of subgrid mG from the 1D index
subroutine get_SubgridPHI_from_1Dindex(n_in,mG)
integer,intent(in)::n_in
real(dp),intent(out)::mG

mG=int(mod(n_in,NUM_TOT_PHI+1)/NUM_PHI)

end subroutine get_SubgridPHI_from_1Dindex

!!!! returns the coordinates (x1,x2,x3) from 1D index
subroutine get_X123_from_1Dindex(s,x1,x2,x3)
integer,intent(in)::s
real(dp),intent(out)::x1,x2,x3
real(dp)::rho, phi

call get_RhoP_from_1Dindex(s,rho,phi)

call X123_fromRhoP(rho,phi,x1,x2,x3)

end subroutine get_X123_from_1Dindex

!!!! find the number of subgrid for phi
!!!! it is assumed that phi in [-6,12] (i.e. I do not do the check for full circle only one extra circles)
pure function getSubgrid_PHI(phi_in)
real(dp),intent(in)::phi_in
integer::getSubgrid_PHI

if(phi_in>=6) then
    getSubgrid_PHI=int(phi_in-6)
else if(phi_in<0) then
    getSubgrid_PHI=int(phi_in+6)
else
    getSubgrid_PHI=int(phi_in)
end if

end function getSubgrid_PHI

!!!!-----------------------------------------------------------------------------------
!!!!----------------------------------------- INTERPOLATION ---------------------------
!!!!-----------------------------------------------------------------------------------
!!!! computes the grid from the given function
function GETgrid(func)
procedure(func_tw3)::func
real(dp),dimension(0:NUM_TOT)::GETgrid

integer::s
real(dp)::x1,x2,x3

do s=0,NUM_TOT
    call get_X123_from_1Dindex(s,x1,x2,x3)
    GETgrid(s)=func(x1,x2)
end do
end function GETgrid

!!!! computes the interpolation to the point x1,x2 from the 1D grid
!!!! It is not importnat is the point at the edge or not, because the interpolator at the edge is delta-function
function GETinterpolation(x1,x2,grid)
real(dp),intent(in)::x1,x2
!real(dp),dimension(0:NUM_TOT),intent(in)::grid
real(dp),dimension(0:),intent(in)::grid
real(dp)::GETinterpolation

real(dp)::rho,phi
integer::mG

!!!!!! 1) is to get variables in the terms of rho and phi
call RhoP_fromX12(x1,x2,rho,phi)

!!!!!! 2) to understand to which phi subgrid they belong, and if this point at the EDGE
mG=getSubgrid_PHI(phi)
!!!!! inside the grid just compute the interpolation
GETinterpolation=interpolation_in_grid(rho,phi,mG,grid)

! ! !!!!! debug
! ! if(ISNAN(GETinterpolation)) then
! ! write(*,*) "-------- NAN in snowflake GETinterpolation ------ "
! ! write(*,*) "rho,phi,mG",rho,phi,mG
! ! write(*,*) "grid ----->>>>> ",grid
! ! write(*,*) " "
! ! GETinterpolation=interpolation_in_grid_DEBUG(rho,phi,mG,grid)
! ! error stop
! ! end if

end function GETinterpolation


!!!!! compute the interpolation from the grid to a given point rho and phi in a given subgrid
!!!!! it is assumed that subgrid is already computed correctly
function interpolation_in_grid(rho,phi,mG,grid)
real(dp),intent(in)::rho,phi
integer,intent(in)::mG
!real(dp),dimension(0:NUM_TOT),intent(in)::grid
real(dp),dimension(0:),intent(in)::grid
real(dp)::interpolation_in_grid

real(dp)::t,u !!!! these are values of rho and phi in the subgrid space
real(dp)::deltaT(0:NUM_RHO),deltaU(0:NUM_PHI)

integer::i,j,n_here,m_here

!!!!!! 3) get coordinates in this subgrid
t=1-2*rho
u=1-2*(phi-mG)

!!!!!! 4) constract list of bary-centric weights for t, and u
deltaT=t-nodes_inRHO
deltaU=u-nodes_inPHI

!!!!!! 5) check that if the point is too close to the node, we replace the weight by 1, and rest ->0
!!!!!!    also I use this pass to make the proper baricentic formula =fact(i)/(t-t_i)
do i=0,NUM_RHO
    if(abs(deltaT(i))<zero) then
        deltaT=0._dp !!! all -> zero
        deltaT(i)=1._dp !!!! only i->1
        exit
    end if
    deltaT(i)=nodeFactors_inRHO(i)/deltaT(i)
end do

do j=0,NUM_PHI
    if(abs(deltaU(j))<zero) then
        deltaU=0._dp !!! all -> zero
        deltaU(j)=1._dp !!!! only i->1
        exit
    end if
    deltaU(j)=nodeFactors_inPHI(j)/deltaU(j)
end do

!!!!!! 6) Finally, I can interpolate
interpolation_in_grid=0._dp
do i=0,NUM_RHO
    !!!! the global index is computed as G=n*(NUM_TOT_PHI+1)+M
    !!!! where M=mG*NUM_PHI+m
    !!!! so I precompute some parts of index
    n_here=i*(NUM_TOT_PHI+1) !!!!!=n*(NUM_TOT_PHI+1)
    do j=0,NUM_PHI-1
        m_here=n_here+mG*NUM_PHI
        interpolation_in_grid=interpolation_in_grid+deltaT(i)*deltaU(j)*grid(m_here+j)
    end do
    !!!! since the last point is cycled I should take it separately
    if(mG==5) then
        interpolation_in_grid=interpolation_in_grid+deltaT(i)*deltaU(NUM_PHI)*grid(n_here)
    else
        interpolation_in_grid=interpolation_in_grid+deltaT(i)*deltaU(NUM_PHI)*grid(n_here+(mG+1)*NUM_PHI)
    end if
end do

!!!!!! 7) Devide by total weights
interpolation_in_grid=interpolation_in_grid/sum(deltaT)/sum(deltaU)

end function interpolation_in_grid

!!!!!!!!! ------ THIS VERSION OF THE FUNCTION IS USED FOR DEBUGGIN ------
! ! !!!!! compute the interpolation from the grid to a given point rho and phi in a given subgrid
! ! !!!!! it is assumed that subgrid is already computed correctly
! ! function interpolation_in_grid_DEBUG(rho,phi,mG,grid)
! ! real(dp),intent(in)::rho,phi
! ! integer,intent(in)::mG
! ! !real(dp),dimension(0:NUM_TOT),intent(in)::grid
! ! real(dp),dimension(0:),intent(in)::grid
! ! real(dp)::interpolation_in_grid_DEBUG
! !
! ! real(dp)::t,u !!!! these are values of rho and phi in the subgrid space
! ! real(dp)::deltaT(0:NUM_RHO),deltaU(0:NUM_PHI)
! !
! ! integer::i,j,n_here,m_here
! !
! ! !!!!!! 3) get coordinates in this subgrid
! ! t=1-2*rho
! ! u=1-2*(phi-mG)
! !
! ! write(*,*) "3)"
! ! write(*,*) "t=",t
! ! write(*,*) "u=",u
! ! write(*,*) " "
! !
! ! !!!!!! 4) constract list of bary-centric weights for t, and u
! ! deltaT=t-nodes_inRHO
! ! deltaU=u-nodes_inPHI
! !
! ! write(*,*) "4)"
! ! write(*,*) "deltaT=",deltaT
! ! write(*,*) " "
! ! write(*,*) "deltaU=",deltaU
! ! write(*,*) " "
! ! write(*,*) " "
! !
! ! !!!!!! 5) check that if the point is too close to the node, we replace the weight by 1, and rest ->0
! ! !!!!!!    also I use this pass to make the proper baricentic formula =fact(i)/(t-t_i)
! ! do i=0,NUM_RHO
! !     if(abs(deltaT(i))<zero) then
! !         deltaT=0._dp !!! all -> zero
! !         deltaT(i)=1._dp !!!! only i->1
! !         exit
! !     end if
! !     deltaT(i)=nodeFactors_inRHO(i)/deltaT(i)
! ! end do
! !
! ! do j=0,NUM_PHI
! !     if(abs(deltaU(j))<zero) then
! !         deltaU=0._dp !!! all -> zero
! !         deltaU(j)=1._dp !!!! only i->1
! !         exit
! !     end if
! !     deltaU(j)=nodeFactors_inPHI(j)/deltaU(j)
! ! end do
! !
! ! write(*,*) "5)"
! ! write(*,*) "deltaT=",deltaT
! ! write(*,*) " "
! ! write(*,*) "deltaU=",deltaU
! ! write(*,*) " "
! ! write(*,*) " "
! !
! ! !!!!!! 6) Finally, I can interpolate
! ! interpolation_in_grid_DEBUG=0._dp
! ! do i=0,NUM_RHO
! !     !!!! the global index is computed as G=n*(NUM_TOT_PHI+1)+M
! !     !!!! where M=mG*NUM_PHI+m
! !     !!!! so I precompute some parts of index
! !     n_here=i*(NUM_TOT_PHI+1) !!!!!=n*(NUM_TOT_PHI+1)
! !     do j=0,NUM_PHI-1
! !         m_here=n_here+mG*NUM_PHI
! !         interpolation_in_grid_DEBUG=interpolation_in_grid_DEBUG+deltaT(i)*deltaU(j)*grid(m_here+j)
! !         if(isnan(grid(m_here+j))) then
! !             write(*,*) "  :->",i,j,m_here,grid(m_here+j)
! !         end if
! !     end do
! !     !!!! since the last point is cycled I should take it separately
! !     if(mG==5) then
! !         interpolation_in_grid_DEBUG=interpolation_in_grid_DEBUG+deltaT(i)*deltaU(NUM_PHI)*grid(n_here)
! !     else
! !         interpolation_in_grid_DEBUG=interpolation_in_grid_DEBUG+deltaT(i)*deltaU(NUM_PHI)*grid(n_here+(mG+1)*NUM_PHI)
! !     end if
! ! end do
! !
! ! write(*,*) "6)"
! ! write(*,*) "interpolation_in_grid=",interpolation_in_grid_DEBUG
! ! write(*,*) " "
! !
! ! !!!!!! 7) Devide by total weights
! ! interpolation_in_grid_DEBUG=interpolation_in_grid_DEBUG/sum(deltaT)/sum(deltaU)
! !
! ! end function interpolation_in_grid_DEBUG


!!!! computes the interpolator function to the point x1,x2 from the grid which contains points n, relative to this point
!!!! IMPORTANT! if the point n is at the edge, the grid is the semi-sum of neighbour grids
function GETinterpolatorB(n_external,x1,x2,cc_in)
real(dp),intent(in)::x1,x2
integer,intent(in)::n_external
integer,intent(in),optional::cc_in !!!!=0 usual, =1 derivative in rho, =2 derivative in phi
real(dp)::GETinterpolatorB

logical:: atEDGE
real(dp)::rho,phi,inter1,inter2
integer::mG,mmG,nn,mm,m_phi,cc

if(present(cc_in)) then
    cc=cc_in
else
    cc=0
end if

!!!!!! 1) is to get variables in the terms of rho and phi
call RhoP_fromX12(x1,x2,rho,phi)

!!!!!! 2) to understand to which phi subgrid they belong, and if this point at the EDGE
mG=getSubgrid_PHI(phi)

! write(*,*) "INSIDE"
! write(*,*) "--->",x1,x2,rho,phi,mG

!!!!!!!!3) find out the 2D numbers of subgrids, of input n
!!!!!!!! if the subgrid of x is different from subgrid of n, result is zero.
!!!!!!!! however if it is at the edge, check more accurate
!!!!!!!! If the point N in the grid -> the point x also must belong to this grid (but with 1/2 for the interpolator)
!!!!!!!! If the point N is on the edge -> it must be the same edge as for X (and take the semi sum)
call get_2Dindex_from_1Dindex(n_external,nn,m_phi)
call subgridNUM_from_NUM_PHI(m_phi,mmG,mm)

!write(*,*) "--->",nn,mmG,mm,m_phi

if(mm==0) then!!! at the edge
    if((mmG==0 .and. (0==mG .or. 5==mG)) .or. (mmG==mG .or. mmG-1==mG)) then
        atEDGE=.true.
    else
        !!!! The secntor is not adjusted to the edge
        GETinterpolatorB=0._dp
        return
    end if
else
    if(mG/=mmG) then !!! inside the grid
        !!!! if sector is not the same
        GETinterpolatorB=0._dp
        return
    end if

    atEDGE=.false.
end if

if(atEDGE) then

    inter1=interpolatorB_in_grid(rho,phi,nn,mmG,0,cc)

    !write(*,*) "inter1=",inter1

    if(mmG==0) then
        if(abs(phi)<zero) then
            inter2=interpolatorB_in_grid(rho,phi+6,nn,5,NUM_PHI,cc)
        else
            inter2=interpolatorB_in_grid(rho,phi,nn,5,NUM_PHI,cc)
        end if
    else
        inter2=interpolatorB_in_grid(rho,phi,nn,mmG-1,NUM_PHI,cc)
    end if

    GETinterpolatorB=(inter1+inter2)

    if(abs(phi-mmG)<zero) then
        GETinterpolatorB=GETinterpolatorB/2 !!!! the point between the grids should be taken once
    else if(mmG==0 .and. abs(phi-6)<zero) then
        GETinterpolatorB=GETinterpolatorB/2 !!!! the point between the grids should be taken once (special case of cyclic grid)
    end if

else
    !!!! only half of term from one side
    GETinterpolatorB=interpolatorB_in_grid(rho,phi,nn,mmG,mm,cc)
end if


end function GETinterpolatorB

!!!!! compute the interpolator from at rho,phi of the grid that contains the point n,mG,mm
function interpolatorB_in_grid(rho,phi,n,mG,mm,cc)
real(dp),intent(in)::rho,phi
integer,intent(in)::n,mG,mm
integer,intent(in)::cc !!!!=0 usual, =1 derivative in rho, =2 derivative in phi
real(dp)::interpolatorB_in_grid,der_fac

real(dp)::t,u !!!! these are values of rho and phi in the subgrid space
real(dp)::deltaT(0:NUM_RHO),deltaU(0:NUM_PHI)

integer::i,tZ,uZ
logical::zeroT,zeroU

!!!!!! 3) get coordinates in this subgrid
t=1-2*rho
u=1-2*(phi-mG)

!!!! check that the interpolator is in the grid
if(abs(u)>1+zero .or. abs(t)>1+zero) then
    interpolatorB_in_grid=0._dp
    return
end if

!!!!!! 4) constract list of bary-centric weights for t, and u
deltaT=t-nodes_inRHO
deltaU=u-nodes_inPHI

!!!!!! check that non of these elements is zero
zeroT=.false.
zeroU=.false.
do i=0,NUM_RHO
    if(i==n) cycle
    if(abs(deltaT(i))<=zero) then
        tZ=i
        zeroT=.true.
        exit
    end if
end do
do i=0,NUM_PHI
    if(i==mm) cycle
    if(abs(deltaU(i))<=zero) then
        uZ=i
        zeroU=.true.
        exit
    end if
end do

!!!----------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!! Implementation via direct product
!!!! =(-1)^n beta_n/(2n)*prod 2(x-t_i)

!!!! if any of elements is zero, the interpolation is zero
if(cc==0 .and. (zeroT .or. zeroU)) then
    interpolatorB_in_grid=0._dp
    return
end if
if(cc==1 .and. zeroU) then
    interpolatorB_in_grid=0._dp
    return
end if
if(cc==2 .and. zeroT) then
    interpolatorB_in_grid=0._dp
    return
end if

interpolatorB_in_grid=1.d0
if(cc==1 .and. zeroT) then
    do i=0,NUM_RHO
        if(i==n .or. i==tZ) cycle
        interpolatorB_in_grid=interpolatorB_in_grid*deltaT(i)*2.d0
    end do
    !!! multiply by 2 because there factor 2 in the product
    interpolatorB_in_grid=2*interpolatorB_in_grid

    do i=0,NUM_PHI
        if(i==mm) cycle
        interpolatorB_in_grid=interpolatorB_in_grid*deltaU(i)*2.d0
    end do

    interpolatorB_in_grid=interpolatorB_in_grid*nodeFactors_inRHO(n)*nodeFactors_inPHI(mm)/(4*NUM_RHO*NUM_PHI)

else if(cc==2 .and. zeroU) then
    do i=0,NUM_RHO
        if(i==n) cycle
        interpolatorB_in_grid=interpolatorB_in_grid*deltaT(i)*2.d0
    end do

    do i=0,NUM_PHI
        if(i==mm .or. i==uZ) cycle
        interpolatorB_in_grid=interpolatorB_in_grid*deltaU(i)*2.d0
    end do
    !!! multiply by 2 because there factor 2 in the product
    interpolatorB_in_grid=2*interpolatorB_in_grid

    interpolatorB_in_grid=interpolatorB_in_grid*nodeFactors_inRHO(n)*nodeFactors_inPHI(mm)/(4*NUM_RHO*NUM_PHI)
else !!!! general case

    if(cc==0) then
        der_fac=1.d0
    else
        der_fac=0.d0
    end if

    do i=0,NUM_RHO
        if(i==n) cycle
        interpolatorB_in_grid=interpolatorB_in_grid*deltaT(i)*2.d0
        if(cc==1) der_fac=der_fac+1.d0/deltaT(i)
    end do

    do i=0,NUM_PHI
        if(i==mm) cycle
        interpolatorB_in_grid=interpolatorB_in_grid*deltaU(i)*2.d0
        if(cc==2) der_fac=der_fac+1.d0/deltaU(i)
    end do
    interpolatorB_in_grid=der_fac*interpolatorB_in_grid*nodeFactors_inRHO(n)*nodeFactors_inPHI(mm)/(4*NUM_RHO*NUM_PHI)
end if

!!!----------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!! Implementation via baricentric formula
! !!!!!! 5) check that if the point is too close to the node, we replace the weight by 1, and rest ->0
! !!!!!!    also I use this pass to make the proper baricentic formula =fact(i)/(t-t_i)
! do i=0,NUM_RHO
!     if(abs(deltaT(i))<zero) then
!         deltaT=0._dp !!! all -> zero
!         deltaT(i)=1._dp !!!! only i->1
!         exit
!     end if
!     deltaT(i)=nodeFactors_inRHO(i)/deltaT(i)
! end do
!
! do j=0,NUM_PHI
!     if(abs(deltaU(j))<zero) then
!         deltaU=0._dp !!! all -> zero
!         deltaU(j)=1._dp !!!! only i->1
!         exit
!     end if
!     deltaU(j)=nodeFactors_inPHI(j)/deltaU(j)
! end do

!!!!!! 6) Finally, I can interpolate
!interpolatorB_in_grid=deltaT(n)*deltaU(mm)/sum(deltaT)/sum(deltaU)

end function interpolatorB_in_grid

!!!! computes the derivative of interpolator over the x2
function GETinterpolatorB_dx2(n_external,x1,x2)
real(dp),intent(in)::x1,x2
integer,intent(in)::n_external
real(dp)::GETinterpolatorB_dx2

real(dp)::DRHO,DPHI

!!!! derivatives with respec to rho and phi
DRHO=GETinterpolatorB(n_external,x1,x2,1)
DPHI=GETinterpolatorB(n_external,x1,x2,2)

!!!! combine the complete derivative
GETinterpolatorB_dx2=-2*DRHO*dRho_dX2(x1,x2)-2*DPHI*dPHI_dX2(x1,x2)

end function GETinterpolatorB_dx2

!!!! computes the interpolation to the point x from the grid only vs. RHO (such as grids for g2)
function GETinterpolation_RHO(x,grid)
real(dp),intent(in)::x
real(dp),dimension(0:NUM_RHO),intent(in)::grid
real(dp)::GETinterpolation_RHO

real(dp)::rho,t
real(dp),dimension(0:NUM_RHO)::deltaT
integer::i

!!!!!! 1) is to get variables in the terms of rho (i.e. r=x)
rho=rho_fromR(x)

!!!!!! 3) get coordinates in this subgrid
t=1-2*rho

!!!!!! 4) constract list of bary-centric weights for t
deltaT=t-nodes_inRHO

!!!!!! 5) check that if the point is too close to the node, we replace the weight by 1, and rest ->0
!!!!!!    also I use this pass to make the proper baricentic formula =fact(i)/(t-t_i)
do i=0,NUM_RHO
    if(abs(deltaT(i))<zero) then
        deltaT=0._dp !!! all -> zero
        deltaT(i)=1._dp !!!! only i->1
        exit
    end if
    deltaT(i)=nodeFactors_inRHO(i)/deltaT(i)
end do

!!!!!! 6) Finally, I can interpolate
GETinterpolation_RHO=0._dp
do i=0,NUM_RHO
    GETinterpolation_RHO=GETinterpolation_RHO+deltaT(i)*grid(i)
end do

!!!!!! 7) Devide by total weights
GETinterpolation_RHO=GETinterpolation_RHO/sum(deltaT)

end function GETinterpolation_RHO

!!!!! compute the interpolator from at rho,phi of the grid that contains the point n,mG,mm
function GETinterpolatorRHO_X(n,x)
real(dp),intent(in)::x
integer,intent(in)::n
real(dp)::GETinterpolatorRHO_X

real(dp)::t,rho !!!! these are values of rho and phi in the subgrid space
real(dp)::deltaT(0:NUM_RHO)

integer::i

rho=rho_fromR(x)

!!!!!! 3) get coordinates in this subgrid
t=1-2*rho

!!!! check that the interpolator is in the grid
if(abs(t)>1+zero) then
    GETinterpolatorRHO_X=0._dp
    return
end if

!!!!!! 4) constract list of bary-centric weights for t, and u
deltaT=t-nodes_inRHO

!!!----------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!! Implementation via direct product
!!!! =(-1)^n beta_n/(2n)*prod 2(x-t_i)

GETinterpolatorRHO_X=1.d0

do i=0,NUM_RHO
    if(i==n) cycle
    GETinterpolatorRHO_X=GETinterpolatorRHO_X*deltaT(i)*2.d0
end do
GETinterpolatorRHO_X=GETinterpolatorRHO_X*nodeFactors_inRHO(n)/(2*NUM_RHO)

end function GETinterpolatorRHO_X

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Routines for derivatives of grids !!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! Derivative of rho with respect to x2
pure function dRho_dX2(x1,x2)
real(dp),intent(in)::x1,x2
real(dp)::dRho_dX2
real(dp)::x3
x3=-x1-x2
!
if(x1>0 .and. x2>=0 .and. x3<0) then
    dRho_dX2=1._dp
else if(x1<=0 .and. x2>0 .and. x3<0) then
    dRho_dX2=1._dp
else if(x1<0 .and. x2>0 .and. x3>=0) then
    dRho_dX2=0._dp
else if(x1<0 .and. x2<=0 .and. x3>0) then
    dRho_dX2=-1._dp
else if(x1>=0 .and. x2<0 .and. x3>0) then
    dRho_dX2=-1._dp
else if(x1>0 .and. x2<0 .and. x3<=0) then
    dRho_dX2=0._dp
else
    dRho_dX2=1.d32
end if

#if GRIDTYPE==1
dRho_dX2=-dRho_dX2/(grid_param0*R_fromX(x1,x2))
#elif GRIDTYPE==2
dRho_dX2=dRho_dX2*xMin/(1-xMin)/(R_fromX(x1,x2)**2)
#else
!!!==0
dRho_dX2=dRho_dX2/(1-xMin)
#endif

end function dRho_dX2

!!!! Derivative of phi with respect to x2
!!!! it is just x1/r^2
pure function dPHI_dX2(x1,x2)
real(dp),intent(in)::x1,x2
real(dp)::dPHI_dX2

dPHI_dX2=x1/(max(abs(x1),abs(x2),abs(x1+x2)))**2

end function dPHI_dX2

!!!! computes the derivative of of the grid with respect to rho, applying the derivative matrix
!!!! returns the grid
function Dgrid_dRHO(grid)
!real(dp), dimension(0:NUM_TOT),intent(in):: grid
real(dp), dimension(0:),intent(in):: grid
real(dp), dimension(0:NUM_TOT):: Dgrid_dRHO

integer::s,i,n

real(dp)::x1,x2,x3

do s=0,NUM_TOT
    Dgrid_dRHO(s)=0._dp

    call get_X123_from_1Dindex(s,x1,x2,x3)

    !!!! number of subgrid in n
    n=int(s/(NUM_TOT_PHI+1))
    do i=0,NUM_RHO

        Dgrid_dRHO(s)=Dgrid_dRHO(s)+D_der_RHO(n,i)*grid(s+(NUM_TOT_PHI+1)*(i-n))
    end do

    !if(abs(x1)<zero .and. x2>0) write(*,'("{",F12.8,",",F16.8,",",F16.8,"},")') x2, grid(s),Dgrid_dRHO(s)

end do

! integer::s,mm,mG,n,j
!
! do n=0,NUM_RHO
! do mG=0,5
! do mm=0,NUM_PHI-1!!!! important that it is N-1, because the last point belongs to the next grid
!     s=6*NUM_PHI*n+NUM_PHI*mG+mm
!     Dgrid_dRHO(s)=0._dp
!
!     do j=0,NUM_RHO !!!! but still this last point in included into thederivative of this grid
!         Dgrid_dRHO(s)=Dgrid_dRHO(s)+D_der_RHO(n,j)*grid(6*NUM_PHI*j+NUM_PHI*mG+mm)
!     end do
! end do
! end do
! end do

end function Dgrid_dRHO


!!!! computes the derivative of of the grid with respect to phi, applying the derivative matrix
!!!! returns the grid
function Dgrid_dPHI(grid)
!real(dp), dimension(0:NUM_TOT),intent(in):: grid
real(dp), dimension(0:),intent(in):: grid
real(dp), dimension(0:NUM_TOT):: Dgrid_dPHI

integer::s,n,mG,i,j

do n=0,NUM_RHO
do mG=0,5
do i=0,NUM_PHI-1!!!! important that it is N-1, because the last point belongs to the next grid
    s=6*NUM_PHI*n+NUM_PHI*mG+i
    Dgrid_dPHI(s)=0._dp

    do j=0,NUM_PHI !!!! but still this last point in included into thederivative of this grid
        if(mG==5 .and. j==NUM_PHI) then !!!! this is the cycling condition
            Dgrid_dPHI(s)=Dgrid_dPHI(s)+D_der_PHI(i,j)*grid(6*NUM_PHI*n)
        else
            Dgrid_dPHI(s)=Dgrid_dPHI(s)+D_der_PHI(i,j)*grid(s-i+j)
        end if
    end do
end do
end do
end do

end function Dgrid_dPHI

!!!! computes the derivative of the grid with respect to x2
!!!! returns the grid
function Dgrid_dX2(grid)
!real(dp), dimension(0:NUM_TOT),intent(in):: grid
real(dp), dimension(0:),intent(in):: grid
real(dp), dimension(0:NUM_TOT):: Dgrid_dX2

real(dp), dimension(0:NUM_TOT):: DgridRHO,DgridPHI
integer::s
real(dp)::x1,x2,x3

!!!! derivatives with respec to rho and phi
DgridRHO=Dgrid_dRHO(grid)
DgridPHI=Dgrid_dPHI(grid)

!!!! next run though all elements and combine them together
do s=0,NUM_TOT
    call get_X123_from_1Dindex(s,x1,x2,x3)

    Dgrid_dX2(s)=DgridRHO(s)*dRho_dX2(x1,x2)+DgridPHI(s)*dPHI_dX2(x1,x2)

end do

end function Dgrid_dX2

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!! Routines for determination of intersection !!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!! Subroutine checks if the line (x1-v,x3+v) intersects the grid which incorporates point n_in (in total-numbering)
!!!!!! if negativeV, then (x1+v,x3-v)
!!!!!! it returns intersect=T, if it does; set vMin and vMax
!!!!!! if there is no intersection then intersect=F, and vMin vMax are undefined.
subroutine LimitsX2(n_in,x1,x2,intersect,vMin,vMax,negativeV)
integer,intent(in)::n_in
real(dp),intent(in)::x1,x2
real(dp),intent(out)::vMin,vMax
logical,intent(out)::intersect
logical,optional,intent(in)::negativeV
logical::nV

integer::mG,nNUM,m_phi,mm
real(dp)::y1,z1,v1,v2


if(present(negativeV)) then
    nV=negativeV
else
    nV=.false.
end if

!!!! getting subgrids
call get_2Dindex_from_1Dindex(n_in,nNUM,m_phi)
call subgridNUM_from_NUM_PHI(m_phi,mG,mm)


if(x2>=xMin) then
    if(mG==0) then
        y1=1-x2
        z1=0._dp
    else if(mm==0 .and. mG==1) then !!! edge
        y1=1-x2
        z1=-x2
    else if(mG==1) then
        y1=0._dp
        z1=-x2
    else if(mm==0 .and. mG==2) then !!! edge
        y1=0._dp
        z1=-1._dp
    else if(mG==2) then
        y1=-x2
        z1=-1._dp
    else if(mm==0 .and. mG==3) then !!! further edge
        y1=-x2
        z1=-1._dp
    else
        intersect=.false.
        return
    end if
else if(x2>0) then !!!! case than radius is smaller than xMIN
    if(mG==0) then
        y1=1-x2
        z1=xMin-x2
    else if(mG==2) then
        y1=-xMin
        z1=-1._dp
    else if(mm==0 .and. mG==3) then !!! further edge
        y1=-xMin
        z1=-1._dp
    else
        !!! central sector does not have intersection
        intersect=.false.
        return
    end if
else if(x2<=-xMin) then
    if(mG==3) then
        y1=-1-x2
        z1=0._dp
    else if(mm==0 .and. mG==4) then !!! edge
        y1=-1-x2
        z1=-x2
    else if(mG==4) then
        y1=0._dp
        z1=-x2
    else if(mm==0 .and. mG==5) then !!! edge
        y1=0._dp
        z1=1._dp
    else if(mG==5) then
        y1=-x2
        z1=1._dp
    else if(mm==0 .and. mG==0) then !!! further edge
        y1=-x2
        z1=1._dp
    else
        intersect=.false.
        return
    end if
else if(x2<0) then !!!! case than radiouse is smaller than xMIN
    if(mG==3) then
        y1=-1-x2
        z1=-xMin-x2
    else if(mG==5) then
        y1=xMin
        z1=1._dp
    else if(mm==0 .and. mG==0) then !!! further edge
        y1=xMin
        z1=1._dp
    else
        !!! central sector does not have intersection
        intersect=.false.
        return
    end if
else !!!!! case of x2=0
    if(mG==0 .or. mG==5) then
        y1=xMin
        z1=1._dp
    else if(mG==2 .or. mG==3) then
        y1=-xMin
        z1=-1._dp
    else
        !!! central sector does not have intersection
        intersect=.false.
        return
    end if
end if

intersect=.true.

!!!!! These are intersection distances
v1=x1-y1
v2=x1-z1

if(nV) then
    v1=-v1
    v2=-v2
end if

!!!!! order v's
if(v1<v2) then
    vMin=v1
    vMax=v2
else
    vMin=v2
    vMax=v1
end if
!
! write(*,*) "Grid lims    >>>",r0,r1,phi0,phi1
! write(*,*) "intersection >>>",interLEFT,interRIGHT,interBottom,interTOP
! write(*,*) "point  >",x1,x2,-x1-x2
! write(*,*) "point 1>",y1,y2,y3
! write(*,*) "point 2>",z1,z2,z3
! write(*,*) "v1,v2-->",vMin,vMax

end subroutine LimitsX2


!!!!!! Subroutine checks if the line (x2-v,x3+v) intersects the grid which incorporates point n,k
!!!!!! if negativeV, then (x2+v,x3-v)
!!!!!! it returns intersect=T, if it does; set vMin and vMax
!!!!!! if there is no intersection then intersect=F, and vMin vMax are undefined.
subroutine LimitsX1(n_in,x1,x2,intersect,vMin,vMax,negativeV)
integer,intent(in)::n_in
real(dp),intent(in)::x1,x2
real(dp),intent(out)::vMin,vMax
logical,intent(out)::intersect
logical,optional,intent(in)::negativeV
logical::nV

integer::mG,nNUM,m_phi,mm
real(dp)::y2,z2,v1,v2


if(present(negativeV)) then
    nV=negativeV
else
    nV=.false.
end if

!!!! getting subgrids
call get_2Dindex_from_1Dindex(n_in,nNUM,m_phi)
call subgridNUM_from_NUM_PHI(m_phi,mG,mm)

if(x1>=xMin) then
    if(mG==4) then
        y2=-x1
        z2=-1._dp
    else if(mm==0 .and. mG==5) then !edge
        y2=0._dp
        z2=-1._dp
    else if(mG==5) then
        y2=0._dp
        z2=-x1
    else if(mm==0 .and. mG==0) then !edge
        y2=1-x1
        z2=-x1
    else if(mG==0) then
        y2=1-x1
        z2=0._dp
    else if(mm==0 .and. mG==1) then !!further edge
        y2=1-x1
        z2=0._dp
    else
        intersect=.false.
        return
    end if
else if(x1>0) then !!!! case than radiouse is smaller than xMIN
    if(mG==4) then
        y2=-xMin
        z2=-1._dp
    else if(mG==0) then
        y2=1-x1
        z2=xMin-x1
    else if(mm==0 .and. mG==1) then  !!further edge
        y2=1-x1
        z2=xMin-x1
    else
        !!! central sector does not have intersection
        intersect=.false.
        return
    end if
else if(x1<=-xMin) then
    if(mG==1) then
        y2=-x1
        z2=1._dp
    else if(mm==0 .and. mG==2) then !edge
        y2=0._dp
        z2=1._dp
    else if(mG==2) then
        y2=0._dp
        z2=-x1
    else if(mm==0 .and. mG==3) then !edge
        y2=-1-x1
        z2=-x1
    else if(mG==3) then
        y2=-1-x1
        z2=0._dp
    else if(mm==0 .and. mG==4) then
        y2=-1-x1
        z2=0._dp
    else
        intersect=.false.
        return
    end if
else if(x2<0) then !!!! case than radiouse is smaller than xMIN
    if(mG==1) then
        y2=xMin
        z2=1._dp
    else if(mG==3) then
        y2=-1-x1
        z2=-xMin-x1
    else if(mm==0 .and. mG==4) then  !!! further edge
        y2=-1-x1
        z2=-xMin-x1
    else
        !!! central sector does not have intersection
        intersect=.false.
        return
    end if
else !!!!! case of x2=0
    if(mG==0 .or. mG==1) then
        y2=xMin
        z2=1._dp
    else if(mG==3 .or. mG==4) then
        y2=-xMin
        z2=-1._dp
    else
        !!! central sector does not have intersection
        intersect=.false.
        return
    end if
end if


intersect=.true.

!!!!! These are intersection distances (x2-v,x3+v)
v1=x2-y2
v2=x2-z2

if(nV) then
    v1=-v1
    v2=-v2
end if

!!!!! order v's
if(v1<v2) then
    vMin=v1
    vMax=v2
else
    vMin=v2
    vMax=v1
end if
!
! write(*,*) "Grid lims    >>>",r0,r1,phi0,phi1
! write(*,*) "intersection >>>",interLEFT,interRIGHT,interBottom,interTOP
! write(*,*) "point  >",x1,x2,-x1-x2
! write(*,*) "point 1>",y1,y2,y3
! write(*,*) "point 2>",z1,z2,z3
! write(*,*) "v1,v2-->",vMin,vMax

end subroutine LimitsX1

!!!!!! Subroutine checks if the line (x1-v,x2+v) intersects the grid which incorporates point n,k
!!!!!! if negativeV, then (x1+v,x2-v)
!!!!!! it returns intersect=T, if it does; set vMin and vMax
!!!!!! if there is no intersection then intersect=F, and vMin vMax are undefined.
subroutine LimitsX3(n_in,x1,x2,intersect,vMin,vMax,negativeV)
integer,intent(in)::n_in
real(dp),intent(in)::x1,x2
real(dp),intent(out)::vMin,vMax
logical,intent(out)::intersect
logical,optional,intent(in)::negativeV
logical::nV

integer::mG,nNUM,m_phi,mm
real(dp)::y1,z1,v1,v2,x3

x3=-x1-x2
if(present(negativeV)) then
    nV=negativeV
else
    nV=.false.
end if

!!!! getting subgrids
call get_2Dindex_from_1Dindex(n_in,nNUM,m_phi)
call subgridNUM_from_NUM_PHI(m_phi,mG,mm)


if(x3>=xMin) then
    if(mG==2) then
        y1=-x3
        z1=-1._dp
    else if(mm==0 .and. mG==3) then !edge
        y1=0._dp
        z1=-1._dp
    else if(mG==3) then
        y1=0._dp
        z1=-x3
    else if(mm==0 .and. mG==4) then !edge
        y1=1-x3
        z1=-x3
    else if(mG==4) then
        y1=1-x3
        z1=0._dp
    else if(mm==0 .and. mG==5) then !! further edge
        y1=1-x3
        z1=0._dp
    else
        intersect=.false.
        return
    end if
else if(x3>0) then !!!! case than radiouse is smaller than xMIN
    if(mG==2) then
        y1=-xMin
        z1=-1._dp
    else if(mG==4) then
        y1=1-x3
        z1=xMin-x3
    else if(mm==0 .and. mG==5) then  !!! further edge
        y1=1-x3
        z1=xMin-x3
    else
        !!! central sector does not have intersection
        intersect=.false.
        return
    end if
else if(x3<=-xMin) then
    if(mG==5) then
        y1=-x3
        z1=1._dp
    else if(mm==0 .and. mG==0) then !edge
        y1=0._dp
        z1=1._dp
    else if(mG==0) then
        y1=0._dp
        z1=-x3
    else if(mm==0 .and. mG==1) then !edge
        y1=-x3
        z1=-1-x3
    else if(mG==1) then
        y1=-1-x3
        z1=0._dp
    else if(mm==0 .and. mG==2) then !!! further edge
        y1=-1-x3
        z1=0._dp
    else
        intersect=.false.
        return
    end if
else if(x3<0) then !!!! case than radiouse is smaller than xMIN
    if(mG==5) then
        y1=xMin
        z1=1._dp
    else if(mG==1) then
        y1=-1-x3
        z1=-xMin-x3
    else if(mm==0 .and. mG==2) then !!! further edge
        y1=-1-x3
        z1=-xMin-x3
    else
        !!! central sector does not have intersection
        intersect=.false.
        return
    end if
else !!!!! case of x3=0
    if(mG==4 .or. mG==5) then
        y1=xMin
        z1=1._dp
    else if(mG==2 .or. mG==1) then
        y1=-xMin
        z1=-1._dp
    else
        !!! central sector does not have intersection
        intersect=.false.
        return
    end if
end if

intersect=.true.

!!!!! These are intersection distances(x1-v,x2+v)
v1=x1-y1
v2=x1-z1

if(nV) then
    v1=-v1
    v2=-v2
end if

!!!!! order v's
if(v1<v2) then
    vMin=v1
    vMax=v2
else
    vMin=v2
    vMax=v1
end if
!
! write(*,*) "Grid lims    >>>",r0,r1,phi0,phi1
! write(*,*) "intersection >>>",interLEFT,interRIGHT,interBottom,interTOP
! write(*,*) "point  >",x1,x2,-x1-x2
! write(*,*) "point 1>",y1,y2,y3
! write(*,*) "point 2>",z1,z2,z3
! write(*,*) "v1,v2-->",vMin,vMax

end subroutine LimitsX3


end module HexGrid
