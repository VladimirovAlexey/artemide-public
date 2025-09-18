!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                          !!
!! This file is a part of module EvolutionKernels, it contains                              !!
!! expressions for contraction matrix G2                                                    !!
!!                                                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!! the matrix G2matrix is stored globally
!!!!!! it is sparse matrix (0:NUM_RHO)

!!!! This subroutine computes the kernel-matrix by the function Mfunction
!!!! and stores it into the matrix M.
subroutine PreComputeMatrixG2(M)
real(dp), dimension(0:NUM_RHO,0:NUM_TOT)::M
real(dp), dimension(0:NUM_RHO,0:NUM_TOT)::Minit
real(dp), dimension(0:NUM_TOT)::Hinter,Dvector
integer::n,nn
real(dp)::x1,x2,x3
integer::OMP_GET_THREAD_NUM,thread_id

!!!!! original code used mutiplication by the derivative matrix.
!!!!! but it appered to be rather unpresice. Now, the computation is done directly with the derivative of the interpolator


!$OMP PARALLEL DO private(x1,x2,x3,nn)
do n=0,NUM_RHO
    do nn=0,NUM_TOT
    !!!! get value of x for given point
    !!!! this point is (x,0,-x)
    call get_X123_from_1Dindex(6*NUM_PHI*n,x1,x2,x3)

    !Minit(n,nn)=G2_projM(nn,x1)
    !!! minus, because the derivative with respect to dx3=-dx2
    M(n,nn)=-G2_projM(nn,x1)

    if(ISNAN(M(n,nn))) then
        write(*,*) ErrorString("G2-element computed to NAN","")
        write(*,*) "--->",n,nn,x1
        error stop
    end if

    end do
    write(*,*) "G2 part ",n, "/",NUM_RHO," is complete."
end do
!$OMP END PARALLEL DO

end subroutine PreComputeMatrixG2

!!!! multiplies the matrix of G2``vector'' F (0... NUM_TOT)
!!!! and interpolate it to the point x
function G2xF(x,F)
real(dp)::G2xF
real(dp),dimension(0:NUM_TOT),intent(in)::F
real(dp),intent(in)::x
real(dp),dimension(0:NUM_TOT)::g2_grid

!!!! getting the grid of F
g2_grid=MATMUL(G2matrix, F)
G2xF=GETinterpolation_RHO(x,g2_grid)
end function G2xF

!!!!---------------------------------------------------------------------!
!!!!-------------- Computation routines ---------------------------------!
!!!!---------------------------------------------------------------------!


!!!! Each kernel G2_{n}(x) is a function of n (grid-interpolator-number), and x (at which x it is computed)
!!!! -------------
!!!! Definition of G2 is the following
!!!! 1/2*\int I' (1/z1z3+1/z3(z1+z3)+1/z1(z1+z3))
!!!! where each term has its own integration region. I' is derivative of the interpolation function
function G2_projM(n,x)
integer,intent(in)::n
real(dp),intent(in)::x
real(dp)::G2_projM

G2_projM=(G2_0(n,x)+G2_1(n,x)+G2_2(n,x))/2
!G2_projM=G2_0(n,x)+G2_1(n,x)+G2_2(n,x)+G2_3(n,x)+G2_4(n,x)+G2_5(n,x)
!G2_projM=G2_1(n,x)

!
! G2_projM=G2_0(n,x)
! write(*,*) "0... ",G2_projM
! G2_projM=G2_projM+G2_1(n,x)
! write(*,*) "1... ",G2_projM
! G2_projM=G2_projM+G2_2(n,x)
! write(*,*) "2... ",G2_projM
! G2_projM=G2_projM+G2_3(n,x)
! write(*,*) "3... ",G2_projM
! G2_projM=G2_projM+G2_4(n,x)
! write(*,*) "4... ",G2_projM
! G2_projM=G2_projM+G2_5(n,x)
! write(*,*) "5... ",G2_projM

end function G2_projM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!--------------------------------------------------------
!!!!!!!!!!!!!! INTEGRALS FOR X>0
!!!!!!!!!!!-----------------------------------------------------------------------------


!!!!! the part of function which is related to internal region (z1>x,z3<-x)
function G2_0(n,x)
real(dp)::G2_0
integer,intent(in)::n
real(dp),intent(in)::x

integer::mG,nNUM,m_phi,mm

!!!! the simplest (althouhg not universal way) is the identify the integration limits here
!!!! in this case, I should not bother with xMIN, and boundary cases

!!!! getting subgrids
call get_2Dindex_from_1Dindex(n,nNUM,m_phi)
call subgridNUM_from_NUM_PHI(m_phi,mG,mm)
!!! for subgrid 1,2,3,4 the result is 0
!!! boundaries are never touched, because of xMIN
if(mG==0 .or. mG==5 .or. (mm==0 .and. mG==1)) then

    !!!!! Computing the integral over z1
    G2_0=Integrate_GK(f1,x,1._dp)
else
    G2_0=0._dp
    return
end if
contains

function f1(z1)
real(dp)::f1
real(dp),intent(in)::z1

real(dp)::z3MIN,z3MAX

if(mG==5) then
    z3MIN=-z1
    z3MAX=-x
else if(mm==0 .and. mG==0) then
    z3MIN=-1._dp
    z3MAX=-x
else if(mG==0) then
    z3MIN=-1._dp
    z3MAX=-z1
else if(mm==0 .and. mG==1) then
    z3MIN=-1._dp
    z3MAX=-z1
else
    f1=0._dp
    return
end if

!!!!! sending to integral over z3 at given z1 to a dedicated function
f1=G2_0_int3(z1,n,z3MIN,z3MAX)
end function f1
end function G2_0

!!! Integration over z3 in (z3MIN,z3MAX) at fixed z1
function G2_0_int3(z1,n,z3MIN,z3MAX)
real(dp)::G2_0_int3
real(dp),intent(in)::z1,z3MIN,z3MAX
integer,intent(in)::n

G2_0_int3=Integrate_GK(f2,z3Min,z3Max)
contains
function f2(z3)
real(dp)::f2
real(dp),intent(in)::z3

!f2=GETinterpolatorB(n,z1,-z1-z3)/(z1*z3)
f2=GETinterpolatorB_dX2(n,z1,-z1-z3)/(z1*z3)

end function f2
end function G2_0_int3

!!!!!!!!!!!-----------------------------------------------------------------------------
!!!!! the part of function which is related to internal region (z1<x,z3<-x)
function G2_1(n,x)
real(dp)::G2_1
integer,intent(in)::n
real(dp),intent(in)::x

integer::mG,nNUM,m_phi,mm

!!!! the simplest (althouhg not universal way) is the identify the integration limits here
!!!! in this case, I should not bother with xMIN, and boundary cases

!!!! getting subgrids
call get_2Dindex_from_1Dindex(n,nNUM,m_phi)
call subgridNUM_from_NUM_PHI(m_phi,mG,mm)

!!!!! Computing the integral over z1
if(mG==0) then
    G2_1=Integrate_GK(f1,0._dp,x)
else if(mm==0 .and. mG==1) then
    G2_1=Integrate_GK(f1,-1+x,x)
else if(mG==1) then
    G2_1=Integrate_GK(f1,-1+x,0._dp)
else if(mm==0 .and. mG==2) then
    G2_1=Integrate_GK(f1,-1+x,0._dp)
else
    G2_1=0._dp
    return
end if
contains

function f1(z1)
real(dp)::f1
real(dp),intent(in)::z1

!!!!! sending to integral over z3 at given z1 to a dedicated function
f1=G2_1_int3(z1,n,max(-1-z1,-1._dp),-x)
end function f1
end function G2_1

!!! Integration over z3 at fixed z1
function G2_1_int3(z1,n,z3Min,z3Max)
real(dp)::G2_1_int3
real(dp),intent(in)::z1,z3Min,z3Max
integer,intent(in)::n
G2_1_int3=Integrate_GK(f2,z3Min,z3Max)

!if(abs(z1)>zero) then
!    G2_1_int3=G2_1_int3+log(z3Max*(z3Min+z1)/(z3Min*(z3Max+z1)))/z1*GETinterpolatorB(n,z1,0._dp)
!else
!    G2_1_int3=G2_1_int3+(1._dp/z3Min-1._dp/z3Max)*GETinterpolatorB(n,z1,0._dp)
!end if

contains
function f2(z3)
real(dp)::f2
real(dp),intent(in)::z3

!f2=GETinterpolatorB(n,z1,-z1-z3)/((z1+z3)*z3)
f2=GETinterpolatorB_dX2(n,z1,-z1-z3)/((z1+z3)*z3)

if(ISNAN(f2)) write(*,*) n,z1,-z1-z3, GETinterpolatorB_dX2(n,z1,-z1-z3),1/((z1+z3)*z3)

!f2=(GETinterpolatorB(n,z1,-z1-z3)-GETinterpolatorB(n,z1,0._dp))/((z1+z3)*z3)

end function f2
end function G2_1_int3

!!!!!----------------------------------------------------------------------
!!!!! the part of function which is related to internal region (z1>x,z3>-x)
function G2_2(n,x)
real(dp)::G2_2
integer,intent(in)::n
real(dp),intent(in)::x
integer::mG,nNUM,m_phi,mm

!!!! the simplest (althouhg not universal way) is the identify the integration limits here
!!!! in this case, I should not bother with xMIN, and boundary cases

!!!! getting subgrids
call get_2Dindex_from_1Dindex(n,nNUM,m_phi)
call subgridNUM_from_NUM_PHI(m_phi,mG,mm)

if(mG==4 .or. mG==5 .or. (mm==0 .and. mG==0)) then
    !!!!! Computing the integral over z1
    G2_2=Integrate_GK(f1,x,1._dp)
else
    G2_2=0._dp
    return
end if

contains

function f1(z1)
real(dp)::f1
real(dp),intent(in)::z1

real(dp):: z3MIN,z3MAX

if(mG==4) then
    z3MAX=1-z1
    z3MIN=0._dp
else if(mm==0 .and. mG==5) then
    z3MAX=1-z1
    z3MIN=-x
else if(mG==5) then
    z3MAX=0._dp
    z3MIN=-x
else if(mm==0 .and. mG==0) then
    z3MAX=0._dp
    z3MIN=-x
else
    f1=0._dp
    return
end if
!!!!! sending to integral over z3  at given z1
f1=G2_2_int3(z1,n,z3Min,z3Max)

end function f1
end function G2_2

!!! Integration over z3 at fixed z1
function G2_2_int3(z1,n,z3Min,z3Max)
real(dp)::G2_2_int3
real(dp),intent(in)::z1,z3Min,z3Max
integer,intent(in)::n
G2_2_int3=Integrate_GK(f2,z3Min,z3Max)
contains
function f2(z3)
real(dp)::f2
real(dp),intent(in)::z3

!f2=GETinterpolatorB(n,z1,-z1-z3)/((z1+z3)*z1)
f2=GETinterpolatorB_dX2(n,z1,-z1-z3)/((z1+z3)*z1)

end function f2
end function G2_2_int3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!--------------------------------------------------------
!!!!!!!!!!!!!! INTEGRALS FOR X>0
!!!!!!!!!!!-----------------------------------------------------------------------------

!!!!!!!!!!!-----------------------------------------------------------------------------
!!!!! the part of function which is related to internal region (z3>x,z1<-x)
function G2_3(n,x)
real(dp)::G2_3
integer,intent(in)::n
real(dp),intent(in)::x

integer::mG,nNUM,m_phi,mm

!!!! the simplest (althouhg not universal way) is the identify the integration limits here
!!!! in this case, I should not bother with xMIN, and boundary cases

!!!! getting subgrids
call get_2Dindex_from_1Dindex(n,nNUM,m_phi)
call subgridNUM_from_NUM_PHI(m_phi,mG,mm)
!!! for subgrid 4,5,0,1 the result is 0
!!! boundaries are never touched, because of xMIN
if(mG==2 .or. mG==3 .or. (mm==0 .and. mG==4)) then
    !!!!! Computing the integral over z1
    G2_3=Integrate_GK(f1,-1._dp,-x)
else
    G2_3=0._dp
    return
end if

contains

function f1(z1)
real(dp)::f1
real(dp),intent(in)::z1

real(dp)::z3MIN,z3MAX

if(mG==2) then
    z3MIN=x
    z3MAX=-z1
else if(mm==0 .and. mG==3) then
    z3MIN=x
    z3MAX=1._dp
else if(mG==3) then
    z3MIN=-z1
    z3MAX=1._dp
else if(mm==0 .and. mG==4) then
    z3MIN=-z1
    z3MAX=1._dp
else
    f1=0._dp
    return
end if

!!!!! sending to integral over z3 at given z1 to a dedicated function
!!!!! the integral is the for G2_0 but with minus-sign
f1=-G2_0_int3(z1,n,z3MIN,z3MAX)

end function f1
end function G2_3

!!!!!!!!!!!-----------------------------------------------------------------------------
!!!!! the part of function which is related to internal region (z3<x,z1<-x)
function G2_4(n,x)
real(dp)::G2_4
integer,intent(in)::n
real(dp),intent(in)::x

integer::mG,nNUM,m_phi,mm

!!!! the simplest (althouhg not universal way) is the identify the integration limits here
!!!! in this case, I should not bother with xMIN, and boundary cases

!!!! getting subgrids
call get_2Dindex_from_1Dindex(n,nNUM,m_phi)
call subgridNUM_from_NUM_PHI(m_phi,mG,mm)

if(mG==1 .or. mG==2 .or. (mm==0 .and. mG==3)) then
    !!!!! Computing the integral over z1
    G2_4=Integrate_GK(f1,-1._dp,-x)
else
    G2_4=0._dp
    return
end if

contains

function f1(z1)
real(dp)::f1
real(dp),intent(in)::z1

real(dp):: z3MIN,z3MAX

if(mG==1) then
    z3MIN=-1-z1
    z3MAX=0._dp
else if(mm==0 .and. mG==2) then
    z3MIN=-1-z1
    z3MAX=x
else if(mG==2) then
    z3MIN=0._dp
    z3MAX=x
else if(mm==0 .and. mG==3) then
    z3MIN=0._dp
    z3MAX=x
else
    f1=0._dp
    return
end if
!!!!! sending to integral over z3  at given z1
!!!!! the integral is the same as G2_2 (with - sign)
f1=-G2_2_int3(z1,n,z3Min,z3Max)

end function f1
end function G2_4

!!!!!----------------------------------------------------------------------
!!!!! the part of function which is related to internal region (z3>x,z1>-x)
function G2_5(n,x)
real(dp)::G2_5
integer,intent(in)::n
real(dp),intent(in)::x
integer::mG,nNUM,m_phi,mm

!!!! the simplest (althouhg not universal way) is the identify the integration limits here
!!!! in this case, I should not bother with xMIN, and boundary cases

!!!! getting subgrids
call get_2Dindex_from_1Dindex(n,nNUM,m_phi)
call subgridNUM_from_NUM_PHI(m_phi,mG,mm)

!!!!! Computing the integral over z1
if(mG==3) then
    G2_5=Integrate_GK(f1,-x,0._dp)
else if(mm==0 .and. mG==4) then
    G2_5=Integrate_GK(f1,-x,1-x)
else if(mG==4) then
    G2_5=Integrate_GK(f1,0._dp,1-x)
else if(mm==0 .and. mG==5) then
    G2_5=Integrate_GK(f1,0._dp,1-x)
else
    G2_5=0._dp
    return
end if

contains

function f1(z1)
real(dp)::f1
real(dp),intent(in)::z1

!!!!! the integral is the same as G2_1 (with - sign)
f1=-G2_1_int3(z1,n,x,min(1-z1,1._dp))

end function f1
end function G2_5
