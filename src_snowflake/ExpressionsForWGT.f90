!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                          !!
!! This file is a part of module EvolutionKernels, it contains                              !!
!! expressions for contraction matrix G2                                                    !!
!!                                                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!! the matrix WGTmatrix is stored globally
!!!!!! it is normal matrix (0:NUM_RHO)
!!!!!! the WGTmatrix should be contracted with S^+(normal)
!!!!!! there are two matricex WGTmatrixP and WGTmatrixN for case x>0, and x<0

!!!! This subroutine computes the kernel-matrix by the function Mfunction
!!!! and stores it into the matrix M.
subroutine PreComputeMatrixWGT(Mp,Mn)
real(dp), dimension(0:NUM_RHO,0:NUM_TOT)::Mp,Mn
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
    Mp(n,nn)=WGT_projM(nn,x1)
    Mn(n,nn)=WGT_projM(nn,-x1)

    if(ISNAN(Mp(n,nn))) then
        write(*,*) ErrorString("WGTp-element computed to NAN","")
        write(*,*) "--->",n,nn,x1
        error stop
    end if
    if(ISNAN(Mn(n,nn))) then
        write(*,*) ErrorString("WGTn-element computed to NAN","")
        write(*,*) "--->",n,nn,x1
        error stop
    end if

    end do
    write(*,*) "WGT part ",n, "/",NUM_RHO," is complete."
end do
!$OMP END PARALLEL DO

end subroutine PreComputeMatrixWGT

!!!! multiplies the matrix of G2``vector'' F (0... NUM_TOT)
!!!! and interpolate it to the point x
function WGTxF(x,F)
real(dp)::WGTxF
real(dp),dimension(0:NUM_TOT),intent(in)::F
real(dp),intent(in)::x
real(dp),dimension(0:NUM_TOT)::g_grid

if(x>0) then
    !!!! getting the grid of F
    g_grid=MATMUL(WGTmatrixP, F)
    WGTxF=GETinterpolation_RHO(x,g_grid)
else
    !!!! getting the grid of F
    g_grid=MATMUL(WGTmatrixN, F)
    WGTxF=GETinterpolation_RHO(-x,g_grid)
end if
end function WGTxF

!!!!---------------------------------------------------------------------!
!!!!-------------- Computation routines ---------------------------------!
!!!!---------------------------------------------------------------------!


!!!! Each kernel WGT_{n}(x) is a function of n (grid-interpolator-number), and x (at which x it is computed)
!!!! -------------
!!!! Definition of G2 is the following
!!!! -2x*\int I (1/z1**2 z3+1/z2**2 z1+1/z2**2 z3))
!!!! where each term has its own integration region. S is the S^-
!!!! in the last two terms there is also a subtraction at the point S(-x,0,x)
function WGT_projM(n,x)
integer,intent(in)::n
real(dp),intent(in)::x
real(dp)::WGT_projM

if(x>0) then
    WGT_projM=2._dp*x*(WGTp_0(n,x)+WGTp_1(n,x)+WGTp_2(n,x))
else
    WGT_projM=2._dp*x*(WGTn_0(n,x)+WGTn_1(n,x)+WGTn_2(n,x))
end if

end function WGT_projM

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!--------------------------------------------------------
!!!!!!!!!!!!!! INTEGRALS FOR X>0
!!!!!!!!!!!-----------------------------------------------------------------------------


!!!!! the part of function which is related to internal region (z1<-x,z3>x)
function WGTp_0(n,x)
real(dp)::WGTp_0
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
if(mG==3 .or. mG==4 .or. (mm==0 .and. mG==5)) then

    !!!!! Computing the integral over z1
    WGTp_0=Integrate_GK(f1,-1._dp,-x)
else
    WGTp_0=0._dp
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
f1=WGTp_0_int3(z1,n,z3MIN,z3MAX)
end function f1
end function WGTp_0

!!! Integration over z3 in (z3MIN,z3MAX) at fixed z1
function WGTp_0_int3(z1,n,z3MIN,z3MAX)
real(dp)::WGTp_0_int3
real(dp),intent(in)::z1,z3MIN,z3MAX
integer,intent(in)::n

WGTp_0_int3=Integrate_GK(f2,z3Min,z3Max)
contains
function f2(z3)
real(dp)::f2
real(dp),intent(in)::z3

f2=GETinterpolatorB(n,z1,-z1-z3)/(z1*z3**2)

end function f2
end function WGTp_0_int3

!!!!!----------------------------------------------------------------------
!!!!! the part of function which is related to internal region (z1<-x,z3<x)
function WGTp_1(n,x)
real(dp)::WGTp_1
integer,intent(in)::n
real(dp),intent(in)::x
integer::mG,nNUM,m_phi,mm

!!!! the simplest (althouhg not universal way) is the identify the integration limits here
!!!! in this case, I should not bother with xMIN, and boundary cases

!!!! getting subgrids
call get_2Dindex_from_1Dindex(n,nNUM,m_phi)
call subgridNUM_from_NUM_PHI(m_phi,mG,mm)

if(mG==2 .or. mG==3 .or. (mm==0 .and. mG==4)) then
    !!!!! Computing the integral over z1
    WGTp_1=Integrate_GK(f1,-1._dp,-x)
else
    WGTp_1=0._dp
    return
end if

contains

function f1(z1)
real(dp)::f1
real(dp),intent(in)::z1

real(dp):: z3MIN,z3MAX,sub

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
!!!! subtraction term
sub=GETinterpolatorB(n,-x,0._dp)
!!!!! sending to integral over z3  at given z1
f1=WGTp_1_int3(z1,n,z3Min,z3Max,sub)

end function f1
end function WGTp_1

!!! Integration over z3 at fixed z1
function WGTp_1_int3(z1,n,z3Min,z3Max,sub)
real(dp)::WGTp_1_int3
real(dp),intent(in)::z1,z3Min,z3Max,sub
integer,intent(in)::n
WGTp_1_int3=Integrate_GK(f2,z3Min,z3Max)
contains
function f2(z3)
real(dp)::f2
real(dp),intent(in)::z3
real(dp)::inter

inter=GETinterpolatorB(n,z1,-z1-z3)

!f2=(GETinterpolatorB(n,z1,-z1-z3)-sub)/((z1+z3)**2*z1)
f2=(inter-sub)/((z1+z3)**2*z1)+inter/((z1+z3)*z1**2)

end function f2
end function WGTp_1_int3

!!!!!!!!!!!-----------------------------------------------------------------------------
!!!!! the part of function which is related to internal region (z1>-x,z3>x)
function WGTp_2(n,x)
real(dp)::WGTp_2
integer,intent(in)::n
real(dp),intent(in)::x

integer::mG,nNUM,m_phi,mm
real(dp)::sub

!!!! the simplest (althouhg not universal way) is the identify the integration limits here
!!!! in this case, I should not bother with xMIN, and boundary cases

!!!! getting subgrids
call get_2Dindex_from_1Dindex(n,nNUM,m_phi)
call subgridNUM_from_NUM_PHI(m_phi,mG,mm)

!!!! subtraction term
sub=GETinterpolatorB(n,-x,0._dp)
!!!!! Computing the integral over z1
if(mG==3) then
    WGTp_2=Integrate_GK(f1,-x,0._dp)
else if(mm==0 .and. mG==4) then
    WGTp_2=Integrate_GK(f1,-x,1-x)
else if(mG==4) then
    WGTp_2=Integrate_GK(f1,0._dp,1-x)
else if(mm==0 .and. mG==5) then
    WGTp_2=Integrate_GK(f1,0._dp,1-x)
else
    WGTp_2=0._dp
    return
end if
contains

function f1(z1)
real(dp)::f1
real(dp),intent(in)::z1

!!!!! sending to integral over z3 at given z1 to a dedicated function
f1=WGTp_2_int3(z1,n,x,min(1-z1,1._dp),sub)
end function f1
end function WGTp_2

!!! Integration over z3 at fixed z1
function WGTp_2_int3(z1,n,z3Min,z3Max,sub)
real(dp)::WGTp_2_int3
real(dp),intent(in)::z1,z3Min,z3Max,sub
integer,intent(in)::n
WGTp_2_int3=Integrate_GK(f2,z3Min,z3Max)


contains
function f2(z3)
real(dp)::f2
real(dp),intent(in)::z3

f2=(GETinterpolatorB(n,z1,-z1-z3)-sub)/((z1+z3)**2*z3)

if(ISNAN(f2)) write(*,*) n,z1,-z1-z3, (GETinterpolatorB(n,z1,-z1-z3)-sub)/((z1+z3)**2*z3)

end function f2
end function WGTp_2_int3

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!--------------------------------------------------------
!!!!!!!!!!!!!! INTEGRALS FOR X>0
!!!!!!!!!!!-----------------------------------------------------------------------------


!!!!! the part of function which is related to internal region (z1>-x,z3<x)
function WGTn_0(n,x)
real(dp)::WGTn_0
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
if(mG==5 .or. mG==0 .or. (mm==0 .and. mG==1)) then

    !!!!! Computing the integral over z1
    WGTn_0=Integrate_GK(f1,-x,1._dp)
else
    WGTn_0=0._dp
    return
end if
contains

function f1(z1)
real(dp)::f1
real(dp),intent(in)::z1

real(dp)::z3MIN,z3MAX

if(mG==5) then
    z3MIN=-z1
    z3MAX=x
else if(mm==0 .and. mG==0) then
    z3MIN=-1._dp
    z3MAX=x
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
!!!!! the function is the same as in the positive case (just different limits)
f1=WGTp_0_int3(z1,n,z3MIN,z3MAX)
end function f1
end function WGTn_0


!!!!!----------------------------------------------------------------------
!!!!! the part of function which is related to internal region (z1<-x,z3<x)
function WGTn_1(n,x)
real(dp)::WGTn_1
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
    WGTn_1=Integrate_GK(f1,-x,1._dp)
else
    WGTn_1=0._dp
    return
end if

contains

function f1(z1)
real(dp)::f1
real(dp),intent(in)::z1

real(dp):: z3MIN,z3MAX,sub

if(mG==4) then
    z3MIN=0._dp
    z3MAX=1-z1
else if(mm==0 .and. mG==5) then
    z3MIN=x
    z3MAX=1-z1
else if(mG==5) then
    z3MIN=x
    z3MAX=0._dp
else if(mm==0 .and. mG==0) then
    z3MIN=x
    z3MAX=0._dp
else
    f1=0._dp
    return
end if
!!!! subtraction term
sub=GETinterpolatorB(n,-x,0._dp)
!!!!! sending to integral over z3  at given z1
!!!!! the function is the same as in the positive case (just different limits)
f1=WGTp_1_int3(z1,n,z3Min,z3Max,sub)

end function f1
end function WGTn_1


!!!!!!!!!!!-----------------------------------------------------------------------------
!!!!! the part of function which is related to internal region (z1>-x,z3>x)
function WGTn_2(n,x)
real(dp)::WGTn_2
integer,intent(in)::n
real(dp),intent(in)::x

integer::mG,nNUM,m_phi,mm
real(dp)::sub

!!!! the simplest (althouhg not universal way) is the identify the integration limits here
!!!! in this case, I should not bother with xMIN, and boundary cases

!!!! getting subgrids
call get_2Dindex_from_1Dindex(n,nNUM,m_phi)
call subgridNUM_from_NUM_PHI(m_phi,mG,mm)

!!!! subtraction term
sub=GETinterpolatorB(n,-x,0._dp)
!!!!! Computing the integral over z1
if(mG==0) then
    WGTn_2=Integrate_GK(f1,0._dp,-x)
else if(mm==0 .and. mG==1) then
    WGTn_2=Integrate_GK(f1,-1+x,-x)
else if(mG==1) then
    WGTn_2=Integrate_GK(f1,-1+x,0._dp)
else if(mm==0 .and. mG==2) then
    WGTn_2=Integrate_GK(f1,-1+x,0._dp)
else
    WGTn_2=0._dp
    return
end if
contains

function f1(z1)
real(dp)::f1
real(dp),intent(in)::z1

!!!!! sending to integral over z3 at given z1 to a dedicated function
!!!!! the function is the same as in the positive case (just different limits)
f1=WGTp_2_int3(z1,n,max(-1-z1,-1._dp),x,sub)
end function f1
end function WGTn_2

