!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                          !!
!! This file is a part of module EvolutionKernels, it contains                              !!
!! expressions for elementary Kernels                                                       !!
!!                                                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! Each kernel H_{i,j}^{i',j'} is a function of integers i,j, ii=i', jj=j'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- Unity-operator ---------------------------------------------

!!!!! TODO::: It would be nice to update the kernels for a special case of (n,m); i.e. then x is x_M
!!!!!           this case, is most usefull, and for it the code can be somewhat simplified (especially search for limits)

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- \widehat{H}_{12}---------------------------------------------
function H12hat(nn,x1,x2)
real(dp)::H12hat
integer,intent(in)::nn
real(dp),intent(in)::x1,x2
logical::intersect
real(dp)::vMin,vMax

!write(*,*) "--->",nn,x1,x2,-x1-x2
call LimitsX3(nn,x1,x2,intersect,vMin,vMax)

!write(*,*) "--->>",nn,x1,x2,vMin,vMax,x1-vMin,x1

H12hat=0._dp

! write(*,*) "-------- ENTER ------"
! write(*,*) "nn=",nn
! write(*,*) "x1,x2,x3=",x1,x2,-x1-x2
! write(*,*) "v=",vMin,vMax

if(intersect) then
    if(x1>0 .and. vMin<0) H12hat=H12hat+Integrate_GK(f_1,vMin,min(0._dp,vMax))
    if(x1<0 .and. vMax>0) H12hat=H12hat-Integrate_GK(f_1,max(0._dp,vMin),vMax)

    if(x2>0 .and. vMax>0) H12hat=H12hat+Integrate_GK(f_2,max(0._dp,vMin),vMax)
    if(x2<0 .and. vMin<0) H12hat=H12hat-Integrate_GK(f_2,vMin,min(0._dp,vMax))
    !!! adding the rest of ()_+ terms
    if(x1>0 .and. vMin<0) H12hat=H12hat+log(1-x1/vMin)*GETinterpolatorB(nn,x1,x2)
    if(x1<0 .and. vMax>0) H12hat=H12hat+log(1-x1/vMax)*GETinterpolatorB(nn,x1,x2)
    if(x2>0 .and. vMax>0) H12hat=H12hat+log(1+x2/vMax)*GETinterpolatorB(nn,x1,x2)
    if(x2<0 .and. vMin<0) H12hat=H12hat+log(1+x2/vMin)*GETinterpolatorB(nn,x1,x2)
    if(x2==0) H12hat=H12hat+GETinterpolatorB(nn,x1,x2)
else
    H12hat=0._dp
end if


contains
    function f_1(v)
    real(dp)::f_1
    real(dp), intent(in)::v
        f_1=-x1/(x1-v)/v*(GETinterpolatorB(nn,x1,x2)-GETinterpolatorB(nn,x1-v,x2+v))
    end function f_1

    function f_2(v)
    real(dp)::f_2
    real(dp), intent(in)::v
            f_2=x2/v/(x2+v)*(GETinterpolatorB(nn,x1,x2)-x2/(v+x2)*GETinterpolatorB(nn,x1-v,x2+v))
    end function f_2

end function H12hat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- \widehat{H}_{23}---------------------------------------------
function H23hat(nn,x1,x2)
real(dp)::H23hat
integer,intent(in)::nn
real(dp),intent(in)::x1,x2
logical::intersect
real(dp)::vMin,vMax,x3

call LimitsX1(nn,x1,x2,intersect,vMin,vMax,negativeV=.true.)
x3=-x1-x2

H23hat=0._dp
if(intersect) then
    if(x3>0 .and. vMin<0) H23hat=H23hat+Integrate_GK(f_1,vMin,min(0._dp,vMax))
    if(x3<0 .and. vMax>0) H23hat=H23hat-Integrate_GK(f_1,max(0._dp,vMin),vMax)
    if(x2>0 .and. vMax>0) H23hat=H23hat+Integrate_GK(f_2,max(0._dp,vMin),vMax)
    if(x2<0 .and. vMin<0) H23hat=H23hat-Integrate_GK(f_2,vMin,min(0._dp,vMax))
    !!! adding the rest of ()_+ terms
    if(x3>0 .and. vMin<0) H23hat=H23hat+log(1-x3/vMin)*GETinterpolatorB(nn,x1,x2)
    if(x3<0 .and. vMax>0) H23hat=H23hat+log(1-x3/vMax)*GETinterpolatorB(nn,x1,x2)
    if(x2>0 .and. vMax>0) H23hat=H23hat+log(1+x2/vMax)*GETinterpolatorB(nn,x1,x2)
    if(x2<0 .and. vMin<0) H23hat=H23hat+log(1+x2/vMin)*GETinterpolatorB(nn,x1,x2)
    if(x2==0) H23hat=H23hat+GETinterpolatorB(nn,x1,x2)
else
    H23hat=0._dp
end if

contains
    function f_1(v)
    real(dp)::f_1
    real(dp), intent(in)::v
        f_1=-x3/(x3-v)/v*(GETinterpolatorB(nn,x1,x2)-GETinterpolatorB(nn,x1,x2+v))
    end function f_1

    function f_2(v)
    real(dp)::f_2
    real(dp), intent(in)::v
            f_2=x2/v/(x2+v)*(GETinterpolatorB(nn,x1,x2)-x2/(v+x2)*GETinterpolatorB(nn,x1,x2+v))
    end function f_2

end function H23hat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- \widehat{H}_{13}---------------------------------------------
function H13hat(nn,x1,x2)
real(dp)::H13hat
integer,intent(in)::nn
real(dp),intent(in)::x1,x2
logical::intersect
real(dp)::vMin,vMax,x3

call LimitsX2(nn,x1,x2,intersect,vMin,vMax)
x3=-x1-x2

H13hat=0._dp
if(intersect) then
    if(x1>0 .and. vMin<0) H13hat=H13hat+Integrate_GK(f_1,vMin,min(0._dp,vMax))
    if(x1<0 .and. vMax>0) H13hat=H13hat-Integrate_GK(f_1,max(0._dp,vMin),vMax)
    if(x3>0 .and. vMax>0) H13hat=H13hat+Integrate_GK(f_2,max(0._dp,vMin),vMax)
    if(x3<0 .and. vMin<0) H13hat=H13hat-Integrate_GK(f_2,vMin,min(0._dp,vMax))
    !!! adding the rest of ()_+ terms
    if(x1>0 .and. vMin<0) H13hat=H13hat+log(1-x1/vMin)*GETinterpolatorB(nn,x1,x2)
    if(x1<0 .and. vMax>0) H13hat=H13hat+log(1-x1/vMax)*GETinterpolatorB(nn,x1,x2)
    if(x3>0 .and. vMax>0) H13hat=H13hat+log(1+x3/vMax)*GETinterpolatorB(nn,x1,x2)
    if(x3<0 .and. vMin<0) H13hat=H13hat+log(1+x3/vMin)*GETinterpolatorB(nn,x1,x2)
else
    H13hat=0._dp
end if

contains
    function f_1(v)
    real(dp)::f_1
    real(dp), intent(in)::v
        f_1=-x1/(x1-v)/v*(GETinterpolatorB(nn,x1,x2)-GETinterpolatorB(nn,x1-v,x2))
    end function f_1

    function f_2(v)
    real(dp)::f_2
    real(dp), intent(in)::v
            f_2=x3/v/(x3+v)*(GETinterpolatorB(nn,x1,x2)-GETinterpolatorB(nn,x1-v,x2))
    end function f_2
end function H13hat

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- H^+_{12}---------------------------------------------
function H12plus(nn,x1,x2)
real(dp)::H12plus
integer,intent(in)::nn
real(dp),intent(in)::x1,x2
logical::intersect
real(dp)::vMin,vMax

call LimitsX3(nn,x1,x2,intersect,vMin,vMax)

! vMin=-2.d0
! vMax=2.d0

H12plus=0._dp
if(intersect) then
    if(abs(x1+x2)>zero) then
        if(x1>0 .and. vMin<0) H12plus=H12plus+Integrate_GK(f_1,vMin,min(0._dp,vMax))
        if(x1<0 .and. vMax>0) H12plus=H12plus-Integrate_GK(f_1,max(0._dp,vMin),vMax)
        if(x2>0 .and. vMax>0) H12plus=H12plus+Integrate_GK(f_2,max(0._dp,vMin),vMax)
        if(x2<0 .and. vMin<0) H12plus=H12plus-Integrate_GK(f_2,vMin,min(0._dp,vMax))
    else
        if(x1>0 .and. vMin<0) H12plus=H12plus+Integrate_GK(f_3,vMin,min(0._dp,vMax))
        if(x1<0 .and. vMax>0) H12plus=H12plus-Integrate_GK(f_3,max(0._dp,vMin),vMax)
    end if
else
    H12plus=0._dp
end if

contains
    function f_1(v)
    real(dp)::f_1
    real(dp), intent(in)::v
        f_1=x1*(2*x2+x1)/(x1-v)/(x1+x2)**2/2._dp*GETinterpolatorB(nn,x1-v,x2+v)
    end function f_1

    function f_2(v)
    real(dp)::f_2
    real(dp), intent(in)::v
            f_2=(x2/(x2+v)/(-x1-x2))**2*(2*x2+x1+v)/2._dp*GETinterpolatorB(nn,x1-v,x2+v)
    end function f_2

    function f_3(v)
    real(dp)::f_3
    real(dp), intent(in)::v
        f_3=v*(v-2*x1)/(x1-v)**3/2._dp*GETinterpolatorB(nn,x1-v,x2+v)
    end function f_3
end function H12plus

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- H^+_{23}---------------------------------------------
function H23plus(nn,x1,x2)
real(dp)::H23plus
integer,intent(in)::nn
real(dp),intent(in)::x1,x2
logical::intersect
real(dp)::vMin,vMax,x3

x3=-x1-x2
call LimitsX1(nn,x1,x2,intersect,vMin,vMax,negativeV=.true.)

H23plus=0._dp
if(intersect) then
    if(abs(x1)>zero) then
        if(x3>0 .and. vMin<0) H23plus=H23plus+Integrate_GK(f_1,vMin,min(0._dp,vMax))
        if(x3<0 .and. vMax>0) H23plus=H23plus-Integrate_GK(f_1,max(0._dp,vMin),vMax)
        if(x2>0 .and. vMax>0) H23plus=H23plus+Integrate_GK(f_2,max(0._dp,vMin),vMax)
        if(x2<0 .and. vMin<0) H23plus=H23plus-Integrate_GK(f_2,vMin,min(0._dp,vMax))
    else
        if(x3>0 .and. vMin<0) H23plus=H23plus+Integrate_GK(f_3,vMin,min(0._dp,vMax))
        if(x3<0 .and. vMax>0) H23plus=H23plus-Integrate_GK(f_3,max(0._dp,vMin),vMax)
    end if
else
    H23plus=0._dp
end if

contains
    function f_1(v)
    real(dp)::f_1
    real(dp), intent(in)::v
        f_1=x3*(x2-x1)/2/(x3-v)/(x1)**2*GETinterpolatorB(nn,x1,x2+v)
    end function f_1

    function f_2(v)
    real(dp)::f_2
    real(dp), intent(in)::v
            f_2=(x2/(x2+v)/x1)**2*(x2-x1+v)/2*GETinterpolatorB(nn,x1,x2+v)
    end function f_2

    function f_3(v)
    real(dp)::f_3
    real(dp), intent(in)::v
        f_3=v*(v-2*x3)/(x3-v)**3/2._dp*GETinterpolatorB(nn,x1,x2+v)
    end function f_3
end function H23plus

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- H^+_{13}---------------------------------------------
function H13plus(nn,x1,x2)
real(dp)::H13plus
integer,intent(in)::nn
real(dp),intent(in)::x1,x2
logical::intersect
real(dp)::vMin,vMax,x3

x3=-x1-x2
call LimitsX2(nn,x1,x2,intersect,vMin,vMax)

H13plus=0._dp
if(intersect) then
    if(abs(x2)>zero) then
        if(x1>0 .and. vMin<0) H13plus=H13plus+Integrate_GK(f_1,vMin,min(0._dp,vMax))
        if(x1<0 .and. vMax>0) H13plus=H13plus-Integrate_GK(f_1,max(0._dp,vMin),vMax)
        if(x3>0 .and. vMax>0) H13plus=H13plus+Integrate_GK(f_2,max(0._dp,vMin),vMax)
        if(x3<0 .and. vMin<0) H13plus=H13plus-Integrate_GK(f_2,vMin,min(0._dp,vMax))
    else
        if(x1>0 .and. vMin<0) H13plus=H13plus+Integrate_GK(f_3,vMin,min(0._dp,vMax))
        if(x1<0 .and. vMax>0) H13plus=H13plus-Integrate_GK(f_3,max(0._dp,vMin),vMax)
    end if

else
    H13plus=0._dp
end if


contains
    function f_1(v)
    real(dp)::f_1
    real(dp), intent(in)::v

    f_1=x1/x2/(v-x1)*GETinterpolatorB(nn,x1-v,x2)

    end function f_1

    function f_2(v)
    real(dp)::f_2
    real(dp), intent(in)::v

    f_2=-x3/x2/(v+x3)*GETinterpolatorB(nn,x1-v,x2)

    end function f_2

    function f_3(v)
    real(dp)::f_3
    real(dp), intent(in)::v
        f_3=-v/(v-x1)**2*GETinterpolatorB(nn,x1-v,x2)
    end function f_3

end function H13plus

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- H^-_{12}---------------------------------------------
function H12minus(nn,x1,x2)
real(dp)::H12minus
integer,intent(in)::nn
real(dp),intent(in)::x1,x2
logical::intersect
real(dp)::vMin,vMax

!!!! note exchanged arguments
call LimitsX3(nn,x2,x1,intersect,vMin,vMax,negativeV=.true.)

H12minus=0._dp
if(intersect) then
    if(abs(x1+x2)>zero) then
        if(x1>0 .and. vMin<0) H12minus=H12minus+Integrate_GK(f_1,vMin,min(0._dp,vMax))
        if(x1<0 .and. vMax>0) H12minus=H12minus-Integrate_GK(f_1,max(0._dp,vMin),vMax)
        if(x2>0 .and. vMax>0) H12minus=H12minus+Integrate_GK(f_2,max(0._dp,vMin),vMax)
        if(x2<0 .and. vMin<0) H12minus=H12minus-Integrate_GK(f_2,vMin,min(0._dp,vMax))
    else
        if(x1>0 .and. vMin<0) H12minus=H12minus+Integrate_GK(f_3,vMin,min(0._dp,vMax))
        if(x1<0 .and. vMax>0) H12minus=H12minus-Integrate_GK(f_3,max(0._dp,vMin),vMax)
    end if
else
    H12minus=0._dp
end if

contains
    function f_1(v)
    real(dp)::f_1
    real(dp), intent(in)::v
            f_1=x1*(2*x2*(x1-v)-x1*(x2+v))/2/(x1-v)**2/(-x1-x2)**2*GETinterpolatorB(nn,x2+v,x1-v)
    end function f_1

    function f_2(v)
    real(dp)::f_2
    real(dp), intent(in)::v
        f_2=x2**2/2/(-x1-x2)**2/(v+x2)*GETinterpolatorB(nn,x2+v,x1-v)
    end function f_2

    function f_3(v)
    real(dp)::f_3
    real(dp), intent(in)::v
        f_3=v**2/2/(x1-v)**3*GETinterpolatorB(nn,x2+v,x1-v)
    end function f_3
end function H12minus

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- H^-_{23}---------------------------------------------
function H23minus(nn,x1,x2)
real(dp)::H23minus
integer,intent(in)::nn
real(dp),intent(in)::x1,x2
logical::intersect
real(dp)::vMin,vMax,x3

x3=-x1-x2
!!!! note exchanged arguments
call LimitsX1(nn,x1,x3,intersect,vMin,vMax)

H23minus=0._dp
if(intersect) then
    if(abs(x1)>zero) then
        if(x3>0 .and. vMin<0) H23minus=H23minus+Integrate_GK(f_1,vMin,min(0._dp,vMax))
        if(x3<0 .and. vMax>0) H23minus=H23minus-Integrate_GK(f_1,max(0._dp,vMin),vMax)
        if(x2>0 .and. vMax>0) H23minus=H23minus+Integrate_GK(f_2,max(0._dp,vMin),vMax)
        if(x2<0 .and. vMin<0) H23minus=H23minus-Integrate_GK(f_2,vMin,min(0._dp,vMax))
    else
        if(x3>0 .and. vMin<0) H23minus=H23minus+Integrate_GK(f_3,vMin,min(0._dp,vMax))
        if(x3<0 .and. vMax>0) H23minus=H23minus-Integrate_GK(f_3,max(0._dp,vMin),vMax)
    end if
else
    H23minus=0._dp
end if

contains
    function f_1(v)
    real(dp)::f_1
    real(dp), intent(in)::v
            f_1=x3*(2*x2*(x3-v)-x3*(x2+v))/2/(x3-v)**2/(x1)**2*GETinterpolatorB(nn,x1,x3-v)
    end function f_1

    function f_2(v)
    real(dp)::f_2
    real(dp), intent(in)::v
        f_2=x2**2/2/(x1)**2/(v+x2)*GETinterpolatorB(nn,x1,x3-v)
    end function f_2

    function f_3(v)
    real(dp)::f_3
    real(dp), intent(in)::v
        f_3=v**2/2/(x3-v)**3*GETinterpolatorB(nn,x1,x3-v)
    end function f_3
end function H23minus

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- H^e_{23}P---------------------------------------------
function H23eP(nn,x1,x2)
real(dp)::H23eP
integer,intent(in)::nn
real(dp),intent(in)::x1,x2
logical::intersect
real(dp)::vMin,vMax,x3

x3=-x1-x2

!!!! note exchanged arguments
call LimitsX1(nn,x1,x3,intersect,vMin,vMax)

H23eP=0._dp
if(intersect) then
    if(Abs(x3)<zero) then
        H23eP=GETinterpolatorB(nn,x1,0._dp)
    else
        if(x3>0 .and. vMin<0) H23eP=H23eP+Integrate_GK(f_1,vMin,min(0._dp,vMax))
        if(x3<0 .and. vMax>0) H23eP=H23eP-Integrate_GK(f_1,max(0.,vMin),vMax)
    end if
else
    H23eP=0._dp
end if

contains
    function f_1(v)
    real(dp)::f_1
    real(dp), intent(in)::v
            f_1=(x3)/(x3-v)**2*GETinterpolatorB(nn,x1,x3-v)
    end function f_1
end function H23eP


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- H^d_{13}---------------------------------------------
function H13d(nn,x1,x2)
real(dp)::H13d
integer,intent(in)::nn
real(dp),intent(in)::x1,x2
logical::intersect
real(dp)::vMin,vMax,x3



x3=-x1-x2
call LimitsX2(nn,x1,x2,intersect,vMin,vMax)

if(intersect) then
    if(x3>0 .and. x1>0) then
        H13d=-(x1*x3)/x2**3*Integrate_GK(f_1,vMin,vMax)
    else if(x3<0 .and. x1<0) then
        H13d=(x1*x3)/x2**3*Integrate_GK(f_1,vMin,vMax)
    else
        H13d=0._dp
    end if

else
    H13d=0._dp
end if

contains
    function f_1(v)
    real(dp)::f_1
    real(dp), intent(in)::v
        f_1=GETinterpolatorB(nn,x1-v,x2)
    end function f_1
end function H13d

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! KERNELS FOR GLUONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- \widehat{H}_{12} for FFF---------------------------------------------
function H12hat_FFF(nn,x1,x2)
real(dp)::H12hat_FFF
integer,intent(in)::nn
real(dp),intent(in)::x1,x2
logical::intersect
real(dp)::vMin,vMax

call LimitsX3(nn,x1,x2,intersect,vMin,vMax)

H12hat_FFF=0._dp
if(intersect) then
    if(x1>0 .and. vMin<0) H12hat_FFF=H12hat_FFF+Integrate_GK(f_1,vMin,min(0._dp,vMax))
    if(x1<0 .and. vMax>0) H12hat_FFF=H12hat_FFF-Integrate_GK(f_1,max(0._dp,vMin),vMax)

    if(x2>0 .and. vMax>0) H12hat_FFF=H12hat_FFF+Integrate_GK(f_2,max(0._dp,vMin),vMax)
    if(x2<0 .and. vMin<0) H12hat_FFF=H12hat_FFF-Integrate_GK(f_2,vMin,min(0._dp,vMax))
    !!! adding the rest of ()_+ terms
    if(x1>0 .and. vMin<0) H12hat_FFF=H12hat_FFF+log(1-x1/vMin)*GETinterpolatorB(nn,x1,x2)
    if(x1<0 .and. vMax>0) H12hat_FFF=H12hat_FFF+log(1-x1/vMax)*GETinterpolatorB(nn,x1,x2)
    if(x2>0 .and. vMax>0) H12hat_FFF=H12hat_FFF+log(1+x2/vMax)*GETinterpolatorB(nn,x1,x2)
    if(x2<0 .and. vMin<0) H12hat_FFF=H12hat_FFF+log(1+x2/vMin)*GETinterpolatorB(nn,x1,x2)
    if(x1==0) H12hat_FFF=H12hat_FFF+GETinterpolatorB(nn,x1,x2)
    if(x2==0) H12hat_FFF=H12hat_FFF+GETinterpolatorB(nn,x1,x2)
else
    H12hat_FFF=0._dp
end if

contains
    function f_1(v)
    real(dp)::f_1
    real(dp), intent(in)::v
        f_1=-x1/(x1-v)/v*(GETinterpolatorB(nn,x1,x2)-x1/(x1-v)*GETinterpolatorB(nn,x1-v,x2+v))
    end function f_1

    function f_2(v)
    real(dp)::f_2
    real(dp), intent(in)::v
            f_2=x2/(x2+v)/v*(GETinterpolatorB(nn,x1,x2)-x2/(v+x2)*GETinterpolatorB(nn,x1-v,x2+v))
    end function f_2

end function H12hat_FFF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- \widehat{H}_{23} for FFF---------------------------------------------
function H23hat_FFF(nn,x1,x2)
real(dp)::H23hat_FFF
integer,intent(in)::nn
real(dp),intent(in)::x1,x2
logical::intersect
real(dp)::vMin,vMax,x3

call LimitsX1(nn,x1,x2,intersect,vMin,vMax)
x3=-x1-x2

H23hat_FFF=0._dp
if(intersect) then
    if(x3>0 .and. vMax>0) H23hat_FFF=H23hat_FFF+Integrate_GK(f_1,max(0._dp,vMin),vMax)
    if(x3<0 .and. vMin<0) H23hat_FFF=H23hat_FFF-Integrate_GK(f_1,vMin,min(0._dp,vMax))
    if(x2>0 .and. vMin<0) H23hat_FFF=H23hat_FFF+Integrate_GK(f_2,vMin,min(0._dp,vMax))
    if(x2<0 .and. vMax>0) H23hat_FFF=H23hat_FFF-Integrate_GK(f_2,max(0._dp,vMin),vMax)
    !!! adding the rest of ()_+ terms
    if(x3>0 .and. vMax>0) H23hat_FFF=H23hat_FFF+log(1+x3/vMax)*GETinterpolatorB(nn,x1,x2)
    if(x3<0 .and. vMin<0)H23hat_FFF=H23hat_FFF+log(1+x3/vMin)*GETinterpolatorB(nn,x1,x2)
    if(x2>0 .and. vMin<0) H23hat_FFF=H23hat_FFF+log(1-x2/vMin)*GETinterpolatorB(nn,x1,x2)
    if(x2<0 .and. vMax>0) H23hat_FFF=H23hat_FFF+log(1-x2/vMax)*GETinterpolatorB(nn,x1,x2)
    if(x2==0) H23hat_FFF=H23hat_FFF+GETinterpolatorB(nn,x1,x2)
    if(x3==0) H23hat_FFF=H23hat_FFF+GETinterpolatorB(nn,x1,x2)
else
    H23hat_FFF=0._dp
end if

contains
    function f_1(v)
    real(dp)::f_1
    real(dp), intent(in)::v
        f_1=x3/(x3+v)/v*(GETinterpolatorB(nn,x1,x2)-x3/(v+x3)*GETinterpolatorB(nn,x1,x2-v))
    end function f_1

    function f_2(v)
    real(dp)::f_2
    real(dp), intent(in)::v
        f_2=-x2/(x2-v)/v*(GETinterpolatorB(nn,x1,x2)-x2/(x2-v)*GETinterpolatorB(nn,x1,x2-v))
    end function f_2

end function H23hat_FFF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- \widehat{H}_{13} for FFF---------------------------------------------
function H13hat_FFF(nn,x1,x2)
real(dp)::H13hat_FFF
integer,intent(in)::nn
real(dp),intent(in)::x1,x2
logical::intersect
real(dp)::vMin,vMax,x3

call LimitsX2(nn,x1,x2,intersect,vMin,vMax)
x3=-x1-x2

H13hat_FFF=0._dp
if(intersect) then
    if(x1>0 .and. vMin<0) H13hat_FFF=H13hat_FFF+Integrate_GK(f_1,vMin,min(0._dp,vMax))
    if(x1<0 .and. vMax>0) H13hat_FFF=H13hat_FFF-Integrate_GK(f_1,max(0._dp,vMin),vMax)
    if(x3>0 .and. vMax>0) H13hat_FFF=H13hat_FFF+Integrate_GK(f_2,max(0._dp,vMin),vMax)
    if(x3<0 .and. vMin<0) H13hat_FFF=H13hat_FFF-Integrate_GK(f_2,vMin,min(0._dp,vMax))
    !!! adding the rest of ()_+ terms
    if(x1>0 .and. vMin<0) H13hat_FFF=H13hat_FFF+log(1-x1/vMin)*GETinterpolatorB(nn,x1,x2)
    if(x1<0 .and. vMax>0) H13hat_FFF=H13hat_FFF+log(1-x1/vMax)*GETinterpolatorB(nn,x1,x2)
    if(x3>0 .and. vMax>0) H13hat_FFF=H13hat_FFF+log(1+x3/vMax)*GETinterpolatorB(nn,x1,x2)
    if(x3<0 .and. vMin<0) H13hat_FFF=H13hat_FFF+log(1+x3/vMin)*GETinterpolatorB(nn,x1,x2)
    if(x1==0) H13hat_FFF=H13hat_FFF+GETinterpolatorB(nn,x1,x2)
    if(x3==0) H13hat_FFF=H13hat_FFF+GETinterpolatorB(nn,x1,x2)
else
    H13hat_FFF=0._dp
end if

contains
    function f_1(v)
    real(dp)::f_1
    real(dp), intent(in)::v
        f_1=-x1/(x1-v)/v*(GETinterpolatorB(nn,x1,x2)-x1/(x1-v)*GETinterpolatorB(nn,x1-v,x2))
    end function f_1

    function f_2(v)
    real(dp)::f_2
    real(dp), intent(in)::v
            f_2=x3/v/(x3+v)*(GETinterpolatorB(nn,x1,x2)-x3/(v+x3)*GETinterpolatorB(nn,x1-v,x2))
    end function f_2
end function H13hat_FFF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- H^+_{12} for FFF---------------------------------------------
function H12plus_FFF(nn,x1,x2)
real(dp)::H12plus_FFF
integer,intent(in)::nn
real(dp),intent(in)::x1,x2
logical::intersect
real(dp)::vMin,vMax

call LimitsX3(nn,x1,x2,intersect,vMin,vMax)

H12plus_FFF=0._dp
if(intersect) then
    if(abs(x1+x2)>zero) then
        if(x1>0 .and. vMin<0) H12plus_FFF=H12plus_FFF+Integrate_GK(f_1,vMin,min(0._dp,vMax))
        if(x1<0 .and. vMax>0) H12plus_FFF=H12plus_FFF-Integrate_GK(f_1,max(0._dp,vMin),vMax)
        if(x2>0 .and. vMax>0) H12plus_FFF=H12plus_FFF+Integrate_GK(f_2,max(0._dp,vMin),vMax)
        if(x2<0 .and. vMin<0) H12plus_FFF=H12plus_FFF-Integrate_GK(f_2,vMin,min(0._dp,vMax))
    else
        if(x1>0 .and. vMin<0) H12plus_FFF=H12plus_FFF+Integrate_GK(f_3,vMin,min(0._dp,vMax))
        if(x1<0 .and. vMax>0) H12plus_FFF=H12plus_FFF-Integrate_GK(f_3,max(0._dp,vMin),vMax)
    end if
else
    H12plus_FFF=0._dp
end if

contains
    function f_1(v)
    real(dp)::f_1
    real(dp), intent(in)::v
        f_1=(x1/(v-x1))**2*(-v*(x1+3*x2)+8*x1*x2+3*x1**2+3*x2**2)/(x1+x2)**3/6._dp*GETinterpolatorB(nn,x1-v,x2+v)
    end function f_1

    function f_2(v)
    real(dp)::f_2
    real(dp), intent(in)::v
            f_2=(x2/(v+x2))**2*(x2*(8*x1+v)+3*x1*(x1+v)+3*x2**2)/(x1+x2)**3/6._dp*GETinterpolatorB(nn,x1-v,x2+v)
    end function f_2

    function f_3(v)
    real(dp)::f_3
    real(dp), intent(in)::v
        f_3=-v*(-6*v*x1+6*x1*x1+v*v)/(v-x1)**4/6._dp*GETinterpolatorB(nn,x1-v,x2+v)
    end function f_3
end function H12plus_FFF



!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- H^+_{13} for FFF---------------------------------------------
function H13plus_FFF(nn,x1,x2)
real(dp)::H13plus_FFF
integer,intent(in)::nn
real(dp),intent(in)::x1,x2
logical::intersect
real(dp)::vMin,vMax,x3

x3=-x1-x2
call LimitsX2(nn,x1,x2,intersect,vMin,vMax)

H13plus_FFF=0._dp
if(intersect) then
    if(abs(x2)>zero) then
        if(x1>0 .and. vMin<0) H13plus_FFF=H13plus_FFF+Integrate_GK(f_1,vMin,min(0._dp,vMax))
        if(x1<0 .and. vMax>0) H13plus_FFF=H13plus_FFF-Integrate_GK(f_1,max(0._dp,vMin),vMax)
        if(x3>0 .and. vMax>0) H13plus_FFF=H13plus_FFF+Integrate_GK(f_2,max(0._dp,vMin),vMax)
        if(x3<0 .and. vMin<0) H13plus_FFF=H13plus_FFF-Integrate_GK(f_2,vMin,min(0._dp,vMax))
    else
        if(x1>0 .and. vMin<0) H13plus_FFF=H13plus_FFF+Integrate_GK(f_3,vMin,min(0._dp,vMax))
        if(x1<0 .and. vMax>0) H13plus_FFF=H13plus_FFF-Integrate_GK(f_3,max(0._dp,vMin),vMax)
    end if

else
    H13plus_FFF=0._dp
end if

contains
    function f_1(v)
    real(dp)::f_1
    real(dp), intent(in)::v

    f_1=(x1/(v-x1))**2*(-v*(x1+3*x3)+8*x1*x3+3*x1**2+3*x3**2)/(x1+x3)**3/6._dp*GETinterpolatorB(nn,x1-v,x2)

    end function f_1

    function f_2(v)
    real(dp)::f_2
    real(dp), intent(in)::v

    f_2=(x3/(v+x3))**2*(x3*(8*x1+v)+3*x1*(x1+v)+3*x3**2)/(x1+x3)**3/6._dp*GETinterpolatorB(nn,x1-v,x2)

    end function f_2

    function f_3(v)
    real(dp)::f_3
    real(dp), intent(in)::v
        f_3=-v*(-6*v*x1+6*x1*x1+v*v)/(v-x1)**4/6._dp*GETinterpolatorB(nn,x1-v,x2)
    end function f_3

end function H13plus_FFF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- H^-_{12} for FFF---------------------------------------------
function H12minus_FFF(nn,x1,x2)
real(dp)::H12minus_FFF
integer,intent(in)::nn
real(dp),intent(in)::x1,x2
logical::intersect
real(dp)::vMin,vMax

!!!! note exchanged arguments
call LimitsX3(nn,x2,x1,intersect,vMin,vMax,negativeV=.true.)

H12minus_FFF=0._dp
if(intersect) then
    if(abs(x1+x2)>zero) then
        if(x1>0 .and. vMin<0) H12minus_FFF=H12minus_FFF+Integrate_GK(f_1,vMin,min(0._dp,vMax))
        if(x1<0 .and. vMax>0) H12minus_FFF=H12minus_FFF-Integrate_GK(f_1,max(0._dp,vMin),vMax)
        if(x2>0 .and. vMax>0) H12minus_FFF=H12minus_FFF+Integrate_GK(f_2,max(0._dp,vMin),vMax)
        if(x2<0 .and. vMin<0) H12minus_FFF=H12minus_FFF-Integrate_GK(f_2,vMin,min(0._dp,vMax))
    else
        if(x1>0 .and. vMin<0) H12minus_FFF=H12minus_FFF+Integrate_GK(f_3,vMin,min(0._dp,vMax))
        if(x1<0 .and. vMax>0) H12minus_FFF=H12minus_FFF-Integrate_GK(f_3,max(0._dp,vMin),vMax)
    end if
else
    H12minus_FFF=0._dp
end if

contains
    function f_1(v)
    real(dp)::f_1
    real(dp), intent(in)::v
        f_1=-(x1/(x1-v))**2*(v*(x1+3*x2)-2*x1*x2)/(x1+x2)**3/6*GETinterpolatorB(nn,x2+v,x1-v)
    end function f_1

    function f_2(v)
    real(dp)::f_2
    real(dp), intent(in)::v
        f_2=(x2/(x2+v))**2*(v*(3*x1+x2)+2*x1*x2)/(x1+x2)**3/6*GETinterpolatorB(nn,x2+v,x1-v)
    end function f_2

    function f_3(v)
    real(dp)::f_3
    real(dp), intent(in)::v
        f_3=-v**3/6/(v-x1)**4*GETinterpolatorB(nn,x2+v,x1-v)
    end function f_3
end function H12minus_FFF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- H^-_{13} for FFF---------------------------------------------
function H13minus_FFF(nn,x1,x2)
real(dp)::H13minus_FFF
integer,intent(in)::nn
real(dp),intent(in)::x1,x2
logical::intersect
real(dp)::vMin,vMax,x3

x3=-x1-x2

!!!! note exchanged arguments
call LimitsX2(nn,x3,x2,intersect,vMin,vMax,negativeV=.true.)

H13minus_FFF=0._dp
if(intersect) then
    if(abs(x2)>zero) then
        if(x1>0 .and. vMin<0) H13minus_FFF=H13minus_FFF+Integrate_GK(f_1,vMin,min(0._dp,vMax))
        if(x1<0 .and. vMax>0) H13minus_FFF=H13minus_FFF-Integrate_GK(f_1,max(0._dp,vMin),vMax)
        if(x3>0 .and. vMax>0) H13minus_FFF=H13minus_FFF+Integrate_GK(f_2,max(0._dp,vMin),vMax)
        if(x3<0 .and. vMin<0) H13minus_FFF=H13minus_FFF-Integrate_GK(f_2,vMin,min(0._dp,vMax))
    else
        if(x1>0 .and. vMin<0) H13minus_FFF=H13minus_FFF+Integrate_GK(f_3,vMin,min(0._dp,vMax))
        if(x1<0 .and. vMax>0) H13minus_FFF=H13minus_FFF-Integrate_GK(f_3,max(0._dp,vMin),vMax)
    end if
else
    H13minus_FFF=0._dp
end if

contains
    function f_1(v)
    real(dp)::f_1
    real(dp), intent(in)::v
        f_1=-(x1/(x1-v))**2*(v*(x1+3*x3)-2*x1*x3)/(x1+x3)**3/6*GETinterpolatorB(nn,x3+v,x2)
    end function f_1

    function f_2(v)
    real(dp)::f_2
    real(dp), intent(in)::v
        f_2=(x3/(x3+v))**2*(v*(3*x1+x3)+2*x1*x3)/(x1+x3)**3/6*GETinterpolatorB(nn,x3+v,x2)
    end function f_2

    function f_3(v)
    real(dp)::f_3
    real(dp), intent(in)::v
        f_3=-v**3/6/(v-x1)**4*GETinterpolatorB(nn,x3+v,x2)
    end function f_3
end function H13minus_FFF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- tilde-H^+_{12} for FFF---------------------------------------------
function H12tilde_FFF(nn,x1,x2)
real(dp)::H12tilde_FFF
integer,intent(in)::nn
real(dp),intent(in)::x1,x2
logical::intersect
real(dp)::vMin,vMax

call LimitsX3(nn,x1,x2,intersect,vMin,vMax)

H12tilde_FFF=0._dp
if(intersect) then
    if(abs(x1+x2)>zero) then
        if(x1>0 .and. vMin<0) H12tilde_FFF=H12tilde_FFF+Integrate_GK(f_1,vMin,min(0._dp,vMax))
        if(x1<0 .and. vMax>0) H12tilde_FFF=H12tilde_FFF-Integrate_GK(f_1,max(0._dp,vMin),vMax)
        if(x2>0 .and. vMax>0) H12tilde_FFF=H12tilde_FFF+Integrate_GK(f_2,max(0._dp,vMin),vMax)
        if(x2<0 .and. vMin<0) H12tilde_FFF=H12tilde_FFF-Integrate_GK(f_2,vMin,min(0._dp,vMax))
    else
        if(x1>0 .and. vMin<0) H12tilde_FFF=H12tilde_FFF+Integrate_GK(f_3,vMin,min(0._dp,vMax))
        if(x1<0 .and. vMax>0) H12tilde_FFF=H12tilde_FFF-Integrate_GK(f_3,max(0._dp,vMin),vMax)
    end if
else
    H12tilde_FFF=0._dp
end if

contains
    function f_1(v)
    real(dp)::f_1
    real(dp), intent(in)::v
        f_1=-(x1/(x1-v))**2*(v*(x1+3*x2)-2*x1*x2)/(x1+x2)**3/6*GETinterpolatorB(nn,x1-v,x2+v)
    end function f_1

    function f_2(v)
    real(dp)::f_2
    real(dp), intent(in)::v
            f_2=(x2/(x2+v))**2*(v*(3*x1+x2)+2*x1*x2)/(x1+x2)**3/6*GETinterpolatorB(nn,x1-v,x2+v)
    end function f_2

    function f_3(v)
    real(dp)::f_3
    real(dp), intent(in)::v
        f_3=-v**3/6/(v-x1)**4*GETinterpolatorB(nn,x1-v,x2+v)
    end function f_3
end function H12tilde_FFF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- tilde-H^+_{13} for FFF---------------------------------------------
function H13tilde_FFF(nn,x1,x2)
real(dp)::H13tilde_FFF
integer,intent(in)::nn
real(dp),intent(in)::x1,x2
logical::intersect
real(dp)::vMin,vMax,x3

x3=-x1-x2
call LimitsX2(nn,x1,x2,intersect,vMin,vMax)

H13tilde_FFF=0._dp
if(intersect) then
    if(abs(x2)>zero) then
        if(x1>0 .and. vMin<0) H13tilde_FFF=H13tilde_FFF+Integrate_GK(f_1,vMin,min(0._dp,vMax))
        if(x1<0 .and. vMax>0) H13tilde_FFF=H13tilde_FFF-Integrate_GK(f_1,max(0._dp,vMin),vMax)
        if(x3>0 .and. vMax>0) H13tilde_FFF=H13tilde_FFF+Integrate_GK(f_2,max(0._dp,vMin),vMax)
        if(x3<0 .and. vMin<0) H13tilde_FFF=H13tilde_FFF-Integrate_GK(f_2,vMin,min(0._dp,vMax))
    else
        if(x1>0 .and. vMin<0) H13tilde_FFF=H13tilde_FFF+Integrate_GK(f_3,vMin,min(0._dp,vMax))
        if(x1<0 .and. vMax>0) H13tilde_FFF=H13tilde_FFF-Integrate_GK(f_3,max(0._dp,vMin),vMax)
    end if

else
    H13tilde_FFF=0._dp
end if

contains
    function f_1(v)
    real(dp)::f_1
    real(dp), intent(in)::v

    f_1=-(x1/(x1-v))**2*(v*(x1+3*x3)-2*x1*x3)/(x1+x3)**3/6*GETinterpolatorB(nn,x1-v,x2)

    end function f_1

    function f_2(v)
    real(dp)::f_2
    real(dp), intent(in)::v

    f_2=(x3/(x3+v))**2*(v*(3*x1+x3)+2*x1*x3)/(x1+x3)**3/6*GETinterpolatorB(nn,x1-v,x2)

    end function f_2

    function f_3(v)
    real(dp)::f_3
    real(dp), intent(in)::v
        f_3=-v**3/6/(v-x1)**4*GETinterpolatorB(nn,x1-v,x2)
    end function f_3

end function H13tilde_FFF

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! KERNELS FOR QUAKR-GLUON !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- V^+_{13} ---------------------------------------------
function V13plus(nn,x1,x2)
real(dp)::V13plus
integer,intent(in)::nn
real(dp),intent(in)::x1,x2
logical::intersect
real(dp)::vMin,vMax,x3

x3=-x1-x2
call LimitsX2(nn,x1,x2,intersect,vMin,vMax)

V13plus=0._dp
if(intersect) then
    if(abs(x2)>zero) then
        if(x1>0 .and. vMin<0) V13plus=V13plus+Integrate_GK(f_1,vMin,min(0._dp,vMax))
        if(x1<0 .and. vMax>0) V13plus=V13plus-Integrate_GK(f_1,max(0._dp,vMin),vMax)
        if(x1==0) V13plus=V13plus-GETinterpolatorB(nn,x1,x2)/x3
        if(x3>0 .and. vMax>0) V13plus=V13plus+Integrate_GK(f_2,max(0._dp,vMin),vMax)
        if(x3<0 .and. vMin<0) V13plus=V13plus-Integrate_GK(f_2,vMin,min(0._dp,vMax))
        if(x3==0) V13plus=V13plus+GETinterpolatorB(nn,x1,x2)/x1
    else
        if(x1>0 .and. vMin<0) V13plus=V13plus+Integrate_GK(f_3,vMin,min(0._dp,vMax))
        if(x1<0 .and. vMax>0) V13plus=V13plus-Integrate_GK(f_3,max(0._dp,vMin),vMax)
    end if

else
    V13plus=0._dp
end if

contains
    function f_1(v)
    real(dp)::f_1
    real(dp), intent(in)::v

    f_1=-x1*x3/(v-x1)**2*(3*x1+x3-2*v)/(x1+x3)**3*GETinterpolatorB(nn,x1-v,x2)

    end function f_1

    function f_2(v)
    real(dp)::f_2
    real(dp), intent(in)::v

    f_2=x1*x3/(v+x3)**2*(x1+3*x3+2*v)/(x1+x3)**3*GETinterpolatorB(nn,x1-v,x2)

    end function f_2

    function f_3(v)
    real(dp)::f_3
    real(dp), intent(in)::v
        f_3=x1**2/(v-x1)**4*GETinterpolatorB(nn,x1-v,x2)
    end function f_3

end function V13plus

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- V^-_{13} ---------------------------------------------
function V13minus(nn,x1,x2)
real(dp)::V13minus
integer,intent(in)::nn
real(dp),intent(in)::x1,x2
logical::intersect
real(dp)::vMin,vMax,x3

x3=-x1-x2

!!!! note exchanged arguments
call LimitsX2(nn,x3,x2,intersect,vMin,vMax,negativeV=.true.)

V13minus=0._dp
if(intersect) then
    if(abs(x2)>zero) then
        if(x1>0 .and. vMin<0) V13minus=V13minus+Integrate_GK(f_1,vMin,min(0._dp,vMax))
        if(x1<0 .and. vMax>0) V13minus=V13minus-Integrate_GK(f_1,max(0._dp,vMin),vMax)
        if(x3>0 .and. vMax>0) V13minus=V13minus+Integrate_GK(f_2,max(0._dp,vMin),vMax)
        if(x3<0 .and. vMin<0) V13minus=V13minus-Integrate_GK(f_2,vMin,min(0._dp,vMax))
    else
        if(x1>0 .and. vMin<0) V13minus=V13minus+Integrate_GK(f_3,vMin,min(0._dp,vMax))
        if(x1<0 .and. vMax>0) V13minus=V13minus-Integrate_GK(f_3,max(0._dp,vMin),vMax)
    end if
else
    V13minus=0._dp
end if

contains
    function f_1(v)
    real(dp)::f_1
    real(dp), intent(in)::v
        f_1=x1/(x1-v)**2*(-x1*x3+x1*x1+2*v*x3)/(x1+x3)**3*GETinterpolatorB(nn,x3+v,x2)
    end function f_1

    function f_2(v)
    real(dp)::f_2
    real(dp), intent(in)::v
        f_2=x3/(x3+v)**2*(x1*x3-x3*x3+2*v*x1)/(x1+x3)**3*GETinterpolatorB(nn,x3+v,x2)
    end function f_2

    function f_3(v)
    real(dp)::f_3
    real(dp), intent(in)::v
        f_3=v**2/(v-x1)**4*GETinterpolatorB(nn,x3+v,x2)
    end function f_3
end function V13minus


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! KERNELS FOR GLUON-QUARK !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- W^+ ---------------------------------------------
function Wplus(nn,x1,x2)
real(dp)::Wplus
integer,intent(in)::nn
real(dp),intent(in)::x1,x2
logical::intersect
real(dp)::vMin,vMax,x3

x3=-x1-x2
call LimitsX2(nn,x1,x2,intersect,vMin,vMax)

Wplus=0._dp
if(intersect) then
    if(x1>0 .and. vMin<0) Wplus=Wplus+Integrate_GK(f_1,vMin,min(0._dp,vMax))
    if(x1<0 .and. vMax>0) Wplus=Wplus-Integrate_GK(f_1,max(0._dp,vMin),vMax)
    if(x3>0 .and. vMax>0) Wplus=Wplus+Integrate_GK(f_2,max(0._dp,vMin),vMax)
    if(x3<0 .and. vMin<0) Wplus=Wplus-Integrate_GK(f_2,vMin,min(0._dp,vMax))
else
    Wplus=0._dp
end if

contains
    function f_1(v)
    real(dp)::f_1
    real(dp), intent(in)::v

    f_1=-0.5_dp*GETinterpolatorB(nn,x1-v,x2)

    end function f_1

    function f_2(v)
    real(dp)::f_2
    real(dp), intent(in)::v

    f_2=0.5_dp*GETinterpolatorB(nn,x1-v,x2)

    end function f_2

end function Wplus

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- W^+P_23 ---------------------------------------------
function WplusP(nn,x1,x2)
real(dp)::WplusP
integer,intent(in)::nn
real(dp),intent(in)::x1,x2

WplusP=Wplus(nn,x1,-x1-x2)

end function WplusP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- W^- ---------------------------------------------
function Wminus(nn,x1,x2)
real(dp)::Wminus
integer,intent(in)::nn
real(dp),intent(in)::x1,x2
logical::intersect
real(dp)::vMin,vMax,x3

x3=-x1-x2

!!!! note exchanged arguments
call LimitsX2(nn,x3,x2,intersect,vMin,vMax,negativeV=.true.)

Wminus=0._dp
if(intersect) then
    if(abs(x2)>zero) then
        if(x1>0 .and. vMin<0) Wminus=Wminus+Integrate_GK(f_1,vMin,min(0._dp,vMax))
        if(x1<0 .and. vMax>0) Wminus=Wminus-Integrate_GK(f_1,max(0._dp,vMin),vMax)
        if(x3>0 .and. vMax>0) Wminus=Wminus+Integrate_GK(f_2,max(0._dp,vMin),vMax)
        if(x3<0 .and. vMin<0) Wminus=Wminus-Integrate_GK(f_2,vMin,min(0._dp,vMax))
    else
        if(x1>0 .and. vMin<0) Wminus=Wminus+Integrate_GK(f_3,vMin,min(0._dp,vMax))
        if(x1<0 .and. vMax>0) Wminus=Wminus-Integrate_GK(f_3,max(0._dp,vMin),vMax)
    end if
else
    Wminus=0._dp
end if

contains
    function f_1(v)
    real(dp)::f_1
    real(dp), intent(in)::v
        f_1=(x1**2/(v-x1)/(x1+x3)+0.5_dp)*GETinterpolatorB(nn,x3+v,x2)
    end function f_1

    function f_2(v)
    real(dp)::f_2
    real(dp), intent(in)::v
        f_2=(x3**2/(v+x3)/(x1+x3)-0.5_dp)*GETinterpolatorB(nn,x3+v,x2)
    end function f_2

    function f_3(v)
    real(dp)::f_3
    real(dp), intent(in)::v
        f_3=v**2/(v-x1)**2*GETinterpolatorB(nn,x3+v,x2)
    end function f_3
end function Wminus

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- W^-P_23 ---------------------------------------------
function WminusP(nn,x1,x2)
real(dp)::WminusP
integer,intent(in)::nn
real(dp),intent(in)::x1,x2
WminusP=Wminus(nn,x1,-x1-x2)
end function WminusP

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- Delta W ---------------------------------------------
function DeltaW(nn,x1,x2)
real(dp)::DeltaW
integer,intent(in)::nn
real(dp),intent(in)::x1,x2
logical::intersect
real(dp)::vMin,vMax,x3

x3=-x1-x2
call LimitsX2(nn,x1,x2,intersect,vMin,vMax)

DeltaW=0._dp
if(intersect) then
    if(abs(x3)>zero) then
        if(x3>0 .and. vMax>0) DeltaW=DeltaW+Integrate_GK(f_1,max(0._dp,vMin),vMax)
        if(x3<0 .and. vMin<0) DeltaW=DeltaW-Integrate_GK(f_1,vMin,min(0._dp,vMax))
        if(x1>0 .and. x3>0) DeltaW=DeltaW+Integrate_GK(f_2,vMin,vMax)
        if(x1<0 .and. x3<0) DeltaW=DeltaW-Integrate_GK(f_2,vMin,vMax)
    else
        if(x1>0 .and. vMin<0) DeltaW=DeltaW+Integrate_GK(f_3,vMin,min(0._dp,vMax))
        if(x1<0 .and. vMax>0) DeltaW=DeltaW-Integrate_GK(f_3,max(0._dp,vMin),vMax)
    end if
else
    DeltaW=0._dp
end if

contains
    function f_1(v)
    real(dp)::f_1
    real(dp), intent(in)::v

    f_1=GETinterpolatorB(nn,x1-v,x2)

    end function f_1

    function f_2(v)
    real(dp)::f_2
    real(dp), intent(in)::v

    f_2=-x1**2*(3*x3+x1)/(x1+x3)**3*GETinterpolatorB(nn,x1-v,x2)

    end function f_2

    function f_3(v)
    real(dp)::f_3
    real(dp), intent(in)::v
        f_3=-GETinterpolatorB(nn,x1-v,x2)
    end function f_3

end function DeltaW

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!----------------------- Delta W P_23 ---------------------------------------------
function DeltaWP(nn,x1,x2)
real(dp)::DeltaWP
integer,intent(in)::nn
real(dp),intent(in)::x1,x2
DeltaWP=DeltaW(nn,x1,-x1-x2)
end function DeltaWP
