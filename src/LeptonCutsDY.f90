!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!         arTeMiDe 3.0
!
!   Evaluation of the leptonic cuts for the DrellYan
!
!   if you use this module please, quote 1706.01473
!
!   ver 1.0: release (AV, 10.05.2017)
!   ver 1.32 update nessacary for parallelisation (AV,17.09.2018)
!   ver 1.32 CutFactor4 added, asymetric cuts are introduced (AV,03.12.2018)
!   ver.1.4  Deleted old routines. Incapsulated variables (AV. 18.01.2019)
!   ver.3.0  Massive update to the artemide 3.0. Inclusion angular factors + update of integration routines. (AV+SPA. 9.04.2024)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module LeptonCutsDY
use aTMDe_Numerics
use IO_functions
use IntegrationRoutines
implicit none
private

!! Tolerance (absolute)
real(dp)::toleranceINT=0.000001d0
!! Tolerance (comparison)
real(dp)::zero=1.d-8

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! the kinematic variables are passed via an aray
!!! (q_T,Q,y, Q^2 , qt^2/Q^2+qt^2 , Sqrt[Q^2+qT^2] )
!!! The parameters of cut is storred in the array
!!!!! this is array =(/  pT1lim,pT2lim,etaMin,etaMax,Exp(2etaMin),exp(2etaMax) /)

public:: CutFactor,InitializeLeptonCutDY


contains

!!!! tolerance is for the integration tolerance (if adaptive case used)
!!!! zero is for the comparisons
subroutine InitializeLeptonCutDY(tolerance_in,zero_in)
real(dp),intent(in)::tolerance_in,zero_in
toleranceINT=tolerance_in
zero=zero_in
end subroutine InitializeLeptonCutDY

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!   This function survives after many modification of module. For different version of integral evaluation see /history
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Integration is done by more accurate method. Secting the phi-plane, and define the boundaries of eta
!!!! The integal over eta is done exactly
!!!! Unfortunately works only for large enough Q.
!!!! Gives very precise result 10^-5 accuracy
!!!!
!!!! CutParameters = (/ kT1, kT2, etaMin, etaMax /)
!!!! Cut_Type is the selector for the type of cut factor
function CutFactor(qT_in,Q_in,y_in,CutParameters,Cut_Type)
  real(dp),intent(in)::qT_in,Q_in,y_in
  integer,intent(in)::Cut_Type
  real(dp),dimension(1:4),intent(in)::CutParameters
  real(dp):: CutFactor

  real(dp),dimension(1:6)::var,varC
  real(dp)::qT

  
  if(CutParameters(1)>CutParameters(2)) then
    varC(1)=CutParameters(1)**2
    varC(2)=CutParameters(2)**2
  else
    varC(1)=CutParameters(2)**2
    varC(2)=CutParameters(1)**2
  end if
  varC(3)=CutParameters(3)
  varC(4)=CutParameters(4)
  varC(5)=exp(2d0*CutParameters(3))
  varC(6)=exp(2d0*CutParameters(4))

  if(qT_in<0.01d0) then
    qT=0.01d0
  else
    qT=qT_in
  end if
  ! Array containing the kinematic variables, var = [q_T, Q, y, Q^2 , qT^2/(Q^2+qT^2), Sqrt[Q^2+qT^2]
  var=(/qT,Q_in,y_in,Q_in**2,qT**2/(Q_in**2+qT**2),SQRT(Q_in**2+qT**2)/)
  
  if(varC(3)<varC(4)) then
    if(y_in<varC(3) .or. y_in>varC(4)) then
      CutFactor=0d0
    else
      !!!! 7-point estimation
      !CutFactor=Integrate_G7(IntegrandOverPHI,0._dp,pix2)/2
      !!!!! very accurate (adaptive)
      CutFactor=Integrate_GK(IntegrandOverPHI,0._dp,pix2,toleranceINT)
    end if
  else
    CutFactor=0d0
  end if

contains

function IntegrandOverPHI(phi)
real(dp)::IntegrandOverPHI
real(dp),intent(in)::phi
IntegrandOverPHI=IntegralOverEtaFixedPhiEXACT(var,varC,phi,Cut_Type)

end function IntegrandOverPHI
  
end function CutFactor

!it is the same function as IntegralOverEtaFixedPhi, but the integration is done exactly.
! the boundaries are defined as before
!      1     2    3       4          5                       6
! var =(/ qT,  Q_in, y_in,   Q_in**2,  qT**2/(Q_in**2+qT**2), SQRT(Q_in**2+qT**2)/)
! varC=(/ pT1, pT2,  etaMin, etaMax,   EXP(2*etaMin),         EXP(2*etaMax) /)
!!!! Cut_Type is the selector for the type of cut factor
function IntegralOverEtaFixedPhiEXACT(var,varC,phi,Cut_Type)
  real(dp),intent(in),dimension(1:6)::var,varC
  real(dp),intent(in)::phi
  integer,intent(in)::Cut_Type
  real(dp):: IntegralOverEtaFixedPhiEXACT

  real(dp)::eta1,eta2,par1,par2,par3,delta,del2
  
  if(Integrand2THETA(var,varC,var(3),phi)==0) then
    !!!! Here is the point of assumption!! because if Q~qT the integration region is difficult and the boundaries may be inproper
    IntegralOverEtaFixedPhiEXACT=0d0
  else
    !!lower boundary
    eta1=FindBoundary(var,varC,varC(3)-toleranceINT,var(3),phi)
    !!upper boundary
    eta2=FindBoundary(var,varC,var(3),varC(4)+toleranceINT,phi)
    !write(*,*) eta1,eta2

    delta=var(1)/var(2)
    par1=delta*cos(phi)
    del2=delta**2
    par2=var(6)/var(2)
    par3=sin(phi)**2

    IntegralOverEtaFixedPhiEXACT=3/(16*pi)*&
      (integralEtaExactUNDEFINED(phi,par1,par2,par3,var(3)-eta2,del2,Cut_Type)&
      -integralEtaExactUNDEFINED(phi,par1,par2,par3,var(3)-eta1,del2,Cut_Type))
  end if
end function IntegralOverEtaFixedPhiEXACT

!! this is integral over eta exactly evaluated by mathematica. Udefined.
!! Without the factor 3/16 pi
!! it is used in IntegralOverEtaFixedPhiEXACT
!! a = delta cos phi
!! b = Sqrt{1+delta^2}
!! c = sin^2 phi
!! del2=delta**2
!! uu = y-eta
function integralEtaExactUNDEFINED(phi,a,b,c,uu,del2,Cut_Type)
real(dp),intent(in)::phi,a,b,c,uu,del2
integer,intent(in)::Cut_Type
real(dp)::integralEtaExactUNDEFINED

real(dp)::w,R

w=b*b-a*a
R=b*Cosh(uu)-a

SELECT CASE(Cut_Type)
  CASE(-2)
  !!!! usual LP cut factor

  integralEtaExactUNDEFINED=-a*atan((a+b)/sqrt(w)*tanh(uu/2))/w**3.5*(2*w+(6*w-5)*(1+del2))&
           -(w*(R**2*(11 - 6*w) - 2*w + a*R*(-5 + 6*w)) + 3*R**2*(-5 + 6*w)*(1 + del2))*b*sinh(uu)/6/(w*R)**3

  CASE(-1)
  !!!! Puu cut factor
  integralEtaExactUNDEFINED=-a*atan((a+b)/sqrt(w)*tanh(uu/2))/w**2.5*(1+2*w)&
            -(a*R*w - 2*w**2 + R**2*(3 + 5*w + 3*(1 + 2*w)*del2))*sinh(uu)/6/b/w**2/R**3

  CASE(0)
  !!!! P0 factor
  integralEtaExactUNDEFINED=a*atan((a+b)/sqrt(w)*tanh(uu/2))/w**2.5/2*(3-2*w)&
            +(a*R*w - 2*w**2 + R**2*(3 - 3*c - 2*w)*del2)*sinh(uu)/4/b/R**3/w**2

  CASE(1)
  !!!! P1 factor

  integralEtaExactUNDEFINED=-(2*cos(phi)-3*R*sqrt(del2))/3/b**2/R**3

  CASE(2)
  !!!! P2 factor
  integralEtaExactUNDEFINED=-a*atan((a+b)/sqrt(w)*tanh(uu/2))/2/w**3.5*(10 + w*(-11 + 2*w) + 2*c*(-5 + 2*w))&
            -(-30 - 90*del2 + 33*w - 2*w**2*(1 + 2*w) + 2*c*del2*R**2*(-15 + 11*w) + b**6*(30 - 33*w + 6*w**2) &
            - 3*del2**3*(10 + w*(-11 + 2*w)) + 3*del2**2*(-3 + R**2)*(10 + w*(-11 + 2*w)) &
            + del2*(a*R*(10 - 10*c - 11*w)*w + w*(99 - 2*w*(7 + w)) + R**2*(60 + w*(-85 + 29*w)))&
            )*sinh(uu)/12/b/(R*w)**3/del2

  CASE(3)
  !!!! P3 factor
  integralEtaExactUNDEFINED=-(a*b*atan((a+b)/sqrt(w)*tanh(uu/2))/w**2.5*(2*w-3)&
            +(a**3 - a*(1 + del2) + del2*R*(-3 + 3*c + 2*w))*sinh(uu)/2/(w*R)**2)/sqrt(del2)

  CASE(4)
  !!!! P4 factor
  integralEtaExactUNDEFINED=1._dp/2/b/R**2

  CASE(5)
  !!!! P5 factor
  integralEtaExactUNDEFINED=sin(phi)/sqrt(del2)*(&
            2*del2*b*atan((a+b)/sqrt(w)*tanh(uu/2))/w**3.5_dp*(-5+3*w+c*(5-2*w))&
            +(-2*a**5 + 4*a**3*(1 + del2) - a*(1 + del2)*(2 + 2*del2 + 15*R**2) &
            + a*(11 + 9*del2)*R**2*w + del2*R*w*(-5 + 5*c + 3*w))&
            *sinh(uu)/3/(w*R)**3)

  CASE(6)
  !!!! P6 factor
  integralEtaExactUNDEFINED=-2._dp*sin(phi)/3/b/R**3

  CASE(7)
  !!!! P3A=P8 factor
  integralEtaExactUNDEFINED=-sin(phi)*(&
            atan((a+b)/sqrt(w)*tanh(uu/2))/w**2.5_dp*(3-2*w+3*del2)&
            +(3*a*R+w)*b*sinh(uu)/2/(w*R)**2)

  CASE DEFAULT
  !!!! Puu factor
  integralEtaExactUNDEFINED=-a*atan((a+b)/sqrt(w)*tanh(uu/2))/w**2.5*(1+2*w)&
            -(a*R*w - 2*w**2 + R**2*(3 + 5*w + 3*(1 + 2*w)*del2))*sinh(uu)/6/b/w**2/R**3
END SELECT


end function integralEtaExactUNDEFINED


!!!!the theta (0 or 1) function of the integrand in the coordinates rapidity-angle
!! Search for the boundary of integration between t1 and t2 at fixed phi
!! by devision on 2
!      1     2    3       4          5                       6
! var =(/ qT,  Q_in, y_in,   Q_in**2,  qT**2/(Q_in**2+qT**2), SQRT(Q_in**2+qT**2)/)
! varC=(/ pT1, pT2,  etaMin, etaMax,   EXP(2*etaMin),         EXP(2*etaMax) /)
function Integrand2THETA(var,varC,h1,p1)
  real(dp),dimension(1:6),intent(in)::var,varC
  real(dp),intent(in):: h1,p1
  integer::Integrand2THETA

  real(dp)::cosp1,chhy,l1square,l2square,exp2h2,l1
  
  cosp1=COS(p1)
  chhy=COSH(h1-var(3))
  
  l1=var(4)/2d0/(var(6)*chhy-var(1)*cosp1)
  l1square=l1**2
  
  l2square=var(1)**2+l1square-2d0*var(1)*l1*cosp1
  exp2h2=(EXP(2d0*var(3)+h1)*var(6)-EXP(var(3)+2d0*h1)*l1)/(EXP(h1)*var(6)-EXP(var(3))*l1)
  
    
  If((l1square>varC(1)).and.(l2square>varC(2)).and.&
    (varC(4)>h1).and.(varC(6)>exp2h2).and.&
    (h1>varC(3)).and.(exp2h2>varC(5))) then
    Integrand2THETA=1
    else
    Integrand2THETA=0
   end if  
end function Integrand2THETA

!! Search for the boundary of integration between t1 and t2 at fixed phi
!! by devision on 2
!      1     2    3       4          5                       6
! var =(/ qT,  Q_in, y_in,   Q_in**2,  qT**2/(Q_in**2+qT**2), SQRT(Q_in**2+qT**2)/)
! varC=(/ pT1, pT2,  etaMin, etaMax,   EXP(2*etaMin),         EXP(2*etaMax) /)
function FindBoundary(var,varC,eta1_in,eta2_in,phi)
  real(dp),dimension(1:6),intent(in)::var,varC
  real(dp),intent(in)::eta1_in,eta2_in
  real(dp)::FindBoundary

  real(dp)::eta1,eta2,phi,eta3
  integer:: v1,v2,v3,i
  
  eta1=eta1_in
  eta2=eta2_in
  v1=Integrand2THETA(var,varC,eta1,phi)
  v2=Integrand2THETA(var,varC,eta2,phi)
  
  if(v1==v2) then
    write(*,*) 'ERROR LeptonCuts: problem with boundary evaluation. Probably Q is too close to qT. EVALUATION STOP'
    write(*,*) 'Problematic values:'
    write(*,*) 'Kinematic:  Q=',var(2), 'qT=',var(1),'y=',var(3)
    write(*,*) 'Cut params: (pT1,pT2)=(',SQRT(varC(1)),SQRT(varC(2)),'), eta =(',varC(3),varC(4),')'
    stop
    FindBoundary=(eta1+eta2)/2d0
  else
    do 
      eta3=(eta1+eta2)/2d0
      v3=Integrand2THETA(var,varC,eta3,phi)
      if(v3==v1) then 
        eta1=eta3
      else
        eta2=eta3
      end if
      i=i+1
      if (ABS(eta1-eta2)<toleranceINT) exit
    end do 
    FindBoundary=(eta1+eta2)/2d0
  end if
end function FindBoundary

end module LeptonCutsDY
