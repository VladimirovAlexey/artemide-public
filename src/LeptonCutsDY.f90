!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 1.4
!
!	Evaluation of the leptonic cuts for the DrellYan
!	
!	if you use this module please, quote 1706.01473
!
!	ver 1.0: release (AV, 10.05.2017)
!	ver 1.32 update nessacary for parallelisation (AV,17.09.2018)
!	ver 1.32 CutFactor4 added, asymetric cuts are introduced (AV,03.12.2018)
!	ver.1.4  Deleted old routines. Incapsulated variables (AV. 18.01.2019(
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module LeptonCutsDY
use aTMDe_Numerics
use IO_functions
implicit none
private
! public

!!!Parameters of cut and cross-section
!!!!! this is array =(/  pT1lim,pT2lim,etaMin,etaMax,Exp(2etaMin),exp(2etaMax) /)
real(dp)::cutParam_global(1:6)

!!! number of divisions in Simpsons
integer,parameter::num=64!
!! Tolerance (absolute)
real(dp),parameter::tolerance=0.000001d0

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! the kinematic variables are passed via an aray
!!! (q_T,Q,y, Q^2 , qt^2/Q^2+qt^2 , Sqrt[Q^2+qT^2] )

public:: SetCutParameters
public:: CutFactor4,CutFactorA


  interface SetCutParameters
    module procedure SetCutParameters_sym,SetCutParameters_asym
  end interface

contains


!SEt parameters of cut 
subroutine SetCutParameters_sym(pT_in,etaMin_in,etaMax_in)
  real(dp)::pT_in,etaMax_in,etaMin_in
  
  cutParam_global=(/ pT_in**2,&
	      pT_in**2,&
	      etaMin_in,&
	      etaMax_in,&
	      EXP(2*etaMin_in),&
	      EXP(2*etaMax_in) /)
  
end subroutine SetCutParameters_sym

!SEt parameters of cut 
! with asymetric cuts for pT
subroutine SetCutParameters_asym(pT1_in,pT2_in,etaMin_in,etaMax_in)
  real(dp)::pT1_in,pT2_in,etaMax_in,etaMin_in
  
  !! for definetines we order pt1>pt2
  if(pT1_in>=pT2_in) then 
  cutParam_global=(/ pT1_in**2,&
	      pT2_in**2,&
	      etaMin_in,&
	      etaMax_in,&
	      EXP(2d0*etaMin_in),&
	      EXP(2d0*etaMax_in) /)
  else
  cutParam_global=(/ pT2_in**2,&
	      pT1_in**2,&
	      etaMin_in-tolerance,&
	      etaMax_in+tolerance,&
	      EXP(2d0*etaMin_in),&
	      EXP(2d0*etaMax_in) /)
  end if
  
end subroutine SetCutParameters_asym

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!
!!	This function survives after many modification of module. For different version of integral evaluation see /history
!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! Integration is done by more accurate method. Secting the phi-plane, and define the boundaries of eta
!!!! The integal over eta is done exactly
!!!! Unfortunately works only for large enough Q.
!!!! Gives very precise result 10^-5 accuracy
!!!!
!!!! CutParameters = (/ kT1, kT2, etaMin, etaMax /)
function CutFactor4(qT,Q_in,y_in,CutParameters)
  real(dp):: qT,CutFactor4,Q_in,y_in
  real(dp):: dphi
  integer :: i
  real(dp),dimension(1:6)::var,varC
  real(dp),dimension(1:4),optional,intent(in)::CutParameters
  
  if(present(CutParameters)) then
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
  else 
    varC=cutParam_global
  end if
  
  var=(/qT,Q_in,y_in,Q_in**2,qT**2/(Q_in**2+qT**2),SQRT(Q_in**2+qT**2)/)
  
  if(varC(3)<varC(4)) then
  if(y_in<varC(3) .or. y_in>varC(4)) then
    CutFactor4=0d0
  else
    CutFactor4=0d0
    dphi=pix2/num
    do i=0,num
      CutFactor4=CutFactor4+Simp(i,num)*IntegralOverEtaFixedPhiEXACT(var,varC,i*dphi)
    end do
      CutFactor4=CutFactor4*dphi/3d0/2d0
  end if
  else
    CutFactor4=0d0
  end if
  
end function CutFactor4

!!!! integral for anti-symmetric part of lepton tensor
function CutFactorA(qT,Q_in,y_in,CutParameters)
  real(dp):: qT,CutFactorA,Q_in,y_in
  real(dp):: dphi
  integer :: i
  real(dp),dimension(1:6)::var,varC
  real(dp),dimension(1:4),optional,intent(in)::CutParameters
  
  if(present(CutParameters)) then
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
  else 
    varC=cutParam_global
  end if
  
  var=(/qT,Q_in,y_in,Q_in**2,qT**2/(Q_in**2+qT**2),SQRT(Q_in**2+qT**2)/)
  
  if(varC(3)<varC(4)) then
  if(y_in<varC(3) .or. y_in>varC(4)) then
    CutFactorA=0d0
  else
    CutFactorA=0d0
    dphi=pix2/num
    do i=0,num
      CutFactorA=CutFactorA+Simp(i,num)*IntegralOverEtaFixedPhiEXACT_A(var,varC,i*dphi)
    end do
      CutFactorA=CutFactorA*dphi/3d0/2d0
  end if
  else
    CutFactorA=0d0
  end if
  
end function CutFactorA



!it is the same function as IntegralOverEtaFixedPhi, but the integration is done exactly.
! the boundaries are defined as before
!	   1     2    3       4          5                       6
! var =(/ qT,  Q_in, y_in,   Q_in**2,  qT**2/(Q_in**2+qT**2), SQRT(Q_in**2+qT**2)/)
! varC=(/ pT1, pT2,  etaMin, etaMax,   EXP(2*etaMin),         EXP(2*etaMax) /)
function IntegralOverEtaFixedPhiEXACT(var,varC,phi)
  real(dp),dimension(1:6)::var,varC
  real(dp):: IntegralOverEtaFixedPhiEXACT,phi
  real(dp)::eta1,eta2
  
  if(Integrand2THETA(var,varC,var(3),phi)==0) then
  !!!! Here is the point of assumption!! because if Q~qT the integration region is difficult and the boundaries may be inproper
  IntegralOverEtaFixedPhiEXACT=0d0
  else
  !!lower boundary
  eta1=FindBoundary(var,varC,varC(3)-0.0001d0,var(3),phi)
  !!upper boundary
  eta2=FindBoundary(var,varC,var(3),varC(4)+0.0001d0,phi)
  !write(*,*) eta1,eta2
  
  IntegralOverEtaFixedPhiEXACT=var(4)/(16d0*pi)*(&
	integralEtaExactUNDEFINED(var(1)*cos(phi),var(6),var(3)-eta2,var(4))&
	-integralEtaExactUNDEFINED(var(1)*cos(phi),var(6),var(3)-eta1,var(4)))
	
  end if
end function IntegralOverEtaFixedPhiEXACT

!it is the same function as IntegralOverEtaFixedPhi, but the integration is done exactly.
! the boundaries are defined as before
!	   1     2    3       4          5                       6
! var =(/ qT,  Q_in, y_in,   Q_in**2,  qT**2/(Q_in**2+qT**2), SQRT(Q_in**2+qT**2)/)
! varC=(/ pT1, pT2,  etaMin, etaMax,   EXP(2*etaMin),         EXP(2*etaMax) /)
!!! ANTI SYMMETRIC ONE!
function IntegralOverEtaFixedPhiEXACT_A(var,varC,phi)
  real(dp),dimension(1:6)::var,varC
  real(dp):: IntegralOverEtaFixedPhiEXACT_A,phi
  real(dp)::eta1,eta2
  
  if(Integrand2THETA(var,varC,var(3),phi)==0) then
  !!!! Here is the point of assumption!! because if Q~qT the integration region is difficult and the boundaries may be inproper
  IntegralOverEtaFixedPhiEXACT_A=0d0
  else
  !!lower boundary
  eta1=FindBoundary(var,varC,varC(3)-0.0001d0,var(3),phi)
  !!upper boundary
  eta2=FindBoundary(var,varC,var(3),varC(4)+0.0001d0,phi)
  !write(*,*) eta1,eta2
  
  IntegralOverEtaFixedPhiEXACT_A=var(4)/(16d0*pi)*(&
	integralEtaExactUNDEFINED_A(var(1)*cos(phi),var(6),var(3)-eta2,var(4))&
	-integralEtaExactUNDEFINED_A(var(1)*cos(phi),var(6),var(3)-eta1,var(4)))
	
  end if
end function IntegralOverEtaFixedPhiEXACT_A


!! this is integral over eta exactly evaluated by mathematica. Udefined. (without common factor Q^2/16/pi)
!! it is used in IntegralOverEtaFixedPhiEXACT
!! a = qT cos phi
!! b = Sqrt{Q^2+qT^2}
!! uu = y-eta
function integralEtaExactUNDEFINED(a,b,uu,Q2)
  real(dp)::a,b,uu,integralEtaExactUNDEFINED,Q2
  real(dp)::w,R,TT,bb
  
  bb=b*b
  w=bb-a**2
  R=a-b*Cosh(uu)
  TT=-b*Sinh(uu)/w/R
  
  integralEtaExactUNDEFINED=2d0*Q2*TT/R**2+a*(6d0*w-5d0*Q2)*TT/w/R&
      +(6d0-(18d0*bb+11d0*Q2)/w+15d0*bb*Q2/w**2)*TT&
      -6d0*a*(2d0*(3d0*bb+Q2)*w-5d0*bb*Q2)*atan((a+b)/sqrt(w)*tanh(uu/2d0))/w**(3.5d0)
  
end function integralEtaExactUNDEFINED


!! this is integral over eta exactly evaluated by mathematica. Udefined. (without common factor Q^2/16/pi)
!! it is used in IntegralOverEtaFixedPhiEXACT
!! a = qT cos phi
!! b = Sqrt{Q^2+qT^2}
!! uu = y-eta
!!! ANTI SYMMETRIC ONE!
function integralEtaExactUNDEFINED_A(a,b,uu,Q2)
  real(dp)::a,b,uu,integralEtaExactUNDEFINED_A,Q2
  
  
  integralEtaExactUNDEFINED_A=-6d0/(a-b*Cosh(uu))**2
  
end function integralEtaExactUNDEFINED_A


!!!!the theta (0 or 1) funciton of the integrand in the coordinates rapidity-angle
!! Search for the boundary of integration between t1 and t2 at fixed phi
!! by devision on 2
!	   1     2    3       4          5                       6
! var =(/ qT,  Q_in, y_in,   Q_in**2,  qT**2/(Q_in**2+qT**2), SQRT(Q_in**2+qT**2)/)
! varC=(/ pT1, pT2,  etaMin, etaMax,   EXP(2*etaMin),         EXP(2*etaMax) /)
function Integrand2THETA(var,varC,h1,p1)
  real(dp),dimension(1:6)::var,varC
  real(dp):: h1,p1
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

function Simp(i,n)
  integer::i,n
  real(dp)::Simp
  if((i==0).or.(i==n)) then
  Simp=1d0
  else 
    if(MOD(i,2)==1) then
    Simp=4d0
    else
    Simp=2d0
    end if
  end if
end function Simp


!! Search for the boundary of integration between t1 and t2 at fixed phi
!! by devision on 2
!	   1     2    3       4          5                       6
! var =(/ qT,  Q_in, y_in,   Q_in**2,  qT**2/(Q_in**2+qT**2), SQRT(Q_in**2+qT**2)/)
! varC=(/ pT1, pT2,  etaMin, etaMax,   EXP(2*etaMin),         EXP(2*etaMax) /)
function FindBoundary(var,varC,eta1_in,eta2_in,phi)
  real(dp),dimension(1:6)::var,varC
  real(dp)::FindBoundary,eta1,eta2,phi,eta3,eta1_in,eta2_in
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
      if (ABS(eta1-eta2)<tolerance) exit
    end do 
    FindBoundary=(eta1+eta2)/2d0
  end if
end function FindBoundary

end module LeptonCutsDY
