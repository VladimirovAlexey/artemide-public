!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! List of functions to substitute into the cross-section
!!! These functions contain kernels that multiply the TMD distributions in the convolution
!!!
!!! Q2, tau2, qT2 external variables
!!! S, Lambda -- integration variables
!!! sA=sin(alpha)
!!! cT=cos(theta)
!!! process -- integer that specifies process
function DY_KERNEL(Q2,tau2,qT2,S,Lam,sA,cT,process)
real(dp),intent(in)::Q2,tau2,qT2,S,Lam,sA,cT
integer,intent(in)::process
real(dp)::DY_KERNEL


SELECT CASE(process)

  CASE (1,2,3,101,102,103) !unpol. DY, EW-procces

    DY_KERNEL=1._dp

  !---------------------------------
  CASE (20) !A_0 ff-part

    DY_KERNEL=sA*sA!1-tau2*Lam/Q2

  CASE (21) !A_1 ff-part

    DY_KERNEL=tau2*S*sqrt(Lam/Q2/qT2)

  CASE (22) !A_2 ff-part

    DY_KERNEL=(Q2**2+Lam*tau2**2-Q2*tau2*(1-2*S**2+Lam))/Q2/qT2

  CASE (23) !A_3 ff-part

    DY_KERNEL=2*S*sqrt(tau2/qT2)

  CASE (24) !A_4 ff-part

    DY_KERNEL=2*sqrt(Lam*tau2/Q2)

  !---------------------------------
  CASE (29) !Sigma_UU+Sigma_0/2 ff-part (for normalization of nu)

    DY_KERNEL=1.5_dp-tau2*Lam/Q2/2

  CASE (30) !A_0 hh-part

    DY_KERNEL=tau2/(4*M2*Q2)*(Q2*(1-S**2-Lam)-(1+S**2-Lam)*Lam*tau2)

  CASE (31) !A_1 hh-part

    DY_KERNEL=tau2*S*sqrt(Lam)*(2*Q2-(1-S**2+Lam)*tau2)/(4*M2*sqrt(Q2*qT2))

  CASE (32) !A_2 hh-part

    DY_KERNEL=tau2/(4*M2*Q2*qT2)*(Q2**2*(-1+S**2+Lam)+tau2**2*Lam*(1+S**2-Lam)+Q2*tau2*(2*S**4+(Lam-1)**2-3*S**2*(1+Lam)))

  CASE (35) !A_5 hh-part

    DY_KERNEL=-tau2/(4*M2*qT2)*sqrt(Lam*tau2/Q2)*(Q2*(1+S**2-Lam)-(1-S**2-Lam)*tau2)

  CASE (36) !A_6 hh-part

    DY_KERNEL=-tau2*S/(4*M2*Q2)*sqrt(tau2/qT2)*(Q2*(-1+S**2-Lam)+2*Lam*tau2)

  CASE DEFAULT
    write(*,*) ErrorString('undefined process 2 variables: ',moduleName),process
    write(*,*) color('Evaluation stop',c_red_bold)
    stop
 END SELECT


end function DY_KERNEL
