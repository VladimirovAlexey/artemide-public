!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! List of functions to substitute into the cross-section
!!! These functions contain kernels that multiply the TMD distributions in the convolution
!!!
!!! Q2, tau2, qT2 external variables
!!! S, Lambda -- integration variables
!!! cT=cos(theta)
!!! process -- integer that specifies process

function SIDIS_KERNEL(Q2, tau2, qT2, S, Lam, process)
real(dp), intent(in) :: Q2, tau2, qT2, S, Lam
integer, intent(in) :: process
real(dp) :: SIDIS_KERNEL

SELECT CASE(process)


	CASE (2001) !F_UU,T f D-part

		SIDIS_KERNEL = 1+tau2**2/2/Q2**2*((S+sqrt(Lam))**2-1)

	CASE (2011) !F_UU,L f D-part

		SIDIS_KERNEL = tau2**2/2/Q2**2*((S+sqrt(Lam))**2-1)

	CASE (2021) !F_UU^cos(phi_h) f D-part

		SIDIS_KERNEL = tau2/2/sqrt(qT2)/Q2**(5./2)*(2*Q2*(qT2+tau2*S*(S+sqrt(Lam)))+qT2*tau2*((S+sqrt(Lam))**2-1))

	CASE (2031) !F_UU^cos(2phi_h) f D-part

		SIDIS_KERNEL = (4*Q2**2*tau2*(qT2+S**2*tau2)+2*Q2*qT2*tau2*(2*qT2+tau2*(1-Lam+S**2))+qT2*tau2**2*((S+sqrt(Lam))**2-1))/4/Q2**2/qT2

		!------------------------------------------------
	CASE (2004) !F_UU,T h H-part

		SIDIS_KERNEL = tau2**2/8/Q2**2*(2*qT2+tau2*(1-Lam+S**2))*((S+sqrt(Lam))**2-1)/M2

	CASE (2014) !F_UU,L h H-part

		SIDIS_KERNEL = tau2**2/8/Q2**2*(2*qT2+tau2*(1-Lam+S**2))*((S+sqrt(Lam))**2-1)/M2

	CASE (2024) !F_UU^cos(phi_h) h H-part

		SIDIS_KERNEL = tau2**2*(2*Q2*tau2*S*(S-sqrt(Lam))+2*qT2**2+qT2*tau2*(1-Lam+S**2))*((S+sqrt(Lam))**2-1)/8/M2/Q2**(5./2)/sqrt(qT2)

	CASE (2034) !F_UU^cos(2phi_h) h H-part ! TO BE UPDATED
			SIDIS_KERNEL = 1._dp


! 	! To check convolution integral with gaussians
! 	CASE (999)
!
! 		SIDIS_KERNEL = 1


	CASE DEFAULT
		write(*,*) ErrorString('undefined process 2 variables: ',moduleName),process
		write(*,*) color('Evaluation stop',c_red_bold)
		stop

END SELECT

end function SIDIS_KERNEL


