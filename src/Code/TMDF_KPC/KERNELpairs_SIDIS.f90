!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! List of functions to substitute into the cross-section
!!! These functions contain kernels that multiply the TMD distributions in the convolution
!!!
!!! Q2, tau2, qT2 external variables
!!! S, Lambda -- integration variables
!!! cT=cos(theta)
!!! process -- integer that specifies process

function SIDIS_KERNEL(Q2, tau2, qT2, vareps, S, Lam, process)
real(dp), intent(in) :: Q2, tau2, qT2, vareps, S, Lam
integer, intent(in) :: process
real(dp) :: SIDIS_KERNEL

SELECT CASE(process)


    CASE (2001,2002,2003,2101,2102,2103,2104,2105,2106,2107,2108,2109) !F_UU,T f D-part

        SIDIS_KERNEL = 1+tau2**2/2/Q2**2*((S+sqrt(Lam))**2-1)

    CASE (2011,2012,2013,2111,2112,2113,2114,2115,2116,2117,2118,2119) !F_UU,L f D-part

        SIDIS_KERNEL = tau2**2/2/Q2**2*((S+sqrt(Lam))**2-1)

    CASE (2021,2022,2023,2121,2122,2123,2124,2125,2126,2127,2128,2129) !F_UU^cos(phi_h) f D-part

        SIDIS_KERNEL = tau2/2/sqrt(qT2)/Q2**(5./2)*(2*Q2*(qT2+tau2*S*(S+sqrt(Lam)))+qT2*tau2*((S+sqrt(Lam))**2-1))

    CASE (2031,2032,2033,2131,2132,2133,2134,2135,2136,2137,2138,2139) !F_UU^cos(2phi_h) f D-part

        SIDIS_KERNEL = (4*Q2**2*tau2*(qT2+S**2*tau2) &
				+2*Q2*qT2*tau2*(2*qT2+tau2*(1-Lam+S**2))+qT2*tau2**2*((S+sqrt(Lam))**2-1))/4/Q2**2/qT2

    CASE (2041,2042,2043,2141,2142,2143,2144,2145,2146,2147,2148,2149) !F_{UU,T}+\vareps F_{UU,L} f D-part

		SIDIS_KERNEL = 1+(1+2*vareps)*tau2**2/2/Q2**2*((S+sqrt(Lam))**2-1)

        !------------------------------------------------
    CASE (3001,3002,3003,3101,3102,3103,3104,3105,3106,3107,3108,3109) !F_UU,T h H-part

        SIDIS_KERNEL = tau2**2/8/Q2**2*(2*qT2+tau2*(1-Lam+S**2))*((S+sqrt(Lam))**2-1)/M2

    CASE (3011,3012,3013,3111,3112,3113,3114,3115,3116,3117,3118,3119) !F_UU,L h H-part

        SIDIS_KERNEL = tau2**2/8/Q2**2*(2*qT2+tau2*(1-Lam+S**2))*((S+sqrt(Lam))**2-1)/M2

    CASE (3021,3022,3023,3121,3122,3123,3124,3125,3126,3127,3128,3129) !F_UU^cos(phi_h) h H-part

        SIDIS_KERNEL = tau2**2*(2*Q2*tau2*S*(S-sqrt(Lam))&
				+2*qT2**2+qT2*tau2*(1-Lam+S**2))*((S+sqrt(Lam))**2-1)/8/M2/Q2**(5./2)/sqrt(qT2)

    CASE (3031,3032,3033,3131,3132,3133,3134,3135,3136,3137,3138,3139) !F_UU^cos(2phi_h) h H-part ! TO BE UPDATED
        write(*,*) ErrorString('update process 3031: ',moduleName),process
        write(*,*) color('Evaluation stop',c_red_bold)
        stop
            SIDIS_KERNEL = 1._dp

	CASE (3041,3042,3043,3141,3142,3143,3144,3145,3146,3147,3148,3149) !F_UU,T h H-part

        SIDIS_KERNEL = (1+2*vareps)*tau2**2/8/Q2**2*(2*qT2+tau2*(1-Lam+S**2))*((S+sqrt(Lam))**2-1)/M2


!     ! To check convolution integral with gaussians
!     CASE (999)
!
!         SIDIS_KERNEL = 1


    CASE DEFAULT
        write(*,*) ErrorString('undefined process 2 variables: ',moduleName),process
        write(*,*) color('Evaluation stop',c_red_bold)
        stop

END SELECT

end function SIDIS_KERNEL


