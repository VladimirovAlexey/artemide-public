!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!
!!!!!!!!!       part of TMD_AD module for artemide
!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!-------------------------------------------------------------------------------------------------------
!!!--------------------------------------------SECONDARY ADs------------------------------------------------
!!!-------------------------------------------------------------------------------------------------------

!!!! sets the values of logarithm part for  QUARK
!!!! Expression is taken from [1707.07606]
subroutine SetDnkQuark()
    integer::n
    
    !! 1-loop
    do n=NfMIN,NfMAX
        d_nk_Q_internal(1,3,n)=0d0
        d_nk_Q_internal(1,2,n)=0d0
        d_nk_Q_internal(1,1,n)=GammaCusp_q(0,n)/2d0
    end do
    
    !! 2-loop
    do n=NfMIN,NfMAX
        d_nk_Q_internal(2,3,n)=0d0
        d_nk_Q_internal(2,2,n)=GammaCusp_q(0,n)*betaQCD(0,n)/4d0
        d_nk_Q_internal(2,1,n)=GammaCusp_q(1,n)/2d0
    end do
    
    !! 3-loop
    do n=NfMIN,NfMAX
        d_nk_Q_internal(3,3,n)=GammaCusp_q(0,n)*(betaQCD(0,n)**2)/6d0
        d_nk_Q_internal(3,2,n)=(GammaCusp_q(0,n)*betaQCD(1,n)+2d0*GammaCusp_q(1,n)*betaQCD(0,n))/4d0
        d_nk_Q_internal(3,1,n)=(GammaCusp_q(2,n)+4d0*dnk_q(2,0,n)*betaQCD(0,n))/2d0
    end do
    
end subroutine SetDnkQuark

!!!! sets the values of logarithm part for RAD GLUON
!!!! Expression is taken from [1707.07606]
subroutine SetDnkGluon()
    integer::n
    
    !! 1-loop
    do n=NfMIN,NfMAX
        d_nk_G_internal(1,3,n)=0d0
        d_nk_G_internal(1,2,n)=0d0
        d_nk_G_internal(1,1,n)=GammaCusp_g(0,n)/2d0
    end do
    
    !! 2-loop
    do n=NfMIN,NfMAX
        d_nk_G_internal(2,3,n)=0d0
        d_nk_G_internal(2,2,n)=GammaCusp_g(0,n)*betaQCD(0,n)/4d0
        d_nk_G_internal(2,1,n)=GammaCusp_g(1,n)/2d0
    end do
    
    !! 3-loop
    do n=NfMIN,NfMAX
        d_nk_G_internal(3,3,n)=GammaCusp_g(0,n)*(betaQCD(0,n)**2)/6d0
        d_nk_G_internal(3,2,n)=(GammaCusp_g(0,n)*betaQCD(1,n)+2d0*GammaCusp_g(1,n)*betaQCD(0,n))/4d0
        d_nk_G_internal(3,1,n)=(GammaCusp_g(2,n)+4d0*dnk_g(2,0,n)*betaQCD(0,n))/2d0
    end do
    
end subroutine SetDnkGluon

!!!! sets the values of RAD-resummed coefficients for  QUARK
!!!! D=-GAMMA0/2BETA0(LOG(1-X)+as^n/(1-X)^n dnkl X^k Log(1-X)^l
subroutine SetDnklQuark()
    integer::n
    real(dp)::B1,B2,B3,G1,G2,G3
        
    do n=NfMIN,NfMAX
        B1=betaQCD(1,n)/betaQCD(0,n)
        B2=betaQCD(2,n)/betaQCD(0,n)
        B3=betaQCD(3,n)/betaQCD(0,n)
        
        G1=GammaCusp_q(1,n)/GammaCusp_q(0,n)
        G2=GammaCusp_q(2,n)/GammaCusp_q(0,n)
        G3=GammaCusp_q(3,n)/GammaCusp_q(0,n)
        
        !! 1-loop        
        d_nkl_Q_internal(1,0,1,n)=B1
        d_nkl_Q_internal(1,1,0,n)=B1-G1
        d_nkl_Q_internal(1,0,0,n)=0d0
    
        !! 2-loop    
        d_nkl_Q_internal(2,0,2,n)=-(B1**2)/2d0
        d_nkl_Q_internal(2,0,1,n)=B1*G1
        d_nkl_Q_internal(2,2,0,n)=(G2-B1*G1-B2+B1**2)/2d0
        d_nkl_Q_internal(2,1,0,n)=B1*G1-G2
        d_nkl_Q_internal(2,0,0,n)=-2d0*dnk_q(2,0,n)*betaQCD(0,n)/GammaCusp_q(0,n)        
        d_nkl_Q_internal(2,2,2,n)=0d0
        d_nkl_Q_internal(2,1,1,n)=0d0
        
        !! 3-loop    
        d_nkl_Q_internal(3,0,3,n)=(B1**3)/3d0
        d_nkl_Q_internal(3,0,2,n)=-(B1**3)/2d0-(B1**2)*G1
        d_nkl_Q_internal(3,0,1,n)=B1*G2+4d0*dnk_q(2,0,n)*betaQCD(1,n)/GammaCusp_q(0,n)
        d_nkl_Q_internal(3,1,1,n)=B1*B2-B1**3
        d_nkl_Q_internal(3,3,0,n)=1d0/3d0*(B1**3-2d0*B1*B2+B3-(B1**2)*G1+B2*G1+B1*G2-G3)
        d_nkl_Q_internal(3,2,0,n)=-0.5d0*(B1**3)+B1*B2-0.5d0*B3+(B1**2)*G1-B2*G1-B1*G2+G3
        d_nkl_Q_internal(3,1,0,n)=B1*G2-G3
        d_nkl_Q_internal(3,0,0,n)=-2d0*dnk_q(3,0,n)*betaQCD(0,n)/GammaCusp_q(0,n)
        d_nkl_Q_internal(3,1,2,n)=0d0
        d_nkl_Q_internal(3,1,3,n)=0d0
        d_nkl_Q_internal(3,2,1,n)=0d0
        d_nkl_Q_internal(3,2,2,n)=0d0
        d_nkl_Q_internal(3,2,3,n)=0d0
        d_nkl_Q_internal(3,3,1,n)=0d0
        d_nkl_Q_internal(3,3,2,n)=0d0
        d_nkl_Q_internal(3,3,3,n)=0d0
        
    end do
    
end subroutine SetDnklQuark


!!!! sets the values of RAD-resummed coefficients for  QUARK
!!!! D=-GAMMA0/2BETA0(LOG(1-X)+as^n/(1-X)^n dnkl X^k Log(1-X)^l
subroutine SetDnklGluon()
    integer::n
    real(dp)::B1,B2,B3,G1,G2,G3
        
    do n=NfMIN,NfMAX
        B1=betaQCD(1,n)/betaQCD(0,n)
        B2=betaQCD(2,n)/betaQCD(0,n)
        B3=betaQCD(3,n)/betaQCD(0,n)
        
        G1=GammaCusp_g(1,n)/GammaCusp_g(0,n)
        G2=GammaCusp_g(2,n)/GammaCusp_g(0,n)
        G3=GammaCusp_g(3,n)/GammaCusp_g(0,n)
        
        !! 1-loop        
        d_nkl_G_internal(1,0,1,n)=B1
        d_nkl_G_internal(1,1,0,n)=B1-G1
        d_nkl_G_internal(1,0,0,n)=0d0
    
        !! 2-loop    
        d_nkl_G_internal(2,0,2,n)=-(B1**2)/2d0
        d_nkl_G_internal(2,0,1,n)=B1*G1
        d_nkl_G_internal(2,2,0,n)=(G2-B1*G1-B2+B1**2)/2d0
        d_nkl_G_internal(2,1,0,n)=B1*G1-G2
        d_nkl_G_internal(2,0,0,n)=-2d0*dnk_g(2,0,n)*betaQCD(0,n)/GammaCusp_g(0,n)        
        d_nkl_G_internal(2,2,2,n)=0d0
        d_nkl_G_internal(2,1,1,n)=0d0
        
        !! 3-loop    
        d_nkl_G_internal(3,0,3,n)=(B1**3)/3d0
        d_nkl_G_internal(3,0,2,n)=-(B1**3)/2d0-(B1**2)*G1
        d_nkl_G_internal(3,0,1,n)=B1*G2+4d0*dnk_g(2,0,n)*betaQCD(1,n)/GammaCusp_g(0,n)
        d_nkl_G_internal(3,1,1,n)=B1*B2-B1**3
        d_nkl_G_internal(3,3,0,n)=1d0/3d0*(B1**3-2d0*B1*B2+B3-(B1**2)*G1+B2*G1+B1*G2-G3)
        d_nkl_G_internal(3,2,0,n)=-0.5d0*(B1**3)+B1*B2-0.5d0*B3+(B1**2)*G1-B2*G1-B1*G2+G3
        d_nkl_G_internal(3,1,0,n)=B1*G2-G3
        d_nkl_G_internal(3,0,0,n)=-2d0*dnk_g(3,0,n)*betaQCD(0,n)/GammaCusp_g(0,n)
        d_nkl_G_internal(3,1,2,n)=0d0
        d_nkl_G_internal(3,1,3,n)=0d0
        d_nkl_G_internal(3,2,1,n)=0d0
        d_nkl_G_internal(3,2,2,n)=0d0
        d_nkl_G_internal(3,2,3,n)=0d0
        d_nkl_G_internal(3,3,1,n)=0d0
        d_nkl_G_internal(3,3,2,n)=0d0
        d_nkl_G_internal(3,3,3,n)=0d0
        
    end do
    
end subroutine SetDnklGluon

!!!! sets the values of zeta-line PT coefficients for  QUARK
!!!! zeta=C0 mu/b*exp(-v), v=a^n L^k v^{(nk)}
subroutine SetVnkQuark()
    integer::n
    real(dp)::B1,G1,G2,gg1,gg2,gg3,dd2,dd3
        
    do n=NfMIN,NfMAX
        B1=betaQCD(1,n)/betaQCD(0,n)
        
        G1=GammaCusp_q(1,n)/GammaCusp_q(0,n)
        G2=GammaCusp_q(2,n)/GammaCusp_q(0,n)
        
        gg1=gammaV_q(1,n)/GammaCusp_q(0,n)
        gg2=gammaV_q(2,n)/GammaCusp_q(0,n)
        gg3=gammaV_q(3,n)/GammaCusp_q(0,n)
        
        dd2=dnk_q(2,0,n)/GammaCusp_q(0,n)
        dd3=dnk_q(3,0,n)/GammaCusp_q(0,n)
        
        !! 1-loop        
        v_nk_Q_internal(0,1,n)=0d0
        v_nk_Q_internal(0,0,n)=gg1
    
        !! 2-loop    
        v_nk_Q_internal(1,2,n)=betaQCD(0,n)/12d0
        v_nk_Q_internal(1,1,n)=0d0
        v_nk_Q_internal(1,0,n)=dd2-gg1*G1+gg2
        
        !! 3-loop    
        v_nk_Q_internal(2,3,n)=(betaQCD(0,n)**2)/24d0
        v_nk_Q_internal(2,2,n)=betaQCD(0,n)/12d0*(B1+G1)
        v_nk_Q_internal(2,1,n)=betaQCD(0,n)/2d0*(8d0/3d0*dd2-gg1*G1+gg2)
        v_nk_Q_internal(2,0,n)=dd3-dd2*G1+gg1*(G1**2)-gg2*G1-gg1*G2+gg3
        
    end do
    
end subroutine SetVnkQuark

!!!! sets the values of zeta-line PT coefficients for  QUARK
!!!! zeta=C0 mu/b*exp(-v), v=a^n L^k v^{(nk)}
subroutine SetVnkGluon()
    integer::n
    real(dp)::B1,G1,G2,gg1,gg2,gg3,dd2,dd3
        
    do n=NfMIN,NfMAX
        B1=betaQCD(1,n)/betaQCD(0,n)
        
        G1=GammaCusp_g(1,n)/GammaCusp_g(0,n)
        G2=GammaCusp_g(2,n)/GammaCusp_g(0,n)
        
        gg1=gammaV_g(1,n)/GammaCusp_g(0,n)
        gg2=gammaV_g(2,n)/GammaCusp_g(0,n)
        gg3=gammaV_g(3,n)/GammaCusp_g(0,n)
        
        dd2=dnk_g(2,0,n)/GammaCusp_g(0,n)
        dd3=dnk_g(3,0,n)/GammaCusp_g(0,n)
        
        !! 1-loop        
        v_nk_G_internal(0,1,n)=0d0
        v_nk_G_internal(0,0,n)=gg1
    
        !! 2-loop    
        v_nk_G_internal(1,2,n)=betaQCD(0,n)/12d0
        v_nk_G_internal(1,1,n)=0d0
        v_nk_G_internal(1,0,n)=dd2-gg1*G1+gg2
        
        !! 3-loop    
        v_nk_G_internal(2,3,n)=(betaQCD(0,n)**2)/24d0
        v_nk_G_internal(2,2,n)=betaQCD(0,n)/12d0*(B1+G1)
        v_nk_G_internal(2,1,n)=betaQCD(0,n)/2d0*(8d0/3d0*dd2-gg1*G1+gg2)
        v_nk_G_internal(2,0,n)=dd3-dd2*G1+gg1*(G1**2)-gg2*G1-gg1*G2+gg3
        
    end do
    
end subroutine SetVnkGluon

!!!! sets the values of EXACT zeta-line PT coefficients for  QUARK
!!!! zeta=mu^2*exp(-1/beta0/as * OMEGA), OMEGA=a^n OMEGA^{(nk)} ...
!!!! ... in each term different.
subroutine SetOMEGAnkQuark()
    integer::n
    real(dp)::B1,B2,B3,G1,G2,G3,gg1,gg2,gg3,commonF
    
    OMEGA_nk_Q_internal=0d0*OMEGA_nk_Q_internal
        
    do n=NfMIN,NfMAX
        B1=betaQCD(1,n)/betaQCD(0,n)
        B2=betaQCD(2,n)/betaQCD(0,n)
        B3=betaQCD(3,n)/betaQCD(0,n)
        
        G1=GammaCusp_q(1,n)/GammaCusp_q(0,n)
        G2=GammaCusp_q(2,n)/GammaCusp_q(0,n)
        G3=GammaCusp_q(3,n)/GammaCusp_q(0,n)
        
        gg1=gammaV_q(1,n)*betaQCD(0,n)/GammaCusp_q(0,n)
        gg2=gammaV_q(2,n)*betaQCD(0,n)/GammaCusp_q(0,n)
        gg3=gammaV_q(3,n)*betaQCD(0,n)/GammaCusp_q(0,n)
        
        pFACTOR_Q_internal(n)=2d0*betaQCD(0,n)/GammaCusp_q(0,n)
        commonF=1d0/betaQCD(0,n)
        
        !! 0-loop
        !! * z_1
        OMEGA_nk_Q_internal(0,1,n)=1d0*commonF
    
        !! 1-loop
        !! * z_1
        OMEGA_nk_Q_internal(1,1,n)=(B1-G1)*commonF
        !! * 1
        OMEGA_nk_Q_internal(1,2,n)=gg1*commonF
        !! * p
        OMEGA_nk_Q_internal(1,3,n)=-B1/2d0*commonF
        
        !! 2-loop
        !! * z_1
        OMEGA_nk_Q_internal(2,1,n)=(B2-B1*G1+G1**2-G2)/2d0*commonF
        !! * z_{-1}
        OMEGA_nk_Q_internal(2,2,n)=((-B2+B1*G1+G1**2-G2)/2d0-G1*gg1+gg2)*commonF
        !! * 1
        OMEGA_nk_Q_internal(2,3,n)=(-G1*gg1+gg2)*commonF
        
        !! 3-loop
        !! * z_1
        OMEGA_nk_Q_internal(3,1,n)=(2d0*B3+(B1**2)*G1-B2*G1-G1**3+3*G1*G2-B1*B2-B1*G2-2d0*G3)/6d0*commonF
        !! * z_{-1}
        OMEGA_nk_Q_internal(3,2,n)=(-B2+B1*G1+G1**2-G2-2d0*G1*gg1+2d0*gg2)*(B1-G1)/2d0*commonF
        !! * z_{-2}
        OMEGA_nk_Q_internal(3,3,n)=((-B1*B2-B3+(B1**2)*G1+2d0*B2*G1-4d0*(G1**3)-B1*G2+6d0*G1*G2-2d0*G3)/12d0&
                +(2d0*G1-B1)*(G1*gg1-gg2)/2d0-(G2*gg1-gg3)/2d0)*commonF
        !! * 1
        OMEGA_nk_Q_internal(3,4,n)=((G1**2)*gg1-G2*gg1-G1*gg2+gg3)*commonF
    end do
    
end subroutine SetOMEGAnkQuark

!!!! sets the values of EXACT zeta-line PT coefficients for  QUARK
!!!! zeta=mu^2*exp(-1/as * OMEGA), OMEGA=a^n OMEGA^{(nk)} ...
!!!! ... in each term different.
subroutine SetOMEGAnkGluon()
    integer::n
    real(dp)::B1,B2,B3,G1,G2,G3,gg1,gg2,gg3,commonF
    
    OMEGA_nk_G_internal=0d0*OMEGA_nk_G_internal
        
    do n=NfMIN,NfMAX
        B1=betaQCD(1,n)/betaQCD(0,n)
        B2=betaQCD(2,n)/betaQCD(0,n)
        B3=betaQCD(3,n)/betaQCD(0,n)
        
        G1=GammaCusp_g(1,n)/GammaCusp_g(0,n)
        G2=GammaCusp_g(2,n)/GammaCusp_g(0,n)
        G3=GammaCusp_g(3,n)/GammaCusp_g(0,n)
        
        gg1=gammaV_g(1,n)*betaQCD(0,n)/GammaCusp_g(0,n)
        gg2=gammaV_g(2,n)*betaQCD(0,n)/GammaCusp_g(0,n)
        gg3=gammaV_g(3,n)*betaQCD(0,n)/GammaCusp_g(0,n)
        
        pFACTOR_G_internal(n)=2d0*betaQCD(0,n)/GammaCusp_g(0,n)
        commonF=1d0/betaQCD(0,n)
        
        !! 0-loop
        !! * z_1
        OMEGA_nk_G_internal(0,1,n)=1d0*commonF
    
        !! 1-loop
        !! * z_1
        OMEGA_nk_G_internal(1,1,n)=(B1-G1)*commonF
        !! * 1
        OMEGA_nk_G_internal(1,2,n)=gg1*commonF
        !! * p
        OMEGA_nk_G_internal(1,3,n)=-B1/2d0*commonF
        
        !! 2-loop
        !! * z_1
        OMEGA_nk_G_internal(2,1,n)=(B2-B1*G1+G1**2-G2)/2d0*commonF
        !! * z_{-1}
        OMEGA_nk_G_internal(2,2,n)=((-B2+B1*G1+G1**2-G2)/2d0-G1*gg1+gg2)*commonF
        !! * 1
        OMEGA_nk_G_internal(2,3,n)=(-G1*gg1+gg2)*commonF
        
        !! 3-loop
        !! * z_1
        OMEGA_nk_G_internal(3,1,n)=(2d0*B3+(B1**2)*G1-B2*G1-G1**3+3*G1*G2-B1*B2-B1*G2-2d0*G3)/6d0*commonF
        !! * z_{-1}
        OMEGA_nk_G_internal(3,2,n)=(-B2+B1*G1+G1**2-G2-2d0*G1*gg1+2d0*gg2)*(B1-G1)/2d0*commonF
        !! * z_{-2}
        OMEGA_nk_G_internal(3,3,n)=((-B1*B2-B3+(B1**2)*G1+2d0*B2*G1-4d0*(G1**3)-B1*G2+6d0*G1*G2-2d0*G3)/12d0&
                +(2d0*G1-B1)*(G1*gg1-gg2)/2d0-(G2*gg1-gg3)/2d0)*commonF
        !! * 1
        OMEGA_nk_G_internal(3,4,n)=((G1**2)*gg1-G2*gg1-G1*gg2+gg3)*commonF
    end do
    
end subroutine SetOMEGAnkGluon
