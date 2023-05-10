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
    !! *[checked 18.05.22]
    do n=NfMIN,NfMAX
        d_nk_Q_internal(1,4,n)=0d0
        d_nk_Q_internal(1,3,n)=0d0
        d_nk_Q_internal(1,2,n)=0d0
        d_nk_Q_internal(1,1,n)=GammaCusp_q(0,n)/2d0
    end do
    
    !! 2-loop
    !! *[checked 18.05.22]
    do n=NfMIN,NfMAX
        d_nk_Q_internal(2,4,n)=0d0
        d_nk_Q_internal(2,3,n)=0d0
        d_nk_Q_internal(2,2,n)=GammaCusp_q(0,n)*betaQCD(0,n)/4d0
        d_nk_Q_internal(2,1,n)=GammaCusp_q(1,n)/2d0
    end do
    
    !! 3-loop
    !! *[checked 18.05.22]
    do n=NfMIN,NfMAX
        d_nk_Q_internal(3,4,n)=0
        d_nk_Q_internal(3,3,n)=GammaCusp_q(0,n)*(betaQCD(0,n)**2)/6d0
        d_nk_Q_internal(3,2,n)=(GammaCusp_q(0,n)*betaQCD(1,n)+2d0*GammaCusp_q(1,n)*betaQCD(0,n))/4d0
        d_nk_Q_internal(3,1,n)=(GammaCusp_q(2,n)+4d0*dnk_q(2,0,n)*betaQCD(0,n))/2d0
    end do
    
    !! 4-loop
    do n=NfMIN,NfMAX
        d_nk_Q_internal(4,4,n)=GammaCusp_q(0,n)*(betaQCD(0,n)**3)/8d0
        d_nk_Q_internal(4,3,n)=5d0/12d0*GammaCusp_q(0,n)*betaQCD(0,n)*betaQCD(1,n)&
                                +GammaCusp_q(1,n)*(betaQCD(0,n)**2)/2d0
        d_nk_Q_internal(4,2,n)=GammaCusp_q(0,n)*betaQCD(2,n)/4d0+GammaCusp_q(1,n)*betaQCD(1,n)/2d0&
                                +3d0/4d0*GammaCusp_q(2,n)*betaQCD(0,n)&
                                +3d0*dnk_q(2,0,n)*(betaQCD(0,n)**2)
        d_nk_Q_internal(4,1,n)=GammaCusp_q(3,n)/2d0+2d0*dnk_q(2,0,n)*betaQCD(1,n)&
                                +3d0*dnk_q(3,0,n)*betaQCD(0,n)
    end do
    
end subroutine SetDnkQuark

!!!! sets the values of logarithm part for RAD GLUON
!!!! Expression is taken from [1707.07606]
subroutine SetDnkGluon()
    integer::n
    
    !! 1-loop
    do n=NfMIN,NfMAX
        d_nk_G_internal(1,4,n)=0d0
        d_nk_G_internal(1,3,n)=0d0
        d_nk_G_internal(1,2,n)=0d0
        d_nk_G_internal(1,1,n)=GammaCusp_g(0,n)/2d0
    end do
    
    !! 2-loop
    do n=NfMIN,NfMAX
        d_nk_G_internal(2,4,n)=0d0
        d_nk_G_internal(2,3,n)=0d0
        d_nk_G_internal(2,2,n)=GammaCusp_g(0,n)*betaQCD(0,n)/4d0
        d_nk_G_internal(2,1,n)=GammaCusp_g(1,n)/2d0
    end do
    
    !! 3-loop
    do n=NfMIN,NfMAX
        d_nk_G_internal(3,4,n)=0d0
        d_nk_G_internal(3,3,n)=GammaCusp_g(0,n)*(betaQCD(0,n)**2)/6d0
        d_nk_G_internal(3,2,n)=(GammaCusp_g(0,n)*betaQCD(1,n)+2d0*GammaCusp_g(1,n)*betaQCD(0,n))/4d0
        d_nk_G_internal(3,1,n)=(GammaCusp_g(2,n)+4d0*dnk_g(2,0,n)*betaQCD(0,n))/2d0
    end do
    
    !! 4-loop
    do n=NfMIN,NfMAX
        d_nk_G_internal(4,4,n)=GammaCusp_g(0,n)*(betaQCD(0,n)**3)/8d0
        d_nk_G_internal(4,3,n)=5d0/12d0*GammaCusp_q(0,n)*betaQCD(0,n)*betaQCD(1,n)&
                                +GammaCusp_g(1,n)*(betaQCD(0,n)**2)/2d0
        d_nk_G_internal(4,2,n)=GammaCusp_g(0,n)*betaQCD(2,n)/4d0+GammaCusp_g(1,n)*betaQCD(1,n)/2d0&
                                +3d0/4d0*GammaCusp_g(2,n)*betaQCD(0,n)&
                                +3d0*dnk_g(2,0,n)*(betaQCD(0,n)**2)
        d_nk_G_internal(4,1,n)=GammaCusp_g(3,n)/2d0+2d0*dnk_g(2,0,n)*betaQCD(1,n)&
                                +3d0*dnk_g(3,0,n)*betaQCD(0,n)
    end do
    
end subroutine SetDnkGluon

!!!! sets the values of RAD-resummed coefficients for  QUARK
!!!! D=-GAMMA0/2BETA0(LOG(1-X)+as^n/(1-X)^n dnkl X^k Log(1-X)^l
subroutine SetDnklQuark()
    integer::n,i,j
    real(dp)::B1,B2,B3,B4,G1,G2,G3,G4,DD2,DD3,DD4
        
    do n=NfMIN,NfMAX
        B1=betaQCD(1,n)/betaQCD(0,n)
        B2=betaQCD(2,n)/betaQCD(0,n)
        B3=betaQCD(3,n)/betaQCD(0,n)
        B4=betaQCD(4,n)/betaQCD(0,n)
        
        G1=GammaCusp_q(1,n)/GammaCusp_q(0,n)
        G2=GammaCusp_q(2,n)/GammaCusp_q(0,n)
        G3=GammaCusp_q(3,n)/GammaCusp_q(0,n)
        G4=GammaCusp_q(4,n)/GammaCusp_q(0,n)
        
        DD2=-2d0*dnk_q(2,0,n)*betaQCD(0,n)/GammaCusp_q(0,n)
        DD3=-2d0*dnk_q(3,0,n)*betaQCD(0,n)/GammaCusp_q(0,n)
        DD4=-2d0*dnk_q(4,0,n)*betaQCD(0,n)/GammaCusp_q(0,n) 
        
        !! 1-loop
        !! *checked [18.05.22]
        d_nkl_Q_internal(1,0,1,n)=B1
        d_nkl_Q_internal(1,1,0,n)=B1-G1
        d_nkl_Q_internal(1,0,0,n)=0d0
    
        !! 2-loop 
        !! *checked [18.05.22]
        d_nkl_Q_internal(2,0,2,n)=-(B1**2)/2d0
        d_nkl_Q_internal(2,0,1,n)=B1*G1
        d_nkl_Q_internal(2,2,0,n)=(G2-B1*G1-B2+B1**2)/2d0
        d_nkl_Q_internal(2,1,0,n)=B1*G1-G2
        d_nkl_Q_internal(2,0,0,n)=DD2        
        d_nkl_Q_internal(2,2,2,n)=0d0
        d_nkl_Q_internal(2,1,1,n)=0d0
        
        !! 3-loop    
        d_nkl_Q_internal(3,0,3,n)=(B1**3)/3d0
        d_nkl_Q_internal(3,0,2,n)=-(B1**3)/2d0-(B1**2)*G1
        d_nkl_Q_internal(3,0,1,n)=B1*G2-2d0*B1*DD2
        d_nkl_Q_internal(3,1,1,n)=B1*B2-B1**3
        d_nkl_Q_internal(3,3,0,n)=1d0/3d0*(B1**3-2d0*B1*B2+B3-(B1**2)*G1+B2*G1+B1*G2-G3)
        d_nkl_Q_internal(3,2,0,n)=-0.5d0*(B1**3)+B1*B2-0.5d0*B3+(B1**2)*G1-B2*G1-B1*G2+G3
        d_nkl_Q_internal(3,1,0,n)=B1*G2-G3
        d_nkl_Q_internal(3,0,0,n)=DD3
        d_nkl_Q_internal(3,1,2,n)=0d0
        d_nkl_Q_internal(3,1,3,n)=0d0
        d_nkl_Q_internal(3,2,1,n)=0d0
        d_nkl_Q_internal(3,2,2,n)=0d0
        d_nkl_Q_internal(3,2,3,n)=0d0
        d_nkl_Q_internal(3,3,1,n)=0d0
        d_nkl_Q_internal(3,3,2,n)=0d0
        d_nkl_Q_internal(3,3,3,n)=0d0
        
        !! 4-loop    
        do i=0,4
        do j=0,4
            d_nkl_Q_internal(4,i,j,n)=0d0
        end do
        end do
        
        d_nkl_Q_internal(4,0,4,n)=-(B1**4)/4d0
        d_nkl_Q_internal(4,0,3,n)=5d0/6d0*(B1**4)+(B1**3)*G1
        d_nkl_Q_internal(4,0,2,n)=-(B1**2*B2)/2d0 + 3d0*B1**2*DD2 - B1**3*G1 - 3d0*B1**2*G2/2d0
        d_nkl_Q_internal(4,1,2,n)=B1**4 - B1**2*B2
        d_nkl_Q_internal(4,0,1,n)=-2d0*B1**2*DD2 - 3d0*B1*DD3 + B1*G3
        d_nkl_Q_internal(4,1,1,n)=-(B1**2*B2) + B1*B3 - 2d0*B1**3*G1 + 2d0*B1*B2*G1
        d_nkl_Q_internal(4,2,1,n)=-B1**4/2d0 + B1**2*B2 - (B1*B3)/2d0
        d_nkl_Q_internal(4,4,0,n)=(B1**4 - 3d0*B1**2*B2 + B2**2 + 2d0*B1*B3 - B4 &
                - B1**3*G1 + 2d0*B1*B2*G1 - B3*G1 + B1**2*G2 - B2*G2 - B1*G3 + G4)/4d0
        d_nkl_Q_internal(4,3,0,n)=(-2*B1**4 + 6*B1**2*B2 - 2*B2**2 - 4*B1*B3 + 2*B4 + 3*B1**3*G1 &
                - 6*B1*B2*G1 + 3*B3*G1 - 3*B1**2*G2 + 3*B2*G2 + 3*B1*G3 - 3*G4)/3d0
        d_nkl_Q_internal(4,2,0,n)=(-(B1**2*B2) + 2*B1*B3 - B4 - 2*B1**3*G1 + 4*B1*B2*G1 - 2*B3*G1 &
                + 3*B1**2*G2 - 3*B2*G2 - 3*B1*G3 + 3*G4)/2d0
        d_nkl_Q_internal(4,1,0,n)=-2*B1**2*DD2 + 2*B2*DD2 + B1*G3 - G4
        d_nkl_Q_internal(4,0,0,n)=DD4
    end do
    
end subroutine SetDnklQuark


!!!! sets the values of RAD-resummed coefficients for  QUARK
!!!! D=-GAMMA0/2BETA0(LOG(1-X)+as^n/(1-X)^n dnkl X^k Log(1-X)^l
subroutine SetDnklGluon()
    integer::n, i,j
    real(dp)::B1,B2,B3,B4,G1,G2,G3,G4,DD2,DD3,DD4
        
    do n=NfMIN,NfMAX
        B1=betaQCD(1,n)/betaQCD(0,n)
        B2=betaQCD(2,n)/betaQCD(0,n)
        B3=betaQCD(3,n)/betaQCD(0,n)
        B4=betaQCD(4,n)/betaQCD(0,n)
        
        G1=GammaCusp_g(1,n)/GammaCusp_g(0,n)
        G2=GammaCusp_g(2,n)/GammaCusp_g(0,n)
        G3=GammaCusp_g(3,n)/GammaCusp_g(0,n)
        G4=GammaCusp_g(4,n)/GammaCusp_g(0,n)
        
        DD2=-2d0*dnk_g(2,0,n)*betaQCD(0,n)/GammaCusp_g(0,n)
        DD3=-2d0*dnk_g(3,0,n)*betaQCD(0,n)/GammaCusp_g(0,n)
        DD4=-2d0*dnk_g(4,0,n)*betaQCD(0,n)/GammaCusp_g(0,n)
        
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
        
        !! 4-loop    
        do i=0,4
        do j=0,4
            d_nkl_G_internal(4,i,j,n)=0d0
        end do
        end do
        
        d_nkl_G_internal(4,0,4,n)=-(B1**4)/4d0
        d_nkl_G_internal(4,0,3,n)=5d0/6d0*(B1**4)+(B1**3)*G1
        d_nkl_G_internal(4,0,2,n)=-(B1**2*B2)/2d0 + 3d0*B1**2*DD2 - B1**3*G1 - 3d0*B1**2*G2/2d0
        d_nkl_G_internal(4,1,2,n)=B1**4 - B1**2*B2
        d_nkl_G_internal(4,0,1,n)=-2d0*B1**2*DD2 - 3d0*B1*DD3 + B1*G3
        d_nkl_G_internal(4,1,1,n)=-(B1**2*B2) + B1*B3 - 2d0*B1**3*G1 + 2d0*B1*B2*G1
        d_nkl_G_internal(4,2,1,n)=-B1**4/2d0 + B1**2*B2 - (B1*B3)/2d0
        d_nkl_G_internal(4,4,0,n)=(B1**4 - 3d0*B1**2*B2 + B2**2 + 2d0*B1*B3 - B4 &
                - B1**3*G1 + 2d0*B1*B2*G1 - B3*G1 + B1**2*G2 - B2*G2 - B1*G3 + G4)/4d0
        d_nkl_G_internal(4,3,0,n)=(-2*B1**4 + 6*B1**2*B2 - 2*B2**2 - 4*B1*B3 + 2*B4 + 3*B1**3*G1 &
                - 6*B1*B2*G1 + 3*B3*G1 - 3*B1**2*G2 + 3*B2*G2 + 3*B1*G3 - 3*G4)/3d0
        d_nkl_G_internal(4,2,0,n)=(-(B1**2*B2) + 2*B1*B3 - B4 - 2*B1**3*G1 + 4*B1*B2*G1 - 2*B3*G1 &
                + 3*B1**2*G2 - 3*B2*G2 - 3*B1*G3 + 3*G4)/2d0
        d_nkl_G_internal(4,1,0,n)=-2*B1**2*DD2 + 2*B2*DD2 + B1*G3 - G4
        d_nkl_G_internal(4,0,0,n)=DD4
        
    end do
    
end subroutine SetDnklGluon

!!!! sets the values of zeta-line PT coefficients for  QUARK
!!!! zeta=C0 mu/b*exp(-v), v=a^n L^k v^{(nk)}
subroutine SetVnkQuark()
    integer::n
    real(dp)::B1,B2,B3,G1,G2,G3,gg1,gg2,gg3,gg4,dd2,dd3,dd4
        
    v_nk_Q_internal=0d0*v_nk_Q_internal
        
    do n=NfMIN,NfMAX
        B1=betaQCD(1,n)/betaQCD(0,n)
        B2=betaQCD(2,n)/betaQCD(0,n)
        B3=betaQCD(3,n)/betaQCD(0,n)
        
        G1=GammaCusp_q(1,n)/GammaCusp_q(0,n)
        G2=GammaCusp_q(2,n)/GammaCusp_q(0,n)
        G3=GammaCusp_q(3,n)/GammaCusp_q(0,n)
        
        gg1=gammaV_q(1,n)/GammaCusp_q(0,n)
        gg2=(gammaV_q(2,n)+dnk_q(2,0,n))/GammaCusp_q(0,n)
        gg3=(gammaV_q(3,n)+dnk_q(3,0,n))/GammaCusp_q(0,n)
        gg4=(gammaV_q(4,n)+dnk_q(4,0,n))/GammaCusp_q(0,n)
        
        dd2=dnk_q(2,0,n)/GammaCusp_q(0,n)
        dd3=dnk_q(3,0,n)/GammaCusp_q(0,n)
        dd4=dnk_q(4,0,n)/GammaCusp_q(0,n)
        
        !! 1-loop        
        v_nk_Q_internal(0,1,n)=0d0
        v_nk_Q_internal(0,0,n)=gg1
    
        !! 2-loop    
        v_nk_Q_internal(1,2,n)=betaQCD(0,n)/12d0
        v_nk_Q_internal(1,1,n)=0d0
        v_nk_Q_internal(1,0,n)=-gg1*G1+gg2
        
        !! 3-loop    
        v_nk_Q_internal(2,3,n)=(betaQCD(0,n)**2)/24d0
        v_nk_Q_internal(2,2,n)=betaQCD(0,n)/12d0*(B1+G1)
        v_nk_Q_internal(2,1,n)=betaQCD(0,n)/2d0*(5d0/3d0*dd2-gg1*G1+gg2)
        v_nk_Q_internal(2,0,n)=gg1*(G1**2)-gg2*G1-gg1*G2+gg3
        
        !! 4-loop    
        v_nk_Q_internal(3,4,n)=19d0*(betaQCD(0,n)**3)/720d0
        v_nk_Q_internal(3,3,n)=(betaQCD(0,n)**2)/12d0*(G1+3d0/2d0*B1)
        v_nk_Q_internal(3,2,n)=betaQCD(0,n)/12d0*(B2 + B1*G1 - G1**2 + 2d0*G2 &
                + (14d0*dd2 - 5d0*G1*gg1 + 5d0*gg2)*betaQCD(0,n))
        v_nk_Q_internal(3,1,n)=betaQCD(0,n)/6d0*(5d0*B1*dd2 + 8d0*dd3 - 6d0*dd2*G1 &
                - 3d0*B1*G1*gg1 + 6d0*G1**2*gg1 - 6d0*G2*gg1 + 3d0*B1*gg2 &
                - 6d0*G1*gg2 + 6d0*gg3)
        v_nk_Q_internal(3,0,n)=-(G1**3*gg1) + 2d0*G1*G2*gg1 - G3*gg1 + G1**2*gg2 &
                - G2*gg2 - G1*gg3 + gg4 + ((-5d0*dd2**2)/3d0 - dd2*G1*gg1 + dd2*gg2)*betaQCD(0,n)     
                
        !! 5-loop partially
        v_nk_Q_internal(4,5,n)=3d0*(betaQCD(0,n)**4)/160d0
        v_nk_Q_internal(4,4,n)=(betaQCD(0,n)**3)/720d0*(104*B1 + 57*G1)
        v_nk_Q_internal(4,3,n)=betaQCD(0,n)**2/24d0*(2d0*B1**2 + 4d0*B2 + 5d0*B1*G1 - 3d0*G1**2 + 6d0*G2 &
                + ((511d0*dd2)/15d0 - 9d0*G1*gg1 + 9d0*gg2)*betaQCD(0,n))
        v_nk_Q_internal(4,2,n)=betaQCD(0,n)/12d0*(B3 + B2*G1 - B1*G1**2 + G1**3 + 2d0*B1*G2 - 3d0*G1*G2 &
                + 3*G3 + (35d0*B1*dd2 + 33d0*dd3 - 24d0*dd2*G1 - 12d0*B1*G1*gg1 + 13d0*G1**2*gg1 &
                - 14d0*G2*gg1 + 12d0*B1*gg2 - 13d0*G1*gg2 + 14d0*gg3)*betaQCD(0,n))
        v_nk_Q_internal(4,1,n)=betaQCD(0,n)/6d0*(5*B2*dd2 + 8*B1*dd3 + 11*dd4 - 6*B1*dd2*G1 - 9*dd3*G1 &
                + 7*dd2*G1**2 - 7*dd2*G2 - 3*B2*G1*gg1 + 6*B1*G1**2*gg1 - 9*G1**3*gg1 - 6*B1*G2*gg1 &
                +  18*G1*G2*gg1 - 9*G3*gg1 + 3*B2*gg2 - 6*B1*G1*gg2 + 9*G1**2*gg2 - 9*G2*gg2 + 6*B1*gg3 &
                - 9*G1*gg3 + 9*gg4 + (-29*dd2**2 - 16*dd2*G1*gg1 + 16*dd2*gg2)*betaQCD(0,n))
        v_nk_Q_internal(4,0,n)=0d0
        
    end do
    
end subroutine SetVnkQuark

!!!! sets the values of zeta-line PT coefficients for  QUARK
!!!! zeta=C0 mu/b*exp(-v), v=a^n L^k v^{(nk)}
subroutine SetVnkGluon()
    integer::n
    real(dp)::B1,B2,B3,G1,G2,G3,gg1,gg2,gg3,gg4,dd2,dd3,dd4
    
    v_nk_G_internal=0d0*v_nk_G_internal
        
    do n=NfMIN,NfMAX
        B1=betaQCD(1,n)/betaQCD(0,n)
        B2=betaQCD(2,n)/betaQCD(0,n)
        B2=betaQCD(3,n)/betaQCD(0,n)
        
        G1=GammaCusp_g(1,n)/GammaCusp_g(0,n)
        G2=GammaCusp_g(2,n)/GammaCusp_g(0,n)
        G3=GammaCusp_g(3,n)/GammaCusp_g(0,n)
        
        gg1=gammaV_g(1,n)/GammaCusp_g(0,n)
        gg2=(gammaV_g(2,n)+dnk_g(2,0,n))/GammaCusp_g(0,n)
        gg3=(gammaV_g(3,n)+dnk_g(3,0,n))/GammaCusp_g(0,n)
        gg4=(gammaV_g(4,n)+dnk_g(4,0,n))/GammaCusp_g(0,n)
        
        dd2=dnk_g(2,0,n)/GammaCusp_g(0,n)
        dd3=dnk_g(3,0,n)/GammaCusp_g(0,n)
        dd4=dnk_g(4,0,n)/GammaCusp_g(0,n)
        
        !! 1-loop        
        v_nk_G_internal(0,1,n)=0d0
        v_nk_G_internal(0,0,n)=gg1
    
        !! 2-loop    
        v_nk_G_internal(1,2,n)=betaQCD(0,n)/12d0
        v_nk_G_internal(1,1,n)=0d0
        v_nk_G_internal(1,0,n)=-gg1*G1+gg2
        
        !! 3-loop    
        v_nk_G_internal(2,3,n)=(betaQCD(0,n)**2)/24d0
        v_nk_G_internal(2,2,n)=betaQCD(0,n)/12d0*(B1+G1)
        v_nk_G_internal(2,1,n)=betaQCD(0,n)/2d0*(5d0/3d0*dd2-gg1*G1+gg2)
        v_nk_G_internal(2,0,n)=gg1*(G1**2)-gg2*G1-gg1*G2+gg3
        
        !! 4-loop    
        v_nk_G_internal(3,4,n)=19d0*(betaQCD(0,n)**3)/720d0
        v_nk_G_internal(3,3,n)=(betaQCD(0,n)**2)/12d0*(G1+3d0/2d0*B1)
        v_nk_G_internal(3,2,n)=betaQCD(0,n)/12d0*(B2 + B1*G1 - G1**2 + 2d0*G2 &
                + (14d0*dd2 - 5d0*G1*gg1 + 5d0*gg2)*betaQCD(0,n))
        v_nk_G_internal(3,1,n)=betaQCD(0,n)/6d0*(5d0*B1*dd2 + 8d0*dd3 - 6d0*dd2*G1 &
                - 3d0*B1*G1*gg1 + 6d0*G1**2*gg1 - 6d0*G2*gg1 + 3d0*B1*gg2 &
                - 6d0*G1*gg2 + 6d0*gg3)
        v_nk_G_internal(3,0,n)=-(G1**3*gg1) + 2d0*G1*G2*gg1 - G3*gg1 + G1**2*gg2 &
                - G2*gg2 - G1*gg3 + gg4 + ((-5d0*dd2**2)/3d0 - dd2*G1*gg1 + dd2*gg2)*betaQCD(0,n)
                
        !! 5-loop partially
        v_nk_G_internal(4,5,n)=3d0*(betaQCD(0,n)**4)/160d0
        v_nk_G_internal(4,4,n)=(betaQCD(0,n)**3)/720d0*(104*B1 + 57*G1)
        v_nk_G_internal(4,3,n)=betaQCD(0,n)**2/24d0*(2d0*B1**2 + 4d0*B2 + 5d0*B1*G1 - 3d0*G1**2 + 6d0*G2 &
                + ((511d0*dd2)/15d0 - 9d0*G1*gg1 + 9d0*gg2)*betaQCD(0,n))
        v_nk_G_internal(4,2,n)=betaQCD(0,n)/12d0*(B3 + B2*G1 - B1*G1**2 + G1**3 + 2d0*B1*G2 - 3d0*G1*G2 &
                + 3*G3 + (35d0*B1*dd2 + 33d0*dd3 - 24d0*dd2*G1 - 12d0*B1*G1*gg1 + 13d0*G1**2*gg1 &
                - 14d0*G2*gg1 + 12d0*B1*gg2 - 13d0*G1*gg2 + 14d0*gg3)*betaQCD(0,n))
        v_nk_G_internal(4,1,n)=betaQCD(0,n)/6d0*(5*B2*dd2 + 8*B1*dd3 + 11*dd4 - 6*B1*dd2*G1 - 9*dd3*G1 &
                + 7*dd2*G1**2 - 7*dd2*G2 - 3*B2*G1*gg1 + 6*B1*G1**2*gg1 - 9*G1**3*gg1 - 6*B1*G2*gg1 &
                +  18*G1*G2*gg1 - 9*G3*gg1 + 3*B2*gg2 - 6*B1*G1*gg2 + 9*G1**2*gg2 - 9*G2*gg2 + 6*B1*gg3 &
                - 9*G1*gg3 + 9*gg4 + (-29*dd2**2 - 16*dd2*G1*gg1 + 16*dd2*gg2)*betaQCD(0,n))
        v_nk_G_internal(4,0,n)=0d0
        
    end do
    
end subroutine SetVnkGluon

!!!! sets the values of EXACT zeta-line PT coefficients for  QUARK
!!!! zeta=mu^2*exp(-1/beta0/as * OMEGA), OMEGA=a^n OMEGA^{(nk)} ...
!!!! ... in each term different.
subroutine SetOMEGAnkQuark()
    integer::n
    real(dp)::B1,B2,B3,B4,G1,G2,G3,G4,gg1,gg2,gg3,gg4,commonF
    
    OMEGA_nk_Q_internal=0d0*OMEGA_nk_Q_internal
        
    do n=NfMIN,NfMAX
        B1=betaQCD(1,n)/betaQCD(0,n)
        B2=betaQCD(2,n)/betaQCD(0,n)
        B3=betaQCD(3,n)/betaQCD(0,n)
        B4=betaQCD(4,n)/betaQCD(0,n)
        
        G1=GammaCusp_q(1,n)/GammaCusp_q(0,n)
        G2=GammaCusp_q(2,n)/GammaCusp_q(0,n)
        G3=GammaCusp_q(3,n)/GammaCusp_q(0,n)
        G4=GammaCusp_q(4,n)/GammaCusp_q(0,n)
        
        gg1=gammaV_q(1,n)*betaQCD(0,n)/GammaCusp_q(0,n)
        gg2=gammaV_q(2,n)*betaQCD(0,n)/GammaCusp_q(0,n)
        gg3=gammaV_q(3,n)*betaQCD(0,n)/GammaCusp_q(0,n)
        gg4=gammaV_q(4,n)*betaQCD(0,n)/GammaCusp_q(0,n)
        
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
        !!! corrected 17.11.2022
        OMEGA_nk_Q_internal(3,2,n)=(-B2+B1*G1+G1**2-G2-2d0*G1*gg1+2d0*gg2)*(G1-B1)/2d0*commonF
        !! * z_{-2}
        OMEGA_nk_Q_internal(3,3,n)=((-B1*B2-B3+(B1**2)*G1+2d0*B2*G1-4d0*(G1**3)-B1*G2+6d0*G1*G2-2d0*G3)/12d0&
                +(2d0*G1-B1)*(G1*gg1-gg2)/2d0-(G2*gg1-gg3)/2d0)*commonF
        !! * 1
        OMEGA_nk_Q_internal(3,4,n)=((G1**2)*gg1-G2*gg1-G1*gg2+gg3)*commonF
        !! * p (leading term of expansion z[1]+z[-1]+z[-2] parts.)
        !! this is used to cure the groving tail
        OMEGA_nk_Q_internal(3,5,n)=(-G1**3 + 2*G1*G2 - G3 - B1*G1*gg1 + 3*G1**2*gg1 &
                        - 2*G2*gg1 + B1*gg2 - 3*G1*gg2 + 2*gg3)/2d0*commonF
        
        !! 4-loop
        !! * z_1
        OMEGA_nk_Q_internal(4,1,n)=(2d0*B1**2*B2 - 3d0*B2**2 - 4d0*B1*B3 + 6d0*B4 - 2d0*B1**3*G1 &
                    + 6d0*B1*B2*G1 - 2d0*B3*G1 - B1**2*G1**2 - 2d0*B2*G1**2 + 2d0*B1*G1**3 &
                    + G1**4 + 2d0*B1**2*G2 - 2d0*B1*G1*G2 - 6d0*G1**2*G2 + 3d0*G2**2 - 2d0*B1*G3 &
                    + 8d0*G1*G3 - 6d0*G4)/24d0*commonF
        !! * z_{-1}
        OMEGA_nk_Q_internal(4,2,n)=(2d0*B1**2 - B2 - 3d0*B1*G1 + G1**2 + G2)*(-B2 + B1*G1 &
                + G1**2 - G2 - 2d0*G1*gg1 + 2d0*gg2)/4d0*commonF
        !! * z_{-2}
        OMEGA_nk_Q_internal(4,3,n)=(B1*B2 + B3 - B1**2*G1 - 2d0*B2*G1 + 4d0*G1**3 + B1*G2 - 6d0*G1*G2 &
                + 2d0*G3 + 6d0*B1*G1*gg1 - 12d0*G1**2*gg1 + 6d0*G2*gg1 - 6d0*B1*gg2 &
                + 12d0*G1*gg2 - 6d0*gg3)*(B1-G1)/6d0*commonF
        !! * z_{-3}
        OMEGA_nk_Q_internal(4,4,n)=(-2*B1**2*B2 - B2**2 - 4*B1*B3 - 2*B4 + 2*B1**3*G1 + 10*B1*B2*G1 &
                + 6*B3*G1 - 3*B1**2*G1**2 - 6*B2*G1**2 - 18*B1*G1**3 + 27*G1**4 - 2*B1**2*G2 + 30*B1*G1*G2 &
                - 54*G1**2*G2 + 9*G2**2 - 10*B1*G3 + 24*G1*G3 - 6*G4 - 24*B1**2*G1*gg1 - 12*B2*G1*gg1 &
                + 108*B1*G1**2*gg1 - 108*G1**3*gg1 - 48*B1*G2*gg1 + 108*G1*G2*gg1 - 24*G3*gg1 &
                + 24*B1**2*gg2 + 12*B2*gg2 - 108*B1*G1*gg2 + 108*G1**2*gg2 - 36*G2*gg2 + 48*B1*gg3 &
                - 72*G1*gg3 + 24*gg4)/72d0*commonF
        !! * 1
        OMEGA_nk_Q_internal(4,5,n)=(-(G1**3*gg1) + 2*G1*G2*gg1 - G3*gg1 + G1**2*gg2 &
                - G2*gg2 - G1*gg3 + gg4)*commonF
    end do
    
!     write(*,*) "w0"
!     write(*,*) OMEGA_nk_Q_internal(0,1,3),OMEGA_nk_Q_internal(0,1,4),OMEGA_nk_Q_internal(0,1,5)
!     write(*,*) "w1"
!     write(*,*) OMEGA_nk_Q_internal(1,1,3),OMEGA_nk_Q_internal(1,1,4),OMEGA_nk_Q_internal(1,1,5)
!     write(*,*) OMEGA_nk_Q_internal(1,2,3),OMEGA_nk_Q_internal(1,2,4),OMEGA_nk_Q_internal(1,2,5)
!     write(*,*) OMEGA_nk_Q_internal(1,3,3),OMEGA_nk_Q_internal(1,3,4),OMEGA_nk_Q_internal(1,3,5)
!     write(*,*) "w2"
!     write(*,*) OMEGA_nk_Q_internal(2,1,3),OMEGA_nk_Q_internal(2,1,4),OMEGA_nk_Q_internal(2,1,5)
!     write(*,*) OMEGA_nk_Q_internal(2,2,3),OMEGA_nk_Q_internal(2,2,4),OMEGA_nk_Q_internal(2,2,5)
!     write(*,*) OMEGA_nk_Q_internal(2,3,3),OMEGA_nk_Q_internal(2,3,4),OMEGA_nk_Q_internal(2,3,5)
!     write(*,*) "w3"
!     write(*,*) OMEGA_nk_Q_internal(3,1,3),OMEGA_nk_Q_internal(3,1,4),OMEGA_nk_Q_internal(3,1,5)
!     write(*,*) OMEGA_nk_Q_internal(3,2,3),OMEGA_nk_Q_internal(3,2,4),OMEGA_nk_Q_internal(3,2,5)
!     write(*,*) OMEGA_nk_Q_internal(3,3,3),OMEGA_nk_Q_internal(3,3,4),OMEGA_nk_Q_internal(3,3,5)
!     write(*,*) OMEGA_nk_Q_internal(3,4,3),OMEGA_nk_Q_internal(3,4,4),OMEGA_nk_Q_internal(3,4,5)
!     write(*,*) OMEGA_nk_Q_internal(3,5,3),OMEGA_nk_Q_internal(3,5,4),OMEGA_nk_Q_internal(3,5,5)
!     write(*,*) "w4"
!     write(*,*) OMEGA_nk_Q_internal(4,1,3),OMEGA_nk_Q_internal(4,1,4),OMEGA_nk_Q_internal(4,1,5)
!     write(*,*) OMEGA_nk_Q_internal(4,2,3),OMEGA_nk_Q_internal(4,2,4),OMEGA_nk_Q_internal(4,2,5)
!     write(*,*) OMEGA_nk_Q_internal(4,3,3),OMEGA_nk_Q_internal(4,3,4),OMEGA_nk_Q_internal(4,3,5)
!     write(*,*) OMEGA_nk_Q_internal(4,4,3),OMEGA_nk_Q_internal(4,4,4),OMEGA_nk_Q_internal(4,4,5)
!     write(*,*) OMEGA_nk_Q_internal(4,5,3),OMEGA_nk_Q_internal(4,5,4),OMEGA_nk_Q_internal(4,5,5)
    
end subroutine SetOMEGAnkQuark

!!!! sets the values of EXACT zeta-line PT coefficients for  QUARK
!!!! zeta=mu^2*exp(-1/as * OMEGA), OMEGA=a^n OMEGA^{(nk)} ...
!!!! ... in each term different.
subroutine SetOMEGAnkGluon()
    integer::n
    real(dp)::B1,B2,B3,B4,G1,G2,G3,G4,gg1,gg2,gg3,gg4,commonF
    
    OMEGA_nk_G_internal=0d0*OMEGA_nk_G_internal
        
    do n=NfMIN,NfMAX
        B1=betaQCD(1,n)/betaQCD(0,n)
        B2=betaQCD(2,n)/betaQCD(0,n)
        B3=betaQCD(3,n)/betaQCD(0,n)
        B4=betaQCD(4,n)/betaQCD(0,n)
        
        G1=GammaCusp_g(1,n)/GammaCusp_g(0,n)
        G2=GammaCusp_g(2,n)/GammaCusp_g(0,n)
        G3=GammaCusp_g(3,n)/GammaCusp_g(0,n)
        G4=GammaCusp_g(4,n)/GammaCusp_g(0,n)
        
        gg1=gammaV_g(1,n)*betaQCD(0,n)/GammaCusp_g(0,n)
        gg2=gammaV_g(2,n)*betaQCD(0,n)/GammaCusp_g(0,n)
        gg3=gammaV_g(3,n)*betaQCD(0,n)/GammaCusp_g(0,n)
        gg4=gammaV_g(4,n)*betaQCD(0,n)/GammaCusp_g(0,n)
        
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
        !! * p (leading term of expansion z[-1]+z[-2] parts.)
        !! this is used to cure the groving tail
        OMEGA_nk_G_internal(3,5,n)=(B1*B2 - 2*B3 - B1**2*G1 + B2*G1 - 5*G1**3 + B1*G2 + 9*G1*G2 - 4*G3 &
                - 6*B1*G1*gg1 + 18*G1**2*gg1 - 12*G2*gg1 + 6*B1*gg2 - 18*G1*gg2 + 12*gg3)/12d0*commonF
        
        !! 4-loop
        !! * z_1
        OMEGA_nk_G_internal(4,1,n)=(2d0*B1**2*B2 - 3d0*B2**2 - 4d0*B1*B3 + 6d0*B4 - 2d0*B1**3*G1 &
                    + 6d0*B1*B2*G1 - 2d0*B3*G1 - B1**2*G1**2 - 2d0*B2*G1**2 + 2d0*B1*G1**3 &
                    + G1**4 + 2d0*B1**2*G2 - 2d0*B1*G1*G2 - 6d0*G1**2*G2 + 3d0*G2**2 - 2d0*B1*G3 &
                    + 8d0*G1*G3 - 6d0*G4)/24d0*commonF
        !! * z_{-1}
        OMEGA_nk_G_internal(4,2,n)=(2d0*B1**2 - B2 - 3d0*B1*G1 + G1**2 + G2)*(-B2 + B1*G1 &
                + G1**2 - G2 - 2d0*G1*gg1 + 2d0*gg2)/4d0*commonF
        !! * z_{-2}
        OMEGA_nk_G_internal(4,3,n)=(B1*B2 + B3 - B1**2*G1 - 2d0*B2*G1 + 4d0*G1**3 + B1*G2 - 6d0*G1*G2 &
                + 2d0*G3 + 6d0*B1*G1*gg1 - 12d0*G1**2*gg1 + 6d0*G2*gg1 - 6d0*B1*gg2 &
                + 12d0*G1*gg2 - 6d0*gg3)*(B1-G1)/6d0*commonF
        !! * z_{-3}
        OMEGA_nk_G_internal(4,4,n)=(-2*B1**2*B2 - B2**2 - 4*B1*B3 - 2*B4 + 2*B1**3*G1 + 10*B1*B2*G1 &
                + 6*B3*G1 - 3*B1**2*G1**2 - 6*B2*G1**2 - 18*B1*G1**3 + 27*G1**4 - 2*B1**2*G2 + 30*B1*G1*G2 &
                - 54*G1**2*G2 + 9*G2**2 - 10*B1*G3 + 24*G1*G3 - 6*G4 - 24*B1**2*G1*gg1 - 12*B2*G1*gg1 &
                + 108*B1*G1**2*gg1 - 108*G1**3*gg1 - 48*B1*G2*gg1 + 108*G1*G2*gg1 - 24*G3*gg1 &
                + 24*B1**2*gg2 + 12*B2*gg2 - 108*B1*G1*gg2 + 108*G1**2*gg2 - 36*G2*gg2 + 48*B1*gg3 &
                - 72*G1*gg3 + 24*gg4)/72d0*commonF
        !! * 1
        OMEGA_nk_G_internal(4,5,n)=(-(G1**3*gg1) + 2*G1*G2*gg1 - G3*gg1 + G1**2*gg2 &
                - G2*gg2 - G1*gg3 + gg4)*commonF
    end do
    
end subroutine SetOMEGAnkGluon

subroutine SetBetaRoots()
    integer::n,i,j
    real(dp)::B1,B2,B3,B4
    COMPLEX*16 :: a,answ
    
    !!!!! TABLES OF ROOTS ARE COMPUTED BY MATHEMATICA
    !!!!! for pairs of complex roots the root with Im>0 is first
    !!! 3-loop
    COMPLEX*16,dimension(0:6),parameter:: root21=[&
    (-0.03570178508925446272313615680784,0.080160885478465432575800740549553),&
    (-0.038675539164595526336887677383148,0.086322238704594954924135437271609),&
    (-0.042909558259571733587612197093758,0.094759912532946789061512692156804),&
    (-0.049702303908879109500388299249288,0.107277462500385745931517856995518),&
    (-0.06316365127831199015631408649683,0.12852249664419349413270267645244),&
    (-0.106868666188965093663629849524,0.17594897561664481697781266179055),&
    (-0.21268639236122698844346103757164,0.)]
    COMPLEX*16,dimension(0:6),parameter:: root22=[&
    (-0.03570178508925446272313615680784,-0.080160885478465432575800740549553),&
    (-0.038675539164595526336887677383148,-0.086322238704594954924135437271609),&
    (-0.042909558259571733587612197093758,-0.094759912532946789061512692156804),&
    (-0.049702303908879109500388299249288,-0.107277462500385745931517856995518),&
    (-0.06316365127831199015631408649683,-0.12852249664419349413270267645244),&
    (-0.106868666188965093663629849524,-0.17594897561664481697781266179055),&
    (1.0126863923612269884434610375716,0.)]
    
    !!! 4-loop
    COMPLEX*16,dimension(0:6),parameter:: root31=[&
    -0.07241478299428282,-0.07680070354530623,-0.08208483490867254,&
    -0.08849217512975707,-0.09621061277365012,-0.10508908839765044,-0.11362023213019258]
    COMPLEX*16,dimension(0:6),parameter:: root32=[&
    (0.011782712908695248,0.07110324448715928),&
    (0.012965518598492015,0.07588307304102093),&
    (0.014740652117464371,0.08195773114463775),&
    (0.01762023160041607,0.09000832305675134),&
    (0.02281953013830925,0.10128572857488127),&
    (0.0338021542840822,0.11821077095791804),&
    (0.0633829862094608,0.1445765143534799)]
    COMPLEX*16,dimension(0:6),parameter:: root33=[&
    (0.011782712908695248,-0.07110324448715928),&
    (0.012965518598492015,-0.07588307304102093),&
    (0.014740652117464371,-0.08195773114463775),&
    (0.01762023160041607,-0.09000832305675134),&
    (0.02281953013830925,-0.10128572857488127),&
    (0.0338021542840822,-0.11821077095791804),&
    (0.0633829862094608,-0.1445765143534799)]
    
    !!! 5-loop
    COMPLEX*16,dimension(0:6),parameter:: root41=[&
    (-0.05591091756225722,0.04182926526696407),&
    (-0.06136700935148767,0.044883180457591494),&
    (-0.06931025598997154,0.0486907342529324),&
    (-0.08235582183306786,0.05320124554862426),&
    (-0.10910510735854696,0.05533447765266648),&
    (-0.12511060289367348,0.),&
    (-0.11397406563300194,0.)]
    COMPLEX*16,dimension(0:6),parameter:: root42=[&
    (-0.05591091756225722,-0.04182926526696407),&
    (-0.06136700935148767,-0.044883180457591494),&
    (-0.06931025598997154,-0.0486907342529324),&
    (-0.08235582183306786,-0.05320124554862426),&
    (-0.10910510735854696,-0.05533447765266648),&
    (-0.2793998493629861,0),&
    (-9.12266237853665,0)]
     COMPLEX*16,dimension(0:6),parameter:: root43=[&
    (0.028690317392283806,0.05811172924003517),&
    (0.03054701511242339,0.06261303090284312),&
    (0.03290023100445615,0.06851832657672867),&
    (0.03598915425211804,0.07671897555874341),&
    (0.040205175247139464,0.08906642425176058),&
    (0.0462770317488488,0.10970531146206083),&
    (0.06410815653997126,0.1438532871828738)]
    COMPLEX*16,dimension(0:6),parameter:: root44=[&
    (0.028690317392283806,-0.05811172924003517),&
    (0.03054701511242339,-0.06261303090284312),&
    (0.03290023100445615,-0.06851832657672867),&
    (0.03598915425211804,-0.07671897555874341),&
    (0.040205175247139464,-0.08906642425176058),&
    (0.0462770317488488,-0.10970531146206083),&
    (0.06410815653997126,-0.1438532871828738)]
    
    betaRoots_internal=0d0*betaRoots_internal
    do n=NfMIN,NfMAX
        B1=betaQCD(1,n)/betaQCD(0,n)        
        !!! 2-loop
        !!! 1 real root
        betaRoots_internal(1,1,n)=-1d0/B1        
    end do
    !!! rest roots are tabulated
    if(NfMIN>=0 .and. NfMAX<=6) then
        !!! 3-loop
        betaRoots_internal(2,1,NfMIN:NfMAX)=root21(NfMIN:NfMAX)
        betaRoots_internal(2,2,NfMIN:NfMAX)=root22(NfMIN:NfMAX)
        !!! 4-loop
        betaRoots_internal(3,1,NfMIN:NfMAX)=root31(NfMIN:NfMAX)
        betaRoots_internal(3,2,NfMIN:NfMAX)=root32(NfMIN:NfMAX)
        betaRoots_internal(3,3,NfMIN:NfMAX)=root33(NfMIN:NfMAX)        
        !!! 5-loop
        betaRoots_internal(4,1,NfMIN:NfMAX)=root41(NfMIN:NfMAX)
        betaRoots_internal(4,2,NfMIN:NfMAX)=root42(NfMIN:NfMAX)
        betaRoots_internal(4,3,NfMIN:NfMAX)=root43(NfMIN:NfMAX)
        betaRoots_internal(4,4,NfMIN:NfMAX)=root44(NfMIN:NfMAX)
    else
        !!! 3-loop
        betaRoots_internal(2,1,0:6)=root21(0:6)
        betaRoots_internal(2,2,0:6)=root22(0:6)
        !!! 4-loop
        betaRoots_internal(3,1,0:6)=root31(0:6)
        betaRoots_internal(3,2,0:6)=root32(0:6)
        betaRoots_internal(3,3,0:6)=root33(0:6)   
        !!! 5-loop
        betaRoots_internal(4,1,0:6)=root41(0:6)
        betaRoots_internal(4,2,0:6)=root42(0:6)
        betaRoots_internal(4,3,0:6)=root43(0:6)
        betaRoots_internal(4,4,0:6)=root44(0:6)
        write(*,*) ErrorString('Beta-roots are not defined for Nf<0, and Nf>6.',moduleName)
        write(*,*) 'All cases beyond:  [Nf(min),Nf(max)] = [', NfMIN, ', ', NfMAX, '],  are ignored.'
    end if
    
    
    !!! checking procedure
    do n=NfMIN,NfMAX
        B1=betaQCD(1,n)/betaQCD(0,n)
        B2=betaQCD(2,n)/betaQCD(0,n)
        B3=betaQCD(3,n)/betaQCD(0,n)
        B4=betaQCD(4,n)/betaQCD(0,n)
        
        
        do i=1,4
        do j=1,i
            
            a=betaRoots_internal(i,j,n)
            answ=1d0+B1*a
            if(i>1) answ=answ+B2*a**2
            if(i>2) answ=answ+B3*a**3
            if(i>3) answ=answ+B4*a**4
            
            !!! the test of roots must be relative. I.e. for large roots, absolute precision decreses
            if(abs(dreal(answ)/max(1.,abs(a**i)))<0.000001 .and. abs(aimag(answ)/max(1.,abs(a**i)))<0.000001 )then
                !write(*,'(I1,"-loop , root-", I1, " Nf=",I2, ": PASS")') i,j,n
            else
                write(*,*) ErrorString('The beta-root is incorrect:',moduleName)
                write(*,'(I1,"-loop , root-", I1, " Nf=",I2, ": root incorrect")') i,j,n
                write(*,*) "EQN yealds=", answ
            end if
        end do
        end do
    end do

end subroutine SetBetaRoots
