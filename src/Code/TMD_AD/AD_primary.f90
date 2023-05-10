!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!
!!!!!!!!!       part of TMD_AD module for artemide
!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!-------------------------------------------------------------------------------------------------------
!!!--------------------------------------------PRIMARY ADs------------------------------------------------
!!!-------------------------------------------------------------------------------------------------------
!!!! sets the values of beta function
!!!! the expressions are taken from [1701.01404]
subroutine SetBetaFunction()
    integer::n
    
    !! 1-loop
    do n=NfMIN,NfMAX
        beta_internal(0,n)=11d0/3d0*CA-4d0/3d0*TF*n
    end do
    
    !! 2-loop
    do n=NfMIN,NfMAX
        beta_internal(1,n)=34d0/3d0*(CA**2)-20d0/3d0*CA*TF*n-4d0*CF*TF*n
    end do
    
    !! 3-loop
    do n=NfMIN,NfMAX
        beta_internal(2,n)=2857d0/54d0*(CA**3) -1415d0/27d0*(CA**2)*TF*n -205d0/9d0*CF*CA*TF*n +2d0*(CF**2)*TF*n&
                            +44d0/9d0*CF*(TF**2)*(n**2) +158d0/27d0*CA*(TF**2)*(n**2)
    end do
    
    !! 4-loop
    do n=NfMIN,NfMAX
        beta_internal(3,n)=(CA**4)*(150653d0/486d0 -44d0/9d0*zeta3) +dAA*(-80d0/9d0 +704d0/3d0*zeta3)&
                            +(CA**3)*TF*n*(-39143d0/81d0 +136d0/3d0*zeta3) +(CA**2)*CF*TF*n*(7073d0/243d0 -656d0/9d0*zeta3)&
                            +CA*(CF**2)*TF*n*(-4204d0/27d0 +352d0/9d0*zeta3) +dFA*n*(512d0/9d0 -1664d0/3d0*zeta3)&
                            +46d0*(CF**3)*TF*n +(CA**2)*(TF**2)*(n**2)*(7930d0/81d0 +224d0/9d0*zeta3)&
                            +(CF**2)*(TF**2)*(n**2)*(1352d0/27d0 -704d0/9d0*zeta3)&
                            +CA*CF*(TF**2)*(n**2)*(17152d0/243d0 +448d0/9d0*zeta3) +dFF*(n**2)*(-704d0/9d0 +512d0/3d0*zeta3)&
                            +424d0/243d0*CA*(TF**3)*(n**3) +1232d0/243d0*CF*(TF**3)*(n**3)
    end do
    
    !! 5-loop
    do n=NfMIN,NfMAX
        beta_internal(4,n)=(CA**5)*(8296235d0/3888d0 -1630d0/81d0*zeta3 +121d0/6d0*zeta4 -1045d0/9d0*zeta5)&
                            +dAA*CA*(-514d0/3d0 +18716d0/3d0*zeta3 -968d0*zeta4 -15400d0/3d0*zeta5)&
                            +(CA**4)*TF*n*(-5048959d0/972d0 +10505d0/81d0*zeta3 -583d0/3d0*zeta4 +1230d0*zeta5)&
                            +(CA**3)*CF*TF*n*(8141995d0/1944d0 +146d0*zeta3 +902d0/3d0*zeta4 -8720d0/3d0*zeta5)&
                            +(CA**2)*(CF**2)*TF*n*(-548732d0/81d0 -50581d0/27d0*zeta3 -484d0/3d0*zeta4 +12820d0/3d0*zeta5)&
                            +CA*(CF**3)*TF*n*(3717d0 +5696d0/3d0*zeta3 -7480d0/3d0*zeta5)&
                            -(CF**4)*TF*n*(4157d0/6d0 +128d0*zeta3)&
                            +dAA*TF*n*(904d0/9d0 -20752d0/9d0*zeta3 +352d0*zeta4 +4000d0/9d0*zeta5)&
                            +dFA*CA*n*(11312d0/9d0 -127736d0/9d0*zeta3 +2288d0*zeta4 +67520d0/9d0*zeta5)&
                            +dFA*CF*n*(-320d0 +1280d0/3d0*zeta3 +6400d0/3d0*zeta5)&
                            +(CA**3)*(TF**2)*(n**2)*(843067d0/486d0 +18446d0/27d0*zeta3 -104d0/3d0*zeta4-2200d0/3d0*zeta5)&
                            +(CA**2)*CF*(TF**2)*(n**2)*(5701d0/162d0 +26452d0/27d0*zeta3 -944d0/3d0*zeta4 +1600d0/3d0*zeta5)&
                            +(CF**2)*CA*(TF**2)*(n**2)*(31583d0/18d0 -28628d0/27d0*zeta3 +1144d0/3d0*zeta4 -4400d0/3d0*zeta5)&
                            +(CF**3)*(TF**2)*(n**2)*(-5018d0/9d0 -2144d0/3d0*zeta3 +4640d0/3d0*zeta5)&
                            +dFA*TF*(n**2)*(-3680d0/9d0 +40160d0/9d0*zeta3 -832d0*zeta4 -1280d0/9d0*zeta5)&
                            +dFF*CA*(n**2)*(-7184d0/3d0 +40336d0/9d0*zeta3 -704d0*zeta4 +2240d0/9d0*zeta5)&
                            +dFF*CF*(n**2)*(4160d0/3d0 +5120d0/3d0*zeta3 -12800d0/3d0*zeta5)&
                            +(CA**2)*(TF**3)*(n**3)*(-2077d0/27d0 -9736d0/81d0*zeta3 +112d0/3d0*zeta4 +320d0/9d0*zeta5)&
                            +CA*CF*(TF**3)*(n**3)*(-736d0/81d0 -5680d0/27d0*zeta3 +224d0/3d0*zeta4)&
                            +(CF**2)*(TF**3)*(n**3)*(-9922d0/81d0 +7616d0/27d0*zeta3 -352d0/3d0*zeta4)&
                            +dFF*TF*(n**3)*(3520d0/9d0 - 2624d0/3d0*zeta3 +256d0*zeta4 +1280d0/3d0*zeta5)&
                            +CA*(TF**4)*(n**4)*(916d0/243d0 -640d0/81d0*zeta3) - CF*(TF**4)*(n**4)*(856d0/243d0 +128d0/27d0*zeta3)
    end do


end subroutine SetBetaFunction

!!!! sets the values of cusp anomalous dimension
!!!! the 4-loop expression is taken from [2001.11377] (appendix C)
!!!! the 5-loop expression is taken from [1812.11818] formula (13), central value(!!) 
!!!!      only nf=3,4,5 avalible
subroutine SetGammaCuspQuark()
    integer::n
    
    !! 1-loop
    do n=NfMIN,NfMAX
        GammaCuspQ_internal(0,n)=CF*4d0
    end do
    
    !! 2-loop
    do n=NfMIN,NfMAX
        GammaCuspQ_internal(1,n)=4d0*CF*(CA*(67d0/9d0 -2*zeta2) -10d0/9d0*n)
    end do
    
    !! 3-loop
    do n=NfMIN,NfMAX
        GammaCuspQ_internal(2,n)=4d0*CF*(&
            (CA**2)*(245d0/6d0 -268d0/9d0*zeta2 +22d0/3d0*zeta3 +22d0*zeta4)&
            +CA*n*(-209d0/27d0 +40d0/9d0*zeta2 -28d0/3d0*zeta3)&
            +CF*n*(-55d0/6d0 +8d0*zeta3) -(n**2)*4d0/27d0)
    end do
    
    !! 4-loop
    do n=NfMIN,NfMAX
        GammaCuspQ_internal(3,n)=CF*(&
            dFA*(7040d0/3d0*zeta5 +256d0/3d0*zeta3-768d0*(zeta3**2) -256d0*zeta2 -15872d0/35d0*(zeta2**3))&
            +dFF*n*(-2560d0/3d0*zeta5 -512d0/3d0*zeta3 +512d0*zeta2)&
            +(n**3)*(-32d0/81d0 +64d0/27d0*zeta3)&
            +CF*(n**2)*(2392d0/81d0 -640d0/9d0*zeta3 +64d0/5d0*(zeta2**2))&
            +(CF**2)*n*(572d0/9d0 -320d0*zeta5 +592d0/3d0*zeta3)&
            +CA*(n**2)*(923d0/81d0 +2240d0/27d0*zeta3 -608d0/81d0*zeta2 -224d0/15d0*(zeta2**2))&
            +CF*CA*n*(-34066d0/81d0 +160d0*zeta5 +3712d0/9d0*zeta3 +440d0/3d0*zeta2 -128*zeta2*zeta3 -352d0/5d0*(zeta2**2))&
            +(CA**2)*n*(-24137d0/81d0 +2096d0/9d0*zeta5 -23104d0/27d0*zeta3 +20320d0/81d0*zeta2 &
                +448d0/3d0*zeta2*zeta3 -352d0/15d0*(zeta2**2))&
            +(CA**3)*(84278d0/81d0 - 3608d0/9d0*zeta5 +20944d0/27d0*zeta3 -16d0*(zeta3**2) &
                -88400d0/81d0*zeta2 -352d0/3d0*zeta2*zeta3 +3608d0/5d0*(zeta2**2)-20032d0/105d0*(zeta2**3)))
    end do
    
    !! 5-loop
    do n=NfMIN,NfMAX
        if(n<=3) then
            GammaCuspQ_internal(4,n)=4d0*CF*32417.5d0
        else if(n==4) then
            GammaCuspQ_internal(4,n)=4d0*CF*19949.2d0
        else
            GammaCuspQ_internal(4,n)=4d0*CF*12468.3d0
        end if
    end do    
        
end subroutine SetGammaCuspQuark

!!!! sets the values of cusp anomalous dimension
!!!! the values up to 3-loop are the same as for the quark *CA/CF
!!!! the explicit 4-loop expression is not given yet.
!!!! Instead, I use the numeric values for given in [1805.09638]
!!!! the 5-loop expression is taken as CA/CF from [1812.11818] formula (13), central value(!!) 
!!!!      only nf=3,4,5 avalible,
subroutine SetGammaCuspGluon()
    integer::n
    
    !! 1-loop
    do n=NfMIN,NfMAX
        GammaCuspG_internal(0,n)=CA*4d0
    end do
    
    !! 2-loop
    do n=NfMIN,NfMAX
        GammaCuspG_internal(1,n)=4d0*CA*(CA*(67d0/9d0 -2*zeta2) -10d0/9d0*n)
    end do
    
    !! 3-loop
    do n=NfMIN,NfMAX
        GammaCuspG_internal(2,n)=4d0*CA*(&
            (CA**2)*(245d0/6d0 -268d0/9d0*zeta2 +22d0/3d0*zeta3 +22d0*zeta4)&
            +CA*n*(-209d0/27d0 +40d0/9d0*zeta2 -28d0/3d0*zeta3)&
            +CF*n*(-55d0/6d0 +8d0*zeta3) -(n**2)*4d0/27d0)
    end do
    
    !! 4-loop
    do n=NfMIN,NfMAX
        GammaCuspG_internal(3,n)=40880.d0 - 11714.d0*n + 440.0488d0*(n**2) + 7.362774d0*(n**3)
    end do
    
    !! 5-loop
    do n=NfMIN,NfMAX
        if(n<=3) then
            GammaCuspG_internal(4,n)=4d0*CA*32417.5d0
        else if(n==4) then
            GammaCuspG_internal(4,n)=4d0*CA*19949.2d0
        else
            GammaCuspG_internal(4,n)=4d0*CA*12468.3d0
        end if
    end do
    
end subroutine SetGammaCuspGluon

!!!! sets the values of vector FF anomalous dimension
!!!! it is defined with factor 2 (standard for TMD physics)
!!!! 3-loop expression is taken from [1004.3653]
!!!! 4-loop expression is taken from [2202.04660] (aux-file, 2,3-loops confirmed)
subroutine SetGammaVQuark()
    integer::n
    
    !! 1-loop
    do n=NfMIN,NfMAX
        GammaVQ_internal(1,n)=-6d0*CF
    end do
    
    !! 2-loop
    do n=NfMIN,NfMAX
        GammaVQ_internal(2,n)=2d0*CF*(&
            CF*(-3d0/2d0 +12d0*zeta2 -24d0*zeta3)&
            +CA*(-961d0/54d0 -11d0*zeta2 +26d0*zeta3)&
            +n*(65d0/27d0 +2d0*zeta2))
    end do
    
    !! 3-loop
    do n=NfMIN,NfMAX
        GammaVQ_internal(3,n)=2d0*CF*(&
            CF*n*(2953d0/54d0 -26d0/3d0*zeta2 -140d0/3d0*zeta4 +256d0/9d0*zeta3)&
            +(n**2)*(2417d0/729d0 -20d0/9d0*zeta2 -8d0/27d0*zeta3)&
            +CA*n*(-8659d0/729d0 +2594d0/81d0*zeta2 +22d0*zeta4 -964d0/27d0*zeta3)&
            +(CF**2)*(-29d0/2d0 -18d0*zeta2 -144d0*zeta4 -68d0*zeta3 +32d0*zeta2*zeta3 +240d0*zeta5)&
            +CF*CA*(-151d0/4d0 +410d0/3d0*zeta2 +494d0/3d0*zeta4 -844d0/3d0*zeta3 -16d0*zeta2*zeta3 -120d0*zeta5)&
            +(CA**2)*(-139345d0/2916d0 -7163d0/81d0*zeta2 -83d0*zeta4 +3526d0/9d0*zeta3 -88d0/3d0*zeta2*zeta3 -136d0*zeta5))
    end do
    
        !! 4-loop
    do n=NfMIN,NfMAX
        GammaVQ_internal(4,n)=-24503.6d0 + 7857.17d0*n - 414.043d0*n**2 - 6.51577d0*n**3
    end do 
    
end subroutine SetGammaVQuark

!!!! sets the values of vector FF anomalous dimension
!!!! it is defined with factor 2 (standard for TMD physics)
!!!! 3-loop expression is taken from [1004.3653]
!!!! 4-loop expression is taken from [2202.04660] (aux-file, 2,3-loops confirmed)
subroutine SetGammaVGluon()
    integer::n
    
    !! 1-loop
    do n=NfMIN,NfMAX
        GammaVG_internal(1,n)=-22d0/3d0*CA+4d0/3d0*n
    end do
    
    !! 2-loop
    do n=NfMIN,NfMAX
        GammaVG_internal(2,n)=2d0*(&
            (CA**2)*(-692d0/27d0 +11d0/3d0*zeta2+2d0*zeta3) +CA*n*(128d0/27d0 -2d0/3d0*zeta2) +2*CF*n)
    end do
    
    !! 3-loop
    do n=NfMIN,NfMAX
        GammaVG_internal(3,n)=2d0*(&
            (CA**3)*(-97186d0/729d0 +6109d0/81d0*zeta2 -319d0/3d0*zeta4 +122d0/3d0*zeta3 -40d0/3d0*zeta2*zeta3 -16d0*zeta5)&
            +(CA**2)*n*(30715d0/1458d0 -1198d0/81d0*zeta2 +82d0/3d0*zeta4 +356d0/27d0*zeta3)&
            +CF*CA*n*(1217d0/27d0 -2d0*zeta2 -8d0*zeta4 -152d0/9d0*zeta3)&
            +CA*(n**2)*(-269d0/1458d0 +20d0/27d0*zeta2 -56d0/27d0*zeta3)&
            -(CF**2)*n -11d0/9d0*CF*(n**2))
    end do    
    
    !! 4-loop
    do n=NfMIN,NfMAX
        GammaVG_internal(4,n)=-84790.4d0 + 29077.5d0*n - 906.871d0*n**2 - 2.90651d0*n**3
    end do 
    
end subroutine SetGammaVGluon

!!!! sets the values of finite part for RAD
!!!! Expression is taken from [1707.07606]
!!!! 4-loop expression is taken from aux-files in [2205.02249]
subroutine SetDn0Quark()
    integer::n
    
    !! 1-loop
    do n=NfMIN,NfMAX
        d_nk_Q_internal(1,0,n)=0d0
    end do
    
    !! 2-loop
    do n=NfMIN,NfMAX
        d_nk_Q_internal(2,0,n)=CF*CA*(404d0/27d0-14d0*zeta3)-CF*n*56d0/27d0
    end do
    
    !! 3-loop
    do n=NfMIN,NfMAX
        d_nk_Q_internal(3,0,n)=CF*(CA**2)*(297029d0/1458d0 -3196d0/81d0*zeta2 &
                                -6164d0/27d0*zeta3 -77d0/3d0*zeta4 +88d0/3d0*zeta2*zeta3 +96d0*zeta5)&
                            +CF*CA*n*(-31313d0/729d0 +412d0/81d0*zeta2 +452d0/27d0*zeta3 -10d0/3d0*zeta4)&
                            +(CF**2)*n*(-1711d0/54d0 +152d0/9d0*zeta3 +8d0*zeta4)&
                            +CF*(n**2)*(928d0/729d0 +16d0/9d0*zeta3)
    end do    
    
    !! 4-loop
    do n=NfMIN,NfMAX
        d_nk_Q_internal(4,0,n)=350.834d0 - 2428.14d0*n + 378.306d0*n**2 - 8.07192d0*n**3
    end do  
    
end subroutine SetDn0Quark

!!!! sets the values of finite part for RAD
!!!! Expression is taken from [1707.07606]
!!!! at this order the differance between quark and gluon is CA/CF
!!!! 4-loop expression is taken from aux-files in [2205.02249] (there is no CA/CF any more)
subroutine SetDn0Gluon()
    integer::n
    
    !! 1-loop
    do n=NfMIN,NfMAX
        d_nk_G_internal(1,0,n)=0d0
    end do
    
    !! 2-loop
    do n=NfMIN,NfMAX
        d_nk_G_internal(2,0,n)=(CA**2)*(404d0/27d0-14d0*zeta3)-CA*n*56d0/27d0
    end do
    
    !! 3-loop
    do n=NfMIN,NfMAX
        d_nk_G_internal(3,0,n)=(CA**3)*(297029d0/1458d0 -3196d0/81d0*zeta2 &
                                -6164d0/27d0*zeta3 -77d0/3d0*zeta4 +88d0/3d0*zeta2*zeta3 +96d0*zeta5)&
                            +(CA**2)*n*(-31313d0/729d0 +412d0/81d0*zeta2 +452d0/27d0*zeta3 -10d0/3d0*zeta4)&
                            +CF*CA*n*(-1711d0/54d0 +152d0/9d0*zeta3 +8d0*zeta4)&
                            +CA*(n**2)*(928d0/729d0 +16d0/9d0*zeta3)
    end do 
    
    !! 4-loop
    do n=NfMIN,NfMAX
        d_nk_G_internal(3,0,n)=-333.77d0 - 5506.38d0*n + 851.188d0*n**2 - 18.1618d0*n**3
    end do
    
end subroutine SetDn0Gluon
