!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 2.03
!
!	Contains the expressions for coefficients of anomalous dimensions used by TMDR module.
!   Main purpose is to make evaluation of ADs fast and transparent.
!	
!	if you use this module please, quote 1803.11089
!
!				A.Vladimirov (10.02.2020)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module TMD_AD
use aTMDe_Numerics
use IO_functions

implicit none

private

 !Current version of module
 character (len=5),parameter :: version="v2.03"
 character (len=7),parameter :: moduleName="TMD_AD"
 
 logical::started=.false.

!!! Minimum and maximum avalible values of Nf (number of quark flavors)
integer,parameter::NfMIN=0
integer,parameter::NfMAX=6

!!! color factors  for SU(N=3)
real(dp),parameter::TF=1d0/2d0
real(dp),parameter::CF=4d0/3d0
real(dp),parameter::CA=3d0
real(dp),parameter::dAA=135d0/8d0     !d_A^{abcd}d_A^{abcd}/(N^2-1)
real(dp),parameter::dFA=15d0/16d0     !d_F^{abcd}d_A^{abcd}/(N^2-1)
real(dp),parameter::dFF=5d0/96d0      !d_F^{abcd}d_F^{abcd}/(N^2-1)

!!--------------------------------------------------
!! I precompute values of AD on initialization stage
!! and same them to matrices (internal to module).
!! A pack of functions return the values of matrices.
!!
!! There are two classes of anomalous dimensions primary and secondary.
!! Primary:
!! These are set up to the best known order
!! Secondary:
!! These are derived from primary ADs, up to requested number of loops

!!! beta function (up to 5 loops) (numeration starts from 0)
real(dp),dimension(0:4,NfMIN:NfMAX)::beta_internal

!!! gamma cusp for quark/gluon (up to 4 loops) (numeration starts from 0)
real(dp),dimension(0:3,NfMIN:NfMAX)::GammaCuspQ_internal
real(dp),dimension(0:3,NfMIN:NfMAX)::GammaCuspG_internal


!!! gammaV for quark/gluon (up to 3 loops) (numeration starts from 1)
real(dp),dimension(1:3,NfMIN:NfMAX)::GammaVQ_internal
real(dp),dimension(1:3,NfMIN:NfMAX)::GammaVG_internal


!!! d^{(n,k)} for quark/gluon (up to 3 loops) (numeration starts from 1)
real(dp),dimension(1:3,0:3,NfMIN:NfMAX)::d_nk_Q_internal
real(dp),dimension(1:3,0:3,NfMIN:NfMAX)::d_nk_G_internal


!!! d^{(n,k,l)} (resummed coefficients) for quark/Gluon (up to 3 loops/no tree)(numeration starts from 0/)
real(dp),dimension(1:3,0:3,0:3,NfMIN:NfMAX)::d_nkl_Q_internal
real(dp),dimension(1:3,0:3,0:3,NfMIN:NfMAX)::d_nkl_G_internal

!!! v^{(n,k)} for quark/Gluon (up to 3 loops )(numeration starts from 0)
real(dp),dimension(0:3,0:4,NfMIN:NfMAX)::v_nk_Q_internal
real(dp),dimension(0:3,0:4,NfMIN:NfMAX)::v_nk_G_internal

!!! OMEGA^{(n,k)} for quark/Gluon (up to 3 loops )(numeration starts from 0)
!!! Each terms is multiplied by its own functional coefficeint.
real(dp),dimension(0:3,1:4,NfMIN:NfMAX)::OMEGA_nk_Q_internal
real(dp),dimension(0:3,1:4,NfMIN:NfMAX)::OMEGA_nk_G_internal
!!! 2 beta0/Gamma0
real(dp),dimension(NfMIN:NfMAX)::pFACTOR_Q_internal
real(dp),dimension(NfMIN:NfMAX)::pFACTOR_G_internal

public:: TMD_AD_Initialize
public:: betaQCD,GammaCusp_q,GammaCusp_g,gammaV_q, gammaV_g,dnk_q,dnk_g,dnkl_q,dnkl_g
public:: vnk_q,vnk_g,OMEGA_q,OMEGA_g,pFACTOR_q,pFACTOR_g

contains

subroutine TMD_AD_Initialize()
    
    if(started) return

    !----------Primary coefficients
    call SetBetaFunction()
    call SetGammaCuspQuark()
    call SetGammaCuspGluon()
    call SetGammaVQuark()
    call SetGammaVGluon()
    call SetDn0Quark()
    call SetDn0Gluon()
    
    !---------Secondary coefficients
    call SetDnkQuark()
    call SetDnkGluon()
    call SetDnklQuark()
    call SetDnklGluon()
    call SetVnkQuark()
    call SetVnkGluon()
    call SetOMEGAnkQuark()
    call SetOMEGAnkGluon()
    
    started=.true.

end subroutine TMD_AD_Initialize

!!!-------------------------------------------------------------------------------------------------------
!!!--------------------------------------------PUBLIC INTERFACES------------------------------------------
!!!-------------------------------------------------------------------------------------------------------

!!Coefficients of QCD beta function
function betaQCD(n,Nf)
    real(dp)::betaQCD
    integer::n,Nf
    betaQCD=beta_internal(n,Nf)
end function betaQCD

!!Light-like cusp coefficients QUARK
function GammaCusp_q(n,Nf)
    real(dp)::GammaCusp_q
    integer::n,Nf
    GammaCusp_q=GammaCuspQ_internal(n,Nf)
end function GammaCusp_q

!!Light-like cusp coefficients GLUON
function GammaCusp_g(n,Nf)
    real(dp)::GammaCusp_g
    integer::n,Nf
    GammaCusp_g=GammaCuspG_internal(n,Nf)
end function GammaCusp_g

!!gammaV coefficients QUARK
function gammaV_q(n,Nf)
    real(dp)::gammaV_q
    integer::n,Nf
    gammaV_q=GammaVQ_internal(n,Nf)
end function gammaV_q

!!gammaV coefficients GLUON
function gammaV_g(n,Nf)
    real(dp)::gammaV_g
    integer::n,Nf
    gammaV_g=GammaVG_internal(n,Nf)
end function gammaV_g

!! RAD coefficients QUARK
function dnk_q(n,k,Nf)
    real(dp)::dnk_q
    integer::n,k,Nf
    
    dnk_q=d_nk_Q_internal(n,k,Nf)

end function dnk_q

!! RAD coefficients GLUON
function dnk_g(n,k,Nf)
    real(dp)::dnk_g
    integer::n,k,Nf
    
    dnk_g=d_nk_G_internal(n,k,Nf)

end function dnk_g

!! RAD coefficients resummed QUARK
function dnkl_q(n,k,l,Nf)
    real(dp)::dnkl_q
    integer::n,k,Nf,l
    
    dnkl_q=d_nkl_Q_internal(n,k,l,Nf)

end function dnkl_q

!! RAD coefficients resummed GLUON
function dnkl_g(n,k,l,Nf)
    real(dp)::dnkl_g
    integer::n,k,Nf,l
    
    dnkl_g=d_nkl_G_internal(n,k,l,Nf)

end function dnkl_g

!! RAD coefficients QUARK
function vnk_q(n,k,Nf)
    real(dp)::vnk_q
    integer::n,k,Nf
    
    vnk_q=v_nk_Q_internal(n,k,Nf)

end function vnk_q

!! RAD coefficients GLUON
function vnk_g(n,k,Nf)
    real(dp)::vnk_g
    integer::n,k,Nf
    
    vnk_g=v_nk_G_internal(n,k,Nf)

end function vnk_g

!! RAD resummed coefficients QUARK
function OMEGA_q(n,k,Nf)
    real(dp)::OMEGA_q
    integer::n,k,Nf
    
    OMEGA_q=OMEGA_nk_Q_internal(n,k,Nf)

end function OMEGA_q

!! RAD resummed coefficients GLUON
function OMEGA_g(n,k,Nf)
    real(dp)::OMEGA_g
    integer::n,k,Nf
    
    OMEGA_g=OMEGA_nk_G_internal(n,k,Nf)

end function OMEGA_g

!! coefficeint for p QUARK
function pFACTOR_q(Nf)
    real(dp)::pFACTOR_q
    integer::n,k,Nf
    
    pFACTOR_q=pFACTOR_Q_internal(Nf)

end function pFACTOR_q

!! coefficeint for p GLUON
function pFACTOR_g(Nf)
    real(dp)::pFACTOR_g
    integer::n,k,Nf
    
    pFACTOR_g=pFACTOR_G_internal(Nf)

end function pFACTOR_g


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
    
end subroutine SetGammaCuspQuark

!!!! sets the values of cusp anomalous dimension
!!!! the values up to 3-loop are the same as for the quark *CA/CF
!!!! the explicit 4-loop expression is not given yet.
!!!! Instead, I use the numeric values for given in [1805.09638]
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
    
end subroutine SetGammaCuspGluon

!!!! sets the values of vector FF anomalous dimension
!!!! it is defined with factor 2 (standard for TMD physics)
!!!! 3-loop expression is taken from [1004.3653]
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
    
    
end subroutine SetGammaVQuark

!!!! sets the values of vector FF anomalous dimension
!!!! it is defined with factor 2 (standard for TMD physics)
!!!! 3-loop expression is taken from [1004.3653]
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
    
end subroutine SetGammaVGluon

!!!! sets the values of finite part for RAD
!!!! Expression is taken from [1707.07606]
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
    
end subroutine SetDn0Quark

!!!! sets the values of finite part for RAD
!!!! Expression is taken from [1707.07606]
!!!! at this order the differance between quark and gluon is CA/CF
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
    
end subroutine SetDn0Gluon

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

end module TMD_AD
