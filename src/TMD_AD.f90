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
use QCDinput

implicit none

private

 !Current version of module
 character (len=5),parameter :: version="v2.06"
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
!! A pack of get-functions returns the values of matrices.
!!
!! There are two classes of anomalous dimensions primary and secondary.
!! Primary:
!! These are set up to the best known order
!! Secondary:
!! These are derived from primary ADs, up to requested number of loops

!!! beta function (up to 5 loops) (numeration starts from 0)
real(dp),dimension(0:4,NfMIN:NfMAX)::beta_internal

!!! gamma cusp for quark/gluon (up to 5 loops) (numeration starts from 0)
real(dp),dimension(0:4,NfMIN:NfMAX)::GammaCuspQ_internal
real(dp),dimension(0:4,NfMIN:NfMAX)::GammaCuspG_internal


!!! gammaV for quark/gluon (up to 4 loops) (numeration starts from 1)
real(dp),dimension(1:4,NfMIN:NfMAX)::GammaVQ_internal
real(dp),dimension(1:4,NfMIN:NfMAX)::GammaVG_internal

!!!! Roots of beta-function betaRoots(n,k,Nf) \beta(a)=0
!!!! Depending on Nf and k, these can be complex or real. Each case special
COMPLEX*16,dimension(1:4,1:4,NfMIN:NfMAX)::betaRoots_internal
!!!! coefficients of integral Gamma/2
COMPLEX*16,dimension(1:4,0:4,NfMIN:NfMAX)::GammaIntegral_Q_internal
COMPLEX*16,dimension(1:4,0:4,NfMIN:NfMAX)::GammaIntegral_G_internal

!!! d^{(n,k)} for quark/gluon (up to 4 loops) (numeration starts from 1)
real(dp),dimension(1:4,0:4,NfMIN:NfMAX)::d_nk_Q_internal
real(dp),dimension(1:4,0:4,NfMIN:NfMAX)::d_nk_G_internal


!!! d^{(n,k,l)} (resummed coefficients) for quark/Gluon (up to 4 loops/no tree)(numeration starts from 0/)
real(dp),dimension(1:4,0:4,0:4,NfMIN:NfMAX)::d_nkl_Q_internal
real(dp),dimension(1:4,0:4,0:4,NfMIN:NfMAX)::d_nkl_G_internal

!!! v^{(n,k)} for quark/Gluon (up to 4 loops )(numeration starts from 0)
real(dp),dimension(0:4,0:5,NfMIN:NfMAX)::v_nk_Q_internal
real(dp),dimension(0:4,0:5,NfMIN:NfMAX)::v_nk_G_internal

!!! OMEGA^{(n,k)} for quark/Gluon (up to 5 loops )(numeration starts from 0)
!!! Each terms is multiplied by its own functional coefficeint.
real(dp),dimension(0:4,1:5,NfMIN:NfMAX)::OMEGA_nk_Q_internal
real(dp),dimension(0:4,1:5,NfMIN:NfMAX)::OMEGA_nk_G_internal
!!! 2 beta0/Gamma0
real(dp),dimension(NfMIN:NfMAX)::pFACTOR_Q_internal
real(dp),dimension(NfMIN:NfMAX)::pFACTOR_G_internal

!!!! the orders of used perturbative expressions for AD
!!!! the values of these constants are set thorugh initialization routine
integer::orderCusp      !for Cusp AD
integer::orderV         !for gammaV
integer::orderD         !for RAD
integer::orderDresum    !for resummed RAD
integer::orderZETA      !for zeta-line

public:: TMD_AD_Initialize
!public:: betaQCD,GammaCusp_q,GammaCusp_g,gammaV_q, gammaV_g,dnk_q,dnk_g,dnkl_q,dnkl_g
!public:: vnk_q,vnk_g,OMEGA_q,OMEGA_g,pFACTOR_q,pFACTOR_g
public:: zetaMUpert,gammaV,Dpert,Dresum,GammaCusp,zetaSL,RADEvolution,zFUNC

contains

!!! definitions of coefficients of primary ADs (Gamma, beta, gammaV, d)
INCLUDE 'Code/TMD_AD/AD_primary.f90'
!!! definitions of coefficients of secondary ADs (RAD, zeta-lines, resummed versions)
INCLUDE 'Code/TMD_AD/AD_secondary.f90'
!!! definitions of coefficients of secondary ADs at scale mu
INCLUDE 'Code/TMD_AD/AD_atMu.f90'
!!! Routines for the exact zeta-line
INCLUDE 'Code/TMD_AD/exactZetaLine.f90'
!!! Routines for the analytical evaluation of RGE-integrals
INCLUDE 'Code/TMD_AD/AD_Integral.f90'

subroutine TMD_AD_Initialize(oCusp,oV,oD,oDresum,oZETA)
    integer,intent(in)::oCusp,oV,oD,oDresum,oZETA
    
    real(dp)::aa
    
    if(started) return
    
    !setting orders
    orderCusp=oCusp
    orderV=oV
    orderD=oD
    orderDresum=oDresum
    orderZETA=oZETA

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
    
    !---------Set values for integrals
    call SetBetaRoots()
    call SetIntegralCoefficeintsGAMMA_Q()
    call SetIntegralCoefficeintsGAMMA_G()   
    
    started=.true.

end subroutine TMD_AD_Initialize

!!!-------------------------------------------------------------------------------------------------------
!!!--------------------------------------------GET INTERFACES---------------------------------------------
!!!-------------------------------------------------------------------------------------------------------

!!Coefficients of QCD beta function
pure function betaQCD(n,Nf)
    real(dp)::betaQCD
    integer,intent(in)::n,Nf
    betaQCD=beta_internal(n,Nf)
end function betaQCD

!!Light-like cusp coefficients QUARK
pure function GammaCusp_q(n,Nf)
    real(dp)::GammaCusp_q
    integer,intent(in)::n,Nf
    GammaCusp_q=GammaCuspQ_internal(n,Nf)
end function GammaCusp_q

!!Light-like cusp coefficients GLUON
pure function GammaCusp_g(n,Nf)
    real(dp)::GammaCusp_g
    integer,intent(in)::n,Nf
    GammaCusp_g=GammaCuspG_internal(n,Nf)
end function GammaCusp_g

!!gammaV coefficients QUARK
pure function gammaV_q(n,Nf)
    real(dp)::gammaV_q
    integer,intent(in)::n,Nf
    gammaV_q=GammaVQ_internal(n,Nf)
end function gammaV_q

!!gammaV coefficients GLUON
pure function gammaV_g(n,Nf)
    real(dp)::gammaV_g
    integer,intent(in)::n,Nf
    gammaV_g=GammaVG_internal(n,Nf)
end function gammaV_g

!! RAD coefficients QUARK
pure function dnk_q(n,k,Nf)
    real(dp)::dnk_q
    integer,intent(in)::n,k,Nf    
    dnk_q=d_nk_Q_internal(n,k,Nf)
end function dnk_q

!! RAD coefficients GLUON
pure function dnk_g(n,k,Nf)
    real(dp)::dnk_g
    integer,intent(in)::n,k,Nf    
    dnk_g=d_nk_G_internal(n,k,Nf)
end function dnk_g

!! RAD coefficients resummed QUARK
pure function dnkl_q(n,k,l,Nf)
    real(dp)::dnkl_q
    integer,intent(in)::n,k,Nf,l    
    dnkl_q=d_nkl_Q_internal(n,k,l,Nf)
end function dnkl_q

!! RAD coefficients resummed GLUON
pure function dnkl_g(n,k,l,Nf)
    real(dp)::dnkl_g
    integer,intent(in)::n,k,Nf,l    
    dnkl_g=d_nkl_G_internal(n,k,l,Nf)
end function dnkl_g

!! RAD coefficients QUARK
pure function vnk_q(n,k,Nf)
    real(dp)::vnk_q
    integer,intent(in)::n,k,Nf    
    vnk_q=v_nk_Q_internal(n,k,Nf)
end function vnk_q

!! RAD coefficients GLUON
pure function vnk_g(n,k,Nf)
    real(dp)::vnk_g
    integer,intent(in)::n,k,Nf    
    vnk_g=v_nk_G_internal(n,k,Nf)
end function vnk_g

!! RAD resummed coefficients QUARK
pure function OMEGA_q(n,k,Nf)
    real(dp)::OMEGA_q
    integer,intent(in)::n,k,Nf    
    OMEGA_q=OMEGA_nk_Q_internal(n,k,Nf)
end function OMEGA_q

!! RAD resummed coefficients GLUON
pure function OMEGA_g(n,k,Nf)
    real(dp)::OMEGA_g
    integer,intent(in)::n,k,Nf    
    OMEGA_g=OMEGA_nk_G_internal(n,k,Nf)
end function OMEGA_g

!! coefficeint for p QUARK
pure function pFACTOR_q(Nf)
    real(dp)::pFACTOR_q
    integer,intent(in)::Nf    
    pFACTOR_q=pFACTOR_Q_internal(Nf)
end function pFACTOR_q

!! coefficeint for p GLUON
pure function pFACTOR_g(Nf)
    real(dp)::pFACTOR_g
    integer,intent(in)::Nf    
    pFACTOR_g=pFACTOR_G_internal(Nf)
end function pFACTOR_g

end module TMD_AD
