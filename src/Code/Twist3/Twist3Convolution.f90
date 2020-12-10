!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 2.05
!
!	This file contains the part of the code, which is common for all TMD-evaluation modules that operates at twist-2
!	It shares the variables with the module where it is inlucded (as a text)
!	Such idiotic structure is needed since, FORTRAN does not allow inheritance.
!
!	Be AWARE of possible clash of variable names.
!
!	SO FAR THERE IS NO TWIST-3 CONVOLUTION IMPLEMENTATION
!
!				A.Vladimirov (21.05.2020)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Supplimentary functions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!Evaluate Log[mu^2 bT^2/4/Exp(-2Gamma_E)]
!! the b here is b*, this funciton is used only in Coefficeint functions
function LogMuB(mu,bT)
    real(dp)::LogMuB
    real(dp),intent(in)::bT,mu
    LogMuB=2d0*Log(bSTAR(bT,lambdaNP)*mu*C0_inv_const)
end function LogMuB  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Convolutions!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
!- This is the TMD function evaluated for all quarks simultaniously (-5..5) at x,bT,mu
!--    f =  -5, -4, -3,  -2,  -1,0,1,2,3, 4 ,5
!--      = bbar ,cbar sbar,ubar,dbar,g,d,u,s, c ,b
!!---Gluon contribution is undefined
!- The order is accumulative pertrubative order of coefficient =0,1,2 (LO,NLO,NNLO)
!---------------------------------------------------------------------
function Common_lowScale5(x,bT,hadron)
    real(dp),dimension(-5:5)::Common_lowScale5
    real(dp),intent(in) :: x, bT
    integer,intent(in)::hadron

    !!!!SO FAR THERE IS NO TWIST-3 CONVOLUTION IMPLEMENTATION

    Common_lowScale5=FNP(x,bT,hadron,lambdaNP)

 end function Common_lowScale5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !---------------------------------------------------------------------
!- This is the TMD function evaluated for all quarks simultaniously (-5..5) at x,bT,mu the GLUON INCLUDED
!--    f =  -5, -4, -3,  -2,  -1,0,1,2,3, 4 ,5
!--      = bbar ,cbar sbar,ubar,dbar,g,d,u,s, c ,b
!- The order is accumulative pertrubative order of coefficient =0,1,2 (LO,NLO,NNLO)
!---------------------------------------------------------------------
function Common_lowScale50(x,bT,hadron)
    real(dp),dimension(-5:5)::Common_lowScale50
    real(dp),intent(in) :: x, bT
    integer,intent(in)::hadron

    !!!!SO FAR THERE IS NO TWIST-3 CONVOLUTION IMPLEMENTATION

    Common_lowScale50=FNP(x,bT,hadron,lambdaNP)

end function Common_lowScale50
