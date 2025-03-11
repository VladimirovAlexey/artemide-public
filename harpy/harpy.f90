!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			harpy
!
!	An interface for artemide
!
!				A.Vladimirov (13.01.2019)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module harpy
use TMDX_DY
use TMDX_SIDIS
use aTMDe_control
use uTMDPDF
use uTMDPDF_OPE
use uTMDFF
use uTMDFF_OPE
use SiversTMDPDF
use lpTMDPDF
use wgtTMDPDF
use BoerMuldersTMDPDF
use TMDR_model
use TMDR

!!! this flag is requared to guaranty that artemide is not started twice (it lead to the crush)
logical::started=.false.

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!GENERAL
  subroutine Initialize(file)
    character(len=*)::file
    if(started) then
      write(*,*) 'artemide already runs'
    else
      call artemide_Initialize(file)
      started=.true.
    end if
  end subroutine Initialize
  
  
  subroutine ShowStatistics()
   call artemide_ShowStatistics()
  end  subroutine ShowStatistics

  !!!Sets the non-perturbative parameters lambda
  subroutine SetLambda_Main(lambdaIN)
    real*8,intent(in)::lambdaIN(:)
    call artemide_SetNPparameters(lambdaIN)
  end subroutine SetLambda_Main  
  
    !!!Sets the non-perturbative parameters lambda
  subroutine SetLambda_TMDR(lambdaIN)
    real*8,intent(in)::lambdaIN(:)
    call artemide_SetNPparameters_TMDR(lambdaIN)
  end subroutine SetLambda_TMDR
  
      !!!Sets the non-perturbative parameters lambda
  subroutine SetLambda_uTMDPDF(lambdaIN)
    real*8,intent(in)::lambdaIN(:)
    call artemide_SetNPparameters_uTMDPDF(lambdaIN)
  end subroutine SetLambda_uTMDPDF
  
  !!!Sets the non-perturbative parameters lambda
  subroutine SetLambda_uTMDFF(lambdaIN)
    real*8,intent(in)::lambdaIN(:)
    call artemide_SetNPparameters_uTMDFF(lambdaIN)
  end subroutine SetLambda_uTMDFF
  
  !!!Sets the non-perturbative parameters lambda
  subroutine SetLambda_lpTMDPDF(lambdaIN)
    real*8,intent(in)::lambdaIN(:)
    call artemide_SetNPparameters_lpTMDPDF(lambdaIN)
  end subroutine SetLambda_lpTMDPDF
  
  !!!Sets the non-perturbative parameters lambda
  subroutine SetLambda_SiversTMDPDF(lambdaIN)
    real*8,intent(in)::lambdaIN(:)
    call artemide_SetNPparameters_SiversTMDPDF(lambdaIN)
  end subroutine SetLambda_SiversTMDPDF
  
    !!!Sets the non-perturbative parameters lambda
  subroutine SetLambda_wgtTMDPDF(lambdaIN)
    real*8,intent(in)::lambdaIN(:)
    call artemide_SetNPparameters_wgtTMDPDF(lambdaIN)
  end subroutine SetLambda_wgtTMDPDF

      !!!Sets the non-perturbative parameters lambda
  subroutine SetLambda_BoerMuldersTMDPDF(lambdaIN)
    real*8,intent(in)::lambdaIN(:)
    call artemide_SetNPparameters_BoerMuldersTMDPDF(lambdaIN)
  end subroutine SetLambda_BoerMuldersTMDPDF
  
  
  !!!! this routine set the variations of scales
  !!!! it is used for the estimation of errors
  subroutine SetScaleVariation(c1_in,c2_in,c3_in,c4_in)
    real*8::c1_in,c2_in,c3_in,c4_in
    
    call artemide_SetScaleVariations(c1_in,c2_in,c3_in,c4_in)
    
  end subroutine SetScaleVariation
  
  !! reset the number for PDF replica for uTMDPDF
  subroutine SetPDFreplica(rep,hadron)
    integer::rep,hadron
    call uTMDPDF_SetPDFreplica(rep,hadron)
  end subroutine SetPDFreplica
  
    !! reset the number for PDF replica for uTMDFF
  subroutine SetFFreplica(rep,hadron)
    integer::rep,hadron
    call uTMDFF_SetPDFreplica(rep,hadron)
  end subroutine SetFFreplica
  
    !! reset the number for PDF replica for lpTMDPDF
  subroutine SetlpPDFreplica(rep,hadron)
    integer::rep,hadron
    call lpTMDPDF_SetPDFreplica(rep,hadron)
  end subroutine SetlpPDFreplica
  
    !! reset the number for PDF replica for wgtTMDPDF
  subroutine SetwgtPDFreplica(rep,hadron)
    integer::rep,hadron
    call wgtTMDPDF_SetPDFreplica(rep,hadron)
  end subroutine SetwgtPDFreplica
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!! DNP
  function getDNP(b,mu,f)
    real*8::getDNP
    real*8::b,mu
    integer::f
    
    getDNP=CS_kernel(mu,b,f)
  end function getDNP

  !!!!! optimal R
  function getR(b,mu,zeta,f)
    real*8::TMDR_Rzeta
    real*8::b,mu,zeta
    integer::f

    getR=TMDR_Rzeta_harpy(b,mu,zeta,f)
  end function getR
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!uTMDPDF
  !!!!!!!! upolarized TMDPDF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function uTMDPDF_5_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: uTMDPDF_5_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  uTMDPDF_5_Evolved=uTMDPDF_inB(x,bt,muf,zetaf,hadron)
    
  end function uTMDPDF_5_Evolved
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDPDF_50_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: uTMDPDF_50_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  uTMDPDF_50_Evolved=uTMDPDF_inB(x,bt,muf,zetaf,hadron)
    
  end function uTMDPDF_50_Evolved

  !!!!!!!! upolarized TMDPDF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function uTMDPDF_5_Optimal(x,bt,hadron)
    real*8:: uTMDPDF_5_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  uTMDPDF_5_Optimal=uTMDPDF_inB(x,bt,hadron)
    
  end function uTMDPDF_5_Optimal
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDPDF_50_Optimal(x,bt,hadron)
    real*8:: uTMDPDF_50_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  uTMDPDF_50_Optimal=uTMDPDF_inB(x,bt,hadron)
    
  end function uTMDPDF_50_Optimal
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!uTMDFF
  !!!!!!!! upolarized TMDFF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function uTMDFF_5_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: uTMDFF_5_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  uTMDFF_5_Evolved=uTMDFF_inB(x,bt,muf,zetaf,hadron)
    
  end function uTMDFF_5_Evolved
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDFF_50_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: uTMDFF_50_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  uTMDFF_50_Evolved=uTMDFF_inB(x,bt,muf,zetaf,hadron)
    
  end function uTMDFF_50_Evolved

    !!!!!!!! upolarized TMDFF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function uTMDFF_5_Optimal(x,bt,hadron)
    real*8:: uTMDFF_5_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  uTMDFF_5_Optimal=uTMDFF_inB(x,bt,hadron)
    
  end function uTMDFF_5_Optimal
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDFF_50_Optimal(x,bt,hadron)
    real*8:: uTMDFF_50_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  uTMDFF_50_Optimal=uTMDFF_inB(x,bt,hadron)
    
  end function uTMDFF_50_Optimal
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SiversTMDPDF
  !!!!!!!! Sivers TMDFF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function SiversTMDPDF_5_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: SiversTMDPDF_5_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  SiversTMDPDF_5_Evolved=SiversTMDPDF_inB(x,bt,muf,zetaf,hadron)
    
  end function SiversTMDPDF_5_Evolved
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function SiversTMDPDF_50_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: SiversTMDPDF_50_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  SiversTMDPDF_50_Evolved=SiversTMDPDF_inB(x,bt,muf,zetaf,hadron)
    
  end function SiversTMDPDF_50_Evolved

    !!!!!!!! Sivers TMDPDF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function SiversTMDPDF_5_Optimal(x,bt,hadron)
    real*8:: SiversTMDPDF_5_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  SiversTMDPDF_5_Optimal=SiversTMDPDF_inB(x,bt,hadron)
    
  end function SiversTMDPDF_5_Optimal
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function SiversTMDPDF_50_Optimal(x,bt,hadron)
    real*8:: SiversTMDPDF_50_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  SiversTMDPDF_50_Optimal=SiversTMDPDF_inB(x,bt,hadron)
    
  end function SiversTMDPDF_50_Optimal
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!lpTMDPDF
  !!!!!!!! linearly polarized TMDFF
 
  ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function lpTMDPDF_50_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: lpTMDPDF_50_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  lpTMDPDF_50_Evolved=lpTMDPDF_inB(x,bt,muf,zetaf,hadron)
    
  end function lpTMDPDF_50_Evolved
  
  ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function lpTMDPDF_50_Optimal(x,bt,hadron)
    real*8:: lpTMDPDF_50_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  lpTMDPDF_50_Optimal=lpTMDPDF_inB(x,bt,hadron)
    
  end function lpTMDPDF_50_Optimal
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!wgtTMDPDF
  !!!!!!!! wgt TMDFF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function wgtTMDPDF_5_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: wgtTMDPDF_5_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  wgtTMDPDF_5_Evolved=wgtTMDPDF_inB(x,bt,muf,zetaf,hadron)
    
  end function wgtTMDPDF_5_Evolved
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function wgtTMDPDF_50_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: wgtTMDPDF_50_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  wgtTMDPDF_50_Evolved=wgtTMDPDF_inB(x,bt,muf,zetaf,hadron)
    
  end function wgtTMDPDF_50_Evolved

    !!!!!!!! wgt TMDPDF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function wgtTMDPDF_5_Optimal(x,bt,hadron)
    real*8:: wgtTMDPDF_5_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  wgtTMDPDF_5_Optimal=wgtTMDPDF_inB(x,bt,hadron)
    
  end function wgtTMDPDF_5_Optimal
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function wgtTMDPDF_50_Optimal(x,bt,hadron)
    real*8:: wgtTMDPDF_50_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  wgtTMDPDF_50_Optimal=wgtTMDPDF_inB(x,bt,hadron)
    
  end function wgtTMDPDF_50_Optimal

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!BoerMuldersTMDPDF
  !!!!!!!! BoerMulders TMDFF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function BoerMuldersTMDPDF_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: BoerMuldersTMDPDF_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron

  BoerMuldersTMDPDF_Evolved=BoerMuldersTMDPDF_inB(x,bt,muf,zetaf,hadron)

  end function BoerMuldersTMDPDF_Evolved

    !!!!!!!! BoerMulders TMDPDF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function BoerMuldersTMDPDF_Optimal(x,bt,hadron)
    real*8:: BoerMuldersTMDPDF_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron

  BoerMuldersTMDPDF_Optimal=BoerMuldersTMDPDF_inB(x,bt,hadron)

  end function BoerMuldersTMDPDF_Optimal
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! TMDs IN KT
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!! upolarized TMDPDF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function uTMDPDF_kT_5_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: uTMDPDF_kT_5_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  uTMDPDF_kT_5_Evolved=uTMDPDF_inKT(x,bt,muf,zetaf,hadron)
    
  end function uTMDPDF_kT_5_Evolved
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDPDF_kT_50_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: uTMDPDF_kT_50_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  uTMDPDF_kT_50_Evolved=uTMDPDF_inKT(x,bt,muf,zetaf,hadron)
    
  end function uTMDPDF_kT_50_Evolved

  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function uTMDPDF_kT_5_Optimal(x,bt,hadron)
    real*8:: uTMDPDF_kT_5_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  uTMDPDF_kT_5_Optimal=uTMDPDF_inKT(x,bt,hadron)
    
  end function uTMDPDF_kT_5_Optimal
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDPDF_kT_50_Optimal(x,bt,hadron)
    real*8:: uTMDPDF_kT_50_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  uTMDPDF_kT_50_Optimal=uTMDPDF_inKT(x,bt,hadron)
    
  end function uTMDPDF_kT_50_Optimal

  ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDPDF_G0(x,mu,hadron)
    real*8:: uTMDPDF_G0(-5:5)
    real*8:: x,mu
    integer::hadron

  uTMDPDF_G0=uTMDPDF_TMM_G(x,mu,hadron)

  end function uTMDPDF_G0

    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDPDF_X0(x,mu,hadron)
    real*8:: uTMDPDF_X0(-5:5)
    real*8:: x,mu
    integer::hadron

  uTMDPDF_X0=uTMDPDF_TMM_X(x,mu,hadron)

  end function uTMDPDF_X0

      ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDPDF_ASX0(x,mu,hadron)
    real*8:: uTMDPDF_ASX0(-5:5)
    real*8:: x,mu
    integer::hadron

  uTMDPDF_ASX0=uTMDPDF_X0_AS(x,mu,mu,hadron)

  end function uTMDPDF_ASX0

      ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDPDF_PDF(x,mu,hadron)
    real*8:: uTMDPDF_PDF(-5:5)
    real*8:: x,mu
    integer::hadron

  uTMDPDF_PDF=uTMDPDF_OPE_PDF(x,mu,hadron)

  end function uTMDPDF_PDF
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!! upolarized TMDFF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function uTMDFF_kT_5_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: uTMDFF_kT_5_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  uTMDFF_kT_5_Evolved=uTMDFF_inKT(x,bt,muf,zetaf,hadron)
    
  end function uTMDFF_kT_5_Evolved
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDFF_kT_50_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: uTMDFF_kT_50_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  uTMDFF_kT_50_Evolved=uTMDFF_inKT(x,bt,muf,zetaf,hadron)
    
  end function uTMDFF_kT_50_Evolved

  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function uTMDFF_kT_5_Optimal(x,bt,hadron)
    real*8:: uTMDFF_kT_5_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  uTMDFF_kT_5_Optimal=uTMDFF_inKT(x,bt,hadron)
    
  end function uTMDFF_kT_5_Optimal
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDFF_kT_50_Optimal(x,bt,hadron)
    real*8:: uTMDFF_kT_50_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  uTMDFF_kT_50_Optimal=uTMDFF_inKT(x,bt,hadron)
    
  end function uTMDFF_kT_50_Optimal

    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDFF_G0(x,mu,hadron)
    real*8:: uTMDFF_G0(-5:5)
    real*8:: x,mu
    integer::hadron

  uTMDFF_G0=uTMDFF_TMM_G(x,mu,hadron)

  end function uTMDFF_G0

    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDFF_X0(x,mu,hadron)
    real*8:: uTMDFF_X0(-5:5)
    real*8:: x,mu
    integer::hadron

  uTMDFF_X0=uTMDFF_TMM_X(x,mu,hadron)

  end function uTMDFF_X0

      ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDFF_ASX0(x,mu,hadron)
    real*8:: uTMDFF_ASX0(-5:5)
    real*8:: x,mu
    integer::hadron

  uTMDFF_ASX0=uTMDFF_X0_AS(x,mu,mu,hadron)

  end function uTMDFF_ASX0

      ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDFF_FF(x,mu,hadron)
    real*8:: uTMDFF_FF(-5:5)
    real*8:: x,mu
    integer::hadron

  uTMDFF_FF=uTMDFF_OPE_FF(x,mu,hadron)

  end function uTMDFF_FF
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!! linearly polarized gluon TMDPDF

  ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function lpTMDPDF_kT_50_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: lpTMDPDF_kT_50_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  lpTMDPDF_kT_50_Evolved=lpTMDPDF_inKT(x,bt,muf,zetaf,hadron)
    
  end function lpTMDPDF_kT_50_Evolved
  
  ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function lpTMDPDF_kT_50_Optimal(x,bt,hadron)
    real*8:: lpTMDPDF_kT_50_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  lpTMDPDF_kT_50_Optimal=lpTMDPDF_inKT(x,bt,hadron)
    
  end function lpTMDPDF_kT_50_Optimal
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!! Sivers TMDPDF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function SiversTMDPDF_kT_5_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: SiversTMDPDF_kT_5_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  SiversTMDPDF_kT_5_Evolved=SiversTMDPDF_inKT(x,bt,muf,zetaf,hadron)
    
  end function SiversTMDPDF_kT_5_Evolved
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function SiversTMDPDF_kT_50_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: SiversTMDPDF_kT_50_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  SiversTMDPDF_kT_50_Evolved=SiversTMDPDF_inKT(x,bt,muf,zetaf,hadron)
    
  end function SiversTMDPDF_kT_50_Evolved

  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function SiversTMDPDF_kT_5_Optimal(x,bt,hadron)
    real*8:: SiversTMDPDF_kT_5_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  SiversTMDPDF_kT_5_Optimal=SiversTMDPDF_inKT(x,bt,hadron)
    
  end function SiversTMDPDF_kT_5_Optimal
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function SiversTMDPDF_kT_50_Optimal(x,bt,hadron)
    real*8:: SiversTMDPDF_kT_50_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  SiversTMDPDF_kT_50_Optimal=SiversTMDPDF_inKT(x,bt,hadron)
    
  end function SiversTMDPDF_kT_50_Optimal

    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function SiversTMDPDF_G1(x,mu,hadron)
    real*8:: SiversTMDPDF_G1(-5:5)
    real*8:: x,mu
    integer::hadron

  SiversTMDPDF_G1=SiversTMDPDF_TMM_G(x,mu,hadron)

  end function SiversTMDPDF_G1
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!! wgt TMDPDF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function wgtTMDPDF_kT_5_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: wgtTMDPDF_kT_5_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  wgtTMDPDF_kT_5_Evolved=wgtTMDPDF_inKT(x,bt,muf,zetaf,hadron)
    
  end function wgtTMDPDF_kT_5_Evolved
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function wgtTMDPDF_kT_50_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: wgtTMDPDF_kT_50_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  wgtTMDPDF_kT_50_Evolved=wgtTMDPDF_inKT(x,bt,muf,zetaf,hadron)
    
  end function wgtTMDPDF_kT_50_Evolved

  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function wgtTMDPDF_kT_5_Optimal(x,bt,hadron)
    real*8:: wgtTMDPDF_kT_5_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  wgtTMDPDF_kT_5_Optimal=wgtTMDPDF_inKT(x,bt,hadron)
    
  end function wgtTMDPDF_kT_5_Optimal
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function wgtTMDPDF_kT_50_Optimal(x,bt,hadron)
    real*8:: wgtTMDPDF_kT_50_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  wgtTMDPDF_kT_50_Optimal=wgtTMDPDF_inKT(x,bt,hadron)
    
  end function wgtTMDPDF_kT_50_Optimal

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!! BoerMulders TMDPDF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function BoerMuldersTMDPDF_kT_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: BoerMuldersTMDPDF_kT_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron

  BoerMuldersTMDPDF_kT_Evolved=BoerMuldersTMDPDF_inKT(x,bt,muf,zetaf,hadron)

  end function BoerMuldersTMDPDF_kT_Evolved

  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function BoerMuldersTMDPDF_kT_Optimal(x,bt,hadron)
    real*8:: BoerMuldersTMDPDF_kT_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron

  BoerMuldersTMDPDF_kT_Optimal=BoerMuldersTMDPDF_inKT(x,bt,hadron)

  end function BoerMuldersTMDPDF_kT_Optimal
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!DY CROSS-SECTION

  
  function DY_xSec_Single(process,s,qT,Q,y,includeCuts,CutParameters)
    integer,intent(in),dimension(1:4)::process		!the number of process
    real*8,intent(in)::s				!Mandelshtam s
    real*8,intent(in),dimension(1:2)::qT		!(qtMin,qtMax)
    real*8,intent(in),dimension(1:2)::Q			!(Qmin,Qmax)
    real*8,intent(in),dimension(1:2)::y			!(ymin,ymax)
    logical,intent(in)::includeCuts			!include cuts
    real*8,intent(in),dimension(1:4)::CutParameters	!(p1,p2,eta1,eta2)
    real*8::DY_xSec_Single
    real*8::X
    
    call xSec_DY(X,process,s,qT,Q,y,includeCuts,CutParameters)
    DY_xSec_Single=X
  
  end function DY_xSec_Single
  
  function DY_xSec_List(process,s,qT,Q,y,includeCuts,CutParameters,ListLength)
    integer,intent(in)::ListLength
    integer,intent(in),dimension(:,:)::process			!the number of process
    real*8,intent(in),dimension(:)::s				!Mandelshtam s
    real*8,intent(in),dimension(:,:)::qT			!(qtMin,qtMax)
    real*8,intent(in),dimension(:,:)::Q				!(Qmin,Qmax)
    real*8,intent(in),dimension(:,:)::y				!(ymin,ymax)
    logical,intent(in),dimension(:)::includeCuts		!include cuts
    real*8,intent(in),dimension(:,:)::CutParameters	!(p1,p2,eta1,eta2)
    real*8,dimension(1:ListLength)::DY_xSec_List
    
    call xSec_DY_List(DY_xSec_List,process,s,qT,Q,y,includeCuts,CutParameters)
  
  end function DY_xSec_List

  function DY_xSec_List_APPROXIMATE(process,s,qT,Q,y,includeCuts,CutParameters,ListLength)
    integer,intent(in)::ListLength
    integer,intent(in),dimension(:,:)::process			!the number of process
    real*8,intent(in),dimension(:)::s				!Mandelshtam s
    real*8,intent(in),dimension(:,:)::qT			!(qtMin,qtMax)
    real*8,intent(in),dimension(:,:)::Q				!(Qmin,Qmax)
    real*8,intent(in),dimension(:,:)::y				!(ymin,ymax)
    logical,intent(in),dimension(:)::includeCuts		!include cuts
    real*8,intent(in),dimension(:,:)::CutParameters	!(p1,p2,eta1,eta2)
    real*8,dimension(1:ListLength)::DY_xSec_List_APPROXIMATE

    call xSec_DY_List_APPROXIMATE(DY_xSec_List_APPROXIMATE,process,s,qT,Q,y,includeCuts,CutParameters)

  end function DY_xSec_List_APPROXIMATE
  
  function DY_xSec_BINLESS_List(process,s,qT,Q,y,includeCuts,CutParameters,ListLength)
    integer,intent(in)::ListLength
    integer,intent(in),dimension(:,:)::process			!the number of process
    real*8,intent(in),dimension(:)::s				!Mandelshtam s
    real*8,intent(in),dimension(:)::qT			!(qtMin,qtMax)
    real*8,intent(in),dimension(:)::Q				!(Qmin,Qmax)
    real*8,intent(in),dimension(:)::y				!(ymin,ymax)
    logical,intent(in),dimension(:)::includeCuts		!include cuts
    real*8,intent(in),dimension(:,:)::CutParameters	!(p1,p2,eta1,eta2)
    real*8,dimension(1:ListLength)::DY_xSec_BINLESS_List

    call xSec_DY_List_BINLESS(DY_xSec_BINLESS_List,process,s,qT,Q,y,includeCuts,CutParameters)
  
  end function DY_xSec_BINLESS_List
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SIDIS CROSS-SECTION
  
  function SIDIS_xSec_Single(process,s,pT,z,x,Q,doCut,Cuts)
    integer,intent(in),dimension(1:4)::process			!the number of process
    real*8,intent(in)::s					!Mandelshtam s
    real*8,intent(in),dimension(1:2)::pT			!(qtMin,qtMax)
    real*8,intent(in),dimension(1:2)::z				!(zmin,zmax)
    real*8,intent(in),dimension(1:2)::x				!(xmin,xmax)
    real*8,intent(in),dimension(1:2)::Q				!(Qmin,Qmax)    
    logical,intent(in)::doCut					!triger cuts
    real*8,intent(in),dimension(1:4)::Cuts			!(ymin,yMax,W2min,W2max)
    real*8::SIDIS_xSec_Single
    
    call xSec_SIDIS(SIDIS_xSec_Single,process,s,pT,z,x,Q,doCut,Cuts)
  
  end function SIDIS_xSec_Single
  
  function SIDIS_xSec_Single_withMasses(process,s,pT,z,x,Q,doCut,Cuts,masses)
    integer,intent(in),dimension(1:4)::process			!the number of process
    real*8,intent(in)::s					!Mandelshtam s
    real*8,intent(in),dimension(1:2)::pT			!(qtMin,qtMax)
    real*8,intent(in),dimension(1:2)::z				!(zmin,zmax)
    real*8,intent(in),dimension(1:2)::x				!(xmin,xmax)
    real*8,intent(in),dimension(1:2)::Q				!(Qmin,Qmax)    
    logical,intent(in)::doCut					!triger cuts
    real*8,intent(in),dimension(1:4)::Cuts			!(ymin,yMax,W2min,W2max)
    real*8,intent(in),dimension(1:2)::masses			!(mTARGET,mPRODUCT)
    real*8::SIDIS_xSec_Single_withMasses
    
    call xSec_SIDIS(SIDIS_xSec_Single_withMasses,process,s,pT,z,x,Q,doCut,Cuts,masses)
  
  end function SIDIS_xSec_Single_withMasses
  
  function SIDIS_xSec_List(process,s,pT,z,x,Q,doCut,Cuts,masses,ListLength)
    integer,intent(in)::ListLength
    integer,intent(in),dimension(:,:)::process			!the number of process
    real*8,intent(in),dimension(:)::s				!Mandelshtam s
    real*8,intent(in),dimension(:,:)::pT			!(qtMin,qtMax)
    real*8,intent(in),dimension(:,:)::z				!(zmin,zmax)
    real*8,intent(in),dimension(:,:)::x				!(xmin,xmax)
    real*8,intent(in),dimension(:,:)::Q				!(Qmin,Qmax)        
    logical,intent(in),dimension(:)::doCut			!triger cuts
    real*8,intent(in),dimension(:,:)::Cuts			!(ymin,yMax,W2min,W2max)
    real*8,intent(in),dimension(:,:)::masses			!(mTARGET,mPRODUCT)
    real*8,dimension(1:ListLength)::SIDIS_xSec_List
    
    call xSec_SIDIS_List_forharpy(SIDIS_xSec_List,process,s,pT,z,x,Q,doCut,Cuts,masses)
  
  end function SIDIS_xSec_List
  
  function SIDIS_xSec_BINLESS_List(process,s,pT,z,x,Q,masses,ListLength)
    integer,intent(in)::ListLength
    integer,intent(in),dimension(:,:)::process			!the number of process
    real*8,intent(in),dimension(:)::s				!Mandelshtam s
    real*8,intent(in),dimension(:)::pT			        !(qt)
    real*8,intent(in),dimension(:)::z				!(z)
    real*8,intent(in),dimension(:)::x				!(x)
    real*8,intent(in),dimension(:)::Q				!(Q)        
    real*8,intent(in),dimension(:,:)::masses			!(mTARGET,mPRODUCT)
    real*8,dimension(1:ListLength)::SIDIS_xSec_BINLESS_List
    
    call xSec_SIDIS_BINLESS_List_forharpy(SIDIS_xSec_BINLESS_List,process,s,pT,z,x,Q,masses)
  
  end function SIDIS_xSec_BINLESS_List
  
  
end module harpy
