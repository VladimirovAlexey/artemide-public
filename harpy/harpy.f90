!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			harpy
!
!	An interface for artemide
!
!				A.Vladimirov (13.01.2019)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module harpy
use TMDs
use TMDs_inKT
use TMDX_DY
use TMDX_SIDIS

!!! this flag is requared to guaranty that artemide is not started twice (it lead to the crush)
logical::started=.false.

contains

  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!GENERAL
  subroutine Initialize(orderMain)
    character(len=*)::orderMain 
    if(started) then
      write(*,*) 'artemide already runs'
    else
      call TMDX_DY_Initialize(orderMain)
      call TMDX_SIDIS_Initialize(orderMain)
      call TMDs_inKT_Initialize(orderMain)
      started=.true.
    end if
  end subroutine Initialize
  
    !! call for parameters from the model
  subroutine SetLambda_ByReplica(num)
  integer:: num
  
  call TMDX_DY_SetNPParameters(num)
  
  end subroutine SetLambda_ByReplica
  
  !!!Sets the non-pertrubative parameters lambda
  subroutine SetLambda(lambdaIN)
    real*8,intent(in)::lambdaIN(:)
    
    call TMDX_DY_SetNPParameters(lambdaIN)
    
  end subroutine SetLambda
  
  !!!! this routine set the variations of scales
  !!!! it is used for the estimation of errors
  subroutine SetScaleVariation(c1_in,c2_in,c3_in,c4_in)
    real*8::c1_in,c2_in,c3_in,c4_in
    
    call TMDX_DY_SetScaleVariations(c1_in,c2_in,c3_in,c4_in)
    
  end subroutine SetScaleVariation
  
  !! reset the number for PDF replica for uTMDPDF
  subroutine SetPDFreplica(rep)
    integer::rep
    call TMDs_SetPDFreplica(rep)
  end subroutine SetPDFreplica
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!uTMD
  !!!!!!!! upolarized TMDPDF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function uTMDPDF_5_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: uTMDPDF_5_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  uTMDPDF_5_Evolved=uTMDPDF_5(x,bt,muf,zetaf,hadron)
    
  end function uTMDPDF_5_Evolved
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDPDF_50_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: uTMDPDF_50_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  uTMDPDF_50_Evolved=uTMDPDF_5(x,bt,muf,zetaf,hadron)
    
  end function uTMDPDF_50_Evolved

    !!!!!!!! upolarized TMDPDF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function uTMDPDF_5_Optimal(x,bt,hadron)
    real*8:: uTMDPDF_5_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  uTMDPDF_5_Optimal=uTMDPDF_5(x,bt,hadron)
    
  end function uTMDPDF_5_Optimal
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDPDF_50_Optimal(x,bt,hadron)
    real*8:: uTMDPDF_50_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  uTMDPDF_50_Optimal=uTMDPDF_5(x,bt,hadron)
    
  end function uTMDPDF_50_Optimal
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! UTMD IN KT
    !!!!!!!! upolarized TMDPDF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function uTMDPDF_kT_5_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: uTMDPDF_kT_5_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  uTMDPDF_kT_5_Evolved=uTMDPDF_kT_5(x,bt,muf,zetaf,hadron)
    
  end function uTMDPDF_kT_5_Evolved
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDPDF_kT_50_Evolved(x,bt,muf,zetaf,hadron)
    real*8:: uTMDPDF_kT_50_Evolved(-5:5)
    real*8:: x,bt,muf,zetaf
    integer::hadron
  
  uTMDPDF_kT_50_Evolved=uTMDPDF_kT_5(x,bt,muf,zetaf,hadron)
    
  end function uTMDPDF_kT_50_Evolved

    !!!!!!!! upolarized TMDPDF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function uTMDPDF_kT_5_Optimal(x,bt,hadron)
    real*8:: uTMDPDF_kT_5_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  uTMDPDF_kT_5_Optimal=uTMDPDF_kT_5(x,bt,hadron)
    
  end function uTMDPDF_kT_5_Optimal
  
    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDPDF_kT_50_Optimal(x,bt,hadron)
    real*8:: uTMDPDF_kT_50_Optimal(-5:5)
    real*8:: x,bt
    integer::hadron
  
  uTMDPDF_kT_50_Optimal=uTMDPDF_kT_5(x,bt,hadron)
    
  end function uTMDPDF_kT_50_Optimal
  
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!DY CROSS-SECTION
  
    function DY_xSec_SingleN(process,s,qT,Q,y,includeCuts,CutParameters,Num)
    integer,intent(in),dimension(1:3)::process		!the number of process
    real*8,intent(in)::s				!Mandelshtam s
    real*8,intent(in),dimension(1:2)::qT		!(qtMin,qtMax)
    real*8,intent(in),dimension(1:2)::Q			!(Qmin,Qmax)
    real*8,intent(in),dimension(1:2)::y			!(ymin,ymax)
    logical,intent(in)::includeCuts			!include cuts
    real*8,intent(in),dimension(1:4)::CutParameters	!(p1,p2,eta1,eta2)
    integer,intent(in)::Num				!number of sections
    real*8::DY_xSec_Single
    real*8::X
    
    call xSec_DY(X,process,s,qT,Q,y,includeCuts,CutParameters,Num)
    DY_xSec_SingleN=X
  
  end function DY_xSec_SingleN
  
  function DY_xSec_Single(process,s,qT,Q,y,includeCuts,CutParameters)
    integer,intent(in),dimension(1:3)::process		!the number of process
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
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!SIDIS CROSS-SECTION
  
  function SIDIS_xSec_Single(process,s,pT,z,x,Q,doCut,Cuts)
    integer,intent(in),dimension(1:3)::process			!the number of process
    real*8,intent(in)::s					!Mandelshtam s
    real*8,intent(in),dimension(1:2)::pT			!(qtMin,qtMax)
    real*8,intent(in),dimension(1:2)::z				!(zmin,zmax)
    real*8,intent(in),dimension(1:2)::x				!(xmin,xmax)
    real*8,intent(in),dimension(1:2)::Q				!(Qmin,Qmax)    
    logical,intent(in)::doCut					!triger cuts
    real*8,intent(in),dimension(1:3)::Cuts			!(ymin,yMax,W2)
    real*8::SIDIS_xSec_Single
    
    call xSec_SIDIS(SIDIS_xSec_Single,process,s,pT,z,x,Q,doCut,Cuts)
  
  end function SIDIS_xSec_Single
  
  function SIDIS_xSec_Single_withMasses(process,s,pT,z,x,Q,doCut,Cuts,masses)
    integer,intent(in),dimension(1:3)::process			!the number of process
    real*8,intent(in)::s					!Mandelshtam s
    real*8,intent(in),dimension(1:2)::pT			!(qtMin,qtMax)
    real*8,intent(in),dimension(1:2)::z				!(zmin,zmax)
    real*8,intent(in),dimension(1:2)::x				!(xmin,xmax)
    real*8,intent(in),dimension(1:2)::Q				!(Qmin,Qmax)    
    logical,intent(in)::doCut					!triger cuts
    real*8,intent(in),dimension(1:3)::Cuts			!(ymin,yMax,W2)
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
    real*8,intent(in),dimension(:,:)::Cuts			!(ymin,yMax,W2)
    real*8,intent(in),dimension(:,:)::masses			!(mTARGET,mPRODUCT)
    real*8,dimension(1:ListLength)::SIDIS_xSec_List
    
    call xSec_SIDIS_List_forharpy(SIDIS_xSec_List,process,s,pT,z,x,Q,doCut,Cuts,masses)
  
  end function SIDIS_xSec_List
  
  
end module harpy