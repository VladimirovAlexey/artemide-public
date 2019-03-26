!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 1.2
!
! Interface module to the user defined alpha-s, PDF, FF, etc.
! Could be interfaced to LHAPDF
!
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module QCDinput

implicit none

private

logical:: started

public::QCDinput_Initialize,As,QCDinput_IsInitialized
public::xPDF,xFF
public:: QCDinput_SetPDFreplica

 contains 
 
 
 function QCDinput_IsInitialized()
  logical::QCDinput_IsInitialized
  QCDinput_IsInitialized=started 
 end function QCDinput_IsInitialized
 
 !------------------------------Functions below are to be changed by user (if needed)
 
 !!Initialization with order
 !! order 	='LO'
 !!		='LO+'
 !!		='NLO'
 !!		='NLO+'
 !!		='NNLO'
 !!		='NNLO+'
 !! see details in manual
 subroutine QCDinput_Initialize(order,uPDF,uFF)
  character*64:: name,namePIp,namePIm,namePI,nameKAp,nameKAm,nameKA
  character(len=*)::order
  character(len=*),optional::uPDF,uFF
  character*64:: uPDF_order,uFF_order
  
  if(started) then
    return
  else
    
    !!! wrap the input
    if(present(uPDF)) then
      uPDF_order=uPDF
    else
      uPDF_order=order
    end if
    
    if(present(uFF)) then
      uFF_order=uFF
    else
      uFF_order=order
    end if
    
    !---------------------------------------------------------------------------!
    !--------------- write the initialization code here ------------------------!
    !---------------------------------------------------------------------------!
    
    !-----------------------uPDF set------------------------------------
    if(uPDF_order == "NONE") then
      write(*,*) 'uPDF initialisation skiped.'
    else
    SELECT CASE(uPDF_order)
     CASE ("LO")
! 	name='MMHT2014lo68cl'
	name='NNPDF31_lo_as_0118'
      CASE ("LO+")
! 	name='MMHT2014lo68cl'
	name='NNPDF31_lo_as_0118'
      CASE ("NLO")
! 	name='MMHT2014nlo68cl'
	name='NNPDF31_nlo_as_0118'
      CASE ("NLO+")
! 	name='MMHT2014nlo68cl'
        name='NNPDF31_nlo_as_0118'
      CASE ("NNLO")
! 	name='MMHT2014nnlo68cl'
	name='NNPDF31_nnlo_as_0118'!_1000'
      CASE ("NNLO+")
! 	name='MMHT2014nnlo68cl'
	name='NNPDF31_nnlo_as_0118'
      CASE DEFAULT
	name='MMHT2014nlo68cl'
     END SELECT
    
    call InitPDFsetByNameM(1,name)
    call InitPDFM(1,0)!central
    end if
    
    
    !-----------------------uFF set------------------------------------    
!     if(uFF_order == "NONE") then
!       write(*,*) 'uFF initialisation skiped.'
!     else
!     SELECT CASE(uFF_order)
!      CASE ("LO","LO+","NLO","NLO+","NNLO","NNLO+")
!      !!! NNFF
! ! 	namePIp='NNFF10_PIp_nnlo'
! ! 	namePIm='NNFF10_PIm_nnlo'
! ! 	namePI='NNFF10_PIsum_nnlo'
! ! 	nameKAp='NNFF10_KAp_nnlo'
! ! 	nameKAm='NNFF10_KAm_nnlo'
! ! 	nameKA='NNFF10_KAsum_nnlo'
!     !!! DSS NLO
!      	namePIp='dsspipNLO'
!  	namePIm='dsspimNLO'
!  	namePI='dsspipNLO'
!  	nameKAp='dssKpNLO'
!  	nameKAm='dssKmNLO'
!  	nameKA='dssKpNLO'
!       CASE DEFAULT
! 	name='MMHT2014nlo68cl'
!      END SELECT
!     
!     !pi+
!     call InitPDFsetByNameM(2,namePIp)
!     call InitPDFM(2,0)!central
!     !pi-
!     call InitPDFsetByNameM(3,namePIm)
!     call InitPDFM(3,0)!central
!     !pi+ + pi-
!     call InitPDFsetByNameM(4,namePIp)
!     call InitPDFM(4,0)!central
!     !K+
!     call InitPDFsetByNameM(5,nameKAp)
!     call InitPDFM(5,0)!central
!     !K-
!     call InitPDFsetByNameM(6,nameKAm)
!     call InitPDFM(6,0)!central
!     !K+ + K-
!     call InitPDFsetByNameM(7,nameKA)
!     call InitPDFM(7,0)!central
!     
!     end if
!     
!     
   end if
  
  started=.true.
 
 end subroutine QCDinput_Initialize
 
 !!! set a different replica number for PDF.
 subroutine QCDinput_SetPDFreplica(rep)
 integer:: rep
  call InitPDF(rep)
 end subroutine QCDinput_SetPDFreplica
 
 !!!!alphas(Q)/4pi
 !!! NOT FORGET 4 PI !!!
 function As(Q)
 real*8::as,Q,alphasPDF
 real*8,parameter::Pi4=12.566370614359172d0
 
 As=alphasPDF(Q)/Pi4
 
 end function As

 !!!!array of x times PDF(x,Q) for hadron 'hadron'
 !!!! array is (-5:5) (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
 function xPDF(x,Q,hadron)
      real*8 :: x,Q
      integer:: hadron
      real*8, dimension(-5:5):: xPDF
      real*8, dimension (-6:6)::inputPDF
      
      call evolvePDFM(1,x,Q,inputPDF)
      
      xPDF=inputPDF(-5:5)
      
  end function xPDF
  
  
    !!!! return x*F(x,mu)
  !!!! enumeration of flavors
  !!!!  f = -5,-4, -3,  -2,  -1,0,1,2,3,5,4
  !!!!    = bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b
  !!!! enumeration of hadrons 
  !!!! 1=pion_+, 2=pion_0, 3=pion_-
  !!!! 4=kaon_+, 5=kaon_0, 6=kaon_-
  function xFF(x,Q,hadron)
      integer :: hadron
      real*8 :: x,Q
      real*8,dimension(-5:5):: xFF
      real*8, dimension (-6:6)::inputFF
      
      call evolvePDFM(hadron+1,x,Q,inputFF)
      
      xFF=inputFF(-5:5)
!       real*8::U, UB, D, DB, S, SB,C,B,GL
!       SELECT CASE(hadron)
! 	CASE (1)!pi+
! 	  call fDSS(2,1,1, 1, x, Q**2, U, UB, D, DB, S, SB,C,B,GL)
! 	CASE (2)!pi-
! 	  call fDSS(2,1,-1, 1, x, Q**2, U, UB, D, DB, S, SB,C,B,GL)
! 	CASE (3)!pi0
! 	  call fDSS(2,1,0,1, x, Q**2, U, UB, D, DB, S, SB,C,B,GL)
! 	CASE (4)!K+
! 	  call fDSS(2,2,1, 1, x, Q**2, U, UB, D, DB, S, SB,C,B,GL)
! 	CASE (5)!K-
! 	  call fDSS(2,2,-1, 1, x, Q**2, U, UB, D, DB, S, SB,C,B,GL)
! 	CASE (6)!K0
! 	  call fDSS(2,2,0,1, x, Q**2, U, UB, D, DB, S, SB,C,B,GL)
! 	CASE DEFAULT
! 	  write(*,*) 'WARNING: arTeMiDe.uTFF_xFF: unknown hadron (=',hadron,'). Evaluation stop'
! 	  stop
! 	END SELECT
! 	
!      xFF=(/B,C,SB,UB,DB,GL,D,U,S,C,B/)
  end function xFF
 
end module QCDinput
