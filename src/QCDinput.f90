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
 subroutine QCDinput_Initialize(order)
  character*64:: name
  character(len=*)::order
  
  if(started) then
    return
  else
    
    !---------------------------------------------------------------------------!
    !--------------- write the initialization code here ------------------------!
    !---------------------------------------------------------------------------!
    
    SELECT CASE(order)
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
	name='NNPDF31_nnlo_as_0118'
      CASE ("NNLO+")
! 	name='MMHT2014nnlo68cl'
	name='NNPDF31_nnlo_as_0118'
      CASE DEFAULT
	name='MMHT2014nlo68cl'
     END SELECT
    
    call InitPDFsetByName(name)
    call InitPDF(0)!central
    
  end if
  
  started=.true.
 
 end subroutine QCDinput_Initialize
 
 
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
      
      call evolvePDF(x,Q,inputPDF)
      
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
      real*8::U, UB, D, DB, S, SB,C,B,GL
      SELECT CASE(hadron)
	CASE (1)!pi+
	  call fDSS(2,1,1, 1, x, Q**2, U, UB, D, DB, S, SB,C,B,GL)
	CASE (2)!pi0
	  call fDSS(2,1,0, 1, x, Q**2, U, UB, D, DB, S, SB,C,B,GL)
	CASE (3)!pi-
	  call fDSS(2,1,-1,1, x, Q**2, U, UB, D, DB, S, SB,C,B,GL)
	CASE (4)!K+
	  call fDSS(2,2,1, 1, x, Q**2, U, UB, D, DB, S, SB,C,B,GL)
	CASE (5)!K0
	  call fDSS(2,2,0, 1, x, Q**2, U, UB, D, DB, S, SB,C,B,GL)
	CASE (6)!K-
	  call fDSS(2,2,-1,1, x, Q**2, U, UB, D, DB, S, SB,C,B,GL)
	CASE DEFAULT
	  write(*,*) 'WARNING: arTeMiDe.uTFF_xFF: unknown hadron (=',hadron,'). Evaluation stop'
	  stop
	END SELECT
	
     xFF=(/B,C,SB,UB,DB,GL,D,U,S,C,B/)
  end function xFF
 
end module QCDinput