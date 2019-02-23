!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 1.4
!
! The module defining various QED and electro weak parameters
!	
!	14.02.2019 values of ckm matrix are added AV.
!
!						AV.  10.06.2018
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module EWinput

implicit none

private

logical:: started

real*8::Zmass,Wmass,alphaZ,sW2,cW2
real*8::Vckm_UD,Vckm_US,Vckm_CD,Vckm_CS,Vckm_CB,Vckm_UB

public::alphaEM,EWinput_Initialize

!!-Z-gamma DY
real*8,public::paramU,paramD,paramS,paramC,paramB,paramL
real*8,public::paramMIXU,paramMIXD,paramMIXS,paramMIXC,paramMIXB,paramMIXL
!!-W DY
real*8,public::paramW_UD,paramW_US,paramW_UB,paramW_CD,paramW_CS,paramW_CB,paramW_L

!!-EW-boson parameters
real*8,public::GammaZ2,MZ2
real*8,public::GammaW2,MW2


contains

 function EWinput_IsInitialized()
  logical::EWinput_IsInitialized
  
  EWinput_IsInitialized=started 
 end function EWinput_IsInitialized
 
 
 subroutine EWinput_Initialize(order)
  character(len=*)::order
  character(256)::line
  real*8::dummy
  
  if(started)  return
  
    !---------------------------------------------------------------------------!
    !--------------- Eat the numerical values from constants--------------------!
    !---------------------------------------------------------------------------!
  OPEN(UNIT=51, FILE='constants', ACTION="read", STATUS="old")    
    !!!! Physical constants
    do
    read(51,'(A)') line    
    if(line(1:3)=='*1 ') exit
    end do    
    do
    read(51,'(A)') line
    if(line(1:3)=='*B ') exit
    end do
    
    do
    read(51,'(A)') line
    if(line(1:3)=='*1)') exit
    end do
    read(51,'(A)') line
    read(51,*) dummy	!!!!!!!!!!sin^2 theta_W
    sW2=dummy
    cw2=1d0-sw2		!!!!!!!!!!cos^2 theta_W
    
    do
    read(51,'(A)') line
    if(line(1:3)=='*2)') exit
    end do
    read(51,'(A)') line
    read(51,*) Zmass     !!!!!!!!!!Z mass
    MZ2=Zmass**2
    read(51,'(A)') line
    read(51,*) dummy
    GammaZ2=dummy**2     !!!!!!!!!!Z width
    
    do
    read(51,'(A)') line
    if(line(1:3)=='*3)') exit
    end do
    read(51,'(A)') line
    read(51,*) Wmass     !!!!!!!!!!W mass
    MW2=Wmass**2
    read(51,'(A)') line
    read(51,*) dummy
    GammaW2=dummy**2     !!!!!!!!!!W width
    
    
    do
    read(51,'(A)') line
    if(line(1:3)=='*D ') exit
    end do
    read(51,'(A)') line
    read(51,*) dummy
    alphaZ=1d0/dummy
    
    
    CLOSE (51, STATUS='KEEP')
  
    !---------------------------------------------------------------------------!
    !--------------- write the initialization code here ------------------------!
    !---------------------------------------------------------------------------!
    
    SELECT CASE(order)
     CASE ("LO")
      CASE ("LO+")
      CASE ("NLO")
      CASE ("NLO+")
      CASE ("NNLO")
      CASE ("NNLO+")
      CASE DEFAULT
     END SELECT
    
  call Set_EWconstants()
  
  started=.true.
 
 end subroutine EWinput_Initialize
 
 
   !!!!alpha EM (normalized at MZ as 127^{-1}
   !!!! with 1-loop run
  function alphaEM(mu)
  real*8::mu,alphaEM
  real*8,parameter::beta0=-0.1061032953945969d0*(2d0+3d0*11d0/9d0) !!!! =-1/3pi  * (NUMBER OF LEPTONS+Nc*sum(e^_q))
  alphaEM=alphaZ/(1+alphaZ*beta0*(2d0*LOG(mu/Zmass)-5d0/3d0))
  end function alphaEM

 subroutine Set_EWconstants()
 real*8::ef,t3
 !!!! param is given by
 !!!! ((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2)
 !!!!  it is 2(gV^2+gA^2) for Z boson.
 
 !!!! paramMIX is given by
 !!!! eq(t2-2ef sW^2)/(2sw cW)
 !!!! eq*gV  for Z boson
 
 !-------------------------------------------------------
 !---  Z-boson interaction
 !---------------U quark
 ef=2d0/3d0
 t3=+0.5d0
 paramU=((1d0-2d0*Abs(ef)*sW2)**2+4d0*ef**2*sW2**2)/(8d0*sW2*cW2)
 paramMIXU=ef*(t3-2d0*ef*sW2)/(2d0*Sqrt(sw2*cw2))
 
 !---------------D-quark 
 ef=-1d0/3d0
 t3=-0.5d0
 paramD=((1d0-2d0*Abs(ef)*sW2)**2+4d0*ef**2*sW2**2)/(8d0*sW2*cW2)
 paramMIXD=ef*(t3-2d0*ef*sW2)/(2d0*Sqrt(sw2*cw2))
 
 !---------------S-quark
 ef=-1d0/3d0
 t3=-0.5d0
 paramS=((1d0-2d0*Abs(ef)*sW2)**2+4d0*ef**2*sW2**2)/(8d0*sW2*cW2)
 paramMIXC=ef*(t3-2d0*ef*sW2)/(2d0*Sqrt(sw2*cw2))
 
 !---------------C-quark
 ef=2d0/3d0
 t3=+0.5d0
 paramC=((1d0-2d0*Abs(ef)*sW2)**2+4d0*ef**2*sW2**2)/(8d0*sW2*cW2)
 paramMIXC=ef*(t3-2d0*ef*sW2)/(2d0*Sqrt(sw2*cw2))
 
 !---------------B-quark
 ef=-1d0/3d0
 t3=-0.5d0
 paramB=((1d0-2d0*Abs(ef)*sW2)**2+4d0*ef**2*sW2**2)/(8d0*sW2*cW2)
 paramMIXB=ef*(t3-2d0*ef*sW2)/(2d0*Sqrt(sw2*cw2))
 
  !---------------Lepton
 ef=-1d0
 t3=-0.5d0
 paramL=((1d0-2d0*Abs(ef)*sW2)**2+4d0*ef**2*sW2**2)/(8d0*sW2*cW2)
 paramMIXL=ef*(t3-2d0*ef*sW2)/(2d0*Sqrt(sw2*cw2))
 
 !-------------------------------------------------------
 !---  W-boson interaction
 !----elements of CKM matrix from PDG2018
 Vckm_UD=0.97420d0
 Vckm_US=0.2243d0
 Vckm_CD=0.218d0
 Vckm_CS=0.997d0
 Vckm_CB=0.0422d0
 Vckm_UB=0.0394d0
 
 paramW_UD=abs(Vckm_UD)**2/(4d0*sW2)
 paramW_US=abs(Vckm_US)**2/(4d0*sW2)
 paramW_UB=abs(Vckm_UB)**2/(4d0*sW2)
 paramW_CD=abs(Vckm_CD)**2/(4d0*sW2)
 paramW_CS=abs(Vckm_CS)**2/(4d0*sW2)
 paramW_CB=abs(Vckm_CB)**2/(4d0*sW2)
 
 paramW_L=1d0/(4d0*sW2)
 
 end subroutine Set_EWconstants

end module EWinput


 