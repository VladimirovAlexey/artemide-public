!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 1.4
!
! The module defining various QED and electro weak parameters
!
!    14.02.2019 values of ckm matrix are added AV.
!
!                        AV.  10.06.2018
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module EWinput
use aTMDe_Numerics
use IO_functions
implicit none

private

logical:: started=.false.
integer::outputLevel
character (len=7),parameter :: moduleName="EWinput"
character (len=5),parameter :: version="v3.00"
!Last appropriate verion of constants-file
integer,parameter::inputver=30

real(dp)::massZ,massW,massHIGGS,sW2,cW2
real(dp)::Vckm_UD,Vckm_US,Vckm_CD,Vckm_CS,Vckm_CB,Vckm_UB

public::alphaEM,EWinput_Initialize,EWinput_IsInitialized

!!-Z-gamma DY
real(dp),public::paramU,paramD,paramS,paramC,paramB,paramL,paramL_A
real(dp),public::paramMIXU,paramMIXD,paramMIXS,paramMIXC,paramMIXB,paramMIXL,paramMIXL_A

real(dp),public::zP_gg_U, zP_gg_D, zP_gg_S, zP_gg_C, zP_gg_B, zP_gg_L
real(dp),public::zP_gZ_U, zP_gZ_D, zP_gZ_S, zP_gZ_C, zP_gZ_B, zP_gZ_L
real(dp),public::zP_ZZ_U, zP_ZZ_D, zP_ZZ_S, zP_ZZ_C, zP_ZZ_B, zP_ZZ_L
real(dp),public::zM_gZ_U, zM_gZ_D, zM_gZ_S, zM_gZ_C, zM_gZ_B, zM_gZ_L
real(dp),public::zM_ZZ_U, zM_ZZ_D, zM_ZZ_S, zM_ZZ_C, zM_ZZ_B, zM_ZZ_L
real(dp),public::rP_ZZ_U, rP_ZZ_D, rP_ZZ_S, rP_ZZ_C, rP_ZZ_B, rP_ZZ_L
real(dp),public::rM_gZ_U, rM_gZ_D, rM_gZ_S, rM_gZ_C, rM_gZ_B, rM_gZ_L

!!-W DY
real(dp),public::paramW_UD,paramW_US,paramW_UB,paramW_CD,paramW_CS,paramW_CB,paramW_L

!!-EW-boson parameters
real(dp),public::GammaZ2,MZ2,MZ
real(dp),public::GammaW2,MW2,MW

!!-Higgs-boson parameters
real(dp),public::MH2,GammaH2,VEVH

!!-Lepton parameters
real(dp),public::massELECTRON,massMUON,massTAU
!!-quark masses parameters (used to compute threashold of running alpha)
real(dp)::massCHARM,massBOTTOM,massTOP

!!alphaQED parameters
real(dp)::alphaZ,alphaZinv,alphaTAUinv !!! these are from constans-file
real(dp)::alphaELECTRONinv,alphaMUONinv,alphaCHARMinv,alphaBOTTOMinv,alphaTOPinv!!! these are calculated in Set_betaQED()
!! running parameters
real(dp)::betaQED !!! this is corrected beta function of qed matched to mZ and mTau threasholds with 1-loop evolution

contains

 function EWinput_IsInitialized()
  logical::EWinput_IsInitialized
  
  EWinput_IsInitialized=started 
 end function EWinput_IsInitialized
 
 subroutine EWinput_Initialize(file,prefix)
  character(len=*)::file
  character(len=*),optional::prefix
  character(len=300)::path
  logical::initRequired
  real(dp)::dummy
  integer::FILEver
  
  if(started) return
  
  if(present(prefix)) then
    path=trim(adjustl(prefix))//trim(adjustr(file))
  else
    path=trim(adjustr(file))
  end if
  
  OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")
    !!! Search for output level
    call MoveTO(51,'*0   ')
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) FILEver
    if(FILEver<inputver) then
      write(*,*) 'artemide.'//trim(moduleName)//': const-file version is too old.'
      write(*,*) '             Update the const-file with artemide.setup'
      write(*,*) '  '
      stop
    end if
    call MoveTO(51,'*p2  ')
    read(51,*) outputLevel    
    if(outputLevel>1) write(*,*) '--------------------------------------------- '
    if(outputLevel>1) write(*,*) 'artemide.EWinput: initialization started ... '
  
    !!! check do we need initialisation?
    call MoveTO(51,'*2   ')
    call MoveTO(51,'*p1  ')
    read(51,*) initRequired
    if(.not.initRequired) then
      if(outputLevel>1) write(*,*)'artemide.EWinput: initialization is not required. '
      started=.false.
      return
    end if
    
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) alphaZinv
    alphaZ=1d0/alphaZinv
    
    call MoveTO(51,'*p2   ')
    read(51,*) sW2    !!!!!!!!!!sin^2 theta_W
    cw2=1d0-sw2        !!!!!!!!!!cos^2 theta_W
    
    call MoveTO(51,'*p3   ')!!!!CKM matrix
    read(51,*) Vckm_UD,Vckm_US,Vckm_UB
    read(51,*) Vckm_CD,Vckm_CS,Vckm_CB
    
    call MoveTO(51,'*p4  ')
    read(51,*) alphaTAUinv
    
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) massZ     !!!!!!!!!!Z mass
    MZ=massZ
    MZ2=massZ**2
    call MoveTO(51,'*p2  ')
    read(51,*) dummy
    GammaZ2=dummy**2     !!!!!!!!!!Z width
    
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) massW     !!!!!!!!!!W mass
    MW=massW
    MW2=massW**2
    call MoveTO(51,'*p2  ')
    read(51,*) dummy
    GammaW2=dummy**2     !!!!!!!!!!W width
    
    call MoveTO(51,'*D   ')
    call MoveTO(51,'*p1  ')
    read(51,*) massHIGGS     !!!!!!!!!!Higgs mass
    MH2=massHIGGS**2
    call MoveTO(51,'*p2  ')
    read(51,*) dummy
    GammaH2=dummy**2     !!!!!!!!!!Higgs width
    call MoveTO(51,'*p3  ')
    read(51,*) VEVH      !!!!!!!!!!Higgs VEV
    
    call MoveTO(51,'*E   ')
    call MoveTO(51,'*p1  ')
    read(51,*) massELECTRON     !!!!!!!!!!electron mass in GEV
    call MoveTO(51,'*p2  ')
    read(51,*) massMUON         !!!!!!!!!!muon mass in GEV
    call MoveTO(51,'*p3  ')
    read(51,*) massTAU         !!!!!!!!!!tau-lepton mass in GEV
    
    CLOSE (51, STATUS='KEEP')
    
    !!!!! read quark thresholds from the QCD section
    OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")
    call MoveTO(51,'*1   ')!!! QCD section
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) massCHARM
    call MoveTO(51,'*p2  ')
    read(51,*) massBOTTOM
    call MoveTO(51,'*p3  ')
    read(51,*) massTOP
    CLOSE (51, STATUS='KEEP')
  
    
  call Set_EWconstants()
  call Set_betaQED()
  
   if(outputLevel>0) write(*,*) color('----- arTeMiDe.EWinput '//trim(version)//': .... initialized',c_green)
   if(outputLevel>1) write(*,*) ' '
  started=.true.
 
 end subroutine EWinput_Initialize
 
 
   !!!!alpha EM (normalized at MZ as 127^{-1}
   !!!! with 1-loop run
  function alphaEM(mu)
  real(dp)::mu,alphaEM
  
    if(mu>massTOP) then
      alphaEM=1d0/(alphaTOPinv+16d0*betaQED*log(mu/massTOP))
    else if(mu>massBOTTOM) then
      alphaEM=1d0/(alphaZinv+40d0/3d0*betaQED*log(mu/massZ))
    else if(mu>massTAU) then
      alphaEM=1d0/(alphaTAUinv+38d0/3d0*betaQED*log(mu/massTAU))
    else if(mu>massCHARM) then
      alphaEM=1d0/(alphaTAUinv+32d0/3d0*betaQED*log(mu/massTAU))
    else if(mu>massMUON) then
      alphaEM=1d0/(alphaCHARMinv+8d0*betaQED*log(mu/massCHARM))
    else if(mu>massELECTRON) then
      alphaEM=1d0/(alphaMUONinv+16d0/3d0*betaQED*log(mu/massMUON))
    else
      alphaEM=1d0/alphaELECTRONinv
    end if
  
  end function alphaEM

 subroutine Set_EWconstants()
 real(dp)::ef,t3
 !!!! param is given by
 !!!! ((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2)
 !!!!  it is 2(gV^2+gA^2) for Z boson.
 
 !!!! paramMIX is given by
 !!!! eq(t2-2ef sW^2)/(2sw cW)
 !!!! eq*gV  for Z boson

 !!!! parameters z and r
 !!!! zP_gg = ef^2
 !!!! zP_gZ = eq(t2-2ef sW^2)/(2sw cW) = paramMIX
 !!!! zP_ZZ = ((1-2|eq|sw^2)^2+4eq^2sw^4)/(8sw^2cw^2) = param
 !!!! rP_ZZ = ef(ef sW^2-T3)/cW^2
 !!!! zM_gZ =-abs(ef)/4sWcW
 !!!! zM_ZZ=(4 abs(ef)sW^2-1)/8sw^2cw^2
 !!!! rM_gZ=abs(ef)/4sWcW
 
 !-------------------------------------------------------
 !---  Z-boson interaction
 !---------------U quark
 ef=2._dp/3._dp
 t3=+0.5_dp
 paramU=((1-2*Abs(ef)*sW2)**2+4*ef**2*sW2**2)/(8*sW2*cW2)
 paramMIXU=ef*(t3-2*ef*sW2)/(2*Sqrt(sw2*cw2))

 zP_gg_U=ef**2
 zP_gZ_U=ef*(t3-2*ef*sW2)/(2*Sqrt(sw2*cw2))
 zP_ZZ_U=((1-2*Abs(ef)*sW2)**2+4*ef**2*sW2**2)/(8*sW2*cW2)
 zM_gZ_U=-abs(ef)/4/sqrt(sW2*cW2)
 zM_ZZ_U=(4*abs(ef)*sW2-1)/(8*sW2*cW2)
 rP_ZZ_U=ef*(ef*sW2-t3)/cW2
 rM_gZ_U=abs(ef)/(4*sqrt(sW2*cW2))
 
 !---------------D-quark 
 ef=-1._dp/3._dp
 t3=-0.5_dp
 paramD=((1-2*Abs(ef)*sW2)**2+4*ef**2*sW2**2)/(8*sW2*cW2)
 paramMIXD=ef*(t3-2*ef*sW2)/(2*Sqrt(sw2*cw2))

 zP_gg_D=ef**2
 zP_gZ_D=ef*(t3-2*ef*sW2)/(2*Sqrt(sw2*cw2))
 zP_ZZ_D=((1-2*Abs(ef)*sW2)**2+4*ef**2*sW2**2)/(8*sW2*cW2)
 zM_gZ_D=-abs(ef)/4/sqrt(sW2*cW2)
 zM_ZZ_D=(4*abs(ef)*sW2-1)/(8*sW2*cW2)
 rP_ZZ_D=ef*(ef*sW2-t3)/cW2
 rM_gZ_D=abs(ef)/(4*sqrt(sW2*cW2))
 
 !---------------S-quark
 ef=-1._dp/3._dp
 t3=-0.5_dp
 paramS=((1-2*Abs(ef)*sW2)**2+4*ef**2*sW2**2)/(8*sW2*cW2)
 paramMIXS=ef*(t3-2*ef*sW2)/(2*Sqrt(sw2*cw2))
 
 zP_gg_S=ef**2
 zP_gZ_S=ef*(t3-2*ef*sW2)/(2*Sqrt(sw2*cw2))
 zP_ZZ_S=((1-2*Abs(ef)*sW2)**2+4*ef**2*sW2**2)/(8*sW2*cW2)
 zM_gZ_S=-abs(ef)/4/sqrt(sW2*cW2)
 zM_ZZ_S=(4*abs(ef)*sW2-1)/(8*sW2*cW2)
 rP_ZZ_S=ef*(ef*sW2-t3)/cW2
 rM_gZ_S=abs(ef)/(4*sqrt(sW2*cW2))

 !---------------C-quark
 ef=2._dp/3._dp
 t3=+0.5_dp
 paramC=((1-2*Abs(ef)*sW2)**2+4*ef**2*sW2**2)/(8*sW2*cW2)
 paramMIXC=ef*(t3-2*ef*sW2)/(2*Sqrt(sw2*cw2))

 zP_gg_C=ef**2
 zP_gZ_C=ef*(t3-2*ef*sW2)/(2*Sqrt(sw2*cw2))
 zP_ZZ_C=((1-2*Abs(ef)*sW2)**2+4*ef**2*sW2**2)/(8*sW2*cW2)
 zM_gZ_C=-abs(ef)/4/sqrt(sW2*cW2)
 zM_ZZ_C=(4*abs(ef)*sW2-1)/(8*sW2*cW2)
 rP_ZZ_C=ef*(ef*sW2-t3)/cW2
 rM_gZ_C=abs(ef)/(4*sqrt(sW2*cW2))

 !---------------B-quark
 ef=-1._dp/3._dp
 t3=-0.5_dp
 paramB=((1-2*Abs(ef)*sW2)**2+4*ef**2*sW2**2)/(8*sW2*cW2)
 paramMIXB=ef*(t3-2*ef*sW2)/(2*Sqrt(sw2*cw2))

 zP_gg_B=ef**2
 zP_gZ_B=ef*(t3-2*ef*sW2)/(2*Sqrt(sw2*cw2))
 zP_ZZ_B=((1-2*Abs(ef)*sW2)**2+4*ef**2*sW2**2)/(8*sW2*cW2)
 zM_gZ_B=-abs(ef)/4/sqrt(sW2*cW2)
 zM_ZZ_B=(4*abs(ef)*sW2-1)/(8*sW2*cW2)
 rP_ZZ_B=ef*(ef*sW2-t3)/cW2
 rM_gZ_B=abs(ef)/(4*sqrt(sW2*cW2))

  !---------------Lepton
 ef=-1._dp
 t3=-0.5_dp
 paramL=((1-2*Abs(ef)*sW2)**2+4*ef**2*sW2**2)/(8*sW2*cW2)
 paramMIXL=ef*(t3-2*ef*sW2)/(2*Sqrt(sw2*cw2))
 !!! asymetric combinations
 paramL_A=(4*Abs(ef)*sW2**2-1)/(8*sW2*cW2)
 paramMIXL_A=-t3*ef/(2*Sqrt(sw2*cw2))

 zP_gg_L=ef**2
 zP_gZ_L=ef*(t3-2*ef*sW2)/(2*Sqrt(sw2*cw2))
 zP_ZZ_L=((1-2*Abs(ef)*sW2)**2+4*ef**2*sW2**2)/(8*sW2*cW2)
 zM_gZ_L=-abs(ef)/4/sqrt(sW2*cW2)
 zM_ZZ_L=(4*abs(ef)*sW2-1)/(8*sW2*cW2)
 rP_ZZ_L=ef*(ef*sW2-t3)/cW2
 rM_gZ_L=abs(ef)/(4*sqrt(sW2*cW2))

 !-------------------------------------------------------
 !---  W-boson interaction
 
 paramW_UD=abs(Vckm_UD)**2/(4d0*sW2)
 paramW_US=abs(Vckm_US)**2/(4d0*sW2)
 paramW_UB=abs(Vckm_UB)**2/(4d0*sW2)
 paramW_CD=abs(Vckm_CD)**2/(4d0*sW2)
 paramW_CS=abs(Vckm_CS)**2/(4d0*sW2)
 paramW_CB=abs(Vckm_CB)**2/(4d0*sW2)
 
 paramW_L=1d0/(4d0*sW2)
 
 end subroutine Set_EWconstants

 !!! compute the matched value of QED beta function
 !!! It uses alphaZ and alphaTau to determine the beta function of QED with threasholds
 !!! the threasholds used
 !!! e(+u+d); muon(+s); c; t; b; t;
 !!! betaEFF=Neff*beta
 !!! beta is almost -1/3pi
 !!! Neff=N_leptons+sum_q e_q^2 N_c
 !!! Neff= 8/3; 4; 16/3;, 19/3; 20/3; 8
 !!!
 !!! So fixing at z and tau gives beta=(alpha^{-1}(tau)-alpha^{-1}(Z))/2/(Neff(top<->b) Log[b/MZ]+Neff(b<->tau) Log[tau/b])
 !!! this number should be close to -1/3pi
 !!! 
 !!! It also compute threashold values of alpha QED
 subroutine Set_betaQED()
  real(dp),parameter::betaQED_1loop=-0.1061032953945969d0  
  real(dp)::deltaB
  betaQED=(alphaTAUinv-alphaZinv)/2d0/(20d0/3d0*Log(massBOTTOM/massZ)+19d0/3d0*Log(massTAU/massBOTTOM))

  deltaB=betaQED/betaQED_1loop
  
  if(outputLevel>2) write(*,*) "    Effective QED beta fuction in fractions of LO:",deltaB
  if(abs(deltaB-1d0)>0.05) &
    write(*,*)  WarningString(' Effective QED beta function 5% deviate from LO. Check boundary setup.',modulename)
  
  alphaTOPinv=alphaZinv+2d0*betaQED*20d0/3d0*log(massTOP/massZ)
  alphaBOTTOMinv=alphaZinv+2d0*betaQED*20d0/3d0*log(massBOTTOM/massZ)
  !alphaTAUinv  !! exact
  alphaCHARMinv=alphaTAUinv+2d0*betaQED*16d0/3d0*log(massCHARM/massTAU)
  alphaMUONinv=alphaCHARMinv+2d0*betaQED*4d0*log(massMUON/massCHARM)
  alphaELECTRONinv=alphaMUONinv+2d0*betaQED*8d0/3d0*log(massELECTRON/massMUON)
  
  if(outputLevel>2) write(*,*) "    Theashold values of alpha^(-1) QED (mE,mMU,mC,mTAU,mB,mT):"
  if(outputLevel>2) write(*,"('    ',F7.3,',',F7.3,',',F7.3,',',F7.3,',',F7.3,',',F7.3)") &
          alphaELECTRONinv,alphaMUONinv,alphaCHARMinv,alphaTAUinv,alphaBOTTOMinv,alphaTOPinv

 end subroutine Set_betaQED
 
end module EWinput


 
