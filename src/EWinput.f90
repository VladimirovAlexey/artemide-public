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

logical:: started=.false.
integer::outputLevel

real*8::Zmass,Wmass,alphaZ,sW2,cW2
real*8::Vckm_UD,Vckm_US,Vckm_CD,Vckm_CS,Vckm_CB,Vckm_UB

public::alphaEM,EWinput_Initialize,EWinput_IsInitialized

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
 
  !!! move CURRET in streem to the next line that starts from pos (5 char)
 subroutine MoveTO(streem,pos)
 integer,intent(in)::streem
 character(len=5)::pos
 character(len=300)::line
    do
    read(streem,'(A)') line    
    if(line(1:5)==pos) exit
    end do
 end subroutine MoveTO
 
 subroutine EWinput_Initialize(file,prefix)
  character(len=*)::file
  character(len=*),optional::prefix
  character(len=300)::path
  logical::initRequared
  real*8::dummy
  
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
    call MoveTO(51,'*p2  ')
    read(51,*) outputLevel    
    if(outputLevel>2) write(*,*) '--------------------------------------------- '
    if(outputLevel>2) write(*,*) 'artemide.EWinput: initialization started ... '
  
    !!! check do we need initialisation?
    call MoveTO(51,'*2   ')
    call MoveTO(51,'*p1  ')
    read(51,*) initRequared
    if(.not.initRequared) then
      if(outputLevel>2) write(*,*)'artemide.EWinput: initialization is not requared. '
      started=.false.
      return
    end if
    
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) dummy
    alphaZ=1d0/dummy
    
    call MoveTO(51,'*p2   ')
    read(51,*) sW2	!!!!!!!!!!sin^2 theta_W
    cw2=1d0-sw2		!!!!!!!!!!cos^2 theta_W
    
    call MoveTO(51,'*p3   ')
    read(51,*) Vckm_UD,Vckm_US,Vckm_UB
    read(51,*) Vckm_CD,Vckm_CS,Vckm_CB
    
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) Zmass     !!!!!!!!!!Z mass
    MZ2=Zmass**2
    call MoveTO(51,'*p2  ')
    read(51,*) dummy
    GammaZ2=dummy**2     !!!!!!!!!!Z width
    
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) Wmass     !!!!!!!!!!W mass
    MW2=Wmass**2
    call MoveTO(51,'*p2  ')
    read(51,*) dummy
    GammaW2=dummy**2     !!!!!!!!!!W width
    
    
    
    CLOSE (51, STATUS='KEEP')
  
    
  call Set_EWconstants()
  if(outputLevel>2)	write(*,*)'EWinput succesfully initialized.'
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
 
 paramW_UD=abs(Vckm_UD)**2/(4d0*sW2)
 paramW_US=abs(Vckm_US)**2/(4d0*sW2)
 paramW_UB=abs(Vckm_UB)**2/(4d0*sW2)
 paramW_CD=abs(Vckm_CD)**2/(4d0*sW2)
 paramW_CS=abs(Vckm_CS)**2/(4d0*sW2)
 paramW_CB=abs(Vckm_CB)**2/(4d0*sW2)
 
 paramW_L=1d0/(4d0*sW2)
 
 end subroutine Set_EWconstants

end module EWinput


 