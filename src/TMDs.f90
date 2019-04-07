!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 2.0
!
!	Evaluation of the TMDs
!	
!	if you use this module please, quote 1706.01473
!
!	ver 1.0: release (AV, 10.05.2017)
!	ver 1.4: release (AV, 23.12.2018)
!	ver 1.41:release (AV, 23.12.2018)
!	ver 2.00:release (AV, 29.03.2019)
!
!				A.Vladimirov (23.12.2018)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDs
  use QCDinput
  use TMDR
  use uTMDPDF
  use uTMDFF
  implicit none

  private
!   public
 
 character (len=7),parameter :: moduleName="TMDs"
 character (len=5),parameter :: version="v2.00"
 
!------------------------------------------Physical and mathematical constants------------------------------------------
!------------------------------------------Working variables------------------------------------------------------------
  
  logical::started=.false.
  integer::outputLevel=2
  integer::messageTrigger=5
  
  !!!! The TMD evolution can be tritted differently
  integer::EvolutionType
  logical::include_uTMDPDF
  logical::include_uTMDFF
  
  
  !!!parameters for the uncertanty estimation
  real*8::c1_global,c3_global
  
!-----------------------------------------Public interface--------------------------------------------------------------
  public::TMDs_SetScaleVariations,TMDs_Initialize,TMDs_IsInitialized
  real*8,dimension(-5:5),public:: uTMDPDF_5,uTMDPDF_50
  real*8,dimension(-5:5),public:: uTMDFF_5,uTMDFF_50
  
  public::uPDF_uPDF,uPDF_anti_uPDF
  
  interface uTMDPDF_5
    module procedure uTMDPDF_5_Ev,uTMDPDF_5_optimal
  end interface
  
  interface uTMDPDF_50
    module procedure uTMDPDF_50_Ev,uTMDPDF_50_optimal
  end interface
  
  interface uTMDFF_5
    module procedure uTMDFF_5_Ev,uTMDFF_5_optimal
  end interface
  
  interface uTMDFF_50
    module procedure uTMDFF_50_Ev,uTMDFF_50_optimal
  end interface
  
  contains 
  
   INCLUDE 'Model/TMDs_model.f90'
   
  function TMDs_IsInitialized()
  logical::TMDs_IsInitialized
  TMDs_IsInitialized=started
  end function TMDs_IsInitialized
  
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
   
!!! Initializing routing
!!! Filles the prebuiled arrays
!!! orderAD, is order of anomalous dimension for evolution
!!!! order zeta is the order of pertrubative zeta expression, used only in the "universal TMD"
 subroutine TMDs_Initialize(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=300)::path,line
    logical::initRequared
    
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
    if(outputLevel>2) write(*,*) 'artemide.',moduleName,version,': initialization started ... '
    call MoveTO(51,'*p3  ')
    read(51,*) messageTrigger
    
    
    call MoveTO(51,'*6   ')
    call MoveTO(51,'*p1  ')
    read(51,*) initRequared
    if(.not.initRequared) then
      if(outputLevel>2) write(*,*)'artemide.',moduleName,': initialization is not requared. '
      started=.false.
      return
    end if
    CLOSE (51, STATUS='KEEP') 
    !!! then we read it again from the beginning to fill parameters of other modules
    OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")
    
    !! TMDR
    call MoveTO(51,'*3   ')
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p2  ')
    read(51,*) EvolutionType
    
    !! uTMDPDF
    call MoveTO(51,'*4   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_uTMDPDF
    
    !! uTMDFF
    call MoveTO(51,'*5   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_uTMDFF
    
    CLOSE (51, STATUS='KEEP')
    
    if(.not.QCDinput_IsInitialized()) then
      if(outputLevel>1) write(*,*) '.. initializing QCDinput (from ',moduleName,')'
      if(present(prefix)) then
	call QCDinput_Initialize(file,prefix)
      else
	call QCDinput_Initialize(file)
      end if
    end if
    
    if(.not.TMDR_IsInitialized()) then
      if(outputLevel>1) write(*,*) '.. initializing TMDR (from ',moduleName,')'
      if(present(prefix)) then
	call TMDR_Initialize(file,prefix)
      else
	call TMDR_Initialize(file)
      end if
    end if
    
    if(include_uTMDPDF .and. (.not.uTMDPDF_IsInitialized())) then
      if(outputLevel>1) write(*,*) '.. initializing uTMDPDF (from ',moduleName,')'
      if(present(prefix)) then
	call uTMDPDF_Initialize(file,prefix)
      else
	call uTMDPDF_Initialize(file)
      end if
    end if
    
    if(include_uTMDFF .and. (.not.uTMDFF_IsInitialized())) then
      if(outputLevel>1) write(*,*) '.. initializing uTMDFF (from ',moduleName,')'
      if(present(prefix)) then
	call uTMDFF_Initialize(file,prefix)
      else
	call uTMDFF_Initialize(file)
      end if
    end if
  
      c1_global=1d0
      c3_global=1d0
      
      
      started=.true.
     if(outputLevel>0) write(*,*) '----- arTeMiDe.TMDs ',version,': .... initialized'
     if(outputLevel>1) write(*,*) 
  end subroutine TMDs_Initialize

  !!!! this routine set the variations of scales
  !!!! it is used for the estimation of errors
  subroutine TMDs_SetScaleVariations(c1_in,c3_in)
    real*8::c1_in,c3_in
    
    If(EvolutionType==1) then
    if(c1_in<0.1d0 .or. c1_in>10.d0) then
    if(outputLevel>0) write(*,*) 'WARNING: arTeMiDe.TMDs: variation in c1 is enourmous. c1 is set to 2'
    c1_global=2d0
    else
    c1_global=c1_in
    end if
    else
      if(ABS(c1_in-1d0)>0.00001d0) then
        if(outputLevel>0) write(*,*) "WARNING: arTeMiDe:TMDs: variation of c1 is sensless. &
	    There is no such constant within current evolution scheme. Nothing is done"
      end if
    end if
    
    if(EvolutionType==3) then
      if(ABS(c3_in-1d0)>0.00001d0) then
        if(outputLevel>0) write(*,*) "WARNING: arTeMiDe:TMDs: variation of c3 is sensless. &
	    There is no such constant within current evolution scheme. Nothing is done"
      end if
    else
    if(c3_in<0.1d0 .or. c3_in>10.d0) then
    if(outputLevel>0) write(*,*) 'WARNING: arTeMiDe.TMDs: variation in c3 is enourmous. c3 is set to 2'
    c3_global=2d0
    else
    c3_global=c3_in
    end if
    end if
    
    if(outputLevel>1) write(*,*) 'TMDs: set scale variations constants (c1,c3) as:',c1_global,c3_global
  end subroutine TMDs_SetScaleVariations
  
  !!!!!!!!!!!!!!!!!!!uTMDPDF!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !!!!!!!! upolarized TMDPDF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function uTMDPDF_5_Ev(x,bt,muf,zetaf,hadron)
    real*8::uTMDPDF_5_Ev(-5:5)
    real*8:: x,bt,muf,zetaf
    real*8:: mui,Rkernel
    integer::hadron
    
   SELECT CASE(EvolutionType)
   CASE(1)!!!! improved D
       mui=mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,c3_global*mui,c1_global*mu0(bt),1)
   CASE(2)!!!! improved gamma
       mui=mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,c3_global*mui,1)
   CASE(3)!!!! fixed mu
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
   END SELECT
    uTMDPDF_5_Ev=Rkernel*uTMDPDF_lowScale5(x,bT,hadron)
    
    !!! forcefully set =0 below threshold
    if(muf<mBOTTOM) then
    uTMDPDF_5_Ev(5)=0d0
    uTMDPDF_5_Ev(-5)=0d0
    end if
    if(muf<mCHARM) then
    uTMDPDF_5_Ev(4)=0d0
    uTMDPDF_5_Ev(-4)=0d0
    end if

    
  end function uTMDPDF_5_Ev
  
  
      !!!!!!!! upolarized TMDPDF
  ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDPDF_50_Ev(x,bt,muf,zetaf,hadron)
    real*8::uTMDPDF_50_Ev(-5:5)
    real*8:: x,bt,muf,zetaf
    real*8:: mui,Rkernel ,RkernelG   
    integer::hadron
    
   SELECT CASE(EvolutionType)
   CASE(1)!!!! improved D
      mui=mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,c3_global*mui,c1_global*mu0(bt),1)
      RkernelG=TMDR_Rzeta(bt,muf,zetaf,c3_global*mui,c1_global*mu0(bt),0)
   CASE(2)!!!! improved gamma
       mui=mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,c3_global*mui,1)
      RkernelG=TMDR_Rzeta(bt,muf,zetaf,c3_global*mui,0)
   CASE(3)!!!! fixed mu
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
      RkernelG=TMDR_Rzeta(bt,muf,zetaf,0)
   END SELECT
    uTMDPDF_50_Ev=Rkernel*uTMDPDF_lowScale50(x,bT,hadron)
    uTMDPDF_50_Ev(0)=uTMDPDF_50_Ev(0)*RkernelG/Rkernel
    
    !!! forcefully set =0 below threshold
    if(muf<mBOTTOM) then
    uTMDPDF_50_Ev(5)=0d0
    uTMDPDF_50_Ev(-5)=0d0
    end if
    if(muf<mCHARM) then
    uTMDPDF_50_Ev(4)=0d0
    uTMDPDF_50_Ev(-4)=0d0
    end if
    
  end function uTMDPDF_50_Ev
  
    !!!!!!!! upolarized TMDPDF OPTIMAL
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function uTMDPDF_5_optimal(x,bt,hadron)
    real*8::uTMDPDF_5_optimal(-5:5)
    real*8:: x,bt
    integer::hadron
    
    uTMDPDF_5_optimal=uTMDPDF_lowScale5(x,bT,hadron)
  end function uTMDPDF_5_optimal
  
  
      !!!!!!!! upolarized TMDPDF
  ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDPDF_50_optimal(x,bt,hadron)
    real*8::uTMDPDF_50_optimal(-5:5)
    real*8:: x,bt
    integer::hadron
    
    uTMDPDF_50_optimal=uTMDPDF_lowScale50(x,bT,hadron)
    
  end function uTMDPDF_50_optimal
  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PRODUCTS FOR DY!!!!!!!!!!!!!!!!!!!!!!
    !!! Product of quark*antiquark uTMDPDFs. Slightly faster then just product
    ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function uPDF_uPDF(x1,x2,bt,muf,zetaf,hadron1,hadron2)
    real*8:: x1,x2,bt,muf,zetaf
    real*8:: mui,Rkernel
    integer::hadron1,hadron2
    real*8,dimension(-5:5)::tmd1,tmd2,uPDF_uPDF
    
   SELECT CASE(EvolutionType)
   CASE(1)!!!! improved D
       mui=mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,c3_global*mui,c1_global*mu0(bt),1)
   CASE(2)!!!! improved gamma
       mui=mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,c3_global*mui,1)
   CASE(3)!!!! fixed mu
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
   END SELECT
   tmd1=uTMDPDF_lowScale5(x1,bT,hadron1)
   tmd2=uTMDPDF_lowScale5(x2,bT,hadron2)
   uPDF_uPDF=(Rkernel**2)*tmd1*(tmd2(5:-5:-1))
    
    !!! forcefully set =0 below threshold
    if(muf<mBOTTOM) then
    uPDF_uPDF(5)=0d0
    uPDF_uPDF(-5)=0d0
    end if
    if(muf<mCHARM) then
    uPDF_uPDF(4)=0d0
    uPDF_uPDF(-4)=0d0
    end if

    
  end function uPDF_uPDF
  
    !!! Product of quark*quark uTMDPDFs. Slightly faster then just product
    ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function uPDF_anti_uPDF(x1,x2,bt,muf,zetaf,hadron1,hadron2)
    real*8:: x1,x2,bt,muf,zetaf
    real*8:: mui,Rkernel
    integer::hadron1,hadron2
    real*8,dimension(-5:5)::tmd1,tmd2,uPDF_anti_uPDF
    
   SELECT CASE(EvolutionType)
   CASE(1)!!!! improved D
       mui=mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,c3_global*mui,c1_global*mu0(bt),1)
   CASE(2)!!!! improved gamma
       mui=mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,c3_global*mui,1)
   CASE(3)!!!! fixed mu
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
   END SELECT
   tmd1=uTMDPDF_lowScale5(x1,bT,hadron1)
   tmd2=uTMDPDF_lowScale5(x2,bT,hadron2)
   uPDF_anti_uPDF=(Rkernel**2)*tmd1*tmd2
    
    !!! forcefully set =0 below threshold
    if(muf<mBOTTOM) then
    uPDF_anti_uPDF(5)=0d0
    uPDF_anti_uPDF(-5)=0d0
    end if
    if(muf<mCHARM) then
    uPDF_anti_uPDF(4)=0d0
    uPDF_anti_uPDF(-4)=0d0
    end if

    
  end function uPDF_anti_uPDF
  
  !!!!!!!!!!!!!!!!!!!uTMDFF!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
    !!!!!!!! upolarized TMDFF
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function uTMDFF_5_Ev(x,bt,muf,zetaf,hadron)
    real*8::uTMDFF_5_Ev(-5:5)
    real*8:: x,bt,muf,zetaf
    real*8:: mui,Rkernel
    integer::hadron
    
   SELECT CASE(EvolutionType)
   CASE(1)!!!! improved D
       mui=mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,c3_global*mui,c1_global*mu0(bt),1)
   CASE(2)!!!! improved gamma
       mui=mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,c3_global*mui,1)
   CASE(3)!!!! fixed mu
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
   END SELECT
    uTMDFF_5_Ev=Rkernel*uTMDFF_lowScale5(x,bT,hadron)
    
    !!! forcefully set =0 below threshold
    if(muf<mBOTTOM) then
    uTMDFF_5_Ev(5)=0d0
    uTMDFF_5_Ev(-5)=0d0
    end if
    if(muf<mCHARM) then
    uTMDFF_5_Ev(4)=0d0
    uTMDFF_5_Ev(-4)=0d0
    end if

    
  end function uTMDFF_5_Ev
  
  
      !!!!!!!! upolarized TMDFF
  ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDFF_50_Ev(x,bt,muf,zetaf,hadron)
    real*8::uTMDFF_50_Ev(-5:5)
    real*8:: x,bt,muf,zetaf
    real*8:: mui,Rkernel ,RkernelG   
    integer::hadron
    
   SELECT CASE(EvolutionType)
   CASE(1)!!!! improved D
      mui=mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,c3_global*mui,c1_global*mu0(bt),1)
      RkernelG=TMDR_Rzeta(bt,muf,zetaf,c3_global*mui,c1_global*mu0(bt),0)
   CASE(2)!!!! improved gamma
       mui=mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,c3_global*mui,1)
      RkernelG=TMDR_Rzeta(bt,muf,zetaf,c3_global*mui,0)
   CASE(3)!!!! fixed mu
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
      RkernelG=TMDR_Rzeta(bt,muf,zetaf,0)
   END SELECT
    uTMDFF_50_Ev=Rkernel*uTMDFF_lowScale50(x,bT,hadron)
    uTMDFF_50_Ev(0)=uTMDFF_50_Ev(0)*RkernelG/Rkernel
    
    !!! forcefully set =0 below threshold
    if(muf<mBOTTOM) then
    uTMDFF_50_Ev(5)=0d0
    uTMDFF_50_Ev(-5)=0d0
    end if
    if(muf<mCHARM) then
    uTMDFF_50_Ev(4)=0d0
    uTMDFF_50_Ev(-4)=0d0
    end if
    
  end function uTMDFF_50_Ev
  
  ! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
  function uTMDFF_5_optimal(x,bt,hadron)
    real*8::uTMDFF_5_optimal(-5:5)
    real*8:: x,bt
    integer::hadron
    
    uTMDFF_5_optimal=uTMDFF_lowScale5(x,bT,hadron)
    
  end function uTMDFF_5_optimal

    ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  function uTMDFF_50_optimal(x,bt,hadron)
    real*8::uTMDFF_50_optimal(-5:5)
    real*8:: x,bt
    integer::hadron
    
    uTMDFF_50_optimal=uTMDFF_lowScale50(x,bT,hadron)
    
  end function uTMDFF_50_optimal
  
end module TMDs