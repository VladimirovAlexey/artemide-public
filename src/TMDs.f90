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
!	ver 2.01:release (AV, 08.06.2019)
!
!				A.Vladimirov (23.12.2018)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDs
  use IO_functions
  use QCDinput
  use TMDR
  use uTMDPDF
  use uTMDFF
  use lpTMDPDF
  implicit none

  private
!   public
 
 character (len=7),parameter :: moduleName="TMDs"
 character (len=5),parameter :: version="v2.01"
 !Last appropriate verion of constants-file
  integer,parameter::inputver=3
 
!------------------------------------------Physical and mathematical constants------------------------------------------
!------------------------------------------Working variables------------------------------------------------------------
  
  logical::started=.false.
  integer::outputLevel=2
  integer::messageTrigger=5
  
  !!!! The TMD evolution can be tritted differently
  integer::EvolutionType
  logical::include_uTMDPDF
  logical::include_uTMDFF
  logical::include_lpTMDPDF
  
  
  !!!parameters for the uncertanty estimation
  real*8::c1_global,c3_global
  
!-----------------------------------------Public interface--------------------------------------------------------------
  public::TMDs_SetScaleVariations,TMDs_Initialize,TMDs_IsInitialized
  real*8,dimension(-5:5),public:: uTMDPDF_5,uTMDPDF_50
  real*8,dimension(-5:5),public:: uTMDFF_5,uTMDFF_50
  real*8,dimension(-5:5),public:: lpTMDPDF_50
  
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
  
  interface lpTMDPDF_50
    module procedure lpTMDPDF_50_Ev,lpTMDPDF_50_optimal
  end interface
  
  contains 
  
   INCLUDE 'Model/TMDs_model.f90'
   
  function TMDs_IsInitialized()
  logical::TMDs_IsInitialized
  TMDs_IsInitialized=started
  end function TMDs_IsInitialized

!!! Initializing routing
!!! Filles the prebuiled arrays
!!! orderAD, is order of anomalous dimension for evolution
!!!! order zeta is the order of pertrubative zeta expression, used only in the "universal TMD"
 subroutine TMDs_Initialize(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=300)::path,line
    logical::initRequared
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
      write(*,*) '		     Update the const-file with artemide.setup'
      write(*,*) '  '
      stop
    end if
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
    if(outputLevel>2) write(*,'(A,I3)') ' Evolution type =',EvolutionType
    
    !! uTMDPDF
    call MoveTO(51,'*4   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_uTMDPDF
    
    !! uTMDFF
    call MoveTO(51,'*5   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_uTMDFF
    
    !! uTMDFF
    call MoveTO(51,'*11   ')
    call MoveTO(51,'*p1  ')
    read(51,*) include_lpTMDPDF
    
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
    
     if(include_lpTMDPDF .and. (.not.lpTMDPDF_IsInitialized())) then
      if(outputLevel>1) write(*,*) '.. initializing lpTMDPDF (from ',moduleName,')'
      if(present(prefix)) then
	call lpTMDPDF_Initialize(file,prefix)
      else
	call lpTMDPDF_Initialize(file)
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
      mui=c3_global*mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,c1_global*mu0(bt),1)
   CASE(2)!!!! improved gamma
       mui=c3_global*mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
   CASE(3)!!!! fixed mu
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
   CASE(4)!!!! exact solution via zeta-line
      mui=c3_global*mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
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
      mui=c3_global*mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,c1_global*mu0(bt),1)
      RkernelG=TMDR_Rzeta(bt,muf,zetaf,mui,c1_global*mu0(bt),0)
   CASE(2)!!!! improved gamma
       mui=c3_global*mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
      RkernelG=TMDR_Rzeta(bt,muf,zetaf,mui,0)
   CASE(3)!!!! fixed mu
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
      RkernelG=TMDR_Rzeta(bt,muf,zetaf,0)
   CASE(4)!!!! exact solution via zeta-line
      mui=c3_global*mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bT,muf,zetaf,mui,1)
      RkernelG=TMDR_Rzeta(bt,muf,zetaf,mui,0)
   END SELECT
    !uTMDPDF_50_Ev=Rkernel*uTMDPDF_lowScale50(x,bT,hadron)
    !uTMDPDF_50_Ev(0)=uTMDPDF_50_Ev(0)*RkernelG/Rkernel
    uTMDPDF_50_Ev=uTMDPDF_lowScale50(x,bT,hadron)*&
      (/Rkernel,Rkernel,Rkernel,Rkernel,Rkernel,RkernelG,Rkernel,Rkernel,Rkernel,Rkernel,Rkernel/)
    
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
       mui=c3_global*mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,c1_global*mu0(bt),1)
   CASE(2)!!!! improved gamma
       mui=c3_global*mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
   CASE(3)!!!! fixed mu
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
   CASE(4)!!!! exact solution via zeta-line
      mui=c3_global*mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
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
       mui=c3_global*mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,c1_global*mu0(bt),1)
   CASE(2)!!!! improved gamma
       mui=c3_global*mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
   CASE(3)!!!! fixed mu
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
   CASE(4)!!!! exact solution via zeta-line
      mui=c3_global*mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
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
       mui=c3_global*mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,c1_global*mu0(bt),1)
   CASE(2)!!!! improved gamma
       mui=c3_global*mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
   CASE(3)!!!! fixed mu
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
   CASE(4)!!!! exact solution via zeta-line
      mui=c3_global*mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
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
      mui=c3_global*mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,c1_global*mu0(bt),1)
      RkernelG=TMDR_Rzeta(bt,muf,zetaf,mui,c1_global*mu0(bt),0)
   CASE(2)!!!! improved gamma
       mui=c3_global*mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
      RkernelG=TMDR_Rzeta(bt,muf,zetaf,mui,0)
   CASE(3)!!!! fixed mu
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
      RkernelG=TMDR_Rzeta(bt,muf,zetaf,0)
   CASE(4)!!!! exact solution via zeta-line
      mui=c3_global*mu_LOW(bt)
      Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
      RkernelG=TMDR_Rzeta(bt,muf,zetaf,mui,0)
   END SELECT
    !uTMDFF_50_Ev=Rkernel*uTMDFF_lowScale50(x,bT,hadron)
    !uTMDFF_50_Ev(0)=uTMDFF_50_Ev(0)*RkernelG/Rkernel
    uTMDFF_50_Ev=uTMDFF_lowScale50(x,bT,hadron)*&
      (/Rkernel,Rkernel,Rkernel,Rkernel,Rkernel,RkernelG,Rkernel,Rkernel,Rkernel,Rkernel,Rkernel/)
    
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
  
  
  !!!!!!!!!!!!!!!!!!!lpTMDPDF!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  !!!!!!!! linearly polarized gluon TMDPDF
  ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  ! all quark terms are zero!
  function lpTMDPDF_50_Ev(x,bt,muf,zetaf,hadron)
    real*8::lpTMDPDF_50_Ev(-5:5)
    real*8:: x,bt,muf,zetaf
    real*8:: mui,RkernelG   
    integer::hadron
    
   SELECT CASE(EvolutionType)
   CASE(1)!!!! improved D
      mui=c3_global*mu_LOW(bt)
      RkernelG=TMDR_Rzeta(bt,muf,zetaf,mui,c1_global*mu0(bt),0)
   CASE(2)!!!! improved gamma
      mui=c3_global*mu_LOW(bt)      
      RkernelG=TMDR_Rzeta(bt,muf,zetaf,mui,0)
   CASE(3)!!!! fixed mu
      RkernelG=TMDR_Rzeta(bt,muf,zetaf,0)
   CASE(4)!!!! exact solution via zeta-line
      mui=c3_global*mu_LOW(bt)
      RkernelG=TMDR_Rzeta(bt,muf,zetaf,mui,0)
   END SELECT
    lpTMDPDF_50_Ev=RkernelG*lpTMDPDF_lowScale50(x,bT,hadron)
    
    !!! forcefully set=0 all quarks
!     glTMDPDF_50_Ev(5)=0d0
!     glTMDPDF_50_Ev(4)=0d0
!     glTMDPDF_50_Ev(3)=0d0
!     glTMDPDF_50_Ev(2)=0d0
!     glTMDPDF_50_Ev(1)=0d0
!     glTMDPDF_50_Ev(-1)=0d0
!     glTMDPDF_50_Ev(-2)=0d0
!     glTMDPDF_50_Ev(-3)=0d0
!     glTMDPDF_50_Ev(-4)=0d0
!     glTMDPDF_50_Ev(-5)=0d0
    
  end function lpTMDPDF_50_Ev
  
    !!!!!!!! linearly polarized gluon TMDPDF OPTIMAL
  ! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
  ! all quark terms are zero!
  function lpTMDPDF_50_optimal(x,bt,hadron)
    real*8::lpTMDPDF_50_optimal(-5:5)
    real*8:: x,bt
    integer::hadron
    
    lpTMDPDF_50_optimal=lpTMDPDF_lowScale50(x,bT,hadron)
    
  end function lpTMDPDF_50_optimal
  
  
  
end module TMDs