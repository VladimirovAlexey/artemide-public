!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 1.4
!
!	Evaluation of the TMDs
!	
!	if you use this module please, quote 1706.01473
!
!	ver 1.0: release (AV, 10.05.2017)
!	ver 1.4: release (AV, 23.12.2018)
!	ver 1.41:release (AV, 23.12.2018)
!
!				A.Vladimirov (23.12.2018)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDs
  use TMDR
  use uTMDPDF
  use uTMDFF
  implicit none

  private
!   public
 
 !Current version of module
 character (len=5),parameter :: version="v1.40"
 
!------------------------------------------Physical and mathematical constants------------------------------------------

!!!Threashold parameters
! Set in the constants file
  real*8 :: mCHARM=1.4d0
  real*8 :: mBOTTOM=4.75d0
  
!------------------------------------------Working variables------------------------------------------------------------
  
  logical::started=.false.
  
  !!!! The TMD evolution can be tritted differently
  !!!! Originally we considered evolution along the path 1  (=1)
  !!!! Currently we adapt the solution independent evolution (=2)
  !!!! The universal TMD definition also requares extra constant (=3)
  integer::EvolutionType
  !! Level of output
  !! 0=only critical
  !! 1=initialization details
  !! 2=WARNINGS
  integer::outputLevel=2
  
  !!!parameters for the uncertanty estimation
  real*8::c1_global,c3_global,c4_global
  
  logical::MakeGrid,includeGluon!!! These parameters are for the grid precalculation
  
  !!! this is lenght of array of non-pertrubative parameters
  !!! its length is the number of TMDs modules
  !!! if entry=0 this module is not called
  integer,dimension(0:2)::parameterLength
  !!! total number of NP parameters, to check the input
  integer::totalNumberOfPar=0
  
  logical:: convergenceLost=.false.
  
!-----------------------------------------Public interface--------------------------------------------------------------
  public::TMDs_SetNPParameters,TMDs_SetScaleVariations,TMDs_Initialize,TMDs_SetPDFreplica
  real*8,dimension(-5:5),public:: uTMDPDF_5,uTMDPDF_50
  real*8,dimension(-5:5),public:: uTMDFF_5,uTMDFF_50
  
  public::uPDF_uPDF,uPDF_anti_uPDF
  public::TMDs_convergenceISlost,TMDs_IsconvergenceLost
  
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
  
  interface TMDs_SetNPParameters
    module procedure TMDs_SetNPParameters,TMDs_SetNPParameters_rep
  end interface
  
  contains 
  
   INCLUDE 'Model/TMDs_model.f90'
  
     !!! This subroutine can be called only ones per programm run in the very beginning
  !!! It initializes TMDRm abd TMDconvolution subpackages
  !!! It set the pertrubative orders(!), so pertrubative orders cannot be changed afterwards.
  !!! It also set the paths forPDF and As greeds according to orderPDF (NLO or NNLO).
  subroutine TMDs_Initialize(orderMain)
    character(len=*)::orderMain
    character(256)::line
    integer:: j
    
    OPEN(UNIT=51, FILE='constants', ACTION="read", STATUS="old")    
    !!! Search for output level
    do
    read(51,'(A)') line    
    if(line(1:3)=='*0 ') exit
    end do    
    do
    read(51,'(A)') line
    if(line(1:3)=='*A ') exit
    end do
    read(51,'(A)') line
    read(51,*) outputLevel
    
    if(outputLevel>1) write(*,*) '----- arTeMiDe.TMDs ',version,': .... initialization'
    
    if(started) then
    if(outputLevel>1) write(*,*) 'arTeMiDe.TMDs already initialized'
    else
    
    !!!!search for physic constants entry
    do
    read(51,'(A)') line    
    if(line(1:3)=='*1 ') exit
    end do    
    !!!!search for masses entry entry
    do
    read(51,'(A)') line
    if(line(1:3)=='*A ') exit
    end do
    read(51,'(A)') line
    read(51,*) mCHARM
    read(51,'(A)') line
    read(51,*) mBOTTOM   
    
    !!!!search for grid constants entry
    do
    read(51,'(A)') line    
    if(line(1:3)=='*3 ') exit
    end do    
    !!!!search for uTMDPDF entry
    do
    read(51,'(A)') line
    if(line(1:3)=='*B ') exit
    end do
    read(51,'(A)') line
    read(51,*) MakeGrid
    read(51,'(A)') line
    read(51,*) includeGluon 
    
    if(outputLevel>1) write(*,*) ' Griding = ',MakeGrid,' includeGluon = ',includeGluon
    
    if(outputLevel>1) write(*,*) 'Number of NP parameters: '
    !!!!search for length of parameters
    do
    read(51,'(A)') line    
    if(line(1:3)=='*4 ') exit
    end do
    do
    read(51,'(A)') line    
    if(line(1:3)=='*A ') exit
    end do
    !TMDevolution
    do
    read(51,'(A)') line    
    if(line(1:3)=='*0)') exit
    end do
    read(51,*) j
    parameterLength(0)=j
    if(outputLevel>1) write(*,'(A,I2)') ' | TMD evolution	: ',parameterLength(0)
    !uTMDPDF
    do
    read(51,'(A)') line    
    if(line(1:3)=='*1)') exit
    end do
    read(51,*) j
    parameterLength(1)=j
    if(outputLevel>1) write(*,'(A,I2)') ' | unpolarized TMDPDF	: ',parameterLength(1)
    !uTMDFF
    do
    read(51,'(A)') line    
    if(line(1:3)=='*2)') exit
    end do
    read(51,*) j
    parameterLength(2)=j
    if(outputLevel>1) write(*,'(A,I2)') ' | unpolarized TMDFF	: ',parameterLength(2)
    
    !!!! Selecting the evolution type
    do
    read(51,'(A)') line    
    if(line(1:3)=='*B ') exit
    end do    
    do
    read(51,'(A)') line    
    if(line(1:3)=='*1)') exit
    end do    
    read(51,*) EvolutionType
    if(outputLevel>1) then
      if(EvolutionType==1) write(*,*) ' Selected evolution type: improved D (CSS-like)'
      if(EvolutionType==2) write(*,*) ' Selected evolution type: improved gamma'
      if(EvolutionType==3) write(*,*) ' Selected evolution type: fixed mu (optimal)'
    end if
    
    CLOSE (51, STATUS='KEEP') 
    
    if(outputLevel>1) write(*,*) 'Initialization of submodules ...'
    call TMDR_Initialize(orderMain)
    if(parameterLength(1)>0) then
      call uTMDPDF_Initialize(orderMain)
    else
      if(outputLevel>1) write(*,*)  'uTMDPDF is not used.'
    end if
    if(parameterLength(2)>0) then
      call uTMDFF_Initialize(orderMain)
    else
      if(outputLevel>1) write(*,*)  'uTMDFF is not used.'
    end if
    
    totalNumberOfPar=sum(parameterLength)
    
    if(MakeGrid) then
       if(outputLevel>1) write(*,*) 'TMDs:TMD griding is ON'
     end if
    
      c1_global=1d0
      c3_global=1d0
      c4_global=1d0
      
      
      started=.true.
     if(outputLevel>0) write(*,*) '----- arTeMiDe.TMDs ',version,': .... initialized'
     if(outputLevel>1) write(*,*) ' '
    end if
  end subroutine TMDs_Initialize

  !!!!!!!Functions which carry the trigger on convergences.... Its used in xSec, and probably in other places.
  function TMDs_IsconvergenceLost()
  logical::TMDs_IsconvergenceLost
  TMDs_IsconvergenceLost=convergenceLost
  end function TMDs_IsconvergenceLost
  
  subroutine TMDs_convergenceISlost()
  convergenceLost=.true.
  if(outputLevel>1) write(*,*) 'arTeMiDe.TMDs: convergence triger set to be lost.'
  end subroutine TMDs_convergenceISlost
  
    !sets the NP parameters of TMD model (bmax,gB,etc)
  subroutine TMDs_SetNPParameters(lambda)
    real*8,intent(in):: lambda(:)
    
    if(size(lambda)<totalNumberOfPar) then
      if(outputLevel>0) write(*,*) 'arTeMiDe.TMDs: WARNING! Length of NP parameters array (',size(lambda),&
	  ') is less then expected (',totalNumberOfPar,')'
      if(outputLevel>0) write(*,*) 'NOTHING IS DONE'
    else
    
    convergenceLost=.false.
    
   call TMDR_setNPparameters(lambda(1:parameterLength(0)))
   if(parameterLength(1)>0)&
	call uTMDPDF_SetLambdaNP(lambda(parameterLength(0)+1:parameterLength(0)+parameterLength(1)) ,MakeGrid,includeGluon)
   if(parameterLength(2)>0)&
      call uTMDFF_SetLambdaNP(lambda(parameterLength(0)+parameterLength(1)+1:&
				      parameterLength(0)+parameterLength(1)+parameterLength(2)),&
				      MakeGrid,includeGluon)
   end if

  end subroutine TMDs_SetNPParameters
  
  !sets the NP parameters of TMD model for replica
  subroutine TMDs_SetNPParameters_rep(num)
    integer::num
    
    convergenceLost=.false.
    
   call TMDR_setNPparameters(num)
   
   if(parameterLength(1)>0) call uTMDPDF_SetLambdaNP(num,MakeGrid,includeGluon)
   if(parameterLength(2)>0) call uTMDFF_SetLambdaNP(num,MakeGrid,includeGluon)
  end subroutine TMDs_SetNPParameters_rep
  
  !!!! this routine set the variations of scales
  !!!! it is used for the estimation of errors
  subroutine TMDs_SetScaleVariations(c1_in,c3_in,c4_in)
    real*8::c1_in,c3_in,c4_in,c4_old
    
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
    
    c4_old=c4_global
    if(c4_in<0.1d0 .or. c4_in>10.d0) then
       if(outputLevel>0) write(*,*) 'WARNING: arTeMiDe.TMDs: variation in c4 is enourmous. c4 is set to 2'
      c4_global=2d0
    else
    c4_global=c4_in
    end if
    
    if(outputLevel>1) write(*,*) 'TMDs: set scale variation as:',c1_global,c3_global,c4_global
!     write(*,*) 'ADD INDEPENDENT VARIATIONS!! for TMDs'
    !!! we reset c4 only if it different from the old one. Otherwise it can lead to the grid recalculation 
    if(c4_old/=c4_in) then 
     if(parameterLength(1)>0) then
     call uTMDPDF_SetScaleVariation(c4_in)
     call uTMDPDF_resetGrid(MakeGrid,includeGluon)
     end if
     if(parameterLength(2)>0) then
     call uTMDFF_SetScaleVariation(c4_in)
     call uTMDFF_resetGrid(MakeGrid,includeGluon)
     end if
    end if
  end subroutine TMDs_SetScaleVariations
  
  !!!Pass variable to SetPDFreplica, reset grid
  subroutine TMDs_SetPDFreplica(rep)
  integer::rep
    convergenceLost=.false.
    call uTMDPDF_SetPDFreplica(rep)  
    call uTMDPDF_resetGrid(MakeGrid,includeGluon)
  end subroutine TMDs_SetPDFreplica
  
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

! program example
! use TMDs
!   implicit none
!   
!   call TMDs_Initialize('LO')
!   
!   call TMDs_SetNPParameters((/0.189d0,0.2d0,0.425d0,0.189d0,0.12d0,0.425d0/))
!   
!   write(*,*) uTMDPDF_5(0.1d0,1d0,911d0,91d0**2,1)
! !   write(*,*) uTMDPDF_50(0.1d0,1d0,911d0,91d0**2,1)
!   
!   write(*,*) uTMDFF_5(0.1d0,1d0,911d0,91d0**2,1)
! !   write(*,*) uTMDFF_50(0.1d0,1d0,911d0,91d0**2,1)
! 
! end program example 