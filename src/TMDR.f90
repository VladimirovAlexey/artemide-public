!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!            arTeMiDe 2.00
!
!    Evaluation of the TMD evolution kernel
!    Here we use the improved gamma-solution, and the universal TMD definition.
!
!    if you use this module please, quote 1803.11089
!
!                A.Vladimirov (17.04.2018)
!            v1.32   A.Vladimirov (30.08.2018)
!                b-freeze at 1d-6 A.Vladimirov (16.09.2018)
!            v1.41   transpose-issue fixed A.Vladimirov (11.03.2019)
!                29.03.2019  Update to version 2.00 (AV).
!            v2.01     Added zeta-line with non-pertrubative D A.Vladimirov (06.06.2019)
!                Added gluon evolution A.Vladimirov (12.06.2019)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDR
use aTMDe_Numerics
use IO_functions
use QCDinput
use TMDR_model
use TMD_AD

implicit none

private

!!!!!! 0=STANDARD MODE
!!!!!! 1=CSS mode - former type 1 (ONLY FOR TESTING)
#define EVOLUTION_MODE 0

!Current version of module
character (len=5),parameter :: version="v3.00"
character (len=7),parameter :: moduleName="TMDR"
!Last appropriate verion of constants-file
integer,parameter::inputver=30

!------------------------------------------Working variables------------------------------------------------------------

!!!! the orders of used perturbative expressions for AD
!! Orders of anomalous dimensions as^n
integer::orderCusp      !for Cusp AD
integer::orderV         !for gammaV
integer::orderD         !for RAD
integer::orderDresum    !for resummed RAD
integer::orderZETA      !for zeta-line

!! Level of output
!! 0=only critical
!! 1=initialization details
!! 2=WARNINGS
integer::outputLevel=2
integer::messageTrigger=5
logical::started=.false.

!! Scale variation of OPE
real(dp)::c1_global=1._dp
!! Precision tolerance used in various routines
real(dp)::tolerance=0.0001d0
!! The lowerst value of b (below the expression is freezed)
real(dp)::bFREEZE=1d-6
!! Parameter that interpolaes the small-b part of zeta-line
real(dp)::smoothingParameter=0.01d0

integer::counter,messageCounter

!------------------------------------------Non-pertrubative parameters--------------------------------------------------

!!Number of non-pertrubative parameters
integer::NPlength=0
!! array of non-pertrubative parameters
real(dp),allocatable,dimension(:):: NPparam

public::TMDR_R,TMDR_Rzeta,TMDR_Rzeta_harpy
public:: TMDR_Initialize,TMDR_setNPparameters,LowestQ,TMDR_IsInitialized,TMDR_CurrentNPparameters,TMDR_SetScaleVariation
public:: CS_kernel,zetaNP

  
contains

#if INTEGRATION_MODE==1
    !!! THE EVOLUTION BY CSS-formula
    INCLUDE 'Code/TMDR/typeCSS.f90'
#endif


function TMDR_IsInitialized()
    logical::TMDR_IsInitialized
    TMDR_IsInitialized=started
end function TMDR_IsInitialized

   
!!! Initializing routing
!!! Filles the prebuiled arrays
!!! orderAD, is order of anomalous dimension for evolution
!!!! order zeta is the order of pertrubative zeta expression, used only in the "universal TMD"
subroutine TMDR_Initialize(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=300)::path
    logical::initRequired
    character(len=8)::orderMain,orderI
    logical::overrideORDER=.false.
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
      CLOSE (51, STATUS='KEEP')
    end if
    call MoveTO(51,'*p2  ')
    read(51,*) outputLevel    
    if(outputLevel>2) write(*,*) '--------------------------------------------- '
    if(outputLevel>2) write(*,*) 'artemide.',moduleName,version,': initialization started ... '
    call MoveTO(51,'*p3  ')
    read(51,*) messageTrigger
    
    
    call MoveTO(51,'*3   ')
    call MoveTO(51,'*p1  ')
    read(51,*) initRequired
    if(.not.initRequired) then
      if(outputLevel>2) write(*,*)'artemide.',moduleName,': initialization is not required. '
      started=.false.
      return
    end if
    
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) orderMain
    call MoveTO(51,'*p2  ')
    read(51,*) overrideORDER
    
    if(.not.overrideORDER) then !!! usual definition
    !!!!! IMPORTANT
    !!!!! Gamma cusp start from Gamma0, i.e. orderCusp=0 = as^1    
    SELECT CASE(orderMain)
      CASE ("LO")
    if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: LO'
    orderCusp=0
    orderV=0
    orderD=0
    orderDresum=0
    orderZETA=0
      CASE ("NLO")
    if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NLO'
    orderCusp=1
    orderV=1
    orderD=1
    orderDresum=1
    orderZETA=1
      CASE ("NNLO")
    if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NNLO'
    orderCusp=2
    orderV=2
    orderD=2
    orderDresum=2
    orderZETA=2
      CASE ("N2LO")
    if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NNLO'
    orderCusp=2
    orderV=2
    orderD=2
    orderDresum=2
    orderZETA=2
      CASE ("NNNLO")
    if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NNNLO'
    orderCusp=3
    orderV=3
    orderD=3
    orderDresum=3
    orderZETA=3
      CASE ("N3LO") !!! same as NNNLO
    if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: N3LO'
    orderCusp=3
    orderV=3
    orderD=3
    orderDresum=3
    orderZETA=3
      CASE ("N4LO")
    if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: N4LO'
    orderCusp=4
    orderV=4
    orderD=4
    orderDresum=4
    orderZETA=4
      CASE DEFAULT
    if(outputLevel>0) write(*,*) &
                WarningString(' Initialize:try to set unknown ADs orders. Switch to NLO.',modulename)
    if(outputLevel>1) write(*,*)  trim(moduleName)//' Order set: NLO'
    orderCusp=1
    orderV=1
    orderD=1
    orderDresum=1
    orderZETA=1
    END SELECT

      if(outputLevel>2) then
      write(*,'(A,I1)') ' |  GammaCusp     =as^',orderCusp+1
      write(*,'(A,I1)') ' |  gammaV        =as^',orderV
      write(*,'(A,I1)') ' |  D             =as^',orderD
      write(*,'(A,I1)') ' |  Dresum        =as^',orderDresum
      write(*,'(A,I1)') ' |  zeta_mu       =as^',orderZETA
     end if

    else !!! overrided definition
      if(outputLevel>1) write(*,*) trim(moduleName)//' Starndard order-definition is overrided'

      call MoveTO(51,'*p3  ')
      read(51,*) orderI
      orderCusp=OrderSelector(orderI)
      if(orderCusp<0 .and. outputLevel>0)  then
        write(*,*) WarningString(&
        ' Initialize:try to set unknown CUSP order. Switch to NLO.',modulename)
        orderCusp=1
      end if

      call MoveTO(51,'*p4  ')
      read(51,*) orderI
      orderV=OrderSelector(orderI)
      if(orderV<0 .and. outputLevel>0)  then
        write(*,*) WarningString(&
        ' Initialize:try to set unknown gammaV order. Switch to NLO.',modulename)
        orderV=1
      end if

      call MoveTO(51,'*p5  ')
      read(51,*) orderI
      orderD=OrderSelector(orderI)
      if(orderD<0 .and. outputLevel>0)  then
        write(*,*) WarningString(&
        ' Initialize:try to set unknown RAD order. Switch to NLO.',modulename)
        orderD=1
      end if

      orderDresum=orderD

      call MoveTO(51,'*p6  ')
      read(51,*) orderI
      orderZETA=OrderSelector(orderI)
      if(orderZETA<0 .and. outputLevel>0)  then
        write(*,*) WarningString(&
        ' Initialize:try to set unknown Zeta-line order. Switch to NLO.',modulename)
        orderZETA=1
      end if

      if(outputLevel>1) then
      write(*,'(A,I1)') ' |  GammaCusp     =as^',orderCusp+1
      write(*,'(A,I1)') ' |  gammaV        =as^',orderV
      write(*,'(A,I1)') ' |  D             =as^',orderD
      !write(*,'(A,I1)') ' |  Dresum        =as^',orderDresum
      write(*,'(A,I1)') ' |  zeta_mu       =as^',orderZETA
     end if
    end if

    !-------------------------
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) NPlength
    
    if(outputLevel>2) write(*,'(A,I3)') ' Number of NP parameters =',NPlength
    
    !--------------------------
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) tolerance
    call MoveTO(51,'*p2  ')
    read(51,*) bFREEZE
    call MoveTO(51,'*p3  ')
    read(51,*) smoothingParameter
    
    if(outputLevel>2) write(*,'(A,ES10.3)') ' |   tolerance=',tolerance
    
    CLOSE (51, STATUS='KEEP') 
    
    call TMD_AD_Initialize(orderCusp,orderV,orderD,orderDresum,orderZETA)
    
    
    if(.not.QCDinput_IsInitialized()) then
      if(outputLevel>1) write(*,*) '.. initializing QCDinput (from ',moduleName,')'
      if(present(prefix)) then
          call QCDinput_Initialize(file,prefix)
      else
    call QCDinput_Initialize(file)
      end if
    end if
    
    if(outputLevel>2) write(*,*) 'Model initialization..'
    allocate(NPparam(1:NPlength))
    call ModelInitialization(NPlength)

    c1_global=1._dp

#if INTEGRATION_MODE==1
    if(outputLevel>0) write(*,*) color('  CSS-like evolution path is selected',c_red)
    if(outputLevel>0) write(*,*) color('  !!! IT IS NOT A DEFAULT OPTION !!! ',c_red)
#endif
    
    started=.true.
    counter=0
    messageCounter=0
    if(outputLevel>0) write(*,*) color('----- arTeMiDe.TMDR '//trim(version)//': .... initialized',c_green)
    if(outputLevel>1) write(*,*) ' '
end subroutine TMDR_Initialize

!!! function that associates LO=0, NLO=1, etc.
pure function OrderSelector(orderSTR)
  character*8,intent(in)::orderSTR
  integer::OrderSelector
      SELECT CASE(orderSTR)
      CASE ("LO")
      OrderSelector=0
      CASE ("NLO")
      OrderSelector=1
      CASE ("NNLO")
      OrderSelector=2
      CASE ("N2LO")
      OrderSelector=2
      CASE ("NNNLO")
      OrderSelector=3
      CASE ("N3LO")
      OrderSelector=3
      CASE ("N4LO")
      OrderSelector=4
      CASE DEFAULT
      OrderSelector=-1
    END SELECT
end function OrderSelector

!Reset the NP parameters
subroutine TMDR_setNPparameters(lambdaIN)
    real(dp),intent(in)::lambdaIN(:)
    integer::ll
    
    ll=size(lambdaIN)
    if(ll<NPlength) then 
      if(outputLevel>0) write(*,"(A,I3,A,I3,')')")&
          WarningString('length of lambdaNP(',moduleName),ll,&
          color(') is less then requred (',c_red),NPlength
      if(outputLevel>0) write(*,*)color('                Rest parameters are replaced by zeros!',c_red)
      NPparam=0d0*NPparam
      NPparam(1:ll)=lambdaIN(1:ll)
    else if (ll>NPlength) then
      if(outputLevel>0) write(*,"(A,I3,A,I3,')')")&
          WarningString('length of lambdaNP(',moduleName),ll,&
          color(') is greater then requred (',c_red),NPlength
      if(outputLevel>0) write(*,*)color('                Array is truncated!',c_red)
      NPparam(1:NPlength)=lambdaIN(1:NPlength)
     else
      NPparam=lambdaIN
     end if
    
    if(outputLevel>2) write(*,*) 'arTeMiDe.TMDR: NPparameters reset = (',NPparam,')'
    
    call ModelUpdate(NPparam)
    
end subroutine TMDR_setNPparameters

subroutine TMDR_CurrentNPparameters(var)
    real(dp),dimension(1:NPlength)::var
    var=NPparam
 
end subroutine TMDR_CurrentNPparameters


!!!-------------------------------------------------------------------------------------------------------
!!!-----------------------------Routines for the evolution of type 3--------------------------------------
!!!----------------------------------(zeta-prescription)--------------------------------------------------
!!!-------------------------------------------------------------------------------------------------------

!!! the value of CS-kernel at mu,b for flavor f
function CS_kernel(mu,b,f)
  real(dp)::CS_kernel
  real(dp),intent(in)::mu,b
  integer,intent(in)::f

  real(dp)::muValue

  muValue=muOPE(b,c1_global)

  CS_kernel=Dpert(muValue,bSTAR(b),f)+RADEvolution(muValue,mu,f)+DNP(b,f)

end function CS_kernel

!! This is the non-pertrubative shape of zeta_mu line.
!! It MUST follow the equipotential line in perturbative regime (at small-b), at the level pf PT accuracy.
!! Otherwice, your evolution is completely broken.
!! Use zetaMUpert for perturbative values, use zetaSL for exact values
function zetaNP(mu,b,f)
    real(dp)::zetaNP
    real(dp),intent(in)::mu,b
    integer,intent(in)::f
    real(dp)::rad

    rad=CS_kernel(mu,b,f)
    zetaNP=zetaNP_rad(mu,rad,b,f)
end function zetaNP

function zetaNP_rad(mu,rad,b,f)
    real(dp)::zetaNP_rad
    real(dp),intent(in)::mu,b,rad
    integer,intent(in)::f
    real(dp)::zz

    !! this ofset is required to guaranty a good numerical behavior at b->0.
    !! In principle, zz=0 also works
    zz=Exp(-b**2/smoothingParameter)

    zetaNP_rad=zetaMUpert(mu,b,f)*zz+zetaSL(mu,rad,f)*(1d0-zz)

end function zetaNP_rad

 !!! general interface for evolutions from (muf,zetaf)->(mui,zetai)
function TMDR_R(b,muf,zetaf,mui,zetai,f)
  real(dp)::TMDR_R
  real(dp),intent(in)::b,muf,zetaf,mui,zetai
  integer,intent(in)::f
  TMDR_R=TMDR_Rzeta(b,muf,zetaf,f)/TMDR_Rzeta(b,mui,zetai,f)
end function TMDR_R


!!! Evolution exponent in the improved gamma-picture to zeta-line (defined by zetaNP)
function TMDR_Rzeta(b,muf,zetaf,f)
  real(dp), intent(in)::b,muf,zetaf
  integer,intent(in)::f
  real(dp)::TMDR_Rzeta,bLocal

  real(dp)::Dvalue,zetaP,muLOW

  if(b<bFREEZE) then
    bLocal=bFREEZE
  else
    bLocal=b
  end if

#if INTEGRATION_MODE==0
  !!!!! THIS IS DEFAULT VERSION
    Dvalue=CS_kernel(muf,b,f)

    zetaP=zetaNP_rad(muf,Dvalue,bLocal,f)

    TMDR_Rzeta=EXP(-Dvalue*Log(zetaf/zetaP))
#elif INTEGRATION_MODE==1
    !!!!! THIS IS FOR CSS-like evolution !!!!!
    muLOW=muOPE(b,c1_global)
    TMDR_Rzeta=TMDR_Rzeta_type1(b,muf,zetaf,muLOW,muLOW,f)
#endif

  if(TMDR_Rzeta>1d6) then
    write(*,*) ErrorString('Evolution factor is TOO HUGE check the formula',moduleName)
    write(*,*) 'b=',bLocal,'zetaf=',zetaf,'muf=',muf,'zetaEXACT=',zetaP
    write(*,*) 'c1=',c1_global,'log(zeta/zetamu)=',Log(zetaf/zetaP)
    write(*,*) 'NPparameters= (',NPparam,')'
    write(*,*) 'muOPE=',muOPE(b,c1_global), 'b*=',bSTAR(bLocal)
    write(*,*) 'Dpert=',Dpert(muOPE(b,c1_global),bSTAR(bLocal),f),&
    '\int G =',RADEvolution(muOPE(b,c1_global),muf,f), 'DNP=',DNP(bLocal,f)
    write(*,*) 'Evaluation stop'
    stop
  end if

end function TMDR_Rzeta

function TMDR_Rzeta_harpy(b,muf,zetaf,f)
  real(dp)::TMDR_Rzeta_harpy,b,muf,zetaf
  integer::f

  TMDR_Rzeta_harpy=TMDR_Rzeta(b,muf,zetaf,f)
end function TMDR_Rzeta_harpy

!!!!!!!!!! ------------------------ SUPPORINTG ROUTINES --------------------------------------

!!!! this routine set the variations of scales
!!!! it is used for the estimation of errors
subroutine TMDR_SetScaleVariation(c1_in)
    real(dp),intent(in)::c1_in
    if(c1_in<0.1d0 .or. c1_in>10.d0) then
        if(outputLevel>0) write(*,*) WarningString('variation in c1 is enourmous. c1 is set to 2',moduleName)
        c1_global=2d0
    else if(abs(c1_in-c1_global)<tolerance) then
        if(outputLevel>1) write(*,*) color('TMDR: c1-variation is ignored. c1='//real8ToStr(c1_global),c_yellow)
    else
        c1_global=c1_in
        if(outputLevel>1) write(*,*) color('TMDR: set scale variations c1 as:'//real8ToStr(c1_global),c_yellow)
    end if
end subroutine TMDR_SetScaleVariation

 !! we search for the lowest avalible Q, under which the evolution is inverted (for quark)
  !! it is defined by Q^2<zeta_Q(b->infty)
  !! with variabtions it is Q^2<zeta_{c Q}(b->infty)
  !! the function return the values for c={0.5,1,2} at b=25 (let it be assimptotic value)
  !! The function returns 1 if Q<1
 function LowestQ()
  real(dp),dimension(1:3)::LowestQ
  
  real(dp),parameter::b=25d0
  integer,parameter::f=1
  real(dp)::Qs1,Qs2,Qs3,V1,V2,V3
  
  !c=1
  Qs1=1d0
  V1=Qs1**2-zetaNP(Qs1,b,f)
  
  Qs3=25d0
  V3=Qs3**2-zetaNP(Qs3,b,f)
  
  if(V1<0d0 .and. V3>0d0) then
    do
      Qs2=(Qs3+Qs1)/2d0
      V2=Qs2**2-zetaNP(Qs2,b,f)
      if(V2>0d0) then 
    Qs3=Qs2
    V3=V2
      else
    Qs1=Qs2
    V1=V2
      end if
      if(Qs3-Qs1<tolerance) exit
    end do
    LowestQ(2)=(Qs3+Qs1)/2d0
  else if(V1>0d0 .and. V3>0d0) then
    LowestQ(2)=1d0
  else
    write(*,*) ErrorString('LowestQ: Negative value at large Q (c=1).',moduleName)
    stop
  end if
  
    !c=0.5
  Qs1=1d0
  V1=Qs1**2-zetaNP(0.5d0*Qs1,b,f)
  
  Qs3=25d0
  V3=Qs3**2-zetaNP(0.5d0*Qs3,b,f)
  
  if(V1<0d0 .and. V3>0d0) then
    do
      Qs2=(Qs3+Qs1)/2d0
      V2=Qs2**2-zetaNP(0.5d0*Qs2,b,f)
      if(V2>0d0) then 
    Qs3=Qs2
    V3=V2
      else
    Qs1=Qs2
    V1=V2
      end if
      if(Qs3-Qs1<tolerance) exit
    end do
    LowestQ(1)=(Qs3+Qs1)/2d0
  else if(V1>0d0 .and. V3>0d0) then
    LowestQ(1)=1d0
  else
    write(*,*) ErrorString('LowestQ: Negative value at large Q (c=1).',moduleName)
    stop
  end if
  
    !c=2
  Qs1=1d0
  V1=Qs1**2-zetaNP(2d0*Qs1,b,f)
  
  Qs3=25d0
  V3=Qs3**2-zetaNP(2d0*Qs3,b,f)
  
  if(V1<0d0 .and. V3>0d0) then
    do
      Qs2=(Qs3+Qs1)/2d0
      V2=Qs2**2-zetaNP(2d0*Qs2,b,f)
      if(V2>0d0) then 
    Qs3=Qs2
    V3=V2
      else
    Qs1=Qs2
    V1=V2
      end if
      if(Qs3-Qs1<tolerance) exit
    end do
    LowestQ(3)=(Qs3+Qs1)/2d0
  else if(V1>0d0 .and. V3>0d0) then
    LowestQ(3)=1d0
  else
    write(*,*) ErrorString('LowestQ: Negative value at large Q (c=1).',moduleName)
    stop
  end if
 
 end function LowestQ


end module TMDR
