!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 2.00
!
!	Evaluation of the TMD evolution kernel
!	Here we use the improved gamma-solution, and the universal TMD definition.
!	
!	if you use this module please, quote 1803.11089
!
!				A.Vladimirov (17.04.2018)
!			v1.32   A.Vladimirov (30.08.2018)
!				b-freeze at 1d-6 A.Vladimirov (16.09.2018)
!			v1.41   transpose-issue fixed A.Vladimirov (11.03.2019)
!				29.03.2019  Update to version 2.00 (AV).
!			v2.01 	Added zeta-line with non-pertrubative D A.Vladimirov (06.06.2019)
!				Added gluon evolution A.Vladimirov (12.06.2019)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDR
use aTMDe_Numerics
use IO_functions
use QCDinput
use TMDR_model
use TMD_AD

implicit none

private

!Current version of module
character (len=5),parameter :: version="v2.06"
character (len=7),parameter :: moduleName="TMDR"
!Last appropriate verion of constants-file
integer,parameter::inputver=10

character(256)::replicaFILE
logical::usereplicaFILE=.false.
character(50)::name

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

!! Precision tolerance used in various routines
real(dp)::tolerance=0.0001d0

!! Evolution type 
integer:: EvolutionType=3

integer::counter,messageCounter

!!!!!!!
!!model parameters
real(dp),allocatable::Stable(:,:)

!------------------------------------------Non-pertrubative parameters--------------------------------------------------

!!Number of non-pertrubative parameters
integer::NPlength=0
!! array of non-pertrubative parameters
real(dp),allocatable,dimension(:):: NPparam

public::TMDR_R,TMDR_Rzeta,TMDR_Rzeta_harpy
public:: TMDR_Initialize,TMDR_setNPparameters,LowestQ,TMDR_IsInitialized,TMDR_CurrentNPparameters

!public::DNP!,GammaCusp,gammaV



interface TMDR_setNPparameters
    module procedure TMDR_setNPparameters, TMDR_SetReplica
end interface 

interface TMDR_R
    module procedure TMDR_R_typeZ0,TMDR_R_typeZ
end interface

interface TMDR_Rzeta
    module procedure TMDR_Rzeta_type1,TMDR_Rzeta_type2,TMDR_Rzeta_type3
end interface
  
 contains

!!! Routines for evolution type1
INCLUDE 'Code/TMDR/type1.f90'
!!! Routines for evolution type2
INCLUDE 'Code/TMDR/type2.f90'
!!! Routines for evolution type3
INCLUDE 'Code/TMDR/type3.f90'

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
    character(len=300)::path,line
    logical::initRequired
    character(len=8)::orderMain
    integer::i,FILEver
    
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
      CASE ("LO+")
	if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: LO+'
	orderCusp=0
	orderV=1
	orderD=1
	orderDresum=0
	orderZETA=0
      CASE ("NLO")
	if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NLO'
	orderCusp=1
	orderV=1
	orderD=1
	orderDresum=1
	orderZETA=1
      CASE ("NLO+")
	if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NLO+'
	orderCusp=1
	orderV=2
	orderD=2
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
      CASE ("NNLO+")
	if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NNLO+'
	orderCusp=2
	orderV=3
	orderD=3
	orderDresum=3!2
	orderZETA=3!2
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
      CASE ("N3LO+")
	if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: N3LO+'
	orderCusp=3
	orderV=4
	orderD=4
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

    call MoveTO(51,'*p2  ')
    read(51,*) EvolutionType
    
    if(EvolutionType<1 .or. EvolutionType>3) then
        if(outputLevel>0) write(*,*) &
                WarningString(' Initialize: Inappropriate evolution type (not 1,2,3). Switch to 3',modulename)
        EvolutionType=3
    end if
    
    if(outputLevel>2) write(*,'(A,I3)') ' Evolution type =',EvolutionType
    
    
    !-------------------------
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) NPlength
    
    if(outputLevel>2) write(*,'(A,I3)') ' Number of NP parameters =',NPlength
    
    !!allocating number of NP input
    allocate(NPparam(1:NPlength))
    call MoveTO(51,'*p2  ')
    do i=1,NPlength
      read(51,*) NPparam(i)
    end do
    
    !--------------------------
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) tolerance
    
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
    call ModelInitialization(NPparam)
    
    !!! we also initialize the variables
    call TMDR_setNPparameters(NPparam)
    
    started=.true.
    counter=0
    messageCounter=0
    if(outputLevel>0) write(*,*) color('----- arTeMiDe.TMDR '//trim(version)//': .... initialized',c_green)
    if(outputLevel>1) write(*,*) ' '
end subroutine TMDR_Initialize

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

subroutine TMDR_SetReplica(rep)
    integer,intent(in)::rep    
    real(dp),allocatable::NParray(:)

    call GetReplicaParameters(rep,NParray)
    
    call TMDR_setNPparameters(NParray)
    
end subroutine TMDR_SetReplica 
 
subroutine TMDR_CurrentNPparameters(var)
    real(dp),dimension(1:NPlength)::var
    var=NPparam
 
end subroutine TMDR_CurrentNPparameters

 !!!--------------------------------------------------------------------------
 !!!--------------------------- interfaces -----------------------------------
 
 !!! general interface for evolutions from (muf,zetaf)->(mui,zetai)
 !!! contains parameter mu0 which is used in type-1 evolution.
 !!! for other cases mu0 is ignored
 !!! these functions are the part of TMDR_R definition
 function TMDR_R_typeZ0(b,muf,zetaf,mui,zetai,mu0,f)
  real(dp)::TMDR_R_typeZ0
  real(dp),intent(in)::b,muf,zetaf,mui,zetai,mu0
  integer,intent(in)::f
  
  SELECT CASE(EvolutionType)
    CASE(1)      
      TMDR_R_typeZ0=TMDR_R_type1(b,muf,zetaf,mui,zetai,mu0,f)
    CASE(2)
      TMDR_R_typeZ0=TMDR_R_type2(b,muf,zetaf,mui,zetai,f)
    CASE(3)
      TMDR_R_typeZ0=TMDR_Rzeta_type3(b,muf,zetaf,f)/TMDR_Rzeta_type3(b,mui,zetai,f)
    CASE DEFAULT
      TMDR_R_typeZ0=TMDR_Rzeta_type3(b,muf,zetaf,f)/TMDR_Rzeta_type3(b,mui,zetai,f)
    end SELECT
 
 end function TMDR_R_typeZ0
 
 !!! general interface for evolutions from (muf,zetaf)->(mui,zetai) for types 2,3,4
 !!! for the type 1 the parameter mu0 = mui
 !!! these functions are the part of TMDR_R definition
 function TMDR_R_typeZ(b,muf,zetaf,mui,zetai,f)
  real(dp)::TMDR_R_typeZ
  real(dp),intent(in)::b,muf,zetaf,mui,zetai
  integer,intent(in)::f
  
  SELECT CASE(EvolutionType)
    CASE(1)
      !!!!! parameter mu0 = mui !!!
      TMDR_R_typeZ=TMDR_R_type1(b,muf,zetaf,mui,zetai,mui,f)
    CASE(2)
      TMDR_R_typeZ=TMDR_R_type2(b,muf,zetaf,mui,zetai,f)
    CASE(3)
      TMDR_R_typeZ=TMDR_Rzeta_type3(b,muf,zetaf,f)/TMDR_Rzeta_type3(b,mui,zetai,f)
    CASE DEFAULT
      TMDR_R_typeZ=TMDR_Rzeta_type3(b,muf,zetaf,f)/TMDR_Rzeta_type3(b,mui,zetai,f)
    end SELECT
 
 end function TMDR_R_typeZ
 
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
