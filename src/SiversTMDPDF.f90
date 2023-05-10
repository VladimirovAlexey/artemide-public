!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 2.0
!
!	Evaluation of the Sivers TMD PDF f_{1T}^\perp at low normalization point in zeta-prescription.
!	
!	if you use this module please, quote 20??.????
!
!
!				A.Vladimirov (21.05.2020)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module SiversTMDPDF
use aTMDe_Numerics
use IO_functions
use QCDinput
use SiversTMDPDF_model

implicit none
!------------------------LOCALs -----------------------------------------------

private 

!Current version of module
character (len=5),parameter :: version="v2.05"
character (len=12),parameter :: moduleName="SiversTMDPDF"
!Last appropriate verion of constants-file
integer,parameter::inputver=12

!--------------------------------Working variables-----------------------------------------------
!--- general
logical:: started=.false.
!! Level of output
!! 0=only critical
!! 1=initialization details
!! 2=WARNINGS
integer::outputLevel=2
!! variable that count number of WRNING mesagges. In order not to spam too much
integer::messageTrigger=6

!!! The global order which is used in programm
! it is set by SetOrderForC subroutine
! 0=LO, 1=NLO, 2=NNLO
integer :: order_global

integer::lambdaNPlength
real(dp),dimension(:),allocatable::lambdaNP
real(dp),dimension(:),allocatable::lambdaNP_grid !!! this is the value of lambda on which the grid is build

real(dp)::c4_global!!!this is the variation constant for mu_OPE

!!!Parameters of numerics
real(dp) :: tolerance=0.0001d0!!! relative tolerance of the integration
integer :: maxIteration=5000

integer :: counter,messageCounter

INCLUDE 'Code/Twist3/Twist3Convolution-VAR.f90'
INCLUDE 'Code/Grids/TMDGrid-B-VAR.f90'

!!--------------------------------- variables for the griding the TMD.---------------------------------------------
logical :: gridReady!!!!indicator that grid is ready to use. If it is .true., the TMD calculated from the grid
logical :: prepareGrid!!!idicator that grid must be prepared
logical :: withGluon!!!indicator the gluon is needed in the grid
logical :: IsFnpZdependent !!! indicator that the grid must be recalculated with the change of Lambda

!!--------------------------------- variables for hadron composition---------------------------------------------
integer::numberOfHadrons				        !!!number of hadrons/components
integer,dimension(:),allocatable::hadronsInGRID	!!!list of hadron to be pre-grid
logical::IsComposite=.false.					!!!flag to use the composite TMD

!!-----------------------------------------------Public interface---------------------------------------------------

public::SiversTMDPDF_Initialize,SiversTMDPDF_SetLambdaNP,SiversTMDPDF_SetScaleVariation,SiversTMDPDF_resetGrid,&
    SiversTMDPDF_SetPDFreplica
public::SiversTMDPDF_IsInitialized,SiversTMDPDF_CurrentNPparameters
public::SiversTMDPDF_lowScale5,SiversTMDPDF_lowScale50
!   public::CheckCoefficient  
!   public::mu_OPE

interface SiversTMDPDF_SetLambdaNP
    module procedure SiversTMDPDF_SetLambdaNP_usual,SiversTMDPDF_SetReplica_optional
end interface
  
contains

INCLUDE 'Code/Twist3/Twist3Convolution.f90'
INCLUDE 'Code/Grids/TMDGrid-B-2.f90'


!! Testing the model
INCLUDE 'Code/SiversTMDPDF/modelTest.f90'
!! The convolution interfaces
INCLUDE 'Code/SiversTMDPDF/convolutions.f90'

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Interface subroutines!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

function SiversTMDPDF_IsInitialized()
    logical::SiversTMDPDF_IsInitialized
    SiversTMDPDF_IsInitialized=started
end function SiversTMDPDF_IsInitialized

!! Initialization of the package
subroutine SiversTMDPDF_Initialize(file,prefix)
    character(len=*)::file
    character(len=*),optional::prefix
    character(len=300)::path,line
    logical::initRequired
    character(len=8)::orderMain
    logical::bSTAR_lambdaDependent
    integer::i,FILEver

    if(started) return

    if(.not.QCDinput_IsInitialized()) then
        if(outputLevel>1) write(*,*) '.. initializing QCDinput (from ',moduleName,')'
        if(present(prefix)) then
        call QCDinput_Initialize(file,prefix)
        else
        call QCDinput_Initialize(file)
        end if
    end if

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
    if(outputLevel>1) write(*,*) '--------------------------------------------- '
    if(outputLevel>1) write(*,*) 'artemide.',moduleName,version,': initialization started ... '
    call MoveTO(51,'*p3  ')
    read(51,*) messageTrigger


    call MoveTO(51,'*12   ')
    call MoveTO(51,'*p1  ')
    read(51,*) initRequired
    if(.not.initRequired) then
        if(outputLevel>1) write(*,*)'artemide.',moduleName,': initialization is not required. '
        started=.false.
        return
    end if

    !----- ORDER
    call MoveTO(51,'*A   ')
    call MoveTO(51,'*p1  ')
    read(51,*) orderMain

    SELECT CASE(trim(orderMain))
        CASE ("NA")
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: NA',color(" (TMD=fNP)",c_yellow)
            order_global=-50
        CASE ("LO")
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: LO'
            order_global=0
        CASE ("LO+")
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: LO+'
            order_global=0        
        CASE DEFAULT
            if(outputLevel>0) write(*,*) &
                WarningString('Initialize: unknown order for coefficient function. Switch to LO.',moduleName)
            if(outputLevel>1) write(*,*) trim(moduleName)//' Order set: LO'
            order_global=0
    END SELECT

    if(outputLevel>2 .and. order_global>-1) write(*,'(A,I1)') ' |  Coef.func.    =as^',order_global
        
    !------ Compositeness
    call MoveTO(51,'*p2  ')
    read(51,*) IsComposite

    if(outputLevel>2) then
        if(IsComposite) then
        write(*,'(A,I1)') ' |  Use compsite  =TRUE'
        else
        write(*,'(A,I1)') ' |  Use compsite  =FALSE'
        end if
        end if

    !-------------parameters of NP model
    call MoveTO(51,'*B   ')
    call MoveTO(51,'*p1  ')
    read(51,*) lambdaNPlength

    if(outputLevel>2) write(*,'(A,I3)') ' Number of NP parameters =',lambdaNPlength

    if(lambdaNPlength<=0) then
    write(*,*) ErrorString('Initialize: number of non-pertrubative &
            parameters should be >=1. Check the constants-file. Evaluation STOP',moduleName)
    stop
    end if
    allocate(lambdaNP(1:lambdaNPlength))
    
    call MoveTO(51,'*p2  ')
    do i=1,lambdaNPlength
        read(51,*) lambdaNP(i)
    end do


    !-------------Numeric parameters
    call MoveTO(51,'*C   ')
    call MoveTO(51,'*p1  ')
    read(51,*) tolerance
    call MoveTO(51,'*p2  ')
    read(51,*) maxIteration

    if(outputLevel>2) then
        write(*,'(A,ES10.3)') ' |  tolerance     =',tolerance
        write(*,'(A,ES10.3)') ' |  max iteration =',REAL(maxIteration)
        end if
        
    !-------------Make grid options
    call MoveTO(51,'*D   ')
    call MoveTO(51,'*p1  ')
    read(51,*) prepareGrid    
    call MoveTO(51,'*p2  ')
    read(51,*) withGluon
    call MoveTO(51,'*p3  ')
    read(51,*) numberOfHadrons
    allocate(hadronsInGRID(1:numberOfHadrons))
    call MoveTO(51,'*p4  ')
    read(51,*) hadronsInGRID

    !-------------Parameters of grid
    call MoveTO(51,'*E   ')
    call MoveTO(51,'*p1  ')
    read(51,*) xGrid_Min
    call MoveTO(51,'*p2  ')
    read(51,*) bGrid_Max
    call MoveTO(51,'*p3  ')
    read(51,*) GridSizeX
    call MoveTO(51,'*p4  ')
    read(51,*) GridSizeB
    call MoveTO(51,'*p5  ')
    read(51,*) slope

        if(outputLevel>2) then
        write(*,*) 'Grid options:'
        write(*,'(A,ES10.3)') ' |  xGrid_Min                 =',xGrid_Min
        write(*,'(A,ES10.3)') ' |  bGrid_Max                 =',bGrid_Max 
        write(*,'(A,I6,A,I6,A)') ' |  (GridSizeX,GridSizeB)     =(',GridSizeX,',',GridSizeB,')'
        write(*,'(A,F6.3)') ' |  slope                     =',slope 
        write(*,'(A,I3)')   ' |  hadrons to grid           =',numberOfHadrons
        write(*,*)   ' | list of hadrons in grid    =(',hadronsInGRID,')'
        end if

    CLOSE (51, STATUS='KEEP') 

    allocate(gridMain(0:GridSizeX,0:GridSizeB,-5:5,1:numberOfHadrons))
    allocate(boundaryValues(0:GridSizeX,-5:5,1:numberOfHadrons))

    allocate(lambdaNP_grid(1:lambdaNPlength))

    c4_global=1d0

    call ModelInitialization(lambdaNP)

    !!!!!!!Checking the lambda-dependance of bSTAR
    bSTAR_lambdaDependent=testbSTAR()

    if(bSTAR_lambdaDependent) then
        if(outputLevel>2) write(*,*) 'arTeMiDe.SiversTMDPDF: bSTAR is dependent on lambda'
    else
        if(outputLevel>2) write(*,*) 'arTeMiDe.SiversTMDPDF: bSTAR is independent on lambda'
    end if

    !!! if fnp depende on z or bSTAR depeds on lambda
    !!! grid must be recalculate ech time. It canbe saved to single IsFnpZdependent
    if(bSTAR_lambdaDependent) then
        IsFnpZdependent=.true.
        if(outputLevel>2) write(*,*) 'arTeMiDe.SiversTMDPDF: ............. convolution is lambda sensitive.'
    end if

    started=.true.
    messageCounter=0

    if(outputLevel>0) write(*,*) color('----- arTeMiDe.SiversTMDPDF '//trim(version)//': .... initialized',c_green)
    if(outputLevel>1) write(*,*) ' '
end subroutine SiversTMDPDF_Initialize

!! call for parameters from the model
subroutine SiversTMDPDF_SetReplica_optional(num,buildGrid, gluonRequired)
    integer,intent(in):: num
    logical,optional,intent(in):: buildGrid,gluonRequired
    real(dp),allocatable::NParray(:)

    call GetReplicaParameters(num,NParray)

    if(present(buildGrid)) then
    if(present(gluonRequired)) then
        call SiversTMDPDF_SetLambdaNP_usual(NParray,buildGrid=buildGrid,gluonRequired=gluonRequired)
    else
        call SiversTMDPDF_SetLambdaNP_usual(NParray,buildGrid=buildGrid)
    end if
    else
    if(present(gluonRequired)) then
        call SiversTMDPDF_SetLambdaNP_usual(NParray,gluonRequired=gluonRequired)
    else
        call SiversTMDPDF_SetLambdaNP_usual(NParray)
    end if
    end if

end subroutine SiversTMDPDF_SetReplica_optional

!! call QCDinput to change the PDF replica number
!! unset the grid, since it should be recalculated fro different PDF replica.
subroutine SiversTMDPDF_SetPDFreplica(rep,hadron)
    integer,intent(in):: rep,hadron
    
    if(outputLevel>0) write(*,*) &
                WarningString('Cannot change PDF replica. No PDFs for Sivers. NOTHING IS DONE!',moduleName)
    
end subroutine SiversTMDPDF_SetPDFreplica

!!!Sets the non-pertrubative parameters lambda
!!! carries additionl option to build the grid
!!! if need to build grid, specify the gluon required directive.
subroutine SiversTMDPDF_SetLambdaNP_usual(lambdaIN,buildGrid, gluonRequired)
    real(dp),intent(in)::lambdaIN(:)
    logical,optional,intent(in) :: buildGrid,gluonRequired
    real(dp),dimension(1:lambdaNPlength)::lambdaOLD
    logical::IsNewValues
    integer::i,ll
    messageCounter=0

    if(present(buildGrid)) prepareGrid=buildGrid
    if(present(gluonRequired)) withGluon=gluonRequired

    ll=size(lambdaIN)
    if(ll<lambdaNPlength) then 
        if(outputLevel>0) write(*,"(A,I3,A,I3,')')")&
                WarningString('SetLambdaNP:length of lambdaNP(',moduleName),&
            ll,color(') is less then requred (',c_red),lambdaNPlength
        if(outputLevel>0) write(*,*)color('                Rest parameters are replaced by zeros!',c_red)
        lambdaNP=0d0*lambdaNP
        lambdaNP(1:ll)=lambdaIN(1:ll)
    else if (ll>lambdaNPlength) then
        if(outputLevel>0) write(*,"(A,I3,A,I3,')')")&
                WarningString('SetLambdaNP:length of lambdaNP(',moduleName),&
            ll,color(') is greater then requred (',c_red),lambdaNPlength
        if(outputLevel>0) write(*,*)color('                Array is truncated!',c_red)
        lambdaNP(1:lambdaNPlength)=lambdaIN(1:lambdaNPlength)
        else
        lambdaOLD=lambdaNP
        lambdaNP=lambdaIN
        end if


    lambdaOLD=lambdaNP
    lambdaNP=lambdaIN
    IsNewValues=.false.
    do i=1,lambdaNPlength
        if(ABS(lambdaNP(i)-lambdaOLD(i))>10d-10) then
        IsNewValues=.true.
        exit
        end if
    end do    


    if(IsNewValues.and.(outputLevel>2)) write(*,*) 'arTeMiDe.',moduleName,': NPparameters reset = (',lambdaNP,')'
    call ModelUpdate(lambdaNP)

    !! further if's are only for griding

    if(prepareGrid) then !!!grid is requred

        if(IsNewValues) then  !! values are new

            if(gridReady) then  !!! grid is already build
        
                if(IsFnpZdependent) then !!! check the z-dependance of FNP
                    !! if it is z- dependent, rebuild the grid
                    gridReady=.false.
                    call MakeGrid()
                    gridReady=.true.
        
                else !!! if z-Independent just do nothing.
                    if(outputLevel>2) write(*,*) 'arTeMiDe.',moduleName,': the values are to be restored from the initial grid'
                end if
        
            else !!! grid is not ready (how comes?)
                call MakeGrid()
                gridReady=.true.
            end if

        else  !! values are old

            if(gridReady) then  !!! grid is already build
                !!!nothing to do
            else!!rare option then parameters are not new but grit is not build
                if(outputLevel>2) write(*,*) 'arTeMiDe.',moduleName,': parameters are not reset. But grid is not ready.'
                call MakeGrid()
                gridReady=.true.
            end if
        
        end if

    else
        gridReady=.false.
    end if 

end subroutine SiversTMDPDF_SetLambdaNP_usual

!!! returns current value of NP parameters
subroutine SiversTMDPDF_CurrentNPparameters(var)
    real(dp),dimension(1:lambdaNPlength),intent(out)::var
    var=lambdaNP
end subroutine SiversTMDPDF_CurrentNPparameters

!!! This subroutine ask for the grid reconstruction (or destruction)
subroutine SiversTMDPDF_resetGrid(buildGrid,gluonRequired)
    logical,optional,intent(in)::buildGrid,gluonRequired
    logical::previousState

    if(present(buildGrid)) prepareGrid=buildGrid
    if(present(gluonRequired)) withGluon=gluonRequired

    previousState=gridReady
    gridReady=.false.

    !! we recalculate grid only if it was already calculated!
    if(prepareGrid .and. previousState) then
        if(outputLevel>1) write(*,*) 'arTeMiDe ',moduleName,':  Grid Reset. with c4=',c4_global
        call MakeGrid()
        gridReady=.true.
    end if
end subroutine SiversTMDPDF_resetGrid


!!!! this routine set the variations of scales
!!!! it is used for the estimation of errors
subroutine SiversTMDPDF_SetScaleVariation(c4_in)
    real(dp),intent(in)::c4_in
    if(c4_in<0.1d0 .or. c4_in>10.d0) then
        if(outputLevel>0) write(*,*) WarningString('variation in c4 is enourmous. c4 is set to 2',moduleName)
        c4_global=2d0
        call SiversTMDPDF_resetGrid()
    else if(abs(c4_in-c4_global)<tolerance) then
        if(outputLevel>1) write(*,*) color('SiversTMDPDF: c4-variation is ignored. c4='//real8ToStr(c4_global),c_yellow)
    else
        c4_global=c4_in
        if(outputLevel>1) write(*,*) color('SiversTMDPDF: set scale variations constant c4 as:'//real8ToStr(c4_global),c_yellow)
        call SiversTMDPDF_resetGrid()
    end if
end subroutine SiversTMDPDF_SetScaleVariation

end module SiversTMDPDF
