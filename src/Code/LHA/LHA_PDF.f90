!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           arTeMiDe 3.01
!
! The module that handle the PDF from the LHAPDF-file. 1412.7420
! Reads info-file, interpolate, etc.
! This module is contains the definition of class
!
!           A.Vladimirov (16.10.2025)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module LHA_PDF
use aTMDe_Numerics
use aTMDe_IO
implicit none

!!!!!! DEBUGMODE 1 switches many messages
#define DEBUGMODE 0

private
!Current version of module
character (len=5),parameter :: version="v3.04"

type, public :: LHAPDFgridReader
    private
    character(:),allocatable :: moduleName
    integer::outputLevel=2
    !!!!! global variables to
    character(len=:), allocatable::ReferenceString !!! contains description of reference
    character(len=:), allocatable::PDFname          !!! name of PDF for messages
    character(len=:), allocatable::MainPath         !!! the path to the PDF-set without ".000?"
    integer:: NumMembers
    integer,allocatable,dimension(:):: FlavorArray
    real(dp)::XMin,XMax,QMin,QMax
    !!!!! from [1412.7420]
    !0, 1, or 2 to, respectively, indicate no forcing, forcing negative values to 0, or forcing negative-or-zero values to a very small positive constan
    integer::ForcePositive_Global=0

    real(dp)::tolerance=10.d-8
    real(dp)::tolerancePDF=10.d-12   !!!! tolerance in interpretation of PDF values

    !!!!!! Main grid (x,Q,f)
    real(dp),dimension(:,:,:),allocatable::MainGrid
    real(dp),dimension(:),allocatable::Xnodes,Qnodes,PDF_logQs,PDF_logXs
    real(dp),allocatable,dimension(:,:):: PDF_WX,PDF_WQ!! barycentric wieghts
    integer,allocatable,dimension(:,:):: PDF_indices!! indices

    real(dp),dimension(:,:),allocatable::extrapolationGAMMA !!! parameter of extrapolation

    logical::TablesAreReady=.false. !!! flag that shows that module is ready

contains
    procedure,public::xPDF=>xPDF_this
    procedure,public::SetReplica=>SetReplica_this
end type

interface LHAPDFgridReader
    procedure :: ReadInfo
end interface LHAPDFgridReader


contains

!!! split a line by ":"
!!! attempt to identify the entry, and parse it
subroutine ParseInfoLine(this,line)
type(LHAPDFgridReader),intent(inout)::this
character(len=*), intent(in)::line

character(len=:), allocatable::linePart1, linePart2, linePart3
integer::position_i,position_j
integer::i,j

this%ForcePositive_Global=0

position_i=SCAN(line,":")
!!! First part of the line contrain the name
linePart1=trim(line(:position_i-1))
!!! second part of the line contains the value
linePart2=trim(line(position_i+1:))

!!!! some times there are leading space-terms remove them
if(linePart1(1:1)==" ") linePart1=linePart1(2:)


SELECT CASE(linePart1)

CASE("SetDesc")
    !!!!! Set Description
    if(this%outputLevel>2) then
        write(*,*) "Description of LHAPDF set:"
        write(*,*) linePart2
    end if

CASE("SetIndex","Authors","FlavorScheme","OrderQCD","ErrorType","DataVersion", &
    "Format", "NumFlavors", "AlphaS_OrderQCD", "MZ","MUp","MDown","MStrange",  &
    "MCharm", "MBottom", "MTop", "AlphaS_MZ", "AlphaS_Type","AlphaS_Qs", "AlphaS_Vals", &
    "Particle","AlphaS_Lambda3","AlphaS_Lambda4","AlphaS_Lambda5","ErrorConfLevel",&
    "ThresholdUP","ThresholdDown","ThresholdStrange","ThresholdCharm","ThresholdBottom","ThresholdTop")
    !!!!! Unused values of LHA-PDF file
    if(this%outputLevel>3) then
        write(*,*) "Value of "//trim(linePart1)//" in LHAPDF-info :"
        write(*,*) linePart2
    end if

CASE("Reference")
    this%ReferenceString=linePart2
    if(this%outputLevel>2) then
        write(*,*) "Reference :",this%ReferenceString
    end if
CASE("NumMembers")
    read(linePart2,*) this%NumMembers
    if(this%outputLevel>2) then
        write(*,*) "Number of Members of PDF-set :", this%NumMembers
    end if
CASE("Flavors")
    !!!! flavors have structure [,,,]
    position_i=SCAN(linePart2,"[")
    position_j=SCAN(linePart2,"]")
    linePart3=linePart2(position_i+1:position_j-1)
    !!! determine lengh of array by counting ","
    j=0
    do i=0,len(linePart3)
        !write(*,*)linePart3(i:i)
        if(trim(linePart3(i:i)) .eq. ",") j=j+1
    end do
    if(allocated(this%FlavorArray) .and. this%outputLevel>1) then
    write(*,*) WarningString(&
        "Double entry of the flavor array in the info-file of "//trim(this%PDFname)//". Second entry ignored.",this%moduleName)
    else
        allocate(this%FlavorArray(1:j+1))
        read(linePart3,*) this%FlavorArray
        if(this%outputLevel>2) then
            write(*,'(A, I4, 999(", ",I4))') "Flavor-enumeration Array :", this%FlavorArray
        end if
    end if

CASE("XMin")
    read(linePart2,*) this%XMin
    if(this%outputLevel>2) then
        write(*,'(A,F10.8)') " XMin     :", this%XMin
    end if
CASE("XMax")
    read(linePart2,*) this%XMax
    if(this%outputLevel>2) then
        write(*,'(A,F10.8)') " XMax     :", this%XMax
    end if
CASE("QMin")
    read(linePart2,*) this%QMin
    if(this%outputLevel>2) then
        write(*,'(A,F10.2)') " QMin     :", this%QMin
    end if
CASE("QMax")
    read(linePart2,*) this%QMax
    if(this%outputLevel>2) then
        write(*,'(A,F10.2)') " QMax     :", this%QMax
    end if

CASE("ForcePositive")
    read(linePart2,*) this%ForcePositive_Global
    if(this%outputLevel>2) then
        write(*,'(A,F10.2)') " ForcePositive :", this%ForcePositive_Global
    end if

CASE DEFAULT
    if(this%outputLevel>1) then
    write(*,*) color("Unknown entry in LHAPDF-info file :"//linePart1,c_red)//" ... ignore"
    end if

END SELECT
    
end subroutine ParseInfoLine

!!! opens the info-file and read lines.
function ReadInfo(name,directory,outL) result(this)
    type(LHAPDFgridReader)::this
    character(len=*),intent(in)::name
    character(len=*),intent(in)::directory
    integer,intent(in)::outL
    character(len=300)::path
    character(len=1024)::line !!!! 1024 line size!
    integer::ios,i,j
    
    this%TablesAreReady=.false.
    this%outputLevel=outL

    this%PDFname=trim(name)
    if(len(this%PDFname)>9) then
        this%moduleName=this%PDFname(1:9)//".reader"
    else
        this%moduleName=trim(this%PDFname)//".reader"
    end if

    this%MainPath=trim(adjustl(directory))//trim(adjustr(name))//"/"//trim(adjustr(name))
    path=trim(this%MainPath)//".info"

    if(this%outputLevel>1) write(*,'(A)') color("----- Loading PDF from "//trim(name),c_yellow)

    OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old", IOSTAT=ios)
    if(ios /= 0) ERROR STOP ErrorString('The info-file is not found at '//trim(path),this%moduleName)
    
    !!! reading all lines
    do
        read(51,'(A)', iostat=ios) line
        if(ios /= 0) exit !!! end of file
        
        !!!! there could be empty lines.
        if(len(trim(line))>0) call ParseInfoLine(this,line)
    end do
    
    CLOSE (51, STATUS='KEEP')

    if(this%outputLevel>1) write(*,*) "Info-file for "//trim(this%PDFname)//" parsed succesfully"

    call this%SetReplica(0)

    if(this%outputLevel>1) write(*,'(A)') color("----- PDF "//trim(this%PDFname)//" loaded. ",c_yellow)&
        //color("["//trim(this%ReferenceString)//"]",c_blue)
    
end function ReadInfo

!!! read the table for PDF-replica, and prepare all the Tables
!!! requires only the num. The PDF is searched in the MainPath, which is set by info-file
subroutine SetReplica_this(this,num)
    class(LHAPDFgridReader),intent(inout)::this
    integer,intent(in)::num !!! the number of replica to set
    character(len=500)::path
    character(len=4096)::line !!!! 1024 line size!
    character(len=:),allocatable:: lineToParse
    real(dp),dimension(:),allocatable::listToCheck,PDFentry,Qentry
    integer,dimension(:),allocatable::listOfF
    integer,dimension(-5:5)::FlavorPermutationIndex !!! list that contain the enumeration change variable
    logical::XisSet
    integer::Xsize,Qsize,QsizeSUBGRID
    integer::ios,i,j,k,iX,iQ,iF,indexQ,dummyI

    this%TablesAreReady=.false.

    if(allocated(this%Xnodes)) deallocate(this%Xnodes)
    if(allocated(this%Qnodes)) deallocate(this%Qnodes)
    if(allocated(this%PDF_logQs)) deallocate(this%PDF_logQs)
    if(allocated(this%PDF_logXs)) deallocate(this%PDF_logXs)
    if(allocated(this%MainGrid)) deallocate(this%MainGrid)
    if(allocated(this%PDF_WQ)) deallocate(this%PDF_WQ)
    if(allocated(this%PDF_WX)) deallocate(this%PDF_WX)
    if(allocated(this%PDF_indices)) deallocate(this%PDF_indices)
    if(allocated(this%extrapolationGAMMA)) deallocate(this%extrapolationGAMMA)

    XisSet=.false.!!!! trigger fro the first set of X's
    Qsize=0
    allocate(listOfF(0:size(this%FlavorArray)-1))

    if(num<0 .or. num>this%NumMembers) &
        ERROR STOP ErrorString('Attempt to call for replica '//trim(numToStr(num))//&
            '. The maximum number of replicas (according to info) is '//trim(numToStr(this%NumMembers)),this%moduleName)

    !!!! create the name of replica-file MainPath_000n.dat
    write(path,'(A,"_",I4.4,".dat")') trim(this%MainPath),num

#if DEBUGMODE==1
    write(*,*) "------------ LHA-READER DEBUG MODE ----------------"
    write(*,*) "-------------- "//trim(this%PDFname)//" ----------------"
#endif


    OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old", IOSTAT=ios)
    if(ios /= 0) ERROR STOP ErrorString('The file is not found at '//trim(path),this%moduleName)

    !!!! the file should be structure as follows [1412.7420]
    !!!PdfType: centra
    !!!Format: lhagrid1
    !!!---
    !!! List of X
    !!! List of Q
    !!! List of flavors
    !!! values per flavor for x1,Q1
    !!! values per flavor for x1,Q2
    !!! ...
    !!! values per flavor for x1,Qn
    !!! values per flavor for x2,Q1
    !!! ...
    !!! ...
    !!! values per flavor for xn,Qn
    !!!---
    !!! List of X
    !!! List of Q
    !!! List of flavors
    !!! values per flavor for x1,Q1
    !!! etc.
    !!!---
    !!!

    !!! To read such structure I do two runs of reading

    !!!! the first run is the check of the consistency
    !!!! 1) all lists of flavors are the same
    !!!! 2) all lists of X are the same
    !!!! 3) Count the number of Q's

    !!! then I can allocate the table for PDF, and other tables
    !!! The second run
    !!! 1) Fill the list of Q's
    !!! 2) Fill the list of PDF's

#if DEBUGMODE==1
    write(*,*) "First read of the table :"
#endif
    !!! reading all lines (first run)
    do
        read(51,'(A)', iostat=ios) line
        if(ios /= 0) exit !!! end of file
        !write(*,*) "---->>>   ",trim(line)

        !!!! search for "---"
        if(trim(line)=="---") then
            read(51,'(A)', iostat=ios) line
            if(ios /= 0) exit !!! end of file

            !!!!!!!!-------------This should be the line of X's
            !!!!!!!! parse and check the correctness
            lineToParse=trim(line)

            !!! the list is unformated (!!!) so I count the number of ".", assuming that each number has some decimals.
            j=0 !! (starts from 0 because its number of ".")
            do i=1,len(lineToParse)
                if(lineToParse(i:i) .eq. ".") j=j+1
            end do

#if DEBUGMODE==1
            write(*,*) 'Subtable: counted number of Xs--->',j
#endif

            if(XisSet) then
                if(allocated(listToCheck)) deallocate(listToCheck)
                allocate(listToCheck(0:j-1)) !!! number of "."=number of terms
                read(lineToParse,*) listToCheck
                if(size(listToCheck)/=Xsize) then
                    CLOSE (51, STATUS='KEEP')
                    ERROR STOP ErrorString("Different sizes of x-grids in LHAPDF file.",this%moduleName)
                end if
                if(sum(abs(this%Xnodes/listToCheck)-1)>this%tolerance) then
                    CLOSE (51, STATUS='KEEP')
                    ERROR STOP ErrorString("X-grids for subgrids in LHAPDF file do not coincide.",this%moduleName)
                end if
            else
                XisSet=.true.
                !!! set up the first list of X's
                allocate(this%Xnodes(0:j-1)) !!! number of "."=number of terms
                Xsize=j
                read(lineToParse,*) this%Xnodes
            end if

            !!!!!!!!-------------This should be the line of Q's
            !!!!!!!! count number of entries
            read(51,'(A)', iostat=ios) line
            lineToParse=trim(line)
            !!! the list is space separated (count the number of spaces)
            j=0 !! (starts from 0 because its number of ".")
            do i=1,len(lineToParse)
                if(lineToParse(i:i) .eq. ".") j=j+1
            end do
            QsizeSUBGRID=j
            Qsize=Qsize+QsizeSUBGRID

#if DEBUGMODE==1
            write(*,*) 'Subtable: counted number of Qs--->',QsizeSUBGRID
            write(*,*) 'Current total Q               --->',Qsize
#endif

            !!!!!!!!-------------This should be the line of flavors
            !!!!!!!! check that it coincides with the line of flavors from the info-file
            read(51,'(A)', iostat=ios) line
            lineToParse=trim(line)
            read(lineToParse,*) listOfF
            if(sum(listOfF-this%FlavorArray)/=0) then
                CLOSE (51, STATUS='KEEP')
                ERROR STOP ErrorString("Flavor-grids for subgrids in LHAPDF file do not coincide.",this%moduleName)
            end if

            !!! then it should be j * Xsize lines of grids
            do i=1, QsizeSUBGRID*Xsize
                read(51,'(A)', iostat=ios) line
                if(ios /= 0) then
                    CLOSE (51, STATUS='KEEP')
                    ERROR STOP ErrorString("Unexpected END-OF-FILE in "//trim(path),this%moduleName)
                end if

                if(trim(line)=="---") then
                    CLOSE (51, STATUS='KEEP')
                    ERROR STOP ErrorString("Incorrect number of entries subgrid of "//trim(path),this%moduleName)
                end if
            end do
        else
#if DEBUGMODE==1
            write(*,*) "--->>> EXPECTED --- LINE <<<---"
#endif
        end if
    end do

    CLOSE (51, STATUS='KEEP')

#if DEBUGMODE==1
    write(*,*) "First read complete"
    write(*,*) "X size: ( 0 :",Xsize-1,")"
    write(*,*) "Q size: ( 0 :",Qsize-1,")"
    write(*,*) "FlavorArray", this%FlavorArray
#endif

    !!!!! now the Qsize, Xsize, and Xnodes are set up

    !!! Allocate varriables
    allocate(this%Qnodes(0:Qsize-1))
    allocate(this%MainGrid(0:Xsize-1,0:Qsize-1,-5:5))
    allocate(PDFentry(1:size(this%FlavorArray)))

    !!!! the list of parcing the flavor to (-5:5)
    k=0
    do i=-5,5
        do j=1,size(this%FlavorArray)
            if(i==0 .and. this%FlavorArray(j)==21) then
                FlavorPermutationIndex(i)=j
                k=k+1
                exit
            !else if(i/=0 .and. i>=-5 .and. i<=5 .and. FlavorArray(j)==i) then
            else if(this%FlavorArray(j)==i) then
                FlavorPermutationIndex(i)=j
                k=k+1
                exit
            else !!!!! if there was no entry, I set the corresponding index to -50
                FlavorPermutationIndex(i)=-50
            end if
        end do
    end do

#if DEBUGMODE==1
    write(*,*) "FlavorPermutationIndex", FlavorPermutationIndex
#endif

    if(k/=11) then
        write(*,*) WarningString(&
            "The flavor numbering list in "//trim(this%PDFname)//" does not map to standard (-5:5).",this%moduleName)
#if DEBUGMODE==1
        write(*,*) "The interpretaion of permutation is the following (-50=empty)"
        write(*,*) "Entry :", this%FlavorArray
        write(*,*) "Interpretation :", FlavorPermutationIndex
#endif
        !!ERROR STOP ErrorString("The flavor numbering list does not map to (-5:5)",this%moduleName)
    end if

#if DEBUGMODE==1
    write(*,*) " -- "
    write(*,*) "Second read of the table :"
#endif

    !!! now read the second time and write all grids
    !k=1!!!! counter of PDF-entry
    !j=1!!!! counter of Q-entry
    indexQ=0 !!!! index of current position in Qnodes

    OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old", IOSTAT=ios)
    do
        read(51,'(A)', iostat=ios) line
        if(ios /= 0) exit !!! end of file
        !write(*,*) "---->>>   ",trim(line)

        !!!! search for "---"
        if(trim(line)=="---") then
            read(51,'(A)', iostat=ios) line !!!! this is X-line (skip) or EOF
            !write(*,*) "---->>>",trim(line)
            if(ios /= 0) exit !!! end of file

#if DEBUGMODE==1
            write(*,*) "Enter block :",indexQ
#endif

            !!!!!!!!-------------This should be the line of Q's
            !!!!!!!! re-read entries and add to the Qline
            !!!!!!!! count number of entries
            read(51,'(A)', iostat=ios) line
            lineToParse=trim(line)
            !!! the list is space separated (count the number of spaces)
            j=0 !! (starts from 0 because its number of ".")
            do i=1,len(lineToParse)
                if(lineToParse(i:i) .eq. ".") j=j+1
            end do
            QsizeSUBGRID=j

            allocate(Qentry(0:QsizeSUBGRID-1))
            read(lineToParse,*) Qentry
            !!! add it to the list of Qnodes
            this%Qnodes(indexQ:indexQ+QsizeSUBGRID-1)=Qentry(0:QsizeSUBGRID-1)

#if DEBUGMODE==1
            write(*,*) 'Subtable: counted number of Qs--->',QsizeSUBGRID
            write(*,*) 'Q nodes from ',indexQ, " to ", indexQ+QsizeSUBGRID-1 ," set."
#endif

            !!!!!!!------------This is line of flavors (skip)
            read(51,'(A)', iostat=ios) line

            !!!!!! reading the PDF-entries
            !!!!!! and insert it into the corresponding Grid-elements
            dummyI=0 !!!! counter of lines
            do i=0,Xsize-1
            do k=0,QsizeSUBGRID-1
                read(51,*) PDFentry
                !!!write(*,*) i,indexQ+k
                do j=-5,5
                    if(FlavorPermutationIndex(j)==-50) then !!!! absent entry replaced by 0
                        this%MainGrid(i,indexQ+k,j)=0._dp
                    else
                        this%MainGrid(i,indexQ+k,j)=PDFentry(FlavorPermutationIndex(j))
                    end if
                end do
                dummyI=dummyI+1
            end do
            end do
#if DEBUGMODE==1
            write(*,*) 'Have read ', dummyI,' lines.'
#endif

            deallocate(Qentry)
            indexQ=indexQ+QsizeSUBGRID  !!!! not forget to update indexQ
        else
#if DEBUGMODE==1
            write(*,*) "--->>> EXPECTED --- LINE <<<---"
#endif
        end if

    end do
    CLOSE (51, STATUS='KEEP')

#if DEBUGMODE==1
    write(*,*) "Second read complete."
#endif

    !!!!! checking the consistency of the read table
    !!!!! and implicating the ForcePositive statement
    do i=0,Xsize-1
    do k=0,Qsize-1
        do j=-5,5
            if(this%ForcePositive_Global==1 .and. this%MainGrid(i,k,j)<0) this%MainGrid(i,k,j)=0._dp
            if(this%ForcePositive_Global==2 .and. this%MainGrid(i,k,j)<=this%tolerancePDF) this%MainGrid(i,k,j)=this%tolerancePDF

            if(ISNAN(this%MainGrid(i,k,j)) .or. this%MainGrid(i,k,j)>10d5) then
                write(*,*) this%PDFname," SOMETHING STRANGE IN THE PDF ENTRY (x,Q,f) ",i,k,j
                write(*,*) this%MainGrid(i,k,:)
              error stop ErrorString("An entry in LHA tables has been interpreted as NaN or >10^5.",this%moduleName)
            end if
        end do
    end do
    end do

    !!!!!------------ prepare the table of interpolation
    allocate(this%PDF_logQs(0:Qsize-1))
    allocate(this%PDF_logXs(0:Xsize-1))
    allocate(this%PDF_indices(1:4,0:Qsize-2))
    allocate(this%PDF_WQ(1:4,0:Qsize-2))
    allocate(this%PDF_WX(1:4,0:Xsize-2))

    this%PDF_LogQs=log(this%Qnodes)
    this%PDF_LogXs=log(this%Xnodes)

    !!! compute the index table for cubic-interpolation
    !!! if the x belongs to (x_{j-1}<x<x_j) the interpolation is done by indices{j}
    do i=0,size(this%PDF_LogQs)-2
        if(i==0) then
            this%PDF_indices(1:4,i)=(/0,1,2,3/)
        else if(i==size(this%PDF_LogQs)-2) then
            this%PDF_indices(1:4,i)=i-1-(/1,0,-1,-2/)
        else
            this%PDF_indices(1:4,i)=i-(/1,0,-1,-2/)
            !!! if the value of alpha is repeated these are ends of subgrids.
            !!! I shift the indicing
            if(abs(this%PDF_LogQs(i+1)-this%PDF_LogQs(i+2))<10d-4) then
                this%PDF_indices(1:4,i)=this%PDF_indices(1:4,i)-1
            else if(abs(this%PDF_LogQs(i)-this%PDF_LogQs(i+1))<10d-4) then
                !!! zero-size bin
                !!! this should be never reached. but I assign it to previous one
                this%PDF_indices(1:4,i)=this%PDF_indices(1:4,i)-2
            else if(abs(this%PDF_LogQs(i-1)-this%PDF_LogQs(i))<10d-4) then
                this%PDF_indices(1:4,i)=this%PDF_indices(1:4,i)+1
            end if
        end if
    end do

#if DEBUGMODE==1
        write(*,*) " "
        write(*,*) "----- TABLE FOR INDICES OF Q-INTERPOLATION -----"
        do i = 0, Qsize-2
            write(*,*) i, " --> ",this%PDF_indices(1:4,i)
        end do
        write(*,*) " "
#endif

    !!!! table for X is very simple,  there is no repeating entries.
    !!!! so index is i+(-1,0,1,2), except for i=0 -> (0,1,2,3), and i=S-1 (S-3,S-2,S-1,S)

    !!! precompute barycentric weights for Q
    do i=0,size(this%PDF_LogQs)-2
        this%PDF_WQ(1,i)=1.d0/(&
        (this%PDF_LogQs(this%PDF_indices(1,i))-this%PDF_LogQs(this%PDF_indices(2,i)))*&
        (this%PDF_LogQs(this%PDF_indices(1,i))-this%PDF_LogQs(this%PDF_indices(3,i)))*&
        (this%PDF_LogQs(this%PDF_indices(1,i))-this%PDF_LogQs(this%PDF_indices(4,i))))

        this%PDF_WQ(2,i)=1.d0/(&
        (this%PDF_LogQs(this%PDF_indices(2,i))-this%PDF_LogQs(this%PDF_indices(1,i)))*&
        (this%PDF_LogQs(this%PDF_indices(2,i))-this%PDF_LogQs(this%PDF_indices(3,i)))*&
        (this%PDF_LogQs(this%PDF_indices(2,i))-this%PDF_LogQs(this%PDF_indices(4,i))))

        this%PDF_WQ(3,i)=1.d0/(&
        (this%PDF_LogQs(this%PDF_indices(3,i))-this%PDF_LogQs(this%PDF_indices(1,i)))*&
        (this%PDF_LogQs(this%PDF_indices(3,i))-this%PDF_LogQs(this%PDF_indices(2,i)))*&
        (this%PDF_LogQs(this%PDF_indices(3,i))-this%PDF_LogQs(this%PDF_indices(4,i))))

        this%PDF_WQ(4,i)=1.d0/(&
        (this%PDF_LogQs(this%PDF_indices(4,i))-this%PDF_LogQs(this%PDF_indices(1,i)))*&
        (this%PDF_LogQs(this%PDF_indices(4,i))-this%PDF_LogQs(this%PDF_indices(2,i)))*&
        (this%PDF_LogQs(this%PDF_indices(4,i))-this%PDF_LogQs(this%PDF_indices(3,i))))
    end do

    !!! precompute barycentric weights for X
    do i=0,size(this%PDF_LogXs)-2
        if(i==0) then
            k=1
        else if(i==size(this%PDF_LogXs)-2) then
            k=i-1
        else
            k=i
        end if

        this%PDF_WX(1,i)=1.d0/(&
        (this%PDF_LogXs(k-1)-this%PDF_LogXs(k))*(this%PDF_LogXs(k-1)-this%PDF_LogXs(k+1))*(this%PDF_LogXs(k-1)-this%PDF_LogXs(k+2)))

        this%PDF_WX(2,i)=1.d0/(&
        (this%PDF_LogXs(k)-this%PDF_LogXs(k-1))*(this%PDF_LogXs(k)-this%PDF_LogXs(k+1))*(this%PDF_LogXs(k)-this%PDF_LogXs(k+2)))

        this%PDF_WX(3,i)=1.d0/(&
        (this%PDF_LogXs(k+1)-this%PDF_LogXs(k-1))*(this%PDF_LogXs(k+1)-this%PDF_LogXs(k))*(this%PDF_LogXs(k+1)-this%PDF_LogXs(k+2)))

        this%PDF_WX(4,i)=1.d0/(&
        (this%PDF_LogXs(k+2)-this%PDF_LogXs(k-1))*(this%PDF_LogXs(k+2)-this%PDF_LogXs(k))*(this%PDF_LogXs(k+2)-this%PDF_LogXs(k+1)))
    end do

    !!! extrapolation is done by formula F(x)*(Q/Q0)^[gamma(x) Q +2]
    !!! where gamma(x)=(gg-2)/Q0, where gg=D[logF]/dlogQ
    !!! see (4,5) in 1412.7420
    allocate(this%extrapolationGAMMA(0:Xsize-1,-5:5))
    !!!!! it is possible that the F(x)=0, in this case there is nothing to extrapolate, set gamma=0
    !!!!! since it is an exceptional case, one should check for it one-by-one
    do i=0,Xsize-1
    do j=-5,5

        if(this%MainGrid(i,0,j)*this%MainGrid(i,1,j)>this%tolerancePDF) then
        !!! this is derivative (the ratio is important to handle the both negative signs (YES, it happen sometimes)
            this%extrapolationGAMMA(i,j)=(log(this%MainGrid(i,1,j)/this%MainGrid(i,0,j)))/(this%PDF_LogQs(1)-this%PDF_LogQs(0))
            this%extrapolationGAMMA(i,j)=(this%extrapolationGAMMA(i,j)-2)/sqrt(this%Qnodes(0))
        else if(this%MainGrid(i,0,j)*this%MainGrid(i,1,j)<-this%tolerancePDF) then
            if(this%outputLevel>2) then
            write(*,*) WarningString("Low-Q extrapolation formula fails, because PDF cross the zero-value.",this%moduleName)
            write(*,*) "Extrapolation parameters: x-node =",i," flavor=",j
            write(*,*) "Nodes of to extrapolate [f(Q0),f(Q1)]: = [",this%MainGrid(i,0,j)," ,",this%MainGrid(i,1,j),"]"
            write(*,*) WarningString("Extrapolation parameter set to zero.",this%moduleName)
            end if
            this%extrapolationGAMMA(i,j)=0.d0
        else
            this%extrapolationGAMMA(i,j)=0.d0
        end if
    end do
    end do


    do i=0,Xsize-1
    do j=-5,5
            if(ISNAN(this%extrapolationGAMMA(i,j)) .or. this%extrapolationGAMMA(i,j)>10d5) then
                write(*,*) this%PDFname," SOMETHING STRANGE IN THE PDF Extrapolation part (x,f) ",i,j
                write(*,*) this%extrapolationGAMMA(i,:)
                write(*,*) "the points to extrapolate", this%MainGrid(i,0,j), this%MainGrid(i,1,j)
              error stop ErrorString("An extrapolation parameter has been computed NaN or >10^5.",this%moduleName)
            end if
    end do
    end do

    if(this%outputLevel>1) write(*,"(A,I4,A)") " "//trim(this%moduleName)//": replica ",num," of "//trim(this%PDFname)//" loaded."

end subroutine SetReplica_this

!!! returns the value of xPDF extrapolated from the grid
#if DEBUGMODE==1
function xPDF_this(this,x,Q)
#else
pure function xPDF_this(this,x,Q)
#endif
class(LHAPDFgridReader),intent(in)::this
real(dp),intent(in)::x,Q
real(dp),dimension(-5:5):: xPDF_this

real(dp)::logQ,dd,deltaQ(1:4),subQ(1:4,-5:5),logX,deltaX(1:4),gg(-5:5)
integer::iQ,iX,j
logical::flag

if(X<this%Xmin .or. X>this%Xmax) ERROR STOP ErrorString('X ='//numTostr(X)//' is outside of X-range',this%moduleName)

if(Q>this%Qmax) ERROR STOP ErrorString('Q ='//numTostr(this%Qmax)//' is outside of QMax-range',this%moduleName)
if(Q<0.4d0) ERROR STOP ErrorString('Q ='//numTostr(this%Qmax)//' is smaller that 0.4 GeV of Q-range',this%moduleName)

if(Q<this%Qmin) then !!!! extrapolate!
    logX=log(X)

    !!! search for the index of X
    if(logX<this%PDF_LogXs(1)) then
        iX=1
    else if(logX>this%PDF_LogXs(size(this%PDF_LogXs)-2)) then
        iX=size(this%PDF_LogXs)-3
    else
        do iX=1,size(this%PDF_LogXs)-3
            if(logX<this%PDF_LogXs(iX+1)) exit
        end do
    end if

    !!!! first interpolate gamma to the x
    subQ(1:4,-5:5)=this%extrapolationGAMMA(iX-1:iX+2,-5:5)
    do j=1,4
    dd=logX-this%PDF_LogXs(iX+j-2)
    !!! I check if the point is close to the node
    if(abs(dd)<this%tolerance) then
        gg(-5:5)=subQ(j,-5:5)
        goto 21
    end if
    deltaX(j)=this%PDF_WX(j,iX)/dd
    end do

    gg=(deltaX(1)*subQ(1,-5:5)+deltaX(2)*subQ(2,-5:5)+deltaX(3)*subQ(3,-5:5)+deltaX(4)*subQ(4,-5:5))/sum(deltaX)

    !!!! second interpolate PDF at Q0 to the x
21  subQ(1:4,-5:5)=this%MainGrid(iX-1:iX+2,0,-5:5)
    do j=1,4
    dd=logX-this%PDF_LogXs(iX+j-2)
    !!! I check if the point is close to the node
    if(abs(dd)<this%tolerance) then
        xPDF_this(-5:5)=subQ(j,-5:5)
        goto 22
    end if
    deltaX(j)=this%PDF_WX(j,iX)/dd
    end do

    xPDF_this(-5:5)=(deltaX(1)*subQ(1,-5:5)+deltaX(2)*subQ(2,-5:5)+deltaX(3)*subQ(3,-5:5)+deltaX(4)*subQ(4,-5:5))/sum(deltaX)

22  xPDF_this(-5:5)=xPDF_this(-5:5)*(Q/this%Qnodes(0))**(gg(-5:5)*sqrt(Q) +2)

    return
end if

logQ=log(Q)
logX=log(X)

!!! search for the index of Q
do iQ=0,size(this%PDF_LogQs)-2
     if(logQ<this%PDF_LogQs(iQ+1)) exit
end do
!!! search for the index of X
if(logX<this%PDF_LogXs(1)) then
    iX=1
else if(logX>this%PDF_LogXs(size(this%PDF_LogXs)-2)) then
    iX=size(this%PDF_LogXs)-3
else
    do iX=1,size(this%PDF_LogXs)-3
        if(logX<this%PDF_LogXs(iX+1)) exit
    end do
end if

!!!! Q-interpolation
do j=1,4
    dd=logQ-this%PDF_LogQs(this%PDF_indices(j,iQ))
    !!! I check if the point is close to the node
    if(abs(dd)<this%tolerance) then
        subQ(1:4,-5:5)=this%MainGrid(iX-1:iX+2,this%PDF_indices(j,iQ),-5:5)
        goto 10
    end if
    deltaQ(j)=this%PDF_WQ(j,iQ)/dd
end do

subQ=(deltaQ(1)*this%MainGrid(iX-1:iX+2,this%PDF_indices(1,iQ),-5:5)+&
deltaQ(2)*this%MainGrid(iX-1:iX+2,this%PDF_indices(2,iQ),-5:5)+&
deltaQ(3)*this%MainGrid(iX-1:iX+2,this%PDF_indices(3,iQ),-5:5)+&
deltaQ(4)*this%MainGrid(iX-1:iX+2,this%PDF_indices(4,iQ),-5:5))/sum(deltaQ)

!!! X-interpolation
10 do j=1,4
    dd=logX-this%PDF_LogXs(iX+j-2)
    !!! I check if the point is close to the node
    if(abs(dd)<this%tolerance) then
        xPDF_this(-5:5)=subQ(j,-5:5)
        return
    end if
    deltaX(j)=this%PDF_WX(j,iX)/dd
end do

xPDF_this=(deltaX(1)*subQ(1,-5:5)+deltaX(2)*subQ(2,-5:5)+deltaX(3)*subQ(3,-5:5)+deltaX(4)*subQ(4,-5:5))/sum(deltaX)

#if DEBUGMODE==1
  do j=-5,5
   if(ISNAN(xPDF_this(j))) then

    write(*,*) ErrorString('a values of PDF '//trim(this%PDFname)//' computed to NAN',this%moduleName)
    write(*,*) '----- information on last call -----'
    write(*,*) 'x=', x, 'Q=',Q,' j=',j,' result=',xPDF_this(j)
    write(*,*) iX, iQ
    write(*,*) this%PDF_indices(:,iQ)
    write(*,*) " -- ", size(this%MainGrid(:,1,1)),size(this%MainGrid(1,:,1)),size(this%MainGrid(1,1,:))
    write(*,*) ' -------- part of grid that has been used ------  '
    write(*,*) "Q1->", this%MainGrid(iX-1:iX+2,this%PDF_indices(1,iQ),j)
    write(*,*) "Q2->", this%MainGrid(iX-1:iX+2,this%PDF_indices(2,iQ),j)
    write(*,*) "Q3->", this%MainGrid(iX-1:iX+2,this%PDF_indices(3,iQ),j)
    write(*,*) "Q4->", this%MainGrid(iX-1:iX+2,this%PDF_indices(4,iQ),j)
    error stop
   end if
  end do
#endif

end function xPDF_this

end module LHA_PDF
