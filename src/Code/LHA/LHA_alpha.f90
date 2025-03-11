!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!           arTeMiDe 3.01
!
! The module that reads the LHAPDF-info-file and extract the table for values of alpha_s
!  Interpolates and returns its value.
! This module is unique, because all artemide-submodules use the same alpha_s for consitency
! For a moment only "ipol" version of the alpha_s is possible.
!
!           A.Vladimirov
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module LHA_alpha
use IO_functions
implicit none

private
!Current version of module
character (len=5),parameter :: version="v3.01"
character (len=9),parameter :: moduleName="LHA_alpha"
integer::outputLevel=2

!!!!! global variables to 
character(len=:), allocatable::ReferenceString !!! contains description of reference
real(dp):: QMin,QMax
real(dp):: MZ, MUp, MDown, MStrange, MCharm, MBottom, MTop
real(dp):: AlphaS_MZ
real(dp),allocatable,dimension(:):: AlphaS_Qs,AlphaS_Vals,AlphaS_LogQs
real(dp),allocatable,dimension(:,:):: AlphaS_W!! barycentric wieghts
integer,allocatable,dimension(:,:):: AlphaS_indices!! indices

real(dp)::extrapolA,extrapolB1,extrapolB2,LambdaEFF !!! precomputed parameters of extrapolation

real(dp)::tolerance=10.d-6

logical::AlphaStype_IsRecognized=.false. !!!! TRUE= correct input for alphaS (used during load)
                                

public:: ReadInfo,AlphaS

contains

!!! split a line by ":"
!!! attempt to identify the entry, and parse it
subroutine ParseInfoLine(line)
    character(len=*), intent(in)::line
    character(len=:), allocatable::linePart1, linePart2, linePart3
    integer::position_i,position_j
    integer::i,j
    
    position_i=SCAN(line,":")
    !!! First part of the line contrain the name
    linePart1=trim(line(:position_i-1))
    !!! second part of the line contains the value
    linePart2=trim(line(position_i+1:))
    
    SELECT CASE(linePart1)
    
    CASE("SetDesc")
        !!!!! Set Description
        if(outputLevel>2) then
            write(*,*) "Description of LHAPDF set:"
            write(*,*) linePart2
        end if
        
    CASE("SetIndex","Authors","FlavorScheme","OrderQCD","ErrorType","DataVersion", &
        "Format", "NumFlavors", "AlphaS_OrderQCD", "NumMembers", "Flavors", "XMin", &
        "XMax","Particle","AlphaS_Lambda3","AlphaS_Lambda4","AlphaS_Lambda5",&
        "ThresholdUP","ThresholdDown","ThresholdStrange","ThresholdCharm","ThresholdBottom","ThresholdTop")
        !!!!! Unused values of LHA-PDF file
        if(outputLevel>3) then
            write(*,*) "Value of "//trim(linePart1)//" in LHAPDF-info :"
            write(*,*) linePart2
        end if
        
    CASE("Reference")
        ReferenceString=linePart2
        if(outputLevel>2) then
            write(*,*) "Reference :",ReferenceString
        end if

    CASE("QMin")
        read(linePart2,*) QMin
        if(outputLevel>2) then
            write(*,'(A,F10.2)') " QMin     :", QMin
        end if
    CASE("QMax")
        read(linePart2,*) QMax
        if(outputLevel>2) then
            write(*,'(A,F10.2)') " QMax     :", QMax
        end if
        
    CASE("MZ")
        read(linePart2,*) MZ
        if(outputLevel>2) then
            write(*,'(A,F10.6)') " MZ       :", MZ
        end if
    CASE("MUp")
        read(linePart2,*) MUp
        if(outputLevel>2) then
            write(*,'(A,F10.8)') " MUp      :", MUp
        end if
    CASE("MDown")
        read(linePart2,*) MDown 
        if(outputLevel>2) then
            write(*,'(A,F10.8)') " MDown    :", MDown
        end if
    CASE("MStrange")
        read(linePart2,*) MStrange
        if(outputLevel>2) then
            write(*,'(A,F10.6)') " MStrange :", MStrange
        end if
    CASE("MCharm")
        read(linePart2,*) MCharm
        if(outputLevel>2) then
            write(*,'(A,F10.6)') " MCharm   :", MCharm
        end if
    CASE("MBottom")
        read(linePart2,*) MBottom
        if(outputLevel>2) then
            write(*,'(A,F10.6)') " MBottom  :", MBottom
        end if
    CASE("MTop")
        read(linePart2,*) MTop
        if(outputLevel>2) then
            write(*,*) "MTop     :", MTop
        end if
    CASE("AlphaS_MZ")
        read(linePart2,*) AlphaS_MZ
        if(outputLevel>2) then
            write(*,'(A,F10.4)') " AlphaS at MZ :", AlphaS_MZ
        end if
    
    !!!! The reading of alpha_s part only if it is required
    CASE("AlphaS_Type")
        if(index(linePart2,"ipol")==0) then
            ERROR STOP ErrorString('Only "ipol" type of alphaS can be used. Current input is '//trim(linePart2),moduleName)
        else
            AlphaStype_IsRecognized=.true.!!! check passed
        end if
    
    CASE("AlphaS_Qs")
        !!!! list of Qs has the structure [,,,]
        position_i=SCAN(linePart2,"[")
        position_j=SCAN(linePart2,"]")
        linePart3=linePart2(position_i+1:position_j-1)
        !!! determine lengh of array by counting ","
        j=0
        do i=0,len(linePart3)
            !write(*,*)linePart3(i:i)
            if(trim(linePart3(i:i)) .eq. ",") j=j+1
        end do
        allocate(AlphaS_Qs(0:j))
        read(linePart3,*) AlphaS_Qs
        if(outputLevel>1) write(*,*) "AlphaS_Qs values parsed succesfully"
        
    CASE("AlphaS_Vals")
        !!!! list of vals has the structure [,,,]
        position_i=SCAN(linePart2,"[")
        position_j=SCAN(linePart2,"]")
        linePart3=linePart2(position_i+1:position_j-1)
        !!! determine lengh of array by counting ","
        j=0
        do i=0,len(linePart3)
            !write(*,*)linePart3(i:i)
            if(trim(linePart3(i:i)) .eq. ",") j=j+1
        end do
        allocate(AlphaS_Vals(0:j))
        read(linePart3,*) AlphaS_Vals
        if(outputLevel>1) write(*,*) "AlphaS_Vals values parsed succesfully"
    
    CASE DEFAULT
        if(outputLevel>1) then
        write(*,*) color("Unknown entry in LHAPDF-info file :"//linePart1,c_red)//" ... ignore"
        end if
    
    END SELECT
    
end subroutine ParseInfoLine

!!! opens the info-file and read lines.
subroutine ReadInfo(name,directory,outP)
    character(len=*),intent(in)::name
    character(len=*),intent(in)::directory
    integer,intent(in)::outP
    character(len=300)::path
    character(len=4096)::line !!!! Long line to guaranty the input
    integer::ios,i,j
    
    outputLevel=outP
    AlphaStype_IsRecognized=.false.
    
    path=trim(adjustl(directory))//trim(adjustr(name))//"/"//trim(adjustr(name))//".info"

    if(outputLevel>1) write(*,'(A)') color("----- Loading Alpha_s from "//trim(name),c_yellow)

    OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old", IOSTAT=ios)
    if(ios /= 0) ERROR STOP ErrorString('The info-file is not found at '//trim(path),moduleName)
    
    !!! reading all lines
    do
        read(51,'(A)', iostat=ios) line
        if(ios /= 0) exit !!! end of file
        
        !!!! there could be empty lines.
        if(len(trim(line))>0) call ParseInfoLine(line)
    end do
    
    CLOSE (51, STATUS='KEEP')
    
    
    !!! check that alphas-variables are loaded succesfully.
    if(.not.AlphaStype_IsRecognized) &
        ERROR STOP ErrorString('Type of alphaS is not recognized in  '//trim(path),moduleName)
    if(.not.allocated(AlphaS_Qs)) &
        ERROR STOP ErrorString('AlphaS_Qs variables are missed in '//trim(path),moduleName)
    if(.not.allocated(AlphaS_Vals)) &
        ERROR STOP ErrorString('AlphaS_Vals variables are missed in '//trim(path),moduleName)

    !!!!!------------ prepare the table of alphaS-interpolation
    allocate(AlphaS_LogQs(0:size(AlphaS_Qs)-1))
    allocate(AlphaS_indices(1:4,0:size(AlphaS_Qs)-2))
    allocate(AlphaS_W(1:4,0:size(AlphaS_Qs)-2))

    AlphaS_LogQs=log(AlphaS_Qs)

    !!! compute the index table for cubic-interpolation
    !!! if the x belongs to (x_{j-1}<x<x_j) the interpolation is done by indices{j}
    do i=0,size(AlphaS_Qs)-2
        if(i==0) then
            AlphaS_indices(1:4,i)=(/0,1,2,3/)
        else if(i==size(AlphaS_Qs)-2) then
            AlphaS_indices(1:4,i)=i-1-(/1,0,-1,-2/)
        else
            AlphaS_indices(1:4,i)=i-(/1,0,-1,-2/)
            !!! if the value of alpha is repeated these are ends of subgrids.
            !!! I shift the indicing
            if(abs(AlphaS_Qs(i+1)-AlphaS_Qs(i+2))<10d-4) then
                AlphaS_indices(1:4,i)=AlphaS_indices(1:4,i)-1
            else if(abs(AlphaS_Qs(i)-AlphaS_Qs(i+1))<10d-4) then
                !!! zero-size bin
                !!! this should be never reached. but I assign it to previous one
                AlphaS_indices(1:4,i)=AlphaS_indices(1:4,i)-2
            else if(abs(AlphaS_Qs(i-1)-AlphaS_Qs(i))<10d-4) then
                AlphaS_indices(1:4,i)=AlphaS_indices(1:4,i)+1
            end if
        end if
    end do
!
!         do i=0,size(AlphaS_Qs)-2
!         write(*,*) AlphaS_indices(1:4,i),AlphaS_Qs(i),i
!         end do
!         do i=0,size(AlphaS_Qs)-1
!             write(*,*) "{",AlphaS_Qs(i),",",AlphaS_Vals(i),"},"
!         end do

        !!! now compute the barycentric weights
    do i=0,size(AlphaS_Qs)-2
        AlphaS_W(1,i)=1.d0/(&
        (AlphaS_LogQs(AlphaS_indices(1,i))-AlphaS_LogQs(AlphaS_indices(2,i)))*&
        (AlphaS_LogQs(AlphaS_indices(1,i))-AlphaS_LogQs(AlphaS_indices(3,i)))*&
        (AlphaS_LogQs(AlphaS_indices(1,i))-AlphaS_LogQs(AlphaS_indices(4,i))))

        AlphaS_W(2,i)=1.d0/(&
        (AlphaS_LogQs(AlphaS_indices(2,i))-AlphaS_LogQs(AlphaS_indices(1,i)))*&
        (AlphaS_LogQs(AlphaS_indices(2,i))-AlphaS_LogQs(AlphaS_indices(3,i)))*&
        (AlphaS_LogQs(AlphaS_indices(2,i))-AlphaS_LogQs(AlphaS_indices(4,i))))

        AlphaS_W(3,i)=1.d0/(&
        (AlphaS_LogQs(AlphaS_indices(3,i))-AlphaS_LogQs(AlphaS_indices(1,i)))*&
        (AlphaS_LogQs(AlphaS_indices(3,i))-AlphaS_LogQs(AlphaS_indices(2,i)))*&
        (AlphaS_LogQs(AlphaS_indices(3,i))-AlphaS_LogQs(AlphaS_indices(4,i))))

        AlphaS_W(4,i)=1.d0/(&
        (AlphaS_LogQs(AlphaS_indices(4,i))-AlphaS_LogQs(AlphaS_indices(1,i)))*&
        (AlphaS_LogQs(AlphaS_indices(4,i))-AlphaS_LogQs(AlphaS_indices(2,i)))*&
        (AlphaS_LogQs(AlphaS_indices(4,i))-AlphaS_LogQs(AlphaS_indices(3,i))))
    end do

    !!!! linear extrapolation in the inverse alpha and log-scale
    !!!! alpha=A/(B1 log(Q/Q0)+B2)
    extrapolA=AlphaS_Vals(0)*AlphaS_Vals(1)*Log(AlphaS_Qs(1)/AlphaS_Qs(0))
    extrapolB1=(AlphaS_Vals(0)-AlphaS_Vals(1))
    extrapolB2=AlphaS_Vals(1)*Log(AlphaS_Qs(1)/AlphaS_Qs(0))

    LambdaEFF=exp(-extrapolB2/extrapolB1)*AlphaS_Qs(0)

    if(LambdaEFF+0.1>Qmin) &
    ERROR STOP ErrorString('Effective LambdaQCD computed as '//real8Tostr(LambdaEFF)//' It is too high...',moduleName)

    if(outputLevel>1) write(*,'("AlphaS prepared with Effective LambdaQCD = ",F10.6)') LambdaEFF
    LambdaEFF=max(LambdaEFF+0.1d0,0.4d0)

    if(outputLevel>1) write(*,'(A)') color("----- Alpha_s from "//trim(name)//" loaded. ",c_yellow)&
        //color("["//trim(ReferenceString)//"]",c_blue)

end subroutine ReadInfo


!!! returns the value of AlphaS interpolated from the table.
!!! the interpolation is log-linear.
function AlphaS(Q)
real(dp),intent(in)::Q
real(dp):: AlphaS

real(dp)::logQ,deltas(1:4)
integer::i,j
if(Q<Qmin) then !!! logarith log-extrapolation
    if(Q<LambdaEFF) ERROR STOP ErrorString('Q ='//real8Tostr(Q)//' is smaller than Effective LambdaQCD',moduleName)

    AlphaS=extrapolA/(extrapolB1*Log(Q/AlphaS_Qs(0))+extrapolB2)
else if(Q>Qmax) then !!! constant
    AlphaS=AlphaS_Vals(size(AlphaS_Qs)-1)
else
    logQ=log(Q)

    do i=0,size(AlphaS_Qs)-2
        if(logQ<AlphaS_LogQs(i+1)) exit
    end do

    ! !!! this is index of the Q-box
    do j=1,4
        deltas(j)=logQ-AlphaS_LogQs(AlphaS_indices(j,i))
        !!! I check if the point is close to the node
        if(abs(deltas(j))<tolerance) then
            AlphaS=AlphaS_Vals(AlphaS_indices(j,i))
            return
        end if
        deltas(j)=AlphaS_W(j,i)/deltas(j)
    end do

    AlphaS=sum(deltas*AlphaS_Vals(AlphaS_indices(1,i):AlphaS_indices(4,i)))/sum(deltas)
end if

end function AlphaS

end module LHA_alpha
