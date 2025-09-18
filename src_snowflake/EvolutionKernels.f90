!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                          !!
!! This module computes and saves the matrices of evolution kernels for the tw-3 evolution  !!
!!                                                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module EvolutionKernels
use IO_snowflake
use HexGrid
implicit none

private

INCLUDE 'commonVariables.f90'

!!!! GK_ORDER=1 = G7K15
!!!! GK_ORDER=1 = G20K41
#define GK_ORDER 2

#if GK_ORDER==1
INCLUDE 'Tables/G7K15.f90'
#elif GK_ORDER==2
INCLUDE 'Tables/G20K41.f90'
#endif

!!! QCD constants
!!!! color coefficient CA
real(dp),parameter::CA=3._dp
!!!! color coefficient CF
real(dp),parameter::CF=4._dp/3._dp
!!!! color coefficient CF-CA/2
real(dp),parameter::CFminusCA2=-1/6._dp
!!!! color coefficient ColorQG = - (n^2-4)/n
real(dp),parameter::ColorQG=-5/3._dp

!!!!! the sparse structure for convolution matrices
!!!!! it is the same for kernel matrices and for G2-convolution matrices
type SparseM2
    integer, allocatable :: i(:)  !!!! contains the list of 1D index which is in 1-to-1 asociation with the value
    real(dp), allocatable :: val(:) !!!! value for list of indices
end type SparseM2

!!!! Evolution kernels for the quark. Hqq_J is the mixture contribution
type(SparseM2),allocatable,dimension(:):: Hqq, Hqq_J, Hqq_CO
!!!! Evolution kernels for the quark/gluons. Plus and minus components
type(SparseM2),allocatable,dimension(:):: Hqg_PLUS, Hgq_PLUS, Hgg_PLUS
type(SparseM2),allocatable,dimension(:):: Hqg_MINUS, Hgq_MINUS, Hgg_MINUS

public:: EvolutionKernels_Initialize,EvSingletPLUS,EvSingletMINUS,EvNonSinglet,EvChiralOdd
public:: SaveKernels,ReadKernels

public:: H12hat,kernelHqq,TEST

!!!!!! matrix of G2 is usual matrix because it has all elements non-zero (except boundary, but there are a few of them)
public::G2_projM,G2xF,D2xF
real(dp),allocatable,dimension(:,:)::G2matrix
real(dp),allocatable,dimension(:)::D2matrix

!!!!!! matrix for WGT function
public::WGT_projM,WGTxF
real(dp),allocatable,dimension(:,:)::WGTmatrixP,WGTmatrixN


contains

INCLUDE 'ExpressionsForKernels.f90'
INCLUDE 'ExpressionsForG2.f90'
INCLUDE 'ExpressionsForD2.f90'
INCLUDE 'ExpressionsForWGT.f90'

subroutine EvolutionKernels_Initialize(path)
real*8::t1,t2,tt1,tt2
!$ real*8::omp_get_wtime
integer::sparsity
character(len=*)::path
logical::readKERNEL
character(len=300)::pathToKernels

!----------------- reading ini-file --------------------------------------
OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")
!!! Number of grid-nodes in the sector-angle
call MoveTO(51,'0   :')
read(51,*) showINI

!!! Minimal value of X
call MoveTO(51,'A.1 :')
read(51,*) xMIN

!!! Number of nodes in each subgrid in R
call MoveTO(51,'A.2 :')
read(51,*) NUM_RHO

!!! Number of nodes in each subgrid in PHI
call MoveTO(51,'A.3 :')
read(51,*) NUM_PHI

!!!! the last integer for grid in RHO and PHI, 0... N (i.e. total size is +1)
NUM_TOT_RHO=NUM_RHO
NUM_TOT_PHI=6*NUM_PHI-1

!!!! the last integer for 1D grid in RHO and PHI, 0... N (i.e. total size is +1)
NUM_TOT=(NUM_TOT_RHO+1)*(NUM_TOT_PHI+1)-1

!!!!------------------ numerical parameters
call MoveTO(51,'B.1 :')
read(51,*) zero

!!! Integration tolerance
call MoveTO(51,'B.2 :')
read(51,*) toleranceINT

!!! Maximal size of step in the Runge-Kutta iteraction (for t=log[mu])
call MoveTO(51,'B.3 :')
read(51,*) RGh_max

!!! Maximum number of processors to use
call MoveTO(51,'C.1 :')
read(51,*) allowedNumProcessor

!!! Initialize Chiral-Even kernels
call MoveTO(51,'D.1 :')
read(51,*) IncludeChiralEvenEvolution
!!! Initialize Chiral-Odd kernels
call MoveTO(51,'D.2 :')
read(51,*) IncludeChiralOddEvolution
!!! Take into account gluon and flavor mixing (singlet evolution)
call MoveTO(51,'D.3 :')
read(51,*) useSingletEvolution
!!! Mass of CHARM threshold [GeV]
call MoveTO(51,'D.4 :')
read(51,*) massCHARM
!!! Mass of BOTTOM threshold [GeV]
call MoveTO(51,'D.5 :')
read(51,*) massBOTTOM
!!! Initialize G2 matrix
call MoveTO(51,'D.6 :')
read(51,*) IncludeG2Matrix
!!! Initialize WGT matrix
call MoveTO(51,'D.7 :')
read(51,*) IncludeWGTMatrix


!!! Read kernels from files
call MoveTO(51,'X.1 :')
read(51,*) readKERNEL
if(readKERNEL) then
    !!! Path to the kernels
    call MoveTO(51,'X.2 :')
    read(51,'(A)') pathToKernels
end if

CLOSE (51, STATUS='KEEP')

!$ call OMP_set_num_threads(allowedNumProcessor)

if(showINI) then
    write(*,*) "Estimated zero:",zero
    write(*,*) "Tolerance of integration:",toleranceINT
    write(*,*) "Maximal step of Runge-Kutta:",RGh_max
#if GK_ORDER==1
write(*,*) ">>> Using GK7/15 for integration"
#elif GK_ORDER==2
write(*,*) ">>> Using GK20/41 for integration"
#else
write(*,*) ">>> Using GK7/15 for integration"
#endif
!$ write(*,'("Parallel computation is ",A," with ",I3," processors.")') color("ON",c_yellow),allowedNumProcessor
end if

if(readKERNEL) then
    call ReadKernels(pathToKernels)
else

    call cpu_time(t1)
    !$ t1=omp_get_wtime()


    if(IncludeChiralEvenEvolution) then
            !---------------------------- Hqq
        call cpu_time(tt1)
        !$ tt1=omp_get_wtime()

        write(*,'(A)',advance="no") "Preparing Hqq kernels ..."
        allocate(Hqq(0:NUM_TOT))
        call PreComputeMatrix(kernelHqq, Hqq, sparsity)

        call cpu_time(tt2)
        !$ tt2=omp_get_wtime()

        write(*,*) " Complete. Computation time =",tt2-tt1

!         !!-----------------
!         return
!         !!------------------

        if(showINI) write(*,*) "Sparsity of Hqq =",real(sparsity,dp)/NUM_TOT/NUM_TOT

        if(useSingletEvolution) then
        !!!!! in thecase of Singlet evolution mixing with gluon is accounted

        !---------------------------- Hqq - J
        call cpu_time(tt1)
        !$ tt1=omp_get_wtime()

        write(*,'(A)',advance="no") "Preparing Hqq' kernels ..."
        allocate(Hqq_J(0:NUM_TOT))
        call PreComputeMatrix(kernelHqq_J, Hqq_J, sparsity)

        call cpu_time(tt2)
        !$ tt2=omp_get_wtime()

        write(*,*) " Complete. Computation time =",tt2-tt1
        if(showINI) write(*,*) "Sparsity of Hqq =",real(sparsity,dp)/NUM_TOT/NUM_TOT


        !---------------------------- Hqg plus
        call cpu_time(tt1)
        !$ tt1=omp_get_wtime()

        write(*,'(A)',advance="no") "Preparing Hqg+ kernels ..."
        allocate(Hqg_PLUS(0:NUM_TOT))
        call PreComputeMatrix(kernelHqg_plus, Hqg_PLUS, sparsity)

        call cpu_time(tt2)
        !$ tt2=omp_get_wtime()

        write(*,*) " Complete. Computation time =",tt2-tt1
        if(showINI) write(*,*) "Sparsity of Hqg+ =",real(sparsity,dp)/NUM_TOT/NUM_TOT

        !---------------------------- Hgq plus
        call cpu_time(tt1)
        !$ tt1=omp_get_wtime()

        write(*,'(A)',advance="no") "Preparing Hgq+ kernels ..."
        allocate(Hgq_PLUS(0:NUM_TOT))
        call PreComputeMatrix(kernelHgq_plus, Hgq_PLUS, sparsity)

        call cpu_time(tt2)
        !$ tt2=omp_get_wtime()

        write(*,*) " Complete. Computation time =",tt2-tt1
        if(showINI) write(*,*) "Sparsity of Hgq+ =",real(sparsity,dp)/NUM_TOT/NUM_TOT

        !---------------------------- Hgg plus
        call cpu_time(tt1)
        !$ tt1=omp_get_wtime()

        write(*,'(A)',advance="no") "Preparing Hgg+ kernels ..."
        allocate(Hgg_PLUS(0:NUM_TOT))
        call PreComputeMatrix(kernelHgg_plus, Hgg_PLUS, sparsity)

        call cpu_time(tt2)
        !$ tt2=omp_get_wtime()

        write(*,*) " Complete. Computation time =",tt2-tt1
        if(showINI) write(*,*) "Sparsity of Hgg+ =",real(sparsity,dp)/NUM_TOT/NUM_TOT

        !---------------------------- Hqg plus
        call cpu_time(tt1)
        !$ tt1=omp_get_wtime()

        write(*,'(A)',advance="no") "Preparing Hqg- kernels ..."
        allocate(Hqg_MINUS(0:NUM_TOT))
        call PreComputeMatrix(kernelHqg_minus, Hqg_MINUS, sparsity)

        call cpu_time(tt2)
        !$ tt2=omp_get_wtime()

        write(*,*) " Complete. Computation time =",tt2-tt1
        if(showINI) write(*,*) "Sparsity of Hqg- =",real(sparsity,dp)/NUM_TOT/NUM_TOT

        !---------------------------- Hgq plus
        call cpu_time(tt1)
        !$ tt1=omp_get_wtime()

        write(*,'(A)',advance="no") "Preparing Hgq- kernels ..."
        allocate(Hgq_MINUS(0:NUM_TOT))
        call PreComputeMatrix(kernelHgq_minus, Hgq_MINUS, sparsity)

        call cpu_time(tt2)
        !$ tt2=omp_get_wtime()

        write(*,*) " Complete. Computation time =",tt2-tt1
        if(showINI) write(*,*) "Sparsity of Hgq- =",real(sparsity,dp)/NUM_TOT/NUM_TOT

        !---------------------------- Hgg plus
        call cpu_time(tt1)
        !$ tt1=omp_get_wtime()

        write(*,'(A)',advance="no") "Preparing Hgg- kernels ..."
        allocate(Hgg_MINUS(0:NUM_TOT))
        call PreComputeMatrix(kernelHgg_minus, Hgg_MINUS, sparsity)

        call cpu_time(tt2)
        !$ tt2=omp_get_wtime()

        write(*,*) " Complete. Computation time =",tt2-tt1
        if(showINI) write(*,*) "Sparsity of Hgg- =",real(sparsity,dp)/NUM_TOT/NUM_TOT
        end if
    end if

    if(IncludeChiralOddEvolution) then
        !---------------------------- Hqq - J
        call cpu_time(tt1)
        !$ tt1=omp_get_wtime()

        write(*,'(A)',advance="no") "Preparing Hqq (chiral odd)  kernels ..."
        allocate(Hqq_CO(0:NUM_TOT))
        call PreComputeMatrix(kernelHqq_CO, Hqq_CO, sparsity)

        call cpu_time(tt2)
        !$ tt2=omp_get_wtime()

        write(*,*) " Complete. Computation time =",tt2-tt1
        if(showINI) write(*,*) "Sparsity of Hqq_CO =",real(sparsity,dp)/NUM_TOT/NUM_TOT
    end if

    if(IncludeG2Matrix) then
        call cpu_time(tt1)
        !$ tt1=omp_get_wtime()

        write(*,'(A)',advance="no") "Preparing G2 projection matrix ..."
        allocate(G2matrix(0:NUM_RHO,0:NUM_TOT))
        call PreComputeMatrixG2(G2matrix)

        call cpu_time(tt2)
        !$ tt2=omp_get_wtime()

        write(*,*) " Complete. Computation time =",tt2-tt1

        call cpu_time(tt1)
        !$ tt1=omp_get_wtime()

        write(*,'(A)',advance="no") "Preparing D2 projection matrix ..."
        allocate(D2matrix(0:NUM_TOT))
        call PreComputeMatrixD2(D2matrix,G2matrix)

        call cpu_time(tt2)
        !$ tt2=omp_get_wtime()

        write(*,*) " Complete. Computation time =",tt2-tt1
    end if

    if(IncludeWGTMatrix) then
        call cpu_time(tt1)
        !$ tt1=omp_get_wtime()

        write(*,'(A)',advance="no") "Preparing WGT projection matrices ..."
        allocate(WGTmatrixP(0:NUM_RHO,0:NUM_TOT))
        allocate(WGTmatrixN(0:NUM_RHO,0:NUM_TOT))
        call PreComputeMatrixWGT(WGTmatrixP,WGTmatrixN)

        call cpu_time(tt2)
        !$ tt2=omp_get_wtime()

        write(*,*) " Complete. Computation time =",tt2-tt1

        call cpu_time(tt1)
        !$ tt1=omp_get_wtime()
    end if

    call cpu_time(t2)
    !$ t2=omp_get_wtime()

    if(showINI) write(*,*) "--------- All necesary kernels computed. Total time =",t2-t1
end if

end subroutine EvolutionKernels_Initialize

!!!! This subroutine computes the kernel-matrix by the function Mfunction
!!!! and stores it into the COMPACT-matrix M
!!!! the function must have the interface (n,x1,x2) where n is any integer
!!!! the result is the matrix [n,m(x1,x2)]
subroutine PreComputeMatrix(Mfunction,M,sparsity)
type(SparseM2), dimension(0:)::M
real(dp)::Mfunction
integer::n,nn,counter,sparsity,NUM_M
real(dp)::x1,x2,x3
integer::OMP_GET_THREAD_NUM,thread_id

real(dp),dimension(0:NUM_TOT):: Hinter

sparsity=0

NUM_M=size(M)-1

!allocate(M(0:NUM_TOT))

!$OMP PARALLEL DO private(counter,thread_id,x1,x2,x3,nn,Hinter)
do n=0,NUM_M
    thread_id = OMP_GET_THREAD_NUM()

    !!!! first: compute the matrix for element (n,k)
    counter=0
    do nn=0,NUM_TOT

        call get_X123_from_1Dindex(n,x1,x2,x3)
        Hinter(nn)=Mfunction(nn,x1,x2)

        !write(*,*) n,nn,Hinter(nn)

        if(abs(Hinter(nn))>zero) counter=counter+1
    end do

    !!!! second: allocate the value of M
    allocate(M(n)%i(1:counter))
    allocate(M(n)%val(1:counter))

    !!!! third: fill the elements of M
    counter=1
    do nn=0,NUM_TOT
        if(abs(Hinter(nn))>zero) then
            M(n)%i(counter)=nn
            M(n)%val(counter)=Hinter(nn)
            counter=counter+1
        end if
    end do

end do
!$OMP END PARALLEL DO

sparsity=0
do n=0,NUM_M
    sparsity=sparsity+size(M(n)%i)
end do

end subroutine PreComputeMatrix

!!!!!----------------------------------------------------------------------------
!!!!!-------------------- Routines associated with saving and reading -----------
!!!!!------------------------- previously stored kernels ------------------------
!!!!!----------------------------------------------------------------------------

!!!!!! Save Kernel into a directory
subroutine SaveKernels(path)
character(len=*)::path

write(*,'(A)',advance="no") "Saving Hqq ..."
call SaveKernelMatrix(Hqq,trim(path)//trim("Hqq.ker"))
write(*,'(A)',advance="no") "Saving Hqq_J ..."
call SaveKernelMatrix(Hqq_J,trim(path)//trim("Hqq_J.ker"))
write(*,'(A)',advance="no") "Saving Hqq_CO ..."
call SaveKernelMatrix(Hqq_CO,trim(path)//trim("Hqq_CO.ker"))
write(*,'(A)',advance="no") "Saving Hqg_PLUS ..."
call SaveKernelMatrix(Hqg_PLUS,trim(path)//trim("Hqg_PLUS.ker"))
write(*,'(A)',advance="no") "Saving Hgq_PLUS ..."
call SaveKernelMatrix(Hgq_PLUS,trim(path)//trim("Hgq_PLUS.ker"))
write(*,'(A)',advance="no") "Saving Hgg_PLUS ..."
call SaveKernelMatrix(Hgg_PLUS,trim(path)//trim("Hgg_PLUS.ker"))
write(*,'(A)',advance="no") "Saving Hqg_MINUS ..."
call SaveKernelMatrix(Hqg_MINUS,trim(path)//trim("Hqg_MINUS.ker"))
write(*,'(A)',advance="no") "Saving Hgq_MINUS ..."
call SaveKernelMatrix(Hgq_MINUS,trim(path)//trim("Hgq_MINUS.ker"))
write(*,'(A)',advance="no") "Saving Hgg_MINUS ..."
call SaveKernelMatrix(Hgg_MINUS,trim(path)//trim("Hgg_MINUS.ker"))
if(IncludeG2Matrix) then
    write(*,'(A)',advance="no") "Saving G2 projector ..."
    call SaveMatrix(G2matrix,trim(path)//trim("G2.ker"))
    !!!call PreComputeMatrixD2(D2matrix,G2matrix)
    write(*,'(A)',advance="no") "Saving D2 projector ..."
    call Save1DMatrix(D2matrix,trim(path)//trim("D2.ker"))
end if
if(IncludeWGTMatrix) then
    !call PreComputeMatrixWGT(WGTmatrixP,WGTmatrixN)
    write(*,'(A)',advance="no") "Saving WGT(positive) projector ..."
    call SaveMatrix(WGTmatrixP,trim(path)//trim("WGTp.ker"))
    write(*,'(A)',advance="no") "Saving WGT(negative) projector ..."
    call SaveMatrix(WGTmatrixN,trim(path)//trim("WGTn.ker"))
end if
write(*,'(A)') "All kernels are stored in "//trim(path)

end subroutine SaveKernels

!!!!!! Save matrix For the Kernel into a file
!!!!!! i..e it stores the sparse object to the file
subroutine SaveKernelMatrix(M,path)
type(SparseM2),dimension(0:),intent(in)::M
character(len=*)::path
integer::i,j,num_M

write(*,'(A)') trim(path)

num_M=size(M)-1 !!!! -1, because hte counting is from 0

OPEN(UNIT=51, FILE=path, ACTION="write", STATUS="replace")
write(51,*) num_M
write(51,*) NUM_RHO
write(51,*) NUM_PHI

do i=0,num_M
write(51,*) i, size(M(i)%i)
end do
do i=0,num_M
do j=1,size(M(i)%i)
write(51,*) M(i)%i(j), M(i)%val(j)
end do
end do

CLOSE(51, STATUS='KEEP')

end subroutine SaveKernelMatrix

!!!!!! Save matrix For the into a file
!!!!!! it stores the fixed-size matrix to the file
subroutine SaveMatrix(M,path)
real(dp),dimension(0:,0:),intent(in)::M
character(len=*)::path
integer::i,j,num1,num2

write(*,'(A)') trim(path)

num1=size(M(:,0))-1 !!!! -1, because the counting is from 0
num2=size(M(0,:))-1 !!!! -1, because the counting is from 0

OPEN(UNIT=51, FILE=path, ACTION="write", STATUS="replace")
write(51,*) num1
write(51,*) num2

do i=0,num1
do j=0,num2
write(51,*) M(i,j)
end do
end do

CLOSE(51, STATUS='KEEP')

end subroutine SaveMatrix

!!!!!! Save 1D matrix For the into a file
!!!!!! it stores the fixed-size matrix to the file
subroutine Save1DMatrix(M,path)
real(dp),dimension(0:),intent(in)::M
character(len=*)::path
integer::i,num1

write(*,'(A)') trim(path)

num1=size(M(:))-1 !!!! -1, because the counting is from 0

OPEN(UNIT=51, FILE=path, ACTION="write", STATUS="replace")
write(51,*) num1

do i=0,num1
write(51,*) M(i)
end do

CLOSE(51, STATUS='KEEP')

end subroutine Save1DMatrix

!!!!!! Read all Kernels from the directory
subroutine ReadKernels(path)
character(len=*)::path

if(showINI) write(*,*) "Reading Hqq ..."
if(allocated(Hqq)) deallocate(Hqq)
allocate(Hqq(0:NUM_TOT))
call ReadKernelMatrix(Hqq,trim(path)//trim("Hqq.ker"))

if(showINI) write(*,*) "Reading Hqq_J ..."
if(allocated(Hqq_J)) deallocate(Hqq_J)
allocate(Hqq_J(0:NUM_TOT))
call ReadKernelMatrix(Hqq_J,trim(path)//trim("Hqq_J.ker"))

if(showINI) write(*,*) "Reading Hqq_CO ..."
if(allocated(Hqq_CO)) deallocate(Hqq_CO)
allocate(Hqq_CO(0:NUM_TOT))
call ReadKernelMatrix(Hqq_CO,trim(path)//trim("Hqq_CO.ker"))

if(showINI) write(*,*) "Reading Hqg_PLUS ..."
if(allocated(Hqg_PLUS)) deallocate(Hqg_PLUS)
allocate(Hqg_PLUS(0:NUM_TOT))
call ReadKernelMatrix(Hqg_PLUS,trim(path)//trim("Hqg_PLUS.ker"))

if(showINI) write(*,*) "Reading Hgq_PLUS ..."
if(allocated(Hgq_PLUS)) deallocate(Hgq_PLUS)
allocate(Hgq_PLUS(0:NUM_TOT))
call ReadKernelMatrix(Hgq_PLUS,trim(path)//trim("Hgq_PLUS.ker"))

if(showINI) write(*,*) "Reading Hgg_PLUS ..."
if(allocated(Hgg_PLUS)) deallocate(Hgg_PLUS)
allocate(Hgg_PLUS(0:NUM_TOT))
call ReadKernelMatrix(Hgg_PLUS,trim(path)//trim("Hgg_PLUS.ker"))

if(showINI) write(*,*) "Reading Hqg_MINUS ..."
if(allocated(Hqg_MINUS)) deallocate(Hqg_MINUS)
allocate(Hqg_MINUS(0:NUM_TOT))
call ReadKernelMatrix(Hqg_MINUS,trim(path)//trim("Hqg_MINUS.ker"))

if(showINI) write(*,*) "Reading Hgq_MINUS ..."
if(allocated(Hgq_MINUS)) deallocate(Hgq_MINUS)
allocate(Hgq_MINUS(0:NUM_TOT))
call ReadKernelMatrix(Hgq_MINUS,trim(path)//trim("Hgq_MINUS.ker"))

if(showINI) write(*,*) "Reading Hgg_MINUS ..."
if(allocated(Hgg_MINUS)) deallocate(Hgg_MINUS)
allocate(Hgg_MINUS(0:NUM_TOT))
call ReadKernelMatrix(Hgg_MINUS,trim(path)//trim("Hgg_MINUS.ker"))

if(IncludeG2Matrix) then
    if(showINI) write(*,*) "Reading G2 ..."
    if(allocated(G2matrix)) deallocate(G2matrix)
    allocate(G2matrix(0:NUM_RHO,0:NUM_TOT))
    call ReadMatrix(G2matrix,trim(path)//trim("G2.ker"))

    if(showINI) write(*,*) "Reading D2 ..."
    if(allocated(D2matrix)) deallocate(D2matrix)
    allocate(D2matrix(0:NUM_TOT))
    call Read1DMatrix(D2matrix,trim(path)//trim("D2.ker"))

    !call PreComputeMatrixD2(D2matrix,G2matrix)
end if

if(IncludeWGTMatrix) then
    if(showINI) write(*,*) "Reading WGT(positive) ..."
    if(allocated(WGTmatrixP)) deallocate(WGTmatrixP)
    allocate(WGTmatrixP(0:NUM_RHO,0:NUM_TOT))
    call ReadMatrix(WGTmatrixP,trim(path)//trim("WGTp.ker"))

    if(showINI) write(*,*) "Reading WGT(negative) ..."
    if(allocated(WGTmatrixN)) deallocate(WGTmatrixN)
    allocate(WGTmatrixN(0:NUM_RHO,0:NUM_TOT))
    call ReadMatrix(WGTmatrixN,trim(path)//trim("WGTn.ker"))

    !call PreComputeMatrixWGT(WGTmatrixP,WGTmatrixN)
end if

end subroutine ReadKernels

!!!!!! Read matrix For the Kernel from a file
!!!!!! i.e. it reads the preveously stored sparse matrix
subroutine ReadKernelMatrix(M,path)
type(SparseM2),dimension(0:)::M
character(len=*)::path
integer::i,j,sizeM
integer::int1,int2,int0
real(dp)::real1
logical::fileExist

inquire( file=trim(path), exist=fileExist )
if(.not.fileExist) then
    write(*,*) ErrorString("file for kernel is not found. Check the path."," ")
    write(*,*) trim(path)
    stop
end if

sizeM=size(M)-1

OPEN(UNIT=51, FILE=trim(path), ACTION="read", STATUS="old")
!!! read size
read(51,*) int0
read(51,*) int1
read(51,*) int2

if(sizeM/=int0) then
    write(*,*) ErrorString("Saved kernels use different total-size than declared in INI-file"," ")
    write(*,'("Grid-file = ",I4 ," vs. INI-file = ",I4)') sizeM,int0
    write(*,*) "Setup INI-file appropriately. Evaluation STOP"
    CLOSE(51, STATUS='KEEP')
    stop
end if

if(NUM_RHO/=int1 .or. NUM_PHI/=int2) then
    write(*,*) ErrorString("Saved kernels use different grid-sizes than declared in INI-file"," ")
    write(*,'("Files = (",I4 ," x ",I4,")  vs. INI-file (",I4 ," x ",I4,")")') int1,int2,NUM_RHO,NUM_PHI
    write(*,*) "Setup INI-file appropriately. Evaluation STOP"
    CLOSE(51, STATUS='KEEP')
    stop
end if

!!! read sizes of individual SparseM2, and allocate them
do i=0,sizeM
    if(allocated(M(i)%i)) deallocate(M(i)%i)
    if(allocated(M(i)%val)) deallocate(M(i)%val)
    read(51,*) int1, int2
    if(int1/=i) then
        write(*,*) ErrorString("Kernel file corrupted. ERROR1"," ")
        write(*,*) trim(path)
        write(*,*) "Evaluation STOP"
        CLOSE(51, STATUS='KEEP')
        stop
    end if

    allocate(M(i)%i(1:int2),M(i)%val(1:int2))
end do

!!! Finally read enties
do i=0,sizeM
   do j=1,size(M(i)%i)
    read(51,*) int1,real1
    M(i)%i(j)=int1
    M(i)%val(j)=real1

    if(abs(real1)>100000) write(*,*) i,int1,real1
    if(ISNAN(real1)) write(*,*) "file currepted",i,int1,real1
   end do
end do

CLOSE(51, STATUS='KEEP')

end subroutine ReadKernelMatrix

!!!!!! Read matrix For the Kernel from a file
!!!!!! i.e. it reads the preveously stored normal matrix
!!!!!! in this matrix the numbering must strat from zero (0:,0:)
subroutine ReadMatrix(M,path)
real(dp),dimension(0:,0:)::M
character(len=*)::path
integer::i,j,num1,num2
integer::int1,int2
real(dp)::real1
logical::fileExist

inquire( file=trim(path), exist=fileExist )
if(.not.fileExist) then
    write(*,*) ErrorString("file for kernel is not found. Check the path."," ")
    write(*,*) trim(path)
    stop
end if

num1=size(M(:,0))-1 !!!! -1, because the counting is from 0
num2=size(M(0,:))-1 !!!! -1, because the counting is from 0

OPEN(UNIT=51, FILE=trim(path), ACTION="read", STATUS="old")
!!! read size
read(51,*) int1
read(51,*) int2


if(num1/=int1 .or. num2/=int2) then
    write(*,*) ErrorString("Saved matrix uses different of different shape than required"," ")
    write(*,'("Files = (",I4 ," x ",I4,")  vs. requested (",I4 ," x ",I4,")")') int1,int2,num1,num2
    write(*,*) "Setup appropriately. Evaluation STOP"
    CLOSE(51, STATUS='KEEP')
    stop
end if

!!! Finally read enties
do i=0,num1
   do j=0,num2
    read(51,*) real1
    M(i,j)=real1

    if(abs(real1)>100000) write(*,*) i,int1,real1
    if(ISNAN(real1)) write(*,*) "file currepted",i,int1,real1
   end do
end do

CLOSE(51, STATUS='KEEP')

end subroutine ReadMatrix

!!!!!! Read 1D matrix For the Kernel from a file
!!!!!! i.e. it reads the preveously stored normal matrix
!!!!!! in this matrix the numbering must strat from zero (0:)
subroutine Read1DMatrix(M,path)
real(dp),dimension(0:)::M
character(len=*)::path
integer::i,num1
integer::int1
real(dp)::real1
logical::fileExist

inquire( file=trim(path), exist=fileExist )
if(.not.fileExist) then
    write(*,*) ErrorString("file for kernel is not found. Check the path."," ")
    write(*,*) trim(path)
    stop
end if

num1=size(M(:))-1 !!!! -1, because the counting is from 0

OPEN(UNIT=51, FILE=trim(path), ACTION="read", STATUS="old")
!!! read size
read(51,*) int1


if(num1/=int1) then
    write(*,*) ErrorString("Saved 1D matrix uses different of different shape than required"," ")
    write(*,'("Files = ",I4 ,"  vs. requested ",I4 )') int1,num1
    write(*,*) "Setup appropriately. Evaluation STOP"
    CLOSE(51, STATUS='KEEP')
    stop
end if

!!! Finally read enties
do i=0,num1
    read(51,*) real1
    M(i)=real1
    if(abs(real1)>100000) write(*,*) i,int1,real1
    if(ISNAN(real1)) write(*,*) "file currepted",i,int1,real1
end do

CLOSE(51, STATUS='KEEP')

end subroutine Read1DMatrix

!!!!!----------------------------------------------------------------------------
!!!!!-------------- Routines associated with evolution kernels ------------------
!!!!!----------------------------------------------------------------------------

function TEST(F)
real(dp),intent(in),dimension(0:NUM_TOT)::F
real(dp),dimension(0:NUM_TOT)::TEST
TEST=HxF(Hqq,F)
end function

!!!!!!! composition of elementary kernels
!!!----- quark-quark part
function kernelHqq(n,x1,x2)
real(dp)::kernelHqq
integer,intent(in)::n
real(dp),intent(in)::x1,x2
       !kernelHqq=H13plus(n,x1,x2)
! !
     kernelHqq=CA*(H12hat(n,x1,x2)+H23hat(n,x1,x2)-2._dp*H12plus(n,x1,x2))&
        -1._dp/CA*(H13hat(n,x1,x2)-H13plus(n,x1,x2)-H23eP(n,x1,x2)+2*H12minus(n,x1,x2))&
       -3._dp*CF*GETinterpolatorB(n,x1,x2)
end function kernelHqq

!!!----- quark-quark mixture part
function kernelHqq_J(n,x1,x2)
real(dp)::kernelHqq_J
integer,intent(in)::n
real(dp),intent(in)::x1,x2
    kernelHqq_J=4*H13d(n,x1,x2)
end function kernelHqq_J

!!!----- quark-quark part (chiral odd)
function kernelHqq_CO(n,x1,x2)
real(dp)::kernelHqq_CO
integer,intent(in)::n
real(dp),intent(in)::x1,x2
    kernelHqq_CO=CA*(H12hat(n,x1,x2)+H23hat(n,x1,x2)-2*H12plus(n,x1,x2)-2*H23plus(n,x1,x2))&
       -1._dp/CA*(H13hat(n,x1,x2)+2*H12minus(n,x1,x2)+2*H23minus(n,x1,x2))&
       -3._dp*CF*GETinterpolatorB(n,x1,x2)
end function kernelHqq_CO

!!!----- gluon-gluon mixture part [WITHOUT BETA0 PART]
function kernelHgg_plus(n,x1,x2)
real(dp)::kernelHgg_plus
integer,intent(in)::n
real(dp),intent(in)::x1,x2
    kernelHgg_plus=CA*(H12hat_FFF(n,x1,x2)+H23hat_FFF(n,x1,x2)+H13hat_FFF(n,x1,x2) &
       -4*(H12plus_FFF(n,x1,x2)+H13plus_FFF(n,x1,x2)) &
       -2*(H12tilde_FFF(n,x1,x2)+H13tilde_FFF(n,x1,x2)) &
       +6*(H12minus_FFF(n,x1,x2)+H13minus_FFF(n,x1,x2)))
end function kernelHgg_plus

!!!----- gluon-gluon mixture part [WITHOUT BETA0 PART]
function kernelHgg_minus(n,x1,x2)
real(dp)::kernelHgg_minus
integer,intent(in)::n
real(dp),intent(in)::x1,x2
    kernelHgg_minus=CA*(H12hat_FFF(n,x1,x2)+H23hat_FFF(n,x1,x2)+H13hat_FFF(n,x1,x2) &
        -4*(H12plus_FFF(n,x1,x2)+H13plus_FFF(n,x1,x2)) &
        -2*(H12tilde_FFF(n,x1,x2)+H13tilde_FFF(n,x1,x2)) &
        -6*(H12minus_FFF(n,x1,x2)+H13minus_FFF(n,x1,x2)))
end function kernelHgg_minus

!!!----- quark-gluon mixture part
function kernelHqg_plus(n,x1,x2)
real(dp)::kernelHqg_plus
integer,intent(in)::n
real(dp),intent(in)::x1,x2
    kernelHqg_plus=V13plus(n,x1,x2)-V13minus(n,x1,x2)
end function kernelHqg_plus

!!!----- quark-gluon mixture part
function kernelHqg_minus(n,x1,x2)
real(dp)::kernelHqg_minus
integer,intent(in)::n
real(dp),intent(in)::x1,x2
    kernelHqg_minus=V13plus(n,x1,x2)+V13minus(n,x1,x2)
end function kernelHqg_minus

!!!----- gluon-quark mixture part
function kernelHgq_plus(n,x1,x2)
real(dp)::kernelHgq_plus
integer,intent(in)::n
real(dp),intent(in)::x1,x2
    kernelHgq_plus=CA*(Wplus(n,x1,x2)+Wminus(n,x1,x2)-2*DeltaW(n,x1,x2)&
                       - WplusP(n,x1,x2)-WminusP(n,x1,x2)+2*DeltaWP(n,x1,x2))
end function kernelHgq_plus

!!!----- gluon-quark mixture part
function kernelHgq_minus(n,x1,x2)
real(dp)::kernelHgq_minus
integer,intent(in)::n
real(dp),intent(in)::x1,x2
    kernelHgq_minus=ColorQG*(Wplus(n,x1,x2)+Wminus(n,x1,x2)+ WplusP(n,x1,x2)+WminusP(n,x1,x2))
end function kernelHgq_minus

!!!! multiplies the sparse ``H''-matrix by ``vector'' F (0... NUM_TOT)
!!!! the dimension of the multiplication is defined by the H which is (w x NUM_TOT), where w=0...?
!!!! so for evolution kernels w=0...NUM_TOT
!!!! so for G2 matrix w=0...NUM_RHO
function HxF(H,F)
    real(dp),dimension(0:NUM_TOT),intent(in)::F
    type(SparseM2),dimension(0:),intent(in)::H
    real(dp),dimension(0:size(H)-1)::HxF

    integer::c,n

    !$OMP PARALLEL DO private(n)
    do c=0,size(H)-1
        HxF(c)=0._dp
        do n=1,size(H(c)%i)
            HxF(c)=HxF(c)+H(c)%val(n)*F(H(c)%i(n))
        end do
    end do
    !$OMP END PARALLEL DO
end function

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Runge Kuta-4 within vector space !!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! various combinations !!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! All functions have the same interface
!!!! F is the initial boundary condition as the vector in S-space
!!!! alpha is the a(t) is alpha_s/4pi(log mu^2)
!!!! t0, t1 are initial and final scales of evolution
!!!! nf is the active number of flavors
!!!! output is placed into F0

subroutine EvNonSinglet(F,alpha,t0,t1)
    real(dp),dimension(0:NUM_TOT)::F
    real(dp)::alpha
    real(dp),intent(in)::t0,t1

    real(dp)::tn,RGh
    real(dp),dimension(0:NUM_TOT)::k1,k2,k3,k4
    integer::i,numSteps

   !!! number of steps, and the step size
    numSteps=int((t1-t0)/RGh_max)!+1
    RGh=(t1-t0)/numSteps

    tn=t0
    do i=1,numSteps
        !!! the factor -1 is the factor in the evolution equation

        k1=-RGh*alpha(tn)*HxF(Hqq,F)
        k2=-RGh*alpha(tn+RGh/2)*HxF(Hqq,F+k1/2)
        k3=-RGh*alpha(tn+RGh/2)*HxF(Hqq,F+k2/2)
        k4=-RGh*alpha(tn+RGh)*HxF(Hqq,F+k3)

        F=F+(k1+2*k2+2*k3+k4)/6
        tn=tn+RGh
    end do
end subroutine EvNonSinglet

subroutine EvSingletPLUS(F,Fg,alpha,t0,t1,nf)
    real(dp),dimension(0:NUM_TOT)::F,Fg
    real(dp)::alpha
    real(dp),intent(in)::t0,t1
    integer,intent(in)::nf

    real(dp)::tn,RGh
    real(dp),dimension(0:NUM_TOT)::k1,k2,k3,k4,k1g,k2g,k3g,k4g,F_dum,Fg_dum
    integer::i,numSteps


    !!! number of steps, and the step size
    numSteps=int((t1-t0)/RGh_max)+1
    RGh=(t1-t0)/numSteps
    tn=t0
    do i=1,numSteps
        !!! the factor -1 is the factor in the evolution equation
        k1=-RGh*alpha(tn)*(HxF(Hqq,F)+Nf*(HxF(Hqq_J,F)+HxF(Hqg_PLUS,Fg)))
        k1g=-RGh*alpha(tn)*(HxF(Hgq_PLUS,F)+HxF(Hgg_PLUS,Fg)-(11._dp-2._dp*nf/3._dp)*Fg)

        F_dum=F+k1/2
        Fg_dum=Fg+k1g/2

        k2=-RGh*alpha(tn+RGh/2)*(HxF(Hqq,F_dum)+Nf*(HxF(Hqq_J,F_dum)+HxF(Hqg_PLUS,Fg_dum)))
        k2g=-RGh*alpha(tn+RGh/2)*(HxF(Hgq_PLUS,F_dum)+HxF(Hgg_PLUS,Fg_dum)-(11._dp-2._dp*nf/3._dp)*Fg_dum)

        F_dum=F+k2/2
        Fg_dum=Fg+k2g/2

        k3=-RGh*alpha(tn+RGh/2)*(HxF(Hqq,F_dum)+Nf*(HxF(Hqq_J,F_dum)+HxF(Hqg_PLUS,Fg_dum)))
        k3g=-RGh*alpha(tn+RGh/2)*(HxF(Hgq_PLUS,F_dum)+HxF(Hgg_PLUS,Fg_dum)-(11._dp-2._dp*nf/3._dp)*Fg_dum)

        F_dum=F+k3
        Fg_dum=Fg+k3g

        k4=-RGh*alpha(tn+RGh)*(HxF(Hqq,F_dum)+Nf*(HxF(Hqq_J,F_dum)+HxF(Hqg_PLUS,Fg_dum)))
        k4g=-RGh*alpha(tn+RGh)*(HxF(Hgq_PLUS,F_dum)+HxF(Hgg_PLUS,Fg_dum)-(11._dp-2._dp*nf/3._dp)*Fg_dum)

        F=F+(k1+2*k2+2*k3+k4)/6
        Fg=Fg+(k1g+2*k2g+2*k3g+k4g)/6
        tn=tn+RGh



     end do

end subroutine EvSingletPLUS

subroutine EvSingletMINUS(F,Fg,alpha,t0,t1,nf)
    real(dp),dimension(0:NUM_TOT)::F,Fg
    real(dp)::alpha
    real(dp),intent(in)::t0,t1

    real(dp)::tn,RGh
    real(dp),dimension(0:NUM_TOT)::k1,k2,k3,k4,k1g,k2g,k3g,k4g,F_dum,Fg_dum
    integer::nf
    integer::i,numSteps

    !!! number of steps, and the step size
    numSteps=int((t1-t0)/RGh_max)+1
    RGh=(t1-t0)/numSteps
    tn=t0

   do i=1,numSteps
        !!! the factor -1 is the factor in the evolution equation

        k1=-RGh*alpha(tn)*(HxF(Hqq,F)+Nf*HxF(Hqg_MINUS,Fg))
        k1g=-RGh*alpha(tn)*(HxF(Hgq_MINUS,F)+HxF(Hgg_MINUS,Fg)-(11._dp-2._dp*nf/3._dp)*Fg)

        F_dum=F+k1/2
        Fg_dum=Fg+k1g/2

        k2=-RGh*alpha(tn+RGh/2)*(HxF(Hqq,F_dum)+Nf*HxF(Hqg_MINUS,Fg_dum))
        k2g=-RGh*alpha(tn+RGh/2)*(HxF(Hgq_MINUS,F_dum)+HxF(Hgg_MINUS,Fg_dum)-(11._dp-2._dp*nf/3._dp)*Fg_dum)

        F_dum=F+k2/2
        Fg_dum=Fg+k2g/2

        k3=-RGh*alpha(tn+RGh/2)*(HxF(Hqq,F_dum)+Nf*HxF(Hqg_MINUS,Fg_dum))
        k3g=-RGh*alpha(tn+RGh/2)*(HxF(Hgq_MINUS,F_dum)+HxF(Hgg_MINUS,Fg_dum)-(11._dp-2._dp*nf/3._dp)*Fg_dum)

        F_dum=F+k3
        Fg_dum=Fg+k3g

        k4=-RGh*alpha(tn+RGh)*(HxF(Hqq,F_dum)+Nf*HxF(Hqg_MINUS,Fg_dum))
        k4g=-RGh*alpha(tn+RGh)*(HxF(Hgq_MINUS,F_dum)+HxF(Hgg_MINUS,Fg_dum)-(11._dp-2._dp*nf/3._dp)*Fg_dum)

        F=F+(k1+2*k2+2*k3+k4)/6
        Fg=Fg+(k1g+2*k2g+2*k3g+k4g)/6

        tn=tn+RGh

    end do

end subroutine EvSingletMINUS

subroutine EvChiralOdd(F,alpha,t0,t1)
    real(dp),dimension(0:NUM_TOT)::F
    real(dp)::alpha
    real(dp),intent(in)::t0,t1

    real(dp)::tn,RGh
    real(dp),dimension(0:NUM_TOT)::k1,k2,k3,k4
    integer::i,numSteps

    !!! number of steps, and the step size
    numSteps=int((t1-t0)/RGh_max)+1
    RGh=(t1-t0)/numSteps
    tn=t0

    !!! all steps except the last
    do i=1,numSteps
        k1=-RGh*alpha(tn)*HxF(Hqq_CO,F)
        k2=-RGh*alpha(tn+RGh/2)*HxF(Hqq_CO,F+k1/2)
        k3=-RGh*alpha(tn+RGh/2)*HxF(Hqq_CO,F+k2/2)
        k4=-RGh*alpha(tn+RGh)*HxF(Hqq_CO,F+k3)

        F=F+(k1+2*k2+2*k3+k4)/6
        tn=tn+RGh
    end do

end subroutine EvChiralOdd

!!!!-------------------------------------------------------------------------------------------------------------
!!! Gauss-Kronrod 7/15 adaptive
!!! f::  function of 1 variable
!!! xMin, and xMax boundaries of the integral. xMax>xMin !!
!!! tolerance is relative (it is weighted by approximate value of integral)
function Integrate_GK(f,xMin,xMax)
    real(dp)::f,Integrate_GK
    real(dp),intent(in)::xMin,xMax
    real(dp)::delta,av,g7,k15,eps,fI
    integer::i

    delta=(xMax-xMin)/2._dp
    av=(xMax+xMin)/2._dp

    if(delta<zero) then
        Integrate_GK=0._dp
        return
    end if

    g7=0._dp
    k15=0._dp
#if GK_ORDER==1
    do i=1,15
        fI=f(Xi_k15(i)*delta+av)
        g7=g7+Wi_g7(i)*fI
        k15=k15+Wi_k15(i)*fI
    end do
#elif  GK_ORDER==2
    do i=1,41
        fI=f(Xi_k41(i)*delta+av)
        g7=g7+Wi_g20(i)*fI
        k15=k15+Wi_k41(i)*fI
    end do
#endif

    if(delta*abs(k15)>1.) then
        eps=delta*abs(k15)*toleranceINT
    else
        eps=toleranceINT
    end if

!    write(*,*) "-->",delta*abs(k15-g7),eps,delta*k15
    if(abs(delta*k15)<zero) then
        Integrate_GK=delta*k15
    else if(delta*abs(k15-g7)>eps) then
        Integrate_GK=GK_Rec(f,xMin,av,eps)+GK_Rec(f,av,xMax,eps)
    else
        Integrate_GK=delta*k15
    end if

end function Integrate_GK

recursive function GK_Rec(f,xMin,xMax,eps) result(res)
    real(dp)::f,res
    real(dp),intent(in)::xMin,xMax
    real(dp)::delta,av,g7,k15,eps,fI
    integer::i

    delta=(xMax-xMin)/2._dp
    av=(xMax+xMin)/2._dp

    if(delta<zero) then
        res=0._dp
        return
    end if

    g7=0._dp
    k15=0._dp
#if GK_ORDER==1
    do i=1,15
        fI=f(Xi_k15(i)*delta+av)
        g7=g7+Wi_g7(i)*fI
        k15=k15+Wi_k15(i)*fI
    end do
#elif  GK_ORDER==2
    do i=1,41
        fI=f(Xi_k41(i)*delta+av)
        g7=g7+Wi_g20(i)*fI
        k15=k15+Wi_k41(i)*fI
    end do
#endif

    if(abs(delta*k15)<zero) then
        res=delta*k15
    else if(delta*abs(k15-g7)>eps) then
        res=GK_Rec(f,xMin,av,eps)+GK_Rec(f,av,xMax,eps)
    else
        res=delta*k15
    end if

end function GK_Rec
!!!!-------------------------------------------------------------------------------------------------------------

end module EvolutionKernels
