!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                          !!
!! This module is the central module for SnowFlake library.                                 !!
!! It decompose the input and call correponsing evolution routines                          !!
!!                                                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module SnowFlake
use IO_snowflake
use HexGrid
use EvolutionKernels
use SnowFlake_Model
implicit none

private

INCLUDE 'commonVariables.f90'

!!!!! These global variables save the result of current evolution
!!!!! first index is the Q-index
!!!!! These are C-odd functions
real(dp),allocatable:: Uplus(:,:),Dplus(:,:),Splus(:,:),Cplus(:,:),Bplus(:,:),Gplus(:,:)
real(dp),allocatable:: Uminus(:,:),Dminus(:,:),Sminus(:,:),Cminus(:,:),Bminus(:,:),Gminus(:,:)
!!!!! These are C-even functions
real(dp),allocatable:: Uodd(:,:),Dodd(:,:),Sodd(:,:),Codd(:,:),Bodd(:,:)

!!!!! Global variables for saving the Q-grid
!!!!! the grid is saved with Q=Qmin*exp(t/2), t=2 log(Q/Qmin)
!!!!! i.e. t=0,...,maxT
!!!!! step for saving the grids
real(dp)::stepT
!!!!! ultimate points of grid
real(dp)::Qmin,Qmax
!!!!! ultimate points of grid
integer::maxT

!!!!! flags to setup if the evolution is already preapread or not
logical:: evolutionEvenPrepared,evolutionOddPrepared

public:: SnowFlake_Initialize,ComputeEvolution, ComputeEvolutionChiralOdd, GetPDF,GetPDFChiralOdd
public:: G2,D2,G2_list,D2_List,WGT,WGT_fList

interface QFromt
    module procedure QFromt_int, QFromt_real
end interface

contains


subroutine SnowFlake_Initialize(file,prefix)
character(len=*)::file
character(len=*),optional::prefix
character(len=300)::path
!$ real*8::omp_get_wtime
real*8::t1,t2

if(present(prefix)) then
    path=trim(adjustl(prefix))//trim(adjustr(file))
else
    path=trim(adjustr(file))
end if

call cpu_time(t1)
!$ t1=omp_get_wtime()

write(*,*) color("----------------------- SNOWFLAKE V2.0 -----------------------",c_cyan)
write(*,*) color("                           /\‾‾‾/\                            ",c_cyan)
write(*,*) color("                          /  \ /  \                           ",c_cyan)
write(*,*) color("                          ----o----                           ",c_cyan)
write(*,*) color("                          \  / \  /                           ",c_cyan)
write(*,*) color("                           \/___\/                            ",c_cyan)
write(*,*) color("- please, cite [S.Rodini, L.Rossi, A.Vladimirov, [2404.01162] -",c_cyan)

write(*,*) "SnowFlake initlization with INI-file:"//trim(path)

!----------------- reading ini-file --------------------------------------
OPEN(UNIT=51, FILE=path, ACTION="read", STATUS="old")

!!! Number of grid-nodes in the sector-angle
call MoveTO(51,'0   :')
read(51,*) showINI

call MoveTO(51,'1   :')
read(51,*) showPROCESS

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

!!! Mass of BOTTOM threshold [GeV]
call MoveTO(51,'E.1 :')
read(51,*) stepT

CLOSE (51, STATUS='KEEP')

!$ call OMP_set_num_threads(allowedNumProcessor)

Qmin=1.
Qmax=2.
maxT=1

evolutionEvenPrepared=.false.
evolutionOddPrepared=.false.

if(showINI) then
    write(*,*) "SnowFlake parameters:"
    write(*,'("Grid size = (",I4," x",I4,") = ",I8, " nodes.")') NUM_RHO,NUM_PHI,NUM_TOT
    write(*,'(A)',advance="no") "Include Chiral-even evolution:"
    if(IncludeChiralEvenEvolution) then
        write(*,*) color(" YES",c_green)
    else
        write(*,*) color(" NO",c_red)
    end if
    write(*,'(A)',advance="no") "Include Chiral-odd evolution:"
    if(IncludeChiralOddEvolution) then
        write(*,*) color(" YES",c_green)
    else
        write(*,*) color(" NO",c_red)
    end if
    write(*,'(A)',advance="no") "Include gluon/flavor mixing :"
    if(useSingletEvolution) then
        write(*,*) color(" YES",c_green)
    else
        write(*,*) color(" NO",c_red)
    end if
    write(*,'(A,F8.3)') "Mass of CHARM threshold [GeV] :",massCHARM
    write(*,'(A,F8.3)') "Mass of BOTTOM threshold [GeV]:",massBOTTOM
end if

if(.not.(IncludeChiralEvenEvolution .or. IncludeChiralOddEvolution)) then
    ERROR STOP ErrorString("Non of evolution types are included. EVALUTION TERMINATED.","main")
end if

call Initialize_HexGrid(path)
call EvolutionKernels_Initialize(path)

call cpu_time(t2)
!$ t2=omp_get_wtime()
write(*,*) "--------- SnowFlake initialization complete. Initialization time =",t2-t1

end subroutine SnowFlake_Initialize


!!!! return the parameter t for the given Q
pure function tFromQ(Q)
real(dp)::tFromQ
real(dp),intent(in)::Q
tFromQ=2.d0*log(Q/Qmin)/stepT
end function tFromQ

!!!! return Q correpsonding to given t
pure function QFromt_real(t)
real(dp)::QFromt_real
real(dp),intent(in)::t
QFromt_real=Qmin*exp(t*stepT/2.d0)
end function QFromt_real

pure function QFromt_int(t)
real(dp)::QFromt_int
integer,intent(in)::t
QFromt_int=QFromt_real(real(t,dp))
end function QFromt_int

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! INTERFACES GLOBAL !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! compute the evolution of chiral-even distributions
!!! mu0, mu1 initial and final scales of evolution (GEV) [i.e. the grid is stored from mu0, till mu1+..]
!!! alpha = alpha_s(mu)
!!! inputQ,inputG defines the type of input functions
!!! --------
!!! inputQ = 'C'  C-definite input for quarks, i.e. Q1=\mathcal{S}+, Q2=\mathcal{S}-
!!! inputQ = 'S'  S-function input for quarks, i.e. Q1=S+, Q2=S-
!!! inputQ = 'T'  T-function input for quarks, i.e. Q1=T, Q2=DeltaT
!!! --------
!!! inputQ = 'C'  C-definite input for gluons, i.e. G1=\mathcal{G}+, G2=\mathcal{G}-
!!! inputQ = 'T'  T-function input for gluons, i.e. G1=T, G2=DeltaT
!!! --------
!!! the result of evolution is stored in the global variable and can be called separately.
subroutine ComputeEvolution(mu0,mu1,alpha,G1,U1,D1,S1,C1,B1,G2,U2,D2,S2,C2,B2,inputQ,inputG)
real(dp),external,optional::G1,U1,D1,S1,C1,B1,G2,U2,D2,S2,C2,B2
real(dp),external::alpha
real(dp),dimension(0:NUM_TOT)::G_p,U_p,D_p,S_p,C_p,B_p
real(dp),dimension(0:NUM_TOT)::G_m,U_m,D_m,S_m,C_m,B_m
real(dp),intent(in)::mu0,mu1
character(len=1),optional,intent(in)::inputQ,inputG
character(len=1)::inQ,inG

integer::n,i
real(dp)::x1,x2,x3,muT,muT1,t1,t2
!$ real*8::omp_get_wtime

if(present(inputQ)) then
    inQ=inputQ
else
    inQ='C'
end if
if(present(inputG)) then
    inG=inputG
else
    inG='C'
end if

!!!! if some function is absent it is interpreted as zero


!!! deshifrating input
SELECT CASE (inQ)
    CASE ('T')
        !!!! \mathcal{S}+=T(123)+T(321)-DeltaT(123)+DeltaT(321)
        !!!! \mathcal{S}-=T(123)-T(321)-DeltaT(123)-DeltaT(321)
        do n=0,NUM_TOT
        call get_X123_from_1Dindex(n,x1,x2,x3)
            !!!! U-quark
            if(present(U1)) then !!!! T
                U_p(n)=U1(x1,x2)+U1(x3,x2)
                U_m(n)=U1(x1,x2)-U1(x3,x2)
            else
                U_p(n)=0._dp
                U_m(n)=0._dp
            end if
            if(present(U2)) then !!!! DeltaT
                U_p(n)=U_p(n)-U2(x1,x2)+U2(x3,x2)
                U_m(n)=U_m(n)-U2(x1,x2)-U2(x3,x2)
            end if

            !!!! D-quark
            if(present(D1)) then !!!! T
                D_p(n)=D1(x1,x2)+D1(x3,x2)
                D_m(n)=D1(x1,x2)-D1(x3,x2)
            else
                D_p(n)=0._dp
                D_m(n)=0._dp
            end if
            if(present(D2)) then !!!! DeltaT
                D_p(n)=D_p(n)-D2(x1,x2)+D2(x3,x2)
                D_m(n)=D_m(n)-D2(x1,x2)-D2(x3,x2)
            end if

            !!!! S-quark
            if(present(S1)) then !!!! T
                S_p(n)=S1(x1,x2)+S1(x3,x2)
                S_m(n)=S1(x1,x2)-S1(x3,x2)
            else
                S_p(n)=0._dp
                S_m(n)=0._dp
            end if
            if(present(S2)) then !!!! DeltaT
                S_p(n)=S_p(n)-S2(x1,x2)+S2(x3,x2)
                S_m(n)=S_m(n)-S2(x1,x2)-S2(x3,x2)
            end if

            !!!! C-quark
            if(present(C1)) then !!!! T
                C_p(n)=C1(x1,x2)+C1(x3,x2)
                C_m(n)=C1(x1,x2)-C1(x3,x2)
            else
                C_p(n)=0._dp
                C_m(n)=0._dp
            end if
            if(present(C2)) then !!!! DeltaT
                C_p(n)=C_p(n)-C2(x1,x2)+C2(x3,x2)
                C_m(n)=C_m(n)-C2(x1,x2)-C2(x3,x2)
            end if

            !!!! B-quark
            if(present(B1)) then !!!! T
                B_p(n)=B1(x1,x2)+B1(x3,x2)
                B_m(n)=B1(x1,x2)-B1(x3,x2)
            else
                B_p(n)=0._dp
                B_m(n)=0._dp
            end if
            if(present(B2)) then !!!! DeltaT
                B_p(n)=B_p(n)-B2(x1,x2)+B2(x3,x2)
                B_m(n)=B_m(n)-B2(x1,x2)-B2(x3,x2)
            end if

        end do

    CASE ('S')
        !!!! \mathcal{S}+=-2(S+(123)+S-(321))
        !!!! \mathcal{S}-=-2(S+(123)-S-(321))
        do n=0,NUM_TOT
        call get_X123_from_1Dindex(n,x1,x2,x3)
            !!!! U-quark
            if(present(U1)) then !!!! S+
                U_p(n)=-2*U1(x1,x2)
                U_m(n)=-2*U1(x1,x2)
            else
                U_p(n)=0._dp
                U_m(n)=0._dp
            end if
            if(present(U2)) then !!!! S-
                U_p(n)=U_p(n)-2*U2(x3,x2)
                U_m(n)=U_m(n)+2*U2(x3,x2)
            end if

            !!!! D-quark
            if(present(D1)) then !!!! S+
                D_p(n)=-2*D1(x1,x2)
                D_m(n)=-2*D1(x1,x2)
            else
                D_p(n)=0._dp
                D_m(n)=0._dp
            end if
            if(present(D2)) then !!!! S-
                D_p(n)=D_p(n)-2*D2(x3,x2)
                D_m(n)=D_m(n)+2*D2(x3,x2)
            end if

            !!!! S-quark
            if(present(S1)) then !!!! S+
                S_p(n)=-2*S1(x1,x2)
                S_m(n)=-2*S1(x1,x2)
            else
                S_p(n)=0._dp
                S_m(n)=0._dp
            end if
            if(present(S2)) then !!!! S-
                S_p(n)=S_p(n)-2*S2(x3,x2)
                S_m(n)=S_m(n)+2*S2(x3,x2)
            end if

            !!!! C-quark
            if(present(C1)) then !!!! S+
                C_p(n)=-2*C1(x1,x2)
                C_m(n)=-2*C1(x1,x2)
            else
                C_p(n)=0._dp
                C_m(n)=0._dp
            end if
            if(present(C2)) then !!!! S-
                C_p(n)=C_p(n)-2*C2(x3,x2)
                C_m(n)=C_m(n)+2*C2(x3,x2)
            end if

            !!!! B-quark
            if(present(B1)) then !!!! S+
                B_p(n)=-2*B1(x1,x2)
                B_m(n)=-2*B1(x1,x2)
            else
                B_p(n)=0._dp
                B_m(n)=0._dp
            end if
            if(present(B2)) then !!!! S-
                B_p(n)=B_p(n)-2*B2(x3,x2)
                B_m(n)=B_m(n)+2*B2(x3,x2)
            end if

        end do

    CASE ('C')  !!! the definition is C-definite
        if(present(U1)) then
            U_p=GETgrid(U1)
        else
            U_p=0._dp
        end if
        if(present(D1)) then
            D_p=GETgrid(D1)
        else
            D_p=0._dp
        end if
        if(present(S1)) then
            S_p=GETgrid(S1)
        else
            S_p=0._dp
        end if
        if(present(C1)) then
            C_p=GETgrid(C1)
        else
            C_p=0._dp
        end if
        if(present(B1)) then
            B_p=GETgrid(B1)
        else
            B_p=0._dp
        end if
        !--------------
        if(present(U2)) then
            U_m=GETgrid(U2)
        else
            U_m=0._dp
        end if
        if(present(D2)) then
            D_m=GETgrid(D2)
        else
            D_m=0._dp
        end if
        if(present(S2)) then
            S_m=GETgrid(S2)
        else
            S_m=0._dp
        end if
        if(present(C2)) then
            C_m=GETgrid(C2)
        else
            C_m=0._dp
        end if
        if(present(B2)) then
            B_m=GETgrid(B2)
        else
            B_m=0._dp
        end if
    CASE DEFAULT
        write(*,*) ErrorString("unknown specification for inputQ. Evaluation STOP"," ")
        stop
END SELECT

!!!! if below thresholds nullify baoundary values of C and B
if(mu0<massCHARM) then
    C_p=0._dp
    C_m=0._dp
end if
if(mu0<massBOTTOM) then
    B_p=0._dp
    B_m=0._dp
end if

!!! deshifrating gluon input
SELECT CASE(inG)
    CASE ('T')
        !!! \mathcal{G}+=T+(123)-T+(132)+T+(213)
        !!! \mathcal{G}-=T-(123)+T-(132)-T-(213)
        do n=0,NUM_TOT
        call get_X123_from_1Dindex(n,x1,x2,x3)
            !!!! G+
            if(present(G1)) then !!!! T
                G_p(n)=G1(x1,x2)-G1(x1,x3)+G1(x2,x1)
            else
                G_p(n)=0._dp
            end if
            !!!! G-
            if(present(G2)) then !!!! T
                G_m(n)=G2(x1,x2)+G2(x1,x3)-G2(x2,x1)
            else
                G_m(n)=0._dp
            end if
        end do
    CASE ('C') !!! the definition is C-definite
        if(present(G1)) then
            G_p=GETgrid(G1)
        else
            G_p=0._dp
        end if
        if(present(G2)) then
            G_m=GETgrid(G2)
        else
            G_m=0._dp
        end if
    CASE DEFAULT
        write(*,*) ErrorString("unknown specification for inputG. Evaluation STOP"," ")
        stop
END SELECT

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!! start the grid computation!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Qmin=mu0
!!!  determining the values of grid
maxT=int(tFromQ(mu1))+1
!!! the upper bound is not exactly mu1, but a rounded value with respect to t
Qmax=QFromt(maxT)

if(showPROCESS) then
    write(*,'("Snowflake: start to compute grid from Q0=",F8.2," till Q1=",F8.2,". In total ", I6," nodes.")') Qmin,Qmax,maxT
end if

call cpu_time(t1)
!$ t1=omp_get_wtime()

!!!!! free memory
if(allocated(Uplus)) deallocate(Uplus)
if(allocated(Dplus)) deallocate(Dplus)
if(allocated(Splus)) deallocate(Splus)
if(allocated(Cplus)) deallocate(Cplus)
if(allocated(Bplus)) deallocate(Bplus)
if(allocated(Gplus)) deallocate(Gplus)

if(allocated(Uminus)) deallocate(Uminus)
if(allocated(Dminus)) deallocate(Dminus)
if(allocated(Sminus)) deallocate(Sminus)
if(allocated(Cminus)) deallocate(Cminus)
if(allocated(Bminus)) deallocate(Bminus)
if(allocated(Gminus)) deallocate(Gminus)

!!!!!!! allocate the space
allocate(Uplus(0:maxT,0:NUM_TOT))
allocate(Dplus(0:maxT,0:NUM_TOT))
allocate(Splus(0:maxT,0:NUM_TOT))
allocate(Cplus(0:maxT,0:NUM_TOT))
allocate(Bplus(0:maxT,0:NUM_TOT))
allocate(Gplus(0:maxT,0:NUM_TOT))

allocate(Uminus(0:maxT,0:NUM_TOT))
allocate(Dminus(0:maxT,0:NUM_TOT))
allocate(Sminus(0:maxT,0:NUM_TOT))
allocate(Cminus(0:maxT,0:NUM_TOT))
allocate(Bminus(0:maxT,0:NUM_TOT))
allocate(Gminus(0:maxT,0:NUM_TOT))

!!!!!!!! save initial values
Uplus(0,:)=U_p
Dplus(0,:)=D_p
Splus(0,:)=S_p
Cplus(0,:)=C_p
Bplus(0,:)=B_p
Gplus(0,:)=G_p

Uminus(0,:)=U_m
Dminus(0,:)=D_m
Sminus(0,:)=S_m
Cminus(0,:)=C_m
Bminus(0,:)=B_m
Gminus(0,:)=G_m

!!! finally call evolution, from step to step, storing at each iteration
!!! the variables G_p,U_p,D_p,S_p,C_p,B_p, are globally update each step
do i=1,maxT
    muT=QFromt(i-1) !!!! previous scale
    muT1=QFromt(i)    !!!! next scale

    call EvolvePLUS(alpha,muT,muT1,G_p,U_p,D_p,S_p,C_p,B_p)
    call EvolveMINUS(alpha,muT,muT1,G_m,U_m,D_m,S_m,C_m,B_m)

    !!!!!!!! save values each step
    Uplus(i,:)=U_p
    Dplus(i,:)=D_p
    Splus(i,:)=S_p
    Cplus(i,:)=C_p
    Bplus(i,:)=B_p
    Gplus(i,:)=G_p

    Uminus(i,:)=U_m
    Dminus(i,:)=D_m
    Sminus(i,:)=S_m
    Cminus(i,:)=C_m
    Bminus(i,:)=B_m
    Gminus(i,:)=G_m
end do

call cpu_time(t2)
!$ t2=omp_get_wtime()

if(showPROCESS) then
    write(*,'("Snowflake: evolution grid computed. Timing = ",F12.6," sec.")') t2-t1
end if

evolutionEvenPrepared=.true.

end subroutine ComputeEvolution

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! compute the evolution of chiral-odd distributions
!!! mu0, mu1 initial and final scales of evolution (GEV)
!!! alpha = alpha_s(mu)
!!! input are functions H, or E
!!! --------------------------
!!! the result of evolution is stored in the global variable and can be called separately.
subroutine ComputeEvolutionChiralOdd(mu0,mu1,alpha,U1,D1,S1,C1,B1)
real(dp),external,optional::U1,D1,S1,C1,B1
real(dp),external::alpha
real(dp),dimension(0:NUM_TOT)::U_p,D_p,S_p,C_p,B_p
real(dp),intent(in)::mu0,mu1

integer::i
real(dp)::t0,t1,muT,muT1,time1,time2
!$ real*8::omp_get_wtime

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!! start the grid computation!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

Qmin=mu0
!!!  determining the values of grid
maxT=int(tFromQ(mu1))+1
!!! the upper bound is not exactly mu1, but a rounded value with respect to t
Qmax=QFromt(maxT)

if(showPROCESS) then
    write(*,'("Snowflake: start to compute grid (chiral Odd) from Q0=",F8.2," till Q1=",F8.2,". In total ", I6," nodes.")') &
        Qmin,Qmax,maxT
end if

call cpu_time(time1)
!$ time1=omp_get_wtime()

!!!!! free memory
if(allocated(Uodd)) deallocate(Uodd)
if(allocated(Dodd)) deallocate(Dodd)
if(allocated(Sodd)) deallocate(Sodd)
if(allocated(Codd)) deallocate(Codd)
if(allocated(Bodd)) deallocate(Bodd)

!!!!!!! allocate the space
allocate(Uodd(0:maxT,0:NUM_TOT))
allocate(Dodd(0:maxT,0:NUM_TOT))
allocate(Sodd(0:maxT,0:NUM_TOT))
allocate(Codd(0:maxT,0:NUM_TOT))
allocate(Bodd(0:maxT,0:NUM_TOT))

!!! the chiral odd functions do not mix in flavor.
!!! Thus if something is zero, it remains zero. And there is no need to compute it

if(present(U1)) then
    U_p=GETgrid(U1)
    Uodd(0,:)=U_p

    do i=1,maxT
        muT=QFromt(i-1) !!!! previous scale
        muT1=QFromt(i)    !!!! next scale

        !!!! t=log[mu^2]
        t0=2*log(muT)
        t1=2*log(muT1)
        call EvChiralOdd(U_p,as,t0,t1)
        !!!!!!!! save values each step
        Uodd(i,:)=U_p
    end do
else

    Uodd(:,:)=0._dp
end if

if(present(D1)) then
    D_p=GETgrid(D1)
    Dodd(0,:)=D_p

    do i=1,maxT
        muT=QFromt(i-1) !!!! previous scale
        muT1=QFromt(i)    !!!! next scale

        !!!! t=log[mu^2]
        t0=2*log(muT)
        t1=2*log(muT1)
        call EvChiralOdd(D_p,as,t0,t1)
        !!!!!!!! save values each step
        Dodd(i,:)=D_p
    end do
else

    Dodd(:,:)=0._dp
end if

if(present(S1)) then
    S_p=GETgrid(S1)
    Sodd(0,:)=S_p

    do i=1,maxT
        muT=QFromt(i-1) !!!! previous scale
        muT1=QFromt(i)    !!!! next scale

        !!!! t=log[mu^2]
        t0=2*log(muT)
        t1=2*log(muT1)
        call EvChiralOdd(S_p,as,t0,t1)
        !!!!!!!! save values each step
        Sodd(i,:)=S_p
    end do
else

    Sodd(:,:)=0._dp
end if

if(present(C1)) then
    C_p=GETgrid(C1)
    Codd(0,:)=C_p

    do i=1,maxT
        muT=QFromt(i-1) !!!! previous scale
        muT1=QFromt(i)    !!!! next scale

        !!!! t=log[mu^2]
        t0=2*log(muT)
        t1=2*log(muT1)
        call EvChiralOdd(C_p,as,t0,t1)
        !!!!!!!! save values each step
        Codd(i,:)=C_p
    end do
else

    Codd(:,:)=0._dp
end if

if(present(B1)) then
    B_p=GETgrid(B1)
    Bodd(0,:)=B_p

    do i=1,maxT
        muT=QFromt(i-1) !!!! previous scale
        muT1=QFromt(i)    !!!! next scale

        !!!! t=log[mu^2]
        t0=2*log(muT)
        t1=2*log(muT1)
        call EvChiralOdd(B_p,as,t0,t1)
        !!!!!!!! save values each step
        Bodd(i,:)=B_p
    end do
else

    Bodd(:,:)=0._dp
end if

call cpu_time(time2)
!$ time2=omp_get_wtime()

if(showPROCESS) then
    write(*,'("Snowflake: grid computed. Timing = ",F12.6," sec.")') time2-time1
end if

evolutionOddPrepared=.true.

contains

    !!!!! as=alpha_s(mu)/(4pi)
    function as(t)
        real(dp)::as,t
        !!! 1/(4 pi)
        real(dp),parameter::pi4minus1=real(0.0795774715459476678844418816862571810,dp)
        as=alpha(Exp(t/2))*pi4minus1
    end function as

end subroutine ComputeEvolutionChiralOdd

!!!!! computed the values of tw3 PDF at point x1,x2, and scale Q
!!!!! The flavor f is defined as
!!!!! f=0=10=gluon
!!!!! f=1=d
!!!!! f=2=u
!!!!! f=3=s
!!!!! f=4=c
!!!!! f=5=b
!!!!! outputT specifies the type of requested function
function GetPDF(x1,x2,Q,f,outputT)
real(dp)::GetPDF
real(dp),intent(in)::x1,x2,Q
integer,intent(in)::f
character(len=1),optional,intent(in)::outputT
character(len=1)::outT
real(dp)::x3,r,t_here,f1,f2
integer::t_low

if(.not.evolutionEvenPrepared) then
    error stop "The evolution tables are not computed. Run ComputeEvolution"
end if

!!!!! processing input x's
x3=-x1-x2
r=max(abs(x1),abs(x2),abs(x1+x2))
if(r<xMIN) then
    write(*,*) ErrorString("requested value with |x|<xMIN"," ")
    write(*,*) "Requested (x1,x2,x3)=(",x1,x2,x3,") with |x|=",r
    write(*,*) "Current xMIN =",xMIN
    write(*,*) "Increase lower limit in INI-file. Evaluation STOP"
    stop
end if
if(abs(x1)>1.or.abs(x2)>1.or.abs(x1+x2)>1) then
    write(*,*) WarningString("requested values at |x|>1. Zero returned."," ")
GetPDF=0._dp
end if

!!!!! processing input Q
if(Q<Qmin) then
    write(*,*) WarningString("requested value with Q<QMIN. Returned linear extrapolation."," ")
end if

if(Q>Qmax) then
    write(*,*) WarningString("requested value with Q>QMAX. Returned linear extrapolation."," ")
end if

    !!!!! processing input type
if(present(outputT)) then
outT=outputT
else
outT='C'
end if

!!!!! these are the nodes
t_here=tFromQ(Q)
t_low=int(t_here)
!!!!! the computation is by linear interpolation
!!! F(t) = (t-t1)f2+(t2-t)f1 (with t2=t1+1

SELECT CASE (outT)
CASE ('C')
    if(f==10.or.f==0) then
        f1=GETinterpolation(x1,x2,Gplus(t_low,:))
        f2=GETinterpolation(x1,x2,Gplus(t_low+1,:))
    else if(f==1) then
        f1=GETinterpolation(x1,x2,Dplus(t_low,:))
        f2=GETinterpolation(x1,x2,Dplus(t_low+1,:))
    else if(f==2) then
        f1=GETinterpolation(x1,x2,Uplus(t_low,:))
        f2=GETinterpolation(x1,x2,Uplus(t_low+1,:))
    else if(f==3) then
        f1=GETinterpolation(x1,x2,Splus(t_low,:))
        f2=GETinterpolation(x1,x2,Splus(t_low+1,:))
    else if(f==4) then
        f1=GETinterpolation(x1,x2,Cplus(t_low,:))
        f2=GETinterpolation(x1,x2,Cplus(t_low+1,:))
    else if(f==5) then
        f1=GETinterpolation(x1,x2,Bplus(t_low,:))
        f2=GETinterpolation(x1,x2,Bplus(t_low+1,:))
    else if(f==-10) then
        f1=GETinterpolation(x1,x2,Gminus(t_low,:))
        f2=GETinterpolation(x1,x2,Gminus(t_low+1,:))
    else if(f==-1) then
        f1=GETinterpolation(x1,x2,Dminus(t_low,:))
        f2=GETinterpolation(x1,x2,Dminus(t_low+1,:))
    else if(f==-2) then
        f1=GETinterpolation(x1,x2,Uminus(t_low,:))
        f2=GETinterpolation(x1,x2,Uminus(t_low+1,:))
    else if(f==-3) then
        f1=GETinterpolation(x1,x2,Sminus(t_low,:))
        f2=GETinterpolation(x1,x2,Sminus(t_low+1,:))
    else if(f==-4) then
        f1=GETinterpolation(x1,x2,Cminus(t_low,:))
        f2=GETinterpolation(x1,x2,Cminus(t_low+1,:))
    else if(f==-5) then
        f1=GETinterpolation(x1,x2,Bminus(t_low,:))
        f2=GETinterpolation(x1,x2,Bminus(t_low+1,:))
    else
        write(*,*) WarningString("the flavor can be only -5,...,5, 10, -10. Zero returned."," ")
        GetPDF=0._dp
        return
    end if
CASE ('T')
    if(f==10.or.f==0) then
        f1=(GETinterpolation(x1,x2,Gplus(t_low,:))-GETinterpolation(x3,x2,Gplus(t_low,:)))/2
        f2=(GETinterpolation(x1,x2,Gplus(t_low+1,:))-GETinterpolation(x3,x2,Gplus(t_low+1,:)))/2
    else if(f==1) then
        f1=(GETinterpolation(x1,x2,Dplus(t_low,:))+GETinterpolation(x3,x2,Dplus(t_low,:)) &
            +GETinterpolation(x1,x2,Dminus(t_low,:))-GETinterpolation(x3,x2,Dminus(t_low,:)))/4
        f2=(GETinterpolation(x1,x2,Dplus(t_low+1,:))+GETinterpolation(x3,x2,Dplus(t_low+1,:))&
            +GETinterpolation(x1,x2,Dminus(t_low+1,:))-GETinterpolation(x3,x2,Dminus(t_low+1,:)))/4
    else if(f==2) then
        f1=(GETinterpolation(x1,x2,Uplus(t_low,:))+GETinterpolation(x3,x2,Uplus(t_low,:))&
            +GETinterpolation(x1,x2,Uminus(t_low,:))-GETinterpolation(x3,x2,Uminus(t_low,:)))/4
        f2=(GETinterpolation(x1,x2,Uplus(t_low+1,:))+GETinterpolation(x3,x2,Uplus(t_low+1,:))&
            +GETinterpolation(x1,x2,Uminus(t_low+1,:))-GETinterpolation(x3,x2,Uminus(t_low+1,:)))/4
    else if(f==3) then
        f1=(GETinterpolation(x1,x2,Splus(t_low,:))+GETinterpolation(x3,x2,Splus(t_low,:))&
            +GETinterpolation(x1,x2,Sminus(t_low,:))-GETinterpolation(x3,x2,Sminus(t_low,:)))/4
        f2=(GETinterpolation(x1,x2,Splus(t_low+1,:))+GETinterpolation(x3,x2,Splus(t_low+1,:))&
            +GETinterpolation(x1,x2,Sminus(t_low+1,:))-GETinterpolation(x3,x2,Sminus(t_low+1,:)))/4
    else if(f==4) then
        f1=(GETinterpolation(x1,x2,Cplus(t_low,:))+GETinterpolation(x3,x2,Cplus(t_low,:))&
            +GETinterpolation(x1,x2,Cminus(t_low,:))-GETinterpolation(x3,x2,Cminus(t_low,:)))/4
        f2=(GETinterpolation(x1,x2,Cplus(t_low+1,:))+GETinterpolation(x3,x2,Cplus(t_low+1,:))&
            +GETinterpolation(x1,x2,Cminus(t_low+1,:))-GETinterpolation(x3,x2,Cminus(t_low+1,:)))/4
    else if(f==5) then
        f1=(GETinterpolation(x1,x2,Bplus(t_low,:))+GETinterpolation(x3,x2,Bplus(t_low,:))&
            +GETinterpolation(x1,x2,Bminus(t_low,:))-GETinterpolation(x3,x2,Bminus(t_low,:)))/4
        f2=(GETinterpolation(x1,x2,Bplus(t_low+1,:))+GETinterpolation(x3,x2,Bplus(t_low+1,:))&
            +GETinterpolation(x1,x2,Bminus(t_low+1,:))-GETinterpolation(x3,x2,Bminus(t_low+1,:)))/4
    else if(f==-10) then
        f1=(GETinterpolation(x1,x2,Gminus(t_low,:))+GETinterpolation(x3,x2,Gminus(t_low,:)))/2
        f2=(GETinterpolation(x1,x2,Gminus(t_low+1,:))+GETinterpolation(x3,x2,Gminus(t_low+1,:)))/2
    else if(f==-1) then
        f1=-(GETinterpolation(x1,x2,Dplus(t_low,:))-GETinterpolation(x3,x2,Dplus(t_low,:))&
            +GETinterpolation(x1,x2,Dminus(t_low,:))+GETinterpolation(x3,x2,Dminus(t_low,:)))/4
        f2=-(GETinterpolation(x1,x2,Dplus(t_low+1,:))-GETinterpolation(x3,x2,Dplus(t_low+1,:))&
            +GETinterpolation(x1,x2,Dminus(t_low+1,:))+GETinterpolation(x3,x2,Dminus(t_low+1,:)))/4
    else if(f==-2) then
        f1=-(GETinterpolation(x1,x2,Uplus(t_low,:))-GETinterpolation(x3,x2,Uplus(t_low,:))&
            +GETinterpolation(x1,x2,Uminus(t_low,:))+GETinterpolation(x3,x2,Uminus(t_low,:)))/4
        f2=-(GETinterpolation(x1,x2,Uplus(t_low+1,:))-GETinterpolation(x3,x2,Uplus(t_low+1,:))&
            +GETinterpolation(x1,x2,Uminus(t_low+1,:))+GETinterpolation(x3,x2,Uminus(t_low+1,:)))/4
    else if(f==-3) then
        f1=-(GETinterpolation(x1,x2,Splus(t_low,:))-GETinterpolation(x3,x2,Splus(t_low,:))&
            +GETinterpolation(x1,x2,Sminus(t_low,:))+GETinterpolation(x3,x2,Sminus(t_low,:)))/4
        f2=-(GETinterpolation(x1,x2,Splus(t_low+1,:))-GETinterpolation(x3,x2,Splus(t_low+1,:))&
            +GETinterpolation(x1,x2,Sminus(t_low+1,:))+GETinterpolation(x3,x2,Sminus(t_low+1,:)))/4
    else if(f==-4) then
        f1=-(GETinterpolation(x1,x2,Cplus(t_low,:))-GETinterpolation(x3,x2,Cplus(t_low,:))&
            +GETinterpolation(x1,x2,Cminus(t_low,:))+GETinterpolation(x3,x2,Cminus(t_low,:)))/4
        f2=-(GETinterpolation(x1,x2,Cplus(t_low+1,:))-GETinterpolation(x3,x2,Cplus(t_low+1,:))&
            +GETinterpolation(x1,x2,Cminus(t_low+1,:))+GETinterpolation(x3,x2,Cminus(t_low+1,:)))/4
    else if(f==-5) then
        f1=-(GETinterpolation(x1,x2,Bplus(t_low,:))-GETinterpolation(x3,x2,Bplus(t_low,:))&
            +GETinterpolation(x1,x2,Bminus(t_low,:))+GETinterpolation(x3,x2,Bminus(t_low,:)))/4
        f2=-(GETinterpolation(x1,x2,Bplus(t_low+1,:))-GETinterpolation(x3,x2,Bplus(t_low+1,:))&
            +GETinterpolation(x1,x2,Bminus(t_low+1,:))+GETinterpolation(x3,x2,Bminus(t_low+1,:)))/4
    else
        write(*,*) WarningString("the flavor can be only -5,...,5, 10, -10. Zero returned."," ")
        GetPDF=0._dp
        return
    end if
    CASE ('S')
    if(f==10.or.f==0) then
        f1=(GETinterpolation(x1,x2,Gplus(t_low,:))-GETinterpolation(x3,x2,Gplus(t_low,:)))/2
        f1=(GETinterpolation(x1,x2,Gplus(t_low+1,:))-GETinterpolation(x3,x2,Gplus(t_low+1,:)))/2
    else if(f==1) then
        f1=-(GETinterpolation(x1,x2,Dplus(t_low,:))+GETinterpolation(x1,x2,Dminus(t_low,:)))/4
        f2=-(GETinterpolation(x1,x2,Dplus(t_low+1,:))+GETinterpolation(x1,x2,Dminus(t_low+1,:)))/4
    else if(f==2) then
        f1=-(GETinterpolation(x1,x2,Uplus(t_low,:))+GETinterpolation(x1,x2,Uminus(t_low,:)))/4
        f2=-(GETinterpolation(x1,x2,Uplus(t_low+1,:))+GETinterpolation(x1,x2,Uminus(t_low+1,:)))/4
    else if(f==3) then
        f1=-(GETinterpolation(x1,x2,Splus(t_low,:))+GETinterpolation(x1,x2,Sminus(t_low,:)))/4
        f2=-(GETinterpolation(x1,x2,Splus(t_low+1,:))+GETinterpolation(x1,x2,Sminus(t_low+1,:)))/4
    else if(f==4) then
        f1=-(GETinterpolation(x1,x2,Cplus(t_low,:))+GETinterpolation(x1,x2,Cminus(t_low,:)))/4
        f2=-(GETinterpolation(x1,x2,Cplus(t_low+1,:))+GETinterpolation(x1,x2,Cminus(t_low+1,:)))/4
    else if(f==5) then
        f1=-(GETinterpolation(x1,x2,Bplus(t_low,:))+GETinterpolation(x1,x2,Bminus(t_low,:)))/4
        f2=-(GETinterpolation(x1,x2,Bplus(t_low+1,:))+GETinterpolation(x1,x2,Bminus(t_low+1,:)))/4
    else if(f==-10) then
        f1=(GETinterpolation(x1,x2,Gminus(t_low,:))+GETinterpolation(x3,x2,Gminus(t_low,:)))/2
        f2=(GETinterpolation(x1,x2,Gminus(t_low+1,:))+GETinterpolation(x3,x2,Gminus(t_low+1,:)))/2
    else if(f==-1) then
        f1=-(GETinterpolation(x3,x2,Dplus(t_low,:))-GETinterpolation(x3,x2,Dminus(t_low,:)))/4
        f2=-(GETinterpolation(x3,x2,Dplus(t_low+1,:))-GETinterpolation(x3,x2,Dminus(t_low+1,:)))/4
    else if(f==-2) then
        f1=-(GETinterpolation(x3,x2,Uplus(t_low,:))-GETinterpolation(x3,x2,Uminus(t_low,:)))/4
        f2=-(GETinterpolation(x3,x2,Uplus(t_low+1,:))-GETinterpolation(x3,x2,Uminus(t_low+1,:)))/4
    else if(f==-3) then
        f1=-(GETinterpolation(x3,x2,Splus(t_low,:))-GETinterpolation(x3,x2,Sminus(t_low,:)))/4
        f2=-(GETinterpolation(x3,x2,Splus(t_low+1,:))-GETinterpolation(x3,x2,Sminus(t_low+1,:)))/4
    else if(f==-4) then
        f1=-(GETinterpolation(x3,x2,Cplus(t_low,:))-GETinterpolation(x3,x2,Cminus(t_low,:)))/4
        f2=-(GETinterpolation(x3,x2,Cplus(t_low+1,:))-GETinterpolation(x3,x2,Cminus(t_low+1,:)))/4
    else if(f==-5) then
        f1=-(GETinterpolation(x3,x2,Bplus(t_low,:))-GETinterpolation(x3,x2,Bminus(t_low,:)))/4
        f2=-(GETinterpolation(x3,x2,Bplus(t_low+1,:))-GETinterpolation(x3,x2,Bminus(t_low+1,:)))/4
    else
        write(*,*) WarningString("the flavor can be only -5,...,5, 10, -10. Zero returned."," ")
        GetPDF=0._dp
        return
    end if
CASE DEFAULT
    write(*,*) ErrorString("unknown outputT. Evaluation STOP."," ")
    stop
END SELECT

!!!!! interpolation
GetPDF=(t_here-t_low)*(f2-f1)+f1

if(ISNAN(GetPDF)) then
write(*,*) "------ NAN INSIDE THE snowflake GETPDF -----"
write(*,*) "x1, x2, Q, f, outputT ----> ",x1,x2,Q,f,outputT
write(*,*) "t_here,t_low, f2,f1 ---->", t_here,t_low, f2,f1
end if

end function GetPDF

!!!!! Computes interpolation flavor f for chiral odd-case
!!!!! f=1=d
!!!!! f=2=u
!!!!! f=3=s
!!!!! f=4=c
!!!!! f=5=b
!!!!!
function GetPDFChiralOdd(x1,x2,Q,f)
real(dp)::GetPDFChiralOdd
real(dp),intent(in)::x1,x2,Q
integer,intent(in)::f
real(dp)::x3,r,t_here,f1,f2
integer::t_low

if(.not.evolutionOddPrepared) then
    error stop "The evolution tables are not computed. Run ComputeEvolutionChiralOdd"
end if

x3=-x1-x2

r=max(abs(x1),abs(x2),abs(x1+x2))
if(r<xMIN) then
    write(*,*) ErrorString("requested value with |x|<xMIN"," ")
    write(*,*) "Requested (x1,x2,x3)=(",x1,x2,x3,") with |x|=",r
    write(*,*) "Current xMIN =",xMIN
    write(*,*) "Increase lower limit in INI-file. Evaluation STOP"
    stop
end if
if(abs(x1)>1.or.abs(x2)>1.or.abs(x1+x2)>1) then
    write(*,*) WarningString("requested values at |x|>1. Zero returned."," ")
GetPDFChiralOdd=0._dp
end if

!!!!! processing input Q
if(Q<Qmin) then
    write(*,*) WarningString("requested value with Q<QMIN. Returned linear extrapolation."," ")
end if

if(Q>Qmax) then
    write(*,*) WarningString("requested value with Q>QMAX. Returned linear extrapolation."," ")
end if

!!!!! these are the nodes
t_here=tFromQ(Q)
t_low=int(t_here)
!!!!! the computation is by linear interpolation
!!! F(t) = (t-t1)f2+(t2-t)f1 (with t2=t1+1



if(f==1) then
    f1=GETinterpolation(x1,x2,Dodd(t_low,:))
    f2=GETinterpolation(x1,x2,Dodd(t_low+1,:))
else if(f==2) then
    f1=GETinterpolation(x1,x2,Uodd(t_low,:))
    f2=GETinterpolation(x1,x2,Uodd(t_low+1,:))
else if(f==3) then
    f1=GETinterpolation(x1,x2,Sodd(t_low,:))
    f2=GETinterpolation(x1,x2,Sodd(t_low+1,:))
else if(f==4) then
    f1=GETinterpolation(x1,x2,Codd(t_low,:))
    f2=GETinterpolation(x1,x2,Codd(t_low+1,:))
else if(f==5) then
    f1=GETinterpolation(x1,x2,Bodd(t_low,:))
    f2=GETinterpolation(x1,x2,Bodd(t_low+1,:))
else
    write(*,*) WarningString("the flavor can be only -5,...,5, 10, -10. Zero returned."," ")
    GetPDFChiralOdd=0._dp
    return
end if

!!!!! interpolation
GetPDFChiralOdd=(t_here-t_low)*(f2-f1)+f1

end function GetPDFChiralOdd

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!! INTERFACES FOR EVOLUTION OF DITRIBUTIONS !!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! Evolution of PLUS components
!!! G,U,... = 1D flavor arrays to evolve
!!! alpha = alpha_s(mu)
!!! mu0, mu1 = initial/final scales in GeV
subroutine EvolvePLUS(alpha,mu0,mu1,G,U,D,S,C,B)
real(dp),external::alpha
real(dp),intent(in)::mu0,mu1
real(dp)::t0,t1
real(dp),dimension(0:NUM_TOT)::G,U,D,S,C,B
real(dp),allocatable,dimension(:)::Singlet,N1,N2,N3,N4!!! singlet, and 3 non-singlet combinations

    if(useSingletEvolution) then

    if(mu0>mu1) then
        write(*,*) ErrorString("evolution from high-to-low scale is not implemented yet"," ")
        write(*,*) "       computation terminated."
        stop
    end if

    !!!! depending on the range of mu one need to allocate different set of variables
    !!!! in anycase singlet Singlet,N1,N2
    allocate(Singlet(0:NUM_TOT))
    allocate(N1(0:NUM_TOT))
    allocate(N2(0:NUM_TOT))
    if(mu1>massCHARM) allocate(N3(0:NUM_TOT))
    if(mu1>massBOTTOM) allocate(N4(0:NUM_TOT))

!         N1=U
!         N2=U
!         call EvNonSinglet(N1,as,0.d0,1.d0)
!         call EvNonSingletOLD(N2,as,0.d0,1.d0)
!         write(*,*) N1-N2
!         stop



    !!!! The boundary condition of the singlet combination depend on the range of mu
    if(mu0<massCHARM) then
        Singlet=U+D+S
    else if(mu0<massBOTTOM) then
        Singlet=U+D+S+C
    else
        Singlet=U+D+S+C+B
    end if
    !!! boundary conditions for non-singlet
    N1=U-D
    N2=U+D-2*S

    !!!!! evolving the singlet+gluon part
    !!!!! depending on range of mu, there are several thresholds
    !!!!! Evolution through the threshold defines boundary for N3 and N4 (if necesary)
    if(mu1<massCHARM) then
        !!!! CASE :  mu0<mu1<massCHARM  Nf=3

        !!!! t=log[mu^2]
        t0=2*log(mu0)
        t1=2*log(mu1)
        call EvSingletPLUS(Singlet,G,as,t0,t1,3)

    else if(mu1<massBOTTOM) then

        if(mu0<massCHARM) then
            !!!! CASE :  mu0<massCharm<mu1<massBOTTOM  Nf=3, then Nf=4
            !!!! t=log[mu^2]
            t0=2*log(mu0)
            t1=2*log(massCHARM)
            call EvSingletPLUS(Singlet,G,as,t0,t1,3)

            N3=Singlet !!! boundary conidition for N3-combination

            t0=t1
            t1=2*log(mu1)
            call EvSingletPLUS(Singlet,G,as,t0,t1,4)
        else
            N3=U+D+S-3*C !!! boundary condition for N3-combination

            !!!! CASE :  massCharm<mu0<mu1<massBOTTOM  Nf=4
            !!!! t=log[mu^2]
            t0=2*log(mu0)
            t1=2*log(mu1)
            call EvSingletPLUS(Singlet,G,as,t0,t1,4)
        end if
    else
        if(mu0<massCHARM) then

            !!!! CASE :  mu0<massCharm<massBOTTOM<mu1  Nf=3, then Nf=4, then Nf=5
            !!!! t=log[mu^2]
            t0=2*log(mu0)
            t1=2*log(massCHARM)
            call EvSingletPLUS(Singlet,G,as,t0,t1,3)

            N3=Singlet !!! boundary conidition for N3-combination

            t0=t1
            t1=2*log(massBOTTOM)
            call EvSingletPLUS(Singlet,G,as,t0,t1,4)

            N4=Singlet !!! boundary conidition for N4-combination

            t0=t1
            t1=2*log(mu1)
            call EvSingletPLUS(Singlet,G,as,t0,t1,5)

        else if(mu0<massBOTTOM) then

            !!!! CASE :  massCharm<mu0<massBOTTOM<mu1  Nf=4 then Nf=5
            !!!! t=log[mu^2]
            t0=2*log(mu0)
            t1=2*log(massBOTTOM)
            call EvSingletPLUS(Singlet,G,as,t0,t1,4)

            N3=U+D+S-3*C !!! boundary conidition for N3-combination
            N4=Singlet !!! boundary conidition for N4-combination

            t0=t1
            t1=2*log(mu1)
            call EvSingletPLUS(Singlet,G,as,t0,t1,5)

        else
            N3=U+D+S-3*C !!! boundary conidition for N3-combination
            N4=U+D+S+C-4*B !!! boundary conidition for N4-combination

            !!!! CASE :  massCharm<massBOTTOM<mu0<mu1  Nf=5
            !!!! t=log[mu^2]
            t0=2*log(mu0)
            t1=2*log(mu1)
            call EvSingletPLUS(Singlet,G,as,t0,t1,5)
        end if
    end if
    !!! non-singletes evolve as they are

    !!!! t=log[mu^2]
    t0=2*log(mu0)
    t1=2*log(mu1)

    call EvNonSinglet(N1,as,t0,t1)
    call EvNonSinglet(N2,as,t0,t1)

    !!!! the evolution for N3, and N4 is used only above thresholds.
    if(mu1>massCHARM) then
        if(mu0<massCHARM) then
            t0=2*log(massCHARM)
        end if
        call EvNonSinglet(N3,as,t0,t1)
    end if

    if(mu1>massBOTTOM) then
        if(mu0<massBOTTOM) then
            t0=2*log(massBOTTOM)
        else
            t0=2*log(mu0)
        end if
        call EvNonSinglet(N4,as,t0,t1)
    end if

    !!! FINALY we uncover the quark contributions
    if(mu1<massCHARM) then
        U=(3*N1+N2+2*Singlet)/6
        D=(-3*N1+N2+2*Singlet)/6
        S=(-N2+Singlet)/3
        C=0._dp
        B=0._dp
    else if(mu1<massBOTTOM) then
        U=(6*N1+2*N2+N3+3*Singlet)/12
        D=(-6*N1+2*N2+N3+3*Singlet)/12
        S=(-4*N2+N3+3*Singlet)/12
        C=(-N3+Singlet)/4
        B=0._dp
    else
        U=(30*N1+10*N2+5*N3+3*n4+12*Singlet)/60
        D=(-30*N1+10*N2+5*N3+3*n4+12*Singlet)/60
        S=(-20*N2+5*N3+3*n4+12*Singlet)/60
        C=(-5*N3+n4+4*Singlet)/20
        B=(-N4+Singlet)/5
    end if

else

    !!! non-singlet evolution is Nf-independent
    !!! thus boundary as they are

    !!!! t=log[mu^2]
    t0=2*log(mu0)
    t1=2*log(mu1)

    call EvNonSinglet(U,as,t0,t1)
    call EvNonSinglet(D,as,t0,t1)
    call EvNonSinglet(S,as,t0,t1)
    call EvNonSinglet(C,as,t0,t1)
    call EvNonSinglet(B,as,t0,t1)

end if

contains

!!!!! as=alpha_s(mu)/(4pi)
function as(t)
    real(dp)::as,t
    !!! 1/(4 pi)
    real(dp),parameter::pi4minus1=real(0.0795774715459476678844418816862571810,dp)
    as=alpha(Exp(t/2))*pi4minus1
end function as

end subroutine EvolvePLUS

!!! Evolution of MINUS components
!!! G,U,... = 1D flavor arrays to evolve
!!! alpha = alpha_s(mu)
!!! mu0, mu1 = initial/final scales in GeV
subroutine EvolveMINUS(alpha,mu0,mu1,G,U,D,S,C,B)
real(dp),external::alpha
real(dp),intent(in)::mu0,mu1
real(dp)::t0,t1
real(dp),dimension(0:NUM_TOT)::G,U,D,S,C,B
real(dp),allocatable,dimension(:)::Singlet,N1,N2,N3,N4!!! singlet, and 3 non-singlet combinations

if(useSingletEvolution) then
    !!! singlet evolution

    if(mu0>mu1) then
        write(*,*) ErrorString("evolution from high-to-low scale is not implemented yet"," ")
        write(*,*) "       computation terminated."
        stop
    end if

    !!!! depending on the range of mu one need to allocate different set of variables
    !!!! in anycase singlet Singlet,N1,N2
    allocate(Singlet(0:NUM_TOT))
    allocate(N1(0:NUM_TOT))
    allocate(N2(0:NUM_TOT))
    if(mu1>massCHARM) allocate(N3(0:NUM_TOT))
    if(mu1>massBOTTOM) allocate(N4(0:NUM_TOT))


    !!!! The boundary condition of the singlet combination depend on the range of mu
    if(mu0<massCHARM) then
        Singlet=U+D+S
    else if(mu0<massBOTTOM) then
        Singlet=U+D+S+C
    else
        Singlet=U+D+S+C+B
    end if
    !!! boundary conditions for non-singlet
    N1=U-D
    N2=U+D-2*S

    !!!!! evolving the singlet+gluon part
    !!!!! depending on range of mu, there are several thresholds
    !!!!! Evolution through the threshold defines boundary for N3 and N4 (if necesary)
    if(mu1<massCHARM) then
        !!!! CASE :  mu0<mu1<massCHARM  Nf=3

        !!!! t=log[mu^2]
        t0=2*log(mu0)
        t1=2*log(mu1)
        call EvSingletMINUS(Singlet,G,as,t0,t1,3)

    else if(mu1<massBOTTOM) then

        if(mu0<massCHARM) then
            !!!! CASE :  mu0<massCharm<mu1<massBOTTOM  Nf=3, then Nf=4
            !!!! t=log[mu^2]
            t0=2*log(mu0)
            t1=2*log(massCHARM)
            call EvSingletMINUS(Singlet,G,as,t0,t1,3)

            N3=Singlet !!! boundary conidition for N3-combination

            t0=t1
            t1=2*log(mu1)
            call EvSingletMINUS(Singlet,G,as,t0,t1,4)
        else
            N3=U+D+S-3*C !!! boundary condition for N3-combination

            !!!! CASE :  massCharm<mu0<mu1<massBOTTOM  Nf=4
            !!!! t=log[mu^2]
            t0=2*log(mu0)
            t1=2*log(mu1)
            call EvSingletMINUS(Singlet,G,as,t0,t1,4)
        end if
    else
        if(mu0<massCHARM) then
            !!!! CASE :  mu0<massCharm<massBOTTOM<mu1  Nf=3, then Nf=4, then Nf=5
            !!!! t=log[mu^2]
            t0=2*log(mu0)
            t1=2*log(massCHARM)
            call EvSingletMINUS(Singlet,G,as,t0,t1,3)

            N3=Singlet !!! boundary conidition for N3-combination

            t0=t1
            t1=2*log(massBOTTOM)
            call EvSingletMINUS(Singlet,G,as,t0,t1,4)

            N4=Singlet !!! boundary conidition for N4-combination

            t0=t1
            t1=2*log(mu1)
            call EvSingletMINUS(Singlet,G,as,t0,t1,5)

        else if(mu0<massBOTTOM) then

            !!!! CASE :  massCharm<mu0<massBOTTOM<mu1  Nf=4 then Nf=5
            !!!! t=log[mu^2]
            t0=2*log(mu0)
            t1=2*log(massBOTTOM)
            call EvSingletMINUS(Singlet,G,as,t0,t1,4)

            N3=U+D+S-3*C !!! boundary conidition for N3-combination
            N4=Singlet !!! boundary conidition for N4-combination

            t0=t1
            t1=2*log(mu1)
            call EvSingletMINUS(Singlet,G,as,t0,t1,5)

        else
            N3=U+D+S-3*C !!! boundary conidition for N3-combination
            N4=U+D+S+C-4*B !!! boundary conidition for N4-combination

            !!!! CASE :  massCharm<massBOTTOM<mu0<mu1  Nf=5
            !!!! t=log[mu^2]
            t0=2*log(mu0)
            t1=2*log(mu1)
            call EvSingletMINUS(Singlet,G,as,t0,t1,5)
        end if
    end if

    !!! non-singletes evolve as they are

    !!!! t=log[mu^2]
    t0=2*log(mu0)
    t1=2*log(mu1)
    call EvNonSinglet(N1,as,t0,t1)
    call EvNonSinglet(N2,as,t0,t1)

    !!!! the evolution for N3, and N4 is used only above thresholds.
    if(mu1>massCHARM) then
        if(mu0<massCHARM) then
            t0=2*log(massCHARM)
        end if

        call EvNonSinglet(N3,as,t0,t1)
    end if

    if(mu1>massBOTTOM) then
        if(mu0<massBOTTOM) then
            t0=2*log(massBOTTOM)
        else
            t0=2*log(mu0)
        end if
        call EvNonSinglet(N4,as,t0,t1)
    end if

    !!! FINALY we uncover the quark contributions
    if(mu1<massCHARM) then
        U=(3*N1+N2+2*Singlet)/6
        D=(-3*N1+N2+2*Singlet)/6
        S=(-N2+Singlet)/3
        C=0._dp
        B=0._dp
    else if(mu1<massBOTTOM) then
        U=(6*N1+2*N2+N3+3*Singlet)/12
        D=(-6*N1+2*N2+N3+3*Singlet)/12
        S=(-4*N2+N3+3*Singlet)/12
        C=(-N3+Singlet)/4
        B=0._dp
    else
        U=(30*N1+10*N2+5*N3+3*N4+12*Singlet)/60
        D=(-30*N1+10*N2+5*N3+3*N4+12*Singlet)/60
        S=(-20*N2+5*N3+3*N4+12*Singlet)/60
        C=(-5*N3+N4+4*Singlet)/20
        B=(-N4+Singlet)/5
    end if

else

    !!! non-singlet evolution is Nf-independent
    !!! thus boundary as they are

    !!!! t=log[mu^2]
    t0=2*log(mu0)
    t1=2*log(mu1)

    call EvNonSinglet(U,as,t0,t1)
    call EvNonSinglet(D,as,t0,t1)
    call EvNonSinglet(S,as,t0,t1)
    call EvNonSinglet(C,as,t0,t1)
    call EvNonSinglet(B,as,t0,t1)

end if


contains

!!!!! as=alpha_s(mu)/(4pi)
function as(t)
    real(dp)::as,t
    !!! 1/(4 pi)
    real(dp),parameter::pi4minus1=real(0.0795774715459476678844418816862571810,dp)
    as=alpha(Exp(t/2))*pi4minus1
end function as

end subroutine EvolveMINUS

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! G2 function !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! compute the G2 function from the stored projection matrix
!!!! The definition of G2 matrix already includes 1/2, so it is should be just applied to \mathfrak{S}^+
function G2(x,Q,f)
real(dp)::G2
real(dp),intent(in)::x,Q
integer,intent(in)::f

real(dp)::t_here
integer::t_low,t_low2
real(dp),dimension(0:NUM_TOT)::functionD1,functionD2,functionD

if(.not.evolutionEvenPrepared) then
    error stop "The evolution tables are not computed. Run ComputeEvolution"
end if

if(x>1.d0 .or. x<xMIN) then
    error stop "ERROR in G2: the value of x is outside of the allowed range [xMin,1]"
end if

!!!!! processing input Q
if(Q<Qmin) then
    write(*,*) WarningString("requested value with Q<QMIN. Returned linear extrapolation."," ")
end if

if(Q>Qmax) then
    write(*,*) WarningString("requested value with Q>QMAX. Returned linear extrapolation."," ")
end if

!!!!! these are the nodes for Q
t_here=tFromQ(Q)
t_low=int(t_here)
t_low2=t_low+1

SELECT CASE(f)
    CASE(1)
        functionD1=Dplus(t_low,:)
        functionD2=Dplus(t_low2,:)
    CASE(2)
        functionD1=Uplus(t_low,:)
        functionD2=Uplus(t_low2,:)
    CASE(3)
        functionD1=Splus(t_low,:)
        functionD2=Splus(t_low2,:)
    CASE(4)
        functionD1=Cplus(t_low,:)
        functionD2=Cplus(t_low2,:)
    CASE(5)
        functionD1=Bplus(t_low,:)
        functionD2=Bplus(t_low2,:)
    CASE(11) !!! u-d
        functionD1=Uplus(t_low,:)-Dplus(t_low,:)
        functionD2=Uplus(t_low2,:)-Dplus(t_low2,:)
    CASE(12) !!! u+d
        functionD1=Uplus(t_low,:)+Dplus(t_low,:)
        functionD2=Uplus(t_low2,:)+Dplus(t_low2,:)
    CASE(100)!!! proton
        !!!!!! the D=eq^2/2*(Sp-Sm)
        functionD1=(Uplus(t_low,:)+Cplus(t_low,:))*4/9 +(Dplus(t_low,:)+Splus(t_low,:)+Bplus(t_low,:))/9
        functionD2=(Uplus(t_low2,:)+Cplus(t_low2,:))*4/9 +(Dplus(t_low2,:)+Splus(t_low2,:)+Bplus(t_low2,:))/9
    CASE(101)!!! neutron = proton[u<->d]
        functionD1=(Dplus(t_low,:)+Cplus(t_low,:))*4/9 +(Uplus(t_low,:)+Splus(t_low,:)+Bplus(t_low,:))/9
        functionD2=(Dplus(t_low2,:)+Cplus(t_low2,:))*4/9 +(Uplus(t_low2,:)+Splus(t_low2,:)+Bplus(t_low2,:))/9
    CASE(102)!!! deutron=(p+n)/2
        functionD1=(Uplus(t_low,:)+Cplus(t_low,:)+Dplus(t_low,:)+Splus(t_low,:)+Bplus(t_low,:))*5/18
        functionD2=(Uplus(t_low2,:)+Cplus(t_low2,:)+Dplus(t_low2,:)+Splus(t_low2,:)+Bplus(t_low2,:))*5/18
    CASE DEFAULT
        write(*,*) ErrorString("D2 routine: unknown flavor"," ")
        write(*,*) "f=",f
        error stop
END SELECT

functionD=(t_low2-t_here)*functionD1+(t_here-t_low)*functionD2
!
! if(t_low==0) then
!     write(*,'("{",F8.6,",",F12.8,"},")') x, GETinterpolation(x,0.0d0,functionD)
! end if

G2=G2xF(x,functionD)

end function G2

!!!!! return the function of G2 computed over the parralel list
subroutine G2_List(res,x,Q,f)
real(dp),dimension(1:),intent(in)::x
real(dp),dimension(1:),intent(in)::Q
integer,dimension(1:),intent(in)::f
real(dp),dimension(1:),intent(out)::res

integer::i,length
length=size(x)
if(size(res)/=length) error stop ErrorString('sizes of output and X lists are not equal.',"Snoflake")
if(size(Q)/=length) error stop ErrorString('sizes of Q and X lists are not equal.',"Snoflake")
if(size(f)/=length) error stop ErrorString('sizes of f and X lists are not equal.',"Snoflake")

!$OMP PARALLEL DO DEFAULT(SHARED)
do i=1,length
    res(i)=G2(x(i),Q(i),f(i))
end do
!$OMP END PARALLEL DO

end subroutine G2_List

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! G2 function !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! compute the D2 Moment from the stored projection matrix
!!!! f=flavor, it. could be 1:5 (for individial flavors; gluon-undefined) or 100=proton, 101=neutron [u<->d]
!!!! The function d2 is 3\int dx x^2 g_2 (the factor 3 is ake ninto account in the matrix)
function D2(Q,f)
real(dp)::D2
integer,intent(in)::f
real(dp),intent(in)::Q

real(dp)::t_here
integer::t_low,t_low2
real(dp),dimension(0:NUM_TOT)::functionD1,functionD2,functionD

if(.not.evolutionEvenPrepared) then
    error stop "The evolution tables are not computed. Run ComputeEvolution"
end if

!!!!! processing input Q
if(Q<Qmin) then
    write(*,*) WarningString("requested value with Q<QMIN. Returned linear extrapolation."," ")
end if

if(Q>Qmax) then
    write(*,*) WarningString("requested value with Q>QMAX. Returned linear extrapolation."," ")
end if

!!!!! these are the nodes for Q
t_here=tFromQ(Q)
t_low=int(t_here)
t_low2=t_low+1


SELECT CASE(f)
    CASE(1)
        functionD1=Dplus(t_low,:)
        functionD2=Dplus(t_low2,:)
    CASE(2)
        functionD1=Uplus(t_low,:)
        functionD2=Uplus(t_low2,:)
    CASE(3)
        functionD1=Splus(t_low,:)
        functionD2=Splus(t_low2,:)
    CASE(4)
        functionD1=Cplus(t_low,:)
        functionD2=Cplus(t_low2,:)
    CASE(5)
        functionD1=Bplus(t_low,:)
        functionD2=Bplus(t_low2,:)
    CASE(11) !!! u-d
        functionD1=Uplus(t_low,:)-Dplus(t_low,:)
        functionD2=Uplus(t_low2,:)-Dplus(t_low2,:)
    CASE(12) !!! u+d
        functionD1=Uplus(t_low,:)+Dplus(t_low,:)
        functionD2=Uplus(t_low2,:)+Dplus(t_low2,:)
    CASE(100)!!! proton
        !!!!!! the D=eq^2/2*(Sp-Sm)
        functionD1=(Uplus(t_low,:)+Cplus(t_low,:))*4/9 +(Dplus(t_low,:)+Splus(t_low,:)+Bplus(t_low,:))/9
        functionD2=(Uplus(t_low2,:)+Cplus(t_low2,:))*4/9 +(Dplus(t_low2,:)+Splus(t_low2,:)+Bplus(t_low2,:))/9
    CASE(101)!!! neutron = proton[u<->d]
        functionD1=(Dplus(t_low,:)+Cplus(t_low,:))*4/9 +(Uplus(t_low,:)+Splus(t_low,:)+Bplus(t_low,:))/9
        functionD2=(Dplus(t_low2,:)+Cplus(t_low2,:))*4/9 +(Uplus(t_low2,:)+Splus(t_low2,:)+Bplus(t_low2,:))/9
    CASE(102)!!! deutron=(p+n)/2
        functionD1=(Uplus(t_low,:)+Cplus(t_low,:)+Dplus(t_low,:)+Splus(t_low,:)+Bplus(t_low,:))*5/18
        functionD2=(Uplus(t_low2,:)+Cplus(t_low2,:)+Dplus(t_low2,:)+Splus(t_low2,:)+Bplus(t_low2,:))*5/18
    CASE DEFAULT
        write(*,*) ErrorString("D2 routine: unknown flavor"," ")
        write(*,*) "f=",f
        error stop
END SELECT

functionD=(t_low2-t_here)*functionD1+(t_here-t_low)*functionD2

D2=D2xF(functionD)

end function D2

!!!!! return the function of G2 computed over the parralel list
subroutine D2_List(res,Q,f)
real(dp),dimension(1:),intent(in)::Q
integer,dimension(1:),intent(in)::f
real(dp),dimension(1:),intent(out)::res

integer::i,length
length=size(Q)
if(size(res)/=length) error stop ErrorString('sizes of output and Q lists are not equal.',"Snoflake.D2")
if(size(f)/=length) error stop ErrorString('sizes of f and Q lists are not equal.',"Snoflake.D2")

!$OMP PARALLEL DO DEFAULT(SHARED)
do i=1,length
    res(i)=D2(Q(i),f(i))
end do
!$OMP END PARALLEL DO

end subroutine D2_List

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! WGT function !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!! compute the WGT function from the stored projection matrix
!!!! The definition of WGT matrix already includes all factor, so it is should be just applied to S^-
!!!! it is defined only for quakr flavors
function WGT(x,Q,f)
real(dp)::WGT
real(dp),intent(in)::x,Q
integer,intent(in)::f

real(dp)::t_here
integer::t_low,t_low2
real(dp),dimension(0:NUM_TOT)::functionD1,functionD2,functionD

if(.not.evolutionEvenPrepared) then
    error stop "The evolution tables are not computed. Run ComputeEvolution"
end if

if(abs(x)>1.d0 .or. abs(x)<xMIN) then
    error stop "ERROR in WGT: the value of x is outside of the allowed range [-1,xMin] [xMin,1]"
end if

!!!!! processing input Q
if(Q<Qmin) then
    write(*,*) WarningString("requested value with Q<QMIN. Returned linear extrapolation."," ")
end if

if(Q>Qmax) then
    write(*,*) WarningString("requested value with Q>QMAX. Returned linear extrapolation."," ")
end if

!!!!! these are the nodes for Q
t_here=tFromQ(Q)
t_low=int(t_here)
t_low2=t_low+1

!!!! S^+=-(\frak{S}^++\frak{S}^-)/4

SELECT CASE(f)
    CASE(1)
        functionD1=-(Dplus(t_low,:)+Dminus(t_low,:))/4
        functionD2=-(Dplus(t_low2,:)+Dminus(t_low2,:))/4
    CASE(2)
        functionD1=-(Uplus(t_low,:)+Uminus(t_low,:))/4
        functionD2=-(Uplus(t_low2,:)+Uminus(t_low2,:))/4
    CASE(3)
        functionD1=-(Splus(t_low,:)+Sminus(t_low,:))/4
        functionD2=-(Splus(t_low2,:)+Sminus(t_low2,:))/4
    CASE(4)
        functionD1=-(Cplus(t_low,:)+Cminus(t_low,:))/4
        functionD2=-(Cplus(t_low2,:)+Cminus(t_low2,:))/4
    CASE(5)
        functionD1=-(Bplus(t_low,:)+Bminus(t_low,:))/4
        functionD2=-(Bplus(t_low2,:)+Bminus(t_low2,:))/4
    CASE DEFAULT
        write(*,*) ErrorString("D2 routine: unknown flavor"," ")
        write(*,*) "f=",f
        error stop
END SELECT

functionD=(t_low2-t_here)*functionD1+(t_here-t_low)*functionD2
!
! if(t_low==0) then
!     write(*,'("{",F8.6,",",F12.8,"},")') x, GETinterpolation(x,0.0d0,functionD)
! end if

WGT=WGTxF(x,functionD)

end function WGT

!!!! same as above, but for flavor list 2-times faster since there is no need to sum matirces 2 times
function WGT_fList(x,Q)
real(dp),dimension(-5:5)::WGT_fList
real(dp),intent(in)::x,Q

real(dp)::t_here
integer::t_low,t_low2
real(dp),dimension(0:NUM_TOT)::functionD1,functionD2,functionD
real(dp)::u,ub,d,db,s,sb,c,cb,b,bb

if(.not.evolutionEvenPrepared) then
    error stop "The evolution tables are not computed. Run ComputeEvolution"
end if

if(x>1.d0 .or. x<xMIN) then
    error stop "ERROR in WGT_fList: the value of x is outside of the allowed range [xMin,1]"
end if

!!!!! processing input Q
if(Q<Qmin) then
    write(*,*) WarningString("requested value with Q<QMIN. Returned linear extrapolation."," ")
end if

if(Q>Qmax) then
    write(*,*) WarningString("requested value with Q>QMAX. Returned linear extrapolation."," ")
end if

!!!!! these are the nodes for Q
t_here=tFromQ(Q)
t_low=int(t_here)
t_low2=t_low+1

!!!! S^+=-(\frak{S}^++\frak{S}^-)/4


functionD1=-(Dplus(t_low,:)+Dminus(t_low,:))/4
functionD2=-(Dplus(t_low2,:)+Dminus(t_low2,:))/4
functionD=(t_low2-t_here)*functionD1+(t_here-t_low)*functionD2
d=WGTxF(x,functionD)
db=WGTxF(-x,functionD)

functionD1=-(Uplus(t_low,:)+Uminus(t_low,:))/4
functionD2=-(Uplus(t_low2,:)+Uminus(t_low2,:))/4
functionD=(t_low2-t_here)*functionD1+(t_here-t_low)*functionD2
u=WGTxF(x,functionD)
ub=WGTxF(-x,functionD)

functionD1=-(Splus(t_low,:)+Sminus(t_low,:))/4
functionD2=-(Splus(t_low2,:)+Sminus(t_low2,:))/4
functionD=(t_low2-t_here)*functionD1+(t_here-t_low)*functionD2
s=WGTxF(x,functionD)
sb=WGTxF(-x,functionD)

functionD1=-(Cplus(t_low,:)+Cminus(t_low,:))/4
functionD2=-(Cplus(t_low2,:)+Cminus(t_low2,:))/4
functionD=(t_low2-t_here)*functionD1+(t_here-t_low)*functionD2
c=WGTxF(x,functionD)
cb=WGTxF(-x,functionD)

functionD1=-(Bplus(t_low,:)+Bminus(t_low,:))/4
functionD2=-(Bplus(t_low2,:)+Bminus(t_low2,:))/4
functionD=(t_low2-t_here)*functionD1+(t_here-t_low)*functionD2
b=WGTxF(x,functionD)
bb=WGTxF(-x,functionD)

!!!! the minus sign is because the TMD WGT function for anti-quark is defined as -f(-x)

WGT_fList=(/-bb,-cb,-sb,-ub,-db,0._dp,d,u,s,c,b/)

end function WGT_fList

end module SnowFlake
