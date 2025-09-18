!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                          !!
!! Code that save the kernels for future use.      A.Vladimirov 01.04.2024                  !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
program SaveKernelsToFile
use SnowFlake
use EvolutionKernels

implicit none

real*8::t1,t2
!$ real*8::omp_get_wtime

!------------------------------------------------------
! ! !!!! initialize the Snowflake with particular ini file
call  SnowFlake_Initialize("kernels16x8.ini","Prog_snowflake/")

call cpu_time(t1)
!$ t1=omp_get_wtime()
call SaveKernels("kernels/16x8/")


! call ReadKernels("kernels/16x8/")

! !------------------------------------------------------
! ! initialize the Snowflake with particular ini file
! call  SnowFlake_Initialize("kernels24x12.ini","Prog_snowflake/")
!
! call cpu_time(t1)
! !$ t1=omp_get_wtime()
! call SaveKernels("kernels/24x12/")
!
!
! call ReadKernels("kernels/24x12/")

!------------------------------------------------------
! ! ! initialize the Snowflake with particular ini file
! call  SnowFlake_Initialize("kernels25x20.ini","Prog_snowflake/")
!
! call cpu_time(t1)
! !$ t1=omp_get_wtime()
! call SaveKernels("kernels/25x20/")
!
!
! call ReadKernels("kernels/25x20/")
! !------------------------------------------------------

call cpu_time(t2)
!$ t2=omp_get_wtime()
write(*,*) "Time for computation of evolution",t2-t1

end program SaveKernelsToFile
