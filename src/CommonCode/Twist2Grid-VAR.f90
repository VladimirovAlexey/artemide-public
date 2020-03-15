!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 1.4
!
!	This file contains the part of the code, which is common for all TMD-evaluation modules that operates at twist-2
!	ONLY DECLARATION OF VARIABLES that are used in twist2Convolution
!	It shares the variables with the module where it is inlucded (as a text)
!	Such idiotic structure is needed since, FORTRAN does not allow inheritance.
!
!	Be AWARE of possible clash of variable names.
!
!	This part is devoted to the Grid evaluation
!
!				A.Vladimirov (08.10.2018)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  real(dp)::xGrid_Min=0.00001d0 !!!!min x in grid. if lower requared remake it
  real(dp)::bGrid_Max=100d0 !!!max b in grid if higher requared, extrapolation is done
  integer::GridSizeX=250
  integer::GridSizeB=750
  real(dp) :: slope=10d0
  real(dp), dimension(:,:,:,:), allocatable :: gridMain !!!! THIS IS HUGE(!) matrix for the grid
  real(dp), dimension(:,:,:), allocatable :: boundaryValues !!!! for b>bGrid_Max we approximate by fNP. With prefactor given by boundary value/fNP(bMax)
  real(dp):: bFactor,slopeFactor !! some often used combinations
