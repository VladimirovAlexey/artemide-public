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
!	This part is devoted to the calculation of Mellin convolution
!
!				A.Vladimirov (08.10.2018)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  !these are globals for the exchange of parameters between functions
  ! it is used only in the calculation without grid!!!
  real*8:: xCurrent,bTcurrent,muCurrent
  real*8::muAt1!!!special case for z-dependant mu. it is = mu(z=1)
  real*8,dimension(-5:5)::FPDFcurrent !!!value of fNP(1)*f(xCurrent)
  real*8,dimension(-5:5)::PDFcurrent !!!value of f(xCurrent)
  real*8,dimension(-5:5)::Fcurrent !!!value of fNP(1)
  real*8,dimension(-5:5) :: integralWeight
  !!!triger for calculation
  logical:: IsMuXdependent, IsFnpZdependent