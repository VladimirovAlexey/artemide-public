!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!                                                                                          !!
!! This file is a part of module EvolutionKernels, it contains                              !!
!! expressions for contraction matrix D2                                                    !!
!!                                                                                          !!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!! the matrix D2matrix is stored globally
!!!!!! it is matrix (0:NUM_RHO)

!!!! This subroutine computes the kernel-matrix for D2-moment
!!!! and stores it into the matrix M.
!!!! The computation is done for the already computed matrix G2matrix
!!!! the definition used here 3*\int_0^1 dx x^2 g2[x]
subroutine PreComputeMatrixD2(M,G2matrix)
real(dp), dimension(0:NUM_TOT)::M
real(dp), dimension(0:NUM_RHO,0:NUM_TOT), intent(in)::G2matrix
real(dp), dimension(0:NUM_RHO)::weights

integer::n

M=0._dp
!!!!! compute weight for each matris
do n=0,NUM_RHO
    weights(n)=Integrate_GK(f_d2,xMin,1._dp)
    M(:)=M(:)+weights(n)*G2matrix(n,:)
end do

contains
function f_d2(z3)
real(dp)::f_d2
real(dp),intent(in)::z3

f_d2=3._dp*GETinterpolatorRHO_X(n,z3)*z3**2

end function f_d2

end subroutine PreComputeMatrixD2

!!!! multiplies the matrix of G2``vector'' F (0... NUM_TOT)
!!!! and interpolate it to the point x
function D2xF(F)
real(dp)::D2xF
real(dp),dimension(0:NUM_TOT),intent(in)::F

!!!! getting the grid of F
D2xF=sum(D2matrix*F)
end function D2xF
