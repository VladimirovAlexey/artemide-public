!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Model for unpolarized TMD evolution for SV19
!
!			corresponds to bb* model
!			DNP=Dpert(b*)+g bb*
!			zeta=zetaPert(b) exp[-b2/BB]+zetaSL(b)(1-exp(-b2/BB)
!
!			Requres two NP parameters (initated by best values)
!
!				A.Vladimirov (11.07.2019)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module TMDR_model
use aTMDe_Numerics
use IO_functions
use TMD_AD, only : Dresum,Dpert,zetaMUpert,zetaSL,RADEvolution
implicit none

private

!!!!!------------------------------------------------------------------------------------
!!!!! These functions MUST defined in module  !!
!!!!!
!!!!! 1) The subroutine is called during the initialization of TMDR
!!!!!    arg=array of initial NP-parameters
public:: ModelInitialization
!!!!! 2) The subroutine that is called on reset of NP-parameters in TMDR
!!!!!    arg=array of new NP-parameters
public:: ModelUpdate
!!!!! 3) Function which returns RAD function
!!!!!    arg=(mu,b,f) with mu=scale(real_dp), b=transverse distance(real_dp), f=flavor(integer)
real(dp),public:: DNP
!!!!! 4) Function which mu_OPE used as the evaluation scale for OPE for RAD
!!!!!    arg=(mu,b,c1) with mu=scale(real_dp), b=transverse distance(real_dp), c1=scale variation parameter(real_dp)
real(dp),public:: muOPE
!!!!! 5) Function which returns b*-value used as argument b, in OPE part of RAD
!!!!!    arg=(b) with b=transverse distance(real_dp)
real(dp),public:: bSTAR
!!!!!------------------------------------------------------------------------------------

real(dp),allocatable::NPparam(:)
  
contains
  
!!!!!! Write nessecery model intitialization.
!!!!!! InitialNPParams is the initial NP-array (in principle, any non-pathological NP-array)
subroutine ModelInitialization(NPlength)
    integer,intent(in):: NPlength

    if(NPlength<4) then
        write(*,*) color('ART23-model: Number NP parameters for TMDR is less then 4',c_red)
        write(*,*) 'Evaluation STOP'
        stop
    end if
    
    allocate(NPparam(1:NPlength))
    NPparam=(/2d0,0.0001d0,0d0,1d0/)
    
    write(*,*) &
    color(">>>  The model for TMD evolution is ART23. Please, cite [1907.10356]&[2305.????]   <<<",c_cyan)

end subroutine ModelInitialization 

!!!!!! Write nessecery model update (e.g. save current NP-parameters)
!!!!!! newNPParams is the new NP-array
subroutine ModelUpdate(newNPParams)  
    real(dp),intent(in):: newNPParams(:)
    
    NPparam=newNPParams !! save new vector of NP-parameters

end subroutine ModelUpdate
  
 
!!! This is the rapidity anomalous dimension non-perturbative model
!!! In your evaluation take care that the saddle point is inside the pertrubative regeme
!!! Use function Dpert(mu,b,f) for D perturbative, use Dresum for D resum
function DNP(b,f)
    real(dp),intent(in)::b
    integer,intent(in)::f
    real(dp)::bS,bSS
    bS=bSTAR(b)

    bSS=b/(1_dp+b**6/20**6)**(1.d0/6.d0)

    DNP=NPparam(2)*bSS*bS+NPparam(3)*bSS*bS*Log(bS/NPparam(1))
end function DNP

!!! the function which is used insted of in the expression for perturbative RAD
!!!! This is the function b* that enters the logarithms of coefficient function
!!!! at small-b it should be ~b to match the collinear regime
!!!! at large-b it is a part of model
pure function bSTAR(b)
    real(dp),intent(in)::b

    bSTAR=b/SQRT(1_dp+b**2/NPparam(1)**2)
end function bSTAR

!!!!This function is the mu(b), which is used inside the OPE for RAD
!!!! c4-- is the scale variation variable
pure function muOPE(bT,c1)
    real(dp),intent(in)::bT,c1

    muOPE=C0_const/bSTAR(bT)*c1

    if(muOPE>1000d0) then
        muOPE=1000d0
    end if
end function muOPE
  
end module TMDR_model
