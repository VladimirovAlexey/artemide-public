!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Model for unpolarized TMD evolution for BSV19.bFIX.noLHC [1902.08474]
!
!			corresponds to bb* model, at fixed bMAX=2.5 GeV^{-1}
!			DNP=Dpert(b*)+g bb*
!			zeta=zetaPert(b*)
!
!			Requres two NP parameters (initated by best values values)
!
!				A.Vladimirov (20.02.2019)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

 module TMDR_model
use aTMDe_Numerics
use IO_functions
use TMD_AD, only : Dresum,zetaMUpert,zetaSL
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
!!!!! 4) Function which returns special equi-potential line, which is used in the evolution solutions
!!!!!    arg=(mu,b,f) with mu=scale(real_dp), b=transverse distance(real_dp), f=flavor(integer)
real(dp),public:: zetaNP
!!!!! 5) Subroutine which returns the array of NP-parameters corresponding to certain integer (replica)
!!!!!    arg=rep input integer,  NParray (real_dp(:), allocatable, intent(out))  returned array
public:: GetReplicaParameters
!!!!!------------------------------------------------------------------------------------

real(dp),allocatable::NPparam(:)
  
contains
  
!!!!!! Write nessecery model intitialization.
!!!!!! InitialNPParams is the initial NP-array (in principle, any non-pathological NP-array)
subroutine ModelInitialization(InitialNPParams)  
    real(dp),intent(in):: InitialNPParams(:)

    if(size(InitialNPParams)<2) then
        write(*,*) color('BSV19-model: Number NP parameters for TMDR is less then 2',c_red)
        write(*,*) 'Evaluation STOP'
        stop
    end if
    
    allocate(NPparam(1:size(InitialNPParams)))
    NPparam=InitialNPParams
    
    write(*,*) color(">>>  The model for TMD evolution is BSV19.bFIX.noLHC. Please, cite [1902.08474]   <<<",c_cyan)

end subroutine ModelInitialization 

!!!!!! Write nessecery model update (e.g. save current NP-parameters)
!!!!!! newNPParams is the new NP-array
subroutine ModelUpdate(newNPParams)  
    real(dp),intent(in):: newNPParams(:)
    
    NPparam=newNPParams !! save new vector of NP-parameters

end subroutine ModelUpdate
  
 
!!! This is the rapidity anomalous dimension non-pertrubative model
!!! In your evaluation take care that the saddle point is inside the pertrubative regeme
!!! Use function Dpert(mu,b,f) for D pertrubative, use Dresum for D resum
function DNP(mu,b,f)
    real(dp),intent(in)::mu,b
    integer,intent(in)::f
    real(dp)::bSTAR

    bSTAR=b/SQRT(1_dp+b**2/NPparam(1)**2)
    DNP=Dresum(mu,bSTAR,1)+NPparam(2)*bSTAR*b  !!!! D*+gK b b*, it smoother turns perturbative to b^2 assimptotic
    
end function DNP
  
!! This is the non-pertrubative shape of zeta_mu line.
!! It MUST follow the equipotential line in pertrubative regime (at small-b), at the level pf PT accuracy.
!! Otherwice, your evolution is completely broken.
!! Use zetaMUpert for perturbative values, use zetaSL for exact values
function zetaNP(mu,b,f)
    real(dp),intent(in)::mu,b
    integer,intent(in)::f
    real(dp)::sqrtBMAX
    
    sqrtBMAX=SQRT(1+b**2/NPparam(1)**2)
    zetaNP=zetaMUpert(mu,b/sqrtBMAX,f)
end function zetaNP
 
 !!! this is the table of replica prameters extracted in fit BSV19.
 !!! -2 is suggested for initialization replica
 !!! -1 is the best fit
 !!! 0 is the mean reaplics
 !!! 1 -- 100 replicas
subroutine GetReplicaParameters(rep,NParray)
    integer,intent(in)::rep
    real(dp),allocatable,intent(out)::NParray(:)
 real,parameter,dimension(1:206)::replicas=(/&
  2.5000,    0.011883,& !! suggested
  2.5000,    0.011883,&
  2.5000,    0.011883,&!! mean replica
  2.5000,    0.0109,&
  2.5000,    0.0100,&
  2.5000,    0.0127,&
  2.5000,    0.0100,&
  2.5000,    0.0100,&
  2.5000,    0.0107,&
  2.5000,    0.0301,&
  2.5000,    0.0103,&
  2.5000,    0.0100,&
  2.5000,    0.0120,&
  2.5000,    0.0105,&
  2.5000,    0.0101,&
  2.5000,    0.0123,&
  2.5000,    0.0100,&
  2.5000,    0.0100,&
  2.5000,    0.0101,&
  2.5000,    0.0110,&
  2.5000,    0.0218,&
  2.5000,    0.0102,&
  2.5000,    0.0102,&
  2.5000,    0.0102,&
  2.5000,    0.0143,&
  2.5000,    0.0100,&
  2.5000,    0.0108,&
  2.5000,    0.0114,&
  2.5000,    0.0101,&
  2.5000,    0.0105,&
  2.5000,    0.0106,&
  2.5000,    0.0100,&
  2.5000,    0.0135,&
  2.5000,    0.0102,&
  2.5000,    0.0100,&
  2.5000,    0.0100,&
  2.5000,    0.0100,&
  2.5000,    0.0107,&
  2.5000,    0.0103,&
  2.5000,    0.0103,&
  2.5000,    0.0109,&
  2.5000,    0.0100,&
  2.5000,    0.0323,&
  2.5000,    0.0100,&
  2.5000,    0.0101,&
  2.5000,    0.0117,&
  2.5000,    0.0115,&
  2.5000,    0.0100,&
  2.5000,    0.0100,&
  2.5000,    0.0100,&
  2.5000,    0.0102,&
  2.5000,    0.0100,&
  2.5000,    0.0101,&
  2.5000,    0.0115,&
  2.5000,    0.0102,&
  2.5000,    0.0100,&
  2.5000,    0.0102,&
  2.5000,    0.0100,&
  2.5000,    0.0103,&
  2.5000,    0.0178,&
  2.5000,    0.0117,&
  2.5000,    0.0102,&
  2.5000,    0.0101,&
  2.5000,    0.0101,&
  2.5000,    0.0103,&
  2.5000,    0.0107,&
  2.5000,    0.0101,&
  2.5000,    0.0201,&
  2.5000,    0.0149,&
  2.5000,    0.0101,&
  2.5000,    0.0197,&
  2.5000,    0.0113,&
  2.5000,    0.0128,&
  2.5000,    0.0100,&
  2.5000,    0.0310,&
  2.5000,    0.0100,&
  2.5000,    0.0112,&
  2.5000,    0.0101,&
  2.5000,    0.0100,&
  2.5000,    0.0110,&
  2.5000,    0.0101,&
  2.5000,    0.0100,&
  2.5000,    0.0177,&
  2.5000,    0.0100,&
  2.5000,    0.0103,&
  2.5000,    0.0109,&
  2.5000,    0.0102,&
  2.5000,    0.0103,&
  2.5000,    0.0178,&
  2.5000,    0.0174,&
  2.5000,    0.0102,&
  2.5000,    0.0170,&
  2.5000,    0.0101,&
  2.5000,    0.0130,&
  2.5000,    0.0107,&
  2.5000,    0.0118,&
  2.5000,    0.0100,&
  2.5000,    0.0101,&
  2.5000,    0.0100,&
  2.5000,    0.0100,&
  2.5000,    0.0100,&
  2.5000,    0.0127,&
  2.5000,    0.0100/)
  
  allocate(NParray(1:2))
  
  if(rep>100) then
   write(*,*) color('ERROR in BSV19.bFIX.noLHC model. It has only 100 replicas. Central replica is set',c_red)
   NParray=1d0*replicas((0+2)*2+1:(0+2)*2+2)
  else 
   NParray=1d0*replicas((rep+2)*2+1:(rep+2)*2+2)
  end if
 
 end subroutine GetReplicaParameters

 end module TMDR_model