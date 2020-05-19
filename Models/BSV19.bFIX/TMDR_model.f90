!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Model for unpolarized TMD evolution for BSV19.bFIX [1902.08474]
!
!			corresponds to bb* model, at b =2.5 GeV^{-1}
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
    
    write(*,*) color(">>>  The model for TMD evolution is BSV19.bFIX. Please, cite [1902.08474]   <<<",c_cyan)

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
  2.5000,    0.037484,& !! suggested
  2.5000,    0.037484,&
  2.5000,    0.037484,&!! mean replica
  2.5000,    0.0263,&
  2.5000,    0.0414,&
  2.5000,    0.0310,&
  2.5000,    0.0307,&
  2.5000,    0.0338,&
  2.5000,    0.0431,&
  2.5000,    0.0398,&
  2.5000,    0.0368,&
  2.5000,    0.0310,&
  2.5000,    0.0273,&
  2.5000,    0.0248,&
  2.5000,    0.0402,&
  2.5000,    0.0410,&
  2.5000,    0.0310,&
  2.5000,    0.0359,&
  2.5000,    0.0348,&
  2.5000,    0.0332,&
  2.5000,    0.0300,&
  2.5000,    0.0248,&
  2.5000,    0.0389,&
  2.5000,    0.0355,&
  2.5000,    0.0411,&
  2.5000,    0.0389,&
  2.5000,    0.0341,&
  2.5000,    0.0316,&
  2.5000,    0.0357,&
  2.5000,    0.0437,&
  2.5000,    0.0464,&
  2.5000,    0.0436,&
  2.5000,    0.0427,&
  2.5000,    0.0420,&
  2.5000,    0.0102,&
  2.5000,    0.0311,&
  2.5000,    0.0393,&
  2.5000,    0.0479,&
  2.5000,    0.0281,&
  2.5000,    0.0388,&
  2.5000,    0.0446,&
  2.5000,    0.0361,&
  2.5000,    0.0360,&
  2.5000,    0.0268,&
  2.5000,    0.0387,&
  2.5000,    0.0315,&
  2.5000,    0.0436,&
  2.5000,    0.0396,&
  2.5000,    0.0468,&
  2.5000,    0.0433,&
  2.5000,    0.0369,&
  2.5000,    0.0434,&
  2.5000,    0.0332,&
  2.5000,    0.0475,&
  2.5000,    0.0418,&
  2.5000,    0.0489,&
  2.5000,    0.0280,&
  2.5000,    0.0324,&
  2.5000,    0.0293,&
  2.5000,    0.0442,&
  2.5000,    0.0389,&
  2.5000,    0.0382,&
  2.5000,    0.0492,&
  2.5000,    0.0214,&
  2.5000,    0.0342,&
  2.5000,    0.0506,&
  2.5000,    0.0380,&
  2.5000,    0.0343,&
  2.5000,    0.0383,&
  2.5000,    0.0370,&
  2.5000,    0.0435,&
  2.5000,    0.0463,&
  2.5000,    0.0423,&
  2.5000,    0.0463,&
  2.5000,    0.0375,&
  2.5000,    0.0372,&
  2.5000,    0.0339,&
  2.5000,    0.0421,&
  2.5000,    0.0419,&
  2.5000,    0.0342,&
  2.5000,    0.0333,&
  2.5000,    0.0307,&
  2.5000,    0.0395,&
  2.5000,    0.0382,&
  2.5000,    0.0428,&
  2.5000,    0.0385,&
  2.5000,    0.0454,&
  2.5000,    0.0395,&
  2.5000,    0.0305,&
  2.5000,    0.0456,&
  2.5000,    0.0420,&
  2.5000,    0.0396,&
  2.5000,    0.0423,&
  2.5000,    0.0335,&
  2.5000,    0.0329,&
  2.5000,    0.0389,&
  2.5000,    0.0455,&
  2.5000,    0.0389,&
  2.5000,    0.0369,&
  2.5000,    0.0368,&
  2.5000,    0.0323,&
  2.5000,    0.0511,&
  2.5000,    0.0293/)
  
  allocate(NParray(1:2))
  
  if(rep>100) then
   write(*,*) color('ERROR in BSV19.bFIX model. It has only 100 replicas. Central replica is set',c_red)
   NParray=1d0*replicas((0+2)*2+1:(0+2)*2+2)
  else 
   NParray=1d0*replicas((rep+2)*2+1:(rep+2)*2+2)
  end if
 
 end subroutine GetReplicaParameters
 
end module TMDR_model
