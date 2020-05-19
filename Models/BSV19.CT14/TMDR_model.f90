!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Model for unpolarized TMD evolution for BSV19
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
    
    write(*,*) color(">>>  The model for TMD evolution is BSV19.CT14. Please, cite [1902.08474]   <<<",c_cyan)

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
    real(dp)::zz,rad
    
    rad=DNP(mu,b,f)
    zz=Exp(-b**2/NPparam(1)**2)
    zetaNP=zetaMUpert(mu,b,f)*zz+zetaSL(mu,rad,f)*(1d0-zz)
end function zetaNP
 
 !!! this is the table of replica prameters extracted in fit BSV19.
 !!! -2 is suggested for initialization replica
 !!! -1 is the best fit
 !!! 0 is the mean reaplics
 !!! 1 -- 100 replicas
subroutine GetReplicaParameters(rep,NParray)
    integer,intent(in)::rep
    real(dp),allocatable,intent(out)::NParray(:)
    
 real,parameter,dimension(1:206)::replicas=(/ &
  2.34665, 0.022735,& !!! mean
  2.5448, 0.0212,&  !!! best
  2.34665, 0.022735,&!!! mean
  2.4906,   0.0194,&
  2.4452,   0.0172,&
  2.4436,   0.0245,&
  2.5633,   0.0139,&
  2.4189,   0.0203,&
  2.4837,   0.0183,&
  1.3065,   0.0485,&
  2.5064,   0.0197,&
  2.2734,   0.0196,&
  2.3795,   0.0184,&
  1.3396,   0.0397,&
  2.4808,   0.0186,&
  1.1596,   0.0644,&
  1.7797,   0.0228,&
  2.7073,   0.0149,&
  1.3184,   0.0473,&
  2.5226,   0.0206,&
  2.2461,   0.0224,&
  2.6879,   0.0155,&
  1.3016,   0.0527,&
  2.5288,   0.0200,&
  2.4319,   0.0189,&
  1.1998,   0.0575,&
  2.5064,   0.0197,&
  2.5553,   0.0216,&
  2.5149,   0.0202,&
  2.2688,   0.0218,&
  2.5023,   0.0195,&
  2.4230,   0.0181,&
  2.2768,   0.0190,&
  1.8900,   0.0223,&
  2.4696,   0.0180,&
  1.5067,   0.0335,&
  2.5298,   0.0184,&
  2.3921,   0.0148,&
  2.4659,   0.0172,&
  2.5118,   0.0186,&
  1.7936,   0.0266,&
  2.5221,   0.0162,&
  2.4951,   0.0204,&
  3.1675,   0.0090,&
  2.4008,   0.0170,&
  2.4448,   0.0185,&
  2.1945,   0.0209,&
  2.4875,   0.0194,&
  1.8707,   0.0253,&
  2.4415,   0.0205,&
  2.2610,   0.0167,&
  2.5263,   0.0207,&
  2.3945,   0.0136,&
  2.5122,   0.0174,&
  2.5423,   0.0198,&
  2.4939,   0.0152,&
  2.4804,   0.0181,&
  2.5885,   0.0162,&
  2.4555,   0.0171,&
  2.3539,   0.0223,&
  1.8714,   0.0163,&
  2.4216,   0.0166,&
  2.5020,   0.0198,&
  1.1408,   0.0654,&
  3.5519,   0.0095,&
  2.6619,   0.0165,&
  1.4011,   0.0369,&
  2.1600,   0.0167,&
  2.1706,   0.0253,&
  2.4922,   0.0194,&
  2.4704,   0.0206,&
  2.7530,   0.0169,&
  2.6040,   0.0159,&
  2.9708,   0.0168,&
  2.2526,   0.0160,&
  2.1865,   0.0232,&
  2.4791,   0.0189,&
  2.4619,   0.0180,&
  2.5163,   0.0202,&
  1.1135,   0.0705,&
  2.3721,   0.0206,&
  1.1212,   0.0699,&
  2.5304,   0.0214,&
  5.4120,   0.0006,&
  1.1379,   0.0652,&
  2.2109,   0.0214,&
  2.4267,   0.0196,&
  1.8546,   0.0233,&
  2.4978,   0.0180,&
  2.4242,   0.0170,&
  2.1816,   0.0156,&
  2.4144,   0.0185,&
  2.3769,   0.0187,&
  2.4529,   0.0191,&
  2.9448,   0.0140,&
  2.5067,   0.0184,&
  2.4534,   0.0189,&
  4.3912,   0.0039,&
  1.3865,   0.0437,&
  2.2920,   0.0166,&
  2.5420,   0.0194,&
  2.5319,   0.0186,&
  3.7664,   0.0030/)
 
 allocate(NParray(1:2))
 
 NParray=1d0*replicas((rep+2)*2+1:(rep+2)*2+2)
 
 end subroutine GetReplicaParameters

end module TMDR_model