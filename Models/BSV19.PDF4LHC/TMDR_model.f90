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
    
    write(*,*) color(">>>  The model for TMD evolution is BSV19.PDF4LHC. Please, cite [1902.08474]   <<<",c_cyan)

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
  1.9266,  0.03655,& !!! mean
  2.3256,  0.0281,&  !!! best
  1.9266,  0.03655,&!!! mean
  2.2206,    0.0257 ,&
  2.0664,    0.0306 ,&
  2.2733,    0.0265 ,&
  1.1022,    0.0782 ,&
  2.1690,    0.0250 ,&
  2.2398,    0.0261 ,&
  1.0238,    0.0987 ,&
  1.1276,    0.0815 ,&
  2.2261,    0.0250 ,&
  2.0831,    0.0199 ,&
  2.2727,    0.0269 ,&
  2.2106,    0.0248 ,&
  2.2261,    0.0250 ,&
  1.1405,    0.0692 ,&
  2.2981,    0.0259 ,&
  1.1077,    0.0809 ,&
  2.2935,    0.0269 ,&
  1.1959,    0.0643 ,&
  1.1312,    0.0765 ,&
  2.2754,    0.0231 ,&
  1.8391,    0.0301 ,&
  1.9944,    0.0273 ,&
  1.2503,    0.0592 ,&
  1.0501,    0.0947 ,&
  2.2828,    0.0266 ,&
  2.3170,    0.0272 ,&
  1.8890,    0.0162 ,&
  3.6206,    0.0100 ,&
  1.0676,    0.0904 ,&
  2.1306,    0.0240 ,&
  1.5473,    0.0386 ,&
  2.2340,    0.0272 ,&
  2.1479,    0.0264 ,&
  1.1096,    0.0847 ,&
  1.9455,    0.0313 ,&
  2.2363,    0.0281 ,&
  2.2252,    0.0237 ,&
  2.1159,    0.0220 ,&
  1.3685,    0.0575 ,&
  1.1773,    0.0740 ,&
  2.2285,    0.0196 ,&
  2.2358,    0.0253 ,&
  2.1282,    0.0252 ,&
  1.2738,    0.0519 ,&
  2.2602,    0.0269 ,&
  1.5561,    0.0360 ,&
  2.2225,    0.0256 ,&
  2.2293,    0.0265 ,&
  1.1156,    0.0785 ,&
  2.2386,    0.0231 ,&
  1.5787,    0.0442 ,&
  2.2268,    0.0274 ,&
  1.6076,    0.0358 ,&
  1.1897,    0.0758 ,&
  1.7967,    0.0313 ,&
  2.2281,    0.0245 ,&
  1.6965,    0.0380 ,&
  2.2820,    0.0213 ,&
  1.6235,    0.0311 ,&
  2.1892,    0.0263 ,&
  2.0927,    0.0271 ,&
  2.2006,    0.0257 ,&
  2.0288,    0.0280 ,&
  2.1704,    0.0216 ,&
  2.2334,    0.0233 ,&
  2.2721,    0.0252 ,&
  1.8104,    0.0382 ,&
  1.2230,    0.0635 ,&
  2.1985,    0.0226 ,&
  2.2909,    0.0255 ,&
  1.9920,    0.0290 ,&
  2.2477,    0.0259 ,&
  2.1319,    0.0282 ,&
  2.1250,    0.0238 ,&
  2.2488,    0.0263 ,&
  1.3524,    0.0494 ,&
  1.9944,    0.0273 ,&
  1.0939,    0.0841 ,&
  2.2596,    0.0267 ,&
  2.2747,    0.0253 ,&
  2.1433,    0.0251 ,&
  2.1859,    0.0250 ,&
  2.2642,    0.0262 ,&
  2.1320,    0.0258 ,&
  2.2978,    0.0261 ,&
  1.3362,    0.0551 ,&
  2.2847,    0.0264 ,&
  1.6045,    0.0281 ,&
  2.0768,    0.0204 ,&
  2.2102,    0.0259 ,&
  2.3621,    0.0223 ,&
  2.0880,    0.0234 ,&
  1.6781,    0.0384 ,&
  1.7566,    0.0320 ,&
  2.3012,    0.0287 ,&
  1.0711,    0.0871 ,&
  2.2432,    0.0273 ,&
  2.2555,    0.0262 ,&
  2.3107,    0.0277 ,&
  1.6754,    0.0167 /)
 allocate(NParray(1:2))
 
 NParray=1d0*replicas((rep+2)*2+1:(rep+2)*2+2)
 
 end subroutine GetReplicaParameters
 
end module TMDR_model
