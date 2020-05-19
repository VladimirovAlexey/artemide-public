!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Model for unpolarized TMD evolution for BSV19.bFIT [1902.08474]
!
!			corresponds to bb* model
!			DNP=Dpert(b*)+g bb*
!			zeta=zetaPert(b*)
!
!			Requres two NP parameters (initated by best values values)
!
!				A.Vladimirov (26.12.2018)
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
    
    write(*,*) color(">>>  The model for TMD evolution is BSV19.bFIT. Please, cite [1902.08474]   <<<",c_cyan)

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
 real,parameter,dimension(1:606)::replicas=(/&
  3.30213,   0.0237183,& !! suggested
  0.,0.,                &
  3.30213,   0.0237183,&!! mean replica
  3.3580,    0.0209,&
  3.4604,    0.0283,&
  3.3426,    0.0214,&
  3.9970,    0.0113,&
  3.3484,    0.0279,&
  3.5535,    0.0180,&
  3.2125,    0.0324,&
  3.3701,    0.0245,&
  3.5259,    0.0290,&
  3.3457,    0.0102,&
  3.3353,    0.0228,&
  3.7684,    0.0215,&
  3.3516,    0.0256,&
  3.4591,    0.0177,&
  3.3242,    0.0249,&
  3.3461,    0.0210,&
  2.7297,    0.0398,&
  3.2607,    0.0229,&
  3.4936,    0.0180,&
  3.2106,    0.0206,&
  3.4415,    0.0168,&
  3.3656,    0.0260,&
  3.3589,    0.0206,&
  3.2515,    0.0228,&
  3.1335,    0.0187,&
  3.2726,    0.0201,&
  3.3500,    0.0221,&
  3.3085,    0.0219,&
  3.3616,    0.0280,&
  3.1933,    0.0239,&
  3.2602,    0.0280,&
  3.2324,    0.0288,&
  3.1244,    0.0236,&
  3.3895,    0.0279,&
  3.0724,    0.0180,&
  3.2062,    0.0149,&
  2.6114,    0.0319,&
  3.3651,    0.0255,&
  3.2676,    0.0240,&
  3.3022,    0.0272,&
  3.1225,    0.0182,&
  2.2804,    0.0387,&
  3.1848,    0.0311,&
  3.4637,    0.0150,&
  3.1035,    0.0286,&
  3.2315,    0.0212,&
  3.3793,    0.0292,&
  3.3564,    0.0265,&
  3.3374,    0.0178,&
  2.6811,    0.0260,&
  3.3924,    0.0280,&
  3.3682,    0.0244,&
  3.4057,    0.0174,&
  3.3016,    0.0277,&
  3.1848,    0.0311,&
  3.9031,    0.0178,&
  3.6019,    0.0215,&
  3.3703,    0.0192,&
  3.0483,    0.0262,&
  3.7700,    0.0112,&
  2.7808,    0.0308,&
  3.3883,    0.0270,&
  3.5692,    0.0164,&
  3.3279,    0.0214,&
  3.3077,    0.0216,&
  3.3513,    0.0270,&
  2.9579,    0.0272,&
  3.1271,    0.0263,&
  3.3541,    0.0221,&
  3.3822,    0.0285,&
  3.3674,    0.0238,&
  3.3799,    0.0276,&
  3.0187,    0.0240,&
  3.2956,    0.0245,&
  3.3564,    0.0265,&
  3.5661,    0.0304,&
  2.9854,    0.0174,&
  3.2335,    0.0233,&
  3.5778,    0.0280,&
  3.5192,    0.0261,&
  2.4829,    0.0197,&
  3.3551,    0.0221,&
  3.3392,    0.0214,&
  3.3757,    0.0265,&
  3.1645,    0.0210,&
  3.3930,    0.0230,&
  3.3461,    0.0199,&
  3.4621,    0.0262,&
  3.3110,    0.0260,&
  3.9456,    0.0124,&
  3.3659,    0.0212,&
  3.1949,    0.0219,&
  3.2106,    0.0206,&
  3.4845,    0.0260,&
  3.3093,    0.0219,&
  3.0940,    0.0183,&
  3.3669,    0.0264,&
  3.1489,    0.0272,&
  3.5535,    0.0180,&
  3.3549,    0.0298,&
  3.3530,    0.0289,&
  3.2915,    0.0163,&
  3.3184,    0.0211,&
  3.3593,    0.0229,&
  3.1902,    0.0208,&
  3.6908,    0.0232,&
  3.4587,    0.0251,&
  3.1748,    0.0207,&
  3.3723,    0.0228,&
  3.1598,    0.0277,&
  3.4216,    0.0259,&
  3.3189,    0.0244,&
  3.5332,    0.0257,&
  3.5938,    0.0195,&
  3.2919,    0.0194,&
  3.4006,    0.0253,&
  3.1035,    0.0286,&
  3.3160,    0.0268,&
  3.3692,    0.0274,&
  3.3644,    0.0266,&
  3.3344,    0.0218,&
  3.0141,    0.0278,&
  3.1575,    0.0265,&
  3.3911,    0.0281,&
  3.4433,    0.0288,&
  3.4216,    0.0259,&
  3.3341,    0.0245,&
  3.7341,    0.0136,&
  3.1250,    0.0224,&
  3.0817,    0.0290,&
  2.8704,    0.0216,&
  3.3426,    0.0214,&
  3.3855,    0.0277,&
  3.2708,    0.0206,&
  3.3723,    0.0228,&
  3.7463,    0.0206,&
  3.3564,    0.0265,&
  3.3495,    0.0268,&
  3.3580,    0.0209,&
  3.3801,    0.0278,&
  3.3072,    0.0246,&
  3.3733,    0.0228,&
  3.8563,    0.0190,&
  3.3926,    0.0241,&
  3.4241,    0.0291,&
  3.4433,    0.0288,&
  2.7841,    0.0102,&
  3.3800,    0.0220,&
  3.5080,    0.0247,&
  3.3699,    0.0225,&
  3.3650,    0.0158,&
  3.3461,    0.0210,&
  3.2754,    0.0196,&
  3.0483,    0.0262,&
  3.3534,    0.0235,&
  3.0906,    0.0247,&
  3.3110,    0.0260,&
  3.3650,    0.0211,&
  3.2236,    0.0239,&
  3.3626,    0.0254,&
  1.9807,    0.0448,&
  2.2804,    0.0387,&
  2.9727,    0.0339,&
  3.3525,    0.0199,&
  3.8067,    0.0250,&
  3.0292,    0.0298,&
  3.3892,    0.0245,&
  3.3084,    0.0254,&
  3.1455,    0.0188,&
  2.4034,    0.0380,&
  2.8004,    0.0283,&
  3.3677,    0.0238,&
  3.1933,    0.0239,&
  3.5535,    0.0180,&
  3.2824,    0.0230,&
  2.8350,    0.0300,&
  3.0316,    0.0275,&
  2.7519,    0.0126,&
  3.1760,    0.0234,&
  3.3398,    0.0187,&
  3.3651,    0.0239,&
  3.2053,    0.0227,&
  3.9864,    0.0176,&
  3.3733,    0.0228,&
  2.5463,    0.0314,&
  3.4092,    0.0278,&
  3.2919,    0.0194,&
  3.2761,    0.0101,&
  2.2119,    0.0507,&
  3.4884,    0.0248,&
  2.7297,    0.0398,&
  3.4115,    0.0223,&
  3.3645,    0.0236,&
  3.4202,    0.0231,&
  3.2159,    0.0245,&
  3.9308,    0.0124,&
  3.3329,    0.0224,&
  2.8671,    0.0262,&
  3.3384,    0.0121,&
  3.1760,    0.0234,&
  3.3114,    0.0238,&
  2.6621,    0.0233,&
  3.3561,    0.0272,&
  3.5157,    0.0281,&
  3.3680,    0.0262,&
  3.1434,    0.0216,&
  3.4644,    0.0116,&
  3.1434,    0.0268,&
  3.2567,    0.0243,&
  2.8350,    0.0300,&
  3.4972,    0.0175,&
  3.0889,    0.0221,&
  3.1875,    0.0245,&
  3.3876,    0.0158,&
  3.9031,    0.0127,&
  3.9951,    0.0158,&
  3.1271,    0.0263,&
  3.3719,    0.0242,&
  3.3812,    0.0277,&
  3.3674,    0.0218,&
  3.9515,    0.0170,&
  3.3430,    0.0201,&
  3.0508,    0.0265,&
  3.1895,    0.0227,&
  3.3677,    0.0301,&
  3.3072,    0.0246,&
  3.2889,    0.0274,&
  3.3170,    0.0222,&
  3.3524,    0.0188,&
  3.2777,    0.0230,&
  3.3424,    0.0238,&
  3.5606,    0.0209,&
  3.4920,    0.0261,&
  3.3651,    0.0239,&
  3.6434,    0.0167,&
  3.3612,    0.0205,&
  3.4762,    0.0296,&
  3.1489,    0.0272,&
  3.9749,    0.0159,&
  3.7286,    0.0191,&
  2.8695,    0.0103,&
  3.9720,    0.0158,&
  3.3426,    0.0248,&
  3.3476,    0.0211,&
  3.2103,    0.0226,&
  3.3384,    0.0224,&
  3.3571,    0.0160,&
  3.4241,    0.0291,&
  2.9818,    0.0260,&
  3.3834,    0.0248,&
  3.3674,    0.0238,&
  3.6240,    0.0260,&
  3.3461,    0.0210,&
  3.3927,    0.0268,&
  2.9630,    0.0188,&
  3.3639,    0.0266,&
  3.5068,    0.0268,&
  3.3233,    0.0228,&
  2.8551,    0.0180,&
  3.5062,    0.0273,&
  3.2266,    0.0206,&
  3.3189,    0.0244,&
  3.7433,    0.0317,&
  3.5152,    0.0144,&
  3.0306,    0.0272,&
  3.3603,    0.0214,&
  3.7649,    0.0219,&
  3.2982,    0.0302,&
  3.0940,    0.0183,&
  3.3753,    0.0276,&
  3.3760,    0.0255,&
  3.5606,    0.0306,&
  3.3760,    0.0271,&
  3.3037,    0.0286,&
  3.6019,    0.0215,&
  3.4774,    0.0247,&
  3.2789,    0.0221,&
  3.3313,    0.0244,&
  3.4996,    0.0280,&
  3.3116,    0.0224,&
  2.8970,    0.0225,&
  3.2463,    0.0185,&
  2.6036,    0.0280,&
  3.2463,    0.0185,&
  3.4366,    0.0283,&
  3.9617,    0.0116,&
  2.9146,    0.0395,&
  3.3144,    0.0235,&
  3.3549,    0.0298,&
  3.3651,    0.0255,&
  3.3757,    0.0265,&
  3.3869,    0.0279,&
  3.2341,    0.0169,&
  3.2179,    0.0217,&
  3.9716,    0.0160,&
  3.2790,    0.0262,&
  3.3624,    0.0228,&
  3.1727,    0.0233,&
  2.6114,    0.0319,&
  3.1401,    0.0220/)
  
  allocate(NParray(1:2))
  
  if(rep>300) then
   write(*,*) color('ERROR in BSV19.bFIT model. It has only 300 replicas. Central replica is set',c_red)
   NParray=1d0*replicas((0+2)*2+1:(0+2)*2+2)
  else 
   NParray=1d0*replicas((rep+2)*2+1:(rep+2)*2+2)
  end if
 
 end subroutine GetReplicaParameters
 
end module TMDR_model
