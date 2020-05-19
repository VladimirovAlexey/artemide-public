!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Model for unpolarized TMD PDF  BSV19  [1902.08474]
!
!				A.Vladimirov (11.07.2019)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
module uTMDPDF_model
use aTMDe_Numerics
use IO_functions
implicit none

private

!!!!!------------------------------------------------------------------------------------
!!!!! These functions MUST defined in module  !!
!!!!!
!!!!! 1) The subroutine is called during the initialization of TMD-module
!!!!!    arg=array of initial NP-parameters
public:: ModelInitialization
!!!!! 2) The subroutine that is called on reset of NP-parameters in TMD-module
!!!!!    arg=array of new NP-parameters
public:: ModelUpdate
!!!!! 3) Function which returns FNP function
!!!!!    arg=(x,z,b,hadron,lambdaNP) with x=x_Bj for TMD (real_dp), z=convolution variable(real_dp), 
!!!!!    b=transverse distance(real_dp), hadron=number of the hadron in grid(integer)
!!!!!    lambdaNP = array of NP parameters (real_dp(:))
real(dp),public,dimension(-5:5):: FNP
!!!!! 4) Function which returns the value of b used as argument of convolution integrals
!!!!!    arg=(b,lambdaNP) with b=transverse distance(real_dp), lambdaNP = array of NP parameters (real_dp(:))
real(dp),public:: bSTAR
!!!!! 5) Function which returns the scale of matching (OPE scale)
!!!!!    arg=(z,bt) with z=convolution variable(real_dp), b=transverse distance(real_dp)
real(dp),public:: mu_OPE
!!!!! 6) Subroutine which returns the array of parameters CA which compose the TMDs into a single one
!!!!!    i.e. the TMD for hardon=h is build as TMD(h)=Sum_c CA(h,c) TMD(c)
!!!!!    it is used only if the option UseComposite TMD is ON,
!!!!!    arg=(h,lambdaNP,includeArray,CA) with h=hadron(integer),lambdaNP = array of NP parameters (real_dp(:))
!!!!!    includeArray=logical array with .true. for terms included in the sum (logical(:),allocatable,intent(out))
!!!!!    CA=coefficient CA (real_dp(:),allocatable,intent(out))
public:: GetCompositionArray
!!!!! 7) Subroutine which returns the array of NP-parameters corresponding to certain integer (replica)
!!!!!    arg=rep input integer,  NParray (real_dp(:), allocatable, intent(out))  returned array
public:: GetReplicaParameters
!!!!!------------------------------------------------------------------------------------

real(dp),allocatable::NPparam(:)

contains  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! USER DEFINED FUNCTIONS   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!! Write nessecery model intitialization.
subroutine ModelInitialization(NPstart)
    real(dp),intent(in)::NPstart(:)
    allocate(NPparam(1:size(NPstart)))
    NPparam=NPstart
    
    write(*,*) color(">>>  The model for uTMDPDF is BSV19.CT14. Please, cite [1902.08474]   <<<",c_cyan)
    
end subroutine ModelInitialization

!!!!!! Write nessecery model update (e.g. save current NP-parameters)
!!!!!! newNPParams is the new NP-array
subroutine ModelUpdate(newNPParams)  
    real(dp),intent(in):: newNPParams(:)
    
    NPparam=newNPParams !! save new vector of NP-parameters

end subroutine ModelUpdate
  
!!! This is  non-pertrubative function
!!! non=pertrubative parameters are lambdaNP()
!!! x-- is the bjorken variable of TMD
!!! z-- is convolution variable
function FNP(x,z,bT,hadron,lambdaNP)
    real(dp),intent(in)::x,z,bT    
    integer,intent(in)::hadron
    real(dp),intent(in)::lambdaNP(:)

    real(dp)::bb,w1,w2,w3,FNP0

    bb=bT**2
    w1=lambdaNP(1)*(1-x)+x*lambdaNP(2)+x*(1-x)*lambdaNP(5)
    w2=lambdaNP(3)*x**lambdaNP(4)+lambdaNP(6)

    if(w2<0d0 .or. w1<0d0) then !!! case of negative power, we return absolutely incorrect expression.
        if(bT<1d0) then
    FNP0=-1d0
        else
    FNP0=0d0
        end if
    else
    FNP0=Exp(-w1*bb/sqrt(1+w2*bb))
    end if

    !    FNP0=Exp(-lambdaNP(1)*bb)

    FNP=FNP0*(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)

end function FNP
  
  !!!! This is the function b* that enter the logarithms of coefficient function
  !!!! at small-b it should be ~b to match the collinear regime
  !!!! at large-b it is a part of model
  !!!! NOTE: if it is lambda-dependent, the grid will be recalculate each reset of lambdaNP
pure function bSTAR(bT,lambdaNP)
    real(dp),intent(in)::bT
    real(dp),intent(in)::lambdaNP(:)

    bSTAR=bT/sqrt(1d0+(bT/500d0)**2)

end function bSTAR
  
    !!!!This function is the mu(x,b), which is used inside the OPE
pure function mu_OPE(z,bt)
    real(dp),intent(in)::z,bt

    mu_OPE=C0_const*1d0/bT+2d0

    if(mu_OPE>1000d0) then
        mu_OPE=1000d0
    end if
end function mu_OPE
  
!!!! if the option UseComposite TMD is OFF, this function is ignored
!!!! If the option UseComposite TMD is ON,
!!!! than the TMD for hardon is build as TMD(hadron)=Sum_c CA(h,c) TMD(c)
!!!! where h=hadron, CA=coefficientArray
!!!! coefficientArray real(dp) list of coefficeints
!!!! includeArray is logical array list (true=TMD(c) is computed, false TMD(c) ignored)
subroutine GetCompositionArray(hadron,lambdaNP,includeArray,coefficientArray)  
    real(dp),intent(in)::lambdaNP(:)
    integer::hadron
    logical,allocatable,intent(out)::includeArray(:)
    real(dp),allocatable,intent(out)::coefficientArray(:)

    allocate(includeArray(1:1))
    allocate(coefficientArray(1:1))
end subroutine GetCompositionArray
  
  
!!! In SV19 model the replica parameters are stored in separate file.
subroutine GetReplicaParameters(rep,NParray)
    integer,intent(in)::rep
    real(dp),allocatable,intent(out)::NParray(:)
 real,parameter,dimension(1:618)::replicas=(/&
   0.27706,  24.8883, 1241.34,      2.66869, -23.7589,   0.1, & !!!! initialisation
    0.2770,   23.6983, 1280.7698,    2.5709,  -22.7613,   0.0,&     !!! best
   0.27706,  24.8883, 1241.34,      2.66869, -23.7589,   0.0,&      !!!!! mean
   0.2621,   24.1856, 1271.3300,    2.6204,  -22.8528,    0.0000,&
   0.2662,   23.9763, 1257.0092,    2.6706,  -22.3719,    0.0000,&
   0.3158,   23.0570, 1159.6480,    2.6773,  -23.7376,    0.0000,&
   0.2730,   25.8214, 1282.8932,    2.6427,  -24.0532,    0.0000,&
   0.2839,   23.7717, 1296.6336,    2.6668,  -22.7662,    0.0000,&
   0.2727,   23.8674, 1295.3886,    2.6496,  -22.5049,    0.0000,&
   0.2606,   22.8846,  848.1077,    2.7590,  -22.0318,    0.0000,&
   0.2743,   23.7405, 1285.2846,    2.6120,  -22.5287,    0.0000,&
   0.2772,   26.1798, 1284.1215,    2.6388,  -25.0891,    0.0000,&
   0.2792,   24.4381, 1287.4303,    2.6491,  -22.6446,    0.0000,&
   0.2688,   39.9039, 2548.7189,    2.9660,  -39.5333,    0.0000,&
   0.2726,   24.2872, 1299.4255,    2.6015,  -22.5568,    0.0000,&
   0.2229,   26.1772, 1179.7598,    2.6953,  -24.9276,    0.0000,&
   0.2548,   24.8849, 1337.3501,    2.7625,  -22.2860,    0.0000,&
   0.3101,   23.6535, 1074.5805,    2.6476,  -23.1973,    0.0000,&
   0.2618,   30.5317, 1378.4227,    2.7619,  -29.8677,    0.0000,&
   0.2942,   23.8081, 1282.1327,    2.6476,  -23.3002,    0.0000,&
   0.2942,   24.6880, 1340.4406,    2.7643,  -24.5073,    0.0000,&
   0.2589,   25.4002, 1329.2479,    2.5498,  -23.4590,    0.0000,&
   0.2387,   27.7380, 1295.4606,    2.6731,  -26.5563,    0.0000,&
   0.2769,   23.7923, 1291.2774,    2.6041,  -22.5376,    0.0000,&
   0.2869,   24.5613, 1318.2916,    2.7191,  -23.7513,    0.0000,&
   0.2547,   25.6850, 1316.9782,    2.9080,  -25.0760,    0.0000,&
   0.2743,   23.7405, 1285.2846,    2.6120,  -22.5287,    0.0000,&
   0.2873,   23.8748, 1286.2306,    2.5776,  -23.2444,    0.0000,&
   0.2746,   23.7292, 1285.4617,    2.6023,  -22.5317,    0.0000,&
   0.2634,   24.0179, 1178.3897,    2.5505,  -22.5841,    0.0000,&
   0.2742,   23.8165, 1291.3831,    2.6193,  -22.5564,    0.0000,&
   0.2698,   23.8111, 1289.3909,    2.6568,  -22.1707,    0.0000,&
   0.2658,   24.0632, 1079.2344,    2.6147,  -22.4124,    0.0000,&
   0.2441,   26.8544, 1429.6968,    2.6153,  -23.5993,    0.0000,&
   0.2695,   23.6683, 1306.9278,    2.6528,  -21.9875,    0.0000,&
   0.2998,   23.6415,  841.2932,    2.7752,  -23.1397,    0.0000,&
   0.2811,   24.0827, 1297.2650,    2.6103,  -22.6773,    0.0000,&
   0.2683,   24.9195, 1099.6117,    2.5620,  -22.9143,    0.0000,&
   0.2704,   24.0754, 1298.6630,    2.6733,  -22.3713,    0.0000,&
   0.3001,   25.1511, 1336.1152,    2.7641,  -25.0119,    0.0000,&
   0.2406,   31.2105, 1397.4399,    2.5138,  -28.4379,    0.0000,&
   0.2788,   24.0007, 1259.0519,    2.6486,  -22.3039,    0.0000,&
   0.2740,   23.8214, 1289.3971,    2.6133,  -22.5154,    0.0000,&
   0.3224,   24.8511, 1407.3178,    2.9154,  -25.0727,    0.0000,&
   0.2560,   24.0035, 1093.2701,    2.5004,  -21.6406,    0.0000,&
   0.2719,   23.8917, 1297.6660,    2.6467,  -22.3307,    0.0000,&
   0.3206,   23.6523, 1088.8532,    2.8017,  -23.7783,    0.0000,&
   0.2720,   23.8624, 1297.3434,    2.6406,  -22.6064,    0.0000,&
   0.2970,   24.7999, 1364.0854,    2.9486,  -24.7155,    0.0000,&
   0.3033,   23.1662,  940.5174,    2.5302,  -22.7305,    0.0000,&
   0.3078,   23.9060,  851.7359,    2.6239,  -23.0871,    0.0000,&
   0.2438,   24.5287, 1537.8895,    2.6265,  -22.6140,    0.0000,&
   0.2887,   23.6687,  972.7202,    2.5910,  -22.4446,    0.0000,&
   0.2704,   24.2495, 1044.3275,    2.5500,  -22.9735,    0.0000,&
   0.2927,   24.1707, 1304.9610,    2.6976,  -23.7820,    0.0000,&
   0.2738,   23.6002, 1147.2858,    2.6907,  -22.3376,    0.0000,&
   0.2425,   26.0435, 1292.3148,    2.5219,  -23.9128,    0.0000,&
   0.2878,   23.8044, 1058.7130,    2.5940,  -22.6906,    0.0000,&
   0.2841,   24.2994, 1304.1583,    2.6837,  -23.0902,    0.0000,&
   0.2851,   23.9108, 1290.8387,    2.7022,  -23.5193,    0.0000,&
   0.2936,   42.9284, 3231.9038,    3.1184,  -42.4205,    0.0000,&
   0.2758,   23.6199, 1034.0246,    2.5989,  -22.3215,    0.0000,&
   0.2758,   23.7773, 1288.5363,    2.6245,  -22.5404,    0.0000,&
   0.1555,   27.9437, 1600.6662,    2.6649,  -24.4804,    0.0000,&
   0.3267,   23.6230,  878.3110,    2.5579,  -23.7963,    0.0000,&
   0.3124,   24.7588, 1212.4103,    2.7320,  -24.7765,    0.0000,&
   0.2698,   31.8448, 1336.8453,    2.6888,  -30.5732,    0.0000,&
   0.3007,   25.2846,  880.8124,    2.5523,  -24.0775,    0.0000,&
   0.2764,   23.6352, 1329.1384,    2.7585,  -23.2288,    0.0000,&
   0.2814,   23.8090, 1302.0337,    2.6408,  -22.5519,    0.0000,&
   0.2918,   24.4383, 1311.2853,    2.7092,  -23.9413,    0.0000,&
   0.2643,   24.2109, 1293.2440,    2.5932,  -22.6417,    0.0000,&
   0.2998,   24.0796, 1099.6067,    2.5982,  -23.2346,    0.0000,&
   0.3255,   20.3317,  716.9703,    2.5067,  -20.5550,    0.0000,&
   0.2935,   24.1213, 1189.2255,    2.7868,  -22.9229,    0.0000,&
   0.3034,   23.2724,  958.2393,    2.6028,  -22.8187,    0.0000,&
   0.2713,   23.9126, 1290.2106,    2.6134,  -22.3628,    0.0000,&
   0.2783,   24.2320, 1309.4283,    2.6833,  -22.9436,    0.0000,&
   0.2734,   23.7272, 1287.0896,    2.6181,  -22.8013,    0.0000,&
   0.2067,   22.3731, 1179.3521,    2.7022,  -20.5138,    0.0000,&
   0.3091,   23.9846, 1157.7125,    2.6617,  -23.6612,    0.0000,&
   0.2049,   23.9036, 1233.8305,    2.6777,  -21.7369,    0.0000,&
   0.3115,   22.3052, 1037.5805,    2.6550,  -22.6350,    0.0000,&
   0.3289,   24.5052, 1246.0654,    2.7235,  -24.2835,    0.0000,&
   0.1753,   27.9768, 1508.1282,    2.5805,  -24.5904,    0.0000,&
   0.2725,   25.8952, 1290.5683,    2.5898,  -24.4704,    0.0000,&
   0.3042,   27.1608, 1325.7357,    2.7098,  -27.3987,    0.0000,&
   0.3097,   25.0608,  969.2283,    2.7628,  -24.9842,    0.0000,&
   0.2803,   23.9887, 1296.3842,    2.6332,  -22.5729,    0.0000,&
   0.2456,   25.2210, 1428.8328,    2.5793,  -22.7584,    0.0000,&
   0.2828,   24.6206,  914.8741,    2.5382,  -22.8197,    0.0000,&
   0.2798,   23.9971, 1287.3429,    2.6256,  -22.5485,    0.0000,&
   0.2428,   24.1522, 1409.2876,    2.7226,  -22.0794,    0.0000,&
   0.2755,   23.8227,  993.4868,    2.5350,  -22.8339,    0.0000,&
   0.2779,   24.0332, 1218.6528,    2.6248,  -22.8505,    0.0000,&
   0.2775,   24.0270, 1298.0773,    2.6137,  -22.5152,    0.0000,&
   0.3146,   24.8984, 1148.5847,    2.6898,  -24.9846,    0.0000,&
   0.3148,   22.8581,  852.8205,    2.5149,  -22.0312,    0.0000,&
   0.2329,   26.1313, 1150.0908,    2.5849,  -23.8841,    0.0000,&
   0.3051,   25.1115, 1298.5293,    2.8407,  -24.6607,    0.0000,&
   0.3079,   24.5819, 1330.7214,    2.7656,  -24.8106,    0.0000,&
   0.2649,   23.8632, 1202.9764,    2.5840,  -22.3515,    0.0000,&
   0.2949,   24.8866,    0.7977,    3.4907,  -24.4749,    0.0000/)
    
    allocate(NParray(1:size(NPparam)))
    if(rep>100) then
        write(*,*) color('ERROR in BSM19.CT14 model. It has only 100 replicas. Central replica is set',c_red)
        NParray=1d0*replicas((0+2)*6+1:(0+2)*6+6)
    else
        NParray=1d0*replicas((rep+2)*6+1:(rep+2)*6+6)
    end if
end subroutine GetReplicaParameters

end module uTMDPDF_model
