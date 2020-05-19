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
    
    write(*,*) color(">>>  The model for uTMDPDF is BSV19.bFIX. Please, cite [1902.08474]   <<<",c_cyan)
    
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
    0.248353, 8.15443,   275.912,   2.52581,  -4.96056,    0.1000,&!!!the suggested fo initialization replica is slightly wider. It prevent the 0-values at large b
    0.248353, 8.15443,   275.912,   2.52581,  -4.96056,    0.0000,&!!! best replica
    0.248353, 8.15443,   275.912,   2.52581,  -4.96056,    0.0000,&!!! mean replica
    0.3183,    9.7583,  349.1763,    2.9585,   -8.1546,    0.0000,&
    0.2157,    7.6595,  288.1890,    2.4484,   -3.9180,    0.0000,&
    0.2980,   11.5897,  374.8948,    2.7785,   -9.7057,    0.0000,&
    0.2717,    7.7718,  272.2432,    2.5266,   -4.2645,    0.0000,&
    0.2485,    8.7815,  279.5729,    2.4171,   -4.9592,    0.0000,&
    0.2360,    7.3568,  250.6606,    2.5625,   -4.4443,    0.0000,&
    0.2378,    7.7165,  306.6578,    2.5476,   -4.2307,    0.0000,&
    0.2502,    7.8545,  243.4586,    2.4840,   -4.8395,    0.0000,&
    0.2732,    7.9270,  250.7517,    2.4977,   -4.6666,    0.0000,&
    0.2414,    8.3407,  292.6242,    2.4010,   -3.8862,    0.0000,&
    0.3059,   12.6697,  420.0725,    2.7933,  -10.5186,    0.0000,&
    0.2317,    7.4379,  244.1309,    2.4200,   -4.0075,    0.0000,&
    0.2427,    7.6630,  250.8475,    2.4497,   -4.4898,    0.0000,&
    0.2980,   11.5897,  374.8948,    2.7785,   -9.7057,    0.0000,&
    0.2512,    7.7993,  235.1599,    2.4716,   -4.7726,    0.0000,&
    0.2467,    7.5285,  274.6018,    2.5640,   -4.2796,    0.0000,&
    0.2559,    8.2205,  259.1454,    2.4725,   -4.8577,    0.0000,&
    0.2544,    8.1088,  229.5658,    2.3621,   -4.4674,    0.0000,&
    0.3059,   12.6697,  420.0725,    2.7933,  -10.5186,    0.0000,&
    0.2414,    7.4459,  233.2453,    2.3806,   -3.8121,    0.0000,&
    0.1885,    8.9611,  434.0463,    2.4783,   -3.4805,    0.0000,&
    0.2372,    7.4225,  253.7932,    2.5144,   -4.2873,    0.0000,&
    0.2615,    7.8520,  301.0166,    2.6823,   -5.1720,    0.0000,&
    0.2146,    9.4106,  370.1657,    2.4035,   -4.6659,    0.0000,&
    0.2385,    7.6679,  237.7240,    2.3976,   -3.8943,    0.0000,&
    0.2722,    7.6868,  240.8871,    2.4954,   -4.6169,    0.0000,&
    0.2199,    7.4339,  251.6986,    2.5388,   -4.4520,    0.0000,&
    0.2478,    7.2748,  254.1880,    2.5611,   -4.6149,    0.0000,&
    0.2600,    7.8033,  327.9953,    2.8921,   -5.8782,    0.0000,&
    0.2678,    7.0731,  225.4841,    2.5731,   -4.6740,    0.0000,&
    0.1859,    8.9088,  396.7617,    2.4443,   -4.0500,    0.0000,&
    0.3252,   16.2530,  290.9761,    2.5294,  -13.3920,    0.0000,&
    0.2465,    8.2372,  224.3463,    2.3187,   -4.3525,    0.0000,&
    0.2299,    6.8383,  204.1714,    2.3330,   -3.2798,    0.0000,&
    0.2347,    7.2584,  251.4713,    2.5473,   -4.5875,    0.0000,&
    0.2678,    9.7140,  330.1632,    2.4641,   -5.6242,    0.0000,&
    0.2350,    7.4074,  281.2593,    2.4868,   -3.7386,    0.0000,&
    0.2356,    7.3930,  245.6756,    2.5338,   -4.4752,    0.0000,&
    0.2681,    7.6264,  285.2058,    2.6199,   -4.6296,    0.0000,&
    0.2495,    7.6647,  234.2579,    2.4486,   -4.3203,    0.0000,&
    0.2096,    9.1329,  342.2426,    2.4186,   -3.9152,    0.0000,&
    0.2623,    7.6822,  262.7294,    2.6817,   -5.2354,    0.0000,&
    0.2448,    8.0460,  266.0813,    2.3944,   -3.9883,    0.0000,&
    0.2323,    7.2566,  245.2923,    2.5506,   -4.3476,    0.0000,&
    0.2469,    7.5844,  282.6983,    2.5963,   -4.5497,    0.0000,&
    0.2597,    7.1448,  242.3148,    2.5577,   -4.6377,    0.0000,&
    0.2628,    8.0269,  284.3552,    2.6831,   -5.8550,    0.0000,&
    0.2585,    8.1473,  272.1278,    2.5387,   -5.1721,    0.0000,&
    0.2300,    7.1898,  255.2655,    2.5667,   -4.3783,    0.0000,&
    0.2559,    8.2205,  259.1454,    2.4725,   -4.8577,    0.0000,&
    0.2553,    7.1695,  210.7529,    2.5140,   -4.8163,    0.0000,&
    0.2333,    7.3747,  261.9473,    2.5807,   -4.4241,    0.0000,&
    0.2347,    7.2429,  251.6140,    2.5492,   -4.6801,    0.0000,&
    0.2537,    7.1147,  192.9449,    2.3069,   -3.3640,    0.0000,&
    0.2430,    8.1775,  252.9467,    2.4352,   -4.5367,    0.0000,&
    0.2769,    8.0360,  215.9608,    2.4903,   -5.1542,    0.0000,&
    0.2335,    7.2984,  252.3044,    2.5801,   -4.4854,    0.0000,&
    0.2590,    7.8859,  281.0910,    2.6250,   -5.1784,    0.0000,&
    0.2612,    7.7825,  251.8535,    2.4535,   -4.6412,    0.0000,&
    0.1804,    9.0012,  430.7154,    2.4849,   -4.5794,    0.0000,&
    0.2353,   10.4947,  373.7626,    2.3859,   -5.2160,    0.0000,&
    0.2479,    7.3098,  204.1023,    2.3002,   -3.5652,    0.0000,&
    0.2497,    6.9438,  234.6510,    2.7669,   -5.2940,    0.0000,&
    0.2509,    7.3139,  255.0434,    2.4782,   -3.9361,    0.0000,&
    0.2481,    7.6391,  280.5731,    2.5415,   -4.3075,    0.0000,&
    0.2441,    7.8022,  234.6851,    2.3651,   -4.3189,    0.0000,&
    0.2462,    8.0317,  254.4989,    2.4469,   -4.6550,    0.0000,&
    0.2558,    7.3333,  201.9142,    2.4448,   -4.7429,    0.0000,&
    0.2338,    7.2654,  251.1804,    2.5567,   -4.4925,    0.0000,&
    0.2451,    7.3565,  299.3438,    2.7158,   -4.7219,    0.0000,&
    0.2338,    7.2654,  251.1804,    2.5567,   -4.4925,    0.0000,&
    0.2581,    7.9864,  242.9988,    2.3513,   -4.5648,    0.0000,&
    0.2408,    8.1678,  276.5880,    2.5328,   -4.9945,    0.0000,&
    0.2124,   10.5676,  409.0769,    2.4887,   -6.1128,    0.0000,&
    0.2399,    7.4061,  277.3885,    2.5615,   -4.2775,    0.0000,&
    0.2489,    7.3563,  250.5803,    2.4653,   -4.1817,    0.0000,&
    0.2541,    7.7330,  255.4031,    2.5063,   -4.4257,    0.0000,&
    0.2513,    8.1819,  224.3622,    2.3832,   -4.8166,    0.0000,&
    0.2502,    8.0749,  260.2824,    2.4316,   -4.3011,    0.0000,&
    0.2432,    7.8728,  297.4122,    2.5867,   -4.7330,    0.0000,&
    0.2399,    7.6936,  245.9075,    2.4497,   -4.1629,    0.0000,&
    0.2189,    7.6667,  320.8838,    2.5352,   -3.9767,    0.0000,&
    0.2215,    7.2639,  276.2159,    2.5280,   -3.5505,    0.0000,&
    0.2278,    7.6120,  268.3801,    2.5816,   -4.6081,    0.0000,&
    0.2316,    7.5237,  236.0710,    2.4017,   -4.1679,    0.0000,&
    0.2752,   10.1842,  370.2353,    2.7026,   -7.3806,    0.0000,&
    0.2312,    7.2176,  249.2958,    2.5441,   -4.3637,    0.0000,&
    0.2542,    7.2176,  227.6251,    2.4744,   -4.3065,    0.0000,&
    0.2469,    7.5844,  282.6983,    2.5963,   -4.5497,    0.0000,&
    0.2161,    7.4562,  330.0610,    2.5320,   -3.5504,    0.0000,&
    0.2745,   10.0278,  322.1786,    2.6621,   -7.5642,    0.0000,&
    0.3018,    9.1025,  273.5308,    2.8034,   -7.4178,    0.0000,&
    0.2791,    8.8633,  280.4576,    2.7187,   -6.9765,    0.0000,&
    0.2424,    7.3207,  250.3308,    2.5019,   -4.4059,    0.0000,&
    0.2457,    7.7875,  249.9160,    2.4735,   -4.4884,    0.0000,&
    0.2481,    8.0122,  272.4213,    2.5365,   -4.6672,    0.0000,&
    0.2490,    7.4172,  239.2965,    2.4462,   -4.0470,    0.0000,&
    0.2374,    8.5941,  247.6610,    2.3264,   -4.6401,    0.0000,&
    0.2229,    7.4685,  291.4087,    2.5604,   -4.4464,    0.0000,&
    0.2769,    8.0360,  215.9608,    2.4903,   -5.1542,    0.0000/)
    
    allocate(NParray(1:size(NPparam)))
    if(rep>100) then
        write(*,*) color('ERROR in BSM19.bFIX model. It has only 100 replicas. Central replica is set',c_red)
        NParray=1d0*replicas((0+2)*6+1:(0+2)*6+6)
    else
        NParray=1d0*replicas((rep+2)*6+1:(rep+2)*6+6)
    end if
end subroutine GetReplicaParameters

end module uTMDPDF_model
