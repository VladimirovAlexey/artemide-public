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
    
    write(*,*) color(">>>  The model for uTMDPDF is BSV19.NNPDF31. Please, cite [1902.08474]   <<<",c_cyan)
    
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
    0.253434,  9.04351, 346.9985,    2.47991,  -5.69988,   0.1, & !!!! initialisation
    0.2512,    7.7572,  334.61,      2.4543,   -4.8203,    0.,&     !!! best
    0.253434,  9.04351, 346.9985,    2.47991,  -5.69988,   0.0,&      !!!!! mean
    0.2578,    7.8679,  340.3738,    2.4602,   -4.7937,    0.0000,&
    0.2578,    7.9417,  331.1954,    2.3978,   -4.6781,    0.0000,&
    0.2586,    7.9660,  344.2272,    2.4611,   -4.8205,    0.0000,&
    0.2289,    7.2611,  258.1676,    2.5300,   -4.4596,    0.0000,&
    0.2478,    8.9980,  354.3787,    2.4011,   -5.4446,    0.0000,&
    0.2572,    8.0662,  359.6624,    2.4991,   -4.7360,    0.0000,&
    0.2534,    9.0003,  314.0943,    2.3367,   -4.8661,    0.0000,&
    0.2495,    7.7463,  335.2263,    2.4628,   -4.7678,    0.0000,&
    0.2597,    8.1363,  346.5956,    2.4724,   -4.7693,    0.0000,&
    0.2698,    8.7243,  186.1754,    2.2774,   -5.2725,    0.0000,&
    0.2429,    8.3242,  359.3842,    2.5081,   -4.5756,    0.0000,&
    0.2366,    8.7685,  234.0162,    2.3468,   -5.5331,    0.0000,&
    0.2943,    7.6210,  216.0706,    2.4311,   -5.0584,    0.0000,&
    0.2545,    7.8254,  339.7985,    2.4825,   -4.6657,    0.0000,&
    0.2566,    7.8651,  335.5309,    2.4516,   -4.7690,    0.0000,&
    0.3022,    8.9906,  235.8987,    2.4546,   -6.3597,    0.0000,&
    0.2924,    7.8132,  280.1988,    2.5360,   -5.4671,    0.0000,&
    0.2842,    8.4124,  243.1678,    2.4213,   -5.8005,    0.0000,&
    0.2509,    7.8577,  340.7542,    2.4789,   -4.6714,    0.0000,&
    0.2552,    8.5353,  360.8098,    2.4516,   -4.5495,    0.0000,&
    0.2583,    9.0322,  324.2326,    2.3339,   -5.3440,    0.0000,&
    0.2534,    7.7721,  335.8151,    2.4534,   -4.8000,    0.0000,&
    0.2806,    7.8008,  253.8902,    2.5131,   -5.4256,    0.0000,&
    0.2080,    8.4711,  431.5971,    2.3778,   -3.8581,    0.0000,&
    0.2927,    7.1473,  189.9989,    2.3611,   -4.8752,    0.0000,&
    0.2004,    9.0091,  400.3620,    2.5560,   -5.2387,    0.0000,&
    0.2486,    8.1903,  311.0446,    2.4497,   -4.5668,    0.0000,&
    0.2474,    8.7945,  358.3973,    2.4468,   -5.2981,    0.0000,&
    0.2562,    8.1847,  348.7013,    2.4944,   -4.7059,    0.0000,&
    0.2441,    8.4432,  352.2815,    2.4740,   -4.5574,    0.0000,&
    0.2377,    8.2127,  350.1083,    2.4668,   -4.1377,    0.0000,&
    0.2490,    8.8928,  245.4520,    2.3562,   -5.5560,    0.0000,&
    0.2206,    9.1749,  374.2333,    2.3923,   -4.6236,    0.0000,&
    0.2971,    7.9618,  279.8325,    2.6157,   -6.2620,    0.0000,&
    0.2382,    9.1020,  429.3304,    2.4859,   -5.0293,    0.0000,&
    0.2571,    7.9541,  339.8995,    2.4530,   -4.7595,    0.0000,&
    0.2376,    8.6221,  370.9768,    2.4076,   -4.7857,    0.0000,&
    0.2559,    8.8900,  259.1895,    2.3447,   -5.1052,    0.0000,&
    0.2546,    8.4668,  360.2347,    2.5016,   -4.9166,    0.0000,&
    0.1296,   10.9171,  669.8408,    2.4851,   -4.6548,    0.0000,&
    0.2022,    9.4991,  432.8895,    2.3838,   -4.4407,    0.0000,&
    0.2351,    9.2695,  317.0773,    2.5276,   -6.2270,    0.0000,&
    0.2828,   10.1168,  310.8260,    2.5185,   -7.9906,    0.0000,&
    0.2297,    9.0575,  392.3475,    2.3915,   -4.5366,    0.0000,&
    0.2615,    8.1249,  340.0968,    2.4651,   -4.9240,    0.0000,&
    0.2855,   16.3130,  551.0105,    2.7064,  -13.2152,    0.0000,&
    0.3010,   12.1228,  430.4172,    2.8512,  -10.4697,    0.0000,&
    0.2521,    7.7834,  334.5582,    2.4489,   -4.8330,    0.0000,&
    0.3045,    9.0159,  258.0490,    2.6019,   -6.5603,    0.0000,&
    0.2925,   13.2559,  462.8200,    2.7862,  -11.3629,    0.0000,&
    0.1835,    7.9779,  438.9639,    2.3970,   -2.4262,    0.0000,&
    0.2645,    7.7938,  253.2093,    2.4161,   -4.6657,    0.0000,&
    0.2595,    8.2350,  346.8596,    2.4570,   -5.0000,    0.0000,&
    0.2041,    8.1963,  479.3535,    2.3780,   -3.0293,    0.0000,&
    0.2513,    7.7859,  344.7852,    2.4961,   -4.4296,    0.0000,&
    0.2550,    8.2813,  345.4170,    2.4364,   -4.6987,    0.0000,&
    0.2272,    7.9506,  330.5669,    2.3294,   -3.7019,    0.0000,&
    0.2498,    7.9273,  343.8844,    2.4893,   -4.5797,    0.0000,&
    0.2544,    7.9418,  339.6894,    2.4547,   -4.7372,    0.0000,&
    0.2471,    8.1384,  350.1716,    2.4928,   -4.7155,    0.0000,&
    0.3001,    9.8060,  270.0585,    2.4990,   -7.3311,    0.0000,&
    0.2782,    8.1525,  356.6410,    2.6478,   -5.9318,    0.0000,&
    0.2543,    8.0337,  362.0907,    2.4635,   -4.7212,    0.0000,&
    0.2676,   10.4287,  275.2772,    2.4603,   -7.2230,    0.0000,&
    0.2389,    8.8335,  345.6981,    2.4183,   -4.8867,    0.0000,&
    0.2282,    9.7983,  358.9101,    2.4495,   -6.1289,    0.0000,&
    0.2554,    8.0037,  343.1728,    2.4677,   -4.6594,    0.0000,&
    0.2410,    8.9606,  287.8230,    2.3621,   -4.7543,    0.0000,&
    0.2519,    7.9599,  343.5066,    2.5134,   -4.5314,    0.0000,&
    0.1858,    7.3644,  279.9578,    2.4554,   -3.4198,    0.0000,&
    0.2988,   10.1213,  303.4911,    2.5806,   -7.5968,    0.0000,&
    0.2286,    8.1397,  312.5855,    2.3386,   -3.9624,    0.0000,&
    0.2747,    8.8846,  234.3163,    2.4776,   -5.6714,    0.0000,&
    0.2537,    8.1837,  352.7137,    2.4894,   -4.7056,    0.0000,&
    0.2436,   10.8767,  453.1025,    2.7338,   -8.4455,    0.0000,&
    0.2423,    9.4442,  272.1055,    2.3025,   -5.8305,    0.0000,&
    0.2195,    8.2628,  406.3514,    2.7204,   -5.3865,    0.0000,&
    0.2699,    8.7392,  259.9923,    2.3772,   -5.2352,    0.0000,&
    0.2550,   12.7960,  441.4030,    2.7031,  -10.0856,    0.0000,&
    0.2553,    8.0415,  345.8693,    2.5498,   -4.7970,    0.0000,&
    0.3830,   35.2749, 1348.8256,    3.3968,  -35.3014,    0.0000,&
    0.2522,    7.8991,  335.8411,    2.4767,   -4.6685,    0.0000,&
    0.2410,    9.5349,  322.5648,    2.5293,   -5.3318,    0.0000,&
    0.2204,    8.5627,  351.9208,    2.3441,   -3.8169,    0.0000,&
    0.3168,   12.7056,  494.5550,    2.9003,  -11.8234,    0.0000,&
    0.2989,    7.0566,  216.3749,    2.4572,   -5.3794,    0.0000,&
    0.2252,    7.9162,  337.7414,    2.3568,   -3.4441,    0.0000,&
    0.2501,    7.8795,  340.0092,    2.4723,   -4.6812,    0.0000,&
    0.2912,    9.5952,  279.6307,    2.5727,   -7.3031,    0.0000,&
    0.2245,    8.3983,  303.5963,    2.3973,   -4.2698,    0.0000,&
    0.2103,    9.5049,  418.9691,    2.3943,   -4.8298,    0.0000,&
    0.2460,    8.2005,  343.4480,    2.4185,   -4.9333,    0.0000,&
    0.2314,    8.1926,  345.9195,    2.3547,   -4.1965,    0.0000,&
    0.2928,    7.7296,  238.4706,    2.4851,   -5.2651,    0.0000,&
    0.2332,   10.8998,  447.9679,    2.5406,   -7.5531,    0.0000,&
    0.2389,   10.2660,  371.2749,    2.6375,   -7.2711,    0.0000,&
    0.2707,    8.8445,  248.3621,    2.3799,   -4.9489,    0.0000,&
    0.2567,    7.8326,  324.3136,    2.4470,   -4.7568,    0.0000,&
    0.2454,    9.2014,  319.3286,    2.3263,   -4.8480,    0.0000,&
    0.2512,    8.4771,  201.5357,    2.2364,   -4.9596,    0.0000/)
    
    allocate(NParray(1:size(NPparam)))
    if(rep>100) then
        write(*,*) color('ERROR in BSM19.NNPDF31 model. It has only 100 replicas. Central replica is set',c_red)
        NParray=1d0*replicas((0+2)*6+1:(0+2)*6+6)
    else
        NParray=1d0*replicas((rep+2)*6+1:(rep+2)*6+6)
    end if
end subroutine GetReplicaParameters

end module uTMDPDF_model
