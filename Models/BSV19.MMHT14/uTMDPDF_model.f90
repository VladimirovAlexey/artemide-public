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
    
    write(*,*) color(">>>  The model for uTMDPDF is BSV19.MMHT14. Please, cite [1902.08474]   <<<",c_cyan)
    
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
    0.1982,   26.49689,2727.9766,    3.00668, -23.54749,   0.1, & !!!! initialisation
    0.2188,   25.8465, 1261.4028,    2.9991,  -23.4391,    0.,&     !!! best
    0.1982,   26.49689,2727.9766,    3.00668, -23.54749,   0.,&   !!!!! mean
    0.1910,   22.7207, 2328.0052,    2.9231,  -19.7054,    0.0000,&
    0.1988,   25.1090, 1843.9161,    2.9226,  -21.7653,    0.0000,&
    0.2359,   25.8717, 1630.1565,    2.9027,  -23.3188,    0.0000,&
    0.2194,   25.2122, 2584.1378,    3.0006,  -22.6334,    0.0000,&
    0.2407,   26.7186, 2046.7754,    3.0317,  -24.4652,    0.0000,&
    0.2153,   27.4082, 2932.6219,    3.2082,  -25.4194,    0.0000,&
    0.2167,   25.4598, 2700.2056,    3.0087,  -22.8429,    0.0000,&
    0.2388,   25.5966, 2586.1933,    2.9646,  -23.5357,    0.0000,&
    0.2762,   28.4352, 2679.6622,    3.1084,  -27.2853,    0.0000,&
    0.2108,   25.8419, 2777.0270,    3.0310,  -22.7908,    0.0000,&
    0.1063,   19.0987, 1763.3384,    2.6379,  -13.7639,    0.0000,&
    0.1774,   22.8804, 2687.9651,    3.2070,  -21.2154,    0.0000,&
    0.1629,   19.4402, 1326.3331,    2.6874,  -15.0715,    0.0000,&
    0.2057,   24.6420, 2829.6567,    3.1751,  -23.0156,    0.0000,&
    0.2170,   24.7671, 1659.1057,    2.8549,  -21.4309,    0.0000,&
    0.2199,   26.0642, 2608.8946,    3.2368,  -23.7862,    0.0000,&
    0.2167,   26.6967, 2610.2571,    3.1597,  -25.0978,    0.0000,&
    0.1714,   27.7112, 2697.1590,    2.9929,  -23.4616,    0.0000,&
    0.1976,   27.8737, 3666.7138,    3.2508,  -25.3151,    0.0000,&
    0.1398,   25.1150, 3876.4571,    3.1880,  -22.4828,    0.0000,&
    0.2169,   25.6809, 1851.3858,    2.8325,  -22.1727,    0.0000,&
    0.2154,   25.7077, 2653.4720,    3.1652,  -23.6333,    0.0000,&
    0.1475,   19.6888, 1660.8649,    2.7763,  -15.0173,    0.0000,&
    0.1725,   19.5812, 1400.0339,    2.8411,  -16.3515,    0.0000,&
    0.1329,   25.5504, 3751.3641,    3.0751,  -20.8883,    0.0000,&
    0.2223,   26.4270, 2860.7185,    3.1820,  -24.5229,    0.0000,&
    0.2085,   26.8597, 2978.8109,    3.3056,  -24.3943,    0.0000,&
    0.1333,   24.2930, 2590.3395,    2.9977,  -20.8059,    0.0000,&
    0.2164,   26.0054, 1662.2964,    2.9656,  -22.9146,    0.0000,&
    0.0685,   28.1664, 5590.1096,    3.1098,  -23.6114,    0.0000,&
    0.2737,   24.0169, 1085.8151,    2.8720,  -21.6129,    0.0000,&
    0.2149,   25.6053, 2943.8848,    3.2346,  -23.4147,    0.0000,&
    0.1115,   38.1959, 7377.9210,    3.2227,  -34.6465,    0.0000,&
    0.1834,   26.6079, 3017.1364,    3.1403,  -24.3965,    0.0000,&
    0.1725,   26.6236, 2623.1498,    2.8798,  -22.2063,    0.0000,&
    0.1929,   25.6000, 2288.6458,    2.8122,  -21.3072,    0.0000,&
    0.2820,   26.5104, 2392.5438,    3.1346,  -25.8282,    0.0000,&
    0.1649,   27.9676, 3591.2961,    2.9948,  -23.6516,    0.0000,&
    0.2099,   25.8422, 2547.2389,    3.1359,  -23.5079,    0.0000,&
    0.1949,   29.4788, 3014.4965,    3.0647,  -25.9257,    0.0000,&
    0.2014,   26.0302, 2815.0881,    3.0210,  -23.4841,    0.0000,&
    0.1878,   29.7663, 3164.5021,    3.1950,  -27.3870,    0.0000,&
    0.1912,   24.6995, 2878.5706,    2.9546,  -21.5680,    0.0000,&
    0.2171,   26.2852, 1326.3382,    2.7431,  -22.7946,    0.0000,&
    0.1389,   23.0565, 2314.5330,    2.8431,  -18.0990,    0.0000,&
    0.1464,   27.3768, 3130.9672,    2.8337,  -22.1941,    0.0000,&
    0.2028,   25.0567, 2386.3284,    2.9112,  -22.0686,    0.0000,&
    0.1946,   23.7625, 1896.2752,    2.7943,  -20.1484,    0.0000,&
    0.1422,   34.8494, 8781.6566,    3.4016,  -32.4372,    0.0000,&
    0.2144,   26.7375, 2597.9936,    2.9740,  -24.0063,    0.0000,&
    0.1862,   28.3604, 2325.9362,    3.0367,  -25.3919,    0.0000,&
    0.2054,   26.7665, 3688.9535,    3.1506,  -24.2071,    0.0000,&
    0.2739,   27.7985, 2361.0847,    3.1106,  -26.8771,    0.0000,&
    0.2076,   26.0467, 2794.9946,    2.9725,  -23.1104,    0.0000,&
    0.2130,   26.0827, 2409.3312,    2.9105,  -23.3462,    0.0000,&
    0.1929,   24.5137, 2144.1505,    2.9634,  -21.8433,    0.0000,&
    0.1641,   20.6289, 1560.7148,    2.8473,  -16.9726,    0.0000,&
    0.2130,   24.7077, 2049.6442,    2.9099,  -21.7549,    0.0000,&
    0.2178,   26.0058, 2482.3028,    3.1655,  -23.8934,    0.0000,&
    0.2085,   24.9061, 2494.1352,    2.9491,  -21.9187,    0.0000,&
    0.1955,   23.7207, 1836.3298,    2.8684,  -20.2776,    0.0000,&
    0.1993,   24.7469, 2043.3158,    2.8001,  -21.2718,    0.0000,&
    0.1747,   25.7916, 2487.7728,    2.8677,  -21.4821,    0.0000,&
    0.2146,   25.2828, 2540.4569,    2.9866,  -22.9429,    0.0000,&
    0.1436,   29.8854, 3108.2759,    2.8813,  -24.4397,    0.0000,&
    0.2131,   25.5520, 2696.1895,    3.0007,  -22.5795,    0.0000,&
    0.2121,   25.5775, 2543.1810,    2.9957,  -22.9725,    0.0000,&
    0.2138,   24.4087, 1841.7114,    2.8515,  -21.0144,    0.0000,&
    0.2530,   24.9004, 2097.8796,    2.9338,  -22.9629,    0.0000,&
    0.1642,   30.1542, 3752.7068,    3.1268,  -27.1872,    0.0000,&
    0.2141,   25.5856, 2334.1935,    2.9294,  -22.3717,    0.0000,&
    0.1420,   19.7515, 1832.6418,    2.9153,  -16.2475,    0.0000,&
    0.2251,   27.2092, 2499.6147,    3.0490,  -25.0283,    0.0000,&
    0.1863,   28.2802, 3862.9637,    3.3101,  -26.2006,    0.0000,&
    0.2111,   24.3060, 2105.4692,    2.8705,  -20.7855,    0.0000,&
    0.2181,   26.4717, 1782.0887,    3.0272,  -23.3240,    0.0000,&
    0.2137,   25.0322, 2510.5623,    2.9535,  -22.1704,    0.0000,&
    0.1480,   20.9132, 1485.8545,    2.7328,  -16.2797,    0.0000,&
    0.1282,   25.3007, 3005.5136,    3.0384,  -21.8455,    0.0000,&
    0.2308,   29.7833, 2902.8612,    3.2787,  -28.7431,    0.0000,&
    0.2309,   45.5344, 6834.5817,    3.4664,  -45.4629,    0.0000,&
    0.2067,   24.9896, 2658.0033,    3.0966,  -22.6161,    0.0000,&
    0.1182,   16.5204, 1507.7212,    2.7072,  -12.1525,    0.0000,&
    0.2164,   26.0054, 1662.2964,    2.9656,  -22.9146,    0.0000,&
    0.2245,   26.1088, 1610.8796,    3.0073,  -22.6779,    0.0000,&
    0.2512,   24.5707, 1217.2099,    2.8863,  -22.3076,    0.0000,&
    0.2084,   26.4168, 2007.7002,    2.8691,  -23.2119,    0.0000,&
    0.2335,   26.6455, 2650.3725,    3.2001,  -25.0382,    0.0000,&
    0.2252,   26.2087, 1521.1289,    2.8355,  -23.6286,    0.0000,&
    0.2590,   26.8327, 1918.7637,    3.1876,  -25.9482,    0.0000,&
    0.0994,   34.8991, 6237.6701,    3.0474,  -28.7698,    0.0000,&
    0.2443,   56.8115, 7754.3657,    3.2999,  -55.9696,    0.0000,&
    0.2462,   25.9317, 1425.6144,    2.9568,  -23.5781,    0.0000,&
    0.2219,   26.2578, 1931.5986,    3.1142,  -23.4358,    0.0000,&
    0.2083,   38.5518, 4525.6925,    3.2750,  -36.7051,    0.0000,&
    0.2251,   25.2892, 1923.3761,    2.9437,  -22.3315,    0.0000,&
    0.2186,   25.3225, 1757.9886,    2.9653,  -22.8389,    0.0000,&
    0.1511,   21.1521, 1638.8923,    2.7002,  -16.9010,    0.0000,&
    0.1667,   26.4678, 2212.3378,    2.8169,  -21.9341,    0.0000,&
    0.2108,   25.8419, 2777.0270,    3.0310,  -22.7908,    0.0000/)
    
    allocate(NParray(1:size(NPparam)))
    if(rep>100) then
        write(*,*) color('ERROR in BSM19.MMHT14 model. It has only 100 replicas. Central replica is set',c_red)
        NParray=1d0*replicas((0+2)*6+1:(0+2)*6+6)
    else
        NParray=1d0*replicas((rep+2)*6+1:(rep+2)*6+6)
    end if
end subroutine GetReplicaParameters

end module uTMDPDF_model
