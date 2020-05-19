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
    
    write(*,*) color(">>>  The model for uTMDPDF is BSV19.bFIT.noLHC. Please, cite [1902.08474]   <<<",c_cyan)
    
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
 real,parameter,dimension(1:588)::replicas=(/&
    0.212239, 12.0979, 316.472, 2.10881, -3.50579,   0.1000,&!!!the suggested fo initialization replica is slightly wider. It prevent the 0-values at large b
    0.2511,    7.9263,  297.3788,    2.4499,   -4.3889,    0.0000,&!!! best replica
    0.212239, 12.0979, 316.472, 2.10881, -3.50579,   0.0000,&!!! mean replica
    0.0067,   14.8533,  464.5722,    2.2732,   -3.7824,    0.0000, &
    0.0416,   12.6174,  492.7889,    2.4393,   -3.3180,    0.0000, &
    0.2701,   10.3854,  157.7149,    2.0810,   -4.1190,    0.0000, &
    0.0000,   12.0892,  698.8302,    2.1145,   -0.4552,    0.0000, &
    0.2362,   17.2951,  412.0402,    2.0178,   -5.2196,    0.0000, &
    0.1101,   11.5238,  326.7139,    2.1229,   -1.9468,    0.0000, &
    0.0558,   16.4995,  577.8282,    2.1661,   -3.1094,    0.0000, &
    0.2037,   10.7452,  199.7868,    2.0716,   -2.8553,    0.0000, &
    0.2142,    9.8104,  268.5611,    2.3835,   -2.5687,    0.0000, &
    0.0821,    4.9439,  223.4503,    1.6957,    6.1199,    0.0000, &
    0.2017,   14.5086,  352.4795,    2.2961,   -6.2614,    0.0000, &
    0.2721,    8.9621,  177.9451,    2.3624,   -3.7332,    0.0000, &
    0.0001,   11.7945,  360.8657,    2.3563,   -2.8348,    0.0000, &
    0.3707,    8.6164,  116.6007,    1.7083,   -0.4140,    0.0000, &
    0.2060,    9.3983,  224.0276,    1.9284,    0.4171,    0.0000, &
    0.0163,   15.2515,  710.0887,    2.1771,   -2.9532,    0.0000, &
    0.1704,   19.9877,  776.5426,    2.9720,  -15.0580,    0.0000, &
    0.4846,   28.1113,  580.8608,    2.4921,  -22.9855,    0.0000, &
    0.2059,   14.5532,  353.8423,    2.2895,   -6.2416,    0.0000, &
    0.4381,   11.8621,  157.1466,    1.8570,   -4.3223,    0.0000, &
    0.2942,    6.8216,  108.9433,    1.7964,    0.4898,    0.0000, &
    0.2207,   12.5298,  229.1676,    1.9690,   -3.5256,    0.0000, &
    0.0968,   12.6638,  343.7903,    2.2932,   -3.7287,    0.0000, &
    0.4647,   14.5239,  180.9930,    2.1998,  -10.4474,    0.0000, &
    0.4775,    8.8428,   81.6143,    1.7338,   -3.6704,    0.0000, &
    0.1963,    7.5824,  146.9405,    1.7797,    1.9244,    0.0000, &
    0.1972,    9.4792,  214.2217,    1.9363,   -1.0219,    0.0000, &
    0.3971,   10.7637,  147.7962,    1.8294,   -2.8185,    0.0000, &
    0.3635,   10.7086,  160.7138,    1.9155,   -3.1188,    0.0000, &
    0.1243,    8.0102,  219.3533,    2.1257,    0.5106,    0.0000, &
    0.1073,    9.2102,  282.3304,    1.8771,    1.5500,    0.0000, &
    0.4373,    8.6551,  103.3083,    1.7065,   -1.3936,    0.0000, &
    0.0263,    8.9893,  243.5721,    1.9648,    2.1485,    0.0000, &
    0.4202,   12.1880,  171.2298,    1.9023,   -4.8085,    0.0000, &
    0.1104,   11.6284,  306.6160,    2.1077,   -1.6050,    0.0000, &
    0.0698,    6.9692,  217.9329,    1.8444,    5.1171,    0.0000, &
    0.3133,   13.7336,  444.5688,    2.2658,   -3.5056,    0.0000, &
    0.0836,   10.0169,  269.5349,    2.2588,   -1.7469,    0.0000, &
    0.2390,   11.9522,  255.1875,    2.0189,   -3.2775,    0.0000, &
    0.0042,    9.7648,  360.1134,    2.0836,    2.5893,    0.0000, &
    0.0379,   18.7177,  889.3145,    2.4971,   -7.8109,    0.0000, &
    0.0483,    9.9669,  329.8418,    2.0594,    2.0999,    0.0000, &
    0.0954,   12.5751,  321.5337,    2.0810,   -1.2180,    0.0000, &
    0.5377,   15.9378,  212.1660,    2.0064,   -9.7418,    0.0000, &
    0.2137,    9.4561,  174.4836,    2.0499,   -2.0848,    0.0000, &
    0.0015,   13.3768,  440.3067,    2.3037,   -2.9843,    0.0000, &
    0.0625,   12.0921,  372.9519,    2.2062,   -1.0567,    0.0000, &
    0.0140,   13.8207,  410.5125,    2.4410,   -5.3891,    0.0000, &
    0.4255,    8.5804,   90.9257,    1.7794,   -2.4737,    0.0000, &
    0.1014,    8.0921,  375.2175,    1.8703,    2.2608,    0.0000, &
    0.0000,    9.2265,  258.0443,    1.9649,    2.4215,    0.0000, &
    0.2198,   14.2446,  358.0585,    2.2791,   -5.3996,    0.0000, &
    0.4659,    9.1023,  104.6867,    1.7482,   -2.6464,    0.0000, &
    0.1998,   13.9786,  344.4110,    2.1050,   -4.4856,    0.0000, &
    0.3500,   10.8110,  144.8562,    1.8629,   -3.4716,    0.0000, &
    0.5025,   10.2509,  116.3675,    1.6956,   -2.9027,    0.0000, &
    0.5167,    8.4327,   70.3555,    1.8710,   -4.9650,    0.0000, &
    0.5102,   10.2035,  109.5501,    1.7919,   -4.0216,    0.0000, &
    0.0316,   13.7774,  882.8922,    2.4899,   -2.7168,    0.0000, &
    0.0057,    9.8799,  311.7494,    2.1020,    0.7235,    0.0000, &
    0.2960,   25.9711,  424.6464,    2.2205,  -17.0888,    0.0000, &
    0.0694,   11.6933,  424.5579,    2.0984,   -0.1476,    0.0000, &
    0.0485,   12.1395,  352.5038,    2.2729,   -2.5808,    0.0000, &
    0.0609,   10.6654,  365.0554,    2.1336,   -1.2117,    0.0000, &
    0.1689,   16.7258,  444.8945,    2.3004,   -7.5742,    0.0000, &
    0.0840,    9.0688,  303.5374,    2.5126,   -1.8657,    0.0000, &
    0.0263,    8.9893,  243.5718,    1.9648,    2.1485,    0.0000, &
    0.1275,    9.1275,  229.8819,    2.5473,   -3.8138,    0.0000, &
    0.1220,   14.5082,  420.0979,    2.2844,   -4.6837,    0.0000, &
    0.2042,   11.6645,  247.6604,    2.0366,   -2.6196,    0.0000, &
    0.1808,   12.8893,  327.0427,    2.1013,   -2.7011,    0.0000, &
    0.5498,    9.8506,   91.5100,    1.7134,   -4.1896,    0.0000, &
    0.2804,   12.2557,  278.3116,    1.9123,   -1.0872,    0.0000, &
    0.1199,    9.5396,  298.2181,    1.9596,    0.3167,    0.0000, &
    0.4030,   14.3694,  255.9111,    2.1225,   -6.9087,    0.0000, &
    0.2105,    9.1985,  152.9791,    2.2422,   -3.9061,    0.0000, &
    0.4660,    7.5412,   86.3213,    1.6188,   -0.4537,    0.0000, &
    0.1907,   10.7971,  231.0001,    1.9682,   -1.0449,    0.0000, &
    0.2856,    9.5708,  164.9142,    1.9514,   -1.6343,    0.0000, &
    0.2870,    7.5709,  109.3652,    1.7737,   -0.1974,    0.0000, &
    0.4086,   16.4230,  239.4412,    2.1459,  -10.4953,    0.0000, &
    0.0023,   10.9363,  869.2940,    2.2051,    1.3586,    0.0000, &
    0.4792,    7.2244,   64.7335,    1.5808,   -1.4534,    0.0000, &
    0.2570,   11.8175,  238.7962,    2.0440,   -3.3061,    0.0000, &
    0.0879,   14.3841,  493.2838,    2.3550,   -4.1936,    0.0000, &
    0.0522,   17.1489,  672.7540,    2.6564,   -9.0915,    0.0000, &
    0.0181,   12.3796,  501.7183,    2.4887,   -4.1720,    0.0000, &
    0.0986,   12.2682,  425.6635,    2.2247,   -2.8220,    0.0000, &
    0.1710,   10.8345,  234.5456,    2.1679,   -3.0833,    0.0000, &
    0.1751,   15.0240,  406.7704,    2.2662,   -5.9840,    0.0000, &
    0.0000,   16.1200,  464.7066,    2.0745,   -3.6718,    0.0000, &
    0.3185,   11.7789,  171.2201,    2.1304,   -6.1484,    0.0000, &
    0.6955,   35.7856,  999.6214,    3.1488,  -37.1290,    0.0000, &
    0.1890,    9.4499,  212.1409,    2.5492,   -4.7331,    0.0000, &
    0.4581,   12.2188,  179.3421,    1.9457,   -5.0383,    0.0000/)
    
    allocate(NParray(1:size(NPparam)))
    if(rep>95) then
        write(*,*) color('ERROR in BSM19.bFIT.noLHC model. It has only 95 replicas. Central replica is set',c_red)
        NParray=1d0*replicas((0+2)*6+1:(0+2)*6+6)
    else
        NParray=1d0*replicas((rep+2)*6+1:(rep+2)*6+6)
    end if
end subroutine GetReplicaParameters

end module uTMDPDF_model
