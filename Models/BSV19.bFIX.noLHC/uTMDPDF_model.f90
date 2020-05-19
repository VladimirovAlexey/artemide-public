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
    
    write(*,*) color(">>>  The model for uTMDPDF is BSV19.bFIX.noLHC. Please, cite [1902.08474]   <<<",c_cyan)
    
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
    0.135821, 11.2293,   413.805,   2.07758,  -2.47575,    0.1000,&!!!the suggested fo initialization replica is slightly wider. It prevent the 0-values at large b
    0.135821, 11.2293,   413.805,   2.07758,  -2.47575,    0.0000,&!!! best replica
    0.135821, 11.2293,   413.805,   2.07758,  -2.47575,    0.0000,&!!! mean replica
    0.3203,   10.2384,  207.7621,    2.3082,   -6.8334,    0.0000,&
    0.0493,   16.8501,  592.1730,    2.1721,   -6.6425,    0.0000,&
    0.1237,   14.4919,  578.9114,    2.4078,   -6.9248,    0.0000,&
    0.0207,   11.0275,  477.3483,    1.9698,    1.4794,    0.0000,&
    0.2031,   10.4123,  276.3608,    1.9837,   -2.6992,    0.0000,&
    0.2445,   10.8924,  304.4318,    2.1114,   -4.0831,    0.0000,&
    0.0029,    6.2163,  434.2731,    1.9461,    5.5859,    0.0000,&
    0.1901,   12.3710,  382.1757,    2.1382,   -4.4974,    0.0000,&
    0.1892,   14.1663,  400.8676,    2.2314,   -7.2234,    0.0000,&
    0.0482,   13.1583,  660.0553,    2.1065,   -0.7813,    0.0000,&
    0.1105,   11.2779,  429.0225,    2.0587,   -1.1421,    0.0000,&
    0.1688,    8.9537,  254.1412,    2.0661,   -1.7291,    0.0000,&
    0.2437,    9.7672,  257.8334,    2.0467,   -3.1494,    0.0000,&
    0.0007,   14.3272,  467.4429,    1.8914,   -1.0972,    0.0000,&
    0.1473,   10.3249,  289.7746,    2.0884,   -2.7868,    0.0000,&
    0.0050,   23.3423, 1297.6310,    2.4191,  -11.4937,    0.0000,&
    0.1779,   10.5071,  323.1526,    2.0777,   -2.5492,    0.0000,&
    0.0376,   12.8898,  682.5326,    2.2903,   -2.8702,    0.0000,&
    0.0119,   14.8465,  676.7602,    2.1010,   -2.0680,    0.0000,&
    0.1517,   13.1372,  366.0090,    2.0989,   -5.1827,    0.0000,&
    0.2339,    7.4022,  136.4038,    2.0113,   -2.5848,    0.0000,&
    0.2273,    8.0676,  288.4444,    2.1092,   -0.8464,    0.0000,&
    0.2479,    9.9857,  232.2811,    2.0188,   -3.4313,    0.0000,&
    0.2901,   12.1923,  338.1887,    2.2768,   -6.7207,    0.0000,&
    0.1337,   12.5233,  527.2268,    2.3712,   -4.2851,    0.0000,&
    0.1484,   10.1013,  238.0258,    2.1709,   -3.9868,    0.0000,&
    0.2689,   10.3890,  266.9041,    2.1753,   -4.9113,    0.0000,&
    0.2381,    8.5216,  192.5438,    1.9689,   -2.3800,    0.0000,&
    0.1812,   10.4234,  352.3761,    2.1685,   -2.5765,    0.0000,&
    0.0432,   19.0076, 1199.9113,    2.4114,   -7.1005,    0.0000,&
    0.1001,    8.9491,  307.8871,    1.9237,    1.7622,    0.0000,&
    0.0773,   16.5154,  649.9063,    2.2125,   -6.3812,    0.0000,&
    0.0531,   11.9432,  536.7471,    2.0720,    0.1084,    0.0000,&
    0.1444,   13.1165,  374.2656,    2.3050,   -6.5316,    0.0000,&
    0.0265,   13.9315,  580.7939,    2.1942,   -3.3127,    0.0000,&
    0.1101,   15.8651,  618.8660,    2.3243,   -7.0425,    0.0000,&
    0.1571,    7.8505,  211.0611,    1.9222,    0.2671,    0.0000,&
    0.2191,    7.8065,  204.2413,    1.8320,    0.0062,    0.0000,&
    0.1968,    9.5825,  238.8714,    1.9334,   -1.8419,    0.0000,&
    0.0171,    2.4996,  282.0199,    1.4856,   11.7653,    0.0000,&
    0.1486,   13.2377,  429.6052,    2.2579,   -5.6716,    0.0000,&
    0.2093,   10.2020,  264.8255,    2.0950,   -3.5823,    0.0000,&
    0.1947,   12.8304,  424.6658,    2.2492,   -5.5852,    0.0000,&
    0.1252,   10.5713,  330.1720,    2.0319,   -1.7631,    0.0000,&
    0.1689,    8.1907,  241.7567,    2.0172,   -0.6706,    0.0000,&
    0.2085,    6.1717,  147.0558,    1.7287,    1.2990,    0.0000,&
    0.0718,   17.5079,  718.9594,    2.2407,   -6.7717,    0.0000,&
    0.1775,   10.6615,  265.7828,    2.0743,   -3.6132,    0.0000,&
    0.1474,   10.4157,  270.5497,    1.9825,   -2.4872,    0.0000,&
    0.0169,   10.0866,  394.3153,    1.8101,    3.9096,    0.0000,&
    0.2119,   10.9424,  286.7716,    2.1365,   -4.4780,    0.0000,&
    0.1868,    7.4688,  140.6990,    1.9544,   -1.8862,    0.0000,&
    0.0840,   14.3941,  517.9648,    2.0989,   -3.4241,    0.0000,&
    0.1654,   10.1378,  285.7086,    2.0961,   -2.8246,    0.0000,&
    0.2130,   11.6147,  370.5032,    2.1986,   -4.3690,    0.0000,&
    0.0646,   13.2442,  556.2879,    2.1842,   -2.8382,    0.0000,&
    0.0553,   10.8935,  520.9317,    2.1090,   -0.1409,    0.0000,&
    0.1893,    7.1859,  174.6980,    1.8656,   -0.1360,    0.0000,&
    0.2051,    6.3954,  166.7643,    1.8001,    1.2241,    0.0000,&
    0.1921,   10.4895,  289.7872,    1.9428,   -2.0777,    0.0000,&
    0.1301,   14.8527,  502.2133,    2.0326,   -4.4746,    0.0000,&
    0.1696,   10.2090,  265.0910,    2.1207,   -3.4064,    0.0000,&
    0.0148,   16.7007,  810.9182,    2.2179,   -4.6163,    0.0000,&
    0.0845,   10.8626,  300.5632,    1.9383,   -1.5587,    0.0000,&
    0.1469,   10.6171,  375.3059,    2.2669,   -3.5992,    0.0000,&
    0.1682,    6.8976,  221.2041,    1.8445,    1.6578,    0.0000,&
    0.1158,   12.8450,  371.3145,    2.1803,   -5.2679,    0.0000,&
    0.0389,    4.5646,  230.3756,    1.5515,    8.0539,    0.0000,&
    0.2215,    9.9380,  267.9917,    2.1177,   -3.4787,    0.0000,&
    0.0005,   22.3728, 1321.8078,    2.3598,   -9.5510,    0.0000,&
    0.1021,   13.4773,  462.5500,    2.1250,   -4.0534,    0.0000,&
    0.0577,   19.6361, 1548.8614,    2.5375,   -9.3576,    0.0000,&
    0.1318,    9.3830,  264.8142,    2.0268,   -1.3226,    0.0000,&
    0.1853,   11.7121,  371.1077,    2.2774,   -5.1697,    0.0000,&
    0.1590,    8.7819,  229.4285,    1.8931,   -0.4093,    0.0000,&
    0.2420,    8.6059,  177.8044,    1.8611,   -2.0407,    0.0000,&
    0.0278,   10.2289,  462.6457,    1.9399,    2.6371,    0.0000,&
    0.0743,   24.5931, 1626.2317,    2.7417,  -15.6620,    0.0000,&
    0.1532,   15.0545,  410.3061,    2.1522,   -7.4095,    0.0000,&
    0.0151,   10.9901,  463.4220,    1.9850,    0.9882,    0.0000,&
    0.0152,   12.8055,  507.5059,    2.0400,   -0.7271,    0.0000,&
    0.0701,    5.7273,  179.1455,    1.6772,    4.1556,    0.0000,&
    0.2215,    7.8693,  209.3558,    1.8526,   -0.0502,    0.0000,&
    0.0723,   13.3372,  523.3329,    2.1085,   -2.8653,    0.0000,&
    0.1855,    8.0360,  205.0965,    1.9787,   -1.0894,    0.0000,&
    0.0270,    5.0911,  285.8459,    1.7217,    7.3022,    0.0000,&
    0.2739,    9.1446,  274.5006,    2.4015,   -4.7552,    0.0000,&
    0.2262,    9.4022,  210.9295,    2.0552,   -3.5645,    0.0000,&
    0.1381,    8.7526,  274.4438,    2.1654,   -1.9353,    0.0000,&
    0.0554,   11.6893,  418.7567,    1.9656,   -0.4025,    0.0000,&
    0.1467,   10.2205,  326.6308,    2.0112,   -1.5362,    0.0000,&
    0.1068,    6.0458,  272.5422,    1.7448,    5.6602,    0.0000,&
    0.0769,    9.2644,  380.3605,    2.0682,    0.5898,    0.0000,&
    0.0373,   15.9498,  825.4843,    2.5184,   -7.1001,    0.0000,&
    0.2810,    7.8986,  173.4535,    1.9945,   -2.4347,    0.0000,&
    0.0202,   12.8427,  588.5290,    2.1794,   -1.6503,    0.0000,&
    0.1398,   14.5789,  479.6105,    2.2915,   -6.8369,    0.0000,&
    0.2317,   10.1883,  284.3366,    1.9017,   -1.3917,    0.0000,&
    0.2050,    6.6016,  162.1046,    1.9601,   -0.5491,    0.0000,&
    0.2147,    6.6506,  181.1879,    1.6776,    2.2456,    0.0000/)

    allocate(NParray(1:size(NPparam)))
    if(rep>100) then
        write(*,*) color('ERROR in BSM19.bFIX.noLHC model. It has only 100 replicas. Central replica is set',c_red)
        NParray=1d0*replicas((0+2)*6+1:(0+2)*6+6)
    else
        NParray=1d0*replicas((rep+2)*6+1:(rep+2)*6+6)
    end if
end subroutine GetReplicaParameters

end module uTMDPDF_model
