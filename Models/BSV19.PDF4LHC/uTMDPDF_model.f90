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
    
    write(*,*) color(">>>  The model for uTMDPDF is BSV19.PDF4LHC. Please, cite [1902.08474]   <<<",c_cyan)
    
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
    0.2181,   17.9138, 926.08,       2.5431,  -15.5469,   0.1, & !!!! initialisation
    0.2348,   15.4548, 664.81,       2.4141,  -13.7209,    0.,&     !!! best
    0.2181,   17.9138, 926.08,       2.5431,  -15.5469,   0.,&   !!!!! mean
    0.2310,   15.6306,  677.1034,    2.4804,  -13.3573,    0.0000,&
    0.2452,   15.2663,  664.0363,    2.5153,  -13.7561,    0.0000,&
    0.2586,   15.2555,  558.1804,    2.4307,  -13.8080,    0.0000,&
    0.1740,   19.1635,  871.5355,    2.5635,  -16.9369,    0.0000,&
    0.2090,   17.1114,  757.2022,    2.4209,  -13.8086,    0.0000,&
    0.2540,   15.6191,  614.0251,    2.4596,  -13.9194,    0.0000,&
    0.0667,   24.7671, 3183.1985,    2.8455,  -20.4245,    0.0000,&
    0.1653,   17.7768, 1047.0810,    2.7615,  -16.1245,    0.0000,&
    0.2274,   15.5853,  681.8057,    2.4746,  -13.1291,    0.0000,&
    0.2127,   17.7983,  694.1201,    2.4559,  -13.9922,    0.0000,&
    0.2259,   15.9064,  684.2758,    2.4509,  -13.5933,    0.0000,&
    0.2293,   15.4590,  677.4583,    2.5027,  -13.1278,    0.0000,&
    0.2274,   15.5853,  681.8057,    2.4746,  -13.1291,    0.0000,&
    0.1675,   34.1575, 3166.8033,    3.0436,  -32.1388,    0.0000,&
    0.2558,   15.1635,  757.6213,    2.6144,  -13.5418,    0.0000,&
    0.1236,   19.1814, 1601.9044,    2.7590,  -15.6336,    0.0000,&
    0.2604,   15.1185,  598.6912,    2.4377,  -13.6276,    0.0000,&
    0.1930,   19.0387,  796.3681,    2.6703,  -16.7432,    0.0000,&
    0.1400,   13.2600,  493.1865,    2.2897,   -9.7206,    0.0000,&
    0.2730,   22.1013,  873.6143,    2.6562,  -21.5464,    0.0000,&
    0.2209,   16.1180,  687.1810,    2.5135,  -13.3441,    0.0000,&
    0.2406,   19.3421,  675.0634,    2.4117,  -17.2837,    0.0000,&
    0.2011,   20.7274,  854.2212,    2.6685,  -18.8176,    0.0000,&
    0.1020,   19.0722, 1786.0109,    2.6968,  -15.0486,    0.0000,&
    0.2429,   15.6456,  672.1641,    2.4394,  -13.6875,    0.0000,&
    0.2313,   15.8899,  765.6632,    2.4752,  -13.5611,    0.0000,&
    0.2927,   18.6389,  704.5741,    2.9349,  -17.4406,    0.0000,&
    0.2892,   33.9528, 3194.0752,    3.0966,  -33.7544,    0.0000,&
    0.1345,   14.1732,  643.2673,    2.4163,  -11.3910,    0.0000,&
    0.2176,   16.4910,  627.1834,    2.4406,  -13.5176,    0.0000,&
    0.2344,   16.3511,  475.7215,    2.4932,  -14.3649,    0.0000,&
    0.2380,   16.1014,  770.2195,    2.4855,  -13.7447,    0.0000,&
    0.2639,   16.9566,  687.1668,    2.5647,  -15.5746,    0.0000,&
    0.1628,   18.6915, 1236.6294,    2.8375,  -16.9045,    0.0000,&
    0.1941,   16.3585,  812.1834,    2.4401,  -13.0769,    0.0000,&
    0.2061,   15.9805,  855.4111,    2.4812,  -13.2638,    0.0000,&
    0.2074,   16.6526,  670.4411,    2.3680,  -13.2717,    0.0000,&
    0.2335,   16.5899,  643.1196,    2.5251,  -13.8091,    0.0000,&
    0.1801,   40.4797, 7144.7521,    3.2914,  -39.2418,    0.0000,&
    0.1750,   16.1663,  629.9645,    2.5433,  -14.4684,    0.0000,&
    0.2168,   14.7552,  572.7328,    2.4104,  -11.2275,    0.0000,&
    0.2266,   15.6123,  681.1548,    2.4873,  -13.2889,    0.0000,&
    0.2207,   16.5433,  695.4979,    2.4342,  -13.4237,    0.0000,&
    0.1751,   26.8690, 1397.0611,    2.7454,  -23.4472,    0.0000,&
    0.2542,   17.3964,  704.7176,    2.5293,  -16.2353,    0.0000,&
    0.1926,   18.3893,  702.1283,    2.3978,  -14.2432,    0.0000,&
    0.2274,   15.6688,  674.7055,    2.4686,  -13.3213,    0.0000,&
    0.2512,   15.8351,  593.7068,    2.3861,  -13.9120,    0.0000,&
    0.1337,   27.5632, 2322.6027,    2.8275,  -24.4908,    0.0000,&
    0.2160,   16.2282,  713.5225,    2.4703,  -13.3426,    0.0000,&
    0.2311,   18.5889,  759.6157,    2.6432,  -17.3098,    0.0000,&
    0.2348,   15.4985,  672.6298,    2.4413,  -13.5010,    0.0000,&
    0.1715,   18.2175, 1107.6706,    2.6242,  -13.6181,    0.0000,&
    0.1666,   17.1156, 1141.4458,    2.7249,  -15.1379,    0.0000,&
    0.2610,   19.2751,  730.8568,    2.6121,  -18.0220,    0.0000,&
    0.2331,   15.7220,  678.0447,    2.4523,  -13.1372,    0.0000,&
    0.2283,   14.3607,  531.6170,    2.5730,  -12.8014,    0.0000,&
    0.1507,   18.8856, 1147.5717,    2.4881,  -13.4925,    0.0000,&
    0.2320,   15.9280,  497.5950,    2.4784,  -13.2243,    0.0000,&
    0.2332,   16.2661,  711.9251,    2.4291,  -13.6386,    0.0000,&
    0.2032,   17.0936,  827.0512,    2.4810,  -13.7942,    0.0000,&
    0.2339,   15.6710,  673.3325,    2.4549,  -13.3057,    0.0000,&
    0.2601,   15.3745,  502.1801,    2.4801,  -13.8283,    0.0000,&
    0.2341,   16.1424,  704.7703,    2.5082,  -13.2223,    0.0000,&
    0.2221,   16.4109,  675.7022,    2.4094,  -13.4056,    0.0000,&
    0.2528,   15.5196,  693.1694,    2.5086,  -13.5595,    0.0000,&
    0.2301,   17.2881,  780.2085,    2.5899,  -15.9584,    0.0000,&
    0.1612,   20.3271,  801.4255,    2.4766,  -17.4896,    0.0000,&
    0.2219,   16.6959,  665.7473,    2.4166,  -13.4919,    0.0000,&
    0.2490,   15.4133,  520.6750,    2.3508,  -13.5440,    0.0000,&
    0.2238,   15.9202,  675.7761,    2.4460,  -13.3692,    0.0000,&
    0.2422,   15.8169,  680.7375,    2.4446,  -13.6335,    0.0000,&
    0.2300,   15.4991,  659.0534,    2.4002,  -13.0756,    0.0000,&
    0.2286,   16.9948,  691.1089,    2.4700,  -14.0267,    0.0000,&
    0.2388,   15.7691,  679.9312,    2.4518,  -13.6708,    0.0000,&
    0.2092,   20.5862,  643.9759,    2.4708,  -18.4158,    0.0000,&
    0.2406,   19.3421,  675.0634,    2.4117,  -17.2837,    0.0000,&
    0.1542,   17.0930,  773.9340,    2.5352,  -14.8961,    0.0000,&
    0.2395,   16.0252,  732.7911,    2.5242,  -14.0310,    0.0000,&
    0.2345,   16.2107,  706.6905,    2.4525,  -13.6715,    0.0000,&
    0.2205,   15.5146,  681.4560,    2.4785,  -12.7099,    0.0000,&
    0.2124,   17.2085,  908.0841,    2.5774,  -14.3275,    0.0000,&
    0.2645,   15.8323,  786.6830,    2.6474,  -14.5935,    0.0000,&
    0.2103,   16.6139,  776.1693,    2.4880,  -13.5012,    0.0000,&
    0.2352,   15.6121,  672.4533,    2.4406,  -13.5209,    0.0000,&
    0.1806,   34.9405, 3469.0222,    3.0433,  -33.1939,    0.0000,&
    0.2294,   15.8509,  681.3873,    2.4444,  -13.5447,    0.0000,&
    0.2407,   19.7303,  590.2765,    2.4939,  -16.6299,    0.0000,&
    0.2454,   16.2014,  463.3958,    2.3596,  -13.4469,    0.0000,&
    0.2510,   15.1470,  618.7309,    2.5379,  -13.4359,    0.0000,&
    0.2046,   16.8107,  763.7387,    2.4430,  -13.7554,    0.0000,&
    0.2173,   16.5895,  715.5321,    2.5086,  -13.3976,    0.0000,&
    0.2118,   18.2476,  715.3845,    2.3860,  -15.2801,    0.0000,&
    0.2329,   15.9836,  640.9970,    2.5281,  -13.5132,    0.0000,&
    0.2605,   15.5437,  681.0965,    2.5521,  -14.6574,    0.0000,&
    0.1496,   21.0530, 1153.4877,    2.6739,  -19.0440,    0.0000,&
    0.2469,   15.4494,  657.0395,    2.4474,  -13.6476,    0.0000,&
    0.2602,   15.7598,  664.6346,    2.5119,  -14.2762,    0.0000,&
    0.2441,   15.4437,  607.0260,    2.4190,  -13.7670,    0.0000,&
    0.3157,   28.6080,  946.0958,    2.9921,  -28.3380,    0.0000/)
    
    allocate(NParray(1:size(NPparam)))
    if(rep>100) then
        write(*,*) color('ERROR in BSM19.PDF4LHC model. It has only 100 replicas. Central replica is set',c_red)
        NParray=1d0*replicas((0+2)*6+1:(0+2)*6+6)
    else
        NParray=1d0*replicas((rep+2)*6+1:(rep+2)*6+6)
    end if
end subroutine GetReplicaParameters

end module uTMDPDF_model
