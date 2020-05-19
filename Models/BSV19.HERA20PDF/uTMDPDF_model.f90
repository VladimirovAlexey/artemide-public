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
    
    write(*,*) color(">>>  The model for uTMDPDF is BSV19.HERA20PDF. Please, cite [1902.08474]   <<<",c_cyan)
    
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
 real*8::ReplicaParameters(1:6)
 real,parameter,dimension(1:618)::replicas=(/&
    0.324112, 13.1774, 356.124, 2.05101, -10.4468, 0.1, & !!!! initialisation
    0.3204, 11.8747, 298.59, 1.8738, -9.0685, 0.,&     !!! best
    0.324112, 13.1774, 356.124, 2.05101, -10.4468, 0.,&!!!mean
    0.3056, 11.9201, 348.057, 2.0577, -9.1136, 0.,&
    0.3306, 12.0876, 305.328, 1.9312, -9.6543, 0.,&
    0.3194, 12.0654, 288.875, 2.0521, -9.9809, 0.,&
    0.3216, 11.9711, 305.019, 1.8944, -8.8547, 0.,&
    0.3150, 23.3603, 727.776, 2.2669, -20.736, 0.,&
    0.3531, 10.4019, 213.041, 2.0995, -8.0426, 0.,&
    0.2417, 19.0693, 847.890, 2.5759, -16.152, 0.,&
    0.2744, 22.1776, 1031.37, 2.4917, -18.372, 0.,&
    0.3252, 12.2880, 350.894, 2.0772, -10.123, 0.,&
    0.3155, 11.8648, 301.441, 1.9065, -8.8508, 0.,&
    0.3083, 12.1988, 324.811, 1.9527, -8.7756, 0.,&
    0.3645, 11.5761, 242.206, 1.9468, -10.131, 0.,&
    0.3129, 12.2474, 319.149, 1.9750, -8.8321, 0.,&
    0.3457, 12.2355, 372.363, 2.1657, -10.065, 0.,&
    0.3071, 12.2374, 313.247, 1.8752, -8.9658, 0.,&
    0.3274, 12.4544, 318.701, 2.1118, -9.3767, 0.,&
    0.3180, 12.1120, 286.838, 2.0132, -9.1793, 0.,&
    0.3124, 11.9002, 301.738, 1.9006, -8.6354, 0.,&
    0.3206, 11.9195, 304.158, 1.9743, -8.9143, 0.,&
    0.3365, 11.3330, 273.340, 1.8955, -9.1042, 0.,&
    0.2952, 12.3076, 257.622, 1.8337, -8.9577, 0.,&
    0.3548, 12.0521, 265.022, 1.7961, -9.6444, 0.,&
    0.3425, 12.0143, 327.884, 2.1524, -10.833, 0.,&
    0.2995, 13.1419, 304.492, 1.9632, -9.1726, 0.,&
    0.3610, 13.7303, 428.132, 2.2071, -11.931, 0.,&
    0.3072, 12.2468, 323.809, 1.9252, -8.7189, 0.,&
    0.3135, 11.8331, 301.026, 1.9681, -8.7936, 0.,&
    0.3209, 11.8477, 296.495, 1.8603, -9.1456, 0.,&
    0.3035, 12.9447, 332.113, 1.8853, -9.0541, 0.,&
    0.3180, 11.8151, 292.106, 1.8465, -9.1709, 0.,&
    0.2652, 14.6970, 512.200, 2.0375, -9.3486, 0.,&
    0.3120, 11.7803, 302.945, 1.9051, -8.7837, 0.,&
    0.2857, 13.5101, 323.152, 1.9209, -9.2277, 0.,&
    0.3201, 11.9523, 306.970, 1.9123, -8.8815, 0.,&
    0.3274, 12.1267, 309.561, 1.9422, -9.8200, 0.,&
    0.2936, 12.6886, 335.223, 1.9078, -8.9825, 0.,&
    0.3347, 11.9918, 286.491, 1.8207, -9.1901, 0.,&
    0.3137, 12.8152, 284.951, 2.0023, -9.3510, 0.,&
    0.3833, 21.0787, 719.319, 2.5212, -20.891, 0.,&
    0.3041, 12.5797, 343.898, 1.9383, -9.1325, 0.,&
    0.3709, 17.3762, 551.378, 2.8875, -17.267, 0.,&
    0.3091, 13.3842, 265.608, 1.8682, -9.2338, 0.,&
    0.3203, 11.9569, 303.344, 1.9069, -8.8811, 0.,&
    0.3621, 14.1185, 437.309, 2.4384, -13.105, 0.,&
    0.4656, 25.4515, 852.465, 3.3920, -26.936, 0.,&
    0.3494, 11.9962, 328.557, 2.1406, -10.294, 0.,&
    0.3275, 12.1692, 272.119, 1.9942, -9.5009, 0.,&
    0.3500, 19.9742, 555.072, 2.7342, -18.882, 0.,&
    0.3156, 11.7804, 305.988, 1.9222, -9.0334, 0.,&
    0.3661, 11.5017, 250.507, 2.2435, -10.784, 0.,&
    0.2944, 12.2956, 330.651, 2.0439, -8.9949, 0.,&
    0.3394, 12.1588, 331.428, 1.9473, -9.0867, 0.,&
    0.2818, 21.8992, 756.356, 2.5887, -19.961, 0.,&
    0.3111, 12.2000, 313.883, 1.9637, -8.8319, 0.,&
    0.3575, 11.8890, 229.243, 2.0393, -10.214, 0.,&
    0.3314, 11.7750, 306.600, 1.7965, -8.1548, 0.,&
    0.3010, 12.9970, 333.086, 1.8854, -9.0549, 0.,&
    0.3505, 11.4214, 176.802, 1.8440, -9.2912, 0.,&
    0.3338, 12.1888, 308.569, 1.9200, -9.6191, 0.,&
    0.3680, 21.3249, 1007.01, 2.8427, -20.766, 0.,&
    0.3145, 12.1297, 285.255, 2.0228, -8.9848, 0.,&
    0.3101, 12.9072, 326.391, 2.0181, -9.5522, 0.,&
    0.3253, 11.9948, 302.272, 1.8695, -9.1089, 0.,&
    0.3272, 12.1039, 302.617, 1.8424, -9.1648, 0.,&
    0.3251, 11.9482, 330.235, 1.8816, -9.0425, 0.,&
    0.3604, 10.5715, 160.163, 1.9152, -9.2112, 0.,&
    0.3084, 12.0389, 331.970, 2.0403, -8.6332, 0.,&
    0.3197, 11.8339, 302.786, 1.8942, -8.9273, 0.,&
    0.3415, 11.2497, 249.517, 1.8935, -9.3240, 0.,&
    0.3343, 12.2863, 327.052, 1.9890, -9.6423, 0.,&
    0.3527, 12.8976, 199.374, 1.8179, -10.689, 0.,&
    0.3199, 12.1185, 319.480, 1.9647, -8.7978, 0.,&
    0.3328, 11.0540, 266.653, 1.8992, -9.0306, 0.,&
    0.3201, 11.9734, 323.233, 1.9534, -8.9580, 0.,&
    0.3197, 11.9095, 300.856, 1.8729, -9.0765, 0.,&
    0.2883, 13.7655, 332.016, 2.0077, -9.2920, 0.,&
    0.3153, 12.2592, 285.970, 1.9568, -8.9298, 0.,&
    0.3535, 11.1629, 266.002, 2.0248, -9.4822, 0.,&
    0.3217, 12.3084, 301.996, 1.9290, -9.0054, 0.,&
    0.2908, 12.3705, 380.138, 2.0715, -8.9685, 0.,&
    0.3676, 19.7997, 707.163, 2.5258, -18.912, 0.,&
    0.3370, 12.1938, 314.991, 2.0312, -10.498, 0.,&
    0.3095, 12.5172, 355.286, 2.1882, -9.1919, 0.,&
    0.3624, 12.9066, 211.347, 2.4972, -12.339, 0.,&
    0.3079, 13.1568, 359.410, 2.1079, -9.2492, 0.,&
    0.3180, 19.3917, 602.714, 2.2707, -16.717, 0.,&
    0.3578, 12.5400, 248.219, 2.1852, -10.779, 0.,&
    0.3372, 11.9493, 274.622, 1.8642, -9.0967, 0.,&
    0.3145, 12.1297, 285.255, 2.0228, -8.9848, 0.,&
    0.3124, 12.0197, 302.710, 1.9795, -8.9914, 0.,&
    0.3217, 11.7600, 321.951, 2.0112, -9.0425, 0.,&
    0.3175, 11.3191, 286.836, 2.0040, -8.3964, 0.,&
    0.3056, 11.8520, 319.992, 1.9922, -8.3936, 0.,&
    0.3037, 12.2829, 299.734, 2.0041, -8.8125, 0.,&
    0.2662, 12.5199, 411.944, 1.8930, -7.1477, 0.,&
    0.3009, 12.2624, 361.820, 2.0522, -8.9050, 0.,&
    0.3388, 12.0498, 245.751, 1.8361, -9.4326, 0.,&
    0.2812, 13.2875, 347.417, 2.0236, -8.6733, 0.,&
    0.3348, 12.1934, 303.271, 1.8588, -9.6003, 0.,&
    0.3183, 12.2928, 312.331, 1.9794, -8.9033, 0./)
    
    allocate(NParray(1:size(NPparam)))
    if(rep>100) then
        write(*,*) color('ERROR in BSM19.HERA20PDF model. It has only 100 replicas. Central replica is set',c_red)
        NParray=1d0*replicas((0+2)*6+1:(0+2)*6+6)
    else
        NParray=1d0*replicas((rep+2)*6+1:(rep+2)*6+6)
    end if
end subroutine GetReplicaParameters

end module uTMDPDF_model
