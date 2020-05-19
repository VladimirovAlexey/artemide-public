!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Model for unpolarized TMD PDF  Vpion19.bFIT  [1907.10356]
!
!			proton uTMDPDF is from BSV19.HERA set   (h=1)
!			pion uTMDPDF is here			(h=2)
!
!			Requres 6 (proton)+3 (pion)=9 NP parameters
!			Uses HERAPDF20_NNLO_VAR and JAM18PionPDFnlo
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
    
    write(*,*) color(">>>  The model for uTMDPDF is Vpion19 & BSV19. Please, cite [1902.08474]&[1907.10356]   <<<",c_cyan)
    
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
  real*8::FNP0

   real*8::bb,w1,w2,w3
   
   if(hadron==1) then
   
    bb=bT**2
    w1=lambdaNP(1)*(1-x)+x*lambdaNP(2)+x*(1-x)*lambdaNP(5)
    w2=lambdaNP(3)*x**lambdaNP(4)+lambdaNP(6)
   
    if(w2<0d0 .or. w1<0d0) then
    FNP0=-1d0
    else
    FNP0=Exp(-w1*bb/sqrt(1+w2*bb))
    end if
    
  else 
      bb=bT**2
      w1=(lambdaNP(7)+(1-x)**2*lambdaNP(8))
      w2=lambdaNP(9)
      if(w2<0d0 .or. w1<0d0) then
      FNP0=-1d0
      else
      FNP0=Exp(-w1*bb/sqrt(1+w2*bb))
      end if
      
  end if
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
 real*8,parameter,dimension(1:6):: protonNP=(/0.3204d0, 11.8747d0, 298.593d0, 1.8738d0, -9.0685d0, 0.0d0/)
 real,parameter,dimension(1:309)::replicas=(/&
    0.173426, 0.482789, 2.15172, &
    0.093, 0.264, 0.377,&
    0.173426, 0.482789, 2.15172, &
    0.249541, 0.634629, 3.09413, &
    0.105834, 0.156929, 0.385113, &
    0.0978039, 0.281598, 0.449822, &
    0.154557, 0.350553, 1.10694, &
    0.261972, 0.676967, 3.37234, &
    0.183507, 0.520334, 1.87893, &
    0.0785418, 0.356815, 0.41225, &
    0.355825, 0.846914, 5.60176, &
    0.223488, 0.553231, 2.8377, &
    0.629039, 1.61353, 16.4577, &
    0.0818166, 0.360383, 0.435636, &
    0.222895, 0.629612, 2.65102, &
    0.0965146, 0.178813, 0.219953, &
    0.491635, 1.15167, 10.0341, &
    0.153072, 0.319388, 1.10683, &
    0.161597, 0.382618, 1.2612, &
    0.128724, 0.373911, 0.736818, &
    0.152192, 0.290414, 0.823574, &
    0.0954244, 0.278245, 0.356441, &
    0.165523, 0.345776, 1.29734, &
    0.176371, 0.421179, 1.6543, &
    0.198816, 0.340405, 1.68137, &
    0.0894031, 0.322207, 0.387982, &
    0.163753, 0.473674, 1.29232, &
    0.0947285, 0.198516, 0.326766, &
    0.0814235, 0.329594, 0.422357, &
    0.149341, 0.366549, 0.914248, &
    0.0942002, 0.266578, 0.368842, &
    0.133111, 0.572628, 1.31634, &
    0.180704, 0.41721, 1.62999, &
    0.065896, 0.316252, 0.250545, &
    0.10734, 0.247779, 0.362931, &
    0.139521, 0.471966, 1.31441, &
    0.366519, 1.25787, 8.21266, &
    0.0790098, 0.241259, 0.230682, &
    0.581215, 2.27234, 21.0271, &
    0.0954821, 0.261137, 0.374515, &
    0.115915, 0.368228, 0.786806, &
    0.273399, 0.749383, 4.03135, &
    0.465171, 1.07553, 9.80427, &
    0.0903598, 0.263619, 0.406335, &
    0.123613, 0.374445, 0.849558, &
    0.285171, 0.418185, 3.34914, &
    0.269755, 0.553625, 3.96405, &
    0.259095, 1.16033, 4.84876, &
    0.0899398, 0.248281, 0.399757, &
    0.259753, 0.814591, 4.63706, &
    0.0947479, 0.272567, 0.365655, &
    0.108101, 0.256952, 0.452232, &
    0.0914599, 0.304369, 0.38939, &
    0.170683, 0.272946, 1.06934, &
    0.118159, 0.279235, 0.604779, &
    0.264408, 0.762043, 3.82065, &
    0.0784105, 0.316828, 0.458274, &
    0.360117, 1.33631, 9.64109, &
    0.105368, 0.225053, 0.322375, &
    0.0987314, 0.303631, 0.477949, &
    0.150731, 0.437147, 1.11623, &
    0.238012, 0.87718, 2.98115, &
    0.278189, 0.492043, 3.65615, &
    0.0804673, 0.2964, 0.289875, &
    0.0837756, 0.328657, 0.428778, &
    0.100518, 0.276298, 0.456033, &
    0.104566, 0.200711, 0.347386, &
    0.132109, 0.380439, 1.01348, &
    0.113121, 0.188703, 0.36785, &
    0.103887, 0.26594, 0.400361, &
    0.0936283, 0.272979, 0.366824, &
    0.112749, 0.393731, 0.670924, &
    0.12597, 0.491501, 1.02126, &
    0.184632, 0.567039, 1.97799, &
    0.0897044, 0.244245, 0.395551, &
    0.101595, 0.265109, 0.38515, &
    0.247302, 0.471764, 2.98563, &
    0.284248, 0.821081, 4.66352, &
    0.18231, 1.03437, 3.07118, &
    0.108571, 0.375484, 0.727352, &
    0.140538, 0.270434, 0.67072, &
    0.233778, 0.496306, 3.07228, &
    0.120892, 0.378347, 0.696918, &
    0.322058, 0.91204, 6.34466, &
    0.134719, 0.352275, 0.759533, &
    0.157389, 0.4007, 1.20728, &
    0.0814492, 0.37148, 0.442985, &
    0.239761, 0.604956, 2.83285, &
    0.104431, 0.216468, 0.423611, &
    0.113135, 0.307468, 0.522409, &
    0.128644, 0.357123, 0.837743, &
    0.136476, 0.292455, 0.815463, &
    0.143915, 0.468419, 1.26521, &
    0.0938552, 0.272222, 0.374274, &
    0.17918, 0.457854, 1.82332, &
    0.0827782, 0.270842, 0.342522, &
    0.167811, 0.298295, 1.05922, &
    0.170454, 0.315802, 1.18806, &
    0.0885638, 0.321581, 0.444846, &
    0.33685, 1.1168, 6.69006, &
    0.131763, 0.302245, 0.888346, &
    0.117674, 0.38926, 0.906957, &
    0.391747, 0.989056, 7.27382/)
  
  allocate(NParray(1:9))
  if(rep>100) then
   write(*,*) color('ERROR in Vpion19 model. It has only 100 replicas. Central replica is set',c_red)
   NParray=(/protonNP(1),protonNP(2),protonNP(3),protonNP(4),protonNP(5),protonNP(6),&
	  1d0*replicas((0+2)*3+1),1d0*replicas((0+2)*3+2),1d0*replicas((0+2)*3+3)/)
  else 
  NParray=(/protonNP(1),protonNP(2),protonNP(3),protonNP(4),protonNP(5),protonNP(6),&
	  1d0*replicas((rep+2)*3+1),1d0*replicas((rep+2)*3+2),1d0*replicas((rep+2)*3+3)/)
  end if
    
end subroutine GetReplicaParameters
  
end module uTMDPDF_model
