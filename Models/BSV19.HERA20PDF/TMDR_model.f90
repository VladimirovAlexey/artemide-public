!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Model for unpolarized TMD evolution for BSV19
!
!			corresponds to bb* model
!			DNP=Dpert(b*)+g bb*
!			zeta=zetaPert(b) exp[-b2/BB]+zetaSL(b)(1-exp(-b2/BB)
!
!			Requres two NP parameters (initated by best values)
!
!				A.Vladimirov (11.07.2019)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! USER DEFINED FUNCTIONS   !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
 
  !!!!!! Write nessecery model intitialization.
  subroutine ModelInitialization()  
    real*8,dimension(1:NPlength):: InitVar
    integer::i
    if(NPlength<2) then
      write(*,*) 'arTeMiDe.TMDR-model: Number NP parameters for TMDR has less then 2'
      write(*,*) 'Evaluation STOP'
      stop
    end if
    !!!! hard core set of evolution parmeters
    !!!! to secure the model
    EvolutionType=3
    
    InitVar(1:2)=ReplicaParameters(-2)
    if(NPlength>2) then
     do i=3,NPlength
      InitVar(i)=0d0
     end do
    end if
    
    !!! we also initialize the variables
    call TMDR_setNPparameters(InitVar)
    
  end subroutine ModelInitialization 
  !!! This is the rapidity anomalous dimension non-pertrubative model
  !!! In your evaluation take care that the saddle point is inside the pertrubative regeme
  !!! Use function Dpert(mu,b,f) for D pertrubative, use Dresum for D resum
  !!! use non-pertrubative parameters NPparam(1...)
 function DNP(mu,b,f)
 real*8::DNP,mu,b
 integer::f
 real*8::bSTAR
 
  bSTAR=b/SQRT(1+b**2/NPparam(1)**2)
  DNP=Dresum(mu,bSTAR,1)+NPparam(2)*b*bSTAR!!!! D*+gK b b*, it smoother turns perturbative to b^2 assimptotic
  
 end function DNP
  
 !! This is the non-pertrubative shape of zeta_mu line.
 !! It MUST follow the equipotential line in pertrubative regime (at small-b), at the level pf PT accuracy.
 !! Otherwice, your evolution is completely broken.
 !! DO NOT modify it if you do not understand what does it mean!
 !!
 !!! Use function zetaMUpert(mu,b,f) for zetamu pertrubative, use zetaMUresum for zetaMu resumed
 !!! use non-pertrubative parameters NPparam(1...)
 !!
 !! Typical form of it is just zetaMUpert(mu,b,f), if b* is used then zetaMUpert(mu,b^*,f)
 !! The large-b deviation from the "true" line is the part of NP model :)
 function zetaNP(mu,b,f)
 real*8::zetaNP,mu,b
 integer::f
 real*8::zz
 
  zz=Exp(-b**2/NPparam(1)**2)
  zetaNP=zetaMUpert(mu,b,f)*zz+zetaSL(mu,b,f)*(1d0-zz)
 end function zetaNP
 
 !!! this is the table of replica prameters extracted in fit BSV19.
 !!! -2 is suggested for initialization replica
 !!! -1 is the best fit
 !!! 0 is the mean reaplics
 !!! 1 -- 100 replicas
 function ReplicaParameters(rep)
 integer::rep
 real*8::ReplicaParameters(1:2)
 real,parameter,dimension(1:206)::replicas=(/ &
  2.29477, 0.022191,& !!! mean
  2.2824, 0.02500,&  !!! best
  2.29477, 0.022191,&!!! mean
  2.1339, 0.0250,&
  2.3001, 0.0250,&
  2.3812, 0.0211,&
  2.2791, 0.0238,&
  1.6002, 0.0342,&
  2.1429, 0.0156,&
  1.0333, 0.0755,&
  2.0739, 0.0195,&
  2.3324, 0.0250,&
  2.2853, 0.0246,&
  2.2409, 0.0229,&
  2.5955, 0.0203,&
  2.1417, 0.0231,&
  3.2238, 0.0120,&
  2.2678, 0.0245,&
  2.1221, 0.0180,&
  2.2896, 0.0188,&
  2.2922, 0.0225,&
  2.4317, 0.0181,&
  2.3500, 0.0255,&
  1.9926, 0.0266,&
  3.0130, 0.0173,&
  2.3229, 0.0270,&
  2.4220, 0.0134,&
  4.1165, 0.0073,&
  2.2213, 0.0230,&
  2.2684, 0.0227,&
  2.3112, 0.0253,&
  2.2693, 0.0209,&
  2.3751, 0.0255,&
  2.1757, 0.0195,&
  2.2865, 0.0242,&
  2.1061, 0.0193,&
  2.2686, 0.0235,&
  2.3453, 0.0254,&
  2.3224, 0.0220,&
  2.2926, 0.0248,&
  2.2996, 0.0148,&
  3.2594, 0.0095,&
  2.2636, 0.0232,&
  1.3063, 0.0409,&
  2.0656, 0.0193,&
  2.2617, 0.0236,&
  2.5138, 0.0165,&
  1.5081, 0.0080,&
  2.3369, 0.0230,&
  2.2629, 0.0194,&
  1.1405, 0.0430,&
  2.3161, 0.0240,&
  2.2390, 0.0230,&
  2.2160, 0.0214,&
  2.7022, 0.0173,&
  1.0320, 0.0686,&
  2.2233, 0.0224,&
  3.1809, 0.0096,&
  3.3448, 0.0130,&
  2.2642, 0.0228,&
  2.7598, 0.0112,&
  2.2723, 0.0237,&
  2.6383, 0.0104,&
  2.1558, 0.0199,&
  2.2721, 0.0185,&
  2.2465, 0.0248,&
  2.3097, 0.0252,&
  2.3186, 0.0258,&
  2.1205, 0.0237,&
  2.1330, 0.0226,&
  2.2916, 0.0248,&
  2.2871, 0.0254,&
  2.2249, 0.0241,&
  2.8602, 0.0103,&
  2.2338, 0.0214,&
  2.3190, 0.0262,&
  2.3613, 0.0221,&
  2.2992, 0.0245,&
  2.0467, 0.0159,&
  2.2911, 0.0182,&
  2.2855, 0.0251,&
  2.3836, 0.0178,&
  2.3031, 0.0217,&
  3.2501, 0.0080,&
  2.3455, 0.0256,&
  2.1001, 0.0195,&
  1.2865, 0.0391,&
  2.3812, 0.0127,&
  2.1293, 0.0209,&
  2.2443, 0.0171,&
  2.1914, 0.0248,&
  2.1558, 0.0199,&
  2.3263, 0.0205,&
  2.5408, 0.0187,&
  2.6147, 0.0159,&
  2.1962, 0.0219,&
  2.1313, 0.0214,&
  2.1259, 0.0235,&
  2.3300, 0.0201,&
  2.2846, 0.0224,&
  2.3368, 0.0141,&
  2.5143, 0.0228,&
  2.2143, 0.0209/)
  
 ReplicaParameters=1d0*replicas((rep+2)*2+1:(rep+2)*2+2)
 
 end function ReplicaParameters
