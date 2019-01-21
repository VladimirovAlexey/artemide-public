!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Model for unpolarized TMD evolution for BSV19_EXP
!
!			corresponds to bb* model
!			DNP=Dpert(b*)+g bb*
!			zeta=zetaPert(b*)
!
!			Requres two NP parameters (initated by best values values)
!
!				A.Vladimirov (26.12.2018)
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
  DNP=Dresum(mu,bSTAR,f)+NPparam(2)*b*bSTAR!!!! D*+gK b b*, it smoother turns perturbative to b^2 assimptotic
 
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
 real*8::sqrtBMAX

 sqrtBMAX=SQRT(1+b**2/NPparam(1)**2) 
 zetaNP=zetaMUpert(mu,b/sqrtBMAX,f)
 end function zetaNP
 
 subroutine TMDR_SetReplica(rep)
 integer::rep,i
 real*8,dimension(1:NPlength):: InitVar
   InitVar(1:2)=ReplicaParameters(rep)
    if(NPlength>2) then
     do i=3,NPlength
      InitVar(i)=0d0
     end do
    end if
    
    !!! we also initialize the variables
    call TMDR_setNPparameters(InitVar)
    
 end subroutine TMDR_SetReplica 
 
 !!! this is the table of replica prameters extracted in fit BSV19.
 !!! -2 is suggested for initialization replica
 !!! -1 is the best fit
 !!! 0 is the mean reaplics
 !!! 1 -- 100 replicas
 function ReplicaParameters(rep)
 integer::rep
 real*8::ReplicaParameters(1:2)
 real,parameter,dimension(1:206)::replicas=(/&
  3.5217,    0.0293,&  
  3.5217,    0.0293,&  
  3.3533,    0.0241,&  
  3.1451,    0.0240,&  
  3.2388,    0.0299,&  
  2.4845,    0.0361,&  
  2.6653,    0.0389,&  
  3.9255,    0.0152,&  
  2.9230,    0.0256,&  
  2.8664,    0.0390,&  
  3.2558,    0.0213,&  
  3.4764,    0.0201,&  
  3.5011,    0.0160,&  
  3.1475,    0.0242,&  
  3.4982,    0.0209,&  
  2.7269,    0.0388,&  
  3.3026,    0.0221,&  
  3.1705,    0.0255,&  
  3.4420,    0.0274,&  
  3.1017,    0.0283,&  
  3.4555,    0.0243,&  
  3.7440,    0.0116,&  
  3.6258,    0.0220,&  
  3.7072,    0.0209,&  
  3.9738,    0.0213,&  
  3.3446,    0.0294,&  
  3.4639,    0.0217,&  
  3.1956,    0.0128,&  
  3.2426,    0.0242,&  
  2.7158,    0.0345,&  
  3.5384,    0.0214,&  
  3.2757,    0.0322,&  
  2.9425,    0.0248,&  
  3.3808,    0.0235,&  
  3.5904,    0.0208,&  
  3.4407,    0.0224,&  
  3.7059,    0.0110,&  
  3.4218,    0.0141,&  
  3.1140,    0.0245,&  
  3.3913,    0.0313,&  
  3.4533,    0.0251,&  
  3.4280,    0.0241,&  
  2.8601,    0.0264,&  
  3.5166,    0.0112,&  
  3.4134,    0.0302,&  
  3.3620,    0.0315,&  
  3.4386,    0.0235,&  
  3.4398,    0.0327,&  
  3.9841,    0.0169,&  
  3.6241,    0.0278,&  
  3.1365,    0.0354,&  
  3.9630,    0.0164,&  
  3.3998,    0.0279,&  
  3.9954,    0.0173,&  
  3.7072,    0.0209,&  
  3.2018,    0.0302,&  
  3.5562,    0.0193,&  
  3.0926,    0.0343,&  
  3.5741,    0.0276,&  
  3.9468,    0.0126,&  
  3.1675,    0.0280,&  
  3.5222,    0.0256,&  
  3.3620,    0.0315,&  
  3.3110,    0.0245,&  
  3.2661,    0.0210,&  
  2.4950,    0.0338,&  
  3.3806,    0.0244,&  
  3.4358,    0.0283,&  
  3.9782,    0.0157,&  
  3.2039,    0.0271,&  
  3.3624,    0.0246,&  
  3.2844,    0.0332,&  
  3.3524,    0.0177,&  
  3.8853,    0.0119,&  
  3.5093,    0.0249,&  
  3.4966,    0.0180,&  
  3.3200,    0.0310,&  
  3.3029,    0.0264,&  
  3.2923,    0.0285,&  
  3.4991,    0.0254,&  
  3.3630,    0.0273,&  
  3.4101,    0.0299,&  
  3.4855,    0.0270,&  
  2.5081,    0.0144,&  
  3.4358,    0.0283,&  
  3.1861,    0.0119,&  
  3.4958,    0.0159,&  
  3.0031,    0.0311,&  
  3.3775,    0.0291,&  
  3.1735,    0.0308,&  
  2.9230,    0.0256,&  
  3.5723,    0.0177,&  
  3.3125,    0.0179,&  
  3.2662,    0.0175,&  
  3.4280,    0.0169,&  
  3.4237,    0.0212,&  
  3.2472,    0.0171,&  
  3.4855,    0.0270,&  
  3.9755,    0.0189,&  
  3.5300,    0.0271,&  
  3.1259,    0.0207,&  
  3.1424,    0.0301,&  
  3.5824,    0.0227/)
  
  if(rep>100) then
   write(*,*) 'ERROR in BSV19_EXP model. It has only 100 replicas. Central replica is set'
   rep=0
  end if

 ReplicaParameters=1d0*replicas((rep+2)*2+1 : (rep+2)*2+2)
 
 end function ReplicaParameters
