!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Model for unpolarized TMD evolution for BSV19.bFIT.noLHC (without LHC data points)
!
!			corresponds to bb* model
!			DNP=Dpert(b*)+g bb*
!			zeta=zetaPert(b*)
!
!			Requres two NP parameters (initated by best values values)
!
!				A.Vladimirov (15.01.2019)
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
 real,parameter,dimension(1:196)::replicas=(/&
  1.21499, 0.0571695,& !! suggested
  0.,0.,                &
  1.21499, 0.0571695,&!! mean replica
  0.9914,    0.0981,&
  1.1190,    0.0950,&
  0.9960,    0.0770,&
  3.3494,    0.0100,&
  1.0808,    0.0263,&
  1.2037,    0.0695,&
  1.0769,    0.0594,&
  0.9907,    0.0812,&
  1.0745,    0.0806,&
  2.8843,    0.0101,&
  0.9930,    0.0842,&
  0.9970,    0.0916,&
  1.1434,    0.0816,&
  1.0634,    0.0103,&
  1.4059,    0.0105,&
  2.2713,    0.0124,&
  0.9936,    0.1242,&
  0.9929,    0.0104,&
  0.9902,    0.0834,&
  0.9966,    0.0103,&
  1.1235,    0.0473,&
  1.0047,    0.0680,&
  1.0047,    0.0963,&
  0.9900,    0.0406,&
  0.9910,    0.0415,&
  0.9939,    0.0765,&
  1.6961,    0.0140,&
  0.9904,    0.0217,&
  0.9903,    0.0405,&
  1.0038,    0.1190,&
  1.9087,    0.0144,&
  0.9902,    0.0134,&
  0.9936,    0.1094,&
  0.9901,    0.0291,&
  1.1082,    0.0678,&
  0.9932,    0.0998,&
  1.0477,    0.0426,&
  1.0869,    0.0893,&
  1.4649,    0.0100,&
  1.2430,    0.0630,&
  1.3138,    0.0524,&
  1.0224,    0.0977,&
  0.9901,    0.0747,&
  0.9902,    0.0115,&
  1.0122,    0.0840,&
  1.0064,    0.1247,&
  1.1275,    0.0640,&
  0.9983,    0.1203,&
  1.0014,    0.0203,&
  3.2556,    0.0119,&
  0.9937,    0.1106,&
  0.9903,    0.0724,&
  1.0564,    0.0100,&
  1.4326,    0.0159,&
  0.9901,    0.0520,&
  1.0105,    0.0105,&
  0.9910,    0.0523,&
  0.9903,    0.0106,&
  1.5972,    0.0548,&
  1.2215,    0.0694,&
  0.9939,    0.0269,&
  1.5418,    0.0288,&
  0.9990,    0.1101,&
  2.1103,    0.0100,&
  1.0139,    0.0768,&
  0.9920,    0.1401,&
  0.9936,    0.1094,&
  1.0305,    0.1286,&
  0.9913,    0.0869,&
  0.9906,    0.0795,&
  1.0036,    0.0867,&
  0.9926,    0.0100,&
  0.9902,    0.0504,&
  2.0447,    0.0130,&
  0.9913,    0.0264,&
  1.0046,    0.1044,&
  0.9905,    0.0142,&
  1.2726,    0.0159,&
  1.0012,    0.0652,&
  1.0447,    0.0621,&
  1.0089,    0.0286,&
  3.4126,    0.0123,&
  0.9902,    0.0222,&
  0.9903,    0.0783,&
  1.3941,    0.0126,&
  1.0343,    0.1233,&
  1.3443,    0.0767,&
  1.7136,    0.0282,&
  1.0303,    0.0955,&
  1.2167,    0.0429,&
  1.1075,    0.0868,&
  0.9975,    0.0655,&
  0.9900,    0.0102,&
  0.9901,    0.1269,&
  0.9900,    0.0254/)
  
  if(rep>95) then
   write(*,*) 'ERROR in BSV19_EXP model. It has only 95 replicas. Central replica is set'
   rep=0
  end if

 ReplicaParameters=1d0*replicas((rep+2)*2+1 : (rep+2)*2+2)
 
 end function ReplicaParameters
