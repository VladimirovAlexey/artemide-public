!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Model for unpolarized TMD evolution for BSV19.bFIX [1902.08474]
!
!			corresponds to bb* model, at b =2.5 GeV^{-1}
!			DNP=Dpert(b*)+g bb*
!			zeta=zetaPert(b*)
!
!			Requres two NP parameters (initated by best values values)
!
!				A.Vladimirov (20.02.2019)
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
  2.5000,    0.037484,& !! suggested
  2.5000,    0.037484,&
  2.5000,    0.037484,&!! mean replica
  2.5000,    0.0263,&
  2.5000,    0.0414,&
  2.5000,    0.0310,&
  2.5000,    0.0307,&
  2.5000,    0.0338,&
  2.5000,    0.0431,&
  2.5000,    0.0398,&
  2.5000,    0.0368,&
  2.5000,    0.0310,&
  2.5000,    0.0273,&
  2.5000,    0.0248,&
  2.5000,    0.0402,&
  2.5000,    0.0410,&
  2.5000,    0.0310,&
  2.5000,    0.0359,&
  2.5000,    0.0348,&
  2.5000,    0.0332,&
  2.5000,    0.0300,&
  2.5000,    0.0248,&
  2.5000,    0.0389,&
  2.5000,    0.0355,&
  2.5000,    0.0411,&
  2.5000,    0.0389,&
  2.5000,    0.0341,&
  2.5000,    0.0316,&
  2.5000,    0.0357,&
  2.5000,    0.0437,&
  2.5000,    0.0464,&
  2.5000,    0.0436,&
  2.5000,    0.0427,&
  2.5000,    0.0420,&
  2.5000,    0.0102,&
  2.5000,    0.0311,&
  2.5000,    0.0393,&
  2.5000,    0.0479,&
  2.5000,    0.0281,&
  2.5000,    0.0388,&
  2.5000,    0.0446,&
  2.5000,    0.0361,&
  2.5000,    0.0360,&
  2.5000,    0.0268,&
  2.5000,    0.0387,&
  2.5000,    0.0315,&
  2.5000,    0.0436,&
  2.5000,    0.0396,&
  2.5000,    0.0468,&
  2.5000,    0.0433,&
  2.5000,    0.0369,&
  2.5000,    0.0434,&
  2.5000,    0.0332,&
  2.5000,    0.0475,&
  2.5000,    0.0418,&
  2.5000,    0.0489,&
  2.5000,    0.0280,&
  2.5000,    0.0324,&
  2.5000,    0.0293,&
  2.5000,    0.0442,&
  2.5000,    0.0389,&
  2.5000,    0.0382,&
  2.5000,    0.0492,&
  2.5000,    0.0214,&
  2.5000,    0.0342,&
  2.5000,    0.0506,&
  2.5000,    0.0380,&
  2.5000,    0.0343,&
  2.5000,    0.0383,&
  2.5000,    0.0370,&
  2.5000,    0.0435,&
  2.5000,    0.0463,&
  2.5000,    0.0423,&
  2.5000,    0.0463,&
  2.5000,    0.0375,&
  2.5000,    0.0372,&
  2.5000,    0.0339,&
  2.5000,    0.0421,&
  2.5000,    0.0419,&
  2.5000,    0.0342,&
  2.5000,    0.0333,&
  2.5000,    0.0307,&
  2.5000,    0.0395,&
  2.5000,    0.0382,&
  2.5000,    0.0428,&
  2.5000,    0.0385,&
  2.5000,    0.0454,&
  2.5000,    0.0395,&
  2.5000,    0.0305,&
  2.5000,    0.0456,&
  2.5000,    0.0420,&
  2.5000,    0.0396,&
  2.5000,    0.0423,&
  2.5000,    0.0335,&
  2.5000,    0.0329,&
  2.5000,    0.0389,&
  2.5000,    0.0455,&
  2.5000,    0.0389,&
  2.5000,    0.0369,&
  2.5000,    0.0368,&
  2.5000,    0.0323,&
  2.5000,    0.0511,&
  2.5000,    0.0293/)
  
 if(rep>100) then
   write(*,*) 'ERROR in BSV19_EXP model. It has only 100 replicas. Central replica is set'
   rep=0
 end if
 ReplicaParameters=1d0*replicas((rep+2)*2+1 : (rep+2)*2+2)
 
 end function ReplicaParameters
