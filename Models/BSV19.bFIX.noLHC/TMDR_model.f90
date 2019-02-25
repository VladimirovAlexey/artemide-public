!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Model for unpolarized TMD evolution for BSV19.bFIX.noLHC [1902.08474]
!
!			corresponds to bb* model, at fixed bMAX=2.5 GeV^{-1}
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
  2.5000,    0.011883,& !! suggested
  2.5000,    0.011883,&
  2.5000,    0.011883,&!! mean replica
  2.5000,    0.0109,&
  2.5000,    0.0100,&
  2.5000,    0.0127,&
  2.5000,    0.0100,&
  2.5000,    0.0100,&
  2.5000,    0.0107,&
  2.5000,    0.0301,&
  2.5000,    0.0103,&
  2.5000,    0.0100,&
  2.5000,    0.0120,&
  2.5000,    0.0105,&
  2.5000,    0.0101,&
  2.5000,    0.0123,&
  2.5000,    0.0100,&
  2.5000,    0.0100,&
  2.5000,    0.0101,&
  2.5000,    0.0110,&
  2.5000,    0.0218,&
  2.5000,    0.0102,&
  2.5000,    0.0102,&
  2.5000,    0.0102,&
  2.5000,    0.0143,&
  2.5000,    0.0100,&
  2.5000,    0.0108,&
  2.5000,    0.0114,&
  2.5000,    0.0101,&
  2.5000,    0.0105,&
  2.5000,    0.0106,&
  2.5000,    0.0100,&
  2.5000,    0.0135,&
  2.5000,    0.0102,&
  2.5000,    0.0100,&
  2.5000,    0.0100,&
  2.5000,    0.0100,&
  2.5000,    0.0107,&
  2.5000,    0.0103,&
  2.5000,    0.0103,&
  2.5000,    0.0109,&
  2.5000,    0.0100,&
  2.5000,    0.0323,&
  2.5000,    0.0100,&
  2.5000,    0.0101,&
  2.5000,    0.0117,&
  2.5000,    0.0115,&
  2.5000,    0.0100,&
  2.5000,    0.0100,&
  2.5000,    0.0100,&
  2.5000,    0.0102,&
  2.5000,    0.0100,&
  2.5000,    0.0101,&
  2.5000,    0.0115,&
  2.5000,    0.0102,&
  2.5000,    0.0100,&
  2.5000,    0.0102,&
  2.5000,    0.0100,&
  2.5000,    0.0103,&
  2.5000,    0.0178,&
  2.5000,    0.0117,&
  2.5000,    0.0102,&
  2.5000,    0.0101,&
  2.5000,    0.0101,&
  2.5000,    0.0103,&
  2.5000,    0.0107,&
  2.5000,    0.0101,&
  2.5000,    0.0201,&
  2.5000,    0.0149,&
  2.5000,    0.0101,&
  2.5000,    0.0197,&
  2.5000,    0.0113,&
  2.5000,    0.0128,&
  2.5000,    0.0100,&
  2.5000,    0.0310,&
  2.5000,    0.0100,&
  2.5000,    0.0112,&
  2.5000,    0.0101,&
  2.5000,    0.0100,&
  2.5000,    0.0110,&
  2.5000,    0.0101,&
  2.5000,    0.0100,&
  2.5000,    0.0177,&
  2.5000,    0.0100,&
  2.5000,    0.0103,&
  2.5000,    0.0109,&
  2.5000,    0.0102,&
  2.5000,    0.0103,&
  2.5000,    0.0178,&
  2.5000,    0.0174,&
  2.5000,    0.0102,&
  2.5000,    0.0170,&
  2.5000,    0.0101,&
  2.5000,    0.0130,&
  2.5000,    0.0107,&
  2.5000,    0.0118,&
  2.5000,    0.0100,&
  2.5000,    0.0101,&
  2.5000,    0.0100,&
  2.5000,    0.0100,&
  2.5000,    0.0100,&
  2.5000,    0.0127,&
  2.5000,    0.0100/)
  
 if(rep>100) then
   write(*,*) 'ERROR in BSV19_EXP model. It has only 100 replicas. Central replica is set'
   rep=0
 end if
 ReplicaParameters=1d0*replicas((rep+2)*2+1 : (rep+2)*2+2)
 
 end function ReplicaParameters
