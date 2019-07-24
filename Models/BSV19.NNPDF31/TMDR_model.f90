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
  1.8604,   0.02955,& !!! mean
  2.0340,   0.0299,&  !!! best
  1.8604,   0.02955,&!!! mean
  1.6981,   2.0285 ,&
  1.8210,   2.0247 ,&
  1.7391,   2.0247 ,&
  1.8722,   1.1824 ,&
  1.8708,   2.0369 ,&
  1.7661,   1.9747 ,&
  1.8510,   1.9036 ,&
  1.6216,   2.0391 ,&
  1.9575,   1.9871 ,&
  2.1120,   1.4614 ,&
  1.9420,   1.8833 ,&
  1.9390,   1.4066 ,&
  2.3642,   1.7648 ,&
  1.9985,   2.0176 ,&
  1.6984,   2.0042 ,&
  1.9033,   2.1645 ,&
  1.8626,   2.3748 ,&
  2.0927,   1.8945 ,&
  1.7533,   2.0209 ,&
  1.7597,   1.9772 ,&
  1.9074,   2.1423 ,&
  1.7462,   2.0270 ,&
  1.7080,   1.9922 ,&
  1.7922,   1.9176 ,&
  2.1040,   1.9396 ,&
  1.7676,   1.2512 ,&
  2.0155,   1.9249 ,&
  1.9275,   1.9502 ,&
  1.8443,   1.9660 ,&
  1.8832,   1.8658 ,&
  1.9890,   1.9285 ,&
  2.1006,   1.5354 ,&
  2.0470,   2.0343 ,&
  1.9246,   2.1683 ,&
  2.0512,   2.0070 ,&
  1.8118,   1.9783 ,&
  1.8601,   1.8995 ,&
  2.0893,   1.9941 ,&
  2.1061,   1.8743 ,&
  2.2992,   1.1875 ,&
  1.7746,   1.5608 ,&
  1.6241,   1.3902 ,&
  1.9787,   1.9133 ,&
  1.8280,   1.7713 ,&
  1.8525,   1.9898 ,&
  2.1432,   1.8541 ,&
  2.0290,   1.9946 ,&
  1.6983,   2.0552 ,&
  2.1532,   1.8067 ,&
  1.8415,   1.7589 ,&
  2.2617,   1.6545 ,&
  1.9678,   1.8507 ,&
  1.8145,   1.9812 ,&
  1.6644,   2.1305 ,&
  1.6188,   1.9007 ,&
  1.9156,   1.9758 ,&
  1.8882,   1.7821 ,&
  1.6320,   1.9574 ,&
  1.7781,   2.0274 ,&
  1.9149,   2.0164 ,&
  2.0500,   1.8797 ,&
  1.8957,   2.3051 ,&
  1.7193,   2.0366 ,&
  1.9875,   1.5082 ,&
  1.8568,   1.8564 ,&
  1.6777,   1.3322 ,&
  1.8385,   1.9815 ,&
  1.8384,   1.8680 ,&
  1.9586,   1.9244 ,&
  2.0496,   1.1197 ,&
  1.8393,   2.0035 ,&
  1.8190,   2.0502 ,&
  1.9648,   1.7790 ,&
  1.8311,   1.9292 ,&
  1.8126,   1.3224 ,&
  1.9432,   1.4838 ,&
  1.7008,   1.2465 ,&
  1.7807,   2.0508 ,&
  2.0605,   1.3592 ,&
  1.8394,   1.9822 ,&
  2.8054,   1.6373 ,&
  1.7407,   1.9945 ,&
  1.9846,   1.5422 ,&
  1.8309,   1.8809 ,&
  1.9612,   2.0989 ,&
  2.1561,   2.0610 ,&
  1.8885,   1.9383 ,&
  1.9721,   2.0048 ,&
  1.8482,   2.0006 ,&
  2.0745,   1.6064 ,&
  2.0745,   2.1141 ,&
  1.8822,   2.0222 ,&
  1.7990,   2.0719 ,&
  1.9768,   1.8851 ,&
  1.7682,   1.6112 ,&
  1.6623,   1.2944 ,&
  1.8821,   2.0584 ,&
  1.8768,   2.1209 ,&
  1.9822,   2.1215 ,&
  1.8970,   1.8308 /)
  
 ReplicaParameters=1d0*replicas((rep+2)*2+1:(rep+2)*2+2)
 
 end function ReplicaParameters
