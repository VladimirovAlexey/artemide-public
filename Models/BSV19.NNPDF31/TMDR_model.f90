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
  2.0285,    0.0293 ,&
  2.0247,    0.0284 ,&
  2.0247,    0.0291 ,&
  1.1824,    0.0635 ,&
  2.0369,    0.0276 ,&
  1.9747,    0.0290 ,&
  1.9036,    0.0238 ,&
  2.0391,    0.0297 ,&
  1.9871,    0.0269 ,&
  1.4614,    0.0317 ,&
  1.8833,    0.0268 ,&
  1.4066,    0.0427 ,&
  1.7648,    0.0288 ,&
  2.0176,    0.0283 ,&
  2.0042,    0.0292 ,&
  2.1645,    0.0159 ,&
  2.3748,    0.0203 ,&
  1.8945,    0.0269 ,&
  2.0209,    0.0286 ,&
  1.9772,    0.0246 ,&
  2.1423,    0.0220 ,&
  2.0270,    0.0299 ,&
  1.9922,    0.0262 ,&
  1.9176,    0.0303 ,&
  1.9396,    0.0264 ,&
  1.2512,    0.0582 ,&
  1.9249,    0.0249 ,&
  1.9502,    0.0291 ,&
  1.9660,    0.0266 ,&
  1.8658,    0.0270 ,&
  1.9285,    0.0250 ,&
  1.5354,    0.0394 ,&
  2.0343,    0.0230 ,&
  2.1683,    0.0246 ,&
  2.0070,    0.0248 ,&
  1.9783,    0.0289 ,&
  1.8995,    0.0304 ,&
  1.9941,    0.0185 ,&
  1.8743,    0.0273 ,&
  1.1875,    0.0637 ,&
  1.5608,    0.0385 ,&
  1.3902,    0.0455 ,&
  1.9133,    0.0278 ,&
  1.7713,    0.0309 ,&
  1.9898,    0.0280 ,&
  1.8541,    0.0168 ,&
  1.9946,    0.0205 ,&
  2.0552,    0.0303 ,&
  1.8067,    0.0215 ,&
  1.7589,    0.0260 ,&
  1.6545,    0.0348 ,&
  1.8507,    0.0280 ,&
  1.9812,    0.0287 ,&
  2.1305,    0.0257 ,&
  1.9007,    0.0291 ,&
  1.9758,    0.0277 ,&
  1.7821,    0.0323 ,&
  1.9574,    0.0279 ,&
  2.0274,    0.0284 ,&
  2.0164,    0.0256 ,&
  1.8797,    0.0216 ,&
  2.3051,    0.0225 ,&
  2.0366,    0.0294 ,&
  1.5082,    0.0293 ,&
  1.8564,    0.0270 ,&
  1.3322,    0.0503 ,&
  1.9815,    0.0279 ,&
  1.8680,    0.0220 ,&
  1.9244,    0.0271 ,&
  1.1197,    0.0665 ,&
  2.0035,    0.0182 ,&
  2.0502,    0.0234 ,&
  1.7790,    0.0210 ,&
  1.9292,    0.0276 ,&
  1.3224,    0.0526 ,&
  1.4838,    0.0403 ,&
  1.2465,    0.0613 ,&
  2.0508,    0.0173 ,&
  1.3592,    0.0422 ,&
  1.9822,    0.0271 ,&
  1.6373,    0.0051 ,&
  1.9945,    0.0284 ,&
  1.5422,    0.0308 ,&
  1.8809,    0.0264 ,&
  2.0989,    0.0240 ,&
  2.0610,    0.0282 ,&
  1.9383,    0.0270 ,&
  2.0048,    0.0287 ,&
  2.0006,    0.0208 ,&
  1.6064,    0.0327 ,&
  2.1141,    0.0242 ,&
  2.0222,    0.0296 ,&
  2.0719,    0.0256 ,&
  1.8851,    0.0260 ,&
  1.6112,    0.0391 ,&
  1.2944,    0.0490 ,&
  2.0584,    0.0124 ,&
  2.1209,    0.0281 ,&
  2.1215,    0.0178 ,&
  1.8308,    0.0253 /)
  
 ReplicaParameters=1d0*replicas((rep+2)*2+1:(rep+2)*2+2)
 
 end function ReplicaParameters
