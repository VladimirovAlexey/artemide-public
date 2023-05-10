!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!		Part of code that contains hard coefficient functions for Drell-Yan-like x-Sections
!!					is a part of artemide.TMDX_DY
!!
!!			16.06.2019	created A.Vladimirov
!!
!!								A.Vladimirov
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!-----------------------------------------------------------------------------------------------------
!------------------------------------------Drell-Yan -------------------------------------------------
!!! hard coefficeint taken from 1004.3653 up to 3-loop
!!! it is evaluated at mu=Q
function HardCoefficientDY(mu)
    real(dp)::HardCoefficientDY,mu,LQ!=Log[Q^2/mu^2]=-2Log[c1]
    real(dp)::alpha
    
    If(usePiresum) then 
    
    !!! this expression is get by expanding pi-resummedexponent and CV^2 to fixed order
    HardCoefficientDY=1.d0
    if(orderH_global>=1) then
      LQ=-2d0*LOG(c2_global)
      alpha=As(mu*c2_global)
      HardCoefficientDY=HardCoefficientDY+alpha*&
      (-16.946842488404727d0 + 8d0*LQ-2.6666666666666665d0*LQ**2)
    if(orderH_global>=2) then
      HardCoefficientDY=HardCoefficientDY+alpha**2*&
      (-25.413248632430818d0 -208.56946563098907d0*LQ +28.103117631243492d0*LQ**2&
	- 14.518518518518519d0*LQ**3 + 3.5555555555555554d0*LQ**4)
    if(orderH_global>=3) then
      HardCoefficientDY=HardCoefficientDY+alpha**3*&
      (7884.91043450197d0 -3916.5246445256016d0*LQ +58.76188412075794d0*LQ**2&
      +418.58823871161303d0*LQ**3 +13.708854670518122d0*LQ**4 + 10.271604938271604d0*LQ**5& 
      -3.1604938271604937d0*LQ**6)
    end if
    end if
    end if
    
    HardCoefficientDY=HardCoefficientDY*PiResumFactor_q(alpha)
    
    else
    
    HardCoefficientDY=1.d0
    if(orderH_global>=1) then
      LQ=-2d0*LOG(c2_global)
      alpha=As(mu*c2_global)
      HardCoefficientDY=HardCoefficientDY+alpha*&
      (9.372102581166892d0 + 8d0*LQ-2.6666666666666665d0*LQ**2)
    if(orderH_global>=2) then
      HardCoefficientDY=HardCoefficientDY+alpha**2*&
      (359.39087353234015d0 + 1.9820949255839224d0*LQ - 42.08073588761418d0*LQ**2&
	- 14.518518518518519d0*LQ**3 + 3.5555555555555554d0*LQ**4)
    if(orderH_global>=3) then
      HardCoefficientDY=HardCoefficientDY+alpha**3*&
      (8935.66841729192d0 - 2759.2358438992906d0*LQ - 1417.132743244908d0*LQ**2&
      + 36.47614733116575d0*LQ**3 + 107.28732602899498d0*LQ**4 + 10.271604938271604d0*LQ**5& 
      -3.1604938271604937d0*LQ**6)
    if(orderH_global>=4) then 
      HardCoefficientDY=HardCoefficientDY+alpha**4*&
      (135087.2036735284d0 - 150489.22799257038d0*LQ + 64333.03564525264d0*LQ**2 &
      - 8614.870827808947d0*LQ**3 - 1644.002266122397d0*LQ**4 + 403.72511575802685d0*LQ**5 &
      - 57.47461249405413d0*LQ**6 - 3.950617283950617d0*LQ**7 + 1.8436213991769548d0*LQ**8)
    end if
    end if    
    end if
    end if
    
    end if
end function HardCoefficientDY

!!! factor for resummed pi^2 contributions for DY-process
!!! see [0808.3008]
function PiResumFactor_q(alpha)
  real(dp)::PiResumFactor_q
  real(dp)::alpha,aa,ArcT,LogA,UU
  
  ! Nf=5 everywhere
  
  aa=24.08554367752175d0*alpha  ! beta0*pi*alpha
  ArcT=2d0*atan(aa)
  LogA=Log(1+aa**2)
  UU=24d0/529d0/alpha*(aa*ArcT-LogA)
  if(orderH_global>=2) then
    UU=UU+0.05720391222158297d0*(ArcT**2-LogA**2)+0.6063377746307231d0*LogA
  if(orderH_global>=3) then
    UU=UU+alpha/(1+aa**2)*(-3.4962177171202846d0*aa**2 + 3.3966025898820527d0*aa*ArcT &
      + 0.36427776d0*ArcT**2 + 2.9812442965487187d0*LogA &
      - 0.41535829333333335d0*aa**2+LogA - 0.72855552d0*aa*ArcT*LogA - 0.36427776d0*LogA**2)
   end if
   end if
   
   PiResumFactor_q=Exp(UU)
  
  
end function PiResumFactor_q

!-----------------------------------------------------------------------------------------------------
!------------------------------------------Higgs production-------------------------------------------

!!! effective coupling for vertex HFF through the top-quark trinagle
!!! defined to start from 1 (see e.g.[0809.4283] (11)-(12)
!!! mu is scale of coupling
!!! coeff-function is taken at mu=mT (and Nf=5 or 6) and evolved to the scale mu 
function EffCouplingHFF(mu)
  real(dp)::mu,EffCouplingHFF
  real(dp)::alphaT,alphaMU,betaT,betaMU
  
  alphaT=As(mTOP)
  alphaMU=As(c2_global*mu)
  
  !!! we consider only two situations mu<mTOP (Nf=5) and mu >mTOP (Nf=6)
  !!! for mu<mBOTTOM we do nothing (this situation posible never appears)
  if(mu<=mTOP) then
    if(mu<mBOTTOM) call Warning_Raise('no threashold matching for Higgs-DY for mu<mBOTTOM',&
                            messageCounter,messageTrigger,modulename)
  
  
  !!betaT is beta-function at mT  it is normalized to be 1 at LO
  !!betaMU is beta-function at MU it is normalized to be 1 at LO
  betaT=1d0+alphaT*5.0434782608695645d0
  betaMU=1d0+alphaMU*5.0434782608695645d0
  EffCouplingHFF=1d0
  
  if(orderH_global>=1) then
      betaT=betaT+alphaT**2*23.596618357487923d0
      betaMU=betaMU+alphaMU**2*23.596618357487923d0
      EffCouplingHFF=EffCouplingHFF+alphaT*11d0
      
  if(orderH_global>=2) then
      betaT=betaT+alphaT**3*629.4986515814214d0
      betaMU=betaMU+alphaMU**3*629.4986515814214d0
      EffCouplingHFF=EffCouplingHFF+alphaT**2*98.44444444444444d0
      
  if(orderH_global>=3) then
    call Warning_Raise(&
            'no NNNLO implementation of Higgs coefficient function (so far). Continue NNLO',&
            messageCounter,messageTrigger,modulename)
  end if
  end if
  end if
  
  else !Nf=6
  
  !!betaT is beta-function at mT  it is normalized to be 1 at LO
  !!betaMU is beta-function at MU it is normalized to be 1 at LO
  betaT=1d0+alphaT*3.714285714285714d0
  betaMU=1d0+alphaMU*3.714285714285714d0
  EffCouplingHFF=1d0
  
  if(orderH_global>=1) then
      betaT=betaT+alphaT**2*(-4.642857142857142d0)
      betaMU=betaMU+alphaMU**2*(-4.642857142857142d0)
      EffCouplingHFF=EffCouplingHFF+alphaT*11d0
      
  if(orderH_global>=2) then
      betaT=betaT+alphaT**3*353.1833917971027d0
      betaMU=betaMU+alphaMU**3*353.1833917971027d0
      EffCouplingHFF=EffCouplingHFF+alphaT**2*87.27777777777777d0
      
  if(orderH_global>=3) then
      call Warning_Raise(&
            'no NNNLO implementation of Higgs coefficient function (so far). Continue NNLO',&
            messageCounter,messageTrigger,modulename)
  end if
  end if
  end if
  
  end if
  
  !evolved coupling is 
  
  EffCouplingHFF=betaMU/betaT*EffCouplingHFF

end function EffCouplingHFF


!!! hard coefficeint taken from 1004.3653 up to 3-loop
!!! it is evaluated at mu=Q
function HardCoefficientHIGGS(mu)
    real(dp)::HardCoefficientHIGGS,mu,LQ!=Log[Q^2/mu^2]=-2Log[c1]
    real(dp)::alpha
    
   !Nf=5 here!
    
    If(usePiresum) then 
    
    HardCoefficientHIGGS=1.d0
    if(orderH_global>=1) then
      LQ=-2d0*LOG(c2_global)
      alpha=As(mu*c2_global)
      HardCoefficientHIGGS=HardCoefficientHIGGS+alpha*&
      (9.869604401089358d0 -6d0*LQ**2)
    if(orderH_global>=2) then
      HardCoefficientHIGGS=HardCoefficientHIGGS+alpha**2*&
      (-23.720599432600856d0 -56.83020488600443d0*LQ -100.66666666666666d0*LQ**2&
	+ 15.333333333333334d0*LQ**3 + 18d0*LQ**4)
    if(orderH_global>=3) then
        call Warning_Raise(&
            'no NNNLO implementation of Higgs coefficient function (so far). Continue NNLO',&
            messageCounter,messageTrigger,modulename)
    end if
    end if
    end if
    HardCoefficientHIGGS=HardCoefficientHIGGS*PiResumFactor_g(alpha)
    else 
    
    HardCoefficientHIGGS=1.d0
    
    if(orderH_global>=1) then
      LQ=-2d0*LOG(c2_global)
      alpha=As(mu*c2_global)
      HardCoefficientHIGGS=HardCoefficientHIGGS+alpha*&
      (69.0872308076255d0 -6d0*LQ**2)
    if(orderH_global>=2) then
      HardCoefficientHIGGS=HardCoefficientHIGGS+alpha**2*&
      (2723.1832155557718d0 -510.83200733611494d0*LQ - 455.9724251058836d0*LQ**2&
	+ 15.333333333333334d0*LQ**3 + 18d0*LQ**4)
    if(orderH_global>=3) then
     call Warning_Raise(&
            'no NNNLO implementation of Higgs coefficient function (so far). Continue NNLO',&
            messageCounter,messageTrigger,modulename)
    end if
    end if
    end if
    
    end if
end function HardCoefficientHIGGS



!!! factor for resummed pi^2 contributions for HIGGS-process
!!! see [0808.3008]
function PiResumFactor_g(alpha)
  real(dp)::PiResumFactor_g
  real(dp)::alpha,aa,ArcT,LogA,UU
  
  ! Nf=5 everywhere
  
  aa=24.08554367752175d0*alpha  ! beta0*pi*alpha
  ArcT=2d0*atan(aa)
  LogA=Log(1+aa**2)
  UU=54d0/529d0/alpha*(aa*ArcT-LogA)
  if(orderH_global>=2) then
    UU=UU+0.12870880249856168d0*(ArcT**2-LogA**2)+0.19034694944086603d0*LogA
  if(orderH_global>=3) then
    UU=UU+alpha/(1+aa**2)*(1.0895717932098101d0*aa**2 + 0.989555827234617d0*aa*ArcT &
    + 0.81962496d0*ArcT**2 + 0.05499966723461647d0*LogA - 0.93455616d0*aa**2*LogA &
    - 1.63924992d0*aa*ArcT*LogA - 0.81962496*LogA**2)
   end if
   end if
   
   PiResumFactor_g=Exp(UU)
  
  
end function PiResumFactor_g
