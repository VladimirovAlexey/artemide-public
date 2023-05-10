!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 1.4
!
!	This file contains the part of the code, which is common for all TMD-evaluation modules that operates at twist-2
!	It shares the variables with the module where it is inlucded (as a text)
!	Such idiotic structure is needed since, FORTRAN does not allow inheritance.
!
!	Be AWARE of possible clash of variable names.
!
!	This part is devoted to the calculation of Mellin convolution
!	
!	v.2.00 Added b* AV (27.03.2019)
!
!				A.Vladimirov (08.10.2018)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! Supplimentary functions !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------------------------------------------------------------------------
!Used global variables:
! outputlevel, moduleName
! lambdaNP, FNP, bSTAR, mu_OPE
! QCDinput functions
! aTMDe_Numerics, IO_functions functions
! + variables defined in Twist2Convolution-VAR.f90

!!Evaluate Log[mu^2 bT^2/4/Exp(-2Gamma_E)]
!! the b here is b*, this funciton is used only in Coefficeint functions
function LogMuB(mu,bT)
    real(dp)::LogMuB
    real(dp),intent(in)::bT,mu
    LogMuB=2d0*Log(bSTAR(bT,lambdaNP)*mu*C0_inv_const)
end function LogMuB  


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Convolutions!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
!- This is the TMD function evaluated for all quarks simultaniously (-5..5) at x,bT,mu
!--    f =  -5, -4, -3,  -2,  -1,0,1,2,3, 4 ,5
!--      = bbar ,cbar sbar,ubar,dbar,g,d,u,s, c ,b
!!---Gluon contribution is undefined
!- The order is accumulative pertrubative order of coefficient =0,1,2 (LO,NLO,NNLO)
!---------------------------------------------------------------------
!!!---
!!!	It calculates 1/x int_x^1 dz C[z] f[x/z],
!!!		where f is given by xf,  C is given by coeff. function.
!!!---
!---------------------------------------------------------------------
function Common_lowScale5(x,bT,hadron)
  real(dp),dimension(-5:5)::Common_lowScale5
  real(dp),intent(in) :: x, bT
  integer,intent(in)::hadron
  
  real(dp) :: alpha,alphaAt1
  real(dp) :: Lmu,Nf,LmuAt1, NfAt1
  real(dp) :: remnant
  
  real(dp),dimension(-5:5) :: deltaPart
  real(dp),dimension(-5:5) :: convolutionPart
  integer :: j
  
   xCurrent=x
   !! for extrimely small-values of b we freeze its value at b=10^{-6}.
   if(bT>1d-6) then 
    bTcurrent=bT
   else
    bTcurrent=1d-6
   end if
   
   muCurrent=mu_OPE(x,bTcurrent,c4_global)
   alpha=As(muCurrent)
   Nf=real(activeNf(muCurrent),dp)
   Lmu=LogMuB(muCurrent,bTcurrent)
   
  !! in the x-dependent mu we should additionally calculate values at x=1
  if(IsMuXdependent) then
   !!! we first calculate at z=1
    muAt1=mu_OPE(1d0,bTcurrent,c4_global)
    LmuAt1=LogMuB(muAt1,bTcurrent)
    alphaAt1=As(muAt1)
    NfAt1=real(activeNf(muAt1),dp)
  end if
  
  !! boundary value of FNP*PDF
  if(IsMuXdependent) then
   Fcurrent=FNP(xCurrent,1d0,bTcurrent,hadron,lambdaNP)
   PDFcurrent=xf(x,muAt1,hadron)!!!!!!!! This is important, since y enters definition of d via mu(y)   
   FPDFcurrent=Fcurrent*PDFcurrent
  else
   Fcurrent=FNP(xCurrent,1d0,bTcurrent,hadron,lambdaNP)
   PDFcurrent=xf(x,muCurrent,hadron)
   FPDFcurrent=Fcurrent*PDFcurrent  
  end if  
  
 !------------DELTA PART-------------------
 !Leading order is always here!! 
  if(IsMuXdependent) then  
    deltaPart=FPDFcurrent*C_q_q_delta(alphaAt1,NfAt1,LmuAt1)  
  else
    deltaPart=FPDFcurrent*C_q_q_delta(alpha,Nf,Lmu)
  end if
  
   !!!!evaluate coefficients
  if(order_global>=1) then
    call Set_CoeffSing1_q_q(alpha,Nf,Lmu)
    call Set_Coeff_q_q(alpha,Nf,Lmu)
    call Set_Coeff_q_qb(alpha,Nf,Lmu)
    call Set_Coeff_q_qp(alpha,Nf,Lmu)
    call Set_Coeff_q_g(alpha,Nf,Lmu) 
  end if


    
  if(order_global>=1) then
   !!!! evaluating Mellin convolution
    counter=1 !=1 since there was singe call for FPDFcurrent
    !!! The crude estimation of the integral is its tree-value
    !!! It is needed to weight the adaptive integration
    integralWeight=ABS(FPDFcurrent)
    do j=-5,5!!!this is needed since at low energies some of functions =0.
      if(ABS(2d0*integralWeight(j))<tolerance) integralWeight(j)=tolerance
    end do

    !!!! the integral over 1/(1-x)_+ has a remnant part ~delta(1-x)\int_0^x
    !!!! these terms appears order-by-order,
    !!!! since Log[1-x] can be large I compute it order-by-order to avoid dropping of precision (large number x 0d0)
    remnant=CoeffSing1_q_q(1)*LOG(1d0-xCurrent)  !!! remnant of 1/(1-x)_+
    if(order_global>=2) remnant=remnant+CoeffSing1_q_q(2)*LOG(1d0-xCurrent)**2/2d0  !!! remnant of log[1-x]/(1-x)_+
    if(order_global>=3) remnant=remnant+CoeffSing1_q_q(3)*LOG(1d0-xCurrent)**3/3d0  !!! remnant of Log[1-x]^2/(1-x)_+

    convolutionPart=MellinConvolutionVectorPart5(xCurrent,1d0,hadron)&
      +remnant*FPDFcurrent
      !write(*,*) 'counter GK=',counter
  else
    convolutionPart=0d0
  end if
  
  !write(*,*) 'gluonMIXTUREPart =', gluonMIXTURE/x  
  !write(*,*) 'TDhat', (deltaPart+singularPart+regularPart)/x
  Common_lowScale5=(deltaPart+convolutionPart)*(1d0/x)
  
 end function Common_lowScale5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !---------------------------------------------------------------------
!- This is the TMD function evaluated for all quarks simultaniously (-5..5) at x,bT,mu the GLUON INCLUDED
!--    f =  -5, -4, -3,  -2,  -1,0,1,2,3, 4 ,5
!--      = bbar ,cbar sbar,ubar,dbar,g,d,u,s, c ,b
!- The order is accumulative pertrubative order of coefficient =0,1,2 (LO,NLO,NNLO)
!---------------------------------------------------------------------
!!!---
!!!	It calculates 1/x int_x^1 dz C[z] f[x/z],
!!!		where f is given by xf,  C is given by coeff. function.
!!!---
!---------------------------------------------------------------------
function Common_lowScale50(x,bT,hadron)
  real(dp),dimension(-5:5)::Common_lowScale50
  real(dp),intent(in) :: x, bT
  integer,intent(in)::hadron
  
  real(dp) :: alpha,alphaAt1
  real(dp) :: Lmu,Nf,LmuAt1, NfAt1,dummy1,dummy2
  
  real(dp),dimension(-5:5) :: deltaPart
  real(dp),dimension(-5:5) :: convolutionPart
  
  integer :: j
  
   xCurrent=x
   !! for extrimely small-values of b we freeze its value at b=10^{-6}.
   if(bT>1d-6) then 
    bTcurrent=bT
   else
    bTcurrent=1d-6
   end if
   
   muCurrent=mu_OPE(x,bTcurrent,c4_global)
   alpha=As(muCurrent)
   Nf=real(activeNf(muCurrent),dp)
   Lmu=LogMuB(muCurrent,bTcurrent)
   
   
  !! in the x-dependent mu we should additionally calculate values at x=1
  if(IsMuXdependent) then
   !!! we first calculate at z=1
    muAt1=mu_OPE(1d0,bTcurrent,c4_global)
    LmuAt1=LogMuB(muAt1,bTcurrent)
    alphaAt1=As(muAt1)
    NfAt1=real(activeNf(muAt1),dp)
  end if
  
  !! boundary value of FNP*PDF
  if(IsMuXdependent) then
   Fcurrent=FNP(xCurrent,1d0,bTcurrent,hadron,lambdaNP)
   PDFcurrent=xf(x,muAt1,hadron)!!!!!!!! This is important, since y enters definition of d via mu(y)   
   FPDFcurrent=Fcurrent*PDFcurrent
  else
   Fcurrent=FNP(xCurrent,1d0,bTcurrent,hadron,lambdaNP)
   PDFcurrent=xf(x,muCurrent,hadron)
   FPDFcurrent=Fcurrent*PDFcurrent  
  end if

 !------------DELTA PART-------------------
 !Leading order is always here!! 
  if(IsMuXdependent) then  
    deltaPart=FPDFcurrent*C_q_q_delta(alphaAt1,NfAt1,LmuAt1)  
    deltaPart(0)=FPDFcurrent(0)*C_g_g_delta(alphaAt1,NfAt1,LmuAt1)
  else
    deltaPart=FPDFcurrent*C_q_q_delta(alpha,Nf,Lmu)
    deltaPart(0)=FPDFcurrent(0)*C_g_g_delta(alpha,Nf,Lmu)
  end if

   !!!!evaluate coefficients
  if(order_global>=1) then
    call Set_CoeffSing1_q_q(alpha,Nf,Lmu)
    call Set_CoeffSing1_g_g(alpha,Nf,Lmu)
    call Set_Coeff_q_q(alpha,Nf,Lmu)
    call Set_Coeff_q_qb(alpha,Nf,Lmu)
    call Set_Coeff_q_qp(alpha,Nf,Lmu)
    call Set_Coeff_q_g(alpha,Nf,Lmu)  
    call Set_Coeff_g_g(alpha,Nf,Lmu)
    call Set_Coeff_g_q(alpha,Nf,Lmu)
  end if

    
  if(order_global>=1) then
   !!!! evaluating Mellin convolution
    counter=1 !=1 since there was call in FPDFcurrent
    !!! The crude estimation of the integral is its tree-value
    !!! It is needed to weight the adaptive integration
    integralWeight=ABS(FPDFcurrent)
    do j=-5,5!!!this is needed since at low energies some of function =0.
      if(ABS(2d0*integralWeight(j))<tolerance) integralWeight(j)=tolerance
    end do

    !!!! the integral over 1/(1-x)_+ has a remnant part ~delta(1-x)\int_0^x
    !!!! these terms appears order-by-order,
    !!!! since Log[1-x] can be large I compute it order-by-order to avoid dropping of precision (large number x 0d0)

    !!! quark case
    dummy1=CoeffSing1_q_q(1)*LOG(1d0-xCurrent)  !!! remnant of 1/(1-x)_+
    if(order_global>=2) dummy1=dummy1+CoeffSing1_q_q(2)*LOG(1d0-xCurrent)**2/2d0  !!! remnant of log[1-x]/(1-x)_+
    if(order_global>=3) dummy1=dummy1+CoeffSing1_q_q(3)*LOG(1d0-xCurrent)**3/3d0  !!! remnant of Log[1-x]^2/(1-x)_+

    !!! gluon case
    dummy2=CoeffSing1_g_g(1)*LOG(1d0-xCurrent)  !!! remnant of 1/(1-x)_+
    if(order_global>=2) dummy2=dummy2+CoeffSing1_g_g(2)*LOG(1d0-xCurrent)**2/2d0  !!! remnant of log[1-x]/(1-x)_+
    if(order_global>=3) dummy2=dummy2+CoeffSing1_g_g(3)*LOG(1d0-xCurrent)**3/3d0  !!! remnant of Log[1-x]^2/(1-x)_+

    convolutionPart=MellinConvolutionVectorPart50(xCurrent,1d0,hadron)&
      +(/dummy1,dummy1,dummy1,dummy1,dummy1,&
      dummy2,dummy1,dummy1,dummy1,dummy1,dummy1/)*FPDFcurrent
!     write(*,*) 'counter GK=',counter
  else
    convolutionPart=0d0
  end if

  
  !write(*,*) 'gluonMIXTUREPart =', gluonMIXTURE/x  
  !write(*,*) 'TDhat', (deltaPart+singularPart+regularPart)/x
  Common_lowScale50=(deltaPart+convolutionPart)*(1d0/x)
  
end function Common_lowScale50

  
!!Gauss-Kronrod adaptive quadrature, with explicit evaluation at the end point (if converge slow)
!!!---
!!! It calculates int_x0^x1 dz  C[z] f[x/z],
!!!		where f[x/z] is given by function xf, C[z] is given by coeff.
!!!---
recursive function MellinConvolutionVectorPart5(x0,x1,hadron) result(res5)
    integer,intent(in)::hadron
    real(dp),dimension(-5:5)::res5,PDFs,value,eps,epspdf,vg7,vk15
    
    real(dp) :: x0,x1,xm,xr,z,PDFsum,CqMain,CqAnti,CqPrime,CqGluon
    integer :: j,i
    real(dp),dimension(1:parametrizationLength):: var
    real(dp):: alpha,Lmu,Nf,dummy
    real(dp),dimension(-5:5)::F0
    
    xm=0.5d0*(x1+x0)
    xr=0.5d0*(x1-x0)
    
    vg7=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
    vk15=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
    
    Do j=1,15
      
    
      z=xm+xr*Xi_k15(j)
      
    !!!!If mu(x) we have to recalculate Coefficeints every new x!!! This should be very loooong
    if(IsMuXdependent) then
     if(bTcurrent>1d-8) then
	  muCurrent=mu_OPE(z,bTcurrent,c4_global)
	  Lmu=LogMuB(muCurrent,bTcurrent)
	 else
	  muCurrent=mu_OPE(z,10d-8,c4_global)
	  Lmu=LogMuB(muCurrent,10d-8)
     end if
     alpha=As(muCurrent)
     Nf=real(activeNf(muCurrent),dp)

     call Set_CoeffSing1_q_q(alpha,Nf,Lmu)
     call Set_Coeff_q_q(alpha,Nf,Lmu)
     call Set_Coeff_q_qb(alpha,Nf,Lmu)
     call Set_Coeff_q_qp(alpha,Nf,Lmu)
     call Set_Coeff_q_g(alpha,Nf,Lmu)  
    end if
      
      !!! PDFs are together with non-perp func!
      PDFs=xf(xCurrent/z,muCurrent,hadron)
      PDFsum=PDFs(-5)+PDFs(-4)+PDFs(-3)+PDFs(-2)+PDFs(-1)+PDFs(1)+PDFs(2)+PDFs(3)+PDFs(4)+PDFs(5)
      
      counter=counter+1
      
      var=parametrizationString(z)
      
      !! summing regular part
      CqMain=SUM(Coeff_q_q*var)
      CqPrime=SUM(Coeff_q_qp*var)
      CqAnti=SUM(Coeff_q_qb*var)
      CqGluon=SUM(Coeff_q_g*var)
      
      F0=FNP(xCurrent,z,bTcurrent,hadron,lambdaNP)
     
      !!combingin with PDFs
      value=F0*(/&
      CqMain*PDFs(-5)+CqAnti*PDFs(5)+CqPrime*PDFsum+CqGluon*PDFs(0),&
      CqMain*PDFs(-4)+CqAnti*PDFs(4)+CqPrime*PDFsum+CqGluon*PDFs(0),&
      CqMain*PDFs(-3)+CqAnti*PDFs(3)+CqPrime*PDFsum+CqGluon*PDFs(0),&
      CqMain*PDFs(-2)+CqAnti*PDFs(2)+CqPrime*PDFsum+CqGluon*PDFs(0),&
      CqMain*PDFs(-1)+CqAnti*PDFs(1)+CqPrime*PDFsum+CqGluon*PDFs(0),&
      0d0,&
      CqMain*PDFs(1)+CqAnti*PDFs(-1)+CqPrime*PDFsum+CqGluon*PDFs(0),&
      CqMain*PDFs(2)+CqAnti*PDFs(-2)+CqPrime*PDFsum+CqGluon*PDFs(0),&
      CqMain*PDFs(3)+CqAnti*PDFs(-3)+CqPrime*PDFsum+CqGluon*PDFs(0),&
      CqMain*PDFs(4)+CqAnti*PDFs(-4)+CqPrime*PDFsum+CqGluon*PDFs(0),&
      CqMain*PDFs(5)+CqAnti*PDFs(-5)+CqPrime*PDFsum+CqGluon*PDFs(0)/)
      
      !!adding ()_+ part
      dummy=CoeffSing1_q_q(1)/(1d0-z)
      if(order_global>=2) dummy=dummy+CoeffSing1_q_q(2)*LOG(1d0-z)/(1d0-z)
      if(order_global>=3) dummy=dummy+CoeffSing1_q_q(3)*LOG(1d0-z)**2/(1d0-z)

      value=value+dummy*(F0*PDFs-FPDFcurrent)

      vg7=vg7+Wi_g7(j)*value
      vk15=vk15+Wi_k15(j)*value       
    end do
    
    !!!! check the convergance
    eps=ABS(xr*(vg7-vk15)/integralWeight)
    eps(0)=0d0
    !write(*,*) x0,x1, eps
      
    !!!! section for checking integral
      
    if(MAXVAL(eps)>tolerance) then !!!integral not yet convergent
    if(counter>maxIteration) then
    !!!out of counting limit: rise warning, return current result
       
    if(outputLevel>0) call Warning_Raise('Mellin convolution does not converge. Integral evaluation stop after '//&
	  numToStr(maxIteration)//' iterations.',messageCounter,messageTrigger,moduleName)
       
	if(outputLevel>2) then
	  write(*,*) '----- information on last call -----'
	  write(*,*) 'x=',xCurrent,'b=',bTcurrent,'mu=',muCurrent
	  write(*,*) 'iteration=',counter, 'eps=',eps
	  write(*,*) 'weight=',integralWeight
	end if
	  !stop
	res5=xr*vk15
	else
      !!!! we are inside the counting limit
	  if(x1==1d0 .and. 1d0-x0<tolerance) then !!!!small distance to unity
	    !!!! in the case of the integration from x to 1, we check the convergance at 1
	    !!!! if the change of PDF*fNP is small enough we replace the integral, by the exact integral
	    if(IsMuXdependent) then
	      if(bTcurrent>1d-8) then
            muCurrent=mu_OPE(x0,bTcurrent,c4_global)
	      else
            muCurrent=mu_OPE(x0,1d-8,c4_global)
	      end if
	    end if
	    
	    epspdf=ABS((FNP(xCurrent,x0,bTcurrent,hadron,lambdaNP)*xf(xCurrent/x0,muCurrent,hadron)-FPDFcurrent)/integralWeight)
	    epspdf(0)=0      
	    counter=counter+1
	    !here we will add end point integration
	    !write(*,*) x0,x1, epspdf
	    if(MAXVAL(epspdf)<tolerance) then !!! variation is small
        if(IsMuXdependent) then
		if(bTcurrent>1d-8) then
		  muCurrent=mu_OPE(x0,bTcurrent,c4_global)
		  Lmu=LogMuB(muCurrent,bTcurrent)
        else
          muCurrent=mu_OPE(x0,1d-8,c4_global)
		  Lmu=LogMuB(muCurrent,1d-8)
        end if
        alpha=As(muCurrent)
		Nf=real(activeNf(muCurrent),dp)
		call Set_CoeffSing1_q_q(alpha,Nf,Lmu)
		call Set_Coeff_q_q(alpha,Nf,Lmu)
		call Set_Coeff_q_qb(alpha,Nf,Lmu)
		call Set_Coeff_q_qp(alpha,Nf,Lmu)
		call Set_Coeff_q_g(alpha,Nf,Lmu)  
        end if
	    
	    
	     !!!! integrate by usuming that the pdf is flat (+linear)
	      PDFs=(xf(xCurrent/x0,muCurrent,hadron)+PDFcurrent)/2d0
	      PDFsum=PDFs(-5)+PDFs(-4)+PDFs(-3)+PDFs(-2)+PDFs(-1)+PDFs(1)+PDFs(2)+PDFs(3)+PDFs(4)+PDFs(5)
	      counter=counter+1
	      !!!this is integral over vars from x0 to 1, last terms are <10^-7 for x0=1-10^-3
	      var=parametrizationStringAt1(x0)
	      !! summing regular part
	      CqMain=SUM(Coeff_q_q*var)
	      CqPrime=SUM(Coeff_q_qp*var)
	      CqAnti=SUM(Coeff_q_qb*var)
	      CqGluon=SUM(Coeff_q_g*var)
	      
	      F0=FNP(xCurrent,x0,bTcurrent,hadron,lambdaNP)
     
	      !!combingin with PDFs
	      value=F0*(/&
	      CqMain*PDFs(-5)+CqAnti*PDFs(5)+CqPrime*PDFsum+CqGluon*PDFs(0),&
	      CqMain*PDFs(-4)+CqAnti*PDFs(4)+CqPrime*PDFsum+CqGluon*PDFs(0),&
	      CqMain*PDFs(-3)+CqAnti*PDFs(3)+CqPrime*PDFsum+CqGluon*PDFs(0),&
	      CqMain*PDFs(-2)+CqAnti*PDFs(2)+CqPrime*PDFsum+CqGluon*PDFs(0),&
	      CqMain*PDFs(-1)+CqAnti*PDFs(1)+CqPrime*PDFsum+CqGluon*PDFs(0),&
	      0d0,&
	      CqMain*PDFs(1)+CqAnti*PDFs(-1)+CqPrime*PDFsum+CqGluon*PDFs(0),&
	      CqMain*PDFs(2)+CqAnti*PDFs(-2)+CqPrime*PDFsum+CqGluon*PDFs(0),&
	      CqMain*PDFs(3)+CqAnti*PDFs(-3)+CqPrime*PDFsum+CqGluon*PDFs(0),&
	      CqMain*PDFs(4)+CqAnti*PDFs(-4)+CqPrime*PDFsum+CqGluon*PDFs(0),&
	      CqMain*PDFs(5)+CqAnti*PDFs(-5)+CqPrime*PDFsum+CqGluon*PDFs(0)/)
	      !!! AV: 25.06.22
	      !!! the (..)_+ part is tricky.
	      !!! I estimate int_x0^1 (f[xc/x]-f[xc])/(1-x)g[x] as (f[xc/x0]-f[xc])/(1-x0) int_x0^1 g[x]
	      !!! it follows from the expansion of f[xc/x] at x->1, and comparison to the expansion of f[xc/x0] at x0->1
	      !!! it is justified if x0->1, and f[xc] is smooth
          dummy=CoeffSing1_q_q(1)
          if(order_global>=2) dummy=dummy+CoeffSing1_q_q(2)*(LOG(1d0-x0)-1d0)
          if(order_global>=3) dummy=dummy+CoeffSing1_q_q(3)*(LOG(1d0-x0)**2-2d0*LOG(1d0-x0)+2d0)

          value=value+dummy*(F0*PDFs-FPDFcurrent)
	      res5=value
	    else
	      res5=MellinConvolutionVectorPart5(x0,xm,hadron)+MellinConvolutionVectorPart5(xm,x1,hadron)
	    end if
	  else
	    res5=MellinConvolutionVectorPart5(x0,xm,hadron)+MellinConvolutionVectorPart5(xm,x1,hadron)
	  end if
      end if
      else   !!!integral converges
       res5=xr*vk15
      end if
  end function MellinConvolutionVectorPart5
  
!!!Gauss-Kronrod adaptive quadrature, with explicit evaluation at the end point (if converge slow)
!!!---
!!! It calculates int_x0^x1 dz  C[z] f[x/z],
!!!		where f[x/z] is given by function xf, C[z] is given by coeff.
!!!---
recursive function MellinConvolutionVectorPart50(x0,x1,hadron) result(res5)
    integer,intent(in)::hadron
    real(dp),dimension(-5:5)::res5,PDFs,value,vg7,vk15,eps,epspdf,addV
    
    real(dp) :: x0,x1,xm,xr,z,PDFsum,CqMain,CqAnti,CqPrime,CqGluon,CgMain,CgQuark
    integer :: j,i
    real(dp),dimension(1:parametrizationLength):: var
    real(dp):: alpha,Lmu,Nf,dummy
    real(dp),dimension(-5:5)::F0
    
    xm=0.5d0*(x1+x0)
    xr=0.5d0*(x1-x0)
    
    vg7=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
    vk15=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
    
    Do j=1,15
      z=xm+xr*Xi_k15(j)
      
      !!!!If mu(z) we have to recalculate Coefficeints every new z!!! This should be very loooong
      if(IsMuXdependent) then
      if(bTcurrent>1d-8) then
        muCurrent=mu_OPE(z,bTcurrent,c4_global)
        Lmu=LogMuB(muCurrent,bTcurrent)
      else
        muCurrent=mu_OPE(z,10d-8,c4_global)
        Lmu=LogMuB(muCurrent,10d-8)
      end if
      alpha=As(muCurrent)
      Nf=real(activeNf(muCurrent),dp)
      call Set_CoeffSing1_q_q(alpha,Nf,Lmu)
      call Set_Coeff_q_q(alpha,Nf,Lmu)
      call Set_Coeff_q_qb(alpha,Nf,Lmu)
      call Set_Coeff_q_qp(alpha,Nf,Lmu)
      call Set_Coeff_q_g(alpha,Nf,Lmu)
      call Set_Coeff_g_q(alpha,Nf,Lmu)
      call Set_Coeff_g_g(alpha,Nf,Lmu)
      end if
      
      !!! PDFs are together with non-perp func!
      PDFs=xf(xCurrent/z,muCurrent,hadron)
      PDFsum=PDFs(-5)+PDFs(-4)+PDFs(-3)+PDFs(-2)+PDFs(-1)+PDFs(1)+PDFs(2)+PDFs(3)+PDFs(4)+PDFs(5)
      counter=counter+1
      
      var=parametrizationString(z)
      
      !! summing regular part
      CqMain=SUM(Coeff_q_q*var)
      CqPrime=SUM(Coeff_q_qp*var)
      CqAnti=SUM(Coeff_q_qb*var)
      CqGluon=SUM(Coeff_q_g*var)
      
      CgMain=SUM(Coeff_g_g*var)
      CgQuark=SUM(Coeff_g_q*var)
      
     F0=FNP(xCurrent,z,bTcurrent,hadron,lambdaNP)
      !!combingin with PDFs
      value=F0*(/&
      CqMain*PDFs(-5)+CqAnti*PDFs(5)+CqPrime*PDFsum+CqGluon*PDFs(0),&
      CqMain*PDFs(-4)+CqAnti*PDFs(4)+CqPrime*PDFsum+CqGluon*PDFs(0),&
      CqMain*PDFs(-3)+CqAnti*PDFs(3)+CqPrime*PDFsum+CqGluon*PDFs(0),&
      CqMain*PDFs(-2)+CqAnti*PDFs(2)+CqPrime*PDFsum+CqGluon*PDFs(0),&
      CqMain*PDFs(-1)+CqAnti*PDFs(1)+CqPrime*PDFsum+CqGluon*PDFs(0),&
      CgMain*PDFs(0)+CgQuark*PDFsum,&
      CqMain*PDFs(1)+CqAnti*PDFs(-1)+CqPrime*PDFsum+CqGluon*PDFs(0),&
      CqMain*PDFs(2)+CqAnti*PDFs(-2)+CqPrime*PDFsum+CqGluon*PDFs(0),&
      CqMain*PDFs(3)+CqAnti*PDFs(-3)+CqPrime*PDFsum+CqGluon*PDFs(0),&
      CqMain*PDFs(4)+CqAnti*PDFs(-4)+CqPrime*PDFsum+CqGluon*PDFs(0),&
      CqMain*PDFs(5)+CqAnti*PDFs(-5)+CqPrime*PDFsum+CqGluon*PDFs(0)/)
      
      !!adding ()_+ part
      dummy=CoeffSing1_q_q(1)/(1d0-z)
      if(order_global>=2) dummy=dummy+CoeffSing1_q_q(2)*LOG(1d0-z)/(1d0-z)
      if(order_global>=3) dummy=dummy+CoeffSing1_q_q(3)*LOG(1d0-z)**2/(1d0-z)

      addV=dummy*(F0*PDFs-FPDFcurrent)

      dummy=CoeffSing1_g_g(1)/(1d0-z)
      if(order_global>=2) dummy=dummy+CoeffSing1_g_g(2)*LOG(1d0-z)/(1d0-z)
      if(order_global>=3) dummy=dummy+CoeffSing1_g_g(3)*LOG(1d0-z)**2/(1d0-z)

      addV(0)=dummy*(F0(0)*PDFs(0)-FPDFcurrent(0))

      value=value+addV
      
      !write(*,*)CqMain,CqAnti,CqPrime,singCoeff,singCoeffLog      
      vg7=vg7+Wi_g7(j)*value
      vk15=vk15+Wi_k15(j)*value      
    end do
    
      eps=ABS(xr*(vg7-vk15)/integralWeight)
    
      !write(*,*) x0,x1, eps(0), xr*vk15(0)
    
      if(MAXVAL(eps)>tolerance) then
        if(counter>maxIteration) then
        if(outputLevel>0) call Warning_Raise('Mellin convolution does not converge. Integral evaluation stop after '//&
            numToStr(maxIteration)//' iterations.',messageCounter,messageTrigger,moduleName)
        if(outputLevel>1) then
          write(*,*) '----- information on last call -----'
          write(*,*) 'x=',xCurrent,'b=',bTcurrent,'mu=',muCurrent
          write(*,*) 'iteration=',counter, 'eps=',eps
          write(*,*) 'weight=',integralWeight
        end if

        res5=xr*vk15
      else
	  if((1d0-x1)<1d-12 .and. (1d0-x0)<tolerance) then
	    !!!! in the case of the integration from x to 1, we check the convergance at 1
	    !!!! if the change of PDF*fNP is small enough we replace hte integral, by the exact integral
	    
	    if(IsMuXdependent) then
	      if(bTcurrent>1d-8) then
            muCurrent=mu_OPE(x0,bTcurrent,c4_global)
	      else
            muCurrent=mu_OPE(x0,1d-8,c4_global)
	      end if
	    end if
	    
	    epspdf=ABS((FNP(xCurrent,x0,bTcurrent,hadron,lambdaNP)*xf(xCurrent/x0,muCurrent,hadron)-FPDFcurrent)/integralWeight)
	    counter=counter+1
	    !here we will add end point integration
	    !write(*,*) "====",x0,x1, epspdf(0)
	    if(MAXVAL(epspdf)<tolerance) then !!! variation is small
	      if(IsMuXdependent) then
            if(bTcurrent>1d-8) then
              muCurrent=mu_OPE(x0,bTcurrent,c4_global)
              Lmu=LogMuB(muCurrent,bTcurrent)
            else
              muCurrent=mu_OPE(x0,1d-8,c4_global)
              Lmu=LogMuB(muCurrent,1d-8)
            end if
            alpha=As(muCurrent)
            Nf=real(activeNf(muCurrent),dp)
            call Set_CoeffSing1_q_q(alpha,Nf,Lmu)
            call Set_Coeff_q_q(alpha,Nf,Lmu)
            call Set_Coeff_q_qb(alpha,Nf,Lmu)
            call Set_Coeff_q_qp(alpha,Nf,Lmu)
            call Set_Coeff_q_g(alpha,Nf,Lmu)
            call Set_Coeff_g_q(alpha,Nf,Lmu)
            call Set_Coeff_g_g(alpha,Nf,Lmu)
          end if

          !!! I approximate function f(xC/z) by mean value f(xC/x0)-f(xC)
	      PDFs=(xf(xCurrent/x0,muCurrent,hadron)+PDFcurrent)/2d0
	      PDFsum=PDFs(-5)+PDFs(-4)+PDFs(-3)+PDFs(-2)+PDFs(-1)+PDFs(1)+PDFs(2)+PDFs(3)+PDFs(4)+PDFs(5)
	      counter=counter+1
	      !!!this is integral over vars from x0 to 1, last terms are <10^-7 for x0=1-10^-3
	      
	      var=parametrizationStringAt1(x0)
	      !! summing regular part
	      CqMain=SUM(Coeff_q_q*var)
	      CqPrime=SUM(Coeff_q_qp*var)
	      CqAnti=SUM(Coeff_q_qb*var)
	      CqGluon=SUM(Coeff_q_g*var)
	      CgMain=SUM(Coeff_g_g*var)
	      CgQuark=SUM(Coeff_g_q*var)
	      
	      F0=FNP(xCurrent,x0,bTcurrent,hadron,lambdaNP)
     
	      !!combingin with PDFs
	      value=F0*(/&
	      CqMain*PDFs(-5)+CqAnti*PDFs(5)+CqPrime*PDFsum+CqGluon*PDFs(0),&
	      CqMain*PDFs(-4)+CqAnti*PDFs(4)+CqPrime*PDFsum+CqGluon*PDFs(0),&
	      CqMain*PDFs(-3)+CqAnti*PDFs(3)+CqPrime*PDFsum+CqGluon*PDFs(0),&
	      CqMain*PDFs(-2)+CqAnti*PDFs(2)+CqPrime*PDFsum+CqGluon*PDFs(0),&
	      CqMain*PDFs(-1)+CqAnti*PDFs(1)+CqPrime*PDFsum+CqGluon*PDFs(0),&
	      CgMain*PDFs(0)+CgQuark*PDFsum,&
	      CqMain*PDFs(1)+CqAnti*PDFs(-1)+CqPrime*PDFsum+CqGluon*PDFs(0),&
	      CqMain*PDFs(2)+CqAnti*PDFs(-2)+CqPrime*PDFsum+CqGluon*PDFs(0),&
	      CqMain*PDFs(3)+CqAnti*PDFs(-3)+CqPrime*PDFsum+CqGluon*PDFs(0),&
	      CqMain*PDFs(4)+CqAnti*PDFs(-4)+CqPrime*PDFsum+CqGluon*PDFs(0),&
	      CqMain*PDFs(5)+CqAnti*PDFs(-5)+CqPrime*PDFsum+CqGluon*PDFs(0)/)

	      !!! AV: 25.06.22
	      !!! the (..)_+ part is tricky.
	      !!! I estimate int_x0^1 (f[xc/x]-f[xc])/(1-x)g[x] as (f[xc/x0]-f[xc])/(1-x0) int_x0^1 g[x]
	      !!! it follows from the expansion of f[xc/x] at x->1, and comparison to the expansion of f[xc/x0] at x0->1
	      !!! it is justified if x0->1, and f[xc] is smooth
          dummy=CoeffSing1_q_q(1)
          if(order_global>=2) dummy=dummy+CoeffSing1_q_q(2)*(LOG(1d0-x0)-1d0)
          if(order_global>=3) dummy=dummy+CoeffSing1_q_q(3)*(LOG(1d0-x0)**2-2d0*LOG(1d0-x0)+2d0)

          addV=dummy*(F0*PDFs-FPDFcurrent)

          dummy=CoeffSing1_g_g(1)
          if(order_global>=2) dummy=dummy+CoeffSing1_g_g(2)*(LOG(1d0-x0)-1d0)
          if(order_global>=3) dummy=dummy+CoeffSing1_g_g(3)*(LOG(1d0-x0)**2-2d0*LOG(1d0-x0)+2d0)

          addV(0)=dummy*(F0(0)*PDFs(0)-FPDFcurrent(0))

          value=value+addV

	      res5=value
	    else
	      res5=MellinConvolutionVectorPart50(x0,xm,hadron)+MellinConvolutionVectorPart50(xm,x1,hadron)
	    end if
	  else
	    res5=MellinConvolutionVectorPart50(x0,xm,hadron)+MellinConvolutionVectorPart50(xm,x1,hadron)
	  end if
      end if
      else          
       res5=xr*vk15
      end if

  end function MellinConvolutionVectorPart50
  
