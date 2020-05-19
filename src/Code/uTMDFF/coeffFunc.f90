!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!
!!!!!!!!!       part of uTMDFF module for artemide
!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!-------------------------------------------------------------------------------------------------------
!!!-----------------------------Routines for the mathing coefficient--------------------------------------
!!!--------------The order is accumulative pertrubative order of coefficient =0,1,2 (LO,NLO,NNLO)---------
!!!-------------------------------------------------------------------------------------------------------

  !!! the function which contains the functions of parameterizations
  function parametrizationString(z)
  real(dp)::z,lx,l1x
  real(dp),dimension(1:parametrizationLength)::parametrizationString
      lx=Log(z)
      l1x=Log(1d0-z)
      parametrizationString=(/ l1x,l1x**2, l1x**3,&
	1d0/z, lx/z, lx**2/z, lx**3/z, &
	lx, lx**2, lx**3,&
	1d0, z, z**2,&
	z*lx/(1d0-z), z*lx, (z**2)*lx,&
	z*(lx**2)/(1d0-z), z*lx**2,  z*lx**3, &
	(lx/(1d0-z)+1d0)*l1x, lx*l1x,  z*lx*l1x,&
	(1d0-z)/z*l1x, (1d0-z)*l1x, (1d0-z)*l1x**2/)
  
  end function parametrizationString
  
  !!! the function which contains 
    !!! int_z^1 parameterization at values of z -> 1
    !!! it is used to estimate integration error at z~1
  function parametrizationStringAt1(z0)
  real(dp)::z0
  real(dp),dimension(1:parametrizationLength)::parametrizationStringAt1
    parametrizationStringAt1=(/(1d0-z0)*(-1d0+Log(1d0-z0)),(1d0-z0)*(2d0-2d0*Log(1d0-z0)+Log(1d0-z0)**2),&
	      (1d0-z0)*(-6d0+6d0*Log(1d0-z0)-3d0*Log(1d0-z0)**2+Log(1d0-z0)**3),&
	      0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0,&
	      (1d0-z0)**2*Log(1d0-z0),0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0, 0d0/)
  
  end function parametrizationStringAt1

  


!!!!Each coefficient is split to delta, sing x->1, regular
  
  !!!!!coefficient function q<-q delta-part
  function C_q_q_delta(alpha,Nf,Lmu)
  real(dp)::C_q_q_delta,Nf,alpha,Lmu
  
 C_q_q_delta=1d0
 
  if(order_global>=1) then
      !(checked 27.02.19 AV)
      C_q_q_delta=C_q_q_delta+alpha*(-4d0/3d0*zeta2-4d0*Lmu)
  end if
  if(order_global>=2) then
    !(checked 27.02.19 AV)
     C_q_q_delta=C_q_q_delta+alpha*alpha*(&
     -2416d0/81d0 + (Lmu**2)*(-14d0 + 4d0*Nf/3d0 - 128d0*zeta2/9d0) - 134d0*zeta2/3d0 + 448d0*zeta3/9d0 &
     + Nf*(352d0/243d0 + 20d0*zeta2/9d0 + 56d0*zeta3/27d0) +  Lmu*(-14d0 - 140d0*zeta2/3d0 + Nf*(4d0/9d0&
     + 40d0*zeta2/9d0) + 16d0*zeta3/3d0) + 2360d0*zeta4/9d0 )
  end if
  
  end function C_q_q_delta
  
  !!!!!coefficient function g<-g delta-part
  function C_g_g_delta(alpha,Nf,Lmu)
  real(dp)::C_g_g_delta,Nf,alpha,Lmu
  
  C_g_g_delta=1d0
 
  if(order_global>=1) then
      !(checked 27.02.19 AV)
      C_g_g_delta=C_g_g_delta+alpha*(-3d0*zeta2+(-11d0+2d0/3d0*Nf)*Lmu)
  end if
  if(order_global>=2) then
    !(checked 27.02.19 AV)
     C_g_g_delta=C_g_g_delta+alpha*alpha*(&
     -112d0 - 56d0*(Nf**2)/81d0 - 201d0*zeta2/2d0 - 72d0*(Lmu**2)*zeta2 + Lmu*(-96d0 + 32d0*Nf/3d0 &
     - 108d0*zeta3) + Nf*(548d0/27d0 + 5d0*zeta2 - 28d0*zeta3/3d0) + 154d0*zeta3 + 2385d0*zeta4/4d0)
  end if
  end function C_g_g_delta
  
  !!!!!coefficient function q<-q singular-part  (1/(1-x)_+,(Log(1-x)/(1-x))_+)
  subroutine Set_CoeffSing1_q_q(alpha,Nf,Lmu)
  real(dp)::Nf,alpha,Lmu,s1,s2
    
  s1=0d0!!!coeff 1/(1-x)
  s2=0d0!!!coeff log(1-x)/(1-x)
  

  if(order_global>=1) then    
     !(checked 27.02.19 AV)
    s1=s1+alpha*(-16d0/3d0*Lmu)
  end if
  if(order_global>=2) then
    !(checked 27.02.19 AV)
    s1=s1+alpha*alpha*(&
    -3232d0/27d0 + 448d0*Nf/81d0 + (Lmu**2)*(-8d0 + 16d0*Nf/9d0) + &
    Lmu*(-1072d0/9d0 + 160d0*Nf/27d0 + 352d0*zeta2/9d0) + 112d0*zeta3)
     s2=s2+alpha*alpha*(256d0/9d0*Lmu**2)
  end if
  
  CoeffSing1_q_q=(/s1,s2/)
  
  end subroutine Set_CoeffSing1_q_q
  
  !!!!!coefficient function g<-g singular-part  (1/(1-x)_+,(Log(1-x)/(1-x))_+)
  subroutine Set_CoeffSing1_g_g(alpha,Nf,Lmu)
  real(dp)::Nf,alpha,Lmu,s1,s2
    
  s1=0d0!!!coeff 1/(1-x)
  s2=0d0!!!coeff log(1-x)/(1-x)
 
  if(order_global>=1) then  
   !(checked 27.02.19 AV)
    s1=s1+alpha*(-12d0)*Lmu
  end if
  if(order_global>=2) then
    !(checked 27.02.19 AV)
    s1=s1+alpha*alpha*(&
    -808d0/3d0 + (Lmu**2)*(66d0 - 4d0*Nf) + 112d0*Nf/9d0 +&
      Lmu*(-268d0 + 40d0*Nf/3d0 + 108d0*zeta2) + 252d0*zeta3)
!     
     s2=s2+alpha*alpha*144d0*(Lmu**2)
  end if
  
  
  CoeffSing1_g_g=(/s1,s2/)
  
  end subroutine Set_CoeffSing1_g_g
  
  !!!!!coefficient function q->q
  subroutine Set_Coeff_q_q(alpha,Nf,Lmu)  
  real(dp)::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_q_q=(/&
  0d0,0d0,0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
  0d0,0d0,0d0,0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
  0d0,0d0,0d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
  0d0,0d0,0d0,&		!1 (exact), z, z^2
  0d0,0d0,0d0,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
  0d0,0d0,0d0,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
  0d0,0d0,0d0,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
  0d0,0d0,0d0/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
  if(order_global>=1) then
    !(checked 27.02.19 AV)
    Coeff_q_q=Coeff_q_q+alpha*(/&
      0d0,0d0,0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
      0d0,0d0,0d0,0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
      16d0/3d0,0d0,0d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
      8d0/3d0*(1d0+Lmu),8d0/3d0*(Lmu-1d0), 0d0,&		!1 (exact), z, z^2
      32d0/3d0,-16d0/3d0,0d0,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
      0d0,0d0,0d0,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
      0d0,0d0,0d0,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
      0d0,0d0,0d0/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
    

  !------The kernels are calculated in mathematica
  if(order_global>=2) then
      !(checked 27.02.19 AV)
     Coeff_q_q=Coeff_q_q+alpha*alpha*(/&
      -200d0/9d0 + 256d0*Lmu/3d0 - 256d0*(Lmu**2)/9d0, 64d0/9d0, 0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
      0d0,0d0,0d0,0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
      1496d0/9d0-32d0*(Lmu**2)/9d0+Lmu*(248d0/9d0-16d0*Nf/9d0)-8d0*Nf-560d0*zeta2/9d0,&
	  -130d0/9d0+88d0*Lmu/9d0+4d0*Nf/9d0,-140d0/27d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
      (Lmu**2)*(-28d0/9d0 - 8d0*Nf/9d0) - 296d0*Nf/81d0 + Lmu*(64d0/3d0 - 32d0*Nf/27d0 - 176d0*zeta2/9d0),&
	  (Lmu**2)*(100d0/9d0 - 8d0*Nf/9d0) - 152d0*Nf/81d0 + Lmu*(880d0/9d0 - 128d0*Nf/27d0 - 176d0*zeta2/9d0), 0d0,&		!1 (exact), z, z^2
      -128d0*(Lmu**2)/9d0 + Lmu*(-16d0/3d0 - 32d0*Nf/9d0) - 80d0*Nf/9d0, &
	  32d0*(Lmu**2)/3d0 + 8d0*Nf + Lmu*(-24d0 + 16d0*Nf/9d0), 0d0,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
      368d0*Lmu/9d0 + 8d0*Nf/9d0, -280d0*Lmu/9d0 - 4d0*Nf/9d0, 0d0,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
      -256d0*Lmu/3d0, 128d0*Lmu/3d0, 128d0*Lmu/3d0,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
      0d0, -128d0*Lmu/9d0 + 128d0*(Lmu**2)/9d0, 0d0/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
      
      !! adding approximate part
      Coeff_q_q=Coeff_q_q+alpha*alpha*(/&
      0d0,0d0,0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
      0d0,0d0,0d0,0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
      0d0,0d0,0d0,&
      -301.03247439776976d0, -989.2167272286393d0, 82.59056818996338d0,&
      -1063.98482846164d0, 206.28577290245227d0, -18.651271690975136d0,&
      83.00296625888389d0, -70.52319745715631d0, -4.911975080877064d0,&
      -1105.7693500845382d0, 327.376932797376d0, -109.45363459015105d0,&
      -174.6471655693391d0, -112.83673919345797d0, -3.5294575557396084d0/)
      
  end if
    
   !write(*,*) 'regularPart=', regularPart/x
  end if
  end subroutine Set_Coeff_q_q
  
   !!!!!coefficient function q->g
  subroutine Set_Coeff_q_g(alpha,Nf,Lmu)  
  real(dp)::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_q_g=(/&
  0d0,0d0,0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
  0d0,0d0,0d0,0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
  0d0,0d0,0d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
  0d0,0d0,0d0,&		!1 (exact), z,z^2
  0d0,0d0,0d0,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
  0d0,0d0,0d0,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
  0d0,0d0,0d0,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
  0d0,0d0,0d0/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
  if(order_global>=1) then
    !(checked 27.02.19 AV)
    Coeff_q_g=Coeff_q_g+alpha*(/&
      0d0,0d0,0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
      -16d0/3d0*Lmu,32d0/3d0,0d0,0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
      -32d0/3d0,0d0,0d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
      16d0/3d0*Lmu,8d0/3d0*(1d0-Lmu), 0d0,&		!1 (exact), z, z^2
      0d0,16d0/3d0,0d0,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
      0d0,0d0,0d0,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
      0d0,0d0,0d0,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
      0d0,0d0,0d0/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
    

  !------The kernels are calculated in mathematica
    if(order_global>=2) then
      Coeff_q_g=Coeff_q_g+alpha*alpha*(/&
      -40d0/9d0 - 128d0*Lmu/9d0 + 208d0*(Lmu**2)/9d0 - 80d0*zeta2/3d0, -40d0/9d0 + 80d0*Lmu/9d0, 40d0/27d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
      12512d0/27d0 - 248d0*(Lmu**2)/3d0 + Lmu*(112d0/3d0 - 64d0*zeta2) - 352d0*zeta2/3d0 - 352d0*zeta3,&
	  400d0/3d0 + 992d0*Lmu/3d0 - 32d0*Lmu**2, -848d0/3d0 + 128d0*Lmu, -320d0/3d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
      1144d0/3d0 - 64d0*Lmu - 224d0*(Lmu**2)/9d0, 64d0/3d0 + 224d0*Lmu/3d0, -1232d0/27d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
      640d0*(Lmu**2)/9d0, 56d0*(Lmu**2)/9d0, 32d0*(Lmu**2)/3d0,&		!1 (exact), z, z^2
      0d0,-320d0/9d0*Lmu**2, 0d0,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
      0d0, 0d0, 0d0,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
      0d0, 0d0, 0d0,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
      416d0*(Lmu**2)/9d0, -208d0*(Lmu**2)/9d0, 0d0/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
      
      !! adding approximate part
      Coeff_q_g=Coeff_q_g+alpha*alpha*(/&
      0d0,0d0,0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
      0d0,0d0,0d0,0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
      0d0,0d0,0d0,&
      -11519.346897372414d0 - 740.2350001396422d0*Lmu, -34241.466106099186d0 - 1895.8089140328416d0*Lmu,&
		5326.770104932414d0 + 117.62109028406978d0*Lmu,&
      -40601.87518176106d0 - 2457.115336720303d0*Lmu, 4178.463030904903d0 + 197.03426376642062d0*Lmu, &
		-1705.6350033291087d0 - 90.35072880713321d0*Lmu,&
      -966.5754106411847d0 - 65.94510382678313d0*Lmu, 1144.6267544753136d0 + 195.6628827312015d0*Lmu, &
		-84.66732541780037d0 + 0.3444214250305629d0*Lmu,&
      -11393.035115581788d0 - 812.6612648812184d0*Lmu, -25857.49295712562d0 - 1326.6218862880226d0*Lmu, &
		-10601.55795204891d0 - 581.4253894551872d0*Lmu,&
      -12368.214954397781d0 - 628.2350001396422d0*Lmu, -29991.60756399795d0 - 1608.4296229988213d0*Lmu, &
		-5.282535747460972d0 + 9.108108570274817d0*Lmu/)
    end if
    
   !write(*,*) 'regularPart=', regularPart/x
  end if
  end subroutine Set_Coeff_q_g
  
   !!!!!coefficient function g->q
  subroutine Set_Coeff_g_q(alpha,Nf,Lmu)  
  real(dp)::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_g_q=(/&
  0d0,0d0,0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
  0d0,0d0,0d0,0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
  0d0,0d0,0d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
  0d0,0d0,0d0,&		!1 (exact), z, z^2
  0d0,0d0,0d0,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
  0d0,0d0,0d0,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
  0d0,0d0,0d0,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
  0d0,0d0,0d0/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
  if(order_global>=1) then
    !(checked 27.02.19 AV)
    Coeff_g_q=Coeff_g_q+alpha*(/&
      0d0,0d0,0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
      0d0,0d0,0d0,0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
      2d0,0d0,0d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
      -Lmu,2d0*(1d0+Lmu), -2d0*(1d0+Lmu),&		!1 (exact), z,  z^2
      0d0,-4d0,4d0,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
      0d0,0d0,0d0,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
      0d0,0d0,0d0,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
      0d0,0d0,0d0/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
    

  !------The kernels are calculated in mathematica
    if(order_global>=2) then
    Coeff_g_q=Coeff_g_q+alpha*alpha*(/&
      -44d0/3d0 + 26d0*(Lmu**2)/3d0 + 10d0*Nf/9d0 + Lmu*(-14d0 + 4d0*Nf/3d0) + 10d0*zeta2,&
	    -7d0/2d0 - 10d0*Lmu/3d0 + Nf/3d0, -5d0/9d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
      -148d0/9d0 - 4d0*Lmu + 4d0*Lmu**2 + 8d0*zeta2, 8d0 - 16d0*Lmu, 16d0, 0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
      -157d0/3d0 + 10d0*Lmu + 14d0*(Lmu**2)/3d0 - 10d0*Nf/3d0 - 100d0*zeta2/3d0, -73d0/3d0 - 14d0*Lmu + Nf/3d0, 77d0/9d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
      -5042.153939120427d0 - 1531.8031611907763d0*Lmu + 120.72027508333485d0*Lmu**2 &
      - 66.62886623235278d0*Nf + 20.43559788824016d0*Lmu*Nf, &
      303.86339171965324d0 - 12985.247663082215d0*Lmu + 913.3142650671598d0*Lmu**2 &
      - 317.3462434009907d0*Nf + 133.2962117132032d0*Lmu*Nf, &
      3703.409480057778d0 - 412.2734474693585d0*Lmu - 17.358284535835764d0*Lmu**2 &
      + 21.346846003738534d0*Nf + 5.209836558024565d0*Lmu*Nf,&		!1 (exact), z, z^2
      -1059.0639607339292d0 - 14942.360693014281d0*Lmu + 1018.6762556146588*Lmu**2 &
      - 361.41246956998253d0*Nf + 156.7194239372457d0*Lmu*Nf, &
      -8193.44131612091d0 + 6625.5091321556465d0*Lmu - 408.80693588618846d0*Lmu**2 &
      + 116.14537890481989d0*Nf - 66.99593883060648d0*Lmu*Nf, &
      -1046.4745476064809d0 + 142.78460715142586d0*Lmu + 4.0942238266959485d0*Lmu**2 &
      - 18.841069748305667d0*Nf + 1.4503934130147407d0*Lmu*Nf,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
      -1458.9834262434276d0 + 407.23365449967133d0*Lmu - 21.2552050492032d0*Lmu**2 &
      + 5.214301688311585d0*Nf - 3.2700315397596995d0*Lmu*Nf,& 
      1590.9730212100676d0 - 528.3020502418173d0*Lmu + 22.83469849939765d0*Lmu**2 &
      - 6.201059346425027d0*Nf + 3.513030532048041d0*Lmu*Nf, &
      73.07631646309086d0 - 1.3686522702588848d0*Lmu + 0.07440056656025083d0*Lmu**2 &
      - 0.014948523947447324d0*Nf + 0.011446240977035602d0*Lmu*Nf,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
      10420.770382128236d0 - 10241.895319402418d0*Lmu + 638.9211015085228d0*Lmu**2 &
      - 200.7638366136303d0*Nf + 98.29555403916862d0*Lmu*Nf, &
      -20512.88777168424d0 + 2160.2971291961044d0*Lmu - 67.7001647874525d0*Lmu**2 &
      - 49.38978687382705d0*Nf - 10.415409899902482d0*Lmu*Nf, &
      -4683.056024325436d0 - 1584.4478029479442d0*Lmu + 158.259017038288d0*Lmu**2 &
      - 58.30447555159878d0*Nf + 24.347541101116374d0*Lmu*Nf,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
      -4929.425576232434d0 - 1575.4433993058212d0*Lmu + 118.38694175000153d0*Lmu**2 &
      - 65.41307217273041d0*Nf + 18.21337566601794d0*Lmu*Nf, &
      -15198.304463044708d0 - 2938.3075818393627d0*Lmu + 274.35887368634087d0*Lmu**2 &
      - 144.13788267277235d0*Nf + 42.20905754379859d0*Lmu*Nf, &
      -11.384648952243019d0 + 6.516317328981653d0*Lmu + 0.004805761454612681d0*Lmu**2 &
      - 0.6611076391799198d0*Nf + 0.0007393479176680426d0*Lmu*Nf/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
      
    end if
    
   !write(*,*) 'regularPart=', regularPart/x
  end if
  end subroutine Set_Coeff_g_q
  
     !!!!!coefficient function g->g
  subroutine Set_Coeff_g_g(alpha,Nf,Lmu)  
  real(dp)::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_g_g=(/&
  0d0,0d0,0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
  0d0,0d0,0d0,0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
  0d0,0d0,0d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
  0d0,0d0,0d0,&		!1 (exact), z, z^2
  0d0,0d0,0d0,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
  0d0,0d0,0d0,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
  0d0,0d0,0d0,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
  0d0,0d0,0d0/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
  if(order_global>=1) then
    !(checked 27.02.19 AV)
    Coeff_g_g=Coeff_g_g+alpha*(/&
      0d0,0d0,0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
      -12d0*Lmu,24d0,0d0,0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
      -24d0,0d0,0d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
      24d0*Lmu,-12d0*Lmu, 12d0*Lmu,&		!1 (exact), z,  z^2
      24d0,24d0,-24d0,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
      0d0,0d0,0d0,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
      0d0,0d0,0d0,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
      0d0,0d0,0d0/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
    

  !------The kernels are calculated in mathematica
    if(order_global>=2) then
      Coeff_g_g=Coeff_g_g+alpha*alpha*(/&
      -6d0 + 432d0*Lmu - 144d0*Lmu**2 + 2d0*Nf, 36d0, 0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
      3268d0/3d0 + Lmu**2*(-198d0 - 4d0*Nf/9d0) + 278d0*Nf/81d0 - 264d0*zeta2 + Lmu*(412d0*Nf/27d0 + 36d0*zeta2) - 792d0*zeta3,& 
	  536d0 - 72d0*Lmu**2 - 244d0*Nf/9d0 + Lmu*(792d0 + 16d0*Nf/9d0) - 360d0*zeta2,&
	  -660d0 + 288d0*Lmu - 8d0*Nf/9d0, -240d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
      819d0 + Lmu*(-12d0 - 24d0*Nf) - 50d0*Nf + Lmu**2*(-72d0 + 16d0*Nf/3d0) + 360d0*zeta2, &
	  33d0 + Lmu*(216d0 - 16d0*Nf) + 22d0*Nf/3d0, -132d0 + 88d0*Nf/9d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
      -4536.42238065816d0 + 17316.05106770077d0*Lmu - 899.5223166525981d0*Lmu**2 &
	  - 163.95224892156722d0*Nf - 17.333333050953353d0*Lmu*Nf + 10.666666751525197d0*Lmu**2*Nf, &
      18080.936974539778d0 + 118008.31936555808d0*Lmu - 7515.995501919461d0*Lmu**2 &
	  - 625.9193713983447d0*Nf - 6.666666026541869d0*Lmu*Nf - 6.666666410951244d0*Lmu**2*Nf, &
      14139.810571485807d0 - 2400.2117629303043d0*Lmu + 84.6688450985717d0*Lmu**2 &
	  - 184.19497231902477d0*Nf - 4.5925927828634645d0*Lmu*Nf + 0.4444444032182875d0*Lmu**2*Nf,&		!1 (exact), z, z^2
      27937.822349639453d0 + 132629.0291759547d0*Lmu - 8534.848973473487d0*Lmu**2 &
	  - 975.5233827623935d0*Nf + 16.0000007322339d0*Lmu*Nf + 2.993477955391831d-7*Lmu**2*Nf, &
      -35293.57801079519d0 - 52643.48365504077d0*Lmu + 3401.7806883476474d0*Lmu**2 &
	  + 498.8898228158924d0*Nf + 18.66666678376387d0*Lmu*Nf + 5.333333298547101d0*Lmu**2*Nf, &
      -5939.483006875214d0 + 1055.004893054739d0*Lmu - 6.321251200921099d0*Lmu**2 &
	  + 81.15085632871367d0*Nf - 15.999999931023279d0*Lmu*Nf + 1.4357365397679303d-8*Lmu**2*Nf,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
      -6775.737946774227d0 - 4375.194862408091d0*Lmu + 176.58169722761855d0*Lmu**2 &
	  + 30.343871838863492d0*Nf + 5.784855685667105d-8*Lmu*Nf + 5.784855685667105d-9*Lmu**2*Nf, &
      7209.990192354363d0 + 5748.449232000771d0*Lmu - 189.7036422130964d0*Lmu**2 &
	  - 51.50813149633906d0*Nf - 16.00000006398323d0*Lmu*Nf - 6.668502320598985d-9*Lmu**2*Nf, &
      -720.1068548927861d0 + 17.402351899478568d0*Lmu - 0.6180969851274909d0*Lmu**2 &
	  + 9.691339026103627d0*Nf - 2.540110346075181d-10*Lmu*Nf,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
      62244.278348304506d0 + 90595.62066115087d0*Lmu - 5307.959912413799d0*Lmu**2 &
	  - 725.7040599397907d0*Nf - 8.373370761775243d-8*Lmu*Nf + 9.071151658589846d-8*Lmu**2*Nf, &
      -73593.14246325135d0 - 13266.231345168342d0*Lmu + 562.4320363891787d0*Lmu**2 &
	  + 340.92281843959506d0*Nf + 9.023935916755379d-7*Lmu*Nf + 1.7810399835701405d-7*Lmu**2*Nf, &
      -19274.85766344364d0 + 13911.862997386164d0*Lmu - 1314.7672551699827d0*Lmu**2 &
	  + 2.9241023778057107d0*Nf + 3.285034924749133d-7*Lmu*Nf + 7.961647144010051d-8*Lmu**2*Nf,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
      -6776.731748013898d0 + 17595.79207895306d0*Lmu - 839.5223166525981d0*Lmu**2 &
	  - 73.730026699345d0*Nf + 2.823799767667272d-7*Lmu*Nf + 8.485853019485028d-8*Lmu**2*Nf, &
      -54745.95834168829d0 + 28143.140496953594d0*Lmu - 2279.2892078350114d0*Lmu**2 &
	  + 54.675824792462755d0*Nf + 9.417084318872613d-7*Lmu*Nf + 2.2242386654313947d-7*Lmu**2*Nf, &
      -17.510572839038723d0 + 0.9308218455759855d0*Lmu - 0.039924795486343194d0*Lmu**2 - 0.003977517619150036d0*Nf/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
    end if
    
   !write(*,*) 'regularPart=', regularPart/x
  end if
  end subroutine Set_Coeff_g_g
  
     !!!!!coefficient function q->qb
  subroutine Set_Coeff_q_qb(alpha,Nf,Lmu)  
  real(dp)::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_q_qb=(/&
  0d0,0d0,0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
  0d0,0d0,0d0,0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
  0d0,0d0,0d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
  0d0,0d0,0d0,&		!1 (exact), z, z^2
  0d0,0d0,0d0,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
  0d0,0d0,0d0,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
  0d0,0d0,0d0,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
  0d0,0d0,0d0/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
  if(order_global>=2) then
    Coeff_q_qb=Coeff_q_qb+alpha*alpha*(/&
      0d0, 0d0, 0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
      0d0, 0d0, 0d0, 0d0, &	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
      -76d0/9d0 + 16d0*Lmu/9d0, -32d0/9d0 + 8d0*Lmu/9d0, -4d0/3d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
      -448.1996646927744d0 - 244.21167904318d0*Lmu,&
      -5647.490349524448d0 - 1621.5639861911413d0*Lmu,&
      -257.33492226817515d0 + 202.68655426265354d0*Lmu,&		!1 (exact), z, z^2
      -6353.024936485396d0 - 1663.0891109716677d0*Lmu,&
      3394.3630072013384d0 + 551.6934658836014d0*Lmu,&
      84.62666948172131d0 - 66.45557368691802d0*Lmu,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
      483.57763340201967d0 + 101.50807540427206d0*Lmu,&
      -518.3011843768235d0 - 111.7376149098951d0*Lmu,&
      -0.3446236997186985d0 - 0.44982527645318066d0*Lmu,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
      -6484.394475423013d0 - 1245.2499053204756d0*Lmu,&
      3796.7810643424796d0 + 204.97152039767298d0*Lmu,&
      339.26036997443003d0 - 143.991684454412d0*Lmu,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
      -436.86534933603923d0 - 244.8429073687832d0*Lmu,&
      1330.4397468097895d0 - 317.06568392666475d0*Lmu,&
      -0.02218996661410778d0 - 0.021493311323627588d0*Lmu/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
    end if
    
  end subroutine Set_Coeff_q_qb
  
  !!!!!coefficient function q->qp
  subroutine Set_Coeff_q_qp(alpha,Nf,Lmu)  
  real(dp)::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_q_qp=(/&
  0d0,0d0,0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
  0d0,0d0,0d0,0d0,&	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
  0d0,0d0,0d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
  0d0,0d0,0d0,&		!1 (exact), z, z^2
  0d0,0d0,0d0,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
  0d0,0d0,0d0,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
  0d0,0d0,0d0,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
  0d0,0d0,0d0/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
  if(order_global>=2) then
    Coeff_q_qp=Coeff_q_qp+alpha*alpha*(/&
      0d0, 0d0, 0d0,&		!Log[1-z], log[1-z]^2, log[1-z]^3  !exact
      -592d0/81d0 - 16d0*Lmu/9d0 + 16d0*Lmu**2/9d0 + 32d0*zeta2/9d0, 32d0/9d0 - 64d0*Lmu/9d0, 64d0/9d0, 0d0, &	!1/z, log[z]/z, Log[z]^2/z, Log[z]^3/z  !exact
      -48d0 + 8d0*Lmu/3d0 + 8d0*Lmu**2/3d0, -22d0/3d0 - 8d0*Lmu, 44d0/9d0,&		!Log[z], log[z]^2, Log[z]^3 !exact
      56d0*Lmu/3d0 + 4d0*Lmu**2/3d0-28.813571629909163d0,&
      -40d0*Lmu/3d0 - 4d0*Lmu**2/3d0+206.17553030550255d0,&
      -32d0*Lmu/9d0 - 16d0*Lmu**2/9d0+76.0334481394452d0,&		!1 (exact), z, z^2
      251.93541929963473d0, 40d0*Lmu/3d0 + 8d0*Lmu**2/3d0-169.05906218222225d0,&
      64d0*Lmu/9d0-32.29013619101719d0,&		!zLog[z]/(1-z), z Log[z], z^2 Log[z]
      -10.685799944808078d0,-8d0*Lmu+3.7898590887852626d0,4.909581801691148d0,&		!z Log[z]^2/(1-z), z Log[z]^2,  z Log[z]^3
      188.58291934528984d0, -90.34188300607897d0, -1.4634823045099683d0,&		!(Log[z]/(1-z)+1)Log[1-z], Log[z]Log[1-z],  Log[z]Log[1-z],
      18.62607672661471d0, -16.127886782439663d0, 0.0009797967055855182d0/)		!(1-z)Log[z]/z, (1-z)Log[1-z], (1-z) Log[1-z]^2
    end if
    
  end subroutine Set_Coeff_q_qp
  
 subroutine CheckCoefficient(as,Nf,Lmu,z)
 real(dp)::Lmu,as,z,Nf
 real(dp), dimension(1:25)::func
 real(dp), dimension(1:2)::func1
 
  func=(/ Log(1d0-z),log(1d0-z)**2, log(1d0-z)**3,&
     1d0/z, log(z)/z, Log(z)**2/z, Log(z)**3/z, &
     log(z), log(z)**2, Log(z)**3,&
     1d0, z, z**2,&
     z*Log(z)/(1d0-z), z*Log(z), (z**2)*Log(z),&
     z*(Log(z)**2)/(1d0-z), z*Log(z)**2,  z*Log(z)**3, &
     (Log(z)/(1d0-z)+1d0)*Log(1d0-z), Log(z)*Log(1d0-z),  z*Log(z)*Log(1d0-z),&
     (1d0-z)/z*Log(1d0-z), (1d0-z)*Log(1d0-z), (1d0-z)*Log(1d0-z)**2/)
     
  func1=(/1d0/(1d0-z),Log(1d0-z)/(1d0-z)/)
 
  call Set_Coeff_q_g(as,Nf,Lmu)
  
  write(*,*) SUM(Coeff_q_g*func)
 
 end subroutine CheckCoefficient
