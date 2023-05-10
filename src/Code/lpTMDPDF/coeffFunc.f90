!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!
!!!!!!!!!       part of lpTMDPDF module for artemide
!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!-------------------------------------------------------------------------------------------------------
!!!-----------------------------Routines for the mathing coefficient--------------------------------------
!!!--------------The order is accumulative pertrubative order of coefficient =0,1,2 (LO,NLO,NNLO)---------
!!!-------------------------------------------------------------------------------------------------------
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!COEFFICIENT FUNCTIONS!!
  
  !!! the function which contains the functions of parameterizations
  function parametrizationString(z)
  real(dp)::z,lz,llz,zz
  real(dp),dimension(1:parametrizationLength)::parametrizationString
      zz=1d0-z
      lz=Log(z)
      llz=Log(zz)
      parametrizationString=(/&
	    1d0/z,lz/z,&  
	    lz,lz**2,&
	    1d0,zz,z*zz,& 
	    zz*llz/z,&
	    z*lz,z**2*lz,&
	    zz*llz,zz**2*llz,&
	    lz**2*llz,lz*llz**2/)
  
  end function parametrizationString
  
    !!! the function which contains 
    !!! int_z^1 parameterization at values of z -> 1
    !!! it is used to estimate integration error at z~1
  function parametrizationStringAt1(z)
  real(dp)::z
  real(dp),dimension(1:parametrizationLength)::parametrizationStringAt1
  
  parametrizationStringAt1=(/1d0-z, 0d0,0d0,0d0,1d0-z,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
  
  end function parametrizationStringAt1

  
    !!!!Each coefficient is split to delta, sing x->1, regular
    
  !!!!!coefficient function q<-q delta-part
  !!!! NO QUARK HERE!
  function C_q_q_delta(alpha,Nf,Lmu)
  real(dp)::C_q_q_delta,Nf,alpha,Lmu
  
    C_q_q_delta=0d0
  end function C_q_q_delta
  
  !!!!!coefficient function g<-g delta-part
  !!!! NO DELTA-function!
  function C_g_g_delta(alpha,Nf,Lmu)
  real(dp)::C_g_g_delta,Nf,alpha,Lmu
  
  C_g_g_delta=0d0
  end function C_g_g_delta
  
  !!!!!coefficient function q<-q singular-part  (1/(1-x)_+,(Log(1-x)/(1-x))_+)
   !!!! NO QUARK HERE!
  subroutine Set_CoeffSing1_q_q(alpha,Nf,Lmu)
  real(dp)::Nf,alpha,LLambda,Lmu
  
  CoeffSing1_q_q=(/0d0,0d0,0d0/)
  
  end subroutine Set_CoeffSing1_q_q
  
  !!!!!coefficient function g<-g singular-part  (1/(1-x)_+,(Log(1-x)/(1-x))_+)
  !!!! NO SINGULAR PART HERE!
  subroutine Set_CoeffSing1_g_g(alpha,Nf,Lmu)
  real(dp)::Nf,alpha,Lmu
  
  CoeffSing1_g_g=(/0d0,0d0,0d0/)
  
  end subroutine Set_CoeffSing1_g_g
  
  !!!!!coefficient function q<-q regular-part
  !!!! NO QUARK HERE!
  subroutine Set_Coeff_q_q(alpha,Nf,Lmu)  
  real(dp)::alpha,Nf,Lmu
  
  Coeff_q_q=(/0d0, 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
  end subroutine Set_Coeff_q_q
 
  !!!!!coefficient function q<-g regular-part  
  !!!! NO QUARK HERE!
  subroutine Set_Coeff_q_g(alpha,Nf,Lmu)  
  real(dp)::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_q_g=(/0d0, 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
  end subroutine Set_Coeff_q_g
  
    !!!!!coefficient function g<-q regular-part  
  subroutine Set_Coeff_g_q(alpha,Nf,Lmu)  
  real(dp)::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_g_q=(/0d0, 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
  
  if(order_global>=1) then
    Coeff_g_q=Coeff_g_q+alpha*(/&
      -16d0/3d0,0d0,&	!1/x, log[x]/x
      0d0,0d0,&		!log[x], log[x]^2
      16d0/3d0,0d0,0d0,&	!1 ,1-x, x(1-x)
      0d0,&		!(1-x)log[1-x]/x
      0d0,0d0,&		!xlog[x], x^2 log[x]
      0d0,0d0,&		!(1-x)Log[1-x],(1-x)^2Log[1-x]
      0d0,0d0/)		!log[1-x]Log[x]^2,log[1-x]^2log[x]
  
  !------The kernels are calculated in mathematica
    if(order_global>=2) then
     Coeff_g_q=Coeff_g_q+alpha*alpha*(/&
      352d0/9d0+128d0/27d0*Nf+Lmu*(32d0/9d0*Nf-512d0/3d0)+80d0*zeta2-128d0*zeta3,& !1/x,
      64d0*(1d0-Lmu),&	! log[x]/x
      1120d0/9d0-448d0/9d0*Lmu,&	!log[x]
      -224d0/9d0,&			!log[x]^2
      -352d0/9d0-128d0/27d0*Nf+Lmu*(-32d0/9d0*Nf+512d0/3d0)-80d0*zeta2+128d0*zeta3,&	!1
      96.93661014992554d0 - 80d0*Lmu/9d0,&	!1-x
      -12.03616313305288d0, &	!x(1-x)
      -70.94828373998075d0+256d0*Lmu/9d0 + 64d0*Nf/9d0,&		!(1-x)log[1-x]/x
      66.13401288059151d0, -4.368597177618905d0,&		!xlog[x], x^2 log[x]
      29.77672764648373d0, -7.277284712894577d0,&		!(1-x)Log[1-x],(1-x)^2Log[1-x]
      -0.21013073998329043d0, 17.922444586392995d0/)		!log[1-x]Log[x]^2,log[1-x]^2log[x]
    end if
    
    
  end if
  end subroutine Set_Coeff_g_q
  
      !!!!!coefficient function g<-g regular-part  
  subroutine Set_Coeff_g_g(alpha,Nf,Lmu)  
  real(dp)::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_g_g=(/0d0, 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
  if(order_global>=1) then
    Coeff_g_g=Coeff_g_g+alpha*(/&
      -12d0,0d0,&	!1/x, log[x]/x
      0d0,0d0,&		!log[x], log[x]^2
      12d0,0d0,0d0,&	!1 ,1-x, x(1-x)
      0d0,&		!(1-x)log[1-x]/x
      0d0,0d0,&		!xlog[x], x^2 log[x]
      0d0,0d0,&		!(1-x)Log[1-x],(1-x)^2Log[1-x]
      0d0,0d0/)		!log[1-x]Log[x]^2,log[1-x]^2log[x]
  
  !------The kernels are calculated in mathematica
    if(order_global>=2) then
     Coeff_g_g=Coeff_g_g+alpha*alpha*(/&
      220d0+12d0*Nf+Lmu*(64d0/9d0*Nf-408d0)+180d0*zeta2-288d0*zeta3,& !1/x,
      144d0*(1d0-Lmu),&	! log[x]/x
      228d0+8d0*Nf+Lmu*(32d0/3d0*Nf-144d0),&	!log[x]
      -72d0+16d0*Nf/3d0,&			!log[x]^2
      -208d0-16d0*Nf+Lmu*(-64d0/9d0*Nf+408d0)-180d0*zeta2+288d0*zeta3,&	!1
      -1069d0/4d0+80d0/3d0*Nf+Lmu*(64d0/9d0*Nf-48d0)+288d0*zeta2,&	!1-x
      -28d0/3d0*Nf+Lmu*(-32d0/9d0*Nf+24d0)-47.26606760498457d0, &	!x(1-x)
      -0.04368111639126644d0+144d0*Lmu,&		!(1-x)log[1-x]/x
      143.74118698578818d0, -3.936195409223967d0,&		!xlog[x], x^2 log[x]
      18.440046947186577d0, -31.890749138002473d0,&		!(1-x)Log[1-x],(1-x)^2Log[1-x]
      1.080993168990866d0, -1.8412301915911546d0/)		!log[1-x]Log[x]^2,log[1-x]^2log[x]
    end if
    
  end if
  end subroutine Set_Coeff_g_g

  !!!!!coefficient function q<-qb regular-part  
  !!!! NO QUARK HERE!
  subroutine Set_Coeff_q_qb(alpha,Nf,Lmu)  
  real(dp)::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_q_qb=(/0d0, 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
  end subroutine Set_Coeff_q_qb
  
  !!!!!coefficient function q<-qp regular-part
  !!!! NO QUARK HERE!
  subroutine Set_Coeff_q_qp(alpha,Nf,Lmu)  
  real(dp)::alpha,Nf,Lmu
  
  !! the Leading order is always zero, therefore calculation should be done only for order >=1
  Coeff_q_qp=(/0d0, 0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
  end subroutine Set_Coeff_q_qp
  
  !!! This function has been used during debuging
 subroutine CheckCoefficient(as,Nf,Lmu,z)
 real(dp)::Lmu,as,z,Nf
 real(dp), dimension(1:parametrizationLength)::func
 real(dp), dimension(1:3)::func1
 
  func=parametrizationString(z)
     
  func1=(/1d0/(1d0-z),Log(1d0-z)/(1d0-z),Log(1d0-z)/(1d0-z)/)
 
 !!Q->Q
!   call Set_CoeffSing1_q_q(as,Nf,Lmu) 
!   call Set_Coeff_q_q(as,Nf,Lmu)  
!   write(*,*) SUM(Coeff_q_q*func)+SUM(CoeffSing1_q_q*func1)
  
!   !!Q->G
!   call Set_Coeff_q_g(as,Nf,Lmu)  
!   write(*,*) SUM(Coeff_q_g*func)

!   !!Q->Q'
!   call Set_Coeff_q_qp(as,Nf,Lmu)  
!   write(*,*) SUM(Coeff_q_qp*func)

  !!Q->Qbar
!   call Set_Coeff_q_qb(as,Nf,Lmu)  
!   write(*,*) SUM(Coeff_q_qb*func)
 
!  !! G->Q
!   call Set_Coeff_g_q(as,Nf,Lmu)  
!   write(*,*) SUM(Coeff_g_q*func)

!	!!G->G

  call Set_CoeffSing1_g_g(as,Nf,Lmu) 
  call Set_Coeff_g_g(as,Nf,Lmu)  
!   
!   write(*,*) Coeff_g_g
!   write(*,*) '---------'
!   write(*,*) func
  
  write(*,*) SUM(Coeff_g_g*func)+SUM(CoeffSing1_g_g*func1)
 end subroutine CheckCoefficient
