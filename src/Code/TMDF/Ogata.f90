 !!!Prepare tables for Ogata quadrature with given h
 subroutine PrepareTables()
  integer::i,k,j
  real(dp)::hS!=h*hSegmentationWeight
  real(dp)::xi,qqq
   
  do j=1,hSegmentationNumber
  do k=0,3
  do i=1,Nmax
    
    hS=hOGATA*hSegmentationWeight(j)    
    xi=JZero(k,i)
    
!     ww(j,k,i)=BESSEL_JN(k,bb(j,k,i))/xi/(BESSEL_JN(k+1,xi)**2)&
! 	    *(pi*xi*hS*Cosh(xi*hS)+Sinh(pi*Sinh(xi*hS)))/(1d0+Cosh(pi*Sinh(xi*hS)))
    
    !!! if we too far away in xI*hS, the double exponential grow rapidly.
    !!! and for >6, it generates term 10^{300} and exceed the presision

    if(xi*hS>6.d0) then
        bb(j,k,i)=xi*Tanh(piHalf*Sinh(xi*hS))
        ww(j,k,i)=BESSEL_JN(k,bb(j,k,i))/xi/(BESSEL_JN(k+1,xi)**2)
        
    else
        bb(j,k,i)=xi*Tanh(piHalf*Sinh(xi*hS))
        ww(j,k,i)=BESSEL_JN(k,bb(j,k,i))/xi/(BESSEL_JN(k+1,xi)**2)&
        *(pi*xi*hS*Cosh(xi*hS)/(2d0*Cosh(piHalf*Sinh(xi*hS))**2)+Tanh(piHalf*Sinh(xi*hS)))
    end if

  end do
  end do
  end do 
 end subroutine PrepareTables
 
 !!!This is the defining module function
 !!! It evaluates the integral 
 !!!  int_0^infty   b^(n+1) db/2  Jn(b qT) zff F1 F2
 !!!
 function TMDF_F(Q2,qT,x1,x2,mu,zeta1,zeta2,process)
  real(dp)::TMDF_F
  real(dp)::qT,x1,x2,mu,zeta1,zeta2,Q2
  integer::process
  real(dp)::integral,eps,delta
  real(dp)::v1,v2,v3,v4
  integer::k,n,j,Nsegment
    
  CallCounter=CallCounter+1
  integral=0d0
  
  if(qT<0.0000001d0 .or. x1>=1d0 .or. x2>=1d0) then  
  integral=0d0
  else if(TMDF_IsconvergenceLost()) then	
  !!!in the case of lost convergence we return huge number (divergent xSec)
	TMDF_F=1d10		
  else
  
  v1=1d0
  v2=1d0
  v3=1d0
  v4=1d0
  
  !!Here we set the order of Bessel
  if(process<10000) then
  n=0
  else if(process<20000) then
  n=1
  else if(process<30000) then
  n=2
  else
  n=3
  end if
  
  !!! define segment of qT
  do j=1,hSegmentationNumber
    if(qT<qTSegmentationBoundary(j)) exit
  end do
  if(j>hSegmentationNumber) then
    Nsegment=hSegmentationNumber
  else
    Nsegment=j
  end if
  
  !!! sum over OGATA nodes
  do k=1,Nmax!!! maximum of number of bessel roots preevaluated in the head
    eps=ww(Nsegment,n,k)*(bb(Nsegment,n,k)**(n+1))*Integrand(Q2,bb(Nsegment,n,k)/qT,x1,x2,mu,zeta1,zeta2,process)
    
    v4=v3
    v3=v2
    v2=v1
    v1=ABS(eps)
    
    delta=(v1+v2+v3+v4)
    integral=integral+eps
    
    !!! here we check that residual term is smaller than already collected integral
    !!! also checking the zerothness of the integral. If already collected integral is null it is null
    !!! Here is potential bug. If the first 10 points give zero (whereas some later points do not), the integral will be zero
    if((delta<tolerance*abs(integral) .or. abs(integral)<1d-32) .and. k>=10) exit
  end do
  if(k>=Nmax) then	
    if(outputlevel>0) WRITE(*,*) WarningString('OGATA quadrature diverge. TMD decaing too slow? ',moduleName)
      if(outputlevel>1) then
      write(*,*) 'Information over the last call ----------'
      write(*,*) 'bt/qT= ',bb(Nsegment,n,Nmax)/qT, 'qT=',qT, '| segmentation zone=',Nsegment,&
	      ' ogata h=',hOGATA*hSegmentationWeight(Nsegment)
      write(*,*) 'W=',Integrand(Q2,bb(Nsegment,n,Nmax)/qT,x1,x2,mu,zeta1,zeta2,process), 'eps/integral =', eps/integral
      write(*,*) 'residual term=',delta, '>',tolerance
      write(*,*) '(x1,x2)=(',x1,',',x2,')'
      write(*,*) 'process =',process,' it is ',CallCounter,'call.'
      write(*,*) '------------------------------------------'
      end if
    call TMDF_convergenceISlost()
  end if
  
  if(k>MaxCounter) MaxCounter=k-1
!   write(*,*) 'Integral=',integral
  TMDF_F=integral/(qT**(n+2)) 
  end if
  !write(*,*) 'Last call: ',k
  
!    write(*,'("{",F6.2,",",F18.16,"},")') qT,x1*x2*TMDF_F
 end function TMDF_F
