!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!
!!!!!!!!!       part of TMD_AD module for artemide
!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!-------------------------------------------------------------------------------------------------------
!!!-----------------------------Anomalous dimensions at scale mu------------------------------------------
!!!-------------------------------------------------------------------------------------------------------
  
! The TMD rapidity anomalous dimension D, pertrubative expression
function Dpert(mu,bT,f)
    real(dp),intent(in):: bT,mu
    integer,intent(in):: f
    real(dp)::Dpert
    integer:: Nf,n,k
    real(dp)::LL,astrong,inter
    
    LL=2_dp*LOG(bt*mu*C0_inv_const)
    astrong=As(mu)
    Dpert=0_dp
    
    Nf=ActiveNf(mu)
    
    if(f==0) then!gluon
      do n=1,orderD
        inter=0_dp
        do k=0,orderD
         inter=inter+dnk_G(n,k,Nf)*(LL**k)
        end do
        Dpert=Dpert+inter*(astrong**n)
      end do
     else!quark
      do n=1,orderD
        inter=0_dp
        do k=0,orderD
         inter=inter+dnk_Q(n,k,Nf)*(LL**k)
        end do
        Dpert=Dpert+inter*(astrong**n)	
      end do
    end if

end function Dpert
  
  ! The TMD UV anomalous dimension gammaV
function gammaV(mu,f)
    real(dp),intent(in)::mu
    integer,intent(in)::f
    real(dp):: gammaV
    integer:: Nf,i    
    real(dp) ::astrong
    
    astrong=As(mu)
    
    Nf=ActiveNf(mu)
   
    gammaV=0_dp
    if(f==0) then   !!gluon case
     do i=1,orderV
        gammaV=gammaV+gammaV_G(i,Nf)*(astrong**i)
     end do
    else  !! quark case
     do i=1,orderV
        gammaV=gammaV+gammaV_Q(i,Nf)*(astrong**i)
     end do
    end if

end function gammaV

  ! The Gamma Cusp anomalous dimension gammaV
function gammaCUSP(mu,f)
    real(dp),intent(in)::mu
    integer,intent(in)::f
    real(dp) :: gammaCUSP
    integer:: Nf,i
    real(dp) ::astrong
    
    astrong=As(mu)
    
    Nf=activeNf(mu)
    
    gammaCUSP=0_dp
    if(f==0) then   !!gluon case
     do i=1,orderCusp
        gammaCUSP=gammaCUSP+GammaCusp_G(i-1,Nf)*(astrong**i)
     end do
    else  !! quark case
     do i=1,orderCusp
        gammaCUSP=gammaCUSP+GammaCusp_Q(i-1,Nf)*(astrong**i)
     end do
    end if
   
end function gammaCUSP

!!the resummed version of rapidity anomalous dimension
function Dresum(mu,bT,f)
    real(dp),intent(in)::mu,bT
    integer,intent(in)::f
    real(dp)::Dresum
    integer::Nf
    integer::n,k,l
    real(dp):: X,alpha,lX,LL,commulant,inter
    
    alpha=As(mu)
    Nf=ActiveNf(mu)
    LL=2_dp*LOG(bT*mu*C0_inv_const)
    
    X=betaQCD(0,Nf)*alpha*LL
    
    
    !!! in many models D evaluated at X=0, to speed up these computations I explicitely state it
    if(abs(X)<1d-7) then
    Dresum=0_dp
    if(f==0) then !! gluon case
        do n=2,orderDresum
            Dresum=Dresum+alpha**n*dnk_G(n,0,Nf)
        end do
    else !! quark case
        do n=2,orderDresum
            Dresum=Dresum+alpha**n*dnk_Q(n,0,Nf)
        end do
    end if
    
    !!! complete case
    else    
    
    lX=Log(1_dp-X)
    
    if(f==0) then !! gluon case
        commulant=lX
        
        do n=1,orderDresum
        inter=0_dp
            do k=0,n
            do l=0,n
                inter=inter+(X**k)*(lX**l)*dnkl_G(n,k,l,Nf)
            end do
            end do
        commulant=commulant+((alpha/(1_dp-X))**n)*inter
        end do
        
        Dresum=-GammaCusp_G(0,Nf)/betaQCD(0,Nf)*commulant
    else !!! quark case
        commulant=lX
        do n=1,orderDresum
        inter=0_dp
            do k=0,n
            do l=0,n
                inter=inter+(X**k)*(lX**l)*dnkl_Q(n,k,l,Nf)
            end do
            end do
        commulant=commulant+((alpha/(1_dp-X))**n)*inter
        end do
        
        Dresum=-GammaCusp_Q(0,Nf)/betaQCD(0,Nf)/2_dp*commulant  
    end if
    
    if(ISNAN(Dresum)) then
        write(*,*) ErrorString('Dresum is NaN.',moduleName)
        write(*,*) 'At mu=',mu,'b=',bT,'Lmu=',2_dp*Log(mu*bT*C0_inv_const), 'X=',X,'log(1-x)=',lX
        write(*,*) 'Evaluation STOP'
        stop
    end if
    
    end if
end function Dresum

!-------------------zeta-lines -------------------------------------
!the exact zeta-line is in TMDR (because it depends on DNP)

!! the value of zeta_mu in the pertrubation theory with ri=0
function zetaMUpert(mu,bt,f)
  real(dp),intent(in)::mu,bT
  integer,intent(in)::f
  real(dp):: zetaMUpert
  integer::Nf,n,k
  real(dp)::alpha,LL,val,iter
  
  if(orderZETA>0) then !!! NLO,NNLO,...
  
    LL=2_dp*LOG(bt*mu*C0_inv_const)
    alpha=As(mu)
    Nf=ActiveNf(mu)    
    
    !!!!! Important!!!!
    !! the perturbative value for zeta, must be taken 1-order higher
    !! because there double logarithms ~~beta0 L
    !! For 4-loop it requires 5-loop rad... The value which contains it is set to zero v(4,0)=0.
    if(f==0) then !!! gluon
        val=vnk_g(0,0,Nf)
        
        do n=1,orderZETA-1
            iter = 0_dp
            do k=0,n+1           
                iter=iter+(LL**k)*vnk_g(n,k,Nf)            
            end do
            val=val+(alpha**n)*iter
        end do   
        
    else !!!! quark
        val=vnk_q(0,0,Nf)
        
        do n=1,orderZETA
            iter = 0_dp
            do k=0,n+1
                iter=iter+(LL**k)*vnk_q(n,k,Nf)
            end do            
            val=val+(alpha**n)*iter
        end do
    end if
  
    zetaMUpert=mu*C0_const/bT*EXP(-val)
  
  else if (orderZETA==0) then   !!! LO 
    zetaMUpert=mu*C0_const/bT
    
  else !!!! just in case
    zetaMUpert=1_dp
  end if
  
end function zetaMUpert
 
