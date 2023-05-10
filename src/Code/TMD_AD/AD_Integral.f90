!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!
!!!!!!!!!       part of TMD_AD module for artemide
!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!-------------------------------------------------------------------------------------------------------
!!!-----------------------------Integrals with anomalous dimensions---------------------------------------
!!!-------------------------------------------------------------------------------------------------------


!!!! Defines the coefficients for the decomposition for the integral Gamma/(-2 beta) over the roots of beta-function
!!!! The zeroth terms is G0/(-2b0) (common)
!!!! k=number of root. The expression is given by Lagrange decomposition
subroutine SetIntegralCoefficeintsGAMMA_Q()
    integer::n,i,j,k,l
    real(dp)::G1,G2,G3,G4,BB
    COMPLEX(dp)::ps
    
    GammaIntegral_Q_internal=0d0*GammaIntegral_Q_internal
    
    do n=NfMIN,NfMAX
        G1=GammaCusp_q(1,n)/GammaCusp_q(0,n)
        G2=GammaCusp_q(2,n)/GammaCusp_q(0,n)
        G3=GammaCusp_q(3,n)/GammaCusp_q(0,n)
        G4=GammaCusp_q(4,n)/GammaCusp_q(0,n)
        
        !!! loops
        do i=1,4
            !!! special case: common prefactor = Gamma0/(-2beta0)
            GammaIntegral_Q_internal(i,0,n)=GammaCusp_q(0,n)/(-2d0*betaQCD(0,n))
            
            !!! Bi
            BB=betaQCD(i,n)/betaQCD(0,n)            
            !!! roots
            do j=1,i
                !!! 0;th-root=0 
                ps=betaRoots_internal(i,j,n)
                !!! rest roots 
                do k=1,i
                if(k.ne.j) ps=ps*(betaRoots_internal(i,j,n)-betaRoots_internal(i,k,n))
                end do
                
                !!! numerator is Gamma at root
                GammaIntegral_Q_internal(i,j,n)=1d0
                if(i>0) GammaIntegral_Q_internal(i,j,n)=GammaIntegral_Q_internal(i,j,n)+G1*betaRoots_internal(i,j,n)
                if(i>1) GammaIntegral_Q_internal(i,j,n)=GammaIntegral_Q_internal(i,j,n)+G2*betaRoots_internal(i,j,n)**2
                if(i>2) GammaIntegral_Q_internal(i,j,n)=GammaIntegral_Q_internal(i,j,n)+G3*betaRoots_internal(i,j,n)**3
                if(i>3) GammaIntegral_Q_internal(i,j,n)=GammaIntegral_Q_internal(i,j,n)+G4*betaRoots_internal(i,j,n)**4
                
                GammaIntegral_Q_internal(i,j,n)=GammaIntegral_Q_internal(i,j,n)/ps/BB
            end do
        end do

    end do

end subroutine SetIntegralCoefficeintsGAMMA_Q

!!!! Defines the coefficients for the decomposition for the integral Gamma/(-2 beta) over the roots of beta-function
!!!! The zeroth terms is G0/(-2b0) (common)
!!!! k=number of root. The expression is given by Lagrange decomposition
subroutine SetIntegralCoefficeintsGAMMA_G()
    integer::n,i,j,k,l
    real(dp)::G1,G2,G3,G4,BB
    COMPLEX(dp)::ps
    
    GammaIntegral_G_internal=0d0*GammaIntegral_G_internal
    
    do n=NfMIN,NfMAX
        G1=GammaCusp_g(1,n)/GammaCusp_g(0,n)
        G2=GammaCusp_g(2,n)/GammaCusp_g(0,n)
        G3=GammaCusp_g(3,n)/GammaCusp_g(0,n)
        G4=GammaCusp_g(4,n)/GammaCusp_g(0,n)
        
        !!! loops
        do i=1,4
            !!! special case: common prefactor = Gamma0/(-2beta0)
            GammaIntegral_G_internal(i,0,n)=GammaCusp_g(0,n)/(-2d0*betaQCD(0,n))
            
            !!! Bi
            BB=betaQCD(i,n)/betaQCD(0,n)            
            !!! roots
            do j=1,i
                !!! 0;th-root=0 
                ps=betaRoots_internal(i,j,n)
                !!! rest roots 
                do k=1,i
                if(k.ne.j) ps=ps*(betaRoots_internal(i,j,n)-betaRoots_internal(i,k,n))
                end do
                
                !!! numerator is Gamma at root
                GammaIntegral_G_internal(i,j,n)=1d0
                if(i>0) GammaIntegral_G_internal(i,j,n)=GammaIntegral_G_internal(i,j,n)+G1*betaRoots_internal(i,j,n)
                if(i>1) GammaIntegral_G_internal(i,j,n)=GammaIntegral_G_internal(i,j,n)+G2*betaRoots_internal(i,j,n)**2
                if(i>2) GammaIntegral_G_internal(i,j,n)=GammaIntegral_G_internal(i,j,n)+G3*betaRoots_internal(i,j,n)**3
                if(i>3) GammaIntegral_G_internal(i,j,n)=GammaIntegral_G_internal(i,j,n)+G4*betaRoots_internal(i,j,n)**4
                
                GammaIntegral_G_internal(i,j,n)=GammaIntegral_G_internal(i,j,n)/ps/BB
            end do
        end do

    end do

end subroutine SetIntegralCoefficeintsGAMMA_G

!!!! Compute the integral
!!!! Gamma(a)/(-2beta(a)) for [a0,a1]
!!!! loop=loop order (0=1loop), Nf=number of flavors
!!!! Uses exact analytic expression with decomposition over roots of beta-function.
function GammaIntegral_Q(a0,a1,loop,nf)
    integer, intent(in)::loop,nf
    real(dp), intent(in)::a1,a0
    real(dp)::GammaIntegral_Q
    
    COMPLEX(dp):: inter
    integer::i
       
    !!! common term from the 0th-root
    inter=LOG(a1/a0)
    !!! sum over roots
    do i=1,loop
        inter=inter+GammaIntegral_Q_internal(loop,i,nf)*LOG((a1-betaRoots_internal(loop,i,nf))/(a0-betaRoots_internal(loop,i,nf)))        
    end do    
    !!!! common factor (I use index=1, it is not defined for loop=0)
    GammaIntegral_Q=GammaIntegral_Q_internal(1,0,nf)*dreal(inter)
end function GammaIntegral_Q

!!!! Compute the integral
!!!! Gamma(a)/(-2beta(a)) for [a0,a1]
!!!! loop=loop order (0=1loop), Nf=number of flavors
!!!! Uses exact analytic expression with decomposition over roots of beta-function.
function GammaIntegral_G(a0,a1,loop,nf)
    integer, intent(in)::loop,nf
    real(dp), intent(in)::a1,a0
    real(dp)::GammaIntegral_G
    
    COMPLEX(dp):: inter
    integer::i
       
    !!! common term from the 0th-root
    inter=LOG(a1/a0)
    !!! sum over roots
    do i=1,loop
        inter=inter+GammaIntegral_G_internal(loop,i,nf)*LOG((a1-betaRoots_internal(loop,i,nf))/(a0-betaRoots_internal(loop,i,nf)))       
    end do
    GammaIntegral_G=GammaIntegral_G_internal(1,0,nf)*dreal(inter)
end function GammaIntegral_G


!!!! gives the additive factor for the RAD evolution from (mu0) to (mu1)
!!!! It is equal to: Gamma(a)/2 dmu2/mu2 for mu in [mu0,mu1]
!!!! loop= orderCusp, f=flavor
function RADevolution(mu0,mu1,f)
    real(dp),intent(in)::mu0,mu1
    integer,intent(in)::f
    real(dp):: RADevolution
    
    real(dp)::a0,a1,inter,ar1,ar2
    logical::naturalOrder !!! mu0<mu1
    integer::n0,n1
    
    if(mu0<mu1) then
        a0=as(mu0)
        a1=as(mu1)
        naturalOrder=.true.
        n0=activeNf(mu0)
        n1=activeNf(mu1)
    else !!! order of limits inverted
        a0=as(mu1)
        a1=as(mu0)
        naturalOrder=.false.
        n0=activeNf(mu1)
        n1=activeNf(mu0)
    end if
    
    if(n0==n1) then
        if(f==0) then
            inter=GammaIntegral_G(a0,a1,orderCusp,n0)
        else
            inter=GammaIntegral_Q(a0,a1,orderCusp,n0)
        end if
    else if (n0==3 .and. n1==4) then
        ar1=as(mCHARM)
        if(f==0) then
            inter=GammaIntegral_G(a0,ar1,orderCusp,3)+GammaIntegral_G(ar1,a1,orderCusp,4)
        else
            inter=GammaIntegral_Q(a0,ar1,orderCusp,3)+GammaIntegral_Q(ar1,a1,orderCusp,4)
        end if
    else if (n0==4 .and. n1==5) then
        ar1=as(mBOTTOM)
        if(f==0) then
            inter=GammaIntegral_G(a0,ar1,orderCusp,4)+GammaIntegral_G(ar1,a1,orderCusp,5)
        else
            inter=GammaIntegral_Q(a0,ar1,orderCusp,4)+GammaIntegral_Q(ar1,a1,orderCusp,5)
        end if
    else if (n0==3 .and. n1==5) then
        ar1=as(mCHARM)
        ar2=as(mBOTTOM)
        if(f==0) then
            inter=GammaIntegral_G(a0,ar1,orderCusp,3)+GammaIntegral_G(ar1,ar2,orderCusp,4)+GammaIntegral_G(ar2,a1,orderCusp,5)
        else
            inter=GammaIntegral_Q(a0,ar1,orderCusp,3)+GammaIntegral_Q(ar1,ar2,orderCusp,4)+GammaIntegral_Q(ar2,a1,orderCusp,5)
        end if
    else
        write(*,*) WarningString(' RADevolution finds impossible combination of thresholds.',moduleName)
        write(*,*) "Numbers for Nf", n0, " to ",n1, "(compute with ",n0,")"
        inter=GammaIntegral_G(a0,a1,orderCusp,n0)
    end if
    
    if(naturalOrder) then
        RADevolution=inter
    else
        RADevolution=-inter
    end if

end function RADevolution
