!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!
!!!!!!!!!       part of TMDR module for artemide
!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!-------------------------------------------------------------------------------------------------------
!!!-----------------------------The expression for exact zeta-line----------------------------------------
!!!-------------------------------------------------------------------------------------------------------


!!!!!!!!!! exact value of zeta-line at given b,mu expanded over as
!!!!!!!!!! rad = dNP
function zetaSL(mu,rad,f)
  real(dp)::zetaSL
  real(dp),intent(in)::rad,mu
  integer,intent(in)::f
  real(dp)::alpha,GD
  
  if (orderZETA<0) then
    zetaSL=1d0
    return
  end if
  
  alpha=As(mu)
  !!dd=DNP(mu,b,f)
  if(f==0) then
    GD=valueOfGD_type4_G(rad,alpha,mu)
  else
    GD=valueOfGD_type4_Q(rad,alpha,mu)
  end if
  
  !write(*,*) '..........',b,rad,alpha,GD
  
  zetaSL=(mu**2)*exp(-GD)  
 
end function zetaSL
 
!!! expression for G, at given D,alpha, and mu
!!! evaluated for type4 evolution for Quark
 function valueOfGD_type4_Q(dd,alpha,mu)
    real(dp)::valueOfGD_type4_Q,dd,alpha,mu
    integer::Nf
    real(dp)::p,pFF,z1,zm1,zm2,zm3,val
    
    if (orderZETA<0) then
        valueOfGD_type4_Q=0d0
        return
    end if
    
    Nf=ActiveNf(mu)
    pFF=pFACTOR_Q(Nf)
    p=dd*pFF
    
    z1=zFUNC(p,1)
    
    val=z1*OMEGA_q(0,1,Nf)
        
    if(orderZETA>=1) then
        val=val+alpha*(z1*OMEGA_q(1,1,Nf)+OMEGA_q(1,2,Nf)+p*OMEGA_q(1,3,Nf))
        if(orderZETA>=2) then
            zm1=zFUNC(p,-1)
            val=val+(alpha**2)*(z1*OMEGA_q(2,1,Nf)+zm1*OMEGA_q(2,2,Nf)+OMEGA_q(2,3,Nf))
            if(orderZETA>=3) then

                !!! the problem is that the coefficients of z[-1], z[-2], are negative
                !!! therefore at large D the expression for R, diverges exponentially.
                !!! however, z[] \sim as, therefore, these terms contribute to as^4.
                !!! here I kill them. Istead I add their expansion up to p
                !val=val+(alpha**3)*(OMEGA_q(3,4,Nf)+p*OMEGA_q(3,5,Nf))
                !!! killed part
                    !!!+z1*OMEGA_q(3,1,Nf)+zm1*OMEGA_q(3,2,Nf)+zm2*OMEGA_q(3,3,Nf)


                !!! full expression
                zm2=zFUNC(p,-2)
                val=val+(alpha**3)*(OMEGA_q(3,4,Nf)+z1*OMEGA_q(3,1,Nf)+zm1*OMEGA_q(3,2,Nf)+zm2*OMEGA_q(3,3,Nf))

            if(orderZETA>=4) then
                !!! here is the same problem. So I just kill z[-1],z[-2],z[-3]
                !!! since they formally gives as^5
                !!! I do not add anything -- no need
                !val=val+(alpha**4)*(OMEGA_q(4,5,Nf))
                    !!! killed part
                    !!!+z1*OMEGA_q(4,1,Nf)+zm1*OMEGA_q(4,2,Nf)+zm2*OMEGA_q(4,3,Nf)+zm3*OMEGA_q(4,4,Nf)

                !!! full expression
                zm3=zFUNC(p,-3)
                val=val+(alpha**4)*(OMEGA_q(4,5,Nf)+z1*OMEGA_q(4,1,Nf) &
                +zm1*OMEGA_q(4,2,Nf)+zm2*OMEGA_q(4,3,Nf)+zm3*OMEGA_q(4,4,Nf))

            end if
            end if
        end if
    end if
    
    valueOfGD_type4_Q=val/alpha
        
end function valueOfGD_type4_Q

!!! expression for G*D, at given D,alpha, and mu
!!! evaluated for type4 evolution for Quark
 function valueOfGD_type4_G(dd,alpha,mu)
    real(dp)::valueOfGD_type4_G,dd,alpha,mu
    integer::Nf
    real(dp)::p,pFF,z1,zm1,zm2,val
    
    if (orderZETA<0) then
        valueOfGD_type4_G=0d0
        return
    end if
    
    Nf=ActiveNf(mu)
    pFF=pFACTOR_G(Nf)
    p=dd*pFF
    
    z1=zFUNC(p,1)
    
    val=z1*OMEGA_g(0,1,Nf)
        
    if(orderZETA>=1) then
        val=val+alpha*(z1*OMEGA_g(1,1,Nf)+OMEGA_g(1,2,Nf)+p*OMEGA_g(1,3,Nf))
        if(orderZETA>=2) then
            zm1=zFUNC(p,-1)
            val=val+(alpha**2)*(z1*OMEGA_g(2,1,Nf)+zm1*OMEGA_g(2,2,Nf)+OMEGA_g(2,3,Nf))
            if(orderZETA>=3) then
!                 zm2=zFUNC(p,-2)
            !!! the problem is that the coefficients of z[-1], z[-2], are negative
            !!! therefore at large D the expression for R, diverges exponentially.
            !!! however, z[] \sim as, therefore, these terms contribute to as^4.
            !!! here I kill them. Istead I add their expansion up to p
                val=val+(alpha**3)*(OMEGA_g(3,4,Nf)+p*OMEGA_g(3,5,Nf))
                    !!! killed part
                    !!!+z1*OMEGA_g(3,1,Nf)+zm1*OMEGA_g(3,2,Nf)+zm2*OMEGA_g(3,3,Nf)
            if(orderZETA>=4) then
            !!! here is the same problem. So I just kill z[-1],z[-2],z[-3]
            !!! since they formally gives as^5
            !!! I do not add anything -- no need
                val=val+(alpha**4)*(OMEGA_g(4,5,Nf))
                    !!! killed part
                    !!!+z1*OMEGA_g(4,1,Nf)+zm1*OMEGA_g(4,2,Nf)+zm2*OMEGA_g(4,3,Nf)+zm3*OMEGA_g(4,4,Nf)
            end if
            end if
        end if
    end if
    
    valueOfGD_type4_G=val/alpha
        
end function valueOfGD_type4_G

!!! this is z=(exp(- np)-1+np)/p
pure function zFUNC(p,n)
    real(dp)::zFUNC
    real(dp),intent(in)::p
    integer,intent(in)::n
    real(dp):: freezeP
    
    !!!! in the case p~0, the computation requires cancelation of zero, I expand in a series
    if(Abs(p)>0.00001d0) then

        !!!! In the case n<0, and p very large, the exponent blows-up
        !!!! this is due to the nonperurbative effects. Because it happen at b->infty
        !!!! To cure this problem, I freeze the expression for n<0, and large p.
        !!!! For typical CS-kernels, this freze starts at p~2 or b>~5-6 GeV [at 2GeV, and b~2 at 100 GeV], and is suppressed as as^3, or as^4.
        !!!! Therefore, it should not be important,
        if(n<0 .and. p>2.) then
            !!! the tail is 3a/2+ p/(1-p*(3/a**2)), where a is the point at which the effect turns on
            !!! in this form, the function is continious with its first derivative
            freezeP=3d0+p/(1d0-0.75d0*p**2)
            zFUNC=(real(n,dp)*freezeP-1d0+Exp(-real(n,dp)*freezeP))/freezeP
        else

            !! normal case
            zFUNC=(real(n,dp)*p-1d0+Exp(-real(n,dp)*p))/p
        end if
    else

        !!!! in the case p~0, the computation requires cancelation of zero, I expand in a series
        zFUNC=real(n**2,dp)*p/2d0-real(n**3,dp)*(p**2)/6d0
    end if



end function zFUNC 
 
