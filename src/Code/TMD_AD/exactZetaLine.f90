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
  
  zetaSL=(mu**2)*exp(-GD/rad)
 
end function zetaSL
 
!!! expression for G*D, at given D,alpha, and mu
!!! evaluated for type4 evolution for Quark
!!!!
!!!! At NNNLO I use NNLO solution. The reason is the negative values of constants at large-d. It is unphysical
!!!! To be solved...
!!!!
 function valueOfGD_type4_Q(dd,alpha,mu)
    real(dp)::valueOfGD_type4_Q,dd,alpha,mu
    integer::Nf
    real(dp)::p,pFF,z1,zm1,zm2,val
    
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
!             if(orderZETA>=3) then
!                 zm2=zFUNC(p,-2)
!                 val=val+(alpha**3)*(z1*OMEGA_q(3,1,Nf)+zm1*OMEGA_q(3,2,Nf)+zm2*OMEGA_q(3,3,Nf)+OMEGA_q(3,4,Nf))
!             end if
        end if
    end if
    
    valueOfGD_type4_Q=dd/alpha*val
        
end function valueOfGD_type4_Q

!!! expression for G*D, at given D,alpha, and mu
!!! evaluated for type4 evolution for Quark
!!!!
!!!! At NNNLO I use NNLO solution. The reason is the negative values of constants at large-d. It is unphysical
!!!! To be solved...
!!!!
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
!             if(orderZETA>=3) then
!                 zm2=zFUNC(p,-2)
!                 val=val+(alpha**3)*(z1*OMEGA_g(3,1,Nf)+zm1*OMEGA_g(3,2,Nf)+zm2*OMEGA_g(3,3,Nf)+OMEGA_g(3,4,Nf))
!             end if
        end if
    end if
    
    valueOfGD_type4_G=dd/alpha*val
        
end function valueOfGD_type4_G

!!! this is z=(exp(- np)-1+np)/p
pure function zFUNC(p,n)
    real(dp)::zFUNC
    real(dp),intent(in)::p
    integer,intent(in)::n
    
    if(Abs(p)>0.00001d0) then
        zFUNC=(real(n,dp)*p-1d0+Exp(-real(n,dp)*p))/p
    else
        zFUNC=real(n**2,dp)*p/2d0-real(n**3,dp)*(p**2)/6d0
    end if

end function zFUNC 
 
