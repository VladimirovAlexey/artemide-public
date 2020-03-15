!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			Model for TMD defining scales
!			used only in the case of evolution types 1 and 2
!	
!			mu_low=C0/b+2
!			mu0=C0/b+2
!
!				A.Vladimirov (25.04.2018)
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

  
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! 
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!USER DEFINED FUNCTION!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  !!!!This function is the mu(b)
  ! it is used for solutions with variable lower evolution scale (EvolutionTYpe=1,2)
  function mu_LOW(bt)
  real(dp)::bt,mu_LOW
  
    mu_LOW=10d0!C0_const*1d0/bT+2d0
    
    if(mu_LOW>1000d0) then
      mu_LOW=1000d0
    end if
  end function mu_LOW
  
  !!!!This function is the mu0(b)
  ! it is used for solutions within improved D solution (EvolutionTYpe=1)
  function mu0(bt)
  real(dp)::bt,mu0
  
  mu0=mu_LOW(bt)
  end function mu0
