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
  real*8::bt,mu_LOW
  
  real*8, parameter :: C0_const=1.1229189671337703d0  
  
  mu_LOW=C0_const*1d0/bT+2d0!mu_OPE(bt)
  end function mu_LOW
  
  !!!!This function is the mu0(b)
  ! it is used for solutions within improved D solution (EvolutionTYpe=1)
  function mu0(bt)
  real*8::bt,mu0
  
  mu0=mu_LOW(bt)
  end function mu0
