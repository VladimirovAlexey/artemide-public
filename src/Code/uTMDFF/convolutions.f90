!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!
!!!!!!!!!       part of uTMDFF module for artemide
!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!-------------------------------------------------------------------------------------------------------
!!!-----------------------------Routines for evaluation of TMD distributions------------------------------
!!!-------------------------the actual calculation takes a place in twist2 code---------------------------
!!!-------------------------------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Convolutions!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
  function xf(x,Q,hadron)
      real(dp) :: x,Q
      integer:: hadron
      real(dp), dimension(-5:5):: xf
      xf=xFF(x,Q,hadron)
      
  end function xf
 
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Functions for calculation of convolution!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!! We evaluate the integral II=\int_z^1 dy/y  C(z/y) fNP(z,y) d(z)/z^2
!! to make sure that its definition is close to PDF definition we change d(x)->xd(x)/x and evaluate
!! II=1/z**3 \int_z^1 dy  C(y) fNP(z,y) xd(z/y)
!! where  1/z \int_z^1 dy  C(y) fNP(z,y) xd(z/y) is calculated in the common code.
!! so in the end I devide by extra factor z^2.

!---------------------------------------------------------------------
!- This is the TMD function evaluated for all quarks simultaniously (-5..5) at x,bT,mu
!--    f =  -5, -4, -3,  -2,  -1,0,1,2,3, 4 ,5
!--      = bbar ,cbar sbar,ubar,dbar,g,d,u,s, c ,b
!!---Gluon contribution is undefined
!!---Base version: hadron=number of PDF
!---------------------------------------------------------------------
!---------------------------------------------------------------------
function uTMDFF_base5(z,bT,hadron)
    real(dp),dimension(-5:5)::uTMDFF_base5
    real(dp) :: z, bT
    integer::hadron,j

    !!! variables for restoration only
    real(dp),dimension(-5:5) :: fNP_grid,fNP_current

    !!!! if the coefficient function is not used at all. Just return fNP
    if(order_global==-50) then
        uTMDFF_base5=fNP(z,1d0,bT,hadron,lambdaNP)

    !!! in the case the greed has been calculated
    else if(gridReady.and. ANY(hadronsInGRID.eq.hadron)) then 
        uTMDFF_base5=ExtractFromGrid(z,bT,hadron)/z**2

        !!!!!!!!!!This is procedure of restoration of function from the initial grid
        !!! if fNP is z-independent then the value can be obtained by TMDFF(initial) fNP(current)/fNP(initial)
        if(.not.IsFnpZdependent) then
            fNP_grid=FNP(z,0d0,bT,hadron,lambdaNP_grid)
            fNP_current=FNP(z,0d0,bT,hadron,lambdaNP)

            do j=-5,5
                if(fNP_grid(j)==0) then
                if(fNP_current(j)/=0 .and. ((j/=0).or.(.not.withGluon))) then
                if(outputLevel>0) &
                    call Warning_Raise('error in restoration: original value is zero. TMDFF set to zero. b='//numToStr(bT),&
                        messageCounter,messageTrigger,moduleName)
                end if
                uTMDFF_base5(j)=0
                else
                uTMDFF_base5(j)=uTMDFF_base5(j)*fNP_current(j)/fNP_grid(j)
                end if
            end do

        end if

    !!!! finally just calculation
    else
    uTMDFF_base5=Common_lowScale5(z,bT,hadron)/z**2

    end if
end function uTMDFF_base5

!---------------------------------------------------------------------
!- This is the TMD function evaluated for all quarks simultaniously (-5..5) at x,bT,mu
!--    f =  -5, -4, -3,  -2,  -1,0,1,2,3, 4 ,5
!--      = bbar ,cbar sbar,ubar,dbar,g,d,u,s, c ,b
!!---Base version: hadron=number of PDF
!---------------------------------------------------------------------
!---------------------------------------------------------------------
function uTMDFF_base50(z,bT,hadron)
    real(dp),dimension(-5:5)::uTMDFF_base50
    real(dp) :: z, bT
    integer::hadron,j

    !!! variables for restoration only
    real(dp),dimension(-5:5) :: fNP_grid,fNP_current

    !!!! if the coefficient function is not used at all. Just return fNP
    if(order_global==-50) then

        uTMDFF_base50=fNP(z,1d0,bT,hadron,lambdaNP)

    !!! in the case the greed has been calculated
    else if(gridReady .and. ANY(hadronsInGRID.eq.hadron)) then 
        uTMDFF_base50=ExtractFromGrid(z,bT,hadron)/z**2

        if(.not.IsFnpZdependent) then
            fNP_grid=FNP(z,1d0,bT,hadron,lambdaNP_grid)
            fNP_current=FNP(z,1d0,bT,hadron,lambdaNP)

            do j=-5,5
                if(fNP_grid(j)==0) then
                if(uTMDFF_base50(j)/=0) then
                if(outputLevel>0) &
                    call Warning_Raise('error in restoration: original value is zero. TMDFF set to zero. b='//numToStr(bT),&
                        messageCounter,messageTrigger,moduleName)
                end if
                uTMDFF_base50(j)=0
                else
                uTMDFF_base50(j)=uTMDFF_base50(j)*fNP_current(j)/fNP_grid(j)
                end if
            end do

        end if

    !!!!! finally just calculation
    else
        uTMDFF_base50=Common_lowScale50(z,bT,hadron)/z**2
    end if
end function uTMDFF_base50

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!OUTPUT INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!---------------------------------------------------------------------
!- This is the TMD function evaluated for all quarks simultaniously (-5..5) at x,bT,mu
!--    f =  -5, -4, -3,  -2,  -1,0,1,2,3, 4 ,5
!--      = bbar ,cbar sbar,ubar,dbar,??,d,u,s, c ,b
!!---Gluon contribution is undefined
!!---Full version: hadron=number of PDF (if compositeness OFF)
!!---		   hadron=sum components (if compositeness ON)
!---------------------------------------------------------------------
!---------------------------------------------------------------------
function uTMDFF_lowScale5(x,bT,hadron)
  real(dp),dimension(-5:5)::uTMDFF_lowScale5
  real(dp) :: x, bT
  integer::hadron
  
  logical,allocatable::includeInComposition(:)
  real(dp),allocatable::compositionCoefficients(:)
  integer::j,jN
  
  if(x>1d0) then
    uTMDFF_lowScale5=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
    return
  end if
  
  if(IsComposite) then
    
    uTMDFF_lowScale5=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
    
    call GetCompositionArray(hadron,lambdaNP,includeInComposition,compositionCoefficients)
    jN=size(includeInComposition)
    
    do j=1,jN
      if(includeInComposition(j)) then
	uTMDFF_lowScale5=uTMDFF_lowScale5+compositionCoefficients(j)*uTMDFF_base5(x,bT,j)
      end if
    end do
    
    deallocate(includeInComposition,compositionCoefficients)
  
  else
  
   uTMDFF_lowScale5=uTMDFF_base5(x,bT,hadron)  
  end if
 
end function uTMDFF_lowScale5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !---------------------------------------------------------------------
!- This is the TMD function evaluated for all quarks simultaniously (-5..5) at x,bT,mu the GLUON INCLUDED
!--    f =  -5, -4, -3,  -2,  -1,0,1,2,3, 4 ,5
!--      = bbar ,cbar sbar,ubar,dbar,g,d,u,s, c ,b
!!---Full version: hadron=number of PDF (if compositeness OFF)
!!---		   hadron=sum components (if compositeness ON)
!---------------------------------------------------------------------
!---------------------------------------------------------------------
function uTMDFF_lowScale50(x,bT,hadron)
  real(dp),dimension(-5:5)::uTMDFF_lowScale50
  real(dp) :: x, bT
  integer::hadron
  
  logical,allocatable::includeInComposition(:)
  real(dp),allocatable::compositionCoefficients(:)
  integer::j,jN
  
  if(x>1d0) then
    uTMDFF_lowScale50=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
    return
  end if
  
  if(IsComposite) then
    
    uTMDFF_lowScale50=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
    
    call GetCompositionArray(hadron,lambdaNP,includeInComposition,compositionCoefficients)
    jN=size(includeInComposition)
    
    do j=1,jN
      if(includeInComposition(j)) then
	uTMDFF_lowScale50=uTMDFF_lowScale50+compositionCoefficients(j)*uTMDFF_base50(x,bT,j)
      end if
    end do
    
    deallocate(includeInComposition,compositionCoefficients)
  
  else
  
   uTMDFF_lowScale50=uTMDFF_base50(x,bT,hadron)  
  end if
end function uTMDFF_lowScale50

