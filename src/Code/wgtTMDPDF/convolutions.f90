!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!
!!!!!!!!!       part of wgtTMDPDF module for artemide
!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!-------------------------------------------------------------------------------------------------------
!!!-----------------------------Routines for evaluation of TMD distributions------------------------------
!!!-------------------------the actual calculation takes a place in twist2 code---------------------------
!!!-------------------------------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Convolutions!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!-------------------------------------------------
 !!!!array of x times PDF(x,Q) for hadron 'hadron'
 !!!! array is (-5:5) (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
 function xf(x,Q,hadron)
      real(dp) :: x,Q
      integer:: hadron
      real(dp), dimension(-5:5):: xf
      
      !!! helicity times x
      xf=x_hPDF(x,Q,hadron)
      !xf=(/x,x,x,x,x,x,x**2,x,x,x,x/)
      
  end function xf


!---------------------------------------------------------------------
!- This is the TMD function evaluated for all quarks simultaniously (-5..5) at x,bT,mu
!--    f =  -5, -4, -3,  -2,  -1,0,1,2,3, 4 ,5
!--      = bbar ,cbar sbar,ubar,dbar,??,d,u,s, c ,b
!!---Gluon contribution is undefined
!!---Base version: hadron=number of PDF
!---------------------------------------------------------------------
!---------------------------------------------------------------------
function wgtTMDPDF_base5(x,bT,hadron)
    real(dp),dimension(-5:5)::wgtTMDPDF_base5
    real(dp) :: x, bT
    integer::hadron

    !!! variables for restoration only
    real(dp),dimension(-5:5) :: fNP_grid,fNP_current
    integer::j

    !!!! if the coefficient function is not used at all. Just return fNP
    if(order_global<0) then

        wgtTMDPDF_base5=fNP(x,1d0,bT,hadron,lambdaNP)
        
    !!! in the case the greed has been calculated AND the hadron is in the grid
    else if(gridReady .and. ANY(hadronsInGRID.eq.hadron)) then 

        wgtTMDPDF_base5=x*ExtractFromGrid(x,bT,hadron)

        !!!!!!!!!!This is procedure of restoration of function from the initial grid
        !!! if fNP is x-independent then the value can be obtained by TMDPDF(initial) fNP(current)/fNP(initial)
        if(.not.IsFnpZdependent) then
            fNP_grid=FNP(x,0d0,bT,hadron,lambdaNP_grid)
            fNP_current=FNP(x,0d0,bT,hadron,lambdaNP)

            do j=-5,5
                if(fNP_grid(j)==0) then
                if(fNP_current(j)/=0 .and. ((j/=0).or.(.not.withGluon))) then
                if(outputLevel>0) &
                    call Warning_Raise('error in restoration: original value is zero. TMDPDF set to zero. b='&
                        //numToStr(bT)//' x='//numToStr(x)//' f='//numToStr(j),&
                        messageCounter,messageTrigger,moduleName)
                end if
                wgtTMDPDF_base5(j)=0!!!! this is case then 0/0
                else
                wgtTMDPDF_base5(j)=wgtTMDPDF_base5(j)*fNP_current(j)/fNP_grid(j)
                end if
                
                
            end do
            
            !!!! only (tw-2 part)*fNP is saved in grid. Therefore, we add tw3*fNP
            wgtTMDPDF_base5=wgtTMDPDF_base5+g1T_tw3NP(x,hadron,lambdaNP)*FNP(x,0d0,bT,hadron,lambdaNP)

        end if
    
    !!!! finally just calculation
    else  
        !!!! only (tw-2 part)*fNP is computed by convolution. Therefore, we add tw3*fNP
        !!!! the factor x is the part of definition of convolution for WW
        wgtTMDPDF_base5=x*Common_lowScale5(x,bT,hadron)+g1T_tw3NP(x,hadron,lambdaNP)*FNP(x,0d0,bT,hadron,lambdaNP)
    end if

end function wgtTMDPDF_base5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !---------------------------------------------------------------------
!- This is the TMD function evaluated for all quarks simultaniously (-5..5) at x,bT,mu the GLUON INCLUDED
!--    f =  -5, -4, -3,  -2,  -1,0,1,2,3, 4 ,5
!--      = bbar ,cbar sbar,ubar,dbar,g,d,u,s, c ,b
!!---Base version: hadron=number of PDF
!---------------------------------------------------------------------
!---------------------------------------------------------------------
function wgtTMDPDF_base50(x,bT,hadron)
    real(dp),dimension(-5:5)::wgtTMDPDF_base50
    real(dp) :: x, bT
    integer::hadron

    !!! variables for restoration only
    real(dp),dimension(-5:5) :: fNP_grid,fNP_current
    integer::j

    !!!! if the coefficient function is not used at all. Just return fNP
    if(order_global<0) then

        wgtTMDPDF_base50=fNP(x,1d0,bT,hadron,lambdaNP)

    !!! in the case the greed has been calculated FOR THIS HADRON
    else if(gridReady .and. ANY(hadronsInGRID.eq.hadron)) then 
        wgtTMDPDF_base50=x*ExtractFromGrid(x,bT,hadron)

        !!!!!!!!!!This is procedure of restoration of function from the initial grid
        !!! if fNP is x-independent then the value can be obtained by TMDPDF(initial) fNP(current)/fNP(initial)
        if(.not.IsFnpZdependent) then
            fNP_grid=FNP(x,0d0,bT,hadron,lambdaNP_grid)
            fNP_current=FNP(x,0d0,bT,hadron,lambdaNP)

            do j=-5,5
                if(fNP_grid(j)==0) then
                if(wgtTMDPDF_base50(j)/=0.and.j/=0) then
                if(outputLevel>0) &
                    call Warning_Raise('error in restoration: original value is zero. TMDPDF set to zero. b='//numToStr(bT),&
                        messageCounter,messageTrigger,moduleName)
                end if
                wgtTMDPDF_base50(j)=0
                else
                wgtTMDPDF_base50(j)=wgtTMDPDF_base50(j)*fNP_current(j)/fNP_grid(j)
                end if
            end do
            
            !!!! only (tw-2 part)*fNP is saved in grid. Therefore, we add tw3*fNP
            !!!! the factor x is the part of definition of convolution for WW
            wgtTMDPDF_base50=wgtTMDPDF_base50+g1T_tw3NP(x,hadron,lambdaNP)*FNP(x,0d0,bT,hadron,lambdaNP)
        end if
    
    !!!! finally just calculation
    else    
        !!!! only (tw-2 part)*fNP is computed by convolution. Therefore, we add tw3*fNP
        wgtTMDPDF_base50=x*Common_lowScale50(x,bT,hadron)+g1T_tw3NP(x,hadron,lambdaNP)*FNP(x,0d0,bT,hadron,lambdaNP)
    end if
end function wgtTMDPDF_base50

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
function wgtTMDPDF_lowScale5(x,bT,hadron)
  real(dp),dimension(-5:5)::wgtTMDPDF_lowScale5
  real(dp) :: x, bT
  integer::hadron
  
  logical,allocatable::includeInComposition(:)
  real(dp),allocatable::compositionCoefficients(:)
  integer::j,jN
  
  if(x>1d0) then
    wgtTMDPDF_lowScale5=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
    return
  end if
  
  if(IsComposite) then
    
    wgtTMDPDF_lowScale5=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
    
    call GetCompositionArray(hadron,lambdaNP,includeInComposition,compositionCoefficients)
    jN=size(includeInComposition)
    
    do j=1,jN
      if(includeInComposition(j)) then
	wgtTMDPDF_lowScale5=wgtTMDPDF_lowScale5+compositionCoefficients(j)*wgtTMDPDF_base5(x,bT,j)
      end if
    end do
    
    deallocate(includeInComposition,compositionCoefficients)
    
  else
  
   wgtTMDPDF_lowScale5=wgtTMDPDF_base5(x,bT,hadron)  
  end if
  
end function wgtTMDPDF_lowScale5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !---------------------------------------------------------------------
!- This is the TMD function evaluated for all quarks simultaniously (-5..5) at x,bT,mu the GLUON INCLUDED
!--    f =  -5, -4, -3,  -2,  -1,0,1,2,3, 4 ,5
!--      = bbar ,cbar sbar,ubar,dbar,g,d,u,s, c ,b
!!---Full version: hadron=number of PDF (if compositeness OFF)
!!---		   hadron=sum components (if compositeness ON)
!---------------------------------------------------------------------
!---------------------------------------------------------------------
function wgtTMDPDF_lowScale50(x,bT,hadron)
  real(dp),dimension(-5:5)::wgtTMDPDF_lowScale50
  real(dp) :: x, bT
  integer::hadron
  
  logical,allocatable::includeInComposition(:)
  real(dp),allocatable::compositionCoefficients(:)
  integer::j,jN
  
  if(x>1d0) then
    wgtTMDPDF_lowScale50=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
    return
  end if
  
  if(IsComposite) then
    
    wgtTMDPDF_lowScale50=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
    
    call GetCompositionArray(hadron,lambdaNP,includeInComposition,compositionCoefficients)
    jN=size(includeInComposition)
    
    do j=1,jN
      if(includeInComposition(j)) then
	wgtTMDPDF_lowScale50=wgtTMDPDF_lowScale50+compositionCoefficients(j)*wgtTMDPDF_base50(x,bT,j)
      end if
    end do
    
    deallocate(includeInComposition,compositionCoefficients)
  
  else
  
   wgtTMDPDF_lowScale50=wgtTMDPDF_base50(x,bT,hadron)  
  end if
end function wgtTMDPDF_lowScale50
