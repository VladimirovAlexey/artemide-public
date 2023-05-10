!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!
!!!!!!!!!       part of uTMDPDF module for artemide
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
      
      xf=xPDF(x,Q,hadron)
      
  end function xf


!---------------------------------------------------------------------
!- This is the TMD function evaluated for all quarks simultaniously (-5..5) at x,bT,mu
!--    f =  -5, -4, -3,  -2,  -1,0,1,2,3, 4 ,5
!--      = bbar ,cbar sbar,ubar,dbar,??,d,u,s, c ,b
!!---Gluon contribution is undefined
!!---Base version: hadron=number of PDF
!---------------------------------------------------------------------
!---------------------------------------------------------------------
function uTMDPDF_base5(x,bT,hadron)
    real(dp),dimension(-5:5)::uTMDPDF_base5
    real(dp) :: x, bT
    integer::hadron

    !!! variables for restoration only
    real(dp),dimension(-5:5) :: fNP_grid,fNP_current
    integer::j

    !!!! if the coefficient function is not used at all. Just return fNP
    if(order_global==-50) then

        uTMDPDF_base5=fNP(x,1d0,bT,hadron,lambdaNP)
        
    !!! in the case the greed has been calculated AND the hadron is in the grid
    else if(gridReady .and. ANY(hadronsInGRID.eq.hadron)) then 

        uTMDPDF_base5=ExtractFromGrid(x,bT,hadron)

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
                uTMDPDF_base5(j)=0!!!! this is case then 0/0
                else
                uTMDPDF_base5(j)=uTMDPDF_base5(j)*fNP_current(j)/fNP_grid(j)
                end if
            end do

        end if
    
    !!!! finally just calculation
    else
        uTMDPDF_base5=Common_lowScale5(x,bT,hadron)
    end if

end function uTMDPDF_base5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !---------------------------------------------------------------------
!- This is the TMD function evaluated for all quarks simultaniously (-5..5) at x,bT,mu the GLUON INCLUDED
!--    f =  -5, -4, -3,  -2,  -1,0,1,2,3, 4 ,5
!--      = bbar ,cbar sbar,ubar,dbar,g,d,u,s, c ,b
!!---Base version: hadron=number of PDF
!---------------------------------------------------------------------
!---------------------------------------------------------------------
function uTMDPDF_base50(x,bT,hadron)
    real(dp),dimension(-5:5)::uTMDPDF_base50
    real(dp) :: x, bT
    integer::hadron

    !!! variables for restoration only
    real(dp),dimension(-5:5) :: fNP_grid,fNP_current
    integer::j

    !!!! if the coefficient function is not used at all. Just return fNP
    if(order_global==-50) then

        uTMDPDF_base50=fNP(x,1d0,bT,hadron,lambdaNP)

    !!! in the case the greed has been calculated FOR THIS HADRON
    else if(gridReady .and. ANY(hadronsInGRID.eq.hadron)) then 
        uTMDPDF_base50=ExtractFromGrid(x,bT,hadron)

        !!!!!!!!!!This is procedure of restoration of function from the initial grid
        !!! if fNP is x-independent then the value can be obtained by TMDPDF(initial) fNP(current)/fNP(initial)
        if(.not.IsFnpZdependent) then
            fNP_grid=FNP(x,0d0,bT,hadron,lambdaNP_grid)
            fNP_current=FNP(x,0d0,bT,hadron,lambdaNP)

            do j=-5,5
                if(fNP_grid(j)==0) then
                if(uTMDPDF_base50(j)/=0.and.j/=0) then
                if(outputLevel>0) &
                    call Warning_Raise('error in restoration: original value is zero. TMDPDF set to zero. b='//numToStr(bT),&
                        messageCounter,messageTrigger,moduleName)
                end if
                uTMDPDF_base50(j)=0
                else
                uTMDPDF_base50(j)=uTMDPDF_base50(j)*fNP_current(j)/fNP_grid(j)
                end if
            end do
        end if
    
    !!!! finally just calculation
    else
        uTMDPDF_base50=Common_lowScale50(x,bT,hadron)
    end if
end function uTMDPDF_base50

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
function uTMDPDF_lowScale5(x,bT,hadron)
  real(dp),dimension(-5:5)::uTMDPDF_lowScale5
  real(dp) :: x, bT
  integer::hadron
  
  logical,allocatable::includeInComposition(:)
  real(dp),allocatable::compositionCoefficients(:)
  integer::j,jN
  
  if(x>1d0) then
    uTMDPDF_lowScale5=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
    return
  end if
  
  if(IsComposite) then
    
    uTMDPDF_lowScale5=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
    
    call GetCompositionArray(hadron,lambdaNP,includeInComposition,compositionCoefficients)
    jN=size(includeInComposition)
    
    do j=1,jN
      if(includeInComposition(j)) then
	uTMDPDF_lowScale5=uTMDPDF_lowScale5+compositionCoefficients(j)*uTMDPDF_base5(x,bT,j)
      end if
    end do
    
    deallocate(includeInComposition,compositionCoefficients)
    
  else
  
   uTMDPDF_lowScale5=uTMDPDF_base5(x,bT,hadron)  
  end if
  
end function uTMDPDF_lowScale5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !---------------------------------------------------------------------
!- This is the TMD function evaluated for all quarks simultaniously (-5..5) at x,bT,mu the GLUON INCLUDED
!--    f =  -5, -4, -3,  -2,  -1,0,1,2,3, 4 ,5
!--      = bbar ,cbar sbar,ubar,dbar,g,d,u,s, c ,b
!!---Full version: hadron=number of PDF (if compositeness OFF)
!!---		   hadron=sum components (if compositeness ON)
!---------------------------------------------------------------------
!---------------------------------------------------------------------
function uTMDPDF_lowScale50(x,bT,hadron)
  real(dp),dimension(-5:5)::uTMDPDF_lowScale50
  real(dp) :: x, bT
  integer::hadron
  
  logical,allocatable::includeInComposition(:)
  real(dp),allocatable::compositionCoefficients(:)
  integer::j,jN
  
  if(x>1d0) then
    uTMDPDF_lowScale50=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
    return
  end if
  
  if(IsComposite) then
    
    uTMDPDF_lowScale50=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
    
    call GetCompositionArray(hadron,lambdaNP,includeInComposition,compositionCoefficients)
    jN=size(includeInComposition)
    
    do j=1,jN
      if(includeInComposition(j)) then
	uTMDPDF_lowScale50=uTMDPDF_lowScale50+compositionCoefficients(j)*uTMDPDF_base50(x,bT,j)
      end if
    end do
    
    deallocate(includeInComposition,compositionCoefficients)
  
  else
  
   uTMDPDF_lowScale50=uTMDPDF_base50(x,bT,hadron)  
  end if
end function uTMDPDF_lowScale50
