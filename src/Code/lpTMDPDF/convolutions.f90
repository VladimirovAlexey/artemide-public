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
      
      xf=x_lp_PDF(x,Q,hadron)
      
  end function xf


!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Convolutions!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!---------------------NO SUCH FUNCTION IN HERE
! function lpTMDPDF_base5(x,bT,hadron)
!   real(dp),dimension(-5:5)::lpTMDPDF_base5
!   real(dp) :: x, bT
!   integer::hadron
!   
!   lpTMDPDF_base5=0d0
!  
!  end function lpTMDPDF_base5

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !---------------------------------------------------------------------
!- This is the TMD function evaluated for all quarks simultaniously (-5..5) at x,bT,mu the GLUON INCLUDED
!--    f =  -5, -4, -3,  -2,  -1,0,1,2,3, 4 ,5
!--      = bbar ,cbar sbar,ubar,dbar,g,d,u,s, c ,b
!- The order is accumulative pertrubative order of coefficient =0,1,2 (LO,NLO,NNLO)
!!---Base version: hadron=number of PDF
!---------------------------------------------------------------------
!---------------------------------------------------------------------
function lpTMDPDF_base50(x,bT,hadron)
    real(dp),dimension(-5:5)::lpTMDPDF_base50
    real(dp) :: x, bT
    integer::hadron

    !!! variables for restoration only
    real(dp),dimension(-5:5) :: fNP_grid,fNP_current
    integer::j

    !!!! if the coefficient function is not used at all. Just return fNP
    if(order_global==-50) then

        lpTMDPDF_base50=fNP(x,1d0,bT,hadron,lambdaNP)

    !!! in the case the greed has been calculated
    else if(gridReady .and. ANY(hadronsInGRID.eq.hadron)) then 
        lpTMDPDF_base50=ExtractFromGrid(x,bT,hadron)

        !!!!!!!!!!This is procedure of restoration of function from the initial grid
        !!! if fNP is x-independent then the value can be obtained by TMDPDF(initial) fNP(current)/fNP(initial)
        if(.not.IsFnpZdependent) then
            fNP_grid=FNP(x,0d0,bT,hadron,lambdaNP_grid)
            fNP_current=FNP(x,0d0,bT,hadron,lambdaNP)

            do j=-5,5
                if(fNP_grid(j)==0) then
                if(lpTMDPDF_base50(j)/=0.and.j/=0) then
                if(outputlevel>0) &
                    call Warning_Raise('error in restoration: original value is zero. TMDPDF set to zero. b='//numToStr(bT),&
                        messageCounter,messageTrigger,moduleName)
                end if
                lpTMDPDF_base50(j)=0
                else
                lpTMDPDF_base50(j)=lpTMDPDF_base50(j)*fNP_current(j)/fNP_grid(j)
                end if
            end do
        end if

    !!!! Finally, Just  calculation
    else
        lpTMDPDF_base50=Common_lowScale50(x,bT,hadron)
    end if
end function lpTMDPDF_base50

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!OUTPUT INTERFACE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!    NO SUCH FUNCTION HERE
! function lpTMDPDF_lowScale5(x,bT,hadron)
!   real(dp),dimension(-5:5)::lpTMDPDF_lowScale5
!   real(dp) :: x, bT
!   integer::hadron
! end function lpTMDPDF_lowScale5
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
 !---------------------------------------------------------------------
!- This is the TMD function evaluated for all quarks simultaniously (-5..5) at x,bT,mu the GLUON INCLUDED
!--    f =  -5, -4, -3,  -2,  -1,0,1,2,3, 4 ,5
!--      = bbar ,cbar sbar,ubar,dbar,g,d,u,s, c ,b
!!---Full version: hadron=number of PDF (if compositeness OFF)
!!---		   hadron=sum components (if compositeness ON)
!---------------------------------------------------------------------
!---------------------------------------------------------------------
function lpTMDPDF_lowScale50(x,bT,hadron)
  real(dp),dimension(-5:5)::lpTMDPDF_lowScale50
  real(dp) :: x, bT
  integer::hadron
  
  logical,allocatable::includeInComposition(:)
  real(dp),allocatable::compositionCoefficients(:)
  integer::j,jN
  
  if(x>1d0) then
    lpTMDPDF_lowScale50=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
    return
  end if
  
  if(IsComposite) then
    
    lpTMDPDF_lowScale50=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)
    
    call GetCompositionArray(hadron,lambdaNP,includeInComposition,compositionCoefficients)
    jN=size(includeInComposition)
    
    do j=1,jN
      if(includeInComposition(j)) then
	lpTMDPDF_lowScale50=lpTMDPDF_lowScale50+compositionCoefficients(j)*lpTMDPDF_base50(x,bT,j)
      end if
    end do
    
    deallocate(includeInComposition,compositionCoefficients)
  
  else
  
   lpTMDPDF_lowScale50=lpTMDPDF_base50(x,bT,hadron)  
  end if
end function lpTMDPDF_lowScale50
