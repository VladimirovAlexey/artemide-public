!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!
!!!!!!!!!       part of TMDR module for artemide
!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!-------------------------------------------------------------------------------------------------------
!!!-----------------------------Routines for the evolution of type 3--------------------------------------
!!!----------------------------------(zeta-prescription)--------------------------------------------------
!!!-------------------------------------------------------------------------------------------------------

!!! Evolution exponent in the improved gamma-picture to zeta-line (defined by zetaNP)
 function TMDR_Rzeta_type3(b,muf,zetaf,f)
  real(dp)::TMDR_Rzeta_type3,b,muf,zetaf,zetaP
  integer::f
  
  if(b<1d-6) b=1d-6
  
  zetaP=zetaNP(muf,b,f)
  
  TMDR_Rzeta_type3=EXP(-DNP(muf,b,f)*Log(zetaf/zetaP))
 
!   write(*,*) 'HERE'
!   write(*,*) DNP(muf,b,f),zetaf,zetaP
  
  if(TMDR_Rzeta_type3>1d6) then
    write(*,*) ErrorString('Evolution factor(type3) is TOO HUGE check the formula',moduleName)
    write(*,*) 'b=',b,'zetaf=',zetaf,'muf=',muf,'zetaP=',zetaP    
    write(*,*) 'DNP=',DNP(muf,b,f), 'log(zeta/zetamu)=',Log(zetaf/zetaP)
    write(*,*) 'NPparameters= (',NPparam,')'
    write(*,*) 'Evaluation continue with R=10^6'
    TMDR_Rzeta_type3=1d6
    stop
  end if

 end function TMDR_Rzeta_type3

function TMDR_Rzeta_harpy(b,muf,zetaf,f)
  real(dp)::TMDR_Rzeta_harpy,b,muf,zetaf
  integer::f

  TMDR_Rzeta_harpy=TMDR_Rzeta_type3(b,muf,zetaf,f)
end function TMDR_Rzeta_harpy
