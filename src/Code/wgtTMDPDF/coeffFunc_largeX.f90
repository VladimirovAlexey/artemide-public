!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!
!!!!!!!!!       part of wgtTMDPDF module for artemide
!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!-------------------------------------------------------------------------------------------------------
!!!-----------------------------Routines for the mathing coefficient--------------------------------------
!!!--------------The order is accumulative pertrubative order of coefficient =1,2,3 (LO,NLO,NNLO)---------
!!!-------------------------------------------------------------------------------------------------------

!!!! parametrizationString=(/1d0,z,Log(z),Log(1d0-z)/)

!!!!! Some WW-terms have common factor to be multiplied by
!!!!! for exampel wgt=1, wgl=-z
pure function commonFactor_largeX(z)
    real(dp)::commonFactor_largeX
    real(dp), intent(in)::z
    commonFactor_largeX=1.d0
end function commonFactor_largeX

!!!!!coefficient function q<-q regular-part
pure function Coeff_q_q_reg_largeX(alpha,Nf,Lmu)
    real(dp),dimension(1:parametrizationLength)::Coeff_q_q_reg_largeX
    real(dp), intent(in)::alpha,Lmu
    integer,intent(in)::Nf

    !! the Leading order is 1, it is WW-part of worm-gear function
    Coeff_q_q_reg_largeX=(/0d0,0d0,0d0,0d0/) !1
    if(orderMain>=1) then
              
        Coeff_q_q_reg_largeX=Coeff_q_q_reg_largeX+alpha*4d0/3d0*(/&
        2d0*Lmu-2d0, -2d0*Lmu+2d0, 2d0*Lmu-2d0, 0d0/) !
        
    !  write(*,*) 'regularPart=', regularPart/x
    end if
end function Coeff_q_q_reg_largeX

!!!!!coefficient function q<-g regular-part  
pure function Coeff_q_g_reg_largeX(alpha,Nf,Lmu)
    real(dp),dimension(1:parametrizationLength)::Coeff_q_g_reg_largeX
    real(dp), intent(in)::alpha,Lmu
    integer,intent(in)::Nf

    !! the Leading order is always zero, therefore calculation should be done only for order >=1
    Coeff_q_g_reg_largeX=(/0d0,0d0,0d0,0d0/)
    if(orderMain>=1) then
        Coeff_q_g_reg_largeX=Coeff_q_g_reg_largeX+alpha*(/-2d0*Lmu+1d0, 2d0*Lmu-1d0, -Lmu+0.5d0, 0d0/)
    end if
end function Coeff_q_g_reg_largeX

!!!!!coefficient function g<-q regular-part  
pure function Coeff_g_q_reg_largeX(alpha,Nf,Lmu)
    real(dp),dimension(1:parametrizationLength)::Coeff_g_q_reg_largeX
    real(dp), intent(in)::alpha,Lmu
    integer,intent(in)::Nf

    !! the Leading order is always zero, therefore calculation should be done only for order >=1
    Coeff_g_q_reg_largeX=(/0d0,0d0,0d0,0d0/)
    !   if(orderMain>=1) then
    !     Coeff_g_q=Coeff_g_q+alpha*(/0d0/)!
    !     
    !   end if
end function Coeff_g_q_reg_largeX

    !!!!!coefficient function g<-g regular-part  
function Coeff_g_g_reg_largeX(alpha,Nf,Lmu)
    real(dp),dimension(1:parametrizationLength)::Coeff_g_g_reg_largeX
    real(dp), intent(in)::alpha,Lmu
    integer,intent(in)::Nf

    !! the Leading order is always zero, therefore calculation should be done only for order >=1
    Coeff_g_g_reg_largeX=(/1d0,0d0,0d0,0d0/)

    if(outputLevel>0) write(*,*) &
            WarningString('gluon part of worm-gear T function is not known. Set alike quark.',moduleName)

end function Coeff_g_g_reg_largeX

!!!!!coefficient function q<-qb regular-part  
pure function Coeff_q_qb_reg_largeX(alpha,Nf,Lmu)
    real(dp),dimension(1:parametrizationLength)::Coeff_q_qb_reg_largeX
    real(dp), intent(in)::alpha,Lmu
    integer,intent(in)::Nf

    !! the Leading order is always zero, therefore calculation should be done only for order >=1
    Coeff_q_qb_reg_largeX=(/0d0,0d0,0d0,0d0/)!

end function Coeff_q_qb_reg_largeX

    !!!!!coefficient function q<-qp regular-part  
pure function Coeff_q_qp_reg_largeX(alpha,Nf,Lmu)
    real(dp),dimension(1:parametrizationLength)::Coeff_q_qp_reg_largeX
    real(dp), intent(in)::alpha,Lmu
    integer,intent(in)::Nf

    !! the Leading order is always zero, therefore calculation should be done only for order >=1
    Coeff_q_qp_reg_largeX=(/0d0,0d0,0d0,0d0/)
end function Coeff_q_qp_reg_largeX
