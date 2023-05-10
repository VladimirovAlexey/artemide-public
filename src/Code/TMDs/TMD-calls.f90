!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!
!!!!!!!!!       part of TMDs module for artemide
!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!-------------------------------------------------------------------------------------------------------
!!!-----------------------------Routines for composition of TMD distributions-----------------------------
!!!-------------------------------------------------------------------------------------------------------

!!!!!!!!!!!!!!!!!!!uTMDPDF!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!! upolarized TMDPDF
! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
function uTMDPDF_5_Ev(x,bt,muf,zetaf,hadron)
    real(dp)::uTMDPDF_5_Ev(-5:5)
    real(dp),intent(in):: x,bt,muf,zetaf    
    integer,intent(in)::hadron
    real(dp):: mui,Rkernel
    SELECT CASE(EvolutionType)
    CASE(1)!!!! improved D
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,c1_global*mu0(bt),1)
    CASE(2)!!!! improved gamma
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
    CASE(3)!!!! fixed mu
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
    CASE(4)!!!! exact solution via zeta-line
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
    END SELECT
    uTMDPDF_5_Ev=Rkernel*uTMDPDF_lowScale5(x,bT,hadron)

    !!! forcefully set =0 below threshold
    if(muf<mBOTTOM) then
    uTMDPDF_5_Ev(5)=0_dp
    uTMDPDF_5_Ev(-5)=0_dp
    end if
    if(muf<mCHARM) then
    uTMDPDF_5_Ev(4)=0_dp
    uTMDPDF_5_Ev(-4)=0_dp
    end if
end function uTMDPDF_5_Ev


!!!!!!!! upolarized TMDPDF
! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
function uTMDPDF_50_Ev(x,bt,muf,zetaf,hadron)
    real(dp)::uTMDPDF_50_Ev(-5:5)
    real(dp),intent(in):: x,bt,muf,zetaf
    integer,intent(in)::hadron
    real(dp):: mui,Rkernel ,RkernelG       

    SELECT CASE(EvolutionType)
    CASE(1)!!!! improved D
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,c1_global*mu0(bt),1)
        RkernelG=TMDR_Rzeta(bt,muf,zetaf,mui,c1_global*mu0(bt),0)
    CASE(2)!!!! improved gamma
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
        RkernelG=TMDR_Rzeta(bt,muf,zetaf,mui,0)
    CASE(3)!!!! fixed mu
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
        RkernelG=TMDR_Rzeta(bt,muf,zetaf,0)
    CASE(4)!!!! exact solution via zeta-line
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bT,muf,zetaf,mui,1)
        RkernelG=TMDR_Rzeta(bt,muf,zetaf,mui,0)
    END SELECT
    !uTMDPDF_50_Ev=Rkernel*uTMDPDF_lowScale50(x,bT,hadron)
    !uTMDPDF_50_Ev(0)=uTMDPDF_50_Ev(0)*RkernelG/Rkernel
    uTMDPDF_50_Ev=uTMDPDF_lowScale50(x,bT,hadron)*&
        (/Rkernel,Rkernel,Rkernel,Rkernel,Rkernel,RkernelG,Rkernel,Rkernel,Rkernel,Rkernel,Rkernel/)

    !!! forcefully set =0 below threshold
    if(muf<mBOTTOM) then
        uTMDPDF_50_Ev(5)=0_dp
        uTMDPDF_50_Ev(-5)=0_dp
    end if
    if(muf<mCHARM) then
        uTMDPDF_50_Ev(4)=0_dp
        uTMDPDF_50_Ev(-4)=0_dp
    end if

end function uTMDPDF_50_Ev

!!!!!!!! upolarized TMDPDF OPTIMAL
! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
function uTMDPDF_5_optimal(x,bt,hadron)
    real(dp)::uTMDPDF_5_optimal(-5:5)
    real(dp),intent(in):: x,bt
    integer,intent(in)::hadron

    uTMDPDF_5_optimal=uTMDPDF_lowScale5(x,bT,hadron)
end function uTMDPDF_5_optimal

!!!!!!!! upolarized TMDPDF
! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
function uTMDPDF_50_optimal(x,bt,hadron)
    real(dp)::uTMDPDF_50_optimal(-5:5)
    real(dp),intent(in):: x,bt
    integer,intent(in)::hadron

    uTMDPDF_50_optimal=uTMDPDF_lowScale50(x,bT,hadron)
end function uTMDPDF_50_optimal

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!! PRODUCTS FOR DY!!!!!!!!!!!!!!!!!!!!!!
!!! Product of quark*antiquark uTMDPDFs. Slightly faster then just product
! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
function uPDF_uPDF(x1,x2,bt,muf,zetaf,hadron1,hadron2)
    real(dp),dimension(-5:5)::uPDF_uPDF
    real(dp),intent(in):: x1,x2,bt,muf,zetaf    
    integer,intent(in)::hadron1,hadron2
    real(dp):: mui,Rkernel
    real(dp),dimension(-5:5)::tmd1,tmd2
    SELECT CASE(EvolutionType)
    CASE(1)!!!! improved D
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,c1_global*mu0(bt),1)
    CASE(2)!!!! improved gamma
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
    CASE(3)!!!! fixed mu
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
    CASE(4)!!!! exact solution via zeta-line
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
    END SELECT
    tmd1=uTMDPDF_lowScale5(x1,bT,hadron1)
    tmd2=uTMDPDF_lowScale5(x2,bT,hadron2)
    uPDF_uPDF=(Rkernel**2)*tmd1*(tmd2(5:-5:-1))

    !!! forcefully set =0 below threshold
    if(muf<mBOTTOM) then
        uPDF_uPDF(5)=0_dp
        uPDF_uPDF(-5)=0_dp
    end if
    if(muf<mCHARM) then
        uPDF_uPDF(4)=0_dp
        uPDF_uPDF(-4)=0_dp
    end if
end function uPDF_uPDF

!!! Product of quark*quark uTMDPDFs. Slightly faster then just product
! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
function uPDF_anti_uPDF(x1,x2,bt,muf,zetaf,hadron1,hadron2)
    real(dp),dimension(-5:5)::uPDF_anti_uPDF
    real(dp),intent(in):: x1,x2,bt,muf,zetaf
    integer,intent(in)::hadron1,hadron2
    real(dp):: mui,Rkernel
    real(dp),dimension(-5:5)::tmd1,tmd2

    SELECT CASE(EvolutionType)
    CASE(1)!!!! improved D
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,c1_global*mu0(bt),1)
    CASE(2)!!!! improved gamma
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
    CASE(3)!!!! fixed mu
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
    CASE(4)!!!! exact solution via zeta-line
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
    END SELECT
    tmd1=uTMDPDF_lowScale5(x1,bT,hadron1)
    tmd2=uTMDPDF_lowScale5(x2,bT,hadron2)
    uPDF_anti_uPDF=(Rkernel**2)*tmd1*tmd2

    !!! forcefully set =0 below threshold
    if(muf<mBOTTOM) then
        uPDF_anti_uPDF(5)=0_dp
        uPDF_anti_uPDF(-5)=0_dp
    end if
    if(muf<mCHARM) then
        uPDF_anti_uPDF(4)=0_dp
        uPDF_anti_uPDF(-4)=0_dp
    end if
end function uPDF_anti_uPDF

!!!!!!!!!!!!!!!!!!!uTMDFF!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!! upolarized TMDFF
! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
function uTMDFF_5_Ev(x,bt,muf,zetaf,hadron)
    real(dp)::uTMDFF_5_Ev(-5:5)
    real(dp),intent(in):: x,bt,muf,zetaf
    integer,intent(in)::hadron
    real(dp):: mui,Rkernel

    SELECT CASE(EvolutionType)
    CASE(1)!!!! improved D
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,c1_global*mu0(bt),1)
    CASE(2)!!!! improved gamma
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
    CASE(3)!!!! fixed mu
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
    CASE(4)!!!! exact solution via zeta-line
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
    END SELECT
    uTMDFF_5_Ev=Rkernel*uTMDFF_lowScale5(x,bT,hadron)

    !!! forcefully set =0 below threshold
    if(muf<mBOTTOM) then
        uTMDFF_5_Ev(5)=0_dp
        uTMDFF_5_Ev(-5)=0_dp
    end if
    if(muf<mCHARM) then
        uTMDFF_5_Ev(4)=0_dp
        uTMDFF_5_Ev(-4)=0_dp
    end if

end function uTMDFF_5_Ev

!!!!!!!! upolarized TMDFF
! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
function uTMDFF_50_Ev(x,bt,muf,zetaf,hadron)
    real(dp)::uTMDFF_50_Ev(-5:5)
    real(dp),intent(in):: x,bt,muf,zetaf
    integer,intent(in)::hadron
    real(dp):: mui,Rkernel ,RkernelG   

    SELECT CASE(EvolutionType)
    CASE(1)!!!! improved D
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,c1_global*mu0(bt),1)
        RkernelG=TMDR_Rzeta(bt,muf,zetaf,mui,c1_global*mu0(bt),0)
    CASE(2)!!!! improved gamma
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
        RkernelG=TMDR_Rzeta(bt,muf,zetaf,mui,0)
    CASE(3)!!!! fixed mu
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
        RkernelG=TMDR_Rzeta(bt,muf,zetaf,0)
    CASE(4)!!!! exact solution via zeta-line
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
        RkernelG=TMDR_Rzeta(bt,muf,zetaf,mui,0)
    END SELECT
    !uTMDFF_50_Ev=Rkernel*uTMDFF_lowScale50(x,bT,hadron)
    !uTMDFF_50_Ev(0)=uTMDFF_50_Ev(0)*RkernelG/Rkernel
    uTMDFF_50_Ev=uTMDFF_lowScale50(x,bT,hadron)*&
        (/Rkernel,Rkernel,Rkernel,Rkernel,Rkernel,RkernelG,Rkernel,Rkernel,Rkernel,Rkernel,Rkernel/)

    !!! forcefully set =0 below threshold
    if(muf<mBOTTOM) then
        uTMDFF_50_Ev(5)=0_dp
        uTMDFF_50_Ev(-5)=0_dp
    end if
    if(muf<mCHARM) then
        uTMDFF_50_Ev(4)=0_dp
        uTMDFF_50_Ev(-4)=0_dp
    end if

end function uTMDFF_50_Ev

! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
function uTMDFF_5_optimal(x,bt,hadron)
    real(dp)::uTMDFF_5_optimal(-5:5)
    real(dp),intent(in):: x,bt
    integer,intent(in)::hadron

    uTMDFF_5_optimal=uTMDFF_lowScale5(x,bT,hadron)

end function uTMDFF_5_optimal

! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
function uTMDFF_50_optimal(x,bt,hadron)
    real(dp)::uTMDFF_50_optimal(-5:5)
    real(dp),intent(in):: x,bt
    integer,intent(in)::hadron

    uTMDFF_50_optimal=uTMDFF_lowScale50(x,bT,hadron)

end function uTMDFF_50_optimal


!!!!!!!!!!!!!!!!!!!lpTMDPDF!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!! linearly polarized gluon TMDPDF
! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
! all quark terms are zero!
function lpTMDPDF_50_Ev(x,bt,muf,zetaf,hadron)
    real(dp)::lpTMDPDF_50_Ev(-5:5)
    real(dp),intent(in):: x,bt,muf,zetaf
    integer,intent(in)::hadron
    real(dp):: mui,RkernelG   

    SELECT CASE(EvolutionType)
    CASE(1)!!!! improved D
        mui=c3_global*mu_LOW(bt)
        RkernelG=TMDR_Rzeta(bt,muf,zetaf,mui,c1_global*mu0(bt),0)
    CASE(2)!!!! improved gamma
        mui=c3_global*mu_LOW(bt)      
        RkernelG=TMDR_Rzeta(bt,muf,zetaf,mui,0)
    CASE(3)!!!! fixed mu
        RkernelG=TMDR_Rzeta(bt,muf,zetaf,0)
    CASE(4)!!!! exact solution via zeta-line
        mui=c3_global*mu_LOW(bt)
        RkernelG=TMDR_Rzeta(bt,muf,zetaf,mui,0)
    END SELECT
    lpTMDPDF_50_Ev=RkernelG*lpTMDPDF_lowScale50(x,bT,hadron)

    !!! forcefully set=0 all quarks
    !     glTMDPDF_50_Ev(5)=0_dp
    !     glTMDPDF_50_Ev(4)=0_dp
    !     glTMDPDF_50_Ev(3)=0_dp
    !     glTMDPDF_50_Ev(2)=0_dp
    !     glTMDPDF_50_Ev(1)=0_dp
    !     glTMDPDF_50_Ev(-1)=0_dp
    !     glTMDPDF_50_Ev(-2)=0_dp
    !     glTMDPDF_50_Ev(-3)=0_dp
    !     glTMDPDF_50_Ev(-4)=0_dp
    !     glTMDPDF_50_Ev(-5)=0_dp

end function lpTMDPDF_50_Ev

!!!!!!!! linearly polarized gluon TMDPDF OPTIMAL
! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
! all quark terms are zero!
function lpTMDPDF_50_optimal(x,bt,hadron)
    real(dp)::lpTMDPDF_50_optimal(-5:5)
    real(dp),intent(in):: x,bt
    integer,intent(in)::hadron

    lpTMDPDF_50_optimal=lpTMDPDF_lowScale50(x,bT,hadron)

end function lpTMDPDF_50_optimal

!!!!!!!!!!!!!!!!!!!SiversTMDPDF!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!! upolarized TMDFF
! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
function SiversTMDPDF_5_Ev(x,bt,muf,zetaf,hadron)
    real(dp)::SiversTMDPDF_5_Ev(-5:5)
    real(dp),intent(in):: x,bt,muf,zetaf
    integer,intent(in)::hadron
    real(dp):: mui,Rkernel

    SELECT CASE(EvolutionType)
    CASE(1)!!!! improved D
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,c1_global*mu0(bt),1)
    CASE(2)!!!! improved gamma
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
    CASE(3)!!!! fixed mu
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
    CASE(4)!!!! exact solution via zeta-line
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
    END SELECT
    SiversTMDPDF_5_Ev=Rkernel*SiversTMDPDF_lowScale5(x,bT,hadron)

    !!! forcefully set =0 below threshold
    if(muf<mBOTTOM) then
        SiversTMDPDF_5_Ev(5)=0_dp
        SiversTMDPDF_5_Ev(-5)=0_dp
    end if
    if(muf<mCHARM) then
        SiversTMDPDF_5_Ev(4)=0_dp
        SiversTMDPDF_5_Ev(-4)=0_dp
    end if

end function SiversTMDPDF_5_Ev

!!!!!!!! upolarized TMDFF
! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
function SiversTMDPDF_50_Ev(x,bt,muf,zetaf,hadron)
    real(dp)::SiversTMDPDF_50_Ev(-5:5)
    real(dp),intent(in):: x,bt,muf,zetaf
    integer,intent(in)::hadron
    real(dp):: mui,Rkernel ,RkernelG   

    SELECT CASE(EvolutionType)
    CASE(1)!!!! improved D
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,c1_global*mu0(bt),1)
        RkernelG=TMDR_Rzeta(bt,muf,zetaf,mui,c1_global*mu0(bt),0)
    CASE(2)!!!! improved gamma
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
        RkernelG=TMDR_Rzeta(bt,muf,zetaf,mui,0)
    CASE(3)!!!! fixed mu
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
        RkernelG=TMDR_Rzeta(bt,muf,zetaf,0)
    CASE(4)!!!! exact solution via zeta-line
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
        RkernelG=TMDR_Rzeta(bt,muf,zetaf,mui,0)
    END SELECT
    !SiversTMDPDF_50_Ev=Rkernel*SiversTMDPDF_lowScale50(x,bT,hadron)
    !SiversTMDPDF_50_Ev(0)=SiversTMDPDF_50_Ev(0)*RkernelG/Rkernel
    SiversTMDPDF_50_Ev=SiversTMDPDF_lowScale50(x,bT,hadron)*&
        (/Rkernel,Rkernel,Rkernel,Rkernel,Rkernel,RkernelG,Rkernel,Rkernel,Rkernel,Rkernel,Rkernel/)

    !!! forcefully set =0 below threshold
    if(muf<mBOTTOM) then
        SiversTMDPDF_50_Ev(5)=0_dp
        SiversTMDPDF_50_Ev(-5)=0_dp
    end if
    if(muf<mCHARM) then
        SiversTMDPDF_50_Ev(4)=0_dp
        SiversTMDPDF_50_Ev(-4)=0_dp
    end if

end function SiversTMDPDF_50_Ev

! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
function SiversTMDPDF_5_optimal(x,bt,hadron)
    real(dp)::SiversTMDPDF_5_optimal(-5:5)
    real(dp),intent(in):: x,bt
    integer,intent(in)::hadron

    SiversTMDPDF_5_optimal=SiversTMDPDF_lowScale5(x,bT,hadron)

end function SiversTMDPDF_5_optimal

! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
function SiversTMDPDF_50_optimal(x,bt,hadron)
    real(dp)::SiversTMDPDF_50_optimal(-5:5)
    real(dp),intent(in):: x,bt
    integer,intent(in)::hadron

    SiversTMDPDF_50_optimal=SiversTMDPDF_lowScale50(x,bT,hadron)

end function SiversTMDPDF_50_optimal


!!!!!!!!!!!!!!!!!!!wgtTMDPDF!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!! upolarized TMDFF
! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
function wgtTMDPDF_5_Ev(x,bt,muf,zetaf,hadron)
    real(dp)::wgtTMDPDF_5_Ev(-5:5)
    real(dp),intent(in):: x,bt,muf,zetaf
    integer,intent(in)::hadron
    real(dp):: mui,Rkernel

    SELECT CASE(EvolutionType)
    CASE(1)!!!! improved D
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,c1_global*mu0(bt),1)
    CASE(2)!!!! improved gamma
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
    CASE(3)!!!! fixed mu
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
    CASE(4)!!!! exact solution via zeta-line
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
    END SELECT
    wgtTMDPDF_5_Ev=Rkernel*wgtTMDPDF_lowScale5(x,bT,hadron)

    !!! forcefully set =0 below threshold
    if(muf<mBOTTOM) then
        wgtTMDPDF_5_Ev(5)=0_dp
        wgtTMDPDF_5_Ev(-5)=0_dp
    end if
    if(muf<mCHARM) then
        wgtTMDPDF_5_Ev(4)=0_dp
        wgtTMDPDF_5_Ev(-4)=0_dp
    end if

end function wgtTMDPDF_5_Ev

!!!!!!!! upolarized TMDFF
! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
function wgtTMDPDF_50_Ev(x,bt,muf,zetaf,hadron)
    real(dp)::wgtTMDPDF_50_Ev(-5:5)
    real(dp),intent(in):: x,bt,muf,zetaf
    integer,intent(in)::hadron
    real(dp):: mui,Rkernel ,RkernelG   

    SELECT CASE(EvolutionType)
    CASE(1)!!!! improved D
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,c1_global*mu0(bt),1)
        RkernelG=TMDR_Rzeta(bt,muf,zetaf,mui,c1_global*mu0(bt),0)
    CASE(2)!!!! improved gamma
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
        RkernelG=TMDR_Rzeta(bt,muf,zetaf,mui,0)
    CASE(3)!!!! fixed mu
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,1)
        RkernelG=TMDR_Rzeta(bt,muf,zetaf,0)
    CASE(4)!!!! exact solution via zeta-line
        mui=c3_global*mu_LOW(bt)
        Rkernel=TMDR_Rzeta(bt,muf,zetaf,mui,1)
        RkernelG=TMDR_Rzeta(bt,muf,zetaf,mui,0)
    END SELECT
    !wgtTMDPDF_50_Ev=Rkernel*wgtTMDPDF_lowScale50(x,bT,hadron)
    !wgtTMDPDF_50_Ev(0)=wgtTMDPDF_50_Ev(0)*RkernelG/Rkernel
    wgtTMDPDF_50_Ev=wgtTMDPDF_lowScale50(x,bT,hadron)*&
        (/Rkernel,Rkernel,Rkernel,Rkernel,Rkernel,RkernelG,Rkernel,Rkernel,Rkernel,Rkernel,Rkernel/)

    !!! forcefully set =0 below threshold
    if(muf<mBOTTOM) then
        wgtTMDPDF_50_Ev(5)=0_dp
        wgtTMDPDF_50_Ev(-5)=0_dp
    end if
    if(muf<mCHARM) then
        wgtTMDPDF_50_Ev(4)=0_dp
        wgtTMDPDF_50_Ev(-4)=0_dp
    end if

end function wgtTMDPDF_50_Ev

! vector (bbar,cbar,sbar,ubar,dbar,??,d,u,s,c,b)
function wgtTMDPDF_5_optimal(x,bt,hadron)
    real(dp)::wgtTMDPDF_5_optimal(-5:5)
    real(dp),intent(in):: x,bt
    integer,intent(in)::hadron

    wgtTMDPDF_5_optimal=wgtTMDPDF_lowScale5(x,bT,hadron)

end function wgtTMDPDF_5_optimal

! vector (bbar,cbar,sbar,ubar,dbar,g,d,u,s,c,b)
function wgtTMDPDF_50_optimal(x,bt,hadron)
    real(dp)::wgtTMDPDF_50_optimal(-5:5)
    real(dp),intent(in):: x,bt
    integer,intent(in)::hadron

    wgtTMDPDF_50_optimal=wgtTMDPDF_lowScale50(x,bT,hadron)

end function wgtTMDPDF_50_optimal
