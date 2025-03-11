!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 3.0
!
!	This file contains the part of the code, which is common for all TMD-OPE modules
!       that operates at twist-2. It is inclucded (as a text).
!	Such idiotic structure is needed since, FORTRAN does not allow inheritance.
!
!	Be AWARE of possible clash of variable names.
!
!	This part is devoted to the calculation of Mellin convolution between
!       PDF and the coefficient function.
!       CxF_compute(x,b,hadron) = real(-5:5)
!
!	v.3.00 Created (AV, 24.07.2023)
!
!				A.Vladimirov (24.07.2023)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!Evaluate Log[muOPE^2 bT^2/4/Exp(-2Gamma_E)]
!! the b here is b*,
!! Functions muOPE and bSTAR are defined in _OPE_model
!! x is global x, y is the convolution variable
!! this funciton is used only in Coefficeint functions
pure function LogMuB(bT,x,y)
    real(dp)::LogMuB
    real(dp),intent(in)::bT,x,y
    LogMuB=2._dp*Log(bSTAR(bT,x,y)*muOPE(bt,x,y,c4_global)*C0_inv_const)
end function LogMuB

!!! test MU for y-dependance
!!! the test consists in the evaluation of FNP at several random sets and NParray
!!! and comparison of the values.
!!! testMU=true muOPE is dependent on y
function TestMU()
    logical::TestMU
    real(dp)::xR,yR,bR
    real(dp)::test1,test2
    integer::i
    TestMU=.false.
    do i=1,10
        call RANDOM_NUMBER(bR)
        bR=5d0*bR
        if(xR>0.99d0) xR=xR/2d0
        if(xR<0.00001d0) xR=0.0001d0+xR
            !!! generate some random input
        call RANDOM_NUMBER(yR)
        if(yR>0.99d0) yR=yR/2d0
        if(yR<0.00001d0) yR=0.0001d0+yR
        test1=muOPE(bR,xR,yR,1._dp)

        !!! generate some random input
        call RANDOM_NUMBER(yR)
        if(yR>0.99d0) yR=yR/2d0
        if(yR<0.00001d0) yR=0.0001d0+yR
        test2=muOPE(bR,xR,yR,1._dp)

        if(ABS(test1-test2)>1d-8) then
            TestMU=.true.
            exit
        end if
    end do
end function TestMU

!!!! The CxF is defined as Mellin convolution of [C xF]
!!!! The integral is \int_x^1 dy C(y) F[x/y], where F[x]=x PDF(x)
!!!! The ordinary convolution \int_x^1 dy/y C(y) PDF(x/y) = CxF(x)/x
!!!! the parameters x,y for model are defined as in this formula
function CxF_compute(x,bT,hadron,includeGluon)
    real(dp),dimension(-5:5)::CxF_compute
    integer, intent(in)::hadron
    real(dp),intent(in)::x,bT
    logical,intent(in)::includeGluon
    real(dp):: bTcurrent,lx

    real(dp),dimension(-5:5)::deltaPart,PLUSremnant
    real(dp):: muAt1,asAt1,LogAt1,Cqq,Cgg,Csingqq,Csinggg
    integer::NfAt1
    real(dp),dimension(1:3)::CplusAt1_gg,CplusAt1_qq
    real(dp),dimension(-5:5)::PDFat1

    real(dp),dimension(1:parametrizationLength):: Bqq,Bqg,Bgq,Bgg,Bqqb,Bqqp
    integer::i

    real(dp)::yCUT
    real(dp),parameter::xCUT=0.99_dp

    !! for extrimely small-values of b we freeze it.
    if(bT>bMIN) then
        bTcurrent=bT
    else
        bTcurrent=bMIN
    end if

    !! drop the case of x>1
    if(x>1._dp-toleranceGEN) then
        CxF_compute=0._dp
        return
    end if

    !!! values of parameters at y=1
    !!! they are used also later
    muAt1=muOPE(bTcurrent,x,1._dp,c4_global)
    asAt1=As(muAt1)
    LogAt1=LogMuB(bTcurrent,x,1._dp)
    NfAt1=activeNf(muAt1)
    PDFat1=xf(x,muAt1,hadron)
    if(ISNAN(PDFat1(-5))) then
    write(*,*) "---------->",x,muAt1,hadron
    error stop
    end if

    !!!! delta-part
    !! C(y)~delta(1-y)
    Cqq=C_q_q_delta(asAt1,NfAt1,LogAt1)
    if(includeGluon) then
        Cgg=C_g_g_delta(asAt1,NfAt1,LogAt1)
    else
        Cgg=0._dp
    end if
    deltaPart=(/Cqq,Cqq,Cqq,Cqq,Cqq,Cgg,Cqq,Cqq,Cqq,Cqq,Cqq/)*PDFat1

    CxF_compute=deltaPart

    !!!! other parts contribute only if order >LO
    if(orderMain>0) then


    lx=Log(1._dp-x)
    !! this value is used in the integration over 1/(..)_+
    CplusAt1_qq=Coeff_q_q_plus(asAt1,NfAt1,LogAt1)

    Csingqq=sum(CplusAt1_qq*(/lx,lx**2/2._dp,lx**3/3._dp/))
    if(includeGluon) then
        CplusAt1_gg=Coeff_g_g_plus(asAt1,NfAt1,LogAt1)
        Csinggg=sum(CplusAt1_gg*(/lx,lx**2/2._dp,lx**3/3._dp/))
    else
        CplusAt1_gg=0._dp
        Csinggg=0._dp
    end if

    !!! account the remnant of the integration over 1/(..)_+
    !!! it is equal \int_0^x c(y) f(1)
    PLUSremnant=(/Csingqq,Csingqq,Csingqq,Csingqq,Csingqq,&
    Csinggg,Csingqq,Csingqq,Csingqq,Csingqq,Csingqq/)*PDFat1

    !!! if mu is y-independent then one can use the value of coeff at y=1
    !!! and do not update them for each iteration of the integral
    if(.not.IsMuYdependent) then
        Bqq=Coeff_q_q_reg(asAt1,NfAt1,LogAt1)
        Bqg=Coeff_q_g_reg(asAt1,NfAt1,LogAt1)
        if(includeGluon) then
        Bgq=Coeff_g_q_reg(asAt1,NfAt1,LogAt1)
        Bgg=Coeff_g_g_reg(asAt1,NfAt1,LogAt1)
        else
        Bgq=0._dp
        Bgg=0._dp
        end if
        Bqqb=Coeff_q_qb_reg(asAt1,NfAt1,LogAt1)
        Bqqp=Coeff_q_qp_reg(asAt1,NfAt1,LogAt1)
    end if

    !!! for smaller x, the part y~1 can be computed approximately (the xCUT is necesary since if x~y the error grows)
    !!! however if x is large ~1, the  should be closer to 1.
    yCUT=0.9999_dp
    if(x>xCUT) then
       yCUT=(100._dp+x)/101._dp
    end if

    CxF_compute=CxF_compute&
            +Integrate_GK_array5(FFplus,x,yCUT,toleranceINT)+PLUSremnant+Integrate_largey_PLUS()&
            +Integrate_GK_array5(FFreg,x,yCUT,toleranceINT)+Integrate_largey_Reg()

    end if

  do i=-5,5
   if(ISNAN(CxF_compute(i))) then

    write(*,*) ErrorString('convolution computed to NAN. CHECK INTEGRATION',moduleName)
    write(*,*) '----- information on last call -----'
    write(*,*) 'x=', x, 'bT=',bT,' i=',i, 'hadron=',hadron,' result=',CxF_compute(i)
    write(*,*) "--FF-->",PDFat1
    write(*,*) "delta->",(/Csingqq,Csingqq,Csingqq,Csingqq,Csingqq,&
    Csinggg,Csingqq,Csingqq,Csingqq,Csingqq,Csingqq/)
    error stop
   end if
  end do

!     !!!!! these are wrapper functions to pass ther integrand to GK routine
!     !!!!! FORTRAN gets only function of one variable
!     !!!!! using the trick with internal function one can by-pass this limitation
    contains

    !!! regular integrand
    function FFreg(y)
        real(dp),dimension(-5:5)::FFreg
        real(dp),intent(in)::y
        real(dp)::muCurrent,asCurrent,LogCurrent,Aqq,Aqg,Agq,Agg,Aqqp,Aqqb,PDFsum
        integer::NfCurrent
        real(dp),dimension(-5:5)::PDFs
        real(dp),dimension(1:parametrizationLength):: var

!         !!! very rare error, if y~1 (up to machine precision) freeze it!
!         if(y<0.999999999d0) then
            var=parametrizationString(y)
!         else
!             var=parametrizationString(0.999999999d0)
!         end if

        !!! if mu is y-dependent one needs to update the values of parameters for each y
        if(IsMuYdependent) then
            muCurrent=muOPE(bTcurrent,x,y,c4_global)
            asCurrent=As(muCurrent)
            LogCurrent=LogMuB(bTcurrent,x,y)
            NfCurrent=activeNf(muCurrent)

            Aqq=sum(var*Coeff_q_q_reg(asCurrent,NfCurrent,LogCurrent))
            Aqg=sum(var*Coeff_q_g_reg(asCurrent,NfCurrent,LogCurrent))
            if(includeGluon) then
            Agq=sum(var*Coeff_g_q_reg(asCurrent,NfCurrent,LogCurrent))
            Agg=sum(var*Coeff_g_g_reg(asCurrent,NfCurrent,LogCurrent))
            else
            Agq=0._dp
            Agg=0._dp
            end if
            Aqqb=sum(var*Coeff_q_qb_reg(asCurrent,NfCurrent,LogCurrent))
            Aqqp=sum(var*Coeff_q_qp_reg(asCurrent,NfCurrent,LogCurrent))

            PDFs=xf(x/y,muCurrent,hadron)
        else
            Aqq=sum(var*Bqq)
            Aqg=sum(var*Bqg)
            if(includeGluon) then
            Agq=sum(var*Bgq)
            Agg=sum(var*Bgg)
            else
            Agq=0._dp
            Agg=0._dp
            end if
            Aqqb=sum(var*Bqqb)
            Aqqp=sum(var*Bqqp)

            PDFs=xf(x/y,muAt1,hadron)
        end if


        PDFsum=PDFs(-5)+PDFs(-4)+PDFs(-3)+PDFs(-2)+PDFs(-1)+PDFs(1)+PDFs(2)+PDFs(3)+PDFs(4)+PDFs(5)

        FFreg=(/&
        Aqq*PDFs(-5)+Aqqb*PDFs(5)+Aqqp*PDFsum+Aqg*PDFs(0),&
        Aqq*PDFs(-4)+Aqqb*PDFs(4)+Aqqp*PDFsum+Aqg*PDFs(0),&
        Aqq*PDFs(-3)+Aqqb*PDFs(3)+Aqqp*PDFsum+Aqg*PDFs(0),&
        Aqq*PDFs(-2)+Aqqb*PDFs(2)+Aqqp*PDFsum+Aqg*PDFs(0),&
        Aqq*PDFs(-1)+Aqqb*PDFs(1)+Aqqp*PDFsum+Aqg*PDFs(0),&
        Agg*PDFs(0)+Agq*PDFsum,&
        Aqq*PDFs(1)+Aqqb*PDFs(-1)+Aqqp*PDFsum+Aqg*PDFs(0),&
        Aqq*PDFs(2)+Aqqb*PDFs(-2)+Aqqp*PDFsum+Aqg*PDFs(0),&
        Aqq*PDFs(3)+Aqqb*PDFs(-3)+Aqqp*PDFsum+Aqg*PDFs(0),&
        Aqq*PDFs(4)+Aqqb*PDFs(-4)+Aqqp*PDFsum+Aqg*PDFs(0),&
        Aqq*PDFs(5)+Aqqb*PDFs(-5)+Aqqp*PDFsum+Aqg*PDFs(0)/)
    end function FFreg

    !!! (..)_+ integrand
    function FFplus(y)
        real(dp),dimension(-5:5)::FFplus
        real(dp),intent(in)::y
        real(dp),dimension(1:3)::Cplus_qq,Cplus_gg,listLY
        real(dp),dimension(-5:5)::PDF,inter1,inter2
        real(dp)::muCurrent,asCurrent,LogCurrent,dummy1,dummy2,ly
        integer::NfCurrent

        !!! very rare error, if y~1 (up to machine precision) freeze it!
        if(y<0.999999999d0) then
            ly=log(1._dp-y)
        else
            ly=log(1._dp-0.999999999d0)
        end if
        listLY=(/1._dp,ly,ly**2/)

        !!!! if mu is y-dependent one needs to update the values each step
        if(IsMuYdependent) then
            muCurrent=muOPE(bTcurrent,x,y,c4_global)
            asCurrent=As(muCurrent)
            LogCurrent=LogMuB(bTcurrent,x,y)
            NfCurrent=activeNf(muCurrent)

            Cplus_qq=Coeff_q_q_plus(asCurrent,NfCurrent,LogCurrent)
            if(includeGluon) then
                Cplus_gg=Coeff_g_g_plus(asCurrent,NfCurrent,LogCurrent)
            else
                Cplus_gg=0._dp
            end if
            PDF=xf(x/y,muCurrent,hadron)


            dummy1=sum(listLY*Cplus_qq)
            dummy2=sum(listLY*Cplus_gg)
            inter1=(/dummy1,dummy1,dummy1,dummy1,dummy1,&
            dummy2,dummy1,dummy1,dummy1,dummy1,dummy1/)

            dummy1=sum(listLY*CplusAt1_qq)
            dummy2=sum(listLY*CplusAt1_gg)
            inter2=(/dummy1,dummy1,dummy1,dummy1,dummy1,&
            dummy2,dummy1,dummy1,dummy1,dummy1,dummy1/)

            FFplus=(inter1*PDF-inter2*PDFat1)/(1._dp-y)
        else

            PDF=xf(x/y,muAt1,hadron)

            dummy1=sum(listLY*CplusAt1_qq)
            dummy2=sum(listLY*CplusAt1_gg)
            inter2=(/dummy1,dummy1,dummy1,dummy1,dummy1,&
            dummy2,dummy1,dummy1,dummy1,dummy1,dummy1/)

            FFplus=inter2*(PDF-PDFat1)/(1._dp-y)
        end if

    end function FFplus


    !!!!! approximate integration of y~1. It is presice up to 5-6 digits for y~0.9999, and x<0.8
    !!!!! This is a part that result from the (..)_+ distribution
    !!!!! the approximation is like that
    !!!!! \int_y0^1 dy/y [f(x/y)-f(x)]/(1-y) \sim (1-y0) x f'(x) for y0->1
    !!!!! it is the same as [f(x/y0)-f(x)] for y0->1
    !!!!! and similar for ~log-terms
    function Integrate_largey_PLUS()
        real(dp),dimension(-5:5)::Integrate_largey_PLUS
        real(dp),dimension(1:3)::Cplus_qq,Cplus_gg
        real(dp)::muCurrent,asCurrent,LogCurrent
        integer::NfCurrent
        real(dp),dimension(-5:5)::PDFs
        real(dp)::lY,dummy1,dummy2

        lY=log(1._dp-yCUT)

        !!!!The integrals over (..)_+ are
        !!!! int_y^1 dy (f[x/y]-f[x])/(1-y) ~ (f[x/y]-f[x])
        !!!! int_y^1 dy (f[x/y]-f[x])*log[1-y]/(1-y) ~ (f[x/y]-f[x])(Log[1-y]-1)
        !!!! int_y^1 dy (f[x/y]-f[x])*log[1-y]**2/(1-y) ~ (f[x/y]-f[x])(2+Log[1-y]*(log[1-y]-2))
        if(IsMuYdependent) then
            !!! if mu is y-dependent one needs to update the values of parameters for each y
            muCurrent=muOPE(bTcurrent,x,yCUT,c4_global)
            asCurrent=As(muCurrent)
            LogCurrent=LogMuB(bTcurrent,x,yCUT)
            NfCurrent=activeNf(muCurrent)

            PDFs=xf(x/yCUT,muCurrent,hadron)

            Cplus_qq=Coeff_q_q_plus(asCurrent,NfCurrent,LogCurrent)
            if(includeGluon) then
                Cplus_gg=Coeff_g_g_plus(asCurrent,NfCurrent,LogCurrent)
            else
                Cplus_gg=0._dp
            end if

            dummy1=sum((/1._dp,lY-1._dp,2._dp+lY*(lY-2._dp)/)*Cplus_qq)
            dummy2=sum((/1._dp,lY-1._dp,2._dp+lY*(lY-2._dp)/)*Cplus_gg)
            Integrate_largey_PLUS=(/dummy1,dummy1,dummy1,dummy1,dummy1,&
            dummy2,dummy1,dummy1,dummy1,dummy1,dummy1/)*(PDFs-PDFat1)
        else
            PDFs=xf(x/yCUT,muAt1,hadron)

            dummy1=sum((/1._dp,lY-1._dp,2._dp+lY*(lY-2._dp)/)*CplusAt1_qq)
            dummy2=sum((/1._dp,lY-1._dp,2._dp+lY*(lY-2._dp)/)*CplusAt1_gg)
            Integrate_largey_PLUS=(/dummy1,dummy1,dummy1,dummy1,dummy1,&
            dummy2,dummy1,dummy1,dummy1,dummy1,dummy1/)*(PDFs-PDFat1)
        end if

    end function Integrate_largey_PLUS

    !!!!! approximate integration of y~1. It is presice up to 5-6 digits for y~0.9999, and x<0.8
    !!!!! the approximation is done as
    !!!!! \int_y0^1 f(x/y) g[y] \sim f[x]\int_y0^1 g[y] the integrals over g[y] are taken analytically for y0->1
    function Integrate_largey_Reg()
        real(dp),dimension(-5:5)::Integrate_largey_Reg
        real(dp)::Aqq,Aqg,Agq,Agg,Aqqp,Aqqb,PDFsum
        real(dp),dimension(-5:5)::PDFs
        real(dp),dimension(1:parametrizationLength):: var

        !!! the integral over regular part is simpler
        !!!! it is just the value f[x]int_x^1 ...
        var=parametrizationStringAt1(yCUT)
        Aqq=sum(var*Bqq)
        Aqg=sum(var*Bqg)
        if(includeGluon) then
            Agq=sum(var*Bgq)
            Agg=sum(var*Bgg)
        else
            Agq=0._dp
            Agg=0._dp
        end if
        Aqqb=sum(var*Bqqb)
        Aqqp=sum(var*Bqqp)

        PDFs=xf(x,muAt1,hadron)

        PDFsum=PDFs(-5)+PDFs(-4)+PDFs(-3)+PDFs(-2)+PDFs(-1)+PDFs(1)+PDFs(2)+PDFs(3)+PDFs(4)+PDFs(5)

        Integrate_largey_Reg=(/&
        Aqq*PDFs(-5)+Aqqb*PDFs(5)+Aqqp*PDFsum+Aqg*PDFs(0),&
        Aqq*PDFs(-4)+Aqqb*PDFs(4)+Aqqp*PDFsum+Aqg*PDFs(0),&
        Aqq*PDFs(-3)+Aqqb*PDFs(3)+Aqqp*PDFsum+Aqg*PDFs(0),&
        Aqq*PDFs(-2)+Aqqb*PDFs(2)+Aqqp*PDFsum+Aqg*PDFs(0),&
        Aqq*PDFs(-1)+Aqqb*PDFs(1)+Aqqp*PDFsum+Aqg*PDFs(0),&
        Agg*PDFs(0)+Agq*PDFsum,&
        Aqq*PDFs(1)+Aqqb*PDFs(-1)+Aqqp*PDFsum+Aqg*PDFs(0),&
        Aqq*PDFs(2)+Aqqb*PDFs(-2)+Aqqp*PDFsum+Aqg*PDFs(0),&
        Aqq*PDFs(3)+Aqqb*PDFs(-3)+Aqqp*PDFsum+Aqg*PDFs(0),&
        Aqq*PDFs(4)+Aqqb*PDFs(-4)+Aqqp*PDFsum+Aqg*PDFs(0),&
        Aqq*PDFs(5)+Aqqb*PDFs(-5)+Aqqp*PDFsum+Aqg*PDFs(0)/)

    end function Integrate_largey_Reg

end function CxF_compute
