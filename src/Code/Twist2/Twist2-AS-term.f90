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
!       CxF_AS(x,mu,hadron) = real(-5:5)
!   In the second-transverse moment [i.e.  \int d^2 kT kT^2] (see [2402.01836] for details)
!
!	v.3.00 Created (AV, 16.11.2023)
!
!				A.Vladimirov (16.11.2023)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!! The Asumptotic term is equal to the small-b matching with
!!!! the coefficient function -[c1+2c2+6 c3]=C_AS
!!!! where CC[L]=[c0+L c1+L^2 c2+L^3 c3] is the original N3LO coefficeint function
!!!! C_AS can be computed simply as
!!!! C_AS=3/4*CC[0]-2/3*CC[1]-CC[4]/12
!!!! If we compute with muOPE the coefficients contain logarithms


!!!! The following code is canibalized from the Twist2Convolution.f90 and simplified
!!!! it computes [C_AS x f](x,mu)/x
!!!! mu is the scale of tranformation
!!!! mu0 is the scale of OPE

!!!! The CxF is defined as Mellin convolution of [C xF]
!!!! The integral is \int_x^1 dy C(y) F[x/y], where F[x]=x PDF(x)
!!!! The ordinary convolution \int_x^1 dy/y C(y) PDF(x/y) = CxF(x)/x
!!!! the parameters x,y for model are defined as in this formula
function CxF_AS(x,mu,mu0,hadron,includeGluon)
    real(dp),dimension(-5:5)::CxF_AS
    integer, intent(in)::hadron
    real(dp),intent(in)::x,mu,mu0
    logical,intent(in)::includeGluon
    real(dp):: lx

    real(dp),dimension(-5:5)::deltaPart,PLUSremnant
    real(dp):: asAt1,Cqq,Cgg,Csingqq,Csinggg
    integer::NfAt1
    real(dp),dimension(1:3)::CplusAt1_gg,CplusAt1_qq
    real(dp),dimension(-5:5)::PDFat1

    real(dp),dimension(1:parametrizationLength):: Bqq,Bqg,Bgq,Bgg,Bqqb,Bqqp
    integer::i
    real(dp)::a0,a1,a2,a4,l0 !!! coefficients of C(L=0,1,2,4)

    real(dp)::yCUT
    real(dp),parameter::xCUT=0.99_dp

    !! for extremely small-values of mu return ERROR
    if(mu<0.8d0) then
        write(*,*) ErrorString("ERROR in AS-moment computation. mu<0.8",moduleName)
        stop
    end if
    if(mu0<0.8d0) then
        write(*,*) ErrorString("ERROR in AS-moment computation. muOPE<0.8",moduleName)
        stop
    end if

    l0=2._dp*Log(mu/mu0)
    a0=3._dp/4+l0+3._dp*l0**2/8
    a1=-2._dp/3-2*l0-l0**2
    a2=l0+3._dp*l0**2/4
    a4=-1._dp/12-l0**2/8

    !! drop the case of x>1
    if(x>1._dp) then
        CxF_AS=0._dp
        return
    end if

    !!! values of parameters at y=1
    !!! they are used also later
    asAt1=As(mu0)
    NfAt1=activeNf(mu0)
    PDFat1=xf(x,mu0,hadron)

    !!!! delta-part
    !! C(y)~delta(1-y)
    Cqq=a0*C_q_q_delta(asAt1,NfAt1,0._dp)+a1*C_q_q_delta(asAt1,NfAt1,1._dp)&
            +a2*C_q_q_delta(asAt1,NfAt1,2._dp)+a4*C_q_q_delta(asAt1,NfAt1,4._dp)
    if(includeGluon) then
        Cgg=a0*C_g_g_delta(asAt1,NfAt1,0._dp)+a1*C_g_g_delta(asAt1,NfAt1,1._dp)&
            +a2*C_g_g_delta(asAt1,NfAt1,2._dp)+a4*C_g_g_delta(asAt1,NfAt1,4._dp)
    else
        Cgg=0._dp
    end if
    deltaPart=(/Cqq,Cqq,Cqq,Cqq,Cqq,Cgg,Cqq,Cqq,Cqq,Cqq,Cqq/)*PDFat1

    CxF_AS=deltaPart

    !!!! other parts contribute only if order >LO
    if(orderMain>0) then


    lx=Log(1._dp-x)
    !! this value is used in the integration over 1/(..)_+
    CplusAt1_qq=a0*Coeff_q_q_plus(asAt1,NfAt1,0._dp)+a1*Coeff_q_q_plus(asAt1,NfAt1,1._dp)&
                +a2*Coeff_q_q_plus(asAt1,NfAt1,2._dp)+a4*Coeff_q_q_plus(asAt1,NfAt1,4._dp)
    Csingqq=sum(CplusAt1_qq*(/lx,lx**2/2._dp,lx**3/3._dp/))

    if(includeGluon) then
        CplusAt1_gg=&
            a0*Coeff_g_g_plus(asAt1,NfAt1,0._dp)+a1*Coeff_g_g_plus(asAt1,NfAt1,1._dp)&
                +a2*Coeff_g_g_plus(asAt1,NfAt1,2._dp)+a4*Coeff_g_g_plus(asAt1,NfAt1,4._dp)
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
    !!! IT IS Y-INDEPENT
    Bqq=a0*Coeff_q_q_reg(asAt1,NfAt1,0._dp)+a1*Coeff_q_q_reg(asAt1,NfAt1,1._dp)&
            +a2*Coeff_q_q_reg(asAt1,NfAt1,2._dp)+a4*Coeff_q_q_reg(asAt1,NfAt1,4._dp)
    Bqg=a0*Coeff_q_g_reg(asAt1,NfAt1,0._dp)+a1*Coeff_q_g_reg(asAt1,NfAt1,1._dp)&
            +a2*Coeff_q_g_reg(asAt1,NfAt1,2._dp)+a4*Coeff_q_g_reg(asAt1,NfAt1,4._dp)
    if(includeGluon) then
    Bgq=a0*Coeff_g_q_reg(asAt1,NfAt1,0._dp)+a1*Coeff_g_q_reg(asAt1,NfAt1,1._dp)&
            +a2*Coeff_g_q_reg(asAt1,NfAt1,2._dp)+a4*Coeff_g_q_reg(asAt1,NfAt1,4._dp)
    Bgg=a0*Coeff_g_g_reg(asAt1,NfAt1,0._dp)+a1*Coeff_g_g_reg(asAt1,NfAt1,1._dp)&
            +a2*Coeff_g_g_reg(asAt1,NfAt1,2._dp)+a4*Coeff_g_g_reg(asAt1,NfAt1,4._dp)
    else
    Bgq=0._dp
    Bgg=0._dp
    end if
    Bqqb=a0*Coeff_q_qb_reg(asAt1,NfAt1,0._dp)+a1*Coeff_q_qb_reg(asAt1,NfAt1,1._dp)&
            +a2*Coeff_q_qb_reg(asAt1,NfAt1,2._dp)+a4*Coeff_q_qb_reg(asAt1,NfAt1,4._dp)
    Bqqp=a0*Coeff_q_qp_reg(asAt1,NfAt1,0._dp)+a1*Coeff_q_qp_reg(asAt1,NfAt1,1._dp)&
            +a2*Coeff_q_qp_reg(asAt1,NfAt1,2._dp)+a4*Coeff_q_qp_reg(asAt1,NfAt1,4._dp)

    !!! for smaller x, the part y~1 can be computed approximately (the xCUT is necesary since if x~y the error grows)
    !!! however if x is large ~1, the  should be closer to 1.
    yCUT=0.9999_dp
    if(x>xCUT) then
       yCUT=(100._dp+x)/101._dp
    end if

    CxF_AS=CxF_AS+PLUSremnant&
            +Integrate_GK_array5(FFreg,x,yCUT,toleranceINT)&
            +Integrate_GK_array5(FFplus,x,yCUT,toleranceINT)&
            +Integrate_largey()

    end if

  do i=-5,5
   if(ISNAN(CxF_AS(i))) then

    write(*,*) ErrorString('AS convolution computed to NAN. CHECK INTEGRATION',moduleName)
    write(*,*) '----- information on last call -----'
    write(*,*) 'x=', x,' i=',i, 'hadron=',hadron,' result=',CxF_AS(i)

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
        real(dp)::Aqq,Aqg,Agq,Agg,Aqqp,Aqqb,PDFsum
        real(dp),dimension(-5:5)::PDFs
        real(dp),dimension(1:parametrizationLength):: var

        var=parametrizationString(y)

        Aqq=sum(var*Bqq)
        Aqg=sum(var*Bqg)
        Agq=sum(var*Bgq)
        Agg=sum(var*Bgg)
        Aqqb=sum(var*Bqqb)
        Aqqp=sum(var*Bqqp)

        PDFs=xf(x/y,mu0,hadron)


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
        real(dp),dimension(1:3)::listLY
        real(dp),dimension(-5:5)::PDF,inter2
        real(dp)::dummy1,dummy2,ly

        !!! very rare error, if y~1 (up to machine precision) freeze it!
        if(y<0.999999999d0) then
            ly=log(1._dp-y)
        else
            ly=log(1._dp-0.999999999d0)
        end if
        listLY=(/1._dp,ly,ly**2/)



        PDF=xf(x/y,mu0,hadron)

        dummy1=sum(listLY*CplusAt1_qq)
        dummy2=sum(listLY*CplusAt1_gg)
        inter2=(/dummy1,dummy1,dummy1,dummy1,dummy1,&
        dummy2,dummy1,dummy1,dummy1,dummy1,dummy1/)

        FFplus=inter2*(PDF-PDFat1)/(1._dp-y)

    end function FFplus


    !!!!! approximate integration of y~1. It is presice up to 5-6 digits for y~0.9999, and x<0.8
    function Integrate_largey()
        real(dp),dimension(-5:5)::Integrate_largey
        real(dp)::Aqq,Aqg,Agq,Agg,Aqqp,Aqqb,PDFsum
        real(dp),dimension(-5:5)::partPLUS,partReg,PDFs
        real(dp)::lY,dummy1,dummy2
        real(dp),dimension(1:parametrizationLength):: var

        lY=log(1._dp-yCUT)

        !!!!The integrals over (..)_+ are
        !!!! int_y^1 dy (f[x/y]-f[x])/(1-y) ~ (f[x/y]-f[x])
        !!!! int_y^1 dy (f[x/y]-f[x])*log[1-y]/(1-y) ~ (f[x/y]-f[x])(Log[1-y]-1)
        !!!! int_y^1 dy (f[x/y]-f[x])*log[1-y]**2/(1-y) ~ (f[x/y]-f[x])(2+Log[1-y]*(log[1-y]-2))

        PDFs=xf(x/yCUT,mu0,hadron)

        dummy1=sum((/1._dp,lY-1._dp,2._dp+lY*(lY-2._dp)/)*CplusAt1_qq)
        dummy2=sum((/1._dp,lY-1._dp,2._dp+lY*(lY-2._dp)/)*CplusAt1_gg)
        partPLUS=(/dummy1,dummy1,dummy1,dummy1,dummy1,&
        dummy2,dummy1,dummy1,dummy1,dummy1,dummy1/)*(PDFs-PDFat1)

        !!! the integral over regular part is simpler
        !!!! it is just the value f[x]int_x^1 ...
        var=parametrizationStringAt1(yCUT)
        Aqq=sum(var*Bqq)
        Aqg=sum(var*Bqg)
        Agq=sum(var*Bgq)
        Agg=sum(var*Bgg)
        Aqqb=sum(var*Bqqb)
        Aqqp=sum(var*Bqqp)

        PDFs=xf(x,mu0,hadron)

        PDFsum=PDFs(-5)+PDFs(-4)+PDFs(-3)+PDFs(-2)+PDFs(-1)+PDFs(1)+PDFs(2)+PDFs(3)+PDFs(4)+PDFs(5)

        partReg=(/&
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

        Integrate_largey=partReg+partPLUS

    end function Integrate_largey

end function CxF_AS
