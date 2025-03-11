!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 3.0
!
!	This file contains the part of the code, which is common for all TMD-OPE modules
!       that operates at twist-2. It is inclucded (as a text).
!	Such idiotic structure is needed since, FORTRAN does not allow inheritance.
!
!	Be AWARE of possible clash of variable names.
!
!	This part collect the definition of expressions used in the large-X resummation
!
!	v.3.00 Created (AV, 19.07.2024)
!
!				A.Vladimirov (19.07.2024)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!! The large-X resummation leds to the following expression
!!! (\delta(1-x)-A/(1-x)^(1+A)) Exp[EXP]
!!! where A and EXP are perturbative expressions
!!! in this module create these expressions from the known part


!!!!!Argument of the exponent
!!! QUARK channel
pure function LargeX_EXP_q_q(alpha,Nf,Lmu)
  real(dp),intent(in)::alpha,Lmu
  integer,intent(in)::Nf
  real(dp)::LargeX_EXP_q_q
  real(dp)::dd1,dd2

  LargeX_EXP_q_q=0.d0

  !!! the coefficients of delta-function are extracted from coefficient files
  if(orderLX>=1) then
    dd1=C_q_q_delta_1(Nf,Lmu)
    LargeX_EXP_q_q=LargeX_EXP_q_q+alpha*dd1

  if(orderLX>=2) then
    dd2=C_q_q_delta_2(Nf,Lmu)
    LargeX_EXP_q_q=LargeX_EXP_q_q+alpha*alpha*(2*dd2-dd1**2)/2.d0

  if(orderLX>=3) then
    LargeX_EXP_q_q=LargeX_EXP_q_q+alpha**3*(3*C_q_q_delta_3(Nf,Lmu)-3*dd2*dd1+dd1**3)/3.d0
  end if
  end if
  end if
end function LargeX_EXP_q_q

!!! GLUON channel
pure function LargeX_EXP_g_g(alpha,Nf,Lmu)
  real(dp),intent(in)::alpha,Lmu
  integer,intent(in)::Nf
  real(dp)::LargeX_EXP_g_g
  real(dp)::dd1,dd2

  LargeX_EXP_g_g=0.d0

  !!! the coefficients of delta-function are extracted from coefficient files
  if(orderLX>=1) then
    dd1=C_g_g_delta_1(Nf,Lmu)
    LargeX_EXP_g_g=LargeX_EXP_g_g+alpha*dd1

  if(orderLX>=2) then
    dd2=C_g_g_delta_2(Nf,Lmu)
    LargeX_EXP_g_g=LargeX_EXP_g_g+alpha*alpha*(2*dd2-dd1**2)/2.d0

  if(orderLX>=3) then
    LargeX_EXP_g_g=LargeX_EXP_g_g+alpha**3*(3*C_g_g_delta_3(Nf,Lmu)-3*dd2*dd1+dd1**3)/3.d0
  end if
  end if
  end if
end function LargeX_EXP_g_g

!!!!!Argument of the exponentiation = 2 D(b*,muOPE)
!!! QUARK channel
pure function LargeX_A_q_q(alpha,Nf,Lmu)
  real(dp),intent(in)::alpha,Lmu
  integer,intent(in)::Nf
  real(dp)::LargeX_A_q_q

  LargeX_A_q_q=2*Dpert_atL(alpha,Nf,Lmu,orderLX,1)

end function LargeX_A_q_q

!!!!!Argument of the exponentiation = 2 D(b*,muOPE)
!!! GLUON channel
pure function LargeX_A_g_g(alpha,Nf,Lmu)
  real(dp),intent(in)::alpha,Lmu
  integer,intent(in)::Nf
  real(dp)::LargeX_A_g_g

  LargeX_A_g_g=2*Dpert_atL(alpha,Nf,Lmu,orderLX,0)

end function LargeX_A_g_g

!!!!!-------------------------- This function is copied and corrected from Twist2Convolution -------------

!!!! The CxF is defined as Mellin convolution of [C xF]
!!!! The integral is \int_x^1 dy C(y) F[x/y], where F[x]=x PDF(x)
!!!! The ordinary convolution \int_x^1 dy/y C(y) PDF(x/y) = CxF(x)/x
!!!! the parameters x,y for model are defined as in this formula

!!!! the computation with large-X resummation has different structure
!!!! The delta-term exponentiates and mutiply the (..)_+ part
function CxF_largeX_compute(x,bT,hadron,includeGluon)
    real(dp),dimension(-5:5)::CxF_largeX_compute
    integer, intent(in)::hadron
    real(dp),intent(in)::x,bT
    logical,intent(in)::includeGluon
    real(dp):: bTcurrent,lx

    real(dp),dimension(-5:5)::deltaPart,PLUSremnant
    real(dp):: muAt1,asAt1,LogAt1,Cqq,Cgg,Csingqq,Csinggg,alpha_powerQ_at1,alpha_powerG_at1
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
        CxF_largeX_compute=0._dp
        return
    end if

    !!! values of parameters at y=1
    !!! they are used also later
    muAt1=muOPE(bTcurrent,x,1._dp,c4_global)
    asAt1=As(muAt1)
    LogAt1=LogMuB(bTcurrent,x,1._dp)
    NfAt1=activeNf(muAt1)
    PDFat1=xf(x,muAt1,hadron)
    alpha_powerQ_at1=LargeX_A_q_q(asAt1,NfAt1,LogAt1)
    if(includeGluon) then
      alpha_powerG_at1=LargeX_A_g_g(asAt1,NfAt1,LogAt1)
    else
      alpha_powerG_at1=0._dp
    end if

!     write(*,*) "--->",LogAt1,alpha_powerQ_at1,exp(LargeX_EXP_q_q(asAt1,NfAt1,LogAt1))
!     write(*,*) "--->",PDFat1(1)

    if(alpha_powerQ_at1>0.9_dp) then
      write(*,*) ErrorString('The powerQ of large-X resummation is bigger than 0.9',moduleName)
      write(*,*) '----- information on last call -----'
      write(*,*) 'x=', x, 'bT=',bT, 'hadron=',hadron,' power=',alpha_powerQ_at1
      ERROR STOP ErrorString('The powerQ of large-X resummation is bigger than 0.9',moduleName)
    end if
    if(alpha_powerG_at1>0.9_dp) then
      write(*,*) ErrorString('The powerG of large-X resummation is bigger than 0.9',moduleName)
      write(*,*) '----- information on last call -----'
      write(*,*) 'x=', x, 'bT=',bT, 'hadron=',hadron,' power=',alpha_powerG_at1
      ERROR STOP ErrorString('The powerG of large-X resummation is bigger than 0.9',moduleName)
    end if

    !!!! the computation with large-X resummation has different structure
    !!!! if there is no large-X resummation all terms are additive
    !!!! if there is large-X resummation them delta-term exponentiates and mutiply the (..)_+ part

    !!!! delta-part in the resummed case is exp(....)
    Cqq=exp(LargeX_EXP_q_q(asAt1,NfAt1,LogAt1))
    if(includeGluon) then
        Cgg=exp(LargeX_EXP_g_g(asAt1,NfAt1,LogAt1))
    else
        Cgg=0._dp
    end if
    deltaPart=(/Cqq,Cqq,Cqq,Cqq,Cqq,Cgg,Cqq,Cqq,Cqq,Cqq,Cqq/)


    !!! it is expected that orderLX>=orderMain
    !!!! other parts contribute only if order >LO
    if(orderLX>0) then

      !!!! the ''remnant term'' is f(x)/(1-x)^a

      Csingqq=(1-x)**(-alpha_powerQ_at1)
      if(includeGluon) then
        Csinggg=(1-x)**(-alpha_powerG_at1)
      else
        Csinggg=0._dp
      end if

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

    !!!!! minus signs are because of definition alpha-coefficient
    CxF_largeX_compute=deltaPart*(PLUSremnant-Integrate_GK_array5(FFplus,x,yCUT,toleranceINT)-Integrate_largey_PLUS())&
            +Integrate_GK_array5(FFreg,x,yCUT,toleranceINT)+Integrate_largey_Reg()

    else
      !!!!! order =0
      CxF_largeX_compute=PDFat1

    end if

  do i=-5,5
   if(ISNAN(CxF_largeX_compute(i))) then

    write(*,*) ErrorString('convolution computed to NAN. CHECK INTEGRATION',moduleName)
    write(*,*) '----- information on last call -----'
    write(*,*) 'x=', x, 'bT=',bT,' i=',i, 'hadron=',hadron,' result=',CxF_largeX_compute(i)

   end if
  end do

!     !!!!! these are wrapper functions to pass ther integrand to GK routine
!     !!!!! FORTRAN gets only function of one variable
!     !!!!! using the trick with internal function one can by-pass this limitation
    contains

    !!! regular integrand in the same as in the usual convolution
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
            Agq=sum(var*Bgq)
            Agg=sum(var*Bgg)
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
    !!! this part is int [alpha(y)f(x/y)/(1-y)^{1+alpha(y)}-alpha(1)f(x)/(1-y)^{1+alpha(1)}]
    function FFplus(y)
        real(dp),dimension(-5:5)::FFplus
        real(dp),intent(in)::y
        real(dp),dimension(-5:5)::PDF,inter1,inter2
        real(dp)::muCurrent,asCurrent,LogCurrent,dummy1,dummy2,alpha_inter
        integer::NfCurrent

        !!!! if mu is y-dependent one needs to update the values each step
        if(IsMuYdependent) then
            muCurrent=muOPE(bTcurrent,x,y,c4_global)
            asCurrent=As(muCurrent)
            LogCurrent=LogMuB(bTcurrent,x,y)
            NfCurrent=activeNf(muCurrent)

            !!!! expression at y
            alpha_inter=LargeX_A_q_q(asCurrent,NfCurrent,LogCurrent)
            dummy1=alpha_inter/(1-y)**(1+alpha_inter)

            if(includeGluon) then
                alpha_inter=LargeX_A_g_g(asCurrent,NfCurrent,LogCurrent)
                dummy2=alpha_inter/(1-y)**(1+alpha_inter)
            else
                dummy2=0._dp
            end if
            PDF=xf(x/y,muCurrent,hadron)

            inter1=(/dummy1,dummy1,dummy1,dummy1,dummy1,&
            dummy2,dummy1,dummy1,dummy1,dummy1,dummy1/)

            !!!! expression at y=1
            dummy1=alpha_powerQ_at1/(1-y)**(1+alpha_powerQ_at1)

            if(includeGluon) then
                dummy2=alpha_powerG_at1/(1-y)**(1+alpha_powerG_at1)
            else
                dummy2=0._dp
            end if

            inter2=(/dummy1,dummy1,dummy1,dummy1,dummy1,&
            dummy2,dummy1,dummy1,dummy1,dummy1,dummy1/)

            FFplus=(inter1*PDF-inter2*PDFat1)
        else

            !!!! expression at y=1
            dummy1=alpha_powerQ_at1/(1-y)**(1+alpha_powerQ_at1)

            if(includeGluon) then
                dummy2=alpha_powerG_at1/(1-y)**(1+alpha_powerG_at1)
            else
                dummy2=0._dp
            end if

            PDF=xf(x/y,muAt1,hadron)

            inter2=(/dummy1,dummy1,dummy1,dummy1,dummy1,&
            dummy2,dummy1,dummy1,dummy1,dummy1,dummy1/)

            FFplus=inter2*(PDF-PDFat1)
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
        real(dp)::muCurrent,asCurrent,LogCurrent,dummy1,dummy2,alpha_inter
        integer::NfCurrent
        real(dp),dimension(-5:5)::PDFs,inter1,inter2

        !!!!The integrals over (..)_+ is
        if(IsMuYdependent) then
            !!! if mu is y-dependent one needs to update the values of parameters for each y
            muCurrent=muOPE(bTcurrent,x,yCUT,c4_global)
            asCurrent=As(muCurrent)
            LogCurrent=LogMuB(bTcurrent,x,yCUT)
            NfCurrent=activeNf(muCurrent)

            PDFs=xf(x/yCUT,muCurrent,hadron)

            !!!! expression at y
            alpha_inter=LargeX_A_q_q(asCurrent,NfCurrent,LogCurrent)
            dummy1=alpha_inter/(1-yCUT)**(alpha_inter)/(1-alpha_inter)

            if(includeGluon) then
                alpha_inter=LargeX_A_g_g(asCurrent,NfCurrent,LogCurrent)
                dummy2=alpha_inter/(1-yCUT)**(1+alpha_inter)/(1-alpha_inter)
            else
                dummy2=0._dp
            end if
            PDFs=xf(x/yCUT,muCurrent,hadron)

            inter1=(/dummy1,dummy1,dummy1,dummy1,dummy1,&
            dummy2,dummy1,dummy1,dummy1,dummy1,dummy1/)

            !!!! expression at y=1
            dummy1=alpha_powerQ_at1/(1-yCUT)**(alpha_powerQ_at1)/(1-alpha_inter)

            if(includeGluon) then
                dummy2=alpha_powerG_at1/(1-yCUT)**(alpha_powerG_at1)/(1-alpha_inter)
            else
                dummy2=0._dp
            end if

            inter2=(/dummy1,dummy1,dummy1,dummy1,dummy1,&
            dummy2,dummy1,dummy1,dummy1,dummy1,dummy1/)

            Integrate_largey_PLUS=(inter1*PDFs-inter2*PDFat1)
        else
            PDFs=xf(x/yCUT,muAt1,hadron)

            !!!! expression at y=1
            dummy1=alpha_powerQ_at1/(1-yCUT)**(alpha_powerQ_at1)/(1-alpha_inter)

            if(includeGluon) then
                dummy2=alpha_powerG_at1/(1-yCUT)**(alpha_powerG_at1)/(1-alpha_inter)
            else
                dummy2=0._dp
            end if

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
        Agq=sum(var*Bgq)
        Agg=sum(var*Bgg)
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

end function CxF_largeX_compute
