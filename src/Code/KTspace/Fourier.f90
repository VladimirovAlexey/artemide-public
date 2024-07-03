!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 3.0
!
!	This file contains the part of the code, which is common for all TMD modules.
!       It is included (as a text).
!	Such idiotic structure is needed since, FORTRAN does not allow inheritance.
!
!	Be AWARE of possible clash of variable names.
!
!	This part is devoted to the calculation of Fourier transform into kT-space.
!
!	v.3.00 Created (AV, 04.02.2024)
!
!				A.Vladimirov (04.02.2024)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!! the container file must include in the header
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! integer,parameter::TMDtypeN=0 !!!!! this is the order of Bessel-transform (IT IS STRICT FOR TMD)
! real(dp)::kT_FREEZE=0.0001_dp  !!!!! parameter of freezing the low-kT-value
!
! !----Ogata Tables---
! integer,parameter::Nmax=1000
! INCLUDE 'Tables/BesselZero1000.f90'
!
! logical:: convergenceLost=.false.
!
! !!!!! I split the qT over runs qT<qTSegmentationBoundary
! !!!!! In each segment I have the ogata quadrature with h=hOGATA*hSegmentationWeight
! !!!!! It helps to convergen integrals, since h(optimal) ~ qT
! integer,parameter::hSegmentationNumber=7
! real(dp),dimension(1:hSegmentationNumber),parameter::hSegmentationWeight=(/0.0001d0,0.001d0,0.01d0,1d0,2d0,5d0,10d0/)
! real(dp),dimension(1:hSegmentationNumber),parameter::qTSegmentationBoundary=(/0.001d0,0.01d0,0.1d0,10d0,50d0,100d0,200d0/)
!
! real(dp)::hOGATA,toleranceOGATA
! !!!weights of ogata quadrature
! real(dp),dimension(1:hSegmentationNumber,1:Nmax)::ww
! !!!nodes of ogata quadrature
! real(dp),dimension(1:hSegmentationNumber,1:Nmax)::bb
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!Prepare tables for Ogata quadrature with given h
!!! note that the factor 1/(2pi) is taken into ww
!!! the difference between definition in TMDF and here is 1/pi
subroutine PrepareTables()
    integer::i,j
    real(dp)::hS!=h*hSegmentationWeight
    real(dp)::xi
    real(dp)::CommonPrefactor

    if(outputlevel>2) write(*,'(A,A)',advance="no") moduleName,": Preparing Ogata tables for kt-transformation ..."

    !!!! common prefactor is M^{2n}/n! as it is defined for TMD distributions
    SELECT CASE(TMDtypeN)
        CASE(0)
            CommonPrefactor=1._dp
        CASE(1)
            CommonPrefactor=TMDmass**2
        CASE(2)
            CommonPrefactor=TMDmass**4/2
        CASE DEFAULT
            write(*,*) ErrorString('Unknown type ot Momentum transformation for TMD',moduleName),TMDtypeN
            stop
    END SELECT

    do j=1,hSegmentationNumber
    do i=1,Nmax

    hS=hOGATA*hSegmentationWeight(j)
    xi=JZero(TMDtypeN,i)

    !     ww(j,k,i)=BESSEL_JN(k,bb(j,k,i))/xi/(BESSEL_JN(k+1,xi)**2)&
    !         *(pi*xi*hS*Cosh(xi*hS)+Sinh(pi*Sinh(xi*hS)))/(1d0+Cosh(pi*Sinh(xi*hS)))

    !!! if we too far away in xI*hS, the double exponential grow rapidly.
    !!! and for >6, it generates term 10^{300} and exceed the presision

    !!! The CommonPrefactor is the part of TMD transformation definition

    if(xi*hS>6.d0) then
        bb(j,i)=xi*Tanh(piHalf*Sinh(xi*hS))
        ww(j,i)=CommonPrefactor*BESSEL_JN(TMDtypeN,bb(j,i))/xi/(BESSEL_JN(TMDtypeN+1,xi)**2)/pi

    else
        bb(j,i)=xi*Tanh(piHalf*Sinh(xi*hS))
        ww(j,i)=CommonPrefactor*BESSEL_JN(TMDtypeN,bb(j,i))/xi/(BESSEL_JN(TMDtypeN+1,xi)**2)&
        *(pi*xi*hS*Cosh(xi*hS)/(2d0*Cosh(piHalf*Sinh(xi*hS))**2)+Tanh(piHalf*Sinh(xi*hS)))/pi
    end if

    end do
    end do

    if(outputlevel>2) write(*,*) color(' done',c_green)

end subroutine PrepareTables


!------------------------------------------FOURIER--------------------------------
!!!This is the defining module function
!!! It evaluates the integral
!!!  int_0^infty   b db/2pi  J_num(b qT) F1  (b/qT)^num M^{2num}/num!
!!! the prefactor M^{2num}/num! is included into the weights of Ogata ww
function Fourier_ev(x,qT_in,mu,zeta,hadron)
    real(dp),intent(in)::x,mu,zeta,qT_in
    integer,intent(in)::hadron
    real(dp)::integral(-5:5),eps(-5:5),qT
    real(dp)::v1(-5:5),v2(-5:5),v3(-5:5),v4(-5:5),delta(-5:5)
    logical:: partDone(-5:5)
    integer::k,j,Nsegment
    real(dp)::Fourier_ev(-5:5)

    integral=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)

    if(qT_in<kT_FREEZE) then
        qT=kT_FREEZE
    else
        qT=qT_in
    end if

    v1=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
    v2=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
    v3=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
    v4=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
    partDone=(/.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false./)

    !!! define segment of qT
    do j=1,hSegmentationNumber
        if(qT<qTSegmentationBoundary(j)) exit
    end do
    if(j>hSegmentationNumber) then
        Nsegment=hSegmentationNumber
    else
        Nsegment=j
    end if

    do k=1,Nmax!!! maximum of number of bessel roots preevaluated in the head
        eps=ww(Nsegment,k)*(bb(Nsegment,k)**(TMDtypeN+1))*TMD_ev(x,bb(Nsegment,k)/qT,mu,zeta,hadron)

        v4=v3
        v3=v2
        v2=v1
        v1=abs(eps)

        delta=(v1+v2+v3+v4)
        integral=integral+eps

        !!! here we check that residual term is smaller than already collected integral
        !!! also checking the zerothness of the integral. If already collected integral is null it is null
        !!! Here is potential bug. If the first 10 points give zero (whereas some later points do not), the integral will be zero
        !!! I check for each separate flavor
        do j=-5,5
            if((delta(j)<toleranceOGATA*ABS(integral(j)) .or. ABS(integral(j))<1d-32) .and. k>=10) partDone(j)=.true.
        end do
        if(partDone(-5).and.partDone(-4).and.partDone(-3).and.partDone(-2).and.partDone(-1)&
            .and.partDone(0).and.partDone(1).and.partDone(2).and.partDone(3).and.partDone(4).and.partDone(5)) exit

    end do

    if(k>=Nmax) then
        if(outputlevel>0) call Warning_Raise('OGATA quadrature diverge. TMD decaing too slow?',&
            messageCounter,messageTrigger,moduleName)
            if(outputlevel>2) then
            write(*,*) 'Information over the last call ----------'
            write(*,*) partDone
            write(*,*) 'bt/qT= ',bb(Nsegment,Nmax)/qT, 'qT=',qT, '| segmentation zone=',Nsegment,&
                ' ogata h=',hOGATA*hSegmentationWeight(Nsegment)
            write(*,*) 'W=',TMD_ev(x,bb(Nsegment,k)/qT,mu,zeta,hadron), 'eps/integral =', eps/integral
            write(*,*) 'v1+v2+v3+v4=',v1+v2+v3+v4, '>',toleranceOGATA*(ABS(integral(1))+ABS(integral(2)))
            write(*,*) '------------------------------------------'
            end if
    end if
    !!! result is scaled by qT [because the argument of Bessel was scaled bqT-> B]
    !!! the extra factor 1/k^n appears due to the defenition of Fourier-trnaform for TMD-n
    Fourier_ev=integral/(qT**(2*TMDtypeN+2))

end function Fourier_ev

!!! It evaluates the integral
!!!  int_0^infty   b db/2pi  J_num(b qT) F1  (b/qT)^num M^{2num}/num!
!!! the prefactor M^{2num}/num! is included into the weights of Ogata ww
function Fourier_opt(x,qT_in,hadron)
    real(dp),intent(in)::x,qT_in
    integer,intent(in)::hadron
    real(dp)::integral(-5:5),eps(-5:5),qT
    real(dp)::v1(-5:5),v2(-5:5),v3(-5:5),v4(-5:5),delta(-5:5)
    logical:: partDone(-5:5)
    integer::k,j,Nsegment
    real(dp)::Fourier_opt(-5:5)

    integral=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)

    if(qT_in<kT_FREEZE) then
        qT=kT_FREEZE
    else
        qT=qT_in
    end if

    v1=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
    v2=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
    v3=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
    v4=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
    partDone=(/.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false./)

    !!! define segment of qT
    do j=1,hSegmentationNumber
        if(qT<qTSegmentationBoundary(j)) exit
    end do
    if(j>hSegmentationNumber) then
        Nsegment=hSegmentationNumber
    else
        Nsegment=j
    end if

    do k=1,Nmax!!! maximum of number of bessel roots preevaluated in the head
        eps=ww(Nsegment,k)*(bb(Nsegment,k)**(TMDtypeN+1))*TMD_opt(x,bb(Nsegment,k)/qT,hadron)

        v4=v3
        v3=v2
        v2=v1
        v1=abs(eps)

        delta=(v1+v2+v3+v4)
        integral=integral+eps

        !!! here we check that residual term is smaller than already collected integral
        !!! also checking the zerothness of the integral. If already collected integral is null it is null
        !!! Here is potential bug. If the first 10 points give zero (whereas some later points do not), the integral will be zero
        !!! I check for each separate flavor
        do j=-5,5
            if((delta(j)<toleranceOGATA*ABS(integral(j)) .or. ABS(integral(j))<1d-32) .and. k>=10) partDone(j)=.true.
        end do
        if(partDone(-5).and.partDone(-4).and.partDone(-3).and.partDone(-2).and.partDone(-1)&
            .and.partDone(0).and.partDone(1).and.partDone(2).and.partDone(3).and.partDone(4).and.partDone(5)) exit

    end do

    if(k>=Nmax) then
        if(outputlevel>0) call Warning_Raise('OGATA quadrature diverge. TMD decaing too slow?',&
            messageCounter,messageTrigger,moduleName)
            if(outputlevel>2) then
            write(*,*) 'Information over the last call ----------'
            write(*,*) partDone
            write(*,*) 'bt/qT= ',bb(Nsegment,Nmax)/qT, 'qT=',qT, '| segmentation zone=',Nsegment,&
                ' ogata h=',hOGATA*hSegmentationWeight(Nsegment)
            write(*,*) 'W=',TMD_opt(x,bb(Nsegment,k)/qT,hadron), 'eps/integral =', eps/integral
            write(*,*) 'v1+v2+v3+v4=',v1+v2+v3+v4, '>',toleranceOGATA*(ABS(integral(1))+ABS(integral(2)))
            write(*,*) '------------------------------------------'
            end if
    end if
    !!! result is scaled by qT [because the argument of Bessel was scaled bqT-> B]
    !!! the extra factor 1/k^n appears due to the defenition of Fourier-trnaform for TMD-n
    Fourier_opt=integral/(qT**(2*TMDtypeN+2))

end function Fourier_opt
