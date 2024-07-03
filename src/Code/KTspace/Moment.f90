!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 3.0
!
!	This file contains the part of the code, which is common for all TMD modules.
!       It is included (as a text).
!	Such idiotic structure is needed since, FORTRAN does not allow inheritance.
!
!	Be AWARE of possible clash of variable names.
!
!	This part is devoted to the calculation of TMM.
!   If you use this module please cite [2402.01836]
!
!	v.3.00 Created (AV, 11.02.2024)
!
!				A.Vladimirov (11.02.2024)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!! the container file must include in the header
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
! real(dp)::muTMM_min=0.8_dp  !!!!! minimal mu
!
! !!!!! I split the qT over runs qT<qTSegmentationBoundary
! !!!!! For TMM this split is the same as for inKT
!
! real(dp)::hOGATA_TMM,toleranceOGATA_TMM
! !!!weights of ogata quadrature
! real(dp),dimension(1:hSegmentationNumber,0:3,1:Nmax)::ww_TMM
! !!!nodes of ogata quadrature
! real(dp),dimension(1:hSegmentationNumber,0:3,1:Nmax)::bb_TMM
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!!Prepare tables for Ogata quadrature with given h
!!! note that the factor 1/(2pi) is taken into ww
!!! the difference between definition in TMDF and here is 1/pi
subroutine PrepareTablesTMM()
    integer::i,j,k
    real(dp)::hS!=h*hSegmentationWeight
    real(dp)::xi

    if(outputlevel>2) write(*,'(A,A)',advance="no") moduleName,": Preparing Ogata tables for TMM-transformation ..."

    do j=1,hSegmentationNumber
    do k=0,3
    do i=1,Nmax

    hS=hOGATA_TMM*hSegmentationWeight(j)
    xi=JZero(k,i)


    !!! if we too far away in xI*hS, the double exponential grow rapidly.
    !!! and for >6, it generates term 10^{300} and exceed the presision

    if(xi*hS>6.d0) then
        bb_TMM(j,k,i)=xi*Tanh(piHalf*Sinh(xi*hS))
        ww_TMM(j,k,i)=BESSEL_JN(k,bb_TMM(j,k,i))/xi/(BESSEL_JN(k+1,xi)**2)/pi

    else
        bb_TMM(j,k,i)=xi*Tanh(piHalf*Sinh(xi*hS))
        ww_TMM(j,k,i)=BESSEL_JN(k,bb_TMM(j,k,i))/xi/(BESSEL_JN(k+1,xi)**2)&
        *(pi*xi*hS*Cosh(xi*hS)/(2d0*Cosh(piHalf*Sinh(xi*hS))**2)+Tanh(piHalf*Sinh(xi*hS)))/pi
    end if

    end do
    end do
    end do

    if(outputlevel>2) write(*,*) color(' done',c_green)

    if(outputLevel>2) write(*,*) color("If you use TMMs, please, cite [2402.01836]",c_cyan)

end subroutine PrepareTablesTMM


!------------------------------------------General Moment--------------------------------
!!! This is the general function for the computation of the Moment
!!! It evaluates the integral
!!! int_0^infty   b^n db   J_k(b mu) TMD_opt[x,b]
function Moment_Gen(n,k,x,mu,hadron)
    real(dp),intent(in)::x,mu
    integer,intent(in)::n,k,hadron
    real(dp),dimension(-5:5)::Moment_Gen(-5:5)
    real(dp),dimension(-5:5)::integral,eps,v1,v2,v3,v4,delta
    logical:: partDone(-5:5)
    integer::r,j,Nsegment

    integral=(/0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0,0d0/)

    if(mu<muTMM_min) then
        write(*,*) ErrorString("ERROR in KT-moment computation. mu<mu_min:",moduleName),muTMM_min
        stop
    end if

    if(k>3) then
        write(*,*) ErrorString("ERROR in KT-moment computation. Called moment with J_k, k>3",moduleName)
        stop
    end if

    if(k+n<0) then
        write(*,*) ErrorString("ERROR in KT-moment computation. Integral is divergent at 0",moduleName)
        stop
    end if

    v1=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
    v2=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
    v3=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
    v4=(/1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0,1d0/)
    partDone=(/.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false.,.false./)

    !!! define segment of qT
    do j=1,hSegmentationNumber
        if(mu<qTSegmentationBoundary(j)) exit
    end do
    if(j>hSegmentationNumber) then
        Nsegment=hSegmentationNumber
    else
        Nsegment=j
    end if

    do r=1,Nmax!!! maximum of number of bessel roots preevaluated in the head
        eps=ww_TMM(Nsegment,k,r)*(bb_TMM(Nsegment,k,r)**n)*TMD_opt(x,bb_TMM(Nsegment,k,r)/mu,hadron)

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
            if((delta(j)<toleranceOGATA_TMM*ABS(integral(j)) .or. ABS(integral(j))<1d-32) .and. r>=10) partDone(j)=.true.
        end do
        if(partDone(-5).and.partDone(-4).and.partDone(-3).and.partDone(-2).and.partDone(-1)&
            .and.partDone(0).and.partDone(1).and.partDone(2).and.partDone(3).and.partDone(4).and.partDone(5)) exit

    end do

    if(r>=Nmax) then
        if(outputlevel>0) call Warning_Raise('OGATA quadrature diverge for TMM. TMD decaing too slow?',&
            messageCounter,messageTrigger,moduleName)
            if(outputlevel>2) then
            write(*,*) 'Information over the last call ----------'
            write(*,*) partDone
            write(*,*) 'bt/mu= ',bb_TMM(Nsegment,k,Nmax)/mu, 'qT=',mu, '| segmentation zone=',Nsegment,&
                ' ogata h=',hOGATA_TMM*hSegmentationWeight(Nsegment)
            write(*,*) 'W=',TMD_opt(x,bb_TMM(Nsegment,k,Nmax)/mu,hadron), 'eps/integral =', eps/integral
            write(*,*) 'v1+v2+v3+v4=',v1+v2+v3+v4, '>',toleranceOGATA_TMM*(ABS(integral(1))+ABS(integral(2)))
            write(*,*) 'x=',x,'n=',n,'type of J =',k
            write(*,*) '------------------------------------------'
            end if
    end if

    !!! result is multiplied by (2pi) [because the weights are defined for db/2pi] and scaled by mu [by definition of the moment]
    Moment_Gen=integral/(mu**(n+1))*pix2

end function Moment_Gen

!!!! The moment G=G_{n,n} in [2402.01836], defined as
!!!! mu^{n+1}/2^n/n! int b^n J_{n+1} (b mu) F[OPTIMAL]
function Moment_G(x,mu,hadron)
    real(dp),intent(in)::x,mu
    integer,intent(in)::hadron
    real(dp),dimension(-5:5)::Moment_G(-5:5)

    SELECT CASE(TMDtypeN)
        CASE(0)
            Moment_G=mu*Moment_Gen(0,1,x,mu,hadron)
        CASE(1)
            Moment_G=mu**2/2*Moment_Gen(1,2,x,mu,hadron)
        CASE(2)
            Moment_G=mu**3/8*Moment_Gen(2,3,x,mu,hadron)
        CASE DEFAULT
            write(*,*) ErrorString("MOMENT G is defined only for n=0,1,2",moduleName)
            stop
    END SELECT

end function Moment_G

!!!! The moment X=G_{n+1,n} in [2402.01836], defined as
!!!! mu^{n+3}/2^n/n!/(2M^2) int b^n ((n+1)J_{n+1} (b mu)-J_{n+3} (b mu))/(n+2)
!!!! =
!!!! mu^{n+3}/2^n/n!/(2M^2) int b^n (J_{n+1} (b mu)-J_{n+2}(b mu)/(b mu))
!!!! I use the last formula because it utilizes lower order of Bessel-function
function Moment_X(x,mu,hadron)
    real(dp),intent(in)::x,mu
    integer,intent(in)::hadron
    real(dp),dimension(-5:5)::Moment_X(-5:5)

    SELECT CASE(TMDtypeN)
        CASE(0)
            Moment_X=mu**2/(2*TMDmass**2)*(mu*Moment_Gen(0,1,x,mu,hadron)-2*Moment_Gen(-1,2,x,mu,hadron))
            !Moment_X=mu**3/(2*2*TMDmass**2)*(Moment_Gen(0,1,x,mu,hadron)-Moment_Gen(0,3,x,mu,hadron))
        CASE(1)
            Moment_X=mu**3/(4*TMDmass**2)*(mu*Moment_Gen(1,2,x,mu,hadron)-2*Moment_Gen(0,3,x,mu,hadron))
        CASE DEFAULT
            write(*,*) ErrorString("MOMENT X is defined only for n=0,1",moduleName)
            stop
    END SELECT

end function Moment_X
