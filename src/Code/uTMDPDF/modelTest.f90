!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!
!!!!!!!!!       part of uTMDPDF module for artemide
!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!-------------------------------------------------------------------------------------------------------
!!!-----------------------------Routines that test the behaviour of NP model------------------------------
!!!-------------------------------------------------------------------------------------------------------
!!! test FNP for z-dependance
!!! the test consists in the evaluation of FNP at several random sets and NParray
!!! and comparison of the values.
function TestFNP()
    logical::TestFNP
    real(dp)::xR,bR
    real(dp),dimension(-5:5)::test1,test2
    integer::h,i,j
    real(dp),allocatable::NPtest(:)

    allocate(NPtest(1:lambdaNPlength))

    call RANDOM_NUMBER(NPtest)
    NPtest=3d0*NPtest

    TestFNP=.false.
    do h=1,numberOfHadrons
    do i=1,10
        !! generate some random input
        call RANDOM_NUMBER(xR)
        if(xR>0.99d0) xR=xR/2d0
        if(xR<0.00001d0) xR=0.0001d0+xR

        call RANDOM_NUMBER(bR)  
        bR=3d0*bR

        test1=FNP(xR,0.9d0,bR,hadronsInGRID(h),NPtest)
        do j=1,10
            test2=FNP(xR,0.8d0/REAL(j),bR,hadronsInGRID(h),NPtest)
            if(SUM(ABS(test1-test2))>1d-10) then
                TestFNP=.true.
                exit
            end if	
            end do
    end do
    end do

    deallocate(NPtest)
end function TestFNP


  
!!! test MU for x-dependance
!!! the test consists in the evaluation of FNP at several random sets and NParray
!!! and comparison of the values.
function TestMU()
    logical::TestMU
    real(dp)::xR,bR
    real(dp)::test1,test2
    integer::i,j
    TestMU=.false.
    do i=1,10
        call RANDOM_NUMBER(bR)  
        bR=3d0*bR
            !!! generate some random input
        call RANDOM_NUMBER(xR)
        if(xR>0.99d0) xR=xR/2d0
        if(xR<0.00001d0) xR=0.0001d0+xR
        test1=mu_OPE(xR,bR)
            
        !!! generate some random input
        call RANDOM_NUMBER(xR)
        if(xR>0.99d0) xR=xR/2d0
        if(xR<0.00001d0) xR=0.0001d0+xR
        test2=mu_OPE(xR,bR)

        if(ABS(test1-test2)>1d-10) then
            TestMU=.true.
            exit
        end if	
    end do
end function TestMU
  
  
!!! test bSTAR for lambda-dependance
!!! .true. = lambda-dependent
function TestbSTAR()
    integer::i,j,k
    logical::TestbSTAR
    real(dp)::bT,dummy1,dummy2
    real(dp),allocatable::NPtest(:)

    allocate(NPtest(1:lambdaNPlength))
    TestbSTAR=.false.
    do i=1,lambdaNPlength
        NPtest=0.5d0+0d0*NPtest
        do k=0,3
            bT=0.1d0+k
            dummy1=bSTAR(bT,NPtest)

            do j=1,3
                NPtest(i)=0.9d0*j
                dummy2=bSTAR(bT,NPtest)
                if(abs(dummy2-dummy1)>1d-10) then
                    TestbSTAR=.true.
                    exit
                end if
            end do
        end do
    end do

end function TestbSTAR
  
