!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!
!!!!!!!!!       part of SiversTMDPDF module for artemide
!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!     The twist-3 convolution with in Sivers function is difficult enough!     !!!!!!!!!!!
!!!!!!!!!!!!!!!!!     So the NP-functions are made independent on convolution library          !!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
  
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
  
