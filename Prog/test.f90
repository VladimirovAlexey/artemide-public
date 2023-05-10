program example
use aTMDe_control
use TMDX_DY
implicit none

 real*8,dimension(1:3)::xSec,pt,ptmin,ptmax
 real*8,dimension(1:4)::ptPLUS
 integer :: j
 real*8 :: time1, time2,xx
 real*8,dimension(1:2):: p1,p2,p3
 real*8,dimension(1:4)::cc
 integer,dimension(1:3)::pp
 
 real*8,dimension(1:3,1:2)::qtList,Qlist,yList
 real*8,dimension(1:3,1:4)::cutList
 logical,dimension(1:3)::inCutList
 integer,dimension(1:3,1:3)::pList
 integer,dimension(1:3)::nnList
 real*8,dimension(1:3)::X,sList
 !$  real*8::OMP_get_wtime,t1,t2
 
 call cpu_time(time1)
 !$  t1=OMP_get_wtime()
 
  write(*,*) "Initialize artemide at LO. It is fast ~1 sek. NNLO could be much longer ~5 min."
  call artemide_Initialize('const-test')
!   call TMDX_DY_Initialize('LO')
  ! call TMDX_DY_setProcess(1,1,6)
  ! call SetCuts(.true.,20d0,-2.4d0,2.4d0)
  ! call TMDX_DY_XSetup(8000d0**2,91d0,0.01d0)
!   call TMDX_DY_SetNPParameters(0)
  !call artemide_SetReplica_uTMDPDF(0)
  !call artemide_SetReplica_TMDR(0)
 
  do j=1,3
  pt(j)=2*REAL(j)
  ptmin(j)=pt(j)-1d0
  ptPLUS(j)=ptmin(j)
  ptmax(j)=pt(j)+1d0
  end do
  ptPLUS(4)=ptMax(3)
  
  
!  call CalcXsec_DY_PTint_Qint_Yint(xSec,ptmin,ptmax,66d0,126d0,-0.5d0,3.5d0)
!   do j=1,3
!   write(*,*) ptmin(j),'--',ptmax(j),xSec(j)
!   end do 
!  
  write(*,*) "Calculating some values for cross-section one-by-one (DY around Z-boson peak, ATLAS 8TeV kinematics)"
  write(*,*) "ptMin 	--	ptMax		xSec"
   do j=1,3
!     write(*,*) ptmin(j),'--',ptmax(j),xSec_DY(process=(/1,1,6/),s=8000d0**2,&
! 		qT=(/ptMin(j),ptmax(j)/),Q=(/66d0,126d0/),y=(/-0.5d0,3.5d0/),&
! 		includeCuts=.true.,CutParameters=(/20d0,20d0,-2.4d0,2.4d0/),Num=4)
      p1=(/ptMin(j),ptmax(j)/)
      p2=(/66d0,116d0/)
      p3=(/-2.4d0,2.4d0/)
      cc=(/20d0,20d0,-2.4d0,2.4d0/)
      pp=(/1,1,5/)
      call xSec_DY(xx,pp,(8000d0)**2,p1,p2,p3,.true.,CutParameters=cc)
      write(*,*) ptmin(j),'--',ptmax(j),xx
   end do
   
   write(*,*) " "
   write(*,*) "Now the same by list"
   !$ write(*,*) "It must be faster since you use OPENMP"
   do j=1,3
    qtList(j,:)=(/ptMin(j),ptmax(j)/)
    Qlist(j,:)=(/66d0,116d0/)
    yList(j,:)=(/-2.4d0,2.4d0/)
    inCutList(j)=.true.
    cutList(j,:)=(/20d0,20d0,-2.4d0,2.4d0/)
    pList(j,:)=(/1,1,5/)
    nnList(j)=4
    sList(j)=(8000d0)**2
   end do
   
   call xSec_Dy_List(X,pList,sList,qtList,Qlist,yList,inCutList,cutList,nnList)
   
   write(*,*) "result:", X
   
   write(*,*) '-----------------------------------------------------------------------------------------------'
  call cpu_time(time2)
  !$ t2=OMP_get_wtime()
  write(*,*) 'The programm evaluation took', time2-time1, 'sec. (CPU time)'
  !$  write(*,*) 'The programm evaluation took', t2-t1, 'sec. (wallclock time)'
  write(*,*) '  '
  write(*,*) 'If you do not like so many terminal messages check the parameter outputlevel in constants file.'
  write(*,*) 'Not forget to cite artemide [1706.01473]' 
  write(*,*) '-----------------------------------------------------------------------------------------------'

end program example
