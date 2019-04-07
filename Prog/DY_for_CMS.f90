!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! Calculation of DY cross-section for CMS.
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program xSec_DY
use aTMDe_control
use TMDX_DY
implicit none

!!! this program evaluate cross-section for SIDIS

real*8 :: time1, time2
integer :: dummyINT,j
real*8 :: dummyREAL
logical::exist

integer:: process,ss
real*8:: Qmin,Qmax,s,ymin,ymax,ptCut,ptCut2
integer::s1,s2,s3,s4 !!! TMD sizes
real*8,allocatable, dimension(:) :: pt, xSec1,xSec2,xSec3,xSec4
real*8,allocatable, dimension(:) :: xSec11,xSec21,xSec31,xSec41!!!xSec for variation of c2+
real*8,allocatable, dimension(:) :: xSec12,xSec22,xSec32,xSec42!!!xSec for variation of c2-
real*8,allocatable, dimension(:) :: xSec13,xSec23,xSec33,xSec43!!!xSec for variation of c4+
real*8,allocatable, dimension(:) :: xSec14,xSec24,xSec34,xSec44!!!xSec for variation of c4-
! real*8,dimension(1:17)::ptBINS=(/0.1d0,2d0, 4d0, 6d0, 8d0, 10d0, &
!       12d0, 14d0, 18d0,22d0, 28d0, 37d0, 52d0, 85d0, 160d0, 240d0, 1000d0/)

real*8,dimension(1:17)::ptBINS=(/0.1d0,2d0, 4d0, 6d0, 8d0, 10d0, &
       13d0, 16d0, 20d0,25d0, 30d0, 37d0, 45d0, 55d0, 65d0, 75d0, 1000d0/)

call cpu_time(time1)

process=1!!

!This is from the artemide ver.1.4
!   call TMDX_DY_Initialize("NNLO")
!   call TMDX_DY_SetNPParameters((/3.3235d0,0.0380d0,0.2204d0, 7.0808d0,351.7950d0, 2.4632d0,-3.8334d0,0.0001d0, 0.0000d0/))

!This is for artemide ver.2.00
  call artemide_Initialize('constants-file','Models/BSV19.bFIT/')
  call artemide_SetReplica_uTMDPDF(0)
  call artemide_SetReplica_TMDR(0)

     
  s=13000d0**2

  ymin=-2.4d0
  ymax=2.4d0
  ptCut=25d0
  ptCut2=20d0
     
   call TMDX_DY_SetProcess(process)
   !call SetCuts(.true.,ptCut,ptCut2,yMin,yMax)
   call TMDX_DY_SetCuts(.true.,ptCut,ptCut2,yMin,yMax)
   call TMDX_DY_XSetup(s,91d0,0d0)
   
   !---------------------------------------------------------CENTRAL VALUE
   
  Qmin=50d0
  Qmax=76d0
  
   do j=1,17
   if(0.2d0*Qmax<ptBINS(j)) exit
   end do
   s1=j-1
   
   allocate(pt(1:j))
   allocate(xSec1(1:j-1),xSec11(1:j-1),xSec12(1:j-1),xSec13(1:j-1),xSec14(1:j-1))
   pt=ptBINS(1:j)
   ss=j-1
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec1,pt,Qmin,Qmax,ymin,ymax)
   xSec1=xSec1/(pt(2:ss+1)-pt(1:ss))
   
   deallocate(pt)
   
   Qmin=76d0
   Qmax=106d0   
   do j=1,17
   if(0.2d0*Qmax<ptBINS(j)) exit
   end do
   s2=j-1
   
   allocate(pt(1:j))
   allocate(xSec2(1:j-1),xSec21(1:j-1),xSec22(1:j-1),xSec23(1:j-1),xSec24(1:j-1))
   pt=ptBINS(1:j)
   ss=j-1
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec2,pt,Qmin,Qmax,ymin,ymax)
   xSec2=xSec2/(pt(2:ss+1)-pt(1:ss))
   
   deallocate(pt)
   
   Qmin=106d0
   Qmax=170d0
   do j=1,17
   if(0.2d0*Qmax<ptBINS(j)) exit
   end do
   s3=j-1
   
   allocate(pt(1:j))
   allocate(xSec3(1:j-1),xSec31(1:j-1),xSec32(1:j-1),xSec33(1:j-1),xSec34(1:j-1))
   pt=ptBINS(1:j)
   ss=j-1
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec3,pt,Qmin,Qmax,ymin,ymax)
   xSec3=xSec3/(pt(2:ss+1)-pt(1:ss))
   
   deallocate(pt)
   
   
   Qmin=170d0
   Qmax=350d0
   do j=1,17
   if(0.2d0*Qmax<ptBINS(j)) exit
   end do
   s4=j-1
   
   allocate(pt(1:j))
   allocate(xSec4(1:j-1),xSec41(1:j-1),xSec42(1:j-1),xSec43(1:j-1),xSec44(1:j-1))
   pt=ptBINS(1:j)
   ss=j-1
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec4,pt,Qmin,Qmax,ymin,ymax)
   xSec4=xSec4/(pt(2:ss+1)-pt(1:ss))
   
   deallocate(pt)
   
   !-----------------------------------------------------------------
   
  !---------------------------------------------------------Variation of c2+
  
!   call TMDX_DY_SetScaleVariations(1d0,2d0,1d0,1d0)
  call artemide_SetScaleVariations(1d0,2d0,1d0,1d0)
   
  Qmin=50d0
  Qmax=76d0
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec11,ptBINS(1:s1+1),Qmin,Qmax,ymin,ymax)
   xSec11=xSec11/(ptBins(2:s1+1)-ptBins(1:s1))
   
   Qmin=76d0
   Qmax=106d0
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec21,ptBINS(1:s2+1),Qmin,Qmax,ymin,ymax)
   xSec21=xSec21/(ptBins(2:s2+1)-ptBins(1:s2))
   
   Qmin=106d0
   Qmax=170d0
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec31,ptBINS(1:s3+1),Qmin,Qmax,ymin,ymax)
   xSec31=xSec31/(ptBins(2:s3+1)-ptBins(1:s3))
   
   
   Qmin=170d0
   Qmax=350d0
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec41,ptBINS(1:s4+1),Qmin,Qmax,ymin,ymax)
   xSec41=xSec41/(ptBins(2:s4+1)-ptBins(1:s4))
   
   !-----------------------------------------------------------------
     !---------------------------------------------------------Variation of c2-
  
!   call TMDX_DY_SetScaleVariations(1d0,0.5d0,1d0,1d0)
  call artemide_SetScaleVariations(1d0,0.5d0,1d0,1d0)
   
  Qmin=50d0
  Qmax=76d0
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec12,ptBINS(1:s1+1),Qmin,Qmax,ymin,ymax)
   xSec12=xSec12/(ptBins(2:s1+1)-ptBins(1:s1))
   
   Qmin=76d0
   Qmax=106d0
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec22,ptBINS(1:s2+1),Qmin,Qmax,ymin,ymax)
   xSec22=xSec22/(ptBins(2:s2+1)-ptBins(1:s2))
   
   Qmin=106d0
   Qmax=170d0
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec32,ptBINS(1:s3+1),Qmin,Qmax,ymin,ymax)
   xSec32=xSec32/(ptBins(2:s3+1)-ptBins(1:s3))
   
   
   Qmin=170d0
   Qmax=350d0
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec42,ptBINS(1:s4+1),Qmin,Qmax,ymin,ymax)
   xSec42=xSec42/(ptBins(2:s4+1)-ptBins(1:s4))
   
  !-----------------------------------------------------------------
  !---------------------------------------------------------Variation of c4+
  
!   call TMDX_DY_SetScaleVariations(1d0,1d0,1d0,2d0)
  call artemide_SetScaleVariations(1d0,1d0,1d0,2d0)
   
  Qmin=50d0
  Qmax=76d0
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec13,ptBINS(1:s1+1),Qmin,Qmax,ymin,ymax)
   xSec13=xSec13/(ptBins(2:s1+1)-ptBins(1:s1))
   
   Qmin=76d0
   Qmax=106d0
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec23,ptBINS(1:s2+1),Qmin,Qmax,ymin,ymax)
   xSec23=xSec23/(ptBins(2:s2+1)-ptBins(1:s2))
   
   Qmin=106d0
   Qmax=170d0
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec33,ptBINS(1:s3+1),Qmin,Qmax,ymin,ymax)
   xSec33=xSec33/(ptBins(2:s3+1)-ptBins(1:s3))
   
   
   Qmin=170d0
   Qmax=350d0
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec43,ptBINS(1:s4+1),Qmin,Qmax,ymin,ymax)
   xSec43=xSec43/(ptBins(2:s4+1)-ptBins(1:s4))
   
   !-----------------------------------------------------------------
   !---------------------------------------------------------Variation of c4-
  
!   call TMDX_DY_SetScaleVariations(1d0,1d0,1d0,0.5d0)
  call artemide_SetScaleVariations(1d0,1d0,1d0,0.5d0)
   
  Qmin=50d0
  Qmax=76d0
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec14,ptBINS(1:s1+1),Qmin,Qmax,ymin,ymax)
   xSec14=xSec14/(ptBins(2:s1+1)-ptBins(1:s1))
   
   Qmin=76d0
   Qmax=106d0
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec24,ptBINS(1:s2+1),Qmin,Qmax,ymin,ymax)
   xSec24=xSec24/(ptBins(2:s2+1)-ptBins(1:s2))
   
   Qmin=106d0
   Qmax=170d0
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec34,ptBINS(1:s3+1),Qmin,Qmax,ymin,ymax)
   xSec34=xSec34/(ptBins(2:s3+1)-ptBins(1:s3))
   
   
   Qmin=170d0
   Qmax=350d0
   
   call CalcXsec_DY_PTint_Qint_Yint(xSec44,ptBINS(1:s4+1),Qmin,Qmax,ymin,ymax)
   xSec44=xSec44/(ptBins(2:s4+1)-ptBins(1:s4))
   
   !-----------------------------------------------------------------
   
   write(*,*) ' sqrt(s) = ',Sqrt(s)
   write(*,*) ' y in (',ymin,',',ymax,')'
   write(*,*) ' pt(cut1) = ',ptCut
   write(*,*) ' pt(cut2) = ',ptCut2
   
   write(*,*) ' '
   
   
   write(*,*) '50<Q< 76:  number of TMD bins ',size(xSec1)
   write(*,*) 'pt(min),   pt(max),    xSec,    err-,     err+'
   do j=1,size(xSec1)
   write(*,'(F6.2,"     ", F6.2,"     ", ES12.3E3,"     ", ES12.3E3,"     ", ES12.3E3)') ptBINS(j),ptBINS(j+1),xSec1(j),&
    min(xSec1(j),xSec11(j),xSec12(j),xSec13(j),xSec14(j)),&
    max(xSec1(j),xSec11(j),xSec12(j),xSec13(j),xSec14(j))
   end do
   
   write(*,*) "  "
   
   write(*,*) '76<Q< 106:  number of TMD bins ',size(xSec2)
   write(*,*) 'pt(min),   pt(max),    xSec,    err-,     err+'
   do j=1,size(xSec2)
   write(*,'(F6.2,"     ", F6.2,"     ", ES12.3E3,"     ", ES12.3E3,"     ", ES12.3E3)') ptBINS(j),ptBINS(j+1),xSec2(j),&
    min(xSec2(j),xSec21(j),xSec22(j),xSec23(j),xSec24(j)),&
    max(xSec2(j),xSec21(j),xSec22(j),xSec23(j),xSec24(j))
   end do
   
   write(*,*) "  "
   
   write(*,*) '106<Q< 170:  number of TMD bins ',size(xSec3)
   write(*,*) 'pt(min),   pt(max),    xSec,    err-,     err+'
   do j=1,size(xSec3)
   write(*,'(F6.2,"     ", F6.2,"     ", ES12.3E3,"     ", ES12.3E3,"     ", ES12.3E3)') ptBINS(j),ptBINS(j+1),xSec3(j),&
    min(xSec3(j),xSec31(j),xSec32(j),xSec33(j),xSec34(j)),&
    max(xSec3(j),xSec31(j),xSec32(j),xSec33(j),xSec34(j))
   end do
   
   write(*,*) "  "
   
   write(*,*) '170<Q< 350:  number of TMD bins ',size(xSec4)
   write(*,*) 'pt(min),   pt(max),    xSec,    err-,     err+'
   do j=1,size(xSec4)
   write(*,'(F6.2,"     ", F6.2,"     ", ES12.3E3,"     ", ES12.3E3,"     ", ES12.3E3)') ptBINS(j),ptBINS(j+1),xSec4(j),&
    min(xSec4(j),xSec41(j),xSec42(j),xSec43(j),xSec44(j)),&
    max(xSec4(j),xSec41(j),xSec42(j),xSec43(j),xSec44(j))
   end do
   
call cpu_time(time2)
! 
write(*,*) 'Calculation time=',time2-time1
! 
   
 end program xSec_DY