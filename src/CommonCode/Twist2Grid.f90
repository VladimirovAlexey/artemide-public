!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 1.4
!
!	This file contains the part of the code, which is common for all TMD-evaluation modules that operates at twist-2
!	It shares the variables with the module where it is inlucded (as a text)
!	Such idiotic structure is needed since, FORTRAN does not allow inheritance.
!
!	Be AWARE of possible clash of variable names.
!
!	This part is devoted to the Grid evaluation
!
!				A.Vladimirov (08.10.2018)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!Griding functions  !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


!!! this subroutine create grid. It is called only from the uTMDPDF_SetLambdaNP function.
subroutine MakeGrid()
  real*8:: x_local,b_local
  
  !!!size of grid is power of 2, to make it testable and reusable. order is the current power.
  integer:: iX,iB,j,h
  
  real*8,dimension(-5:5)::checkValue
  real*8::maxValue,eps,epsB,epsX,B_BAD,x_BAD,eps_av,time1,time2
  integer::grid_size
  real*8,dimension(0:GridSizeX,1:2,-5:5,1:numberOfHadrons)::smallGrid
  real*8::v1_BAD,v2_BAD
  
  if(numberOfHadrons==1 .and. hadronsInGRID(1)==0) return
  
  call cpu_time(time1)
  if(outputlevel>2) write(*,*) 'arTeMiDe.',moduleName,' starts to build grid. At ',time1
  
  !!!save the parameters on which the grid is build
  lambdaNP_grid=lambdaNP

  bFactor=bGrid_Max/(EXP(slope)-1d0) !!! often used combination
  slopeFactor=slope/REAL(GridSizeB)
!  goto 111 
  !!!We use the following grid
  !! b=(e^{a i/imax}-1)/(e^a-1)*b_Max 
  !! where a=slope. If slope higher more points concentrated at b->0
  !! For x we use the grid
  !! x=(xMin)^(i/imax)
  do h=1,numberOfHadrons
   do iX=0,GridSizeX
    do iB=0,GridSizeB
     x_local=xGrid_Min**(REAL(iX)/REAL(GridSizeX))
     b_local=(EXP(slopeFactor*REAL(iB))-1d0)*bFactor
     
     if(x_local==1d0) then
      gridMain(iX,iB,-5:5,h)=0d0
     else
      if(withGluon) then
       gridMain(iX,iB,-5:5,h)=x_local*Common_lowScale50(x_local,b_local,hadronsInGRID(h))
      else
       gridMain(iX,iB,-5:5,h)=x_local*Common_lowScale5(x_local,b_local,hadronsInGRID(h))
      end if
     end if
    end do    
    end do
   end do
  !!!!!!!!!!!!!**************************************************************************************************************
  !!!!!!!!!!!!!**************************************************************************************************************
  !!!!!!!!!!!!!*****************************EXTRAPOLATION PART***************************************************************
  !!!!!!!!!!!!!**************************************************************************************************************
  
  
  !!!The high B tales are extrapolated by EXP[-a-b b-c b^2](x)*(sing at bGridMax)
111 do h=1,numberOfHadrons
    do iX=1,GridSizeX
    do iB=1,2
     x_local=xGrid_Min**(REAL(iX)/REAL(GridSizeX))
     if(iB==1) then
     b_local=bGrid_Max
     else
     b_local=10*bGrid_Max
     end if
     
      if(withGluon) then
       smallGrid(iX,iB,-5:5,h)=x_local*Common_lowScale50(x_local,b_local,hadronsInGRID(h))
      else
       smallGrid(iX,iB,-5:5,h)=x_local*Common_lowScale5(x_local,b_local,hadronsInGRID(h))
       smallGrid(iX,iB,0,h)=0d0
      end if
      
      !!!save signs
      if(iB==1) then
       do j=-5,5 
        if(ABS(smallGrid(iX,iB,j,h))<1d-150) then
         taleSings(iX,j,h)=0
         smallGrid(iX,iB,j,h)=0d0 !!! sing is 0, does not metter what we save
        else
         taleSings(iX,j,h)=SIGN(1d0,smallGrid(iX,iB,j,h))
         smallGrid(iX,iB,j,h)=LOG(ABS(smallGrid(iX,iB,j,h)))
        end if
       end do
      end if
      if(iB==2) then
      do j=-5,5 
        if(taleSings(iX,j,h)==0) then
         smallGrid(iX,iB,j,h)=smallGrid(iX,1,j,h)-LOG(10d0) !!! sing is 0, does not matter what we save, to prevent errors we put some value
        elseif(ABS(smallGrid(iX,iB,j,h))<1d-150) then
         smallGrid(iX,iB,j,h)=smallGrid(iX,1,j,h)-LOG(1d100) !!! if second is zero, let just expect that it drops 1^10 times
        else
         smallGrid(iX,iB,j,h)=LOG(ABS(smallGrid(iX,iB,j,h)))
        end if
       end do
      end if
    end do
    end do
    
    SELECT CASE(asymptoticBehavior)
    CASE(1)!! exponential
    extrapolateParameters(1:GridSizeX,1,-5:5,h)=-10d0/9d0*smallGrid(1:GridSizeX,1,-5:5,h)&
		  +smallGrid(1:GridSizeX,2,-5:5,h)/9d0
    extrapolateParameters(1:GridSizeX,2,-5:5,h)=(smallGrid(1:GridSizeX,1,-5:5,h)&
		  -smallGrid(1:GridSizeX,2,-5:5,h))/bGrid_Max/9d0
    CASE(2)!! Gaussian
    extrapolateParameters(1:GridSizeX,1,-5:5,h)=-10d0/9d0*smallGrid(1:GridSizeX,1,-5:5,h)&
		  +smallGrid(1:GridSizeX,2,-5:5,h)/9d0
    extrapolateParameters(1:GridSizeX,2,-5:5,h)=(smallGrid(1:GridSizeX,1,-5:5,h)&
		  -smallGrid(1:GridSizeX,2,-5:5,h))/bGrid_Max/bGrid_Max/9d0
    CASE DEFAULT!! exponential
    extrapolateParameters(1:GridSizeX,1,-5:5,h)=-10d0/9d0*smallGrid(1:GridSizeX,1,-5:5,h)&
		  +smallGrid(1:GridSizeX,2,-5:5,h)/9d0
    extrapolateParameters(1:GridSizeX,2,-5:5,h)=(smallGrid(1:GridSizeX,1,-5:5,h)&
		  -smallGrid(1:GridSizeX,2,-5:5,h))/bGrid_Max/9d0
    END SELECT
    
    
    !!!for x=1
    extrapolateParameters(0,1:2,-5:5,h)=extrapolateParameters(1,1:2,-5:5,h)
    do j=-5,5
        taleSings(0,j,h)=0d0
    end do
    
    v1_BAD=1d0
     do iX=1,GridSizeX
     do j=-5,5
      if((.not.(j==0)).or.withGluon) then
       if(extrapolateParameters(iX,2,j,h)<0) then
        !write(*,*) 'arTeMiDe.uTMDPDF: WARNING grid diverge. Extrapolate parameters:', extrapolateParameters(iX,1:2,j)
        if(extrapolateParameters(iX,2,j,h)<v1_BAD .and. (taleSings(iX,j,h)/=0)) then
          v1_BAD=extrapolateParameters(iX,2,j,h)
          v2_BAD=taleSings(iX,j,h)*EXP(-extrapolateParameters(iX,1,j,h))
          X_BAD=xGrid_Min**(REAL(iX)/REAL(GridSizeX))
        end if
       end if
       end if
       
       if(extrapolateParameters(iX,1,j,h)+1==extrapolateParameters(iX,1,j,h) .or. &
        extrapolateParameters(iX,2,j,h)+1==extrapolateParameters(iX,2,j,h)) then !!!infinities
        write(*,*) 'arTeMiDe',moduleName,': CRITICAL. Extrapolate parameters are infinities:', extrapolateParameters(iX,1:2,j,h)&
        ,taleSings(iX,j,h)
        write(*,*) 'for hardron =', hadronsInGRID(h)
        stop
      end if
     end do
    end do
    
    if(v1_BAD<0d0 .and. outputLevel>0 .and. messageTriger<6) then
     SELECT CASE(asymptoticBehavior)
      CASE(1)!! exponential
       write(*,'("arTeMiDe.",A,": WARNING possible TMD diverge. TMD ~",ES12.3," exp(",&
          ES12.3," b) at x=",ES12.3)') moduleName,v2_BAD,-v1_BAD,x_BAD
    CASE(2)!! Gaussian
       write(*,'("arTeMiDe.",A,": WARNING possible TMD diverge. TMD ~",ES12.3," exp(",&
         ES12.3," b^2) at x=",ES12.3)') moduleName,v2_BAD,-v1_BAD,x_BAD
    CASE DEFAULT!! exponential
        write(*,'("arTeMiDe.",A,": WARNING possible TMD diverge. TMD ~",ES12.3," exp(",&
         ES12.3," b) at x=",ES12.3)') moduleName,v2_BAD,-v1_BAD,x_BAD
    END SELECT
    if(outputLevel>1) write(*,*) 'for hardron =', hadronsInGRID(h)
    messageTriger=messageTriger+1
    if(messageTriger>5) write(*,*) 'WARNING: arTeMiDe',moduleName,' number of WARNINGS more then 5. Futher WARNING suppresed'
    end if
    end do
    
    call cpu_time(time2)
    
    if(outputlevel>1) write(*,'(" ",A,": Grid build  (",I5," x",I5,")  calc.time=",F6.2,"s. ")')&
	moduleName, GridSizeX,GridSizeB, time2-time1
    
  !!check the grid
  !! we pass though grid restoring a point from near-by points
  !! in this way we estimate error
    if(outputLevel>2) then
    
    do h=1,numberOfHadrons
    epsB=0d0
    epsX=0d0
    eps_av=0d0
    grid_size=0
    !check b-grid
    do iX=0,GridSizeX
    maxValue=MAXVAL(ABS(gridMain(iX,2:GridSizeB,1:3,h)))
    if(maxValue>0d0) then
    do iB=2,GridSizeB-4
      checkValue=((-gridMain(iX,iB,-5:5,h)+4d0*gridMain(iX,iB+1,-5:5,h)+4d0*gridMain(iX,iB+3,-5:5,h)&
        -gridMain(iX,iB+4,-5:5,h))/6d0-gridMain(iX,iB+2,-5:5,h))
        
      if(.not.withGluon) checkValue(0)=0
      
      eps=MAXVAL(ABS(checkValue))/maxValue
      eps_av=eps_av+checkValue(MAXLOC(ABS(checkValue),1)-6)/maxValue      
      grid_size=grid_size+1
      
      if(eps>epsB) then
       epsB=eps
       B_BAD=(EXP(slopeFactor*REAL(iB+2))-1d0)*bFactor
       x_BAD=xGrid_Min**(REAL(iX)/REAL(GridSizeX))
       end if
    end do    
    end if
    end do
!     write(*,*) 'The worst B reconstraction with ',epsB, 'at (x,b)=',x_BAD,B_BAD
    
!   check x-grid
    do iB=0,GridSizeB
    maxValue=MAXVAL(ABS(gridMain(0:GridSizeX,iB,1:3,h)))
    if(maxValue>0d0) then
    do iX=0,GridSizeX-4
      checkValue=((-gridMain(iX,iB,-5:5,h)+4d0*gridMain(iX+1,iB,-5:5,h)+4d0*gridMain(iX+3,iB,-5:5,h)&
      -gridMain(iX+4,iB,-5:5,h))/6d0-gridMain(iX+2,iB,-5:5,h))
      
      if(.not.withGluon) checkValue(0)=0
      
      eps=MAXVAL(ABS(checkValue))/maxValue
      eps_av=eps_av+eps
      grid_size=grid_size+1
      if(eps>epsX) then
       epsX=eps
       B_BAD=(EXP(slopeFactor*REAL(iB))-1d0)*bFactor
       x_BAD=xGrid_Min**(REAL(iX+2)/REAL(GridSizeX))
       end if
    end do    
    end if
    end do
!     write(*,*) 'The worst X reconstraction with ',epsX, 'at (x,b)=',x_BAD,B_BAD
    
    write(*,'(" ",A,": Grid (for hadron ",I3,") av.badness =",ES12.3)')&
    moduleName,hadronsInGRID(h),eps_av/REAL(grid_size)
    
    end do
    end if
end subroutine MakeGrid
  
function ExtractFromGrid(x,bT,hadron)
  real*8,dimension(-5:5)::ExtractFromGrid
  real*8,dimension(0:3,-5:5):: interI
  real*8::x,bT,indexX,indexB,fX,fB
  integer::i,iX,iB,h,hadron
  real*8::var1,var2,var3,var4 !!dummyvariables
  
  !!!searching for hadron
  h=0
  do i=1,numberOfHadrons
    if(hadronsInGRID(i)==hadron) then
      h=i
      exit
    end if
  end do
  if(h==0) then
    write(*,*) 'arTeMiDe.',moduleName,': CRITICAL MISTAKE:: the hadron ',hadron,' is not found in the grid'
    write(*,*) 'arTeMiDe: evaluation STOP'
    stop
  end if
  
  if(x<xGrid_Min) then
   write(*,*) 'arTeMiDe',moduleName,': CRITICAL MISTAKE:: The TMD with x =',x,'is called. Current grid size is up to '&
   ,xGrid_Min,'. Enlarge boundaries.'
   write(*,*) 'arTeMiDe: evaluation STOP'
   stop
  end if
  if(x>1d0) then
   write(*,*) 'arTeMiDe',moduleName,': CRITICAL MISTAKE:: The TMD with x >1 (',x,') is called.' 
   write(*,*) 'arTeMiDe: evaluation STOP'
   stop
  end if
  if(bT<0d0) then
   write(*,*) 'arTeMiDe',moduleName,': CRITICAL MISTAKE:: The TMD with bT <0 (',bT,') is called.' 
   write(*,*) 'arTeMiDe: evaluation STOP'
   stop
  end if
  
  !!!!!!!!!Finding X index
  indexX=GridSizeX*LOG(x)/LOG(xGrid_Min)
  iX=INT(indexX)
    if(iX>0) then
     if(iX<GridSizeX-1) then !!normal restoration X
      iX=INT(indexX)
     else !!X in the last interval
      iX=INT(indexX)-2      
      end if
    else !! X in the first inteval
      iX=INT(indexX)+1
    end if
    fX=indexX-iX !!! fraction part (outomatically shifted by +- 1 if needed)
  
   if(bT>bGrid_Max) then
    
   !!!!here will code for large b
   SELECT CASE(asymptoticBehavior)
   CASE(1)!! exponential
    interI(0,-5:5)=taleSings(iX-1,-5:5,h)*EXP(-extrapolateParameters(iX-1,1,-5:5,h)&
	  -extrapolateParameters(iX-1,2,-5:5,h)*bT)
    interI(1,-5:5)=taleSings(iX,-5:5,h)*EXP(-extrapolateParameters(iX,1,-5:5,h)&
	  -extrapolateParameters(iX,2,-5:5,h)*bT)
    interI(2,-5:5)=taleSings(iX+1,-5:5,h)*EXP(-extrapolateParameters(iX+1,1,-5:5,h)&
	  -extrapolateParameters(iX+1,2,-5:5,h)*bT) 
    interI(3,-5:5)=taleSings(iX+2,-5:5,h)*EXP(-extrapolateParameters(iX+2,1,-5:5,h)&
	  -extrapolateParameters(iX+2,2,-5:5,h)*bT)
   CASE(2)!! Gaussian
    interI(0,-5:5)=taleSings(iX-1,-5:5,h)*EXP(-extrapolateParameters(iX-1,1,-5:5,h)&
	  -extrapolateParameters(iX-1,2,-5:5,h)*bT*bT)
    interI(1,-5:5)=taleSings(iX,-5:5,h)*EXP(-extrapolateParameters(iX,1,-5:5,h)&
	  -extrapolateParameters(iX,2,-5:5,h)*bT*bT)
    interI(2,-5:5)=taleSings(iX+1,-5:5,h)*EXP(-extrapolateParameters(iX+1,1,-5:5,h)&
	  -extrapolateParameters(iX+1,2,-5:5,h)*bT*bT) 
    interI(3,-5:5)=taleSings(iX+2,-5:5,h)*EXP(-extrapolateParameters(iX+2,1,-5:5,h)&
	  -extrapolateParameters(iX+2,2,-5:5,h)*bT*bT)
   CASE DEFAULT!! exponential
    interI(0,-5:5)=taleSings(iX-1,-5:5,h)*EXP(-extrapolateParameters(iX-1,1,-5:5,h)&
	  -extrapolateParameters(iX-1,2,-5:5,h)*bT)
    interI(1,-5:5)=taleSings(iX,-5:5,h)*EXP(-extrapolateParameters(iX,1,-5:5,h)&
	  -extrapolateParameters(iX,2,-5:5,h)*bT)
    interI(2,-5:5)=taleSings(iX+1,-5:5,h)*EXP(-extrapolateParameters(iX+1,1,-5:5,h)&
	  -extrapolateParameters(iX+1,2,-5:5,h)*bT) 
    interI(3,-5:5)=taleSings(iX+2,-5:5,h)*EXP(-extrapolateParameters(iX+2,1,-5:5,h)&
	  -extrapolateParameters(iX+2,2,-5:5,h)*bT)
   END SELECT
   
   else!!! b inside the main region
    
    indexB=LOG(1d0+bT/bFactor)/slopeFactor
    iB=INT(indexB)
    if(iB>0) then
     if(iB<GridSizeB-1) then !!normal restoration B
      iB=INT(indexB)
     else !!B in the last interval
      iB=INT(indexB)-2      
      end if
    else !! B in the first inteval
      iB=INT(indexB)+1
    end if
    fB=indexB-iB !!! fraction part (automatically shifted by +- 1 if needed)
     
     
     
     !! intepolation procedure
        var1=-fB*(fB-1d0)*(fB-2d0)
        var2=3d0*(fB+1d0)*(fB-1d0)*(fB-2d0)
        var3=-3d0*(fB+1d0)*fB*(fB-2d0)
        var4=(fB+1d0)*fB*(fB-1d0)
        interI(0,-5:5)=(var1*gridMain(iX-1,iB-1,-5:5,h)+var2*gridMain(iX-1,iB,-5:5,h)&
	      +var3*gridMain(iX-1,iB+1,-5:5,h)+var4*gridMain(iX-1,iB+2,-5:5,h))/6d0
	interI(1,-5:5)=(var1*gridMain(iX,iB-1,-5:5,h)+var2*gridMain(iX,iB,-5:5,h)&
	      +var3*gridMain(iX,iB+1,-5:5,h)+var4*gridMain(iX,iB+2,-5:5,h))/6d0
	interI(2,-5:5)=(var1*gridMain(iX+1,iB-1,-5:5,h)+var2*gridMain(iX+1,iB,-5:5,h)&
	      +var3*gridMain(iX+1,iB+1,-5:5,h)+var4*gridMain(iX+1,iB+2,-5:5,h))/6d0
	interI(3,-5:5)=(var1*gridMain(iX+2,iB-1,-5:5,h)+var2*gridMain(iX+2,iB,-5:5,h)&
	      +var3*gridMain(iX+2,iB+1,-5:5,h)+var4*gridMain(iX+2,iB+2,-5:5,h))/6d0
	      
     !!! linear interpolation procedure
! 	interI(0,-5:5)=fB*gridMain(iX-1,iB,-5:5,h)+(1d0-fB)*gridMain(iX-1,iB+1,-5:5,h)
! 	interI(1,-5:5)=fB*gridMain(iX,iB,-5:5,h)+(1d0-fB)*gridMain(iX,iB+1,-5:5,h)
! 	interI(2,-5:5)=fB*gridMain(iX+1,iB,-5:5,h)+(1d0-fB)*gridMain(iX+1,iB+1,-5:5,h)
! 	interI(3,-5:5)=fB*gridMain(iX+2,iB,-5:5,h)+(1d0-fB)*gridMain(iX+2,iB+1,-5:5,h)
	 
     end if
	
    var1=-fX*(fX-1d0)*(fX-2d0)
    var2=3d0*(fX+1d0)*(fX-1d0)*(fX-2d0)
    var3=-3d0*(fX+1d0)*fX*(fX-2d0)
    var4=(fX+1d0)*fX*(fX-1d0)
    ExtractFromGrid=(var1*interI(0,-5:5)+var2*interI(1,-5:5)&
			  +var3*interI(2,-5:5)+var4*interI(3,-5:5))/6d0/x
  do i=-5,5
   if(ISNAN(ExtractFromGrid(i))) then
   
    write(*,*) 'here',bT,i,ExtractFromGrid(i)
    
    write(*,*) 'inter',interI(0:3,i)
    
    write(*,*) 'extrap',taleSings(iX-1,i,h),extrapolateParameters(iX-1,1:2,i,h)
    
    stop
   end if
  end do
			  
end function ExtractFromGrid 