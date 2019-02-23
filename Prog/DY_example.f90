!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!  Program that evaluate the cross-section for DY 
!!  Surves as an example of ardemide-based code
!!  Please, refer to 1902.???? if you use it.
!!
!!  Question, suggestions, etc: vladimirov.aleksey@gmail.com
!!  Extra information & details: https://github.com/VladimirovAlexey/artemide-public/blob/master/doc/Manual.pdf
!!
!!	MADE ON ARTEMIDE 1.4 
!!		(compatibility with other versions in not guarantied)
!!
!!			A.Vladimirov: 13.02.2019
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

program DYtest
!!load the artemide library related to DY cross-section
use TMDX_DY
implicit none

real*8::s,Qmin,Qmax,ymin,ymax,qtMin(1:5),qtMax(1:5),pTcut,etaCutMin,etaCutMax,xSec(1:5),xx
integer::process(1:3),i


!!The first step: initialize the artemide
!!		it is a complex routine that runs thorugh modules, and read all initial inputs and prepare all variables
!!		it could be called only ones (further calls will be ignored)
!!		argument is "LO", "NLO", "NNLO" which specifies the order (see details in the manual)
call TMDX_DY_Initialize("NNLO")

write(*,*) "-----------------------------------------------------------------------------------------------------------------"
write(*,*) "---- if you do not need all these loading messages change the output level in the begining of constants-file ----"
write(*,*) "-----------------------------------------------------------------------------------------------------------------"

!!The second step: state the parameters of the non-perturbative model for TMDs
!!		it can be done directly like -->   call TMDX_DY_SetNPParameters((/3.3235d0,0.0380d0,0.2204d0, 7.0808d0,351.7950d0, 2.4632d0,-3.8334d0,0.0001d0, 0.0000d0/))
!!		or using provided set (this is prefered option if you are not going to make a fit of TMDs)
!!		In this case just state the number of the replica you would like to use (0=central) --> call TMDX_DY_SetNPParameters(0)
!!
!!		The model is defined in the set of files in src/Model/ (see manual for details). 
!!		Some extra models are collected in /Models (change of the model module requares recompilation of the artemide)
!!
!!		IMPORTANT: if the preparation of TMDs grid is ON (check it in "constants"-file option *3*B
!!			then these grids will be calculated after evaluation of this command.
!!			At NNLO it can take significant time (~3-5 min depending on computer), so set "LO" if you would like just to test the installation.
!!			Switch off the grid option if you need to evaluate several points, otherwice using the grids significantly speed up the calculation.
call TMDX_DY_SetNPParameters(0)

!!----------------------------------------------
!!Now the artemide is ready for calculation!
!!----------------------------------------------
!!There are many parameters and specifications for xSec, including different phase space elements, fiducial cuts, bin-integrations, etc.
!!Artemide allows to calculate xSec in various combination/configurations. Check manual (sec on TMDX_DY) for details.
!!Naturally, it does not include all possibilities, but many. If your task requares some not included case, ask me -- I will be happy to include it into the code.

!!---------------------------------------------
!!Example of xSec calculation 
!!---------------------------------------------
!! input: 
!! process pp->gamma/Z 
!! s=(8TeV)^2
!! 66<Q<116 GeV
!! 0.8<|y|<1.2
!! qt-bins: 0-2, 2-4, 4-6, 6-8, 8-10 GeV
!! fiducial cuts on detected lepton pair:
!! pT>20 GeV, -2.4<eta<2.4
!!
!! It corresponds to one of y-bins in ATLAS measurement at 8 TeV, presented in 1512.02192

process=(/1,1,5/)!! this is pp->gamma/Z
s=8000d0**2
Qmin=66d0
Qmax=116d0
ymin=0.8d0
ymax=1.2d0
qtMin=(/0.01d0,2d0,4d0,6d0,8d0/) !! for convergence it is better not to set qt=0 exactly in the presentce of fiducial cuts.Anyway cross-section ->0 at qt->0
qtMax=(/2d0,4d0,6d0,8d0,10d0/)
pTcut=20d0
etaCutMin=-2.4d0
etaCutMax=2.4d0


!!There are two equivalent interfaces to evaluation of DY xsec.
!!	1) Specify the process, kinematics, cuts first and then ask for xSec.
!!	2) The same in a single command with plenty of arguments.
!! Here is example of usage of both commands

!!---------------
!!1) step-by-step
!! set process
call TMDX_DY_SetProcess(process)
!! specify cut parameters (call with .false. if not fiducial cuts are needed)
call SetCuts(.true.,PTcut,etaCutMin,etaCutMax)
!! specify kinematic domain (s,Q,y) (since we going to integrate over Q and y, Q and y arguments play no role)
call TMDX_DY_XSetup(s,Qmin,ymin)
!! calculate xSec
!! here xSec must be array of the same length as qt-list.
call CalcXsec_DY_PTint_Qint_Yint(xSec,qtMin,qtMax,Qmin,Qmax,ymin,ymax)
!! Note that the calculate xSec is in the region 0.8<y<1.2. For the region 0.8<|y|<1.2 multiply it by 2 (since it is symmetric in y->-y)
!! also the data is avaraged in bit (while here is integrated), so we devide by a size of the bin
xSec=2d0*xSec/(qtMax-qtMin)
!! Here is the result!
write(*,*) "evaluation 1            : ",xSec

!!---------------
!!2) This way is simpler (and faster if you have many cores and compile artemde with openmp)
!! there is no need to specify all details independently. There is also a listable version of this function
   do i=1,5
      call xSec_DY(xx,process,s,(/qtMin(i),qtMax(i)/),(/Qmin,Qmax/),(/ymin,ymax/),.true.,(/pTcut,pTcut,etaCutMin,etaCutMax/))
      xSec(i)=2d0*xx/(qtMax(i)-qtMin(i))
   end do
!! Here is the result!
write(*,*) "evaluation 2            : ",xSec

!! results of both evaluations should be equal
!! compare it to actual numbers
write(*,*) "values from [1512.02192]: ", (/3.112d0,6.462d0,6.548d0,5.647d0,4.716d0/)

end program DYtest