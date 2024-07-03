!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!Program makes a smooth plot for AS-term
!
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


module KTtest
INCLUDE '/data/arTeMiDe_Repository/artemide/src/Code/KTspace/Fourier_Levin.f90'
end module KTtest


! module KTtest2
! INCLUDE '/data/arTeMiDe_Repository/artemide/src/Code/KTspace/grid_inKT(new).f90'
! end module KTtest2

program example
use KTtest

implicit none

call Initialize_Fourier_Levin("Prog/Tests/const-files/ART23_BM.atmde",'*14  ','*F   ',"JAJA",3,1)

call TestFourier()

end program example
