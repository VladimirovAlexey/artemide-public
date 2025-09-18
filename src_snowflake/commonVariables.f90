!!!! Working precision of variables
integer, parameter :: dp = selected_real_kind(15, 307)

!!!!! type for twist-3 function
abstract interface
    function func_tw3(x1,x2)
        import::dp
        real(dp) :: func_tw3
        real(dp), intent(in) ::x1,x2
    end function func_tw3
end interface

!!!! Print initialization information
logical::showINI=.false.
!!!! Print initialization information
logical::showPROCESS=.false.

!!!! minimal accesable x (r=0)
real(dp)::xMIN=0.01_dp

!!!! number of nodes in rrho = NR
integer::NUM_RHO=8
!!!! number of nodes in phi = Mphi
integer::NUM_PHI=8
!!!! total size of grid in RHO
integer::NUM_TOT_RHO=8
!!!! total size of grid in PHI
integer::NUM_TOT_PHI=47
!!!! total size of grid in 1D
integer::NUM_TOT=22

!!!! Value of zero to compare with
real(dp)::zero=1d-12
!!!! Integration tolerance
real(dp)::toleranceINT=1d-8
!!!! Runge-Kutta maximal step
real(dp)::RGh_max=0.001_dp
!!!! Number of processesor allows for parralel computation
integer::allowedNumProcessor=16

!!!! Setup chiral even evolution
logical:: IncludeChiralEvenEvolution=.true.
!!!! Setup chiral odd evolution
logical:: IncludeChiralOddEvolution=.true.
!!!! Include mixing with gluon, other flavors etc.
logical:: useSingletEvolution=.true.
!!!! Setup G2-computation
logical:: IncludeG2Matrix=.true.
!!!! Setup WGT-computation
logical:: IncludeWGTMatrix=.true.

!!!! Mass of the CHARM threshold [GeV]
real(dp):: massCHARM=1.27_dp
!!!! Mass of the BOTTOM threshold [GeV]
real(dp):: massBOTTOM=4.18_dp
