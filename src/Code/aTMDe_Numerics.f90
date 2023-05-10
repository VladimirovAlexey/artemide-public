!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!			arTeMiDe 2.03
!
!	Contains the definition of the precision and global numeric constants, used in computation
!   Used in each artemide module.
!
!				A.Vladimirov (10.02.2020)
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

module aTMDe_Numerics

implicit none

public
!!-------------------------Definition of presision kinds----------------
!!-------------------------artemide uses dp presision everywhere -------
!!-------------------------declare as real(dp)
integer, parameter :: sp = selected_real_kind(6, 37)
integer, parameter :: dp = selected_real_kind(15, 307)
integer, parameter :: qp = selected_real_kind(33, 4931)
!!! for the complex nimbers I use COMPLEX*16

!!------------------------Constants
!!------------------------written up to 64 digits (in case double-precision will be rised)
!!! Zeta-functions
real(dp),parameter::zeta2=1.644934066848226436472415166646025189218949901206798437735558229_dp
real(dp),parameter::zeta3=1.202056903159594285399738161511449990764986292340498881792271555_dp
real(dp),parameter::zeta4=1.082323233711138191516003696541167902774750951918726907682976215_dp
real(dp),parameter::zeta5=1.036927755143369926331365486457034168057080919501912811974192678_dp
real(dp),parameter::zeta6=1.017343061984449139714517929790920527901817490032853561842408664_dp

!!! Pi 's
real(dp),parameter::piHalf=1.570796326794896619231321691639751442098584699687552910487472296_dp    !!pi/2
real(dp),parameter::pi=3.141592653589793238462643383279502884197169399375105820974944592_dp        !!pi
real(dp),parameter::pix2=6.283185307179586476925286766559005768394338798750211641949889185_dp      !!2*pi
real(dp),parameter::pix4=12.56637061435917295385057353311801153678867759750042328389977837_dp      !!4*pi

!!! Pi^2 's
real(dp),parameter::pi2=9.869604401089358618834490999876151135313699407240790626413349376_dp        !!pi^2
real(dp),parameter::pi2x2=19.73920880217871723766898199975230227062739881448158125282669875_dp      !!2*pi^2
real(dp),parameter::pi2x4=39.47841760435743447533796399950460454125479762896316250565339750_dp      !!4*pi^2

!!! Pi^4 's
real(dp),parameter::pi4=97.40909103400243723644033268870511124972758567268542169146785939_dp        !!pi^4
real(dp),parameter::pi4x2=194.8181820680048744728806653774102224994551713453708433829357188_dp      !!2*pi^4
real(dp),parameter::pi4x4=389.6363641360097489457613307548204449989103426907416867658714376_dp      !!4*pi^4

!!! C0 -constant (typical TMD constant)
real(dp), parameter :: C0_const=1.122918967133770339648286429581761573531420773850306336308318152_dp        !!=2Exp[-gamma_E]
real(dp), parameter :: C0_inv_const=0.8905362089950989926182520515535897745848226071517151026788329383_dp    !!=Exp[gamma_E]/2

end module aTMDe_Numerics
