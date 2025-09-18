!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
!!! List of functions to substitute into the cross-section
!!! These functions contain product of TMD distributions and all what is neccesary to attach
!!!
function TMD_pair(Q2,x1,x2,k1,k2,mu,process)
real(dp),intent(in)::Q2,x1,x2,k1,k2,mu
integer,dimension(1:3),intent(in)::process
real(dp)::TMD_pair
real(dp),dimension(-5:5)::FA,FB,FAB
real(dp)::param
integer::h1,h2

!increment counter
GlobalCounter=GlobalCounter+1
LocalCounter=LocalCounter+1

h1=process(1)
h2=process(2)


SELECT CASE(process(3))
  !!!test cases
  CASE(0,10000,20000,30000)
    !param=0.1d0
    !TMD_pair=(param/(k1+param**2)**1.5)*(param/(k2+param**2)**1.5)/(pix4)
    param=2.5d0
    TMD_pair=(3*param*(2*param**2-3*k1)/(k1+param**2)**3.5)*(3*param*(2*param**2-3*k2)/(k2+param**2)**3.5)/(pix4)
  CASE(9999,19999,29999,39999)
    TMD_pair=(Exp(-0.2d0*k1)+1/(k1+2.))*(Exp(-0.2d0*k2)+1/(k2+2.))
  CASE(9998,19998,29998,39998)
    TMD_pair=(Exp(-0.2d0*k1)+1/(k1+2.))*(Exp(-0.2d0*k2)+1/(k2+2.))

  CASE (1) !pp->gamma
    ! e_q^2 *F_q(A)*F_qbar(B)
     FA=uTMDPDF_inKT(x1,sqrt(k1),mu,h1)
     FB=uTMDPDF_inKT(x2,sqrt(k2),mu,-h2) !!! -h2, to myltiply quarks by anti-quarks in FAB
     FAB=FA*FB

    TMD_pair=FAB(1)/9.d0&
      +FAB(2)*4.d0/9.d0&
      +FAB(3)/9.d0&
      +FAB(4)*4d0/9.d0&
      +FAB(5)/9d0&
      +FAB(-1)/9.d0&
      +FAB(-2)*4.d0/9.d0&
      +FAB(-3)/9.d0&
      +FAB(-4)*4d0/9.d0&
      +FAB(-5)/9d0

 !--------------------------------------------------------------------------------
  CASE (2) !Delta^{GG'}z_{+l}z_{+f}f1f1 (ONLY Z-BOSON without gamma)

     FA=uTMDPDF_inKT(x1,sqrt(k1),mu,h1)
     FB=uTMDPDF_inKT(x2,sqrt(k2),mu,-h2)!!! -h2, to multiply quarks by anti-quarks in FAB
     FAB=FA*FB

     TMD_pair=XTMD_pairZZ(FAB,Q2)

 !--------------------------------------------------------------------------------
  CASE (3,20,21,22,29) !Delta^{GG'}z_{+l}z_{+f}f1f1

     FA=uTMDPDF_inKT(x1,sqrt(k1),mu,h1)
     FB=uTMDPDF_inKT(x2,sqrt(k2),mu,-h2)!!! -h2, to multiply quarks by anti-quarks in FAB
     FAB=FA*FB

     TMD_pair=XTMD_pairZpZp(FAB,Q2)

  !--------------------------------------------------------------------------------
  CASE (23,24) !Delta^{GG'}z_{-l}z_{-f}{f1f1}_A
     FA=uTMDPDF_inKT(x1,sqrt(k1),mu,h1)
     FB=uTMDPDF_inKT(x2,sqrt(k2),mu,-h2)!!! -h2, to multiply quarks by anti-quarks in FAB
     FAB=FA*FB

    TMD_pair=XTMD_pairZmZm_A(FAB,Q2)

  !--------------------------------------------------------------------------------
  CASE (30,31,32) !Delta^{GG'}z_{+l}r_{+f}h1h1
     FA=BoerMuldersTMDPDF_inKT(x1,sqrt(k1),mu,h1)
     FB=BoerMuldersTMDPDF_inKT(x2,sqrt(k2),mu,-h2)!!! -h2, to myltiply quarks by anti-quarks in FAB
     FAB=FA*FB

    TMD_pair=XTMD_pairZpRp(FAB,Q2)

  !--------------------------------------------------------------------------------
  CASE (35,36) !Delta^{GG'}z_{+l}r_{-f}{h1h1}_A
     FA=BoerMuldersTMDPDF_inKT(x1,sqrt(k1),mu,h1)
     FB=BoerMuldersTMDPDF_inKT(x2,sqrt(k2),mu,-h2)!!! -h2, to myltiply quarks by anti-quarks in FAB
     FAB=FA*FB

    TMD_pair=XTMD_pairZpRm_A(FAB,Q2)

  !--------------------------------------------------------------------------------
  CASE (101) !h1+Cu->gamma* !!this is for E288
    !!!! strictly hadron 1
    FA=uTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    FB=uTMDPDF_inKT(x2,sqrt(k2),mu,1)

    TMD_pair=116d0/567d0*(FA(2)*FB(-2)+FA(-2)*FB(2))+136d0/567d0*(FA(-2)*FB(1)+FA(2)*FB(-1))&
          +34d0/567d0*(FA(-1)*FB(2)+FA(1)*FB(-2))+29d0/567d0*(FA(-1)*FB(1)+FA(1)*FB(-1))&
          +1d0/9d0*(FA(-3)*FB(3)+FA(3)*FB(-3)+4d0*FA(-4)*FB(4)+4d0*FA(4)*FB(-4)+FA(-5)*FB(5)+FA(5)*FB(-5))

  !--------------------------------------------------------------------------------
  CASE (102) !h1+2H->gamma* !!this is for E772
    !!!! strictrly hadron 1
    FA=uTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    FB=uTMDPDF_inKT(x2,sqrt(k2),mu,1)
    TMD_pair=2d0/9d0*(FA(2)*FB(-2)+FA(-2)*FB(2))+2d0/9d0*(FA(-2)*FB(1)+FA(2)*FB(-1))&
          +1d0/18d0*(FA(-1)*FB(2)+FA(1)*FB(-2))+1d0/18d0*(FA(-1)*FB(1)+FA(1)*FB(-1))&
          +1d0/9d0*(FA(-3)*FB(3)+FA(3)*FB(-3)+4d0*FA(-4)*FB(4)+4d0*FA(4)*FB(-4)+FA(-5)*FB(5)+FA(5)*FB(-5))
  !--------------------------------------------------------------------------------
  CASE (103) !h1+W->gamma* !!this is for E537
    !!!! strictrly hadron 1
    !Wolfram has A=183,    Z=74,    N=109
    FA=uTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    FB=uTMDPDF_inKT(x2,sqrt(k2),mu,1)
    TMD_pair=296d0/1647d0*(FA(-2)*FB(2)+FA(2)*FB(-2))+436d0/1647d0*(FA(-2)*FB(1)+FA(2)*FB(-1))&
        +109d0/1647d0*(FA(-1)*FB(2)+FA(1)*FB(-2))+74d0/1647d0*(FA(-1)*FB(1)+FA(1)*FB(-1))&
        +1d0/9d0*(FA(-3)*FB(3)+FA(3)*FB(-3)+4d0*FA(-4)*FB(4)+4d0*FA(4)*FB(-4)+FA(-5)*FB(5)+FA(5)*FB(-5))


  !----------------------------------------------------------------------------------
  !-------------------------SIDIS----------------------------------------------------
  !----------------------------------------------------------------------------------
  CASE (2001,2011,2021,2031) !h1->h2 where !!!! unpolarized SIDIS
    ! e_q^2 *F_q(A)*F_q(B)
    FA=uTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    FB=uTMDFF_inKT(x2,sqrt(k2),mu,h2)
    TMD_pair=FA(1)*FB(1)/9.d0&
      +FA(2)*FB(2)*4.d0/9.d0&
      +FA(3)*FB(3)/9.d0&
      +FA(4)*FB(4)*4d0/9.d0&
      +FA(5)*FB(5)/9d0&
      +FA(-1)*FB(-1)/9.d0&
      +FA(-2)*FB(-2)*4.d0/9.d0&
      +FA(-3)*FB(-3)/9.d0&
      +FA(-4)*FB(-4)*4d0/9.d0&
      +FA(-5)*FB(-5)/9d0

  CASE (2004,2014,2024,2034) !h1->h2 where !!!! unpolarized SIDIS (BM x COLLINS)-part
    ! e_q^2 *F_q(A)*F_q(B)
    FA=uTMDPDF_inKT(x1,sqrt(k1),mu,h1) !!!!!! CHANGE TO BM
    FB=uTMDFF_inKT(x2,sqrt(k2),mu,h2)  !!!!!! CHANGE TO COLLINS
    TMD_pair=FA(1)*FB(1)/9.d0&
      +FA(2)*FB(2)*4.d0/9.d0&
      +FA(3)*FB(3)/9.d0&
      +FA(4)*FB(4)*4d0/9.d0&
      +FA(5)*FB(5)/9d0&
      +FA(-1)*FB(-1)/9.d0&
      +FA(-2)*FB(-2)*4.d0/9.d0&
      +FA(-3)*FB(-3)/9.d0&
      +FA(-4)*FB(-4)*4d0/9.d0&
      +FA(-5)*FB(-5)/9d0

  ! Si tomo estas funciones prueba sí saca numerines, así que algo está pasando con uTMDPDF y uTMDFF ?
  CASE(999)
  !  TMD_pair=(Exp(-0.2d0*k1)+1/(k1+2.))*(Exp(-0.2d0*k2)+1/(k2+2.))
!     TMD_pair=1._dp
  TMD_pair=Exp(-0.1*k1)*Exp(-0.4*k2)


  CASE DEFAULT
    write(*,*) ErrorString('undefined process: ',moduleName),process
    write(*,*) color('Evaluation stop',c_red_bold)
    stop
 END SELECT
 write(44,*) sqrt(k1),FA(1)


end function TMD_pair


!-------------------------------------------------------------------------------------------------------------------------------------
!-------------------------------------------------------------------------------------------------------------------------------------
!!! Combination Delta^{GG'} z_{+l}z_{+f} FF (ONLY FOR Z-BOSON)
function XTMD_pairZZ(FAB,Q2)
     real(dp)::XTMD_pairZZ
     real(dp),intent(in)::Q2
    !!cross-seciton parameters
     real(dp),dimension(-5:5),intent(in):: FAB

     XTMD_pairZZ=zP_ZZ_L*(&
       zP_ZZ_U*FAB(2)&
      +zP_ZZ_D*FAB(1)&
      +zP_ZZ_S*FAB(3)&
      +zP_ZZ_C*FAB(4)&
      +zP_ZZ_B*FAB(5)&
      +zP_ZZ_U*FAB(-2)&
      +zP_ZZ_D*FAB(-1)&
      +zP_ZZ_S*FAB(-3)&
      +zP_ZZ_C*FAB(-4)&
      +zP_ZZ_B*FAB(-5))*&
      Q2*Q2/((Q2-MZ2)**2+GammaZ2*MZ2)

end function XTMD_pairZZ

!!! Combination Delta^{GG'} z_{+l}z_{+f} FF
function XTMD_pairZpZp(FAB,Q2)
     real(dp)::XTMD_pairZpZp
     real(dp),intent(in)::Q2
    !!cross-seciton parameters
     real(dp),dimension(-5:5),intent(in):: FAB

     XTMD_pairZpZp=&
     zP_gg_L*(&!gamma-part
      zP_gg_U*FAB(2)&
      +zP_gg_D*FAB(1)&
      +zP_gg_S*FAB(3)&
      +zP_gg_C*FAB(4)&
      +zP_gg_B*FAB(5)&
      +zP_gg_U*FAB(-2)&
      +zP_gg_D*FAB(-1)&
      +zP_gg_S*FAB(-3)&
      +zP_gg_C*FAB(-4)&
      +zP_gg_B*FAB(-5))&
     +&!gamma-Z interference
     zP_gZ_L*(&
      zP_gZ_U*FAB(2)&
      +zP_gZ_D*FAB(1)&
      +zP_gZ_S*FAB(3)&
      +zP_gZ_C*FAB(4)&
      +zP_gZ_B*FAB(5)&
      +zP_gZ_U*FAB(-2)&
      +zP_gZ_D*FAB(-1)&
      +zP_gZ_S*FAB(-3)&
      +zP_gZ_C*FAB(-4)&
      +zP_gZ_B*FAB(-5))*&
      2d0*Q2*(Q2-MZ2)/((Q2-MZ2)**2+GammaZ2*MZ2)&
     +&!ZZ-contributions
     zP_ZZ_L*(&
       zP_ZZ_U*FAB(2)&
      +zP_ZZ_D*FAB(1)&
      +zP_ZZ_S*FAB(3)&
      +zP_ZZ_C*FAB(4)&
      +zP_ZZ_B*FAB(5)&
      +zP_ZZ_U*FAB(-2)&
      +zP_ZZ_D*FAB(-1)&
      +zP_ZZ_S*FAB(-3)&
      +zP_ZZ_C*FAB(-4)&
      +zP_ZZ_B*FAB(-5))*&
      Q2*Q2/((Q2-MZ2)**2+GammaZ2*MZ2)

end function XTMD_pairZpZp

!!! Combination Delta^{GG'} z_{+l}r_{+f} FF
function XTMD_pairZpRp(FAB,Q2)
     real(dp)::XTMD_pairZpRp
     real(dp),intent(in)::Q2
    !!cross-seciton parameters
     real(dp),dimension(-5:5),intent(in):: FAB


     XTMD_pairZpRp=&
     zP_gg_L*(&!gamma-part   !!!!rP_gg=zP_gg
      zP_gg_U*FAB(2)&
      +zP_gg_D*FAB(1)&
      +zP_gg_S*FAB(3)&
      +zP_gg_C*FAB(4)&
      +zP_gg_B*FAB(5)&
      +zP_gg_U*FAB(-2)&
      +zP_gg_D*FAB(-1)&
      +zP_gg_S*FAB(-3)&
      +zP_gg_C*FAB(-4)&
      +zP_gg_B*FAB(-5))&
     +&!gamma-Z interference  !!!!rP_gZ=zP_gZ
     zP_gZ_L*(&
      zP_gZ_U*FAB(2)&
      +zP_gZ_D*FAB(1)&
      +zP_gZ_S*FAB(3)&
      +zP_gZ_C*FAB(4)&
      +zP_gZ_B*FAB(5)&
      +zP_gZ_U*FAB(-2)&
      +zP_gZ_D*FAB(-1)&
      +zP_gZ_S*FAB(-3)&
      +zP_gZ_C*FAB(-4)&
      +zP_gZ_B*FAB(-5))*&
      2d0*Q2*(Q2-MZ2)/((Q2-MZ2)**2+GammaZ2*MZ2)&
     +&!ZZ-contributions
     zP_ZZ_L*(&
       rP_ZZ_U*FAB(2)&
      +rP_ZZ_D*FAB(1)&
      +rP_ZZ_S*FAB(3)&
      +rP_ZZ_C*FAB(4)&
      +rP_ZZ_B*FAB(5)&
      +rP_ZZ_U*FAB(-2)&
      +rP_ZZ_D*FAB(-1)&
      +rP_ZZ_S*FAB(-3)&
      +rP_ZZ_C*FAB(-4)&
      +rP_ZZ_B*FAB(-5))*&
      Q2*Q2/((Q2-MZ2)**2+GammaZ2*MZ2)

end function XTMD_pairZpRp

!!! Combination Delta^{GG'} z_{-l}z_{-f} {FF}_A
function XTMD_pairZmZm_A(FAB,Q2)
     real(dp)::XTMD_pairZmZm_A
     real(dp),intent(in)::Q2
    !!cross-seciton parameters
     real(dp),dimension(-5:5),intent(in):: FAB

     XTMD_pairZmZm_A=&  !zM_gg =0
     zM_gZ_L*(& !gamma-Z interference
      zM_gZ_U*FAB(2)&
      +zM_gZ_D*FAB(1)&
      +zM_gZ_S*FAB(3)&
      +zM_gZ_C*FAB(4)&
      +zM_gZ_B*FAB(5)&
      -zM_gZ_U*FAB(-2)&
      -zM_gZ_D*FAB(-1)&
      -zM_gZ_S*FAB(-3)&
      -zM_gZ_C*FAB(-4)&
      -zM_gZ_B*FAB(-5))*&
      2d0*Q2*(Q2-MZ2)/((Q2-MZ2)**2+GammaZ2*MZ2)&
     +&!ZZ-contributions
       zM_ZZ_L*(&
       zM_ZZ_U*FAB(2)&
      +zM_ZZ_D*FAB(1)&
      +zM_ZZ_S*FAB(3)&
      +zM_ZZ_C*FAB(4)&
      +zM_ZZ_B*FAB(5)&
      -zM_ZZ_U*FAB(-2)&
      -zM_ZZ_D*FAB(-1)&
      -zM_ZZ_S*FAB(-3)&
      -zM_ZZ_C*FAB(-4)&
      -zM_ZZ_B*FAB(-5))*&
      Q2*Q2/((Q2-MZ2)**2+GammaZ2*MZ2)

end function XTMD_pairZmZm_A

!!! Combination Delta^{GG'} i z_{+l}r_{-f} {FF}_A
function XTMD_pairZpRm_A(FAB,Q2)
     real(dp)::XTMD_pairZpRm_A
     real(dp),intent(in)::Q2
    !!cross-seciton parameters
     real(dp),dimension(-5:5),intent(in):: FAB

     XTMD_pairZpRm_A=&  !zM_gg =0
     zP_gZ_L*(& !gamma-Z interference
      rM_gZ_U*FAB(2)&
      +rM_gZ_D*FAB(1)&
      +rM_gZ_S*FAB(3)&
      +rM_gZ_C*FAB(4)&
      +rM_gZ_B*FAB(5)&
      -rM_gZ_U*FAB(-2)&
      -rM_gZ_D*FAB(-1)&
      -rM_gZ_S*FAB(-3)&
      -rM_gZ_C*FAB(-4)&
      -rM_gZ_B*FAB(-5))*&
      2d0*Q2*sqrt(MZ2*GammaZ2)/((Q2-MZ2)**2+GammaZ2*MZ2)
      !! no ZZ-term


end function XTMD_pairZpRm_A
