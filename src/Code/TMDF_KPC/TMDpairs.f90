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
  CASE (4) !h1+h2-> W+
    FA=uTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    FB=uTMDPDF_inKT(x2,sqrt(k2),mu,h2)

    TMD_pair=paramW_L*(&
    paramW_UD*(FA(2)*FB(-1)+FA(-1)*FB(2))&        !u*dbar+dbar*u
    +paramW_US*(FA(2)*FB(-3)+FA(-3)*FB(2))&        !u*sbar+sbar*u
    +paramW_UB*(FA(2)*FB(-5)+FA(-5)*FB(2))&        !u*bbar+bbar*u
    +paramW_CD*(FA(4)*FB(-1)+FA(-1)*FB(4))&        !c*dbar+dbar*c
    +paramW_CS*(FA(4)*FB(-3)+FA(-3)*FB(4))&        !c*sbar+sbar*c
    +paramW_CB*(FA(4)*FB(-5)+FA(-5)*FB(4))&        !c*bbar+bbar*c
    )*Q2*Q2/((Q2-MW2)**2+GammaW2*MW2)
!--------------------------------------------------------------------------------
  CASE (5) !h1+h2-> W-
    FA=uTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    FB=uTMDPDF_inKT(x2,sqrt(k2),mu,h2)

    TMD_pair=paramW_L*(&
    paramW_UD*(FA(1)*FB(-2)+FA(-2)*FB(1))&        !d*ubar+ubar*d
    +paramW_US*(FA(3)*FB(-2)+FA(-2)*FB(3))&        !s*ubar+ubar*s
    +paramW_UB*(FA(5)*FB(-2)+FA(-2)*FB(5))&        !b*ubar+ubar*b
    +paramW_CD*(FA(1)*FB(-4)+FA(-4)*FB(1))&        !d*cbar+cbar*d
    +paramW_CS*(FA(3)*FB(-4)+FA(-4)*FB(3))&        !s*cbar+cbar*s
    +paramW_CB*(FA(5)*FB(-4)+FA(-4)*FB(5))&        !b*cbar+cbar*b
    )*Q2*Q2/((Q2-MW2)**2+GammaW2*MW2)
!--------------------------------------------------------------------------------
  CASE (6) !h1+h2-> W+ + W-
    FA=uTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    FB=uTMDPDF_inKT(x2,sqrt(k2),mu,h2)

    TMD_pair=paramW_L*(&
    paramW_UD*(FA(2)*FB(-1)+FA(1)*FB(-2)+FA(-2)*FB(1)+FA(-1)*FB(2))&    !u*dbar+d*ubar+ubar*d+dbar*u
    +paramW_US*(FA(2)*FB(-3)+FA(3)*FB(-2)+FA(-2)*FB(3)+FA(-3)*FB(2))&    !u*sbar+s*ubar+ubar*s+sbar*u
    +paramW_UB*(FA(2)*FB(-5)+FA(5)*FB(-2)+FA(-2)*FB(5)+FA(-5)*FB(2))&    !u*bbar+b*ubar+ubar*b+bbar*u
    +paramW_CD*(FA(4)*FB(-1)+FA(1)*FB(-4)+FA(-4)*FB(1)+FA(-1)*FB(4))&    !c*dbar+d*cbar+cbar*d+dbar*c
    +paramW_CS*(FA(4)*FB(-3)+FA(3)*FB(-4)+FA(-4)*FB(3)+FA(-3)*FB(4))&    !c*sbar+s*cbar+cbar*s+sbar*c
    +paramW_CB*(FA(4)*FB(-5)+FA(5)*FB(-4)+FA(-4)*FB(5)+FA(-5)*FB(4))&    !c*bbar+b*cbar+cbar*b+bbar*c
    )*Q2*Q2/((Q2-MW2)**2+GammaW2*MW2)

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
  CASE (2001,2011,2021,2031,2041) !h1->h2 where !!!! unpolarized SIDIS
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

  CASE (2002,2012,2022,2032,2042) !d->h2 where d is deutron prepared from hadron 1 [i.e u->(u+d)/2, d->(u+d)/2]
    ! e_q^2 *F_q(A)*F_q(B)
    FA=uTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    FB=uTMDFF_inKT(x2,sqrt(k2),mu,h2)
    TMD_pair=(FA(1)+FA(2))*(FB(1)+4d0*FB(2))/18d0&
      +FA(3)*FB(3)/9.d0&
      +FA(4)*FB(4)*4d0/9.d0&
      +FA(5)*FB(5)/9d0&
      +(FA(-1)+FA(-2))*(FB(-1)+4d0*FB(-2))/18d0&
      +FA(-3)*FB(-3)/9.d0&
      +FA(-4)*FB(-4)*4d0/9.d0&
      +FA(-5)*FB(-5)/9d0

  CASE (2003,2013,2023,2033,2043) !n->h2 where n=last number (n=neutron=p(u<->d))
    ! e_q^2 *F_q(A)*F_q(B)
    FA=uTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    FB=uTMDFF_inKT(x2,sqrt(k2),mu,h2)
    TMD_pair=FA(2)*FB(1)/9.d0&
      +FA(1)*FB(2)*4.d0/9.d0&
      +FA(3)*FB(3)/9.d0&
      +FA(4)*FB(4)*4d0/9.d0&
      +FA(5)*FB(5)/9d0&
      +FA(-2)*FB(-1)/9.d0&
      +FA(-1)*FB(-2)*4.d0/9.d0&
      +FA(-3)*FB(-3)/9.d0&
      +FA(-4)*FB(-4)*4d0/9.d0&
      +FA(-5)*FB(-5)/9d0

  !--------------------------------------------------------------------------------
   CASE (2101,2111,2121,2131,2141) !p->h? where h?=h1+h2
    ! e_q^2 *F_q(A)*F_q(B)
    FA=uTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    if(h2>0) then
        FB=uTMDFF_inKT(x2,sqrt(k2),mu,1)+uTMDFF_inKT(x2,sqrt(k2),mu,2)
    else
        FB=uTMDFF_inKT(x2,sqrt(k2),mu,-1)+uTMDFF_inKT(x2,sqrt(k2),mu,-2)
    end if
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
!--------------------------------------------------------------------------------
    CASE (2102,2112,2122,2132,2142) !p->h? where h?=h1+h2+h3
    ! e_q^2 *F_q(A)*F_q(B)
    FA=uTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    if(h2>0) then
        FB=uTMDFF_inKT(x2,sqrt(k2),mu,1)+uTMDFF_inKT(x2,sqrt(k2),mu,2)+uTMDFF_inKT(x2,sqrt(k2),mu,3)
    else
        FB=uTMDFF_inKT(x2,sqrt(k2),mu,-1)+uTMDFF_inKT(x2,sqrt(k2),mu,-2)+uTMDFF_inKT(x2,sqrt(k2),mu,-3)
    end if
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
!--------------------------------------------------------------------------------
    CASE (2103,2113,2123,2133,2143) !d->h? where h?=h1+h2 (d=deutron=(p+n)/2)
    ! e_q^2 *F_q(A)*F_q(B)
    FA=uTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    if(h2>0) then
        FB=uTMDFF_inKT(x2,sqrt(k2),mu,1)+uTMDFF_inKT(x2,sqrt(k2),mu,2)
    else
        FB=uTMDFF_inKT(x2,sqrt(k2),mu,-1)+uTMDFF_inKT(x2,sqrt(k2),mu,-2)
    end if
    TMD_pair=(FA(1)+FA(2))*(FB(1)+4d0*FB(2))/18d0&
      +FA(3)*FB(3)/9.d0&
      +FA(4)*FB(4)*4d0/9.d0&
      +FA(5)*FB(5)/9d0&
      +(FA(-1)+FA(-2))*(FB(-1)+4d0*FB(-2))/18d0&
      +FA(-3)*FB(-3)/9.d0&
      +FA(-4)*FB(-4)*4d0/9.d0&
      +FA(-5)*FB(-5)/9d0
!--------------------------------------------------------------------------------
    CASE (2104,2114,2124,2134,2144) !d->h? where h?=h1+h2+h3 (d=deutron=(p+n)/2)
    ! e_q^2 *F_q(A)*F_q(B)
    FA=uTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    if(h2>0) then
        FB=uTMDFF_inKT(x2,sqrt(k2),mu,1)+uTMDFF_inKT(x2,sqrt(k2),mu,2)+uTMDFF_inKT(x2,sqrt(k2),mu,3)
    else
        FB=uTMDFF_inKT(x2,sqrt(k2),mu,-1)+uTMDFF_inKT(x2,sqrt(k2),mu,-2)+uTMDFF_inKT(x2,sqrt(k2),mu,-3)
    end if
    TMD_pair=(FA(1)+FA(2))*(FB(1)+4d0*FB(2))/18d0&
      +FA(3)*FB(3)/9.d0&
      +FA(4)*FB(4)*4d0/9.d0&
      +FA(5)*FB(5)/9d0&
      +(FA(-1)+FA(-2))*(FB(-1)+4d0*FB(-2))/18d0&
      +FA(-3)*FB(-3)/9.d0&
      +FA(-4)*FB(-4)*4d0/9.d0&
      +FA(-5)*FB(-5)/9d0
!--------------------------------------------------------------------------------
   CASE (2105,2115,2125,2135,2145) !n->h? where h?=h1+h2 (n=neutron=p(u<->d))
    ! e_q^2 *F_q(A)*F_q(B)
    FA=uTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    if(h2>0) then
        FB=uTMDFF_inKT(x2,sqrt(k2),mu,1)+uTMDFF_inKT(x2,sqrt(k2),mu,2)
    else
        FB=uTMDFF_inKT(x2,sqrt(k2),mu,-1)+uTMDFF_inKT(x2,sqrt(k2),mu,-2)
    end if
    TMD_pair=FA(2)*FB(1)/9.d0&
      +FA(1)*FB(2)*4.d0/9.d0&
      +FA(3)*FB(3)/9.d0&
      +FA(4)*FB(4)*4d0/9.d0&
      +FA(5)*FB(5)/9d0&
      +FA(-2)*FB(-1)/9.d0&
      +FA(-1)*FB(-2)*4.d0/9.d0&
      +FA(-3)*FB(-3)/9.d0&
      +FA(-4)*FB(-4)*4d0/9.d0&
      +FA(-5)*FB(-5)/9d0
!--------------------------------------------------------------------------------
    CASE (2106,2116,2126,2136,2146) !n->h? where h?=h1+h2+h3 (n=neutron=p(u<->d))
    ! e_q^2 *F_q(A)*F_q(B)
    FA=uTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    if(h2>0) then
        FB=uTMDFF_inKT(x2,sqrt(k2),mu,1)+uTMDFF_inKT(x2,sqrt(k2),mu,2)+uTMDFF_inKT(x2,sqrt(k2),mu,3)
    else
        FB=uTMDFF_inKT(x2,sqrt(k2),mu,-1)+uTMDFF_inKT(x2,sqrt(k2),mu,-2)+uTMDFF_inKT(x2,sqrt(k2),mu,-3)
    end if
    TMD_pair=FA(2)*FB(1)/9.d0&
      +FA(1)*FB(2)*4.d0/9.d0&
      +FA(3)*FB(3)/9.d0&
      +FA(4)*FB(4)*4d0/9.d0&
      +FA(5)*FB(5)/9d0&
      +FA(-2)*FB(-1)/9.d0&
      +FA(-1)*FB(-2)*4.d0/9.d0&
      +FA(-3)*FB(-3)/9.d0&
      +FA(-4)*FB(-4)*4d0/9.d0&
      +FA(-5)*FB(-5)/9d0
!--------------------------------------------------------------------------------
   CASE (2107,2117,2127,2137,2147) !p->h? where h?=h1+h2 [from 3+4]
    ! e_q^2 *F_q(A)*F_q(B)
    FA=uTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    if(h2>0) then
        FB=uTMDFF_inKT(x2,sqrt(k2),mu,3)+uTMDFF_inKT(x2,sqrt(k2),mu,4)
    else
        FB=uTMDFF_inKT(x2,sqrt(k2),mu,-3)+uTMDFF_inKT(x2,sqrt(k2),mu,-4)
    end if
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
!--------------------------------------------------------------------------------
    CASE (2108,2118,2128,2138,2148) !d->h? where h?=h1+h2 (d=deutron=(p+n)/2) [from 3+4]
    ! e_q^2 *F_q(A)*F_q(B)
    FA=uTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    if(h2>0) then
        FB=uTMDFF_inKT(x2,sqrt(k2),mu,3)+uTMDFF_inKT(x2,sqrt(k2),mu,4)
    else
        FB=uTMDFF_inKT(x2,sqrt(k2),mu,-3)+uTMDFF_inKT(x2,sqrt(k2),mu,-4)
    end if
    TMD_pair=(FA(1)+FA(2))*(FB(1)+4d0*FB(2))/18d0&
      +FA(3)*FB(3)/9.d0&
      +FA(4)*FB(4)*4d0/9.d0&
      +FA(5)*FB(5)/9d0&
      +(FA(-1)+FA(-2))*(FB(-1)+4d0*FB(-2))/18d0&
      +FA(-3)*FB(-3)/9.d0&
      +FA(-4)*FB(-4)*4d0/9.d0&
      +FA(-5)*FB(-5)/9d0
!--------------------------------------------------------------------------------
   CASE (2109,2119,2129,2139,2149) !n->h? where h?=h1+h2 (n=neutron=p(u<->d))[from 3+4]
    ! e_q^2 *F_q(A)*F_q(B)
    FA=uTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    if(h2>0) then
        FB=uTMDFF_inKT(x2,sqrt(k2),mu,3)+uTMDFF_inKT(x2,sqrt(k2),mu,4)
    else
        FB=uTMDFF_inKT(x2,sqrt(k2),mu,-3)+uTMDFF_inKT(x2,sqrt(k2),mu,-4)
    end if
    TMD_pair=FA(2)*FB(1)/9.d0&
      +FA(1)*FB(2)*4.d0/9.d0&
      +FA(3)*FB(3)/9.d0&
      +FA(4)*FB(4)*4d0/9.d0&
      +FA(5)*FB(5)/9d0&
      +FA(-2)*FB(-1)/9.d0&
      +FA(-1)*FB(-2)*4.d0/9.d0&
      +FA(-3)*FB(-3)/9.d0&
      +FA(-4)*FB(-4)*4d0/9.d0&
      +FA(-5)*FB(-5)/9d0

!--------------------------------------------------------------------------------
   CASE (3001,3011,3021,3031,3041)  !h1->h2 where !!!! unpolarized SIDIS (BM x COLLINS)-part
    ! e_q^2 *F_q(A)*F_q(B)
    FA=BoerMuldersTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    FB=CollinsTMDFF_inKT(x2,sqrt(k2),mu,h2)
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

  CASE (3002,3012,3022,3032,3042) !h1->h2 where !!!! unpolarized SIDIS (BM x COLLINS)-part
    !d->h2 where d is deutron prepared from hadron 1 [i.e u->(u+d)/2, d->(u+d)/2]
    ! e_q^2 *F_q(A)*F_q(B)
    FA=BoerMuldersTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    FB=CollinsTMDFF_inKT(x2,sqrt(k2),mu,h2)
    TMD_pair=(FA(1)+FA(2))*(FB(1)+4d0*FB(2))/18d0&
      +FA(3)*FB(3)/9.d0&
      +FA(4)*FB(4)*4d0/9.d0&
      +FA(5)*FB(5)/9d0&
      +(FA(-1)+FA(-2))*(FB(-1)+4d0*FB(-2))/18d0&
      +FA(-3)*FB(-3)/9.d0&
      +FA(-4)*FB(-4)*4d0/9.d0&
      +FA(-5)*FB(-5)/9d0

  CASE (3003,3013,3023,3033,3043) !!!! unpolarized SIDIS (BM x COLLINS)-part
    !n->h2 where n=last number (n=neutron=p(u<->d))
    ! e_q^2 *F_q(A)*F_q(B)
    ! e_q^2 *F_q(A)*F_q(B)
    FA=BoerMuldersTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    FB=CollinsTMDFF_inKT(x2,sqrt(k2),mu,h2)
    TMD_pair=FA(2)*FB(1)/9.d0&
      +FA(1)*FB(2)*4.d0/9.d0&
      +FA(3)*FB(3)/9.d0&
      +FA(4)*FB(4)*4d0/9.d0&
      +FA(5)*FB(5)/9d0&
      +FA(-2)*FB(-1)/9.d0&
      +FA(-1)*FB(-2)*4.d0/9.d0&
      +FA(-3)*FB(-3)/9.d0&
      +FA(-4)*FB(-4)*4d0/9.d0&
      +FA(-5)*FB(-5)/9d0
  !--------------------------------------------------------------------------------
   CASE (3101,3111,3121,3131,3141) !p->h? where h?=h1+h2 !!!! unpolarized SIDIS (BM x COLLINS)-part
    ! e_q^2 *F_q(A)*F_q(B)
    FA=BoerMuldersTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    if(h2>0) then
        FB=CollinsTMDFF_inKT(x2,sqrt(k2),mu,1)+CollinsTMDFF_inKT(x2,sqrt(k2),mu,2)
    else
        FB=CollinsTMDFF_inKT(x2,sqrt(k2),mu,-1)+CollinsTMDFF_inKT(x2,sqrt(k2),mu,-2)
    end if
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
!--------------------------------------------------------------------------------
    CASE (3102,3112,3122,3132,3142) !p->h? where h?=h1+h2+h3 !!!! unpolarized SIDIS (BM x COLLINS)-part
    ! e_q^2 *F_q(A)*F_q(B)
    FA=BoerMuldersTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    if(h2>0) then
        FB=CollinsTMDFF_inKT(x2,sqrt(k2),mu,1)+CollinsTMDFF_inKT(x2,sqrt(k2),mu,2)+CollinsTMDFF_inKT(x2,sqrt(k2),mu,3)
    else
        FB=CollinsTMDFF_inKT(x2,sqrt(k2),mu,-1)+CollinsTMDFF_inKT(x2,sqrt(k2),mu,-2)+CollinsTMDFF_inKT(x2,sqrt(k2),mu,-3)
    end if
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
!--------------------------------------------------------------------------------
    CASE (3103,3113,3123,3133,3143) !d->h? where h?=h1+h2 (d=deutron=(p+n)/2) !!!! unpolarized SIDIS (BM x COLLINS)-part
    ! e_q^2 *F_q(A)*F_q(B)
    FA=BoerMuldersTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    if(h2>0) then
        FB=CollinsTMDFF_inKT(x2,sqrt(k2),mu,1)+CollinsTMDFF_inKT(x2,sqrt(k2),mu,2)
    else
        FB=CollinsTMDFF_inKT(x2,sqrt(k2),mu,-1)+CollinsTMDFF_inKT(x2,sqrt(k2),mu,-2)
    end if
    TMD_pair=(FA(1)+FA(2))*(FB(1)+4d0*FB(2))/18d0&
      +FA(3)*FB(3)/9.d0&
      +FA(4)*FB(4)*4d0/9.d0&
      +FA(5)*FB(5)/9d0&
      +(FA(-1)+FA(-2))*(FB(-1)+4d0*FB(-2))/18d0&
      +FA(-3)*FB(-3)/9.d0&
      +FA(-4)*FB(-4)*4d0/9.d0&
      +FA(-5)*FB(-5)/9d0
!--------------------------------------------------------------------------------
    CASE (3104,3114,3124,3134,3144) !d->h? where h?=h1+h2+h3 (d=deutron=(p+n)/2) !!!! unpolarized SIDIS (BM x COLLINS)-part
    ! e_q^2 *F_q(A)*F_q(B)
    FA=BoerMuldersTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    if(h2>0) then
        FB=CollinsTMDFF_inKT(x2,sqrt(k2),mu,1)+CollinsTMDFF_inKT(x2,sqrt(k2),mu,2)+CollinsTMDFF_inKT(x2,sqrt(k2),mu,3)
    else
        FB=CollinsTMDFF_inKT(x2,sqrt(k2),mu,-1)+CollinsTMDFF_inKT(x2,sqrt(k2),mu,-2)+CollinsTMDFF_inKT(x2,sqrt(k2),mu,-3)
    end if
    TMD_pair=(FA(1)+FA(2))*(FB(1)+4d0*FB(2))/18d0&
      +FA(3)*FB(3)/9.d0&
      +FA(4)*FB(4)*4d0/9.d0&
      +FA(5)*FB(5)/9d0&
      +(FA(-1)+FA(-2))*(FB(-1)+4d0*FB(-2))/18d0&
      +FA(-3)*FB(-3)/9.d0&
      +FA(-4)*FB(-4)*4d0/9.d0&
      +FA(-5)*FB(-5)/9d0
!--------------------------------------------------------------------------------
   CASE (3105,3115,3125,3135,3145) !n->h? where h?=h1+h2 (n=neutron=p(u<->d)) !!!! unpolarized SIDIS (BM x COLLINS)-part
    ! e_q^2 *F_q(A)*F_q(B)
    FA=BoerMuldersTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    if(h2>0) then
        FB=CollinsTMDFF_inKT(x2,sqrt(k2),mu,1)+CollinsTMDFF_inKT(x2,sqrt(k2),mu,2)
    else
        FB=CollinsTMDFF_inKT(x2,sqrt(k2),mu,-1)+CollinsTMDFF_inKT(x2,sqrt(k2),mu,-2)
    end if
    TMD_pair=FA(2)*FB(1)/9.d0&
      +FA(1)*FB(2)*4.d0/9.d0&
      +FA(3)*FB(3)/9.d0&
      +FA(4)*FB(4)*4d0/9.d0&
      +FA(5)*FB(5)/9d0&
      +FA(-2)*FB(-1)/9.d0&
      +FA(-1)*FB(-2)*4.d0/9.d0&
      +FA(-3)*FB(-3)/9.d0&
      +FA(-4)*FB(-4)*4d0/9.d0&
      +FA(-5)*FB(-5)/9d0
!--------------------------------------------------------------------------------
    CASE (3106,3116,3126,3136,3146) !n->h? where h?=h1+h2+h3 (n=neutron=p(u<->d)) !!!! unpolarized SIDIS (BM x COLLINS)-part
    ! e_q^2 *F_q(A)*F_q(B)
    FA=BoerMuldersTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    if(h2>0) then
        FB=CollinsTMDFF_inKT(x2,sqrt(k2),mu,1)+CollinsTMDFF_inKT(x2,sqrt(k2),mu,2)+CollinsTMDFF_inKT(x2,sqrt(k2),mu,3)
    else
        FB=CollinsTMDFF_inKT(x2,sqrt(k2),mu,-1)+CollinsTMDFF_inKT(x2,sqrt(k2),mu,-2)+CollinsTMDFF_inKT(x2,sqrt(k2),mu,-3)
    end if
    TMD_pair=FA(2)*FB(1)/9.d0&
      +FA(1)*FB(2)*4.d0/9.d0&
      +FA(3)*FB(3)/9.d0&
      +FA(4)*FB(4)*4d0/9.d0&
      +FA(5)*FB(5)/9d0&
      +FA(-2)*FB(-1)/9.d0&
      +FA(-1)*FB(-2)*4.d0/9.d0&
      +FA(-3)*FB(-3)/9.d0&
      +FA(-4)*FB(-4)*4d0/9.d0&
      +FA(-5)*FB(-5)/9d0
!--------------------------------------------------------------------------------
   CASE (3107,3117,3127,3137,3147) !p->h? where h?=h1+h2 [from 3+4] !!!! unpolarized SIDIS (BM x COLLINS)-part
    ! e_q^2 *F_q(A)*F_q(B)
    FA=BoerMuldersTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    if(h2>0) then
        FB=CollinsTMDFF_inKT(x2,sqrt(k2),mu,3)+CollinsTMDFF_inKT(x2,sqrt(k2),mu,4)
    else
        FB=CollinsTMDFF_inKT(x2,sqrt(k2),mu,-3)+CollinsTMDFF_inKT(x2,sqrt(k2),mu,-4)
    end if
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
!--------------------------------------------------------------------------------
    CASE (3108,3118,3128,3138,3148) !d->h? where h?=h1+h2 (d=deutron=(p+n)/2) [from 3+4] !!!! unpolarized SIDIS (BM x COLLINS)-part
    ! e_q^2 *F_q(A)*F_q(B)
    FA=BoerMuldersTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    if(h2>0) then
        FB=CollinsTMDFF_inKT(x2,sqrt(k2),mu,3)+CollinsTMDFF_inKT(x2,sqrt(k2),mu,4)
    else
        FB=CollinsTMDFF_inKT(x2,sqrt(k2),mu,-3)+CollinsTMDFF_inKT(x2,sqrt(k2),mu,-4)
    end if
    TMD_pair=(FA(1)+FA(2))*(FB(1)+4d0*FB(2))/18d0&
      +FA(3)*FB(3)/9.d0&
      +FA(4)*FB(4)*4d0/9.d0&
      +FA(5)*FB(5)/9d0&
      +(FA(-1)+FA(-2))*(FB(-1)+4d0*FB(-2))/18d0&
      +FA(-3)*FB(-3)/9.d0&
      +FA(-4)*FB(-4)*4d0/9.d0&
      +FA(-5)*FB(-5)/9d0
!--------------------------------------------------------------------------------
   CASE (3109,3119,3129,3139,3149) !n->h? where h?=h1+h2 (n=neutron=p(u<->d))[from 3+4] !!!! unpolarized SIDIS (BM x COLLINS)-part
    ! e_q^2 *F_q(A)*F_q(B)
    FA=BoerMuldersTMDPDF_inKT(x1,sqrt(k1),mu,h1)
    if(h2>0) then
        FB=CollinsTMDFF_inKT(x2,sqrt(k2),mu,3)+CollinsTMDFF_inKT(x2,sqrt(k2),mu,4)
    else
        FB=CollinsTMDFF_inKT(x2,sqrt(k2),mu,-3)+CollinsTMDFF_inKT(x2,sqrt(k2),mu,-4)
    end if
    TMD_pair=FA(2)*FB(1)/9.d0&
      +FA(1)*FB(2)*4.d0/9.d0&
      +FA(3)*FB(3)/9.d0&
      +FA(4)*FB(4)*4d0/9.d0&
      +FA(5)*FB(5)/9d0&
      +FA(-2)*FB(-1)/9.d0&
      +FA(-1)*FB(-2)*4.d0/9.d0&
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
!  write(44,*) sqrt(k1),FA(1)


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
