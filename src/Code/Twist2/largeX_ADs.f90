!!!-------------------------------------------------------------------------------------------------------
!!!--------------The coefficeints of the delta-function of twist-2 matching ------------------------------
!!!--------------Used for large-X resummation ------------------------------------------------------------
!!!-------------------------------------------------------------------------------------------------------

    
!!!!!coefficient function q<-q delta-part
!!!! The individual coefficients are written separately in order to possibility to acces them separately.
pure function C_q_q_delta_1(Nf,Lmu)
  integer,intent(in)::Nf
  real(dp),intent(in)::Lmu
  real(dp)::C_q_q_delta_1
  C_q_q_delta_1=(-4d0/3d0*zeta2-4d0*Lmu)
end function C_q_q_delta_1

pure function C_q_q_delta_2(Nf,Lmu)
  real(dp),intent(in)::Lmu
  integer,intent(in)::Nf
  real(dp)::C_q_q_delta_2
  C_q_q_delta_2=(&
     -2416d0/81d0-134d0/3d0*zeta2+448d0/9d0*zeta3+200d0/9d0*zeta4+&
     Nf*(352d0/243d0+20d0/9d0*zeta2+56d0/27d0*zeta3)+&
     Lmu*(-14d0-140d0/3d0*zeta2+16d0/3d0*zeta3+Nf*(4d0/9d0+40d0/9d0*zeta2))+&
     Lmu**2*(-14d0+4d0/3d0*Nf-128d0/9d0*zeta2))
end function C_q_q_delta_2

pure function C_q_q_delta_3(Nf,Lmu)
  real(dp),intent(in)::Lmu
  integer,intent(in)::Nf
  real(dp)::C_q_q_delta_3
  C_q_q_delta_3=(&
    Nf**2*(2800d0/19683d0 - (496d0*zeta2)/81d0 - (3712d0*zeta3)/729d0 - (88d0*zeta4)/81d0)&
    +Nf*(212644d0/6561d0 + (224116d0*zeta2)/729d0 + (83452d0*zeta3)/729d0 - (5960d0*zeta2*zeta3)/81d0 &
      - (1988d0*zeta4)/81d0 - (8144d0*zeta5)/81d0)&
    -578966d0/2187d0 - (667234d0*zeta2)/243d0 + (13954d0*zeta3)/81d0 + (30212d0*zeta2*zeta3)/27d0 &
    -(244d0*zeta3**2)/3d0 + (12796d0*zeta4)/9d0  - (1576d0*zeta5)/3d0 - (59468d0*zeta6)/81d0&
    +Lmu*(&
        Nf**2*(428d0/729d0 - (560d0*zeta2)/81d0 - (32d0*zeta3)/81d0)&
        + Nf*(6496d0/243d0 + (80536d0*zeta2)/243d0 - (112d0*zeta3)/3d0 - (1064d0*zeta4)/9d0)&
        -5947d0/27d0 - (232612d0*zeta2)/81d0 + (35200d0*zeta3)/27d0 + (15680d0*zeta2*zeta3)/27d0&
        + (3364d0*zeta4)/3d0 - (4960d0*zeta5)/9d0)&
    +Lmu**2*(-302d0 + Nf**2*(-8d0/27d0 - (80d0*zeta2)/27d0) - (26476d0*zeta2)/27d0 &
        + Nf*(340d0/9d0 + (7744d0*zeta2)/81d0 - (32d0*zeta3)/9d0) + (112d0*zeta3)/3d0 + (12800d0*zeta4)/27d0)&
    +Lmu**3*(-84d0 - 16d0*Nf**2/27d0 - 896d0*zeta2/9d0 + Nf*(128d0/9d0 + 256d0*zeta2/27d0) - 4096d0*zeta3/81d0)&
    )
end function C_q_q_delta_3
  
!!!!!coefficient function g<-g delta-part
!!!! The individual coefficients are written separately in order to possibility to acces them separately.
pure function C_g_g_delta_1(Nf,Lmu)
  real(dp),intent(in)::Lmu
  integer,intent(in)::Nf
  real(dp)::C_g_g_delta_1
  C_g_g_delta_1=(-3d0*zeta2-(11d0-2d0/3d0*Nf)*Lmu)
end function C_g_g_delta_1

pure function C_g_g_delta_2(Nf,Lmu)
  real(dp),intent(in)::Lmu
  integer,intent(in)::Nf
  real(dp)::C_g_g_delta_2
  C_g_g_delta_2=(&
    -112d0 - 56d0*Nf**2/81d0 - 201d0*zeta2/2d0 - 72d0*Lmu**2*zeta2 +&
    Lmu*(-96d0 + 32d0*Nf/3d0 - 108d0*zeta3) + Nf*(548d0/27 + 5d0*zeta2 - 28d0*zeta3/3d0) + 154d0*zeta3 + 225d0*zeta4/4d0)
end function C_g_g_delta_2

pure function C_g_g_delta_3(Nf,Lmu)
  real(dp),intent(in)::Lmu
  integer,intent(in)::Nf
  real(dp)::C_g_g_delta_3
  C_g_g_delta_3=(&
    Nf**3*(-752d0/2187d0 + (16d0*zeta3)/27d0)&
    + Nf**2*(-73577d0/2187d0 - (200d0*zeta2)/81d0 - (368d0*zeta3)/81d0 - (20d0*zeta4)/9d0)&
    + Nf*(1033259d0/1458d0 + (27305d0*zeta2)/81d0 - (17762d0*zeta3)/81d0 + (122d0*zeta2*zeta3)/3d0 &
    - (2767d0*zeta4)/18d0 + (556d0*zeta5)/9d0)&
    -698456d0/243d0 - (213865d0*zeta2)/54d0 + (1489d0*zeta3)/9d0 + 429d0*zeta2*zeta3 + 2337d0*zeta3**2&
    + (15395d0*zeta4)/4d0 - 2046d0*zeta5  - (28783d0*zeta6)/16d0&
    +Lmu*((112d0*Nf**3)/243d0&
      +Nf**2*(-4471d0/162d0 - (10d0*zeta2)/3d0 + (56d0*zeta3)/9d0)&
      +Nf*(25175d0/54d0 + (904d0*zeta2)/3d0 + (104d0*zeta3)/3d0 - (15d0*zeta4)/2d0)&
      -4597d0/2d0 - (8855d0*zeta2)/2d0 - 3130d0*zeta3 + 3780d0*zeta2*zeta3 + (495d0*zeta4)/4d0 + 2160d0*zeta5)&
    + Lmu**2*(-561d0 - (38d0*Nf**2)/9d0 - 3216d0*zeta2 + Nf*(311d0/3d0 + 160d0*zeta2) + 2700d0*zeta4)&
    - Lmu**3*576d0*zeta3&
    )
end function C_g_g_delta_3
