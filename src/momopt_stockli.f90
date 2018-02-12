!
! This routine creates the aerodynamcal parameters table using MOMOPT code (Stockli's version)
!
! To compile: 
!
! ifort -r8 -traceback aero.f90 momopt_stockli.f90 -o momopt_stockli.out
! gfortran -freal-4-real-8 -fbacktrace aero.f90 momopt_stockli.f90 -o momopt_stockli.out
!
! To run:
!
! ./momopt_stockli.out
!
subroutine momopt_stockli(z1, z2, zc, chil, zlw, zlen, vcov1, zlt, &
      zs, zmet, zwind, ha, z0, d, g2, g3, rbc, rdc, corb1, corb2)

   use const, only : pi, vkc

   implicit none

   !Dummy arguments
   real, intent(in) :: z1, z2, zc, chil, zlw, zlen, vcov1, zlt, zs, zmet, zwind
   real :: ha, z0, d, g2, g3, rbc, rdc, corb1, corb2

   !Local variables
   integer, parameter :: igcm = 0

   real :: flowb, reyno, cdd, dl, ps, cs, cdg
   real :: zl, u2fac
   real :: z0z2, dz2, haz2, ustaru
   real :: sigma     ! (output parameter in this MOMOPT version)

   !Variables to plot curves
   real :: COEFS(4) ! Storage vector for YA and YB

!      zs     , &     ! Ground roughness length (m)

!      z1  , &        ! Canopy-base height (m)
!      z2  , &        ! Canopy-top height (m)
!      zc  , &        ! Inflection height for leaf-area density (m)
!      chil  , &      ! Leaf area distribution factor
!      leafw  , &     ! Leaf width (m)
!      leafl  , &     ! Leaf length (m)
!      vcover  , &    ! Canopy cover fraction

   ! The lines below were copied from SIBX.F code
   FLOWB = (1.0 - CHIL) / PI
   FLOWB = MAX(FLOWB, 0.0)
   REYNO = (ZLW + ZLEN) /2.0 * 10000.0 / 0.15
   CDD = 1.328 *2.0 / SQRT(REYNO) + 0.45 * FLOWB ** 1.6
   CS = 90.0 * SQRT(ZLW)
   CDG = (VKC / LOG(Z1 / ZS)) ** 2.0

   ! The lines below were copied from SIBX.F code
   DL = ZLT / (Z2 - Z1) / VCOV1
   PS = 1.0 + DL ** 0.6

   CALL MOMOPT (Z2, ZC, Z1, ZS, ZWIND, ZMET, ZLT, &
      CDD, PS, CS, &
      Z0, D, RBC, RDC, &
      HA, G2, G3, CORB1, CORB2, U2FAC, SIGMA, ZL, COEFS)

   ! The lines below were copied from SIBX.F code
   IF (IGCM .EQ. 0) U2FAC = 1.0
   RBC = RBC / SQRT(U2FAC)
   RDC = RDC / U2FAC
   Z0Z2 = Z0 / Z2
   DZ2 = D / Z2
   HAZ2 = HA / Z2
   USTARU = 0.41 / LOG ((ZWIND - D) / Z0)

end subroutine momopt_stockli
