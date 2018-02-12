!!!! Reto Stockli (June 2008): Recoded to F90
!!!
!!!subroutine aero (LAI,fVCover,ChiL,z2,z1,zc,LWidth,LLength,z0,zp_disp,RbC,RdC)
!!!
!!!  ! This routine is the outermost shell for calculating the aerodynamic
!!!  ! variables for input into the Simple Biosphere model (SiB).  It converts
!!!  ! single to double precision, assigns some reference values, and calls
!!!  ! Sibx (which actually calculates the aerodynamic parameters).
!!!  !-----------------------------------------------------------------------
!!!
!!!  implicit none
!!!
!!!  ! begin single precision variables
!!!  
!!!  ! begin independant input variables
!!!  real(kind=4) LAI      ! Leaf-area index
!!!  real(kind=4) fVCover  ! Fractional vegetation cover
!!!  
!!!  ! begin input morphilogical variables
!!!  real(kind=4) z2     ! Canopy top height (m)
!!!  real(kind=4) z1     ! Canopy base height (m)
!!!  real(kind=4) zc       ! Inflection point for leaf density
!!!  real(kind=4) LWidth   ! Leaf width
!!!  real(kind=4) LLength  ! Leaf length
!!!  real(kind=4) ChiL     ! Leaf angle distribution factor, or Xl
!!!
!!!  ! begin output aerodynamic variables
!!!  real(kind=4) RbC      ! RB Coefficient (c1 or cc1)
!!!  real(kind=4) RdC      ! RC Coefficient (c2 or cc2)
!!!  real(kind=4) z0       ! Canopy roughness length
!!!  real(kind=4) zp_disp  ! Zero plane displacement
!!!  real(kind=4) G2       ! Ratio Ra (actual) to Ra (log-linear) for momentum
!!!  real(kind=4) G3       ! Ratio Ra (actual) to Ra (log-linear) for heat transfer
!!!  real(kind=4) CORB1    ! Non-neutral correction for Ra between Ha and z2
!!!  real(kind=4) CORB2    ! Neutral value RdC^2 between Ha and Z2.  Scaled by U2.
!!!  real(kind=4) HA       ! Canopy source height
!!!  real(kind=4) sigma    ! Canopy mixing length (m)
!!!
!!!  ! begin internal constants
!!!  real(kind=4) zwind    ! Reference height for momentum (wind), meters
!!!  
!!!  real(kind=4) zmet     ! Reference height for weather (meteorology), meters
!!!  
!!!  real(kind=4) zs       ! Soil roughness length, meters
!!!  
!!!  ! repeat variable definition with double precision
!!!
!!!  ! begin independant variables
!!!  real(kind=8) dLAI   ! Leaf-area index
!!!  real(kind=8) dfVCover ! Fractional vegetation cover
!!!  
!!!  ! begin morphological variables
!!!  real(kind=8) dz2    ! Canopy top height (m)
!!!  real(kind=8) dz1    ! Canopy base height (m)
!!!  real(kind=8) dzc    ! Inflection point for leaf density
!!!  real(kind=8) dLWidth    ! Leaf width
!!!  real(kind=8) dLLength    ! Leaf length
!!!  real(kind=8) dChiL  ! Leaf angle distribution factor, or Xl
!!!  
!!!  ! begin aerodynamic variables
!!!  real(kind=8) dRbC   ! Rb Coefficient (c1 or cc1)
!!!  real(kind=8) dRdC   ! Rd Coefficient (c2 or cc2)
!!!  real(kind=8) dz0    ! Canopy roughness length
!!!  real(kind=8) dzp_disp     ! Zero plane displacement
!!!  real(kind=8) dG2    ! Ratio Ra (actual) to Ra (log-linear) for momentum
!!!  real(kind=8) dG3    ! Ratio Ra (actual) to Ra (log-linear) for heat transfer
!!!  real(kind=8) dCORB1 ! Non-nuetral correction for Ra between Ha and z2
!!!  real(kind=8) dCORB2 ! neutral value of RBB*U2^2 (RdC^2 for upper canopy)
!!!  real(kind=8) dHA    ! Canopy source height
!!!  real(kind=8) dsigma    ! Canopy source height
!!!  
!!!  ! begin internal constants
!!!  real(kind=8) dzwind ! Reference height for momentum (wind), meters
!!!  ! In practice, zwind is the height of measurement
!!!  real(kind=8) dzmet  ! Reference height for weather (meteorology), meters
!!!  ! In practice, zmet is the height of measurement
!!!  real(kind=8) dzs    ! Soil roughness length, meters
!!!  
!!!  ! Set values for Zwind and zmet. They must exceed the canopy top (Z2).
!!!  ! Set them for local domains, or use values of 80-100 m as for GCM.
!!!  
!!!  zwind=100.
!!!  zmet=100.
!!!
!!!  ! Set Soil roughness length (from Sellers et al Part II pp. 723).
!!!  zs=0.05
!!!
!!!  ! Convert to double precision (required for Sibx calculations)
!!!  dzwind=dble(zwind)
!!!  dzmet=dble(zmet)
!!!  dzs=dble(zs)
!!!  
!!!  dz2=dble(z2)
!!!  dz1=dble(z1)
!!!  dzc=dble(zc)
!!!  dLWidth=dble(LWidth)
!!!  dLLength=dble(LLength)
!!!  dLAI=dble(LAI)
!!!  dfVCover=dble(fVCover)
!!!  dChiL=dble(ChiL)
!!!  
!!!  dRbC=0.d0
!!!  dRdC=0.d0
!!!  dz0=0.d0
!!!  dzp_disp=0.d0
!!!  dG2=0.d0
!!!  dG3=0.d0
!!!  dCORB1=0.d0
!!!  dCORB2=0.d0
!!!
!!!! Calculate aerodynamic variables
!!!  call sibx(dzwind, dzmet, dz2, dzc, dz1, dzs, dLAI, dfVCover, dLLength, dLWidth, dChiL, dRbC, dRdC, dz0, dzp_disp,  &
!!!       dG2, dG3,dCORB1,dCORB2,dHA, dsigma)
!!!  
!!!  ! Convert aerodynamic variables back to single precision for output
!!!  RbC=dRbC
!!!  RdC=dRdC
!!!  z0=dz0
!!!  zp_disp=dzp_disp
!!!  G2=dG2
!!!  G3=dG3
!!!  CORB1=dCORB1
!!!  CORB2=dCORB2
!!!  HA=dHA
!!!  sigma=dsigma
!!!
!!!  return                                                                    
!!!end subroutine aero
!!!
!!!subroutine sibx(ZWIND, ZMET, Z2, ZC, Z1, ZS, LAI, fVCover,LLength, LWidth, ChiL, RBC, RDC, Z0, zp_disp, &
!!!     G2,G3,Corb1,Corb2,HA, sigma)
!!! 
!!!  !======================================================================
!!!  ! sibx calculates the aerodynamic properties of vegetation canopies
!!!  ! assuming a triangular distribution of leaf area density.
!!!  ! References:
!!!  !   SELLERS P.J. , Y. MINTZ, Y.S. SUD, A. DALCHER (1986) 'A SIMPLE
!!!  !         BIOSPHERE MODEL (SIB) FOR USE WITHIN GENERAL CIRCULATION
!!!  !         MODELS', ATMOS. SCI., 43, 6, 505-531.
!!!  !   SELLERS P.J. , W.J. SHUTTLEWORTH, J.L. DORMAN, A. DALCHER,
!!!  !         J.M. ROBERTS (1989) ' CALIBRATING THE SIMPLE BIOSPHERE
!!!  !         MODEL (SIB) FOR AMAZONIAN TROPICAL FOREST USING FIELD AND
!!!  !         REMOTE SENSING DATA. PART 1 , AVERAGE CALIBRATION WITH FIELD
!!!  !         DATA ', ( SEE APPENDIX I ), J. APPL. MET.
!!!  !   SHAW R.H. , A.R. PEREIRA (1982) 'AERODYNAMIC ROUGHNESS OF A
!!!  !         CANOPY: A NUMERICAL EXPERIMENT', AGRIC MET., 26, 51-65.
!!!  !-----------------------------------------------------------------------
!!!  
!!!  Implicit None
!!!
!!!  ! Continuous ecological variables passed in from Mapper
!!!  real(kind=8), intent (in) :: &
!!!       LAI , &   ! Total leaf area index (m^2 m^-2)
!!!       fVCover   ! Fractional vegetation cover
!!!  
!!!  ! Canopy morphology parameters passed in from Mapper
!!!  real(kind=8), intent (in) :: &
!!!       Z2 , &   ! Canopy top height (m)
!!!       Z1 , &   ! Canopy base height (m)
!!!       ZC , &   ! Canopy inflection  height (m)
!!!       LLength , & ! Leaf length
!!!       LWidth  , & ! Leaf width
!!!       ChiL    ! Leaf angle distribution factor
!!!  
!!!  ! Site-specific parameters set in subroutine aero
!!!  real(kind=8), intent (in) :: &
!!!       ZWIND , & ! Measurement height for wind (m)
!!!       ZMET  , & ! Measurement height for weather (m)
!!!       ZS       ! Soil roughness length (m)
!!!
!!!  ! Output aerodynamic variables
!!!  real(kind=8), intent (out) :: &
!!!       RBC , & ! Resistance coefficient between leaves and canopy airspace.
!!!       RDC , & ! Resistance coefficient between ground and canopy airspace.
!!!       Z0  , & ! Canopy roughness length
!!!       zp_disp, & ! Zero plane displacement
!!!       G2  , & ! Ratio Ra (actual) to Ra (log-linear) for momentum
!!!       G3  , & ! Ratio Ra (actual) to Ra (log-linear) for heat transfer
!!!       CORB1 ,&! Non-nuetral correction for Ra between Ha and z2
!!!       CORB2 ,&! neutral value of RBB*U2^2 (RdC^2 for upper canopy)
!!!       HA,   &  ! Canopy source height
!!!       sigma   ! Canopy mixing length
!!!
!!!  ! Control variables
!!!  integer IGCM ! flag for GCM output (= 1 for GCM )
!!!
!!!  ! Local variables
!!!  real(kind=8) :: &
!!!       FLOWB , &! Flow parameter based on leaf angle distribution
!!!       REYNO , &! Reynolds number
!!!       CDD   , &! Leaf drag coefficient
!!!       DL    , &! Mean leaf-area density (m^2 m^-3)
!!!       PS    , &! Leaf shelter factor
!!!       CS    , &! leaf Heat-mass transfer coefficient
!!!       U2FAC   ! Scaling factor for GCM Ratio U2-actual to U2-log-linear
!!!
!!!  ! IGCM=1 for GCM, =0 otherwise (affects resistance values).
!!!  IGCM=1
!!!  
!!!  ! calculate leaf drag coefficient (Sellers et al 1996, Eqn (B7))
!!!  FLOWB=1./3.141592*(1.-ChiL)
!!!  FLOWB=DMAX1(FLOWB, 0.00001D0)
!!!  REYNO=(LWidth+LLength)/2.*10000./0.15
!!!  CDD=1.328*2./DSQRT(REYNO)+0.45*FLOWB**1.6
!!!  
!!!  ! calculate shelter factor and heat-mass transfer coefficient
!!!  ! Cs is from Goudriaan 1977, 'Crop Micrometeorology: A Simulation Study'
!!!  DL=(LAI/(Z2-Z1))/fVCover
!!!  PS=1+DL**0.6
!!!  CS=90.*DSQRT(dble(LWidth))
!!!    
!!!  ! Calculate aerodynamic parameters
!!!  CALL MOMOPT ( Z2, ZC, Z1, ZS, ZWIND, ZMET, LAI, CDD, PS, CS, &
!!!       Z0, zp_disp, RBC, RDC, HA, G2, G3, CORB1, CORB2, U2FAC,sigma )
!!!  
!!!  ! re-scale aerodynamic resistances
!!!  IF (IGCM.EQ.0) U2FAC=1.
!!!  RBC=RBC/DSQRT(U2FAC)
!!!  RDC=RDC/U2FAC
!!!
!!!  RETURN
!!!END subroutine sibx

SUBROUTINE MOMOPT ( Z2, ZC, Z1, ZS, ZWIND, ZMET, LAI, CDD, PS, CS, &
     Z0, zp_disp, RBC, RDC, HA, G2, G3, CORB1, CORB2, U2FAC, sigma, ZL, COEFS )

  !======================================================================
  ! MOMOPT calculates the aerodynamic properties of a vegetation canopy
  ! using the first order closure model of Sellers et al. (1986), updated
  ! in Sellers et al. (1989).
  !-----------------------------------------------------------------------
  !     REFERENCES
  !         SELLERS P.J. , Y. MINTZ, Y.S. SUD, A. DALCHER (1986) 'A SIMPLE
  !         BIOSPHERE MODEL (SIB) FOR USE WITHIN GENERAL CIRCULATION
  !         MODELS', ATMOS. SCI., 43, 6, 505-531.
  !
  !         SELLERS P.J. , W.J. SHUTTLEWORTH, J.L. DORMAN, A. DALCHER,
  !         J.M. ROBERTS (1989) ' CALIBRATING THE SIMPLE BIOSPHERE
  !         MODEL (SIB) FOR AMAZONIAN TROPICAL FOREST USING FIELD AND
  !         REMOTE SENSING DATA. PART 1 , AVERAGE CALIBRATION WITH FIELD
  !         DATA ', ( SEE APPENDIX I ), J. APPL. MET.
  !
  !         SHAW R.H. , A.R. PEREIRA (1982) 'AERODYNAMIC ROUGHNESS OF A
  !         CANOPY: A NUMERICAL EXPERIMENT', AGRIC MET., 26, 51-65.
  !-----------------------------------------------------------------------
  
  implicit none
  
  ! Continuous ecological variables
  real(kind=8), intent (in) :: LAI   ! Total leaf area index
  
  ! Output aerodynamic variables
  real(kind=8), intent (out) :: &
       RBC , & ! Resistance coefficient between leaves and canopy airspace.
       RDC , & ! Resistance coefficient between ground and canopy airspace.
       Z0  , & ! Canopy roughness length
       zp_disp, & ! Zero plane displacement
       G2  , & ! Ratio Ra (actual) to Ra (log-linear) for momentum
       G3  , & ! Ratio Ra (actual) to Ra (log-linear) for heat transfer
       CORB1 , & ! Non-nuetral correction for Ra between Ha and z2
       CORB2 , & ! neutral value of RBB*U2^2 (RdC^2 for upper canopy)
       HA    , & ! CANOPY SOURCE HEIGHT FOR HEAT
       U2FAC, &  ! ratio U2 (actual) to U2 (log-linear). Used for GCM only.
       SIGMA     ! Mixing length (Sellers et al (1989), Eqns (A.3) & (A.5))

  ! Canopy morphology parameters
  real(kind=8), intent (in) :: &
       Z2 ,  &  ! Canopy top height
       Z1 ,  &  ! Canopy base height
       ZC      ! Canopy inflection  height

  ! Canopy flow regime constants
  real(kind=8), intent (in) :: &
       PS  ,  &  ! Leaf shelter factor
       CDD ,  &  ! Leaf drag coefficient
       CS       ! Heat-mass transfer coefficient

  ! Site-specific parameters set in subroutine aero
  real(kind=8), intent (in) :: &
       ZWIND , &! Measurement height for wind
       ZMET  , &! Measurement height for weather
       ZS      ! Soil roughness length

! Constants tuned to the model of Shaw and Pereira
  real(kind=8) G1    ! Ratio of eddy diffusivity extrapolated from a 
  ! log-linear profile (Km) to augmented eddy
  ! diffusivity at canopy top (Km*).
  real(kind=8) ZTZ0  ! Also G4.  Ratio of ZT to Z0.  ZT is the transition
  ! height between the perturbed surface layer and the
  ! inertial boundary layer. Z0 is canopy roughness length.    

! Local variables          
  real(kind=8) VKC        ! von Karman's coefficient
  real(kind=8) CDG        ! Ground drag coefficient
  real(kind=8) ZMAT(4,5)  ! Matrix of boundary and matching conditions
  ! Solve for YA and YB (Sellers et al (1989), Eqn (A.6))
  real(kind=8) WORK(4,5)  ! Working matrix for the non-destructive
  ! Gaussian elimination subroutine GAUSSD
  real(kind=8) COEFS(4) ! Storage vector for YA and YB
  ! See Equation (A.6) of Sellers et al (1989)
  real(kind=8) rl     ! Lower bisection bracket (dependent variable is negative)
  real(kind=8) rh     ! Higher bisection bracket (dependent variable is positive)
  real(kind=8) dxold  ! previous bisection step
  real(kind=8) dx     ! immediate bisection step
  real(kind=8) xacc   ! accuracy threshold
  integer MAXIT ! maximum number of iterations for convergence
  integer j     ! index of iterations
  real(kind=8) A0L    ! value of leaf area density at Z=0 for lower canopy
  real(kind=8) B0L    ! slope of leaf area density vs. Z for lower canopy
  real(kind=8) A0U    ! value of leaf area density at Z=0 for upper canopy
  real(kind=8) B0U    ! slope of leaf area density vs. Z for upper canopy
  real(kind=8) AT     ! Coefficients of Bessel equation     
  real(kind=8) BT     ! See Equation A.5 of Sellers et al (1989)
  real(kind=8) AL     
  real(kind=8) BL
  real(kind=8) YA     ! Bessel function term (first kind)
  real(kind=8) YB   ! Bessel function term (second kind)(Sellers (1989) Eqn (A.6))
  integer N     ! # boundary and matching conditions (#rows in ZMAT)
  integer NP1   ! Number of columns in ZMAT
  real(kind=8) Y      ! Error in D or Z0.  SIGMA and Z0 are roots in Y.
  real(kind=8) zp_disp1     ! Zero plane displacement from external scaling
  real(kind=8) zp_disp2     ! Zero plane displacement from concepts of momentum
  real(kind=8) ALP
  real(kind=8) BET
  real(kind=8) ZL     ! height of transition layer above ground
  real(kind=8) AA
  real(kind=8) BB
  real(kind=8) CC
  real(kind=8) BIGH
  real(kind=8) SQHH
  real(kind=8) V2
  real(kind=8) VX
  real(kind=8) ARG1   ! Three components of ra
  real(kind=8) ARG2   ! Sellers et al (1989), Eqn (A.26)
  real(kind=8) ARG3
  real(kind=8) RASTUB
  real(kind=8) ARGEX
  real(kind=8) Z0F
  real(kind=8) TOPG2
  real(kind=8) BOTG2
  real(kind=8) ZUP
  real(kind=8) TOPG3
  real(kind=8) BOTG3
  integer ISOL

  ! Assign Von Karmon's Constant
  DATA VKC/0.41/
  
  ! Reset Matrix
  ZMAT(:,:)=0.d0

  ! Assign calibration constants which force these results to 
  ! match those of Shaw and Pereira (1982)(Sellers et al (1989), App A)
  G1=1.449
  ZTZ0=11.785

  ! Calculate ground drag coefficient 
  ! (Shaw and Pereira, (1982); Sellers et al (1989), Eqn A.16)
  ! Cdg is weakly dependent on wind speed and varies with height.
  ! Height dependence can be eliminated knowing vertical velocity profile.
  CDG=(VKC/DLOG(Z1/ZS))**2

  ! Calculate coefficients of leaf area density (Sellers et al (1989), eqn A.4).
  !!! Ld(z) = a0l + bol * z,   z1 < z < zc
  !!! Ld(z) = a0u + b0u * z,   zc < z < z2
  CALL DENCAL(Z2,Z1,ZC,LAI,A0L,B0L,A0U,B0U)

  !-----------------------------------------------------------------------
  ! Calculate zero plane displacement height and mixing length
  !-----------------------------------------------------------------------
  ! iterate mixing length (sigma) to converge on displacement height (zp_disp)
  
  ! Set up initial mixing length (SIGMA).
  rl=1.0d0*Z2
  rh=0.001d0*Z2
  SIGMA=0.01d0*Z2
  dxold=rh-rl
  dx=dxold

  ! Set max number iterations and acceptance criteria
  MAXIT=200
  xacc=0.0001d0*Z2

  ! Iterate mixing length to converge on correct mixing length, SIGMA.
  ! Values for SIGMA, D, and Z0 should be of order Z2.
  ! See Shaw and Pereira 1982 for help with scaling.
  do j=1,MAXIT
     !       Calculate coefficients of the Bessel equation Ai and Bi.
     AT=2.*A0U*CDD/(PS*SIGMA)
     BT=2.*B0U*CDD/(PS*SIGMA)
     AL=2.*A0L*CDD/(PS*SIGMA)
     BL=2.*B0L*CDD/(PS*SIGMA)

     !       solve Bessel equation for Ya and Yb
     !       assign matrix of boundary and matching conditions
     !       scale u^2 = 1 at z = z2 as a boundary condition
     CALL WINVAL(Z2, AT, BT, YA, YB)
     ZMAT(1,1)=YA
     ZMAT(1,2)=YB
     ZMAT(1,5)=1.

     !       set u^2 to match at z = zc, inflection in leaf area density
     CALL WINVAL(ZC, AT, BT, YA, YB)
     ZMAT(2,1)=YA
     ZMAT(2,2)=YB

     CALL WINVAL(ZC, AL, BL, YA, YB)
     ZMAT(2,3)=-YA
     ZMAT(2,4)=-YB

     !       solve Bessel equation for dYa/dz and dYb/dz
     !       set d(u^2)/dz to match at z = zc
     CALL GRADVA(ZC, AT, BT, YA, YB)
     ZMAT(3,1)=YA
     ZMAT(3,2)=YB

     CALL GRADVA(ZC, AL, BL, YA, YB)
     ZMAT(3,3)=-YA
     ZMAT(3,4)=-YB
     
     !       set shear stress to match at z = z1 (tak added factor of 0.5)
     !       (Sellers et al (1989) eqn A.16)
     CALL WINVAL(Z1, AL, BL, YA, YB)

     CALL GRADVA(Z1, AL, BL, ALP, BET)
     ZMAT(4,3)=0.5*SIGMA*ALP-CDG*YA
     ZMAT(4,4)=0.5*SIGMA*BET-CDG*YB
     ZMAT(4,5)=0.

     !Continuity of dtau/dz in Z=Z1, that replaces the continuity of tau in Z=Z1 (done by nelsonvn)
!     ZMAT(4,3) = ALP - 2.0 * YA / (Z1 * LOG(Z1 / ZS))
!     ZMAT(4,4) = BET - 2.0 * YB / (Z1 * LOG(Z1 / ZS))

     N=4
     NP1=5
     
     !       Solve matrix of velocity and the vertical gradient of velocity
     !       with boundary conditions for coefficients alpha and beta using
     !       Gaussian elimination.
     CALL GAUSSD(ZMAT, N, NP1, COEFS, WORK)

     !       These two values of d depend indirectly on SIGMA and on u(z).
     !       (Sellers et al (1989), Eqs A.13, A.17)
     CALL DEVAL1(COEFS, AT, BT, Z2, G1, SIGMA, VKC, zp_disp1)
     
     CALL DEVAL2(COEFS, AT, BT, AL, BL, A0U, B0U, A0L,B0L,ZS,Z2, &
          ZC, Z1, PS, CDD, SIGMA, zp_disp2 )


     !       calculate difference between 2 displacement heights
     Y=(zp_disp1-zp_disp2)

     !       determine if error meets convergence criteria 
     if (dabs(Y).lt.xacc) goto 100

     !       bracket root 
     if (Y.lt.0.) then
        rl=SIGMA
     else
        rh=SIGMA
     endif

     !       Take bisection
     dxold=dx
     dx=0.5*(rh-rl)
     SIGMA=rl+dx
     SIGMA=DMAX1(SIGMA, 0.0000001D0)
  enddo ! zero plane displacement


  ! if error less than criteria take average displacement height
100 zp_disp=(zp_disp1+zp_disp2)/2.

  !-----------------------------------------------------------------------
  ! calculate Z0 and G2 as a function of G1
  !-----------------------------------------------------------------------
  ! G1 assumes linear variation of KM between Z2 and ZL such that
  !   KM(ACTUAL)=G1*KM(LOG-LINEAR) AT Z=Z2
  !   KM(ACTUAL)=KM(LOG-LINEAR)    AT Z=ZL
  !
  ! Integration follows eqns A.10 and A.11 in Sellers et al. (1989).
  ! z0 follows from Equations (A.18), (A.19) and (A.20).

  ! solve for Z0 using bisection method
  rl=0.50d0*Z2
  rh=0.01d0*Z2
  Z0=0.05d0*Z2
  dxold=rl-rh
  dx=dxold

  MAXIT=200
  xacc=0.00001d0

  do j=1,MAXIT
     ZL=Z2+ZTZ0*Z0   ! eqn A.8
     AA=1.-G1
     BB=G1*ZL-Z2+zp_disp*(G1-1.)
     CC=-zp_disp*(G1*ZL-Z2)
     BIGH=CC/AA-BB*BB/(4.*AA*AA)  ! eqn A.9
     SQHH=DSQRT(DABS(BIGH))
     V2=Z2+BB/(2.*AA)
     VX=ZWIND+BB/(2.*AA)
     IF(ZL.LT.ZWIND) VX=ZL+BB/(2.*AA)
     IF (BIGH.LT.0.) GO TO 200
     ARG1=DATAN(VX/SQHH)/SQHH
     ARG2=DATAN(V2/SQHH)/SQHH
     ISOL=1
     GO TO 300
200  ARG1=DLOG(DABS(VX-SQHH)/(VX+SQHH))/(2.*SQHH)
     ARG2=DLOG(DABS(V2-SQHH)/(V2+SQHH))/(2.*SQHH)
     ISOL=2
300  ARG3=0.
     IF(ZL.LT.ZWIND) ARG3=DLOG((ZWIND-zp_disp)/(ZL-zp_disp))
     RASTUB=(ZL-Z2)/AA*(ARG1-ARG2)+ARG3
     ARGEX = RASTUB + G1*VKC*VKC*(Z2-zp_disp)/SIGMA
     Z0F = (ZWIND-zp_disp) * DEXP(-ARGEX)
     Z0F = DMAX1( Z0F, 0.00001D0 )
     
     Y=Z0F-Z0
     if (dabs(Y).lt.xacc) goto 330 
     if (Y.lt.0.) then
        rl=Z0F
     else
        rh=Z0F
     endif
     
     dxold=dx
     dx=0.5*(rh-rl)
     Z0=rl+dx      
     Z0=DMAX1(Z0, 0.00001D0)
  enddo

330 Z0=Z0F

  TOPG2=RASTUB-ARG3
  ZUP=ZWIND
  IF (ZL.LT.ZWIND) ZUP=ZL
  BOTG2=DLOG((ZUP-zp_disp)/(Z2-zp_disp))
  G2=TOPG2/BOTG2
  
  !-----------------------------------------------------------------------
  ! Calculate G3, Rbc, Rdc, Corb1, Corb2, Ha
  ! assume linear variation of Kh between z2 and zl as for G2
  !-----------------------------------------------------------------------
  AA=1.-G1
  BB=G1*ZL-Z2+zp_disp*(G1-1.)
  CC=-zp_disp*(G1*ZL-Z2)
  BIGH=CC/AA-BB*BB/(4.*AA*AA)
  SQHH=DSQRT(DABS(BIGH))
  V2=Z2+BB/(2.*AA)
  VX=ZMET+BB/(2.*AA)
  IF(ZL.LT.ZMET) VX=ZL+BB/(2.*AA)
  IF (BIGH.LT.0.) GO TO 400
  ARG1=DATAN(VX/SQHH)/SQHH
  ARG2=DATAN(V2/SQHH)/SQHH
  ISOL=1
  GO TO 500
400 ARG1=DLOG(DABS(VX-SQHH)/(VX+SQHH))/(2.*SQHH)
  ARG2=DLOG(DABS(V2-SQHH)/(V2+SQHH))/(2.*SQHH)
  ISOL=2
500 ARG3=0.
  IF(ZL.LT.ZMET) ARG3=DLOG((ZMET-zp_disp)/(ZL-zp_disp))
  RASTUB=(ZL-Z2)/AA*(ARG1-ARG2)+ARG3
  
  TOPG3=RASTUB-ARG3
  ZUP=ZMET
  IF (ZL.LT.ZMET) ZUP=ZL
  BOTG3=DLOG((ZUP-zp_disp)/(Z2-zp_disp))
  G3=TOPG3/BOTG3

  CALL CALFAC ( Z2, ZC, Z1, CDG, COEFS, AT, BT, AL, BL, A0U, B0U, A0L, B0L, SIGMA, CS, PS, &
       RBC, HA, RDC, CORB1, CORB2 )

  ARG1=DLOG((ZL-zp_disp)/Z0)
  ARG2=DLOG((ZL-zp_disp)/(Z2-zp_disp))
  U2FAC=(ARG1-G2*ARG2)/(ARG1-ARG2)

  RETURN
END SUBROUTINE MOMOPT

SUBROUTINE DENCAL( Z2, Z1, ZC, LAI, A0L, B0L, A0U, B0U)
  ! Calculates the leaf area density in the canopy assuming a triangular 
  ! distribution.  The total area under the curve is the leaf area index.
  ! The inflection point divides the canopy into an upper and lower portion.
  ! DENCAL calculates the slopes and z-intercepts for both the upper and 
  ! lower leaf area density curves.
  
  implicit none
  
  ! input variables
  real(kind=8) Z2    ! height of canopy top
  real(kind=8) Z1    ! height of canopy base
  real(kind=8) ZC    ! height of canopy inflection point
  real(kind=8) LAI   ! Leaf Area Index

  ! output variables
  real(kind=8) A0L   ! value of leaf area density at Z=0 for lower canopy
  real(kind=8) B0L   ! slope of leaf area density vs. Z for lower canopy
  real(kind=8) A0U   ! value of leaf area density at Z=0 for lupper canopy
  real(kind=8) B0U   ! slope of leaf area density vs. Z for upper canopy

  ! internal variables
  real(kind=8) ZLMAX ! maximum value of leaf area density at ZC
  real(kind=8) ZLMIN ! minimum value of leaf area density at Z2 and Z1

  ! Set minimum leaf area density at Z1 and Z2
  ! This assures the leaf area density is never zero.
  ZLMIN=0.001

  ! calculate leaf area density at inflection point, ZC
  ! Set LAI equal to the area of a triangle and solve for its height
  ZLMAX=2.*LAI/(Z2-Z1)-ZLMIN

  !  calculate slopes of upper and lower leaf area density curves
  B0L=(ZLMAX-0.001)/(ZC-Z1)
  B0U=(ZLMAX-0.001)/(ZC-Z2)

  ! calculate z-intercepts (Z=0) for upper and lower leaf area density curves
  A0L=ZLMAX-B0L*ZC
  A0U=ZLMAX-B0U*ZC

  RETURN
END SUBROUTINE DENCAL

SUBROUTINE WINVAL( Z, A, B, YA, YB)

  implicit none
  
  real(kind=8) Z
  real(kind=8) A
  real(kind=8) B
  real(kind=8) YA
  real(kind=8) YB
  real(kind=8) ZNU
  real(kind=8) ARG
  real(kind=8) ORDER
  real(kind=8) BI
  real(kind=8) BK
  real(kind=8) BIP
  real(kind=8) BKP

  ZNU=1./DABS(B)**(2./3.)*(A+B*Z)   !!! SE89 A6, eq 4
  ARG=2./3.*ZNU**(3./2.)            !!! Se89 A6, eq 5
  ORDER=1./3.

  call bessik(ARG,ORDER,BI,BK,BIP,BKP)

  ! The XMIN parameter in bessik is 2.0d0
  ! If argument > 2 then use exponential scaling
  if (ARG.lt.(2.0d0)) then
     YA=DSQRT(ZNU)*BI               !!! SE89 A6, eq 2
     YB=DSQRT(ZNU)*BK               !!! SE89 A6, eq 3
  else if (ARG.ge.(2.0d0)) then
     YA=DSQRT(ZNU)*BI*DEXP(ARG)
     YB=DSQRT(ZNU)*BK*DEXP(-ARG)
  endif
  
  RETURN
END SUBROUTINE WINVAL

SUBROUTINE GRADVA(Z, A, B, YA, YB)

  implicit none

  real(kind=8) Z
  real(kind=8) A
  real(kind=8) B
  real(kind=8) YA
  real(kind=8) YB
  real(kind=8) ZNU
  real(kind=8) ARG
  real(kind=8) ORDER
  real(kind=8) BI13
  real(kind=8) BK13
  real(kind=8) BIP
  real(kind=8) BKP
  real(kind=8) BI43
  real(kind=8) BK43
  real(kind=8) BI23
  real(kind=8) BK23

  ZNU=1./DABS(B)**(2./3.)*(A+B*Z)   !!! SE89 A6, eq 4
  ARG=2./3.*ZNU**(3./2.)            !!! SE89 A6, eq 5

  ORDER=1./3.
  call bessik(ARG,ORDER,BI13,BK13,BIP,BKP)

  ORDER=4./3.
  call bessik(ARG,ORDER,BI43,BK43,BIP,BKP)

  ORDER=2./3.
  call bessik(ARG,ORDER,BI23,BK23,BIP,BKP) 

  if (ARG.lt.(2.0d0)) then
     YA=(1./DSQRT(ZNU))*BI13+ZNU*BI43
     YA=YA/2.*(DABS(B)**(1./3.))*DSIGN(1.D0,B)  !!! SE89 A7, eq 2. Is the factor '2' right?
     YB=-ZNU*BK23
     YB=YB/2.*(DABS(B)**(1./3.))*DSIGN(1.D0,B)  !!! SE89 A7, eq 3. Is this formula right?
  else if (ARG.ge.(2.0d0)) then
     YA=(1./DSQRT(ZNU))*BI13*DEXP(ARG)+ZNU*BI43*DEXP(ARG)
     YA=YA/2.*(DABS(B)**(1./3.))*DSIGN(1.D0,B)
     YB=-ZNU*BK23*DEXP(-ARG)
     YB=YB/2.*(DABS(B)**(1./3.))*DSIGN(1.D0,B)
  endif
  
  RETURN
END SUBROUTINE GRADVA

SUBROUTINE DEVAL1(COEFS, AT, BT, Z2, G1, SIGMA, VKC, zp_disp1)
  ! (Sellers et al (1989), eqn A.13)

  implicit none

  real(kind=8) COEFS(4)
  real(kind=8) AT     
  real(kind=8) BT     
  real(kind=8) Z2    
  real(kind=8) G1
  real(kind=8) SIGMA
  real(kind=8) VKC
  real(kind=8) zp_disp1
  real(kind=8) A1
  real(kind=8) B1
  real(kind=8) Y
  real(kind=8) A2
  real(kind=8) B2
  real(kind=8) YDASH

  CALL WINVAL(Z2, AT, BT, A1, B1)
  Y=COEFS(1)*A1+COEFS(2)*B1   !ERROR! Y=U^2 (commented by nelsonvn)

  CALL GRADVA(Z2,AT,BT,A2,B2)
  
  YDASH=COEFS(1)*A2+COEFS(2)*B2   !ERROR! Y'=(U^2)' and we need U' (commented by nelsonvn)

!   Y = SQRT(Y)   !Correct value of U(Z2) (done by nelsonvn)
!   YDASH = 0.5 * YDASH / Y   !Correct value of U'(Z2) (done by nelsonvn)

  zp_disp1=Z2-1./(G1*VKC)*DSQRT(Y*SIGMA/YDASH)
  
  RETURN
END SUBROUTINE DEVAL1

SUBROUTINE DEVAL2(COEFS, AT, BT, AL, BL, A0U, B0U, A0L, B0L, &
     ZS, Z2, ZC, Z1, PS, CDD, SIGMA, zp_disp2 )
  ! (Sellers et al (1989), eqn A.17)

  implicit none

  real(kind=8) COEFS(4)
  real(kind=8) AT     
  real(kind=8) BT     
  real(kind=8) AL     
  real(kind=8) BL     
  real(kind=8) A0L    
  real(kind=8) B0L    
  real(kind=8) A0U    
  real(kind=8) B0U    
  real(kind=8) ZS
  real(kind=8) Z2     
  real(kind=8) Z1     
  real(kind=8) ZC     
  real(kind=8) PS
  real(kind=8) CDD
  real(kind=8) SIGMA
  real(kind=8) zp_disp2
  real(kind=8) TOP
  real(kind=8) BOT
  real(kind=8) A
  real(kind=8) B
  real(kind=8) ZUP
  real(kind=8) ZLO
  real(kind=8) A0
  real(kind=8) B0
  real(kind=8) YA
  real(kind=8) YB
  integer ILEVEL
  real(kind=8) ZN1
  real(kind=8) ZN2
  real(kind=8) ALP
  real(kind=8) BET
  real(kind=8) VAL
  real(kind=8) TOR1
  real(kind=8) BOTTOR

  TOP=0.
  BOT=0.
  A=AT
  B=BT
  ZUP=Z2
  ZLO=ZC
  A0=A0U
  B0=B0U
  YA=COEFS(1)
  YB=COEFS(2)

  DO ILEVEL=1,2
     CALL YINTNU (A, B, A0, B0,  ZUP, ZLO, YA, YB, ZN1, ZN2)
     TOP=TOP+ZN2
     BOT=BOT+ZN1
     A=AL
     B=BL
     ZUP=ZC
     ZLO=Z1
     A0=A0L
     B0=B0L
     YA=COEFS(3)
     YB=COEFS(4)
  ENDDO

  
  CALL WINVAL(Z1, AL, BL, ALP, BET)
  
  VAL=COEFS(3)*ALP+COEFS(4)*BET   !Y(Z1) VALUE. Y=U^2 (commented by nelsonvn)

!   WRITE(*,*) 'VAL', VAL   !Y(Z1) (nelsonvn)

  !TOR1 is RHS of SE89 A16 (commented by nelsonvn)
  TOR1=((0.41*VAL)/log(Z1/ZS))**2.d0   !ERROR! VAL=U(Z1)^2, SO VAL^2=U(Z1)^4 (commented by nelsonvn)
  BOTTOR=TOR1*PS/CDD*SIGMA   !ERROR! WHY DOES 'SIGMA' APPEAR HERE AND NOT 'LD'? (commented by nelsonvn)

!  TOR1=VAL*(0.41/log(Z1/ZS))**2.d0   !Correct formula of SE89 A16 (done by nelsonvn)
!  YA = (PS - 1.0) ** (1.0 / 0.6)   !Correct DL value (done by nelsonvn)
!  BOTTOR=TOR1*PS/CDD/YA   !Correct formula of SE89 A17 (done by nelsonvn)

  zp_disp2=TOP/(BOT+BOTTOR)
  
  RETURN
END SUBROUTINE DEVAL2

SUBROUTINE YINTNU(A, B, A0, B0, ZUP, ZLO, YA, YB, ZN1, ZN2 )
  ! calculates integrals of U*U*LD and U*U*LD*Z numerically

  implicit none

  real(kind=8) A
  real(kind=8) B
  real(kind=8) A0
  real(kind=8) B0
  real(kind=8) ZUP
  real(kind=8) ZLO
  real(kind=8) YA
  real(kind=8) YB
  real(kind=8) ZN1
  real(kind=8) ZN2
  real(kind=8) ZIN
  real(kind=8) DZ
  real(kind=8) A1
  real(kind=8) B1
  real(kind=8) VAL
  real(kind=8) LAIZ  ! Leaf Area Indez as a function of Z in canopy
  integer I

  ! calculate integrals of U*U*LD and U*U*LD*Z numerically
  ZN1=0.
  ZN2=0.
  DZ=(ZUP-ZLO)/20.
  
  DO I=1,20
     ZIN=I*DZ-DZ*0.5+ZLO
     LAIZ=A0+B0*ZIN
     CALL WINVAL ( ZIN, A, B, A1, B1 )
     VAL=YA*A1+YB*B1
     ZN1=ZN1+VAL*LAIZ*DZ
     ZN2=ZN2+VAL*LAIZ*ZIN*DZ
  ENDDO

  RETURN
END SUBROUTINE YINTNU

SUBROUTINE CALFAC ( Z2, ZC, Z1, CDG, COEFS, AT, BT, AL, BL, A0U, B0U, A0L, B0L, SIGMA, CS, PS, &
     RBC, HA, RDC, CORB1, CORB2 )
  
  implicit none

  real(kind=8) Z2      ! canopy top (m)
  real(kind=8) Z1      ! canopy base (m)
  real(kind=8) ZC      ! leaf area density inflection height (m)
  real(kind=8) CDG     ! drag coeficient for ground
  real(kind=8) COEFS(4)! boundary condition Bessel function coefficients
  real(kind=8) AT      ! 1st Bessel function coefficient for top canopy layer
  real(kind=8) BT      ! 2nd Bessel function coefficient for top canopy layer 
  real(kind=8) AL      ! 1st Bessel function coefficient for bottom canopy layer
  real(kind=8) BL      ! 2nd Bessel function coefficient for bottom canopy layer
  real(kind=8) A0L     ! 1st coefficient for leaf area density bottom canopy layer
  real(kind=8) B0L     ! 2nd coefficient for leaf area density bottom canopy layer
  real(kind=8) A0U     ! 1st coefficient for leaf area density top canopy layer
  real(kind=8) B0U     ! 2nd coefficient for leaf area density top canopy layer
  real(kind=8) SIGMA   !mixing length
  real(kind=8) CS      ! leaf Heat-mass transfer coefficient
  real(kind=8) PS      ! leaf shelter factor
  real(kind=8) RBC     ! rb resistance coefficient 
  real(kind=8) HA      ! canopy source height
  real(kind=8) RDC     ! rd resistance coefficient
  real(kind=8) CORB1   ! TBD
  real(kind=8) CORB2   ! resistance across whole canopy coefficient squared
  real(kind=8) HEIGHT(40) ! integration heights in canopy
  real(kind=8) UU(40)     ! velocity squared
  real(kind=8) GRADU(40)  ! (?) vertical gradient of velocity
  real(kind=8) RBADD(40)  ! cumulative rb in canopy
  real(kind=8) A       ! 1st Bessel function coefficient
  real(kind=8) B       ! 2nd Bessel function coefficient
  real(kind=8) ZUP     ! canopy layer top
  real(kind=8) ZLO     ! canopy layer bottom
  real(kind=8) A0      ! 1st coefficient for leaf area density
  real(kind=8) B0      ! 2nd coefficient for leaf area density
  real(kind=8) YA      ! 1st boundary condition Bessel function coefficient
  real(kind=8) YB      ! 2nd boundary condition Bessel function coefficient
  real(kind=8) DZ      ! integration layer thickness
  real(kind=8) ZIN     ! integration height
  real(kind=8) A1      ! 1st Bessel function value
  real(kind=8) B1      ! 2nd Bessel function value
  real(kind=8) LAIZ    ! leaf area density at height zin
  real(kind=8) TARGET  ! 1/2 point in cumulative rb for canopy
  real(kind=8) A2      ! 1st Bessel function derivative value
  real(kind=8) B2      ! 2nd Bessel function derivative value
  real(kind=8) ZINK    ! numerical integration slice
  real(kind=8) WEIGHT  ! weight factor for exact Ha in integration layer
  real(kind=8) U1      ! velocity at canopy base
  real(kind=8) DZSAVE  ! integration layer thickness for Corb1
  integer ILEVEL ! canopy laer index
  integer INDEX  ! integration level index
  integer IZ     ! integration index
  integer ISTART ! 1st level index for integration

  !--------------------------------------------------------
  ! calculate rb resistance coefficient (RBC)
  !--------------------------------------------------------
  ! Also calculate velocity, velocity gradient, and rb 
  ! as a function of height within canopy
  
  ! initialize RBC
  RBC=0.

  ! set constants for lower canopy layer
  A=AL
  B=BL
  ZUP=ZC
  ZLO=Z1
  A0=A0L
  B0=B0L
  YA=COEFS(3)
  YB=COEFS(4)

  ! Integrate upward through canopy
  DO ILEVEL=1,2
     DZ=(ZUP-ZLO)/20.
     DO IZ=1,20
        INDEX=(ILEVEL-1)*20+IZ
        ZIN=ZLO+DZ*(IZ-0.5)
        HEIGHT(INDEX)=ZIN

        CALL WINVAL( ZIN, A, B, A1, B1)
        UU(INDEX)=YA*A1+YB*B1

        CALL GRADVA (ZIN, A, B, A2, B2)
        GRADU(INDEX)=(YA*A2+YB*B2)/DSQRT(UU(INDEX))

        LAIZ=A0+B0*ZIN
        RBC=RBC+LAIZ*DSQRT(DSQRT(UU(INDEX)))*DZ
        RBADD(INDEX)=RBC
     ENDDO
     A=AT
     B=BT
     ZUP=Z2
     ZLO=ZC
     A0=A0U
     B0=B0U
     YA=COEFS(1)
     YB=COEFS(2)
  ENDDO

  ! calculate 1/2 value of RBC for calculating source height
  TARGET=RBC/2.

  ! final manipulation of RBC
  RBC=1./RBC*CS*PS

  !--------------------------------------------------------
  ! calculate source height, Ha
  !--------------------------------------------------------
  DO IZ=2,40
     IF (TARGET.GT.RBADD(IZ)) GO TO 100
     WEIGHT=1.-(RBADD(IZ)-TARGET)/(RBADD(IZ)-RBADD(IZ-1))
     HA=HEIGHT(IZ-1)+(HEIGHT(IZ)-HEIGHT(IZ-1))*WEIGHT
     GO TO 200
100  CONTINUE
  ENDDO
200 CONTINUE

  !--------------------------------------------------------
  ! Calculate RDC
  !--------------------------------------------------------
  ! Calculate wind speed at canopy base
  CALL WINVAL (Z1, AL, BL, A1, B1)
  U1=DSQRT(COEFS(3)*A1+COEFS(4)*B1)

  ! Calculate total resistance below canopy base
  RDC=1./(U1*CDG)

  ! Calculate resistance from canopy base to level below Ha
  DO IZ=1,40
     DZ=HEIGHT(2)-HEIGHT(1)
     IF (IZ.GT.21) DZ=HEIGHT(40)-HEIGHT(39)
     ZINK=1./(SIGMA*DSQRT(UU(IZ)))*DZ
     IF (HEIGHT(IZ+1).GT.HA) GO TO 300
     RDC=RDC+ZINK
  ENDDO

  ! Calculate resistance from next lowest integration layer to Ha
300 RDC=RDC+WEIGHT*ZINK
  DZSAVE=DZ
      
  !--------------------------------------------------------
  ! Calculate Corb2
  !--------------------------------------------------------
  ! Calculate resistance from Ha to next highest integration level 
  CORB2=ZINK*(1.-WEIGHT)

  ! set starting point at next level above Ha
  ! (previous loop sets value of iz to that level containing Ha)
  ISTART=IZ+1

  ! Integrate resistance from level above Ha to canopy top
  DO IZ=ISTART,39
     DZ=HEIGHT(2)-HEIGHT(1)
     IF (IZ.GT.21) DZ=HEIGHT(40)-HEIGHT(39)
     CORB2=CORB2+1./(SIGMA*DSQRT(UU(IZ)))*DZ
  ENDDO

  ! square the final value, I don't know why
  CORB2=CORB2**2

  !--------------------------------------------------------
  ! Calculate Corb1
  !--------------------------------------------------------
  ! Calculate gradient from Ha to next highest integration level
  CORB1=GRADU(ISTART-1)**2*(1.-WEIGHT)*DZSAVE

  ! Integrate gradient from level above Ha to canopy top
  ! (istart already set to next level above Ha)
  DO IZ=ISTART,39
     DZ=HEIGHT(2)-HEIGHT(1)
     IF (IZ.GT.21) DZ=HEIGHT(40)-HEIGHT(39)
     CORB1=CORB1+GRADU(IZ)**2*DZ
  ENDDO

  ! further manipulate Corb1, I don't know why
  CORB1=CORB1/(Z2-HA)
  CORB1=9.*9.81/(1010.*1.2*CORB1)

  RETURN
END SUBROUTINE CALFAC

SUBROUTINE GAUSSD(A,N,NP1,X,WORK)
  !     SOLVE A LINEAR SYSTEM BY GAUSSIAN ELIMINATION.  DEVELOPED BY
  !     DR. CHIN-HOH MOENG.  A IS THE MATRIX OF COEFFICIENTS, WITH THE
  !     VECTOR OF CONSTANTS APPENDED AS AN EXTRA COLUMN.  X IS THE VECTOR
  !     CONTAINING THE RESULTS.  THE INPUT MATRIX IS NOT DESTROYED.

  implicit none

  real(kind=8) A(4,5)
  real(kind=8) WORK(4,5)
  real(kind=8) X(4)
  real(kind=8) R
  integer N
  integer NP1
  integer I
  integer J
  integer K
  integer L

  DO I=1,N
     DO J=1,NP1
        WORK(I,J)=A(I,J)
     ENDDO
  ENDDO
  
  DO I=2,N
     DO J=I,N
        R=WORK(J,I-1)/WORK(I-1,I-1)   
        DO K=1,NP1
           WORK(J,K)=WORK(J,K)-R*WORK(I-1,K)
        ENDDO
     ENDDO
  ENDDO

  DO I=2,N
     K=N-I+2
     R=WORK(K,NP1)/WORK(K,K)
     DO J=I,N
        L=N-J+1
        WORK(L,NP1)=WORK(L,NP1)-R*WORK(L,K)
     ENDDO
  ENDDO
  
  DO I=1,N
     X(I)=WORK(I,NP1)/WORK(I,I)
  ENDDO
  RETURN
END SUBROUTINE GAUSSD

subroutine bessik(x,xnu,ri,rk,rip,rkp)

  implicit none

  real(kind=8) ri
  real(kind=8) rip
  real(kind=8) rk
  real(kind=8) rkp
  real(kind=8) x
  real(kind=8) xnu
  real(kind=8), parameter :: EPS=1.e-16
  real(kind=8), parameter :: FPMIN=1.e-30
  integer, parameter :: MAXIT=10000
  real(kind=8), parameter :: XMIN=2.
  real(kind=8), parameter :: PI=3.141592653589793d0
       
  Integer i
  Integer l
  Integer nl
  real(kind=8) a
  real(kind=8) a1
  real(kind=8) b
  real(kind=8) c
  real(kind=8) d
  real(kind=8) del
  real(kind=8) del1
  real(kind=8) delh
  real(kind=8) dels
  real(kind=8) e
  real(kind=8) f
  real(kind=8) fact
  real(kind=8) fact2
  real(kind=8) ff
  real(kind=8) gam1
  real(kind=8) gam2
  real(kind=8) gammi
  real(kind=8) gampl
  real(kind=8) h
  real(kind=8) p
  real(kind=8) pimu
  real(kind=8) q
  real(kind=8) q1
  real(kind=8) q2
  real(kind=8) qnew
  real(kind=8) ril
  real(kind=8) ril1
  real(kind=8) rimu
  real(kind=8) rip1
  real(kind=8) ripl
  real(kind=8) ritemp
  real(kind=8) rk1
  real(kind=8) rkmu
  real(kind=8) rkmup
  real(kind=8) rktemp
  real(kind=8) s
  real(kind=8) sum
  real(kind=8) sum1
  real(kind=8) x2
  real(kind=8) xi
  real(kind=8) xi2
  real(kind=8) xmu
  real(kind=8) xmu2

! error check
  if (x.le.0..or.xnu.lt.0) print*, 'bad arguments in bessik'
       
  nl=int(xnu+.5d0)
  xmu=xnu-nl
  xmu2=xmu*xmu
  xi=1.d0/x
  xi2=2.d0*xi
  h=xnu*xi
  if(h.lt.FPMIN) h=FPMIN
  b=xi2*xnu
  d=0.d0
  c=h
  do i=1,MAXIT
     b=b+xi2
     d=1.d0/(b+d)
     c=b+1.d0/c
     del=c*d
     h=del*h
     if(abs(del-1.d0).lt.EPS) goto 1
  enddo

  print*, 'x too large in bessik; try asymptotic expansion'

1 continue

  ril=FPMIN
  ripl=h*ril
  ril1=ril
  rip1=ripl
  fact=xnu*xi
  
  do l=nl,1,-1
     ritemp=fact*ril+ripl
     fact=fact-xi
     ripl=fact*ritemp+ril
     ril=ritemp
  enddo

  f=ripl/ril
  
  if(x.lt.XMIN) then
     x2=.5d0*x
     pimu=PI*xmu
     if(abs(pimu).lt.EPS) then
        fact=1.0d0
     else
        fact=pimu/sin(pimu)
     endif
     d=-log(x2)
     e=xmu*d
     if(abs(e).lt.EPS) then
        fact2=1.d0
     else
        fact2=sinh(e)/e
     endif
     call beschb(xmu,gam1,gam2,gampl,gammi)
     ff=fact*(gam1*cosh(e)+gam2*fact2*d)
     sum=ff
     e=exp(e)
     p=0.5d0*e/gampl
     q=0.5d0/(e*gammi)
     c=1.d0
     d=x2*x2
     sum1=p

     do i=1,MAXIT
        ff=(i*ff+p+q)/(i*i-xmu2)
        c=c*d/i
        p=p/(i-xmu)
        q=q/(i+xmu)
        del=c*ff
        sum=sum+del
        del1=c*(p-i*ff)
        sum1=sum1+del1
        if(abs(del).lt.abs(sum)*EPS) goto 2
     enddo
     print*, 'bessk series failed to converge'
2    continue
     
     rkmu=sum
     rk1=sum1*xi2
     
  else
     b=2.d0*(1.d0+x)
     d=1.d0/b
     delh=d
     h=delh
     q1=0.d0
     q2=1.d0
     a1=.25d0-xmu2
     c=a1
     q=c
     a=-a1
     s=1.d0+q*delh
     
     do i=2,MAXIT
        a=a-2*(i-1)
        c=-a*c/i
        qnew=(q1-b*q2)/a
        q1=q2
        q2=qnew
        q=q+c*qnew
        b=b+2.d0
        d=1.d0/(b+a*d)
        delh=(b*d-1.d0)*delh
        h=h+delh
        dels=q*delh
        s=s+dels
        if(abs(dels/s).lt.EPS) goto 3
     enddo
     
     print*, 'bessik: failure to converge in CF2'
3    continue
     
     h=a1*h
     
     ! Exponentially scaled for all four functions
     rkmu=sqrt(PI/(2.d0*x))/s
     rk1=rkmu*(xmu+x+.5d0-h)*xi
     
  endif
  
  rkmup=xmu*xi*rkmu-rk1
  rimu=xi/(f*rkmu-rkmup)
  ri=(rimu*ril1)/ril
  rip=(rimu*rip1)/ril
  
  do i=1,nl
     rktemp=(xmu+i)*xi2*rk1+rkmu
     rkmu=rk1
     rk1=rktemp
  enddo
  
  rk=rkmu
  rkp=xnu*xi*rkmu-rk1

  return
END subroutine bessik

subroutine beschb(x,gam1,gam2,gampl,gammi)
  
  implicit none

  Integer NUSE1
  Integer NUSE2
  real(kind=8) gam1
  real(kind=8) gam2
  real(kind=8) gammi
  real(kind=8) gampl
  real(kind=8) x
  Parameter (NUSE1=7,NUSE2=8)
  
  real(kind=8) xx
  real(kind=8) c1(7)
  real(kind=8) c2(8)
  real(kind=8) chebev
  SAVE c1
  SAVE c2

  data c1/ -1.142022680371168d0, 6.5165112670737d-3, 3.087090173086d-4, -3.4706269649d-6, 6.9437664d-9, &
       3.67795d-11, -1.356d-13/
  
  data c2/ 1.843740587300905d0, -7.68528408447867d-2, 1.2719271366546d-3, -4.9717367042d-6, -3.31261198d-8, &
       2.423096d-10, -1.702d-13, -1.49d-15 /

  xx=8.d0*x*x-1.d0
  gam1=chebev(-1.0d0,1.0d0,c1,NUSE1,xx)
  gam2=chebev(-1.0d0,1.0d0,c2,NUSE2,xx)
  gampl=gam2-x*gam1
  gammi=gam2+x*gam1

  return
END subroutine beschb

function chebev(a,b,c,m,x)

  implicit none

  Integer m
  real(kind=8) chebev
  real(kind=8) a
  real(kind=8) b
  real(kind=8) x
  real(kind=8) c(m)
  Integer j
  real(kind=8) d
  real(kind=8) dd
  real(kind=8) sv
  real(kind=8) y
  real(kind=8) y2

  ! Error check
  if ((x-a)*(x-b).gt.0) print*, 'x not in range in chebev'

  d=0.
  dd=0.
  y=(2.*x-a-b)/(b-a)
  y2=2.*y

  do j=m,2,-1
     sv=d
     d=y2*d-dd+c(j)
     dd=sv
  enddo

  chebev=y2*d-dd+0.5*c(1)

  return 
END function chebev
