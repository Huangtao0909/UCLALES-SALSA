MODULE mo_salsa_coagulation_kernels
  USE mo_salsa_types, ONLY : aero, cloud, precp, ice, & 
                             iaero, faero, icloud, fcloud, iprecp, fprecp,  &
                             iice, fice
  USE mo_submctl, ONLY : nbins, ncld, nprc, nice, spec, pi6,  &
                         lscgaa, lscgcc, lscgpp, lscgii,      & 
                         lscgca, lscgpa, lscgia,              & 
                         lscgpc, lscgic,                      & 
                         lscgip 
  USE classSection, ONLY : Section
  USE mo_ice_shape, ONLY : t_shape_coeffs, getShapeCoefficients
  IMPLICIT NONE

  CONTAINS
  
    SUBROUTINE update_coagulation_kernels(kbdim,klev,ppres,ptemp,   &
                                          zccaa, zcccc, zccca, zccpc, zccpa,  &
                                          zccpp, zccia, zccic, zccii, zccip) 

      INTEGER, INTENT(in) :: kbdim,klev
      REAL, INTENT(in) :: ppres(kbdim,klev), ptemp(kbdim,klev)
      REAL, INTENT(out) :: zccaa(kbdim,klev,nbins,nbins), zcccc(kbdim,klev,ncld,ncld),  &
                           zccca(kbdim,klev,nbins,ncld), zccpc(kbdim,klev,ncld,nprc),   &
                           zccpa(kbdim,klev,nbins,nprc), zccpp(kbdim,klev,nprc,nprc),   &
                           zccia(kbdim,klev,nbins,nice), zccic(kbdim,klev,ncld,nice),   &
                           zccii(kbdim,klev,nice,nice), zccip(kbdim,klev,nprc,nice)

      ! Aero-aero
      zccaa(:,:,:,:) = 0.
      IF (lscgaa) &
           CALL buildKernelSelf( kbdim,klev,nbins,aero,ptemp,ppres,zccaa )

      ! Cloud-cloud 
      zcccc(:,:,:,:) = 0.
      IF (lscgcc) &
           CALL buildKernelSelf( kbdim,klev,ncld,cloud,ptemp,ppres,zcccc )
      
      ! Precp-precp
      zccpp(:,:,:,:) = 0.
      IF (lscgpp) &
           CALL buildKernelSelf( kbdim,klev,nprc,precp,ptemp,ppres,zccpp )

      ! ice-ice
      zccii(:,:,:,:) = 0.
      IF (lscgii) &
           CALL buildKernelSelf( kbdim,klev,nice,ice,ptemp,ppres,zccii )

      ! Aero-cloud
      zccca(:,:,:,:) = 0.
      IF (lscgca) &
           CALL buildKernel( kbdim,klev,nbins,aero,ncld,cloud,ptemp,ppres,zccca )

      ! Aero-precp
      zccpa(:,:,:,:) = 0.
      IF (lscgpa) &
           CALL buildKernel( kbdim,klev,nbins,aero,nprc,precp,ptemp,ppres,zccpa )

      ! Aero-ice
      zccia(:,:,:,:) = 0.
      IF (lscgia) &
           CALL buildKernel( kbdim,klev,nbins,aero,nice,ice,ptemp,ppres,zccia )

      ! Cloud-precp
      zccpc(:,:,:,:) = 0.
      IF (lscgpc) &
           CALL buildKernel( kbdim,klev,ncld,cloud,nprc,precp,ptemp,ppres,zccpc )

      ! Cloud-ice
      zccic(:,:,:,:) = 0.
      IF (lscgic) &
           CALL buildKernel( kbdim,klev,ncld,cloud,nice,ice,ptemp,ppres,zccic )

      ! Precp-ice
      zccip(:,:,:,:) = 0.
      IF (lscgip) &
           CALL buildKernel( kbdim,klev,nprc,precp,nice,ice,ptemp,ppres,zccip )
            
    END SUBROUTINE update_coagulation_kernels

    ! ----------------------
    SUBROUTINE buildKernelSelf( kbdim,klev,nb1,part1,ptemp,ppres,zcc )
      INTEGER, INTENT(in) :: kbdim,klev,nb1
      TYPE(Section), INTENT(inout) :: part1(kbdim,klev,nb1)
      REAL, INTENT(in) :: ptemp(kbdim,klev),ppres(kbdim,klev)
      REAL, INTENT(inout) :: zcc(kbdim,klev,nb1,nb1)

      REAL :: zdiam_mm, zmass_mm
      
      REAL :: zdiam(nb1), zmass(nb1)

      INTEGER :: mm,nn,ii,jj

      DO jj = 1,klev
         DO ii = 1,kbdim
            
            zdiam(:) = 0.
            zmass(:) = 0.
            
            ASSOCIATE ( pp1 => part1(ii,jj,1:nb1) )
              
              DO mm = 1,nb1
                 CALL pp1(mm)%updateDiameter(limit=.TRUE.,type="all")
                 CALL pp1(mm)%updateRhomean()
              END DO

              IF (pp1(1)%phase == 4) zdiam(1:nb1) = pp1(1:nb1)%dnsp ! for ice use nonspherical shape for the actual diameter
              
              DO mm = 1, nb1         ! smaller colliding particle                 
                 IF (pp1(mm)%numc < pp1(mm)%nlim) CYCLE
                 
                 DO nn = mm, nb1            ! larger colliding particle
                    IF (pp1(nn)%numc < pp1(nn)%nlim) CYCLE
                    
                    zcc(ii,jj,mm,nn) = coagc( pp1(mm),pp1(nn),            &
                                              ptemp(ii,jj),ppres(ii,jj),2 )
                    zcc(ii,jj,nn,mm) = zcc(ii,jj,mm,nn)
                 END DO
              END DO
              
            END ASSOCIATE
            
         END DO
      END DO
      
    END SUBROUTINE buildKernelSelf

    ! --------------------------

    SUBROUTINE buildKernel( kbdim,klev,nb1,part1,nb2,part2,ptemp,ppres,zcc )
      ! Always the "smaller" particle indices first
      INTEGER, INTENT(in) :: kbdim,klev
      INTEGER, INTENT(in) :: nb1, nb2
      TYPE(Section), INTENT(inout) :: part1(kbdim,klev,nb1), part2(kbdim,klev,nb2)
      REAL, INTENT(in)    :: ptemp(kbdim,klev),ppres(kbdim,klev)
      REAL, INTENT(inout)   :: zcc(kbdim,klev,nb1,nb2)
      
      REAL :: zdiam1(nb1), zdiam2(nb2), zmass1(nb1), zmass2(nb2)

      INTEGER :: mm,nn,ii,jj

      DO jj = 1,klev
         DO ii = 1,kbdim
            
            zdiam1(:) = 0.; zdiam2(:) = 0.
            zmass1(:) = 0.; zmass2(:) = 0.
            
            ASSOCIATE ( pp1 => part1(ii,jj,1:nb1), pp2 => part2(ii,jj,1:nb2) )

              ! Update wet diameters
              DO mm = 1,nb1
                 CALL pp1(mm)%updateDiameter(limit=.TRUE.,type="all")
                 CALL pp1(mm)%updateRhomean()
              END DO
              DO nn = 1,nb2
                 CALL pp2(nn)%updateDiameter(limit=.TRUE.,type="all")
                 CALL pp2(nn)%updateRhomean()
              END DO
              
              DO mm = 1,nb1
                 IF (pp1(mm)%numc < pp1(mm)%nlim) CYCLE
                 DO nn = 1,nb2
                    IF (pp2(nn)%numc < pp2(nn)%nlim) CYCLE
                    zcc(ii,jj,mm,nn) = coagc( pp1(mm),pp2(nn),            &
                                              ptemp(ii,jj),ppres(ii,jj),2 )
                 END DO
              END DO
              
            END ASSOCIATE

         END DO
      END DO

    END SUBROUTINE buildKernel


    ! ==========================================

    REAL FUNCTION coagc(pp1,pp2,temp,pres,kernel)

      USE mo_submctl, ONLY : pi, pi6, boltz, pstand, grav, rd 
      USE mo_particle_external_properties, ONLY : terminal_vel

      IMPLICIT NONE
      
      !-- Input variables ----------

      TYPE(Section), INTENT(in) :: pp1,pp2

      REAL, INTENT(IN) :: &
           temp,   &   ! ambient temperature [K]
           pres        ! ambient pressure [fxm]
      
      INTEGER, INTENT(in) :: kernel ! select the type of kernel: 1 - aerosol-aerosol coagulation (the original version)
      !                            2 - hydrometeor-aerosol or hydrometeor-hydrometeor coagulation
           
      !-- Local variables ----------
      REAL ::  &
           visc,     &   ! viscosity of air [kg/(m s)]
           vkin,     &   ! Kinematic viscosity of air [m2 s-1]
           zrhoa,    &   ! Density of air [kg m-3]
           mfp,      &   ! mean free path of air molecules [m]
           mdiam,    &   ! mean diameter of colliding particles [m]
           fmdist,   &   ! distance of flux matching [m]
           eddy_dis, &   ! Eddy dissipation time
           zecoll,   &   ! Collision efficiency for graviational collection
           zev,      &   !
           zea,      &
           zbrown,   &   ! Components for coagulation kernel; Brownian
           zbrconv,  &   !                                    Convective diffusion enhancement
           zgrav,    &   !                                    Gravitational collection
           ztshear,  &   ! turbulent shear
           zturbinert    ! turbulent inertia

      REAL :: diam1, diam2,       & ! Particle diameters; for ice this should be the effective max dimension
              dicesph1, dicesph2, & ! Spherical equivalent diameters for ice
              mass1, mass2,       & ! Masses of particles
              rhop1, rhop2,       & ! Particle densities; For ice this is the effective density (low for non-spherical)
              rhoiceb1, rhoiceb2   ! Bulk ice densities
              
      
      REAL, DIMENSION (2) :: &
           diam,   &   ! diameters of particles [m]
           mpart,  &   ! masses of particles [kg]
           knud,   &   ! particle knudsen number [1]
           beta,   &   ! Cunningham correction factor [1]
           zrhop,  &   ! Particle density [kg m-3]; For ice this will be the effective density (==low for non-spherical)
           dfpart, &   ! particle diffusion coefficient [m2/s]
           mtvel,  &   ! particle mean thermal velocity [m/s]
           termv,  &   ! Particle terminal velocity
           omega,  &   !
           tva,    &   ! temporary variable [m]
           flux        ! flux in continuum and free molec. regime [m/s]
      
      REAL ::  &
           schm(2), &   ! Schmidt nubmer
           reyn(2), &    ! Reynolds number
           stok             ! Stokes number
      INTEGER :: lrg,sml

      INTEGER :: ns ! number of species
      
      TYPE(t_shape_coeffs) :: shape1, shape2 ! Shape coefficients needed for ice
      
      zbrown = 0.
      zbrconv = 0.
      zgrav = 0.
      zev = 0.
      coagc = 0.
      
      !-------------------------------------------------------------------------------

      ns = spec%getNSpec(type="total") ! includes rime
      
      !-- 0) Initializing particle and ambient air variables --------------------
      diam1 = MERGE(pp1%dnsp, pp1%dwet, pp1%phase == 4)
      diam2 = MERGE(pp2%dnsp, pp2%dwet, pp2%phase == 4)
      IF (pp1%phase == 4) dicesph1 = pp1%dwet
      IF (pp2%phase == 4) dicesph2 = pp2%dwet

      mass1 = pp1%rhomean*pi6*pp1%dwet**3
      mass2 = pp2%rhomean*pi6*pp2%dwet**3

      ! If this is for self coagulation, put a minor offset on the particle diameters to account for
      ! the bin width
      IF ( ABS(pp1%dmid - pp2%dmid)/pp1%dmid < 1.e-2) THEN
         diam1 = 0.9*diam2
         mass1 = 0.7290*mass2 ! 0.7290 corresponds to 0.9*diameter ;;  NOTE this might cause some small error for ice!!
      END IF

      
      diam  = (/ diam1, diam2 /)       ! particle diameters [m]
      mpart = (/ mass1, mass2 /)       ! particle masses [kg]
      
      visc = (7.44523e-3*SQRT(temp**3))/(5093.*(temp+110.4)) ! viscosity of air [kg/(m s)]
      
      mfp = (1.656e-10*temp+1.828e-8)*pstand/pres ! mean free path of air [m]
      
      !-- 2) Slip correction factor for small particles -------------------------
      
      knud = 2.*mfp/diam                                    ! Knudsen number
      beta = 1.+knud*(1.142+0.558*exp(-0.999/knud))! Cunningham correction factor
      ! (Allen and Raabe, Aerosol Sci. Tech. 4, 269)
      
      !-- 3) Particle properties ------------------------------------------------
      
      dfpart = beta*boltz*temp/(3.*pi*visc*diam)  ! diffusion coefficient [m2/s]
      mtvel  = sqrt((8.*boltz*temp)/(pi*mpart))    ! mean thermal velocity [m/s]
      omega  = 8.*dfpart/(pi*mtvel)
      
      mdiam = 0.5*(diam(1)+diam(2))               ! mean diameter [m]
      
      !-- 4) Calculation of fluxes and flux matching ----------------------------
      
      flux(1) = 4.*pi*mdiam*( dfpart(1)+dfpart(2) )    ! flux in continuum regime [m3/s]
      flux(2) = pi*sqrt((mtvel(1)**2)+(mtvel(2)**2))*(mdiam**2) !  -"- in free molec. regime [m3/s]
      
      tva(1) = ((mdiam+omega(1))**3 - &              ! temporary variable [m]
           (mdiam**2+omega(1)**2)* &
           sqrt((mdiam**2+omega(1)**2)))/ &
           (3.*mdiam*omega(1)) - mdiam
      
      tva(2) = ((mdiam+omega(2))**3 - &              ! temporary variable [m]
           (mdiam**2+omega(2)**2)* &
           sqrt((mdiam**2+omega(2)**2)))/ &
           (3.*mdiam*omega(2)) - mdiam
      
      fmdist = sqrt(tva(1)**2+tva(2)**2)             ! flux matching distance [m]
      
      SELECT CASE(kernel)
      CASE(1)
         
         ! Aerosol-Aerosol coagulation - like the f version
         !-- 5) Coagulation coefficient [m3/s] -------------------------------------
         coagc = flux(1) / (mdiam/(mdiam+fmdist) + flux(1)/flux(2))
         
      CASE(2)
         
         ! Which particle is larger?
         sml = 1; lrg = 2
         IF (diam(1) >= diam(2)) THEN
            lrg = 1; sml = 2
         END IF

         IF (pp1%phase == 4) &
              CALL getShapeCoefficients( shape1, SUM(pp1%volc(1:ns-1)*spec%rhoic),     &
                                         pp1%volc(ns)*spec%rhori,                      &
                                         pp1%numc                                      )

         IF (pp2%phase == 4) &
              CALL getShapeCoefficients( shape2, SUM(pp2%volc(1:ns-1)*spec%rhoic),     &
                                         pp2%volc(ns)*spec%rhori,                      &
                                         pp2%numc                                      )
         
         
         zrhoa = pres/(rd*temp)       ! Density of air
         zrhop = mpart/(pi6*diam**3)  ! Density of particles; For ice this is the effective density using the non-spherical diameter
         vkin  = visc/zrhoa           ! Kinematic viscosity of air [m2 s-1]

         IF (pp1%phase < 4) THEN         
            termv(1) = terminal_vel(pp1%dwet,pp1%rhomean,zrhoa,visc,beta(1),pp1%phase)
         ELSE
            termv(1) = terminal_vel(pp1%dwet,pp1%rhomean,zrhoa,visc,beta(1),pp1%phase,shape1,pp1%dnsp)
         END IF

         IF (pp2%phase < 4) THEN
            termv(2) = terminal_vel(pp2%dwet,pp2%rhomean,zrhoa,visc,beta(2),pp2%phase)
         ELSE
            termv(2) = terminal_vel(pp2%dwet,pp2%rhomean,zrhoa,visc,beta(2),pp2%phase,shape2,pp2%dnsp)
         END IF
         
         ! Reynolds number
         reyn = diam*termv/vkin
         ! Schmidt number for the smaller particle
         schm = vkin/dfpart
         ! Stokes number
         stok = 2.*termv(sml)*ABS(termv(1) - termv(2))/( diam(lrg)*grav )
         
         !Brownian component
         zbrown = flux(1) / (mdiam/(mdiam+fmdist) + flux(1)/flux(2))
         
         ! Convective enhancement
         IF (reyn(lrg) <= 1.) THEN
            zbrconv = 0.45*zbrown*( (reyn(lrg)*schm(sml))**(1./3.) )
         ELSE IF (reyn(lrg) > 1.) THEN
            zbrconv = 0.45*zbrown*SQRT(reyn(lrg))*( schm(sml)**(1./3.) )
         END IF
         
         ! Turbulent Shear
         eddy_dis = 10.e-4 ! Values suggested by Sami - could be taken from the LES model?
         ztshear = SQRT(8.*pi*eddy_dis/(15.*vkin))*(0.5*(diam(1)+diam(2)))**3
         zturbinert = pi*eddy_dis**(0.75) /(grav* SQRT(SQRT( vkin )))  &
              *(0.5*(diam(1)+diam(2)))**2* ABS(termv(1)-termv(2))
         
         ! gravitational collection
         zea = stok**2/( stok + 0.5 )**2
         IF (stok > 1.214) THEN
            zev = 0.75*LOG(2.*stok)/(stok - 1.214)
            zev = (1. + zev)**(-2)
         ELSE IF (stok <= 1.214) THEN
            zev = 0.
         END IF
         
         zecoll = (60.*zev + zea*reyn(lrg))/(60. + reyn(lrg))
         zgrav = zecoll * pi * mdiam**2
         zgrav = zgrav * ABS(termv(1)-termv(2))
         
         ! Total coagulation kernel
         coagc = zbrown  + zbrconv + SQRT(zgrav**2 + ztshear**2 + zturbinert**2)
         
      END SELECT
      
    END FUNCTION coagc
    


 
END MODULE mo_salsa_coagulation_kernels
