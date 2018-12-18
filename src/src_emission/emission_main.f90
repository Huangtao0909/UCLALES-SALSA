MODULE emission_main
  USE mo_seasalt_emission

  USE mo_submctl, ONLY : pi6, in1a, fn2a, in2b, fn2b, nbins, aerobins, spec, pi6

  USE mo_salsa_types, ONLY : aero

  USE mo_salsa_sizedist, ONLY : size_distribution  ! Could this be packaged somehow differently?

  USE grid, ONLY: deltax, deltay, dzt, zt,  & ! Note dzt is inverse of the distance
  ! Ali, addition of emission type 3     
                  xt, yt, deltaz, dtl, &                  
  ! Ali, addition of emission type 3     
                  nxp,nyp,nzp,              &
                  a_up, a_vp, a_dn,      &
                  a_maerot, a_naerot
  USE util, ONLY: smaller, closest, getMassIndex
  USE exceptionHandling, ONLY: errorMessage
  ! Ali, addition of emission type 3  
  USE mpi_interface, ONLY : myid
  ! Ali, addition of emission type 3
  
  IMPLICIT NONE

  CHARACTER(len=50), PARAMETER :: global_name = "emission_main"
 
  !
  ! ----------------------------------------------
  ! Contains arrays for emitted size distribution
  !
  TYPE EmitSizeDist
     REAL, ALLOCATABLE :: numc(:)
     REAL, ALLOCATABLE :: mass(:)
  END TYPE EmitSizeDist

  !
  ! --------------------------------------------
  ! Type for emission configuration options
  !
  TYPE  EmitConfig
     INTEGER          :: emitType = 1                 ! 1: Natural seasalt emissions, 2: custom artificial emissions
                                                      ! Ali, addition of emission type 3
                                                      ! 3: artificial emission given by a map (moving source of airborne emission)
                                                      ! Ali, addition of emission type 3
     INTEGER          :: regime = 1                   ! Destination bin regime for emitted aerosol. 1: A, 2: B
     REAL             :: start_time = 0.,  &          ! Start time for emission (s)
                         end_time = 86400.            ! End time for emission (s)

     ! Parameters below valid for emitType > 1
     CHARACTER(len=3) :: species = 'SS '              ! Which aerosol species to emit (must conform to names in src_salsa/classSpecies)
     REAL             :: emitHeightMin = -999.,  &    ! Min height of airborne emissions (m)
                         emitHeightMax = -999.        ! Max height (m)
     INTEGER          :: emitLevMin = -999            ! Integer levels corresponding to the heights; If both heights and lev indices are specified,
     INTEGER          :: emitLevMax = -999            ! the level indices are preferred (so use one or the other...).

     INTEGER          :: emitSizeDistType = 1         ! 1: Monochromatic aerosol, 2: modal size disribution (lognormal)
     REAL             :: emitDiam = 10.e-6,    &      ! Assumed (dry )diameter of the particles (mode diameter for emitType=2).
                         emitNum  = 10000.            ! Number consentration of particles emitted per second #/m3/s (mode concentration for emitType=2)
     REAL             :: emitSigma = 2.0              ! Geometric standard deviation for emitSizeDist=2
     ! Ali, addition of emission type 3
     CHARACTER(len=40):: emitMap = ''                 ! Name of the file providing all location of emission (only for emitType = 3)
     REAL             :: scS = -999.                  ! Source speed (m/s) (only for emitType = 3)
     INTEGER          :: z_expan_up = 0               ! Epands the emission map to adjacent cells above the given map
     INTEGER          :: z_expan_dw = 0               ! Epands the emission map to adjacent cells down the given map      
     ! Ali, addition of emission type 3
  END TYPE EmitConfig

  ! Ali, addition of emission type 3
   TYPE EmitType3Config
     INTEGER, ALLOCATABLE :: ix(:)  ! x Index of source location, calculated based on 'emitMap' and computational grid
     INTEGER, ALLOCATABLE :: iy(:)  ! y Index of source location, calculated based on 'emitMap' and computational grid
     INTEGER, ALLOCATABLE :: iz(:)  ! z Index of source location, calculated based on 'emitMap' and computational grid
     REAL, ALLOCATABLE :: t_in(:)   ! Enterance times of emission source in each subdomain/processor
     REAL, ALLOCATABLE :: t_out(:)  ! Exit times of emission source from each subdomain/processor
     REAL, ALLOCATABLE :: t(:)      ! Enterance times of emission source in each cell in each subdomain/processor
     REAL :: t_trac                 ! Tracks the time for the emission source in each subdomain/processor
     INTEGER :: np                  ! Number of points of emission trajectory intersecting with cell boundaries for each subdomain/processor
   END type EmitType3Config
  ! Ali, addition of emission type 3
   
  INTEGER :: nEmissionModes = 0                             ! Number of emission modes (NAMELIST variable, max=5)
  INTEGER, PARAMETER :: maxEmissionModes = 5                ! Max number of emission modes
  TYPE(EmitConfig), TARGET :: emitModes(MaxEmissionModes)   ! Configuration instances for emission modes
  TYPE(EmitSizeDist), TARGET :: emitData(MaxEmissionModes)  ! Emission size distribution data for each setup/mode.
  !Ali, addition of emission type 3
  TYPE(EmitType3Config), TARGET :: emitType3(MaxEmissionModes)   ! Configuration instances for emission type 3
  !Ali, addition of emission type 3

  !Ali, addition of emission type 3
  ! Dynamic allocator for real and integer 1D arrays
  INTERFACE arr_resize
    MODULE PROCEDURE arr_resize_int
    MODULE PROCEDURE arr_resize_rel
  END INTERFACE arr_resize
  !Ali, addition of emission type 3
  
  ! The limit of 5 emission modes can be easily increased if needed. The maximum has to be hardcoded here
  ! as a PARAMETER in order to be able to give the configuration options directly from NAMELIST (the instances must be allocated).  

  CONTAINS

  SUBROUTINE init_emission
    IMPLICIT NONE
    CHARACTER(len=50), PARAMETER :: name = "init_emission"

    INTEGER :: ibin
    INTEGER :: ilev
    INTEGER :: st,en
    INTEGER :: maxv
    INTEGER :: nprof
    INTEGER :: mb1,mb12,mb2,mb22,nc1,nc2
    REAL :: core(nbins), naero(1,1,nbins)
    !Ali, addition of emission type 3
    TYPE(EmitType3Config), POINTER :: emdT3 => NULL()
    REAl, ALLOCATABLE ::  x(:), y(:), z(:)
    !Ali, addition of emission type 3

    DO nprof = 1,nEmissionModes

      ASSOCIATE(emd => emitModes(nprof), edt => emitData(nprof))
      
       ! Use pointers for individual emission configuration and data instances to clean up the code
       IF (myid == 0) THEN
         WRITE(*,*) ''
         WRITE(*,*) '------------------------------------'
         WRITE(*,*) TRIM(name)//': INITIALIZING EMISSION PROFILE NO ',nprof
       END IF
       ! The aerosol species specified for emission MUST also be active in SALSA configuration
       IF ( .NOT. spec%isUsed(emd%species) ) THEN
         CALL errorMessage(global_name, name, &
             'Attempt to emit <'//TRIM(emd%species)//'> but the compound '// &
             'is not set to be used in the SALSA namelist.')
         STOP
       END IF
         
       ! Initialize the emission size distribution arrays
       ALLOCATE( edt%numc(nbins), edt%mass(nbins*spec%getNSpec()) )
       edt%numc(:) = 0.; edt%mass(:) = 0.
         
       IF (emd%emitType == 1) THEN       ! Parameterized sea salt emissions
         IF (myid == 0) THEN   
           WRITE(*,*) 'SEA SALT EMISSIONS; NO INITIALIZATION IMPLEMENTED (NOR NEEDED?)'
         END IF
         
         edt%numc(:) = 0.; edt%mass(:) = 0. ! i.e. do nothing
            
       !ELSE IF (emd%emitType == 2) THEN  ! Artificial emission
       ! Ali, addition of emission type 3
       ELSE IF (emd%emitType >= 2) THEN ! Artificial emission, 2:two height limits, 3:a trajectory (e.g. flight map)
       ! Ali, addition of emission type 3
         
         ! Index limits for the regimes
         CALL regime_limits(emd%regime,st,en)
            
         IF (emd%emitSizeDistType == 1 ) THEN ! Monochromatic aerosol
               
               ! Determine the destination bin of the monochromatic emission from the
               ! specified diameter
               ibin = smaller(aerobins(st:en),emd%emitDiam) - 1
               
               ! Get bin indices for emitted mass tendencies (in addition to specified species, put also
               ! some water due to current limitations in condensation in mo_salsa_dynamics. Review this later!)
               CALL bin_indices(nprof,nbins,st,en,   &
                                nc1, nc2, mb1, mb2,  &
                                mb12, mb22, ibin=ibin)
               
               ! Place the emission number and mass concentrations to the emission data instance
               edt%numc(st+ibin) = emd%emitNum
               edt%mass(mb1) = emd%emitNum *   &
                    (pi6*emd%emitDiam**3)*spec%rholiq(nc1)
               ! Small amount of water due the condensation issue
               edt%mass(mb2) = 0.001 * edt%mass(mb1)
               
            ELSE IF (emd%emitSizeDistType == 2 ) THEN! Lognormal mode
               ! Consider specified emitNum and emitDiam as lognormal mode values. Also emitSigma needs to be specified here
               
               ! First, get the core volume of single particle in each aerosol bin based on the bin mean diameter
               core(:) = pi6 * aero(1,1,:)%dmid**3
               
               ! Get the binned size distribution (brackets used because size_distribution expects arrays for mode values)
               CALL size_distribution(1,1,1,1,[emd%emitNum],    &
                                              [emd%emitDiam],   &
                                              [emd%emitSigma],  &
                                               naero)
               !naero = 1.

               ! -- Get the indices for emitted species and water.
               CALL bin_indices(nprof,nbins,st,en,  &
                                nc1, nc2, mb1, mb2, &
                                mb12, mb22          )
               
               ! Set the emission number and mass concentrations to the emission data instance
               edt%numc(st:en) = naero(1,1,st:en)
               edt%mass(mb1:mb12) = naero(1,1,st:en)*core(st:en)*spec%rholiq(nc1)
               ! Small amount of water
               edt%mass(mb2:mb22) = 0.001 * edt%mass(mb1:mb12)
               
            END IF
            
          END IF
       
         ! Set up the emission levels. This will preferentially use the level indices, except if they're not given by the namelist
         !CALL init_emission_heights(emd)

         ! Ali, addition of emission type 3
          IF (emd%emitType == 2)  THEN
            CALL init_emission_heights(emd)
          ELSE IF (emd%emitType == 3)  THEN
            emdT3 => emitType3(nprof)
            CALL init_emitType3_map(emd%emitMap,x,y,z,emdT3%np)
            CALL lagrangian_tracker(emdT3%ix,emdT3%iy,emdT3%iz,emdT3%t,emdT3%np, &
                                 emdT3%t_in,emdT3%t_out,emdT3%t_trac,&
                                 emd%scS,emd%start_time,emd%end_time,x,y,z)
          END IF
          
          emdT3 => NULL()
          ! Ali, addition of emission type    
       
        END ASSOCIATE
          
        END DO

    IF (myid == 0) THEN
    WRITE(*,*) '-------------------------------------'
    WRITE(*,*) ''
    END IF
  
  END SUBROUTINE init_emission

  !
  ! ----------------------------------------------------------------------------
  ! Subroutine init_emission_heights: Sorts out the level and height information
  !                                   in the emission configuration
  !
  SUBROUTINE init_emission_heights(emd)
    IMPLICIT NONE
    CHARACTER(len=50), PARAMETER :: name = "init_emission_heights"

    TYPE(EmitConfig), INTENT(inout) :: emd
    
    INTEGER :: ilev
    INTEGER :: maxi
    REAL    :: maxf


    ASSOCIATE( emitLevMin => emd%emitLevMin,       &
               emitLevMax => emd%emitLevMax,       &
               emitHeightMin => emd%emitHeightMin, &
               emitHeightMax => emd%emitHeightMax  )
    
      ! At least one of the height or level definitions in the emitConfig must be specified
      IF ( ALL( [emitLevMin, emitLevMax] == -999 ) .AND.  &
           ALL( [emitHeightMin, emitHeightMax] == -999. ) ) THEN
      
          CALL errorMessage(global_name, name, &
               "At least of the the height or level definitions must be"// &
               " specified in emitConfig")
          STOP

      END IF

      IF ( emitLevMin == -999 .AND. emitLevMax == -999) THEN 

         ! EmitConfig must have emitHeightMin or emitHeightMax specified as a positive value
         maxf = MAX(emitHeightMin, emitHeightMax)
         IF (emitHeightMin == -999.) emitHeightMin = maxf
         IF (emitHeightMax == -999.) emitHeightMax = maxf

         ilev = closest(zt,emitHeightMin)
         emitLevMin = ilev
         ilev = closest(zt,emitHeightMax)
         emitLevMax = ilev
         
      ELSE IF (emitLevMin > 0 .OR. emitLevMax > 0) THEN
         
         maxi = MAX(emitLevMin, emitLevMax)
         IF (emitLevMin == -999) emitLevMin = maxi
         IF (emitLevMax == -999) emitLevMax = maxi
         
         ! Update the emission height levels according to the indices
         emitHeightMin = zt(emitLevMin)
         emitHeightMax = zt(emitLevMax)
         
      END IF
      
    END ASSOCIATE

  END SUBROUTINE init_emission_heights

  !
  ! --------------------------------------------------------------------------
  ! Subroutine regime_limits: Returns the limiting indices for A or B regimes
  !                           (reg=1 or reg=2, respectively)
  !
  SUBROUTINE regime_limits(reg,st,en)
    IMPLICIT NONE
    CHARACTER(len=50), PARAMETER :: name = "regime_limits"

    INTEGER, INTENT(in) :: reg
    INTEGER, INTENT(out) :: st,en

    IF ( reg == 1 ) THEN
       st = in1a
       en = fn2a
    ELSE IF ( reg == 2 ) THEN
       st = in2b
       en = fn2b
    END IF

  END SUBROUTINE regime_limits

  !
  ! -----------------------------------------------------------------------------
  ! Subroutine bin_indices: Returns bin and mass indices needed in init_emission
  !
  SUBROUTINE bin_indices(nprof,nbins,st,en,          &
                         nc_emit, nc_h2o,      &
                         mb_emit_1, mb_h2o_1,  &
                         mb_emit_2, mb_h2o_2, ibin)
    IMPLICIT NONE
    CHARACTER(len=50), PARAMETER :: name = "bin_indices"

    INTEGER, INTENT(in) :: nprof,nbins, st, en
    INTEGER, INTENT(in), OPTIONAL :: ibin

    INTEGER, INTENT(out) :: nc_emit, nc_h2o,      &  ! Mass indices for emitted species and water
                            mb_emit_1, mb_h2o_1,  &  ! First index of the regime for emission size distribution 
                                                     ! (if ibin present, this will be the single emission bin for monochromatic emissions)
                            mb_emit_2, mb_h2o_2      ! Last index of the regime for emission size distribution
                                                     ! (will be put out, but not used for monochromatic emissions)

    nc_emit  = spec%getIndex(emitModes(nprof)%species)
    nc_h2o  = spec%getIndex("H2O")
    
    IF (myid == 0) THEN
      WRITE(*,*) name,'---------------------------------'
      WRITE(*,*) nc_emit, nc_h2o
      WRITE(*,*) '--------------------------------------'
    END IF
 
    IF ( PRESENT(ibin) ) THEN
       mb_emit_1 = getMassIndex(nbins,st+ibin,nc_emit)
       mb_h2o_1 = getMassIndex(nbins,st+ibin,nc_h2o)       
    ELSE
       mb_emit_1  = getMassIndex(nbins,st,nc_emit)
       mb_h2o_1  = getMassIndex(nbins,st,nc_h2o)
    END IF

    mb_emit_2 = getMassIndex(nbins,en,nc_emit)
    mb_h2o_2 = getMassIndex(nbins,en,nc_h2o)

  END SUBROUTINE bin_indices

  !
  ! -------------------------------------------------------------------
  ! subroutine aerosol_emission:  calls methods to calculate emitted
  !                               aerosols from ground/sea
  !  
  ! Adapted from the original code by Antti Kukkurainen
  ! Juha Tonttila, FMI, 2017
  !
  SUBROUTINE aerosol_emission(time_in)
    IMPLICIT NONE
    CHARACTER(len=50), PARAMETER :: name = "aerosol_emission"

    REAL, INTENT(in) :: time_in   ! time in seconds
    LOGICAL :: condition
    INTEGER :: pr
    !Ali, addition of emitType 3  
    TYPE(EmitType3Config), POINTER :: emdT3
    INTEGER :: conditionT3 = 0

    emdT3 => NULL()
    !Ali, addition of emitType 3  
    
    ! Loop over all specified emission profiles
    DO pr = 1,nEmissionModes
      ASSOCIATE(emd => emitModes(pr), edt => emitData(pr))
      condition = getCondition(emd,time_in)
      IF (condition .AND. emd%emitType == 2) CALL custom_emission(edt,emd)

      !Ali, addition of emitType 3 
      IF (emd%emitType == 3) THEN
        emdT3 => emitType3(pr)
        conditionT3 = getConditionT3(emdT3,time_in)
        IF (conditionT3 > 0) THEN
          CALL custom_emission_typ3(edt,emd,emdT3,time_in,conditionT3)
        END IF
      END IF
      !Ali, addition of emitType 3
         
    END ASSOCIATE
    ! Ali, addition of emission type
    emdT3 => NULL()
    ! Ali, addition of emission type
      
    END DO

  END SUBROUTINE aerosol_emission

  !
  ! ---------------------------------------------------------------
  ! Simulates the emission of seasalt particles from
  ! an ocean surface as a function of the 10-m wind
  ! speed.
  !
  !SUBROUTINE surface_emission()
  !  IMPLICIT NONE

 !   CHARACTER(len=50), PARAMETER :: name = "surface_emission"
 !   
 !   REAL :: mass_flux(1,nbins) !mass flux at given radius
 !   REAL :: numb_flux(1,nbins)        !number flux at given radius
 !   REAL :: pseaice(1) = 0            !sea ice fraction
 !   REAL :: velo10m_salsa(1,1)        !wind speed

!    INTEGER :: nc, st, en, ii, jj
!    INTEGER :: in, fn
    
!    ! Surface seasalt emissions possible only if sea salt aerosol is used
!    IF (spec%isUsed(spec%nss)) THEN
!       nc=spec%getIndex(spec%nss)
       
!       IF (esrfc%regime == 1) THEN
!          in = in1a
!          fn = fn2a
!       ELSE IF (esrfc%regime == 2) THEN
!          in = in2b
!          fn = fn2b
!       END IF
!       st=(nc-1)*nbins+in
!       en=(nc-1)*nbins+fn
!       
!       DO jj=3,nyp-2
!          DO ii=3,nxp-2
!             
!             velo10m_salsa = SQRT(a_up(2,ii,jj)**2 + a_vp(2,ii,jj)**2)
!             
!             CALL seasalt_emissions_lsce_salsa(1, 1, 1, pseaice, velo10m_salsa, mass_flux, numb_flux)
!             
!             !number of particles + more particles per unit of time * scaling factor [#/kg]
!             a_naerot(2,ii,jj,in:fn) = a_naerot(2,ii,jj,in:fn) + numb_flux(1,in:fn)*(dzt(2)/a_dn(2,ii,jj)) 
!             !mass + more mass per unit of time * scaling factor [kg/kg]
!             a_maerot(2,ii,jj,st:en) = a_maerot(2,ii,jj,st:en) + mass_flux(1,in:fn)*(dzt(2)/a_dn(2,ii,jj)) 
!             
!          END DO
!       END DO
!       
!    END IF
!    
!  END SUBROUTINE surface_emission

  !
  ! ----------------------------------------------------------------
  ! Subroutine custom_emissio: "Customized" emission routine, mainly
  !                            for simulating atnhropogenic emissions,
  !                            such as ship or aircraft emissions etc.
  !                            Support for point sources will be included
  !                            soon. Now only does domain-wide emissions
  !                            at specified altitude and time.
  !
  SUBROUTINE custom_emission(edt,emd)
    IMPLICIT NONE
    
    CHARACTER(len=50), PARAMETER :: name = "cloud_seeding"
    
    TYPE(EmitSizeDist), INTENT(in) :: edt  ! Emission data instance
    TYPE(EmitConfig), INTENT(in) :: emd   ! Emission configuration instance

    INTEGER :: k,i,j,bb,ss,mm
    CHARACTER(len=30) :: emit_spec

    
    IF (myid == 0) THEN
      WRITE(*,*) '========================'
      WRITE(*,*) 'CALCULATING EMISSIONS'
      WRITE(*,*) '========================'
    END IF

    ASSOCIATE( k1 => emd%emitLevMin, k2 => emd%emitLevMax )
    
      DO bb = 1,nbins
         DO j = 1,nyp
            DO i = 1,nxp
               a_naerot(k1:k2,i,j,bb) = a_naerot(k1:k2,i,j,bb) + edt%numc(bb)
               DO ss = 1,spec%getNSpec()
                  mm = getMassIndex(nbins,bb,ss)
                  a_maerot(k1:k2,i,j,mm) = a_maerot(k1:k2,i,j,mm) + edt%mass(mm)
               END DO
            END DO
         END DO
      END DO
      
    END ASSOCIATE

  END SUBROUTINE custom_emission
  
  ! ----------------------------------------------------------

  FUNCTION getCondition(emd,time)
    IMPLICIT NONE
    LOGICAL :: getCondition
    CLASS(EmitConfig), INTENT(in) :: emd
    REAL, INTENT(in) :: time
    CHARACTER(len=50), PARAMETER :: name = "getCondition"

    getCondition = (                                &
                    emd%start_time <= time    .AND. &
                    emd%end_time > time            &
                   )   

  END FUNCTION getCondition
  
  !Ali, addition of emission type 3
  ! -----------------------------------------------------------------------------
  ! Subroutine init_emitType3_map: reads emitMap and calls lagrangian_tracker to return 
  ! cell indices and residence times of the emission source in the cells 
  SUBROUTINE init_emitType3_map(emitMap,x,y,z,np)
    CHARACTER(len=50), PARAMETER :: name = "init_emitType3_map"
    
    CHARACTER(len=40), INTENT(in) :: emitMap
    REAl, ALLOCATABLE, INTENT(out) ::  x(:), y(:), z(:)
    INTEGER, INTENT(out) :: np 

    REAl    :: deltax2, deltay2, xlim1, xlim2, ylim1, ylim2
    INTEGER :: i
    REAl    :: eps = 2e-2


    IF (myid == 0) THEN
      WRITE(*,*) ''
      WRITE(*,*) '------------------------------------'
      WRITE(*,*) TRIM(name)//': READING EMISSION MAP FROM ', TRIM(emitMap)
    END IF
    
    np = 8
    OPEN(1001,file=emitMap,status='old',form='formatted')

    IF (ALLOCATED(x)) DEALLOCATE(x)
    IF (ALLOCATED(y)) DEALLOCATE(y)
    IF (ALLOCATED(z)) DEALLOCATE(z)
    
    ALLOCATE (x(np), y(np), z(np))

    
    deltax2 = deltax/2
    deltay2 = deltay/2
  
    xlim1 = MINVAL(xt) + 3*deltax2
    xlim2 = MAXVAL(xt) - 3*deltax2
    ylim1 = MINVAL(yt) + 3*deltay2
    ylim2 = MAXVAL(yt) - 3*deltay2


    i = 1
    DO WHILE (.TRUE.)
      READ(1001,*,end=101) x(i), y(i), z(i)
      IF (x(i) == xlim1) x(i) = xlim1 - eps
      IF (x(i) == xlim2) x(i) = xlim2 - eps
      IF (y(i) == ylim1) y(i) = ylim1 + eps
      IF (y(i) == ylim2) y(i) = ylim2 + eps
      
      IF ( i >= np) THEN
        np = np*2
        CALL arr_resize(x,np)
        CALL arr_resize(y,np)
        CALL arr_resize(z,np)
      END IF
      i = i + 1
    END DO
    
101 CONTINUE
    CLOSE(1001)
    
    np = i - 1
    CALL arr_resize(x,np)
    CALL arr_resize(y,np)
    CALL arr_resize(z,np)
   
END SUBROUTINE init_emitType3_map
  ! -----------------------------------------------------------------------------
  ! Subroutine init_emitType3_map: reads emitMap and calls lagrangian_tracker to return 
SUBROUTINE lagrangian_tracker(ix,iy,iz,t,np,t_in,t_out,t_trac,scS,start_time,end_time,x,y,z)
  
  CHARACTER(len=50), PARAMETER :: name = "lagrangian_tracker"
  
  REAl, INTENT(in) :: scS
  REAl, INTENT(in) :: start_time  
  REAl, INTENT(inout) :: end_time
  REAl, INTENT(inout) :: t_trac
  REAl, INTENT(in) :: x(:), y(:), z(:)
  INTEGER, INTENT(inout) :: np

  REAl, ALLOCATABLE, INTENT(out) ::  t(:),t_in(:),t_out(:)
  INTEGER, ALLOCATABLE, INTENT(out) ::  ix(:), iy(:), iz(:)

  INTEGER :: i, j, k, io, ix1, iy1, iz1
  INTEGER :: np_trac   
  REAl    :: deltax2, deltay2, deltaz2, d, d_x, d_y, d_z, dx, dy, dz, dx2x1, dy2y1, dz2z1, &
             x_trac, y_trac, z_trac, vx, vy, vz, xlim1, xlim2, ylim1, ylim2, x1, x2, y1, y2, z1, z2  
  REAl    :: t_trac_old, sub_dt, dt, dtx, dty, dtz, min_t, dt2(2) = -999., dt4(4) = -999.
  REAl    :: eps = 1e-6
  LOGICAL :: nexpout = .FALSE.

  deltax2 = deltax/2
  deltay2 = deltay/2
  deltaz2 = deltaz/2
  np_trac = 2
  ALLOCATE(ix(np_trac),iy(np_trac),iz(np_trac),t(np_trac))
  ALLOCATE(t_in(1),t_out(1))
  ix = -999
  iy = -999
  iz = -999
  t  = -999.
  t_in = -999.
  t_out = t_in - 999.
  
  xlim1 = MINVAL(xt) + 3*deltax2
  xlim2 = MAXVAL(xt) - 3*deltax2
  ylim1 = MINVAL(yt) + 3*deltay2
  ylim2 = MAXVAL(yt) - 3*deltay2


  x_trac = xlim1 - 999.
  y_trac = ylim1 - 999.
  z_trac = -999.

  j  = 1
  io = 0
  DO i = 1,np-1
    t_trac = start_time
    d_x = x(i+1) - x(i)
    d_y = y(i+1) - y(i)
    d_z = z(i+1) - z(i)
    
    IF ( (d_x > eps) .AND. (ABS(d_y) > eps) ) THEN
      dz2z1 = 0.
      IF (d_y >  eps) THEN
        x1 = (ylim1-y(i))*(d_x/d_y) + x(i)
        x1 = MAX(x1,MAX(x(i),xlim1))
        y1 = (x1-x(i))*(d_y/d_x) + y(i)
        x2 = (ylim2-y(i))*(d_x/d_y) + x(i)
      END IF
      IF (d_y < -eps) THEN
        x1 = (ylim2-y(i))*(d_x/d_y) + x(i)
        x1 = MAX(x1,MAX(x(i),xlim1))
        y1 = (x1-x(i))*(d_y/d_x) + y(i)
        x2 = (ylim1-y(i))*(d_x/d_y) + x(i)
      END IF
      x2 = MIN(x2,MIN(x(i+1),xlim2))
      y2 = (x2-x(i))*(d_y/d_x) + y(i)
      dx2x1 = x2 - x1
      dy2y1 = y2 - y1
      
    ELSE IF ( (d_x < -eps) .AND. (ABS(d_y) > eps) )  THEN
      dz2z1 = 0.
      IF (d_y >  eps) THEN
        x1 = (ylim1-y(i))*(d_x/d_y) + x(i)
        x1 = MIN(x1,MIN(x(i),xlim2))
        y1 = (x1-x(i))*(d_y/d_x) + y(i)
        y1 = (x1-x(i))*(d_y/d_x) + y(i)
        x2 = (ylim2-y(i))*(d_x/d_y) + x(i)
      END IF
      IF (d_y < -eps) THEN
        x1 = (ylim2-y(i))*(d_x/d_y) + x(i)
        x1 = MIN(x1,MIN(x(i),xlim2))
        y1 = (x1-x(i))*(d_y/d_x) + y(i)
        x2 = (ylim1-y(i))*(d_x/d_y) + x(i)
        END IF
      x2 = MAX(x2,MAX(xlim1,x(i+1)))
      y2 = (x2-x(i))*(d_y/d_x) + y(i)
      dx2x1 = x2 - x1
      dy2y1 = y2 - y1
      
    ELSE IF ( (ABS(d_x) > eps) .AND. (ABS(d_y) <= eps) )  THEN
      dz2z1 = 0.
      IF (d_x >  eps) THEN
        x1 = MAX(x(i),xlim1)
        x2 = MIN(x(i+1),xlim2)
        dx2x1 = x2 - x1
      END IF
            IF (d_x < -eps) THEN
        x1 = MIN(x(i),xlim2)
        x2 = MAX(x(i+1),xlim1)
        dx2x1 = x2 - x1
      END IF
      IF ((y(i) > ylim1).AND.(y(i) < ylim2)) THEN
        dy2y1 = 0.
        d_y = 0.
        y1 = y(i)
        y2 = y(i)
      ELSE
        dy2y1 = -1.
        d_y = 1.
      END IF
      
    ELSE IF ( (ABS(d_x) <= eps) .AND. (ABS(d_y) >  eps) )  THEN
      dz2z1 = 0.
      IF (d_y >  eps) THEN
        y1 = MAX(y(i),ylim1)
        y2 = MIN(y(i+1),ylim2)
        dy2y1 = y2 - y1 
      END IF
      IF (d_y < -eps) THEN
        y1 = MIN(y(i),ylim2)
        y2 = MAX(y(i+1),ylim1)
        dy2y1 = y2 - y1 
      END IF
      IF ((x(i) > xlim1).AND.(x(i) < xlim2)) THEN
        dx2x1 = 0.
        d_x = 0.
        x1 = x(i)
        x2 = x(i)
      ELSE
        dx2x1 = -1.
        d_x = 1.
      END IF
      
    ELSE IF ( (ABS(d_x) <= eps) .AND. (ABS(d_y) <= eps) )  THEN
      z1 = z(i)
      z2 = z(i+1)
      dz2z1 = z2 - z1
      IF ( ((x(i) > xlim1).AND.(x(i) < xlim2)).AND.((y(i) > ylim1).AND.(y(i) < ylim2)) ) THEN
        dx2x1 = 0.
        d_x = 0.
        x1 = x(i)
        x2 = x(i)
        dy2y1 = 0.
        d_y = 0.
        y1 = y(i)
        y2 = y(i)
      ELSE
        dx2x1 = -1.
        d_x = 1.
        dy2y1 = -1.
        d_y = 1.
      END IF
      
    END IF

    IF ((d_x*dx2x1 >= 0.).AND.(d_y*dy2y1 >= 0.).AND.(ANY([abs(dx2x1),abs(dy2y1),abs(dz2z1)] > eps))) THEN
    
      IF ((ABS(d_x) > eps) .AND. (ABS(d_y) < eps)) THEN
        z1 = (d_z/d_x) * (x1-x(i)) + z(i)
        z2 = (d_z/d_x) * (x2-x(i)) + z(i)
        dz2z1 = z2 - z1 
      ELSE IF ((ABS(d_x) < eps) .AND. (ABS(d_y) > eps)) THEN
        z1 = (d_z/d_y) * (y1-y(i)) + z(i)
        z2 = (d_z/d_y) * (y2-y(i)) + z(i)
        dz2z1 = z2 - z1
      ELSE IF ((ABS(d_x) > eps) .AND. (ABS(d_y) > eps)) THEN
        z1 = (d_z/d_y) * (y1-y(i)) + z(i)
        z2 = (d_z/d_y) * (y2-y(i)) + z(i)
        dz2z1 = z2 - z1
      END IF

      IF ((j == 1) .OR. ( ( x1 < (x_trac - eps) ).OR.( x1 > (x_trac + eps) ) )) THEN  
        ix1 = closest(xt,x1)
        ix(j) = ix1
      END IF

      IF ((j == 1) .OR. ( (y1 < (y_trac - eps) ).OR.( y1 > (y_trac + eps) ) )) THEN    
        iy1 = closest(yt,y1)
        iy(j) = iy1
      END IF
      
      IF ((j == 1) .OR. ( (z1 < (z_trac - eps) ).OR.( z1 > (z_trac + eps) ) )) THEN  
        iz1 = closest(zt,z1)
        iz(j) = iz1
      END IF
      
      nexpout = ((x(i+1) < xlim1).OR.(x(i+1) > xlim2).OR.(y(i+1) < ylim1).OR.(y(i+1) > ylim2))
      
      IF (i == 1) THEN
        t_trac = t_trac + sqrt((x1-x(i))**2 + (y1-y(i))**2 + (z1-z(i))**2)/scS
        t_trac = t_trac 
        t(j) = t_trac
        io = io + 1
        t_in(io)  = t_trac
      END IF
      
      IF (i > 1) THEN
        IF ((x_trac >= xlim1 + eps ).AND.(x_trac <= xlim2 - eps ) &
            .AND.(y_trac >= ylim1 + eps ).AND.(y_trac<= ylim2 - eps )) THEN
        
          t_trac = t_trac_old
        ELSE
          t_trac = t_trac + sqrt((x1-x(i))**2 + (y1-y(i))**2 + (z1-z(i))**2)/scS
          IF (i > 2) THEN
            DO k = 1,i-2
              d_x = x(k+1) - x(k)
              d_y = y(k+1) - y(k)
              d_z = z(k+1) - z(k)
              d  = sqrt(d_x**2 + d_y**2 + d_z**2)
              dt = d/scS
              vx = d_x/dt
              vy = d_y/dt
              vz = d_z/dt
              t_trac = t_trac + dt
              t(j)   = t_trac
            END DO
          END IF
          d_x = x(i) - x(i-1)
          d_y = y(i) - y(i-1)
          d_z = z(i) - z(i-1)
          d  = sqrt(d_x**2 + d_y**2 + d_z**2)
          dt = d/scS
          vx = d_x/dt
          vy = d_y/dt

          IF ((x(i) >= xlim1 + eps ).AND.(x(i) <= xlim2 - eps ) &
              .AND.(y(i) >= ylim1 + eps).AND.(y(i) <= ylim2 - eps)) THEN
        
            IF  (d_x > eps) THEN
              dx  = xt(ix(j)) - deltax2 - x(i-1) 
              dtx = dx/vx
            ELSE IF (d_x < -eps) THEN
              dx  = xt(ix(j)) + deltax2 - x(i-1) 
              dtx = dx/vx
            ELSE
              dtx = -999.
            END IF

            IF ( (d_y > eps) .AND. (y(i) >= y1) ) THEN
              dy  = yt(iy(j)) - deltay2 - y(i-1) 
              dty = dy/vy
            ELSE IF ( (d_y < -eps) .AND. (y(i) <= y1) ) THEN
              dy  = yt(iy(j)) + deltay2 - y(i-1) 
              dty = dy/vy
            ELSE 
              dty = -999.
            END IF
     
            IF ((ABS(dx) > eps) .OR. (ABS(dy) > eps)) THEN
              dt2 = [dtx,dty]
              min_t = MAXVAL(dt2, MASK = dt2 > 0.) 
              t_trac = t_trac + min_t
              t(j)   = t_trac
            ELSE
              t(j)   = t_trac
            END IF
            
          ELSE
            t_trac = t_trac + dt
            t(j)   = t_trac
          END IF
          
          IF (nexpout.OR.(x(i) < xlim1).OR.(x(i) > xlim2).OR.(y(i) < ylim1).OR.(y(i) > ylim2)) THEN
            io = io + 1
            CALL arr_resize(t_in,io)
            CALL arr_resize(t_out,io)
            t_in(io)  = t_trac
          END IF
        
        END IF
        
      END IF
     
      d  = sqrt(dx2x1**2 + dy2y1**2 + dz2z1**2)
      dt = d/scS
      vx = dx2x1/dt
      vy = dy2y1/dt
      vz = dz2z1/dt
      sub_dt = 0.
      x_trac = x1
      y_trac = y1
      z_trac = z1
      end_time = t_trac + dt
      
      DO WHILE (sub_dt < dt)
      
        IF (dx2x1 > eps)  THEN
          dx  = xt(ix(j)) + deltax2 - x_trac
          dtx = dx/vx
        ELSE IF (dx2x1 < -eps) THEN
          dx  = xt(ix(j)) - deltax2 - x_trac
          dtx = dx/vx
        ELSE
          dtx = -999.
          dx = 0.
          ix(j+1) = ix(j)
        END IF
      
        IF (dy2y1 > eps)  THEN
          dy  = yt(iy(j)) + deltay2 - y_trac
          dty = dy/vy
        ELSE IF (dy2y1 < -eps) THEN
          dy  = yt(iy(j)) - deltay2 - y_trac
          dty = dy/vy
        ELSE
          dty = -999.
          dy = 0. 
          iy(j+1) = iy(j)
        END IF

        IF (dz2z1 > eps)  THEN
          dz  = zt(iz(j)) + deltaz2 - z_trac
          dtz = dz/vz
        ELSE IF (dz2z1 < -eps) THEN
          dz  = zt(iz(j)) - deltaz2 - z_trac
          dtz = dz/vz
        ELSE
          dtz = -999.
          dz = 0. 
          iz(j+1) = iz(j) 
        END IF
      
        IF (d > eps)  THEN
          dt4 = [dtx,dty,dtz,(dt-sub_dt)]
          min_t = MINVAL(dt4, MASK = dt4 >= 0.)
          
          IF (min_t == dtx) THEN
            IF (vx > 0)  ix(j+1) = ix(j) + 1
            IF (vx < 0)  ix(j+1) = ix(j) - 1
            IF (ABS(dtx-dty) < eps) THEN
              IF(vy > 0) iy(j+1) = iy(j) + 1
              IF(vy < 0) iy(j+1) = iy(j) - 1
            ELSE
              iy(j+1) = iy(j)
            END IF
              
            iz(j+1) = iz(j)
            t_trac = t_trac + min_t
            t(j+1) = t_trac
            t_trac_old = t_trac

          ELSE IF (min_t == dty) THEN
            IF (ABS(dtx-dty) < eps) THEN
              IF(vx > 0) ix(j+1) = ix(j) + 1
              IF(vx < 0) ix(j+1) = ix(j) - 1
            ELSE
              ix(j+1) = ix(j)
            END IF
            IF (vy > 0)  iy(j+1) = iy(j) + 1
            IF (vy < 0)  iy(j+1) = iy(j) - 1
            iz(j+1) = iz(j)
            t_trac = t_trac + min_t
            t(j+1) = t_trac
            t_trac_old = t_trac

        ELSE IF (min_t == dtz) THEN
            ix(j+1) = ix(j) 
            iy(j+1) = iy(j)
            IF (vz > 0)  iz(j+1) = iz(j) + 1
            IF (vz < 0)  iz(j+1) = iz(j) - 1
            t_trac = t_trac + min_t
            t(j+1) = t_trac
            t_trac_old = t_trac           

          ELSE IF (min_t == (dt-sub_dt)) THEN
            t_trac = t_trac + min_t
            t(j+1) = t_trac
            t_trac_old = t_trac
            IF ((ABS(dtx-min_t) < eps).AND.(ABS(dty-min_t) < eps)) THEN
              IF(vx > 0) ix(j+1) = ix(j) + 1
              IF(vx < 0) ix(j+1) = ix(j) - 1
              IF(vy > 0) iy(j+1) = iy(j) + 1
              IF(vy < 0) iy(j+1) = iy(j) - 1
            ELSE IF (ABS(dtx-min_t) < eps) THEN
              IF(vx > 0) ix(j+1) = ix(j) + 1
              IF(vx < 0) ix(j+1) = ix(j) - 1
            ELSE IF (ABS(dty-min_t) < eps) THEN
              IF(vy > 0) iy(j+1) = iy(j) + 1
              IF(vy < 0) iy(j+1) = iy(j) - 1
            ELSE
            j = j - 1
          END IF
            
        END IF
          x_trac = x_trac + vx * min_t
          y_trac = y_trac + vy * min_t
          z_trac = z_trac + vz * min_t
          sub_dt = sub_dt + min_t    
      END IF

      IF ( j+1 == np_trac) THEN
        np_trac = (j+1) * 2
        CALL arr_resize(ix,np_trac)
        CALL arr_resize(iy,np_trac)
        CALL arr_resize(iz,np_trac)
        CALL arr_resize(t,np_trac)
      END IF
      
      IF (t(j+1) == t(j)) THEN  
        ix(j) = ix(j+1)
        iy(j) = iy(j+1)
        iz(j) = iz(j+1)
        j = j - 1
      END IF
       
      j = j + 1
      
    END DO

  IF ( ((io >= 1).AND.nexpout).OR.(i+1 == np)) THEN
      t_out(io) = end_time
    END IF

  END IF

    
  END DO
  

  t_trac = t_in(1)
  np = j
  IF (j > 1) THEN
    CALL arr_resize(ix,np)
    CALL arr_resize(iy,np)
    CALL arr_resize(iz,np)
    CALL arr_resize(t,np+1)
    t(np+1) = t_out(io)
  END IF

END SUBROUTINE lagrangian_tracker
! -----------------------------------------------------------------------------
! Subroutine custom_emission_typ3:
!
SUBROUTINE custom_emission_typ3(edt,emd,emdT3,time,conditionT3)

    REAL, INTENT(in) :: time
    CHARACTER(len=50), PARAMETER :: name = "cloud_seeding_typ3"
    INTEGER, INTENT(in) :: conditionT3
    TYPE(EmitSizeDist), INTENT(in) :: edt       ! Emission data instance
    TYPE(EmitConfig), INTENT(in) :: emd         ! Emission configuration instance
    TYPE(EmitType3Config), INTENT(inout) :: emdT3  ! Emission type 3 configuration instance
    
    INTEGER :: j,bb,ss,mm
    REAL :: dt, t_str,t_end
    INTEGER :: ind, i_str,i_end, di
    CHARACTER(len=30) :: emit_spec

     IF (myid == 0) THEN
       WRITE(*,*) '========================'
       WRITE(*,*) 'CALCULATING EMISSIONS TYPE 3'
       WRITE(*,*) '========================'
     END IF
  
    ASSOCIATE( ix => emdT3%ix, iy => emdT3%iy, iz => emdT3%iz, t => emdT3%t, np => emdT3%np, &
              t_trac => emdT3%t_trac,t_in => emdT3%t_in, t_out => emdT3%t_out, &
              z_expan_up => emd%z_expan_up, z_expan_dw => emd%z_expan_dw)
      
    t_str = MAX(t_in(conditionT3), t_trac)
    t_end = MIN(t_out(conditionT3), (time + dtl) )  
    i_str = MAXLOC(t, 1, mask = t <= t_str)
    i_end = MAXLOC(t, 1, mask = t < t_end)

    di = i_end - i_str + 1
    
    DO bb = 1,nbins   
      DO j = 1, di
        dt  = ( MIN(t_end, t(i_str+j)) - MAX(t_str, t(i_str+j-1)) )/dtl
        ind = i_str+j-1
        a_naerot((iz(ind)-z_expan_dw):(iz(ind)+z_expan_up),ix(ind),iy(ind),bb) = &
                       a_naerot((iz(ind)-z_expan_dw):(iz(ind)+z_expan_up),ix(ind),iy(ind),bb) + edt%numc(bb) * dt
        DO ss = 1,spec%getNSpec()
          mm = getMassIndex(nbins,bb,ss)
          a_maerot((iz(ind)-z_expan_dw):(iz(ind)+z_expan_up),ix(ind),iy(ind),mm) = &
                     a_maerot((iz(ind)-z_expan_dw):(iz(ind)+z_expan_up),ix(ind),iy(ind),mm) + edt%mass(mm) * dt
        END DO
      END DO
    END DO

    t_trac = time + dtl
         
    END ASSOCIATE
  END SUBROUTINE custom_emission_typ3
  
  ! -----------------------------------------------------------------------------
  FUNCTION getConditionT3(emdT3,time)
    IMPLICIT NONE
    INTEGER :: getConditionT3
    INTEGER :: N, i
    CLASS(EmitType3Config), INTENT(in) :: emdT3
    REAL, INTENT(in) :: time
    CHARACTER(len=50), PARAMETER :: name = "getConditiontyp3"
    
    ASSOCIATE(t_in => emdT3%t_in, t_out => emdT3%t_out)
      
      getConditionT3 = 0
      N = SIZE(t_in)

      DO i =1,N
        IF (t_in(i) <= time .AND. t_out(i) > time) getConditionT3 = i
      END DO
       
    END ASSOCIATE
    
  END FUNCTION getConditionT3
! -----------------------------------------------------------------------------
!
SUBROUTINE arr_resize_int(arr,newsize)
  INTEGER, ALLOCATABLE, INTENT(inout) :: arr(:)
  INTEGER, INTENT(inout) :: newsize
  INTEGER, ALLOCATABLE  :: tmp_arr(:)
  INTEGER :: oldsize

  oldsize = SIZE(arr)
  ALLOCATE(tmp_arr(newsize))
  
  IF (oldsize >= newsize) THEN
    tmp_arr = arr(1:newsize)
    DEALLOCATE (arr)
    CAll MOVE_ALLOC(tmp_arr,arr)
  ELSE
    tmp_arr(1:oldsize) = arr
    DEALLOCATE (arr)
    CAll MOVE_ALLOC(tmp_arr,arr)
  END IF
    
END SUBROUTINE arr_resize_int
! -----------------------------------------------------------------------------

SUBROUTINE arr_resize_rel(arr,newsize)
  REAL, ALLOCATABLE, INTENT(inout) :: arr(:)
  INTEGER, INTENT(in) :: newsize
  REAL, ALLOCATABLE  :: tmp_arr(:)
  INTEGER :: oldsize

  oldsize = SIZE(arr)
  ALLOCATE(tmp_arr(newsize))
  
  IF (oldsize >= newsize) THEN
    tmp_arr = arr(1:newsize)
    DEALLOCATE (arr)
    CAll MOVE_ALLOC(tmp_arr,arr)
  ELSE
    tmp_arr(1:oldsize) = arr
    DEALLOCATE (arr)
    CAll MOVE_ALLOC(tmp_arr,arr)
  END IF
    
END SUBROUTINE arr_resize_rel
  
! -----------------------------------------------------------------------------
!Ali, addition of emission type 3 
        
END MODULE emission_main
