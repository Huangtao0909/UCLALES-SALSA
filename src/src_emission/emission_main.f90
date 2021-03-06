MODULE emission_main
  use emission_types, ONLY : EmitConfig, EmitSizeDist, EmitType3Config,   &
                             emitModes, emitData, emitType3, nEmissionModes
  
  USE mo_seasalt_emission

  USE mo_submctl, ONLY : pi6, in1a, fn2a, in2b, fn2b, nbins, spec, pi6

  USE mo_salsa_types, ONLY : aero

  USE mo_salsa_sizedist, ONLY : size_distribution  ! Could this be packaged somehow differently?

  USE mo_aux_state, ONLY : dzt,zt,xt,yt
  USE mo_progn_state, ONLY : a_maerot, a_naerot
  USE mo_diag_state, ONLY : a_dn
  !USE mo_vector_state, ONLY : a_up, a_vp ! needed for the seasalt thing
  USE grid, ONLY: deltax, deltay, deltaz, dtl, &                  
                  nxp,nyp,nzp
    
  USE util, ONLY: smaller, closest, getMassIndex
  USE exceptionHandling, ONLY: errorMessage
  USE mpi_interface, ONLY : myid
  
  IMPLICIT NONE

  CHARACTER(len=50), PARAMETER :: global_name = "emission_main"
   
  CONTAINS

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
    
    ! Loop over all specified emission profiles
    DO pr = 1,nEmissionModes
       ASSOCIATE(emd => emitModes(pr), edt => emitData(pr))
         IF (emd%emitType == 2) THEN
            condition = getCondition(emd,time_in)
            IF (condition) CALL custom_emission(edt,emd)
         END IF
         
         !Ali, addition of emitType 3 
         IF (emd%emitType == 3) THEN
            emdT3 => emitType3(pr)
            conditionT3 = getConditionT3(emdT3,time_in)
            IF (conditionT3 > 0) THEN
               CALL custom_emission_typ3(edt,emd,emdT3,time_in,conditionT3)
            END IF
         END IF
         
       END ASSOCIATE
       ! Ali, addition of emission type
       emdT3 => NULL()
       
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

    INTEGER :: i,j,bb,ss,mm
    
    IF (myid == 0) THEN
      WRITE(*,*) '========================'
      WRITE(*,*) 'CALCULATING EMISSIONS'
      WRITE(*,*) '========================'
    END IF

    ASSOCIATE( k1 => emd%emitLevMin, k2 => emd%emitLevMax )
    
      DO bb = 1,nbins
         DO j = 1,nyp
            DO i = 1,nxp
               a_naerot%d(k1:k2,i,j,bb) = a_naerot%d(k1:k2,i,j,bb) + edt%numc(bb)
               DO ss = 1,spec%getNSpec(type="wet")
                  mm = getMassIndex(nbins,bb,ss)
                  a_maerot%d(k1:k2,i,j,mm) = a_maerot%d(k1:k2,i,j,mm) + edt%mass(mm)
               END DO
            END DO
         END DO
      END DO
      
    END ASSOCIATE

  END SUBROUTINE custom_emission
 
 
! -----------------------------------------------------------------------------
! Subroutine custom_emission_typ3:
!
  SUBROUTINE custom_emission_typ3(edt,emd,emdT3,time,conditionT3)
    IMPLICIT NONE
  
    REAL, INTENT(in) :: time
    CHARACTER(len=50), PARAMETER :: name = "cloud_seeding_typ3"
    INTEGER, INTENT(in) :: conditionT3
    TYPE(EmitSizeDist), INTENT(in) :: edt       ! Emission data instance
    TYPE(EmitConfig), INTENT(in) :: emd         ! Emission configuration instance
    TYPE(EmitType3Config), INTENT(inout) :: emdT3  ! Emission type 3 configuration instance
    
    INTEGER :: j,bb,ss,mm
    REAL :: dt, t_str,t_end
    INTEGER :: ind, i_str,i_end, di
    
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
            a_naerot%d((iz(ind)-z_expan_dw):(iz(ind)+z_expan_up),ix(ind),iy(ind),bb) = &
                 a_naerot%d((iz(ind)-z_expan_dw):(iz(ind)+z_expan_up),ix(ind),iy(ind),bb) + edt%numc(bb) * dt
            DO ss = 1,spec%getNSpec(type="wet")
               mm = getMassIndex(nbins,bb,ss)
               a_maerot%d((iz(ind)-z_expan_dw):(iz(ind)+z_expan_up),ix(ind),iy(ind),mm) = &
                    a_maerot%d((iz(ind)-z_expan_dw):(iz(ind)+z_expan_up),ix(ind),iy(ind),mm) + edt%mass(mm) * dt
            END DO
         END DO
      END DO
      
      t_trac = time + dtl
      
    END ASSOCIATE
  END SUBROUTINE custom_emission_typ3
  
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
        
END MODULE emission_main
