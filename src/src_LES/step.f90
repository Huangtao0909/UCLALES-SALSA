!----------------------------------------------------------------------------
! This file is part of UCLALES.
!
! UCLALES is free software; you can redistribute it and/or modify
! it under the terms of the GNU General Public License as published by
! the Free Software Foundation; either version 3 of the License, or
! (at your option) any later version.
!
! UCLALES is distributed in the hope that it will be useful,
! but WITHOUT ANY WARRANTY; without even the implied warranty of
! MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
! GNU General Public License for more details.
!
! You should have received a copy of the GNU General Public License
! along with this program.  If not, see <http://www.gnu.org/licenses/>.
!
! Copyright 199-2007, Bjorn B. Stevens, Dep't Atmos and Ocean Sci, UCLA
!----------------------------------------------------------------------------
!
module step

  implicit none

  integer :: istpfl = 1
  real    :: timmax = 18000.
  logical :: corflg = .false.

  real    :: frqhis =  9000.
  real    :: frqrst =  3600.
  real    :: frqanl =  3600.
  real    :: anl_start = -1.
  real    :: radfrq =  0.

  real    :: time   =  0.
  real    :: strtim =  0.0
  real    :: cntlat =  31.5 ! 30.0
  logical :: outflg = .true.


contains
  !
  ! ----------------------------------------------------------------------
  ! Subroutine model:  This is the main driver for the model's time
  ! integration.  It calls the routine tstep, which steps through the
  ! physical processes active on a time-step and updates variables.  It
  ! then checks to see whether or not different output options are
  ! satisfied.
  subroutine stepper

    use mpi_interface, only : myid, double_scalar_par_max

    use grid, only : dtl, dzt, zt, zm, nzp, dn0, u0, v0, &
         write_hist, write_anal, close_anal, &
         dtlong, nzp, nyp, nxp, level, &
         ! For mass budged
         a_rp, a_rc, a_srp, a_dn

    use stat, only : sflg, savg_intvl, ssam_intvl, write_ps, close_stat, mcflg, acc_massbudged,  &
         write_massbudged
    use thrm, only : thermo

    real, parameter :: cfl_upper = 0.5

    real    :: t1,t2,tplsdt
    REAL(kind=8) :: cflmax,gcflmax
    integer :: istp, iret
    !
    ! Timestep loop for program
    !
    istp = 0

    call cpu_time(t1)

    do while (time < timmax)
       ! Limit time step based on the Courant-Friedrichs-Lewy condition
       call cfl(cflmax)
       call double_scalar_par_max(cflmax,gcflmax)
       cflmax = gcflmax
       dtl = min(dtlong,dtl*cfl_upper/(cflmax+epsilon(1.)))

       ! Determine when to compute statistics
       !    - When a given output or profile sampling time (n*tstep) will be reached or exceeded for the first time
       !    - After the first call (time=0)
       tplsdt = time + dtl ! Time after t_step
       sflg = (min(mod(tplsdt,ssam_intvl),mod(tplsdt,savg_intvl)) < dtl .or. tplsdt < 1.1*dtl)

       call t_step
       time = time + dtl

       ! Write profiles
       if (mod(tplsdt,savg_intvl) < dtl .or. tplsdt < 1.1*dtl)   &
            call write_ps(nzp,dn0,u0,v0,zm,zt,time)

       ! Write restarts (*.<time>s and *.rst)
       if (mod(tplsdt,frqhis) < dtl .and. outflg)   &
            call write_hist(2, time)

       if (mod(tplsdt,frqrst) < dtl .and. outflg)   &
            call write_hist(1, time)

       ! Write analysis files
       if (mod(tplsdt,frqanl) < dtl .and. outflg .and. time >= anl_start)   &
            call write_anal(time)

       if(myid == 0) then
          istp = istp+1
          if (mod(istp,istpfl) == 0 ) THEN
              call cpu_time(t2) ! t2-t1 is the actual CPU time from the previous output
              print "('   Timestep # ',i6," //     &
                 "'   Model time(sec)=',f10.2,3x,'CPU time(sec)=',f8.3)",     &
                 istp, time, t2-t1
              call cpu_time(t1)
          ENDIF
       endif

    enddo

    IF (mcflg) THEN
       !
       ! Juha:
       ! Get the final statistics of atmospheric water for mass budged
       CALL acc_massbudged(nzp,nxp,nyp,1,dtl,dzt,a_dn,    &
            rv=a_rp,rc=a_rc,prc=a_srp)

       CALL write_massbudged

    END IF ! mcflg

    call write_hist(1, time)
    iret = close_anal()
    iret = close_stat()

  end subroutine stepper
  !
  !----------------------------------------------------------------------
  ! subroutine t_step: Called by driver to timestep through the LES
  ! routines.  Within many subroutines, data is accumulated during
  ! the course of a timestep for the purposes of statistical analysis.
  !
  subroutine t_step()

    use grid, only : level, dtl, Tspinup,                                         &
                     ! Added parameters for interfacing with SALSA
                     nxp, nyp, nzp, a_press, a_temp, a_rsl,                       &
                     a_rc, a_wp, a_rp, a_rt, a_rh,                                  &
                     a_naerop, a_naerot, a_ncloudp, a_ncloudt, a_nprecpp, a_nprecpt,    &
                     a_maerop, a_maerot, a_mcloudp, a_mcloudt, a_mprecpp, a_mprecpt,    &
                     a_nicep,  a_nicet,  a_micep,  a_micet,                             &
                     a_nsnowp, a_nsnowt, a_msnowp, a_msnowt,                            &
                     a_gaerop, a_gaerot, a_dn, a_nactd, a_vactd, prtcl, sst, a_rsi,     &
                     nudge_theta, nudge_rv, nudge_u, nudge_v, nudge_ccn, &
                     coag_ra, coag_na, coag_rc, coag_nc, coag_rr, coag_nr, coag_ri, coag_ni, coag_rs, coag_ns, &
                     cond_ra, cond_rc, cond_rr, cond_ri, cond_rs, auto_rr, auto_nr, auto_rs, auto_ns, &
                     cact_rc, cact_nc, nucl_ri, nucl_ni, melt_ri, melt_ni, melt_rs, melt_ns

    use stat, only : sflg, statistics
    use sgsm, only : diffuse
    use srfc, only : surface
    use thrm, only : thermo
    use mcrp, only : micro
    use prss, only : poisson
    use advf, only : fadvect, newdroplet
    use advl, only : ladvect
    use forc, only : forcings
    USE util, ONLY : maskactiv !Juha: Included for SALSA

    USE mo_salsa_driver, ONLY : run_SALSA
    USE class_ComponentIndex, ONLY : GetNcomp

    real :: xtime

    LOGICAL :: zactmask(nzp,nxp,nyp)
    REAL :: zwp(nzp,nxp,nyp)  !! FOR SINGLE-COLUMN RUNS
    INTEGER :: zrm

    INTEGER :: n4

    zwp = 0.5

    xtime = time/86400. + strtim

    ! The runmode parameter zrm is used by SALSA only
    zrm = 3
    IF ( time < Tspinup ) zrm = 2

    ! Reset ALL tendencies here.
    !----------------------------------------------------------------
    ! "Scalar" timestep
    CALL tend0(.FALSE.)

    ! Put the newly activated to zero
    IF (level >= 4) THEN
       a_vactd = 0.
       a_nactd = 0.
    END IF

    call surface(sst)

    call diffuse

    call sponge(0)

    if (level >= 1) then

       call thermo(level)

       call forcings(xtime,cntlat,sst)

       IF (level >= 4) THEN

          n4 = GetNcomp(prtcl) + 1 ! Aerosol components + water
          CALL tend_constrain(n4)
          call update_sclrs
          CALL tend0(.TRUE.)

          IF ( nxp ==5 .and. nyp == 5 ) THEN
             ! 1D -runs
             CALL run_SALSA(nxp,nyp,nzp,n4,a_press,a_temp,a_rp,a_rt,a_rsl,a_rsi,zwp,a_dn,  &
                  a_naerop,  a_naerot,  a_maerop,  a_maerot,   &
                  a_ncloudp, a_ncloudt, a_mcloudp, a_mcloudt,  &
                  a_nprecpp, a_nprecpt, a_mprecpp, a_mprecpt,  &
                  a_nicep,   a_nicet,   a_micep,   a_micet,    &
                  a_nsnowp,  a_nsnowt,  a_msnowp,  a_msnowt,   &
                  a_nactd,   a_vactd,   a_gaerop,  a_gaerot,   &
                  zrm, prtcl, dtl, time, level,  &
                  coag_ra, coag_na, coag_rc, coag_nc, coag_rr, coag_nr, &
                  coag_ri, coag_ni, coag_rs, coag_ns, &
                  cond_ra, cond_rc, cond_rr, cond_ri, cond_rs, &
                  auto_rr, auto_nr, auto_rs, auto_ns, &
                  cact_rc, cact_nc, nucl_ri, nucl_ni, &
                  melt_ri, melt_ni, melt_rs, melt_ns)
          ELSE
             !! for 2D or 3D runs
             CALL run_SALSA(nxp,nyp,nzp,n4,a_press,a_temp,a_rp,a_rt,a_rsl,a_rsi,a_wp,a_dn,  &
                  a_naerop,  a_naerot,  a_maerop,  a_maerot,   &
                  a_ncloudp, a_ncloudt, a_mcloudp, a_mcloudt,  &
                  a_nprecpp, a_nprecpt, a_mprecpp, a_mprecpt,  &
                  a_nicep,   a_nicet,   a_micep,   a_micet,    &
                  a_nsnowp,  a_nsnowt,  a_msnowp,  a_msnowt,   &
                  a_nactd,   a_vactd,   a_gaerop,  a_gaerot,   &
                  zrm, prtcl, dtl, time, level,  &
                  coag_ra, coag_na, coag_rc, coag_nc, coag_rr, coag_nr, &
                  coag_ri, coag_ni, coag_rs, coag_ns, &
                  cond_ra, cond_rc, cond_rr, cond_ri, cond_rs, &
                  auto_rr, auto_nr, auto_rs, auto_ns, &
                  cact_rc, cact_nc, nucl_ri, nucl_ni, &
                  melt_ri, melt_ni, melt_rs, melt_ns)
          END IF !nxp==5 and nyp == 5

          CALL tend_constrain(n4)
       END IF

    end if ! level

    call update_sclrs

    !-------------------------------------------
    ! "Deposition" timestep
    ! -- Reset only scalar tendencies
    CALL tend0(.TRUE.)

    ! Dont perform sedimentation or level 3 autoconversion during spinup
    IF (zrm == 3) THEN
        CALL micro(level)

        IF (level >= 4) CALL tend_constrain(n4)
        CALL update_sclrs
    END IF

    !-------------------------------------------
    ! "Advection" timestep
    ! -- Reset only scalar tendencies
    call tend0(.TRUE.)

    ! Mask for cloud base activation
    IF (level >= 4)  CALL maskactiv(zactmask,nxp,nyp,nzp,2,a_rh,rc=a_rc,w=a_wp)
    ! Get tendencies from cloud base activation
    IF (level >= 4) CALL newdroplet(zactmask)

    CALL fadvect

    IF (level >= 4)  &
         CALL tend_constrain(n4)

    CALL update_sclrs

    CALL thermo(level)

    ! Nudging
     IF (nudge_theta/=0 .OR. nudge_rv/=0 .OR. nudge_u/=0 .OR. &
            nudge_v/=0 .OR. (level>3 .AND. nudge_ccn/=0) ) THEN

        ! Reset tendencies
        call tend0(.TRUE.)

        ! Update diagnostic tracers
        IF (level >= 4)  THEN
             CALL SALSA_diag_update
             call thermo(level)
        ENDIF

        CALL nudging(time)

        CALL update_sclrs

        CALL thermo(level)
    ENDIF

    IF (level >= 4)  THEN
         CALL SALSA_diagnostics(.true.)
         call thermo(level)
    ENDIF

    call corlos

    call ladvect

    call buoyancy

    call sponge(1)

    call poisson

    CALL thermo(level)

    IF (level >= 4)  THEN
         CALL SALSA_diagnostics(.false.)
         call thermo(level)
    ENDIF

    if (sflg) then
       call statistics (time+dtl)
    end if

  end subroutine t_step
  !
  !----------------------------------------------------------------------
  !
  ! Nudging towards the initial state (temperature, water vapor,
  ! horizontal winds and aerosol and/or cloud droplets).
  !
  ! TR 22.3.2017
  !
  SUBROUTINE nudging(time)

    use grid, only : level, dtl, nxp, nyp, nzp, &
                zt, a_rp, a_rt, a_rc, a_srp, a_ri, a_srs, &
                a_naerop, a_naerot, a_ncloudp, a_nicep, &
                a_tp, a_tt, a_up, a_ut, a_vp, a_vt, &
                !th0, th00, rt0, u0, v0, &
                nudge_theta, nudge_theta_time, nudge_theta_zmin, nudge_theta_zmax, nudge_theta_tau, &
                nudge_rv, nudge_rv_time, nudge_rv_zmin, nudge_rv_zmax, nudge_rv_tau, &
                nudge_u, nudge_u_time, nudge_u_zmin, nudge_u_zmax, nudge_u_tau, &
                nudge_v, nudge_v_time, nudge_v_zmin, nudge_v_zmax, nudge_v_tau, &
                nudge_ccn, nudge_ccn_time, nudge_ccn_zmin, nudge_ccn_zmax, nudge_ccn_tau, &
                theta_ref, rv_ref, u_ref, v_ref, aero_ref, nudge_init
    USE mo_submctl, ONLY : nbins, ncld, nice, in2a, fn2b

    IMPLICIT NONE
    REAL, INTENT(IN) :: time
    REAL :: aero_target(nzp,nxp,nyp,nbins)

    ! Initialization
    IF (nudge_init) THEN
        ! Note: the first temperature and humidity values can include random
        ! perturbations, so could take the target values from soundings (th0, rt0).
        ! There are no wind perturbations, but can still could use u0 and v0.
        !
        ! (Liquid water) potential temperature: nudge towards initial theta
        IF (nudge_theta/=0) THEN
            ALLOCATE(theta_ref(nzp))
            theta_ref(:)=a_tp(:,3,3)
            !theta_ref(:)=th0(:)-th00 ! Initial state from soundings
        ENDIF
        !
        ! Water vapor mixing ratio based on total water
        !   Levels 0-3: total = total water (a_rp)
        !   Levels 4-5: total = water vapor (a_rp) + aerosol and cloud water (a_rc) + rain water (a_srp)
        !                            + ice (a_ri) + snow (a_srs)
        IF (nudge_rv/=0)  THEN
            ALLOCATE(rv_ref(nzp))
            IF (level>3) THEN
                rv_ref(:)=a_rp(:,3,3)+a_rc(:,3,3)+a_srp(:,3,3)+a_ri(:,3,3)+a_srs(:,3,3)
            ELSE ! Levels 0-3
                rv_ref(:)=a_rp(:,3,3) ! This includes all
            ENDIF
            !rv_ref(:)=rt0(:) ! Initial state from soundings
        ENDIF
        !
        ! Horizontal winds
        IF (nudge_u/=0) THEN
            ALLOCATE(u_ref(nzp))
            u_ref(:)=a_up(:,3,3)
            !u_ref(:)=u0(:) ! Initial state from soundings
        ENDIF
        IF (nudge_v/=0) THEN
            ALLOCATE(v_ref(nzp))
            v_ref(:)=a_vp(:,3,3)
            !v_ref(:)=v0(:) ! Initial state from soundings
        ENDIF
        !
        ! Nudge level 4 and 5 aerosol concentration based on total CCN = aerosol + cloud droplets + ice.
        ! Precipitation and snow are not included, because these cannot be related to a specific aerosol bin.
        IF (level>3 .AND. nudge_ccn/=0) THEN
            ! Nudge aerosol based on the total number (aerosol+cloud+ice)
            ALLOCATE(aero_ref(nzp,nbins))
            aero_ref(:,:)=a_naerop(:,3,3,:)
            aero_ref(:,in2a:fn2b)=aero_ref(:,in2a:fn2b)+a_ncloudp(:,3,3,1:ncld)
            IF (level==5) aero_ref(:,in2a:fn2b)=aero_ref(:,in2a:fn2b)+a_nicep(:,3,3,1:nice)
        ENDIF
        !
        ! Initialized
        nudge_init=.FALSE.
    ENDIF

    ! (Liquid water) potential temperature:
    IF (nudge_theta>0) &
        CALL nudge_any(nxp,nyp,nzp,zt,a_tp,a_tt,theta_ref,dtl,time,nudge_theta, &
            nudge_theta_time,nudge_theta_zmin,nudge_theta_zmax,nudge_theta_tau)

    ! Water vapor
    IF (nudge_rv>0) THEN
        IF (level>3) THEN
            ! Nudge water vapor (a_rp) based on total (vapor + cloud + rain [+ ice + snow])
            CALL nudge_any(nxp,nyp,nzp,zt,a_rp+a_rc+a_srp+a_ri+a_srs,a_rt,rv_ref,dtl,time,nudge_rv, &
                nudge_rv_time,nudge_rv_zmin,nudge_rv_zmax,nudge_rv_tau)
        ELSE
            ! Nudge total water (a_rp) based on total water
            CALL nudge_any(nxp,nyp,nzp,zt,a_rp,a_rt,rv_ref,dtl,time,nudge_rv, &
                nudge_rv_time,nudge_rv_zmin,nudge_rv_zmax,nudge_rv_tau)
        ENDIF
    ENDIF

    ! Horizontal winds
    IF (nudge_u>0) &
         CALL nudge_any(nxp,nyp,nzp,zt,a_up,a_ut,u_ref,dtl,time,nudge_u, &
            nudge_u_time,nudge_u_zmin,nudge_u_zmax,nudge_u_tau)
    IF (nudge_v>0) &
        CALL nudge_any(nxp,nyp,nzp,zt,a_vp,a_vt,v_ref,dtl,time,nudge_v, &
            nudge_v_time,nudge_v_zmin,nudge_v_zmax,nudge_v_tau)

    ! Aerosol
    IF (level>3 .AND. nudge_ccn/=0) THEN
        ! Target aerosol concentration = aerosol(t)+cloud(t)+ice(t)
        aero_target(:,:,:,:)=a_naerop(:,:,:,:)
        aero_target(:,:,:,in2a:fn2b)=aero_target(:,:,:,in2a:fn2b)+a_ncloudp(:,:,:,1:ncld)
        IF (level==5) aero_target(:,:,:,in2a:fn2b)=aero_target(:,:,:,in2a:fn2b)-a_nicep(:,:,:,1:nice)
        ! Apply to sectional data
        CALL nudge_any_2d(nxp,nyp,nzp,nbins,zt,aero_target,a_naerot,aero_ref,dtl,time,nudge_ccn, &
            nudge_ccn_time,nudge_ccn_zmin,nudge_ccn_zmax,nudge_ccn_tau)
    ENDIF

  END SUBROUTINE nudging
  !
  ! Nudging for any 3D field based on 1D target
  SUBROUTINE nudge_any(nxp,nyp,nzp,zt,ap,at,trgt,dt,time,iopt,tref,zmin,zmax,tau)
    USE util, ONLY : get_avg3
    IMPLICIT NONE
    INTEGER :: nxp,nyp,nzp
    REAL :: zt(nzp), ap(nzp,nxp,nyp), at(nzp,nxp,nyp)
    REAL :: dt, time
    REAL :: trgt(nzp)
    REAL :: tref,zmin,zmax,tau
    INTEGER :: iopt
    INTEGER :: kk
    REAL :: avg(nzp)
    !
    IF (iopt==1 .AND. time<=tref) THEN
        ! Soft nudging with fixed nudging constant (tau [s]), and the time parameter gives the maximum nudging time
        CALL get_avg3(nzp,nxp,nyp,ap,avg)
        DO kk = 1,nzp
            IF (zmin<=zt(kk) .AND. zt(kk)<=zmax) &
                at(kk,:,:)=at(kk,:,:)-(avg(kk)-trgt(kk))/max(tau,dt)
        ENDDO
    ELSEIF (iopt==2 .AND. time<=tref) THEN
        ! Hard nudging with fixed nudging constant (tau [s]), and the time parameter gives the maximum nudging time
        DO kk = 1,nzp
            IF (zmin<=zt(kk) .AND. zt(kk)<=zmax) &
                at(kk,:,:)=at(kk,:,:)-(ap(kk,:,:)-trgt(kk))/max(tau,dt)
        ENDDO
    ENDIF
    !
  END SUBROUTINE nudge_any
  !
  ! Nudging for any 4D field based on 2D target
  SUBROUTINE nudge_any_2d(nxp,nyp,nzp,nb,zt,ap,at,trgt,dt,time,iopt,tref,zmin,zmax,tau)
    USE util, ONLY : get_avg3
    IMPLICIT NONE
    INTEGER :: nxp,nyp,nzp,nb
    REAL :: zt(nzp), ap(nzp,nxp,nyp,nb), at(nzp,nxp,nyp,nb)
    REAL :: dt, time
    REAL :: trgt(nzp,nb)
    REAL :: tref,zmin,zmax,tau
    INTEGER :: iopt
    INTEGER :: ii, kk
    REAL :: avg(nzp)
    !
    IF (iopt==1 .AND. time<=tref) THEN
        ! Soft nudging with fixed nudging constant (tau [s]), and the time parameter gives the maximum nudging time
        DO ii=1,nb
            CALL get_avg3(nzp,nxp,nyp,ap(:,:,:,ii),avg)
            DO kk = 1,nzp
                IF (zmin<=zt(kk) .AND. zt(kk)<=zmax) &
                    at(kk,:,:,ii)=at(kk,:,:,ii)-(avg(kk)-trgt(kk,ii))/max(tau,dt)
            ENDDO
        ENDDO
    ELSEIF (iopt==2 .AND. time<=tref) THEN
        ! Hard nudging with fixed nudging constant (tau [s]), and the time parameter gives the maximum nudging time
        DO ii=1,nb
            DO kk = 1,nzp
                IF (zmin<=zt(kk) .AND. zt(kk)<=zmax) &
                    at(kk,:,:,ii)=at(kk,:,:,ii)-(ap(kk,:,:,ii)-trgt(kk,ii))/max(tau,dt)
            ENDDO
        ENDDO
    ENDIF
    !
  END SUBROUTINE nudge_any_2d
  !
  !----------------------------------------------------------------------
  ! subroutine tend0: sets all tendency arrays to zero
  !
  subroutine tend0(sclonly)

    use grid, only : a_ut, a_vt, a_wt, nscl, a_st, newsclr

    LOGICAL, INTENT(in) :: sclonly ! If true, only put scalar tendencies to zero

    integer :: n

    IF( .NOT. sclonly) THEN
       a_ut=0.; a_vt=0.; a_wt=0.
    ENDIF
    do n=1,nscl
       call newsclr(n)
       a_st=0.
    end do

  end subroutine tend0
  !
  !----------------------------------------------------------------------
  ! In case of negative tendencies to SALSA arrays, put some constrains
  ! in order to avoid concentrations going negative. This will possibly
  ! slightly affect the conservation of mass - needs testing/revision
  ! Juha Tonttila, FMI, 2014
  !
  SUBROUTINE tend_constrain(nn)

    USE grid, ONLY : a_naerop, a_naerot, a_ncloudp, a_ncloudt, a_nprecpp, a_nprecpt,   &
                     a_maerop, a_maerot, a_mcloudp, a_mcloudt, a_mprecpp, a_mprecpt,   &
                     a_nicep, a_nicet, a_nsnowp, a_nsnowt, a_micep, a_micet, a_msnowp, a_msnowt,  &
                     dtl, nxp,nyp,nzp,level
    USE mo_submctl, ONLY : nbins, ncld, nprc, nice, nsnw

    INTEGER, INTENT(in) :: nn

    INTEGER :: cc, ii,jj,kk,ni

    DO jj = 3,nyp-2

       DO ii = 3,nxp-2

          DO kk = 1,nzp

             ! Aerosols
             DO cc = 1,nbins

                IF ( a_naerop(kk,ii,jj,cc)+a_naerot(kk,ii,jj,cc)*dtl < 0. ) THEN

                   a_naerot(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_naerop(kk,ii,jj,cc))/dtl,a_naerot(kk,ii,jj,cc))
                   DO ni = 1,nn
                      a_maerot(kk,ii,jj,(ni-1)*nbins+cc) = MAX( ((1.e-10-1.0)*a_maerop(kk,ii,jj,(ni-1)*nbins+cc))/dtl,  &
                                                               a_maerot(kk,ii,jj,(ni-1)*nbins+cc) )
                   END DO

                END IF

             END DO

             ! Cloud droplets
             DO cc = 1,ncld

                IF ( a_ncloudp(kk,ii,jj,cc)+a_ncloudt(kk,ii,jj,cc)*dtl < 0. ) THEN

                   a_ncloudt(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_ncloudp(kk,ii,jj,cc))/dtl,a_ncloudt(kk,ii,jj,cc))
                   DO ni = 1,nn
                      a_mcloudt(kk,ii,jj,(ni-1)*ncld+cc) = MAX( ((1.e-10-1.0)*a_mcloudp(kk,ii,jj,(ni-1)*ncld+cc))/dtl,  &
                                                               a_mcloudt(kk,ii,jj,(ni-1)*ncld+cc) )
                   END DO

                END IF

             END DO

             ! Precipitation
             DO cc = 1,nprc

                IF ( a_nprecpp(kk,ii,jj,cc)+a_nprecpt(kk,ii,jj,cc)*dtl < 0. ) THEN

                   a_nprecpt(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_nprecpp(kk,ii,jj,cc))/dtl,a_nprecpt(kk,ii,jj,cc))
                   DO ni = 1,nn
                      a_mprecpt(kk,ii,jj,(ni-1)*nprc+cc) = MAX( ((1.e-10-1.0)*a_mprecpp(kk,ii,jj,(ni-1)*nprc+cc))/dtl,  &
                                                               a_mprecpt(kk,ii,jj,(ni-1)*nprc+cc) )
                   END DO

                END IF

             END DO

             ! ice particles
             IF (level<5) CYCLE
             DO cc = 1,nice

                IF ( a_nicep(kk,ii,jj,cc)+a_nicet(kk,ii,jj,cc)*dtl < 0. ) THEN

                   a_nicet(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_nicep(kk,ii,jj,cc))/dtl,a_nicet(kk,ii,jj,cc))
                   DO ni = 1,nn
                      a_micet(kk,ii,jj,(ni-1)*ncld+cc) = MAX( ((1.e-10-1.0)*a_micep(kk,ii,jj,(ni-1)*nice+cc))/dtl,  &
                                                               a_micet(kk,ii,jj,(ni-1)*nice+cc) )
                   END DO

                END IF

             END DO

             ! Snow
             DO cc = 1,nsnw

                IF ( a_nsnowp(kk,ii,jj,cc)+a_nsnowt(kk,ii,jj,cc)*dtl < 0. ) THEN

                   a_nsnowt(kk,ii,jj,cc) = MAX(((1.e-10-1.0)*a_nsnowp(kk,ii,jj,cc))/dtl,a_nsnowt(kk,ii,jj,cc))
                   DO ni = 1,nn
                      a_msnowt(kk,ii,jj,(ni-1)*nsnw+cc) = MAX( ((1.e-10-1.0)*a_msnowp(kk,ii,jj,(ni-1)*nsnw+cc))/dtl,  &
                                                               a_msnowt(kk,ii,jj,(ni-1)*nsnw+cc) )
                   END DO

                END IF

             END DO

          END DO ! kk

       END DO ! ii

    END DO ! jj

  END SUBROUTINE tend_constrain
  !
  !----------------------------------------------------------------------
  ! Subroutine cfl: Driver for calling CFL computation subroutine
  !
  subroutine cfl(cflmax)

    use grid, only : a_up,a_vp,a_wp,nxp,nyp,nzp,dxi,dyi,dzt,dtl
    use stat, only : fill_scalar

    real(KIND=8), intent (out)   :: cflmax
    real, parameter :: cflnum=0.95

    cflmax =  cfll(nzp,nxp,nyp,a_up,a_vp,a_wp,dxi,dyi,dzt,dtl)

    if (cflmax > cflnum) print *, 'Warning CFL Violation :', cflmax
    call fill_scalar(1,REAL(cflmax))

  end subroutine cfl
  !
  !----------------------------------------------------------------------
  ! Subroutine cfll: Gets the peak CFL number
  !
  real(KIND=8) function cfll(n1,n2,n3,u,v,w,dxi,dyi,dzt,dtlt)

    integer, intent (in) :: n1, n2, n3
    real, dimension (n1,n2,n3), intent (in) :: u, v, w
    real, intent (in)    :: dxi,dyi,dzt(n1),dtlt

    integer :: i, j, k
    cfll=0.
    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1
             cfll=max(cfll, dtlt*2.* max(abs(u(k,i,j)*dxi),             &
                  abs(v(k,i,j)*dyi), abs(w(k,i,j)*dzt(k))))
          end do
       end do
    end do

  end function cfll
  !
  !----------------------------------------------------------------------
  ! subroutine update_sclrs:  Updates scalars by applying tendency and
  ! boundary conditions
  !
  subroutine update_sclrs

    use grid, only : a_sp, a_st, a_qp, nscl, nxyzp, nxp, nyp, nzp, dzt, &
         dtl, newsclr, isgstyp
    use sgsm, only : tkeinit
    use util, only : sclrset

    integer :: n

    do n=1,nscl
       call newsclr(n)
       call update(nzp,nxp,nyp,a_sp,a_st,dtl)
       call sclrset('mixd',nzp,nxp,nyp,a_sp,dzt)
    end do

    if (isgstyp == 2) then
       call tkeinit(nxyzp,a_qp)
    end if

  end subroutine update_sclrs
  !
  ! ----------------------------------------------------------------------
  ! subroutine update:
  !
  subroutine update(n1,n2,n3,a,fa,dt)

    integer, intent(in)   :: n1, n2, n3
    real, intent (in)     :: fa(n1,n2,n3),dt
    real, intent (in out) :: a(n1,n2,n3)
    integer :: i, j, k

    do j=3,n3-2
       do i=3,n2-2
          do k=2,n1-1
             a(k,i,j) = a(k,i,j) + fa(k,i,j)*dt
          end do
       end do
    end do

  end subroutine update
  !
  ! ----------------------------------------------------------------------
  ! subroutine buoyancy:
  !
  subroutine buoyancy

    use grid, only : a_uc, a_vc, a_wc, a_wt, a_rv, a_rc, a_theta, &
         a_rp, a_rpp, a_srp, a_ri, a_srs, nxp, nyp, nzp, dzm, th00, level, pi1
    use stat, only : sflg, comp_tke
    use util, only : ae1mm
    use thrm, only : update_pi1

    real :: awtbar(nzp), a_tmp1(nzp,nxp,nyp), rv(nzp,nxp,nyp), rc(nzp,nxp,nyp)

    IF (level<4) THEN
       rv = a_rv ! Water vapor
       rc = a_rc + a_rpp ! Total condensate (cloud + precipitation)
    ELSE
       rv = a_rp ! Water vapor
       rc = a_rc + a_srp + a_ri + a_srs ! Total condensed water (aerosol+cloud+precipitation+ice+snow)
    END IF
    call boyanc(nzp,nxp,nyp,a_wt,a_theta,rv,th00,a_tmp1,rc)

    call ae1mm(nzp,nxp,nyp,a_wt,awtbar)
    call update_pi1(nzp,awtbar,pi1)

    if (sflg)  call comp_tke(nzp,nxp,nyp,dzm,th00,a_uc,a_vc,a_wc,a_tmp1)

  end subroutine buoyancy
  !
  ! ----------------------------------------------------------------------
  ! subroutine boyanc:
  !
  subroutine boyanc(n1,n2,n3,wt,th,rv,th00,scr,rc)

    use defs, only: g, ep2

    integer, intent(in) :: n1,n2,n3
    real, intent(in)    :: th00,th(n1,n2,n3),  &
                           rv(n1,n2,n3), &  ! Total water vapour mixing ratio
                           rc(n1,n2,n3)     ! Total condensed water (aerosol, cloud, rain, ice and snow) mixing ratio
    real, intent(inout) :: wt(n1,n2,n3)
    real, intent(out)   :: scr(n1,n2,n3)

    integer :: k, i, j
    real :: gover2

    gover2  = 0.5*g

    do j=3,n3-2
       do i=3,n2-2
          do k=1,n1
             scr(k,i,j)=gover2*((th(k,i,j)*(1.+ep2*rv(k,i,j))-th00)/th00-rc(k,i,j))
          end do

          do k=2,n1-2
             wt(k,i,j)=wt(k,i,j)+scr(k,i,j)+scr(k+1,i,j)
          end do
       end do
    end do

  end subroutine boyanc
  !
  ! ----------------------------------------------------------------------
  ! subroutine corlos:  This is the coriolis driver, its purpose is to
  ! from the coriolis accelerations for u and v and add them into the
  ! accumulated tendency arrays of ut and vt.
  !
  subroutine corlos

    use defs, only : omega
    use grid, only : a_uc, a_vc, a_ut, a_vt, nxp, nyp, nzp, u0, v0

    logical, save :: initialized = .False.
    real, save    :: fcor

    integer :: i, j, k

    if (corflg) then
       if (.not.initialized) fcor=2.*omega*sin(cntlat*0.01745329)
       do j=3,nyp-2
          do i=3,nxp-2
             do k=2,nzp
                a_ut(k,i,j)=a_ut(k,i,j) - fcor*(v0(k)-0.25*                   &
                     (a_vc(k,i,j)+a_vc(k,i+1,j)+a_vc(k,i,j-1)+a_vc(k,i+1,j-1)))
                a_vt(k,i,j)=a_vt(k,i,j) + fcor*(u0(k)-0.25*                   &
                     (a_uc(k,i,j)+a_uc(k,i-1,j)+a_uc(k,i,j+1)+a_uc(k,i-1,j+1)))
             end do
          end do
       end do
       initialized = .True.
    end if

  end subroutine corlos
!
! ----------------------------------------------------------------------
! subroutine sponge: does the rayleigh friction for the momentum terms,
! and newtonian damping of thermal term the damping is accumulated with the
! other tendencies
!
  subroutine sponge (isponge)

    use grid, only : u0, v0, a_up, a_vp, a_wp, a_tp, a_ut, a_vt, a_wt, a_tt,&
         nfpt, spng_tfct, spng_wfct, nzp, nxp, nyp, th0, th00

    integer, intent (in) :: isponge

    integer :: i, j, k, kk

    if (maxval(spng_tfct) > epsilon(1.) .and. nfpt > 1) then
       do j=3,nyp-2
          do i=3,nxp-2
             do k=nzp-nfpt,nzp-1
                kk = k+1-(nzp-nfpt)
                if (isponge == 0) then
                   a_tt(k,i,j)=a_tt(k,i,j) - spng_tfct(kk)*                   &
                        (a_tp(k,i,j)-th0(k)+th00)
                else
                   a_ut(k,i,j)=a_ut(k,i,j) - spng_tfct(kk)*(a_up(k,i,j)-u0(k))
                   a_vt(k,i,j)=a_vt(k,i,j) - spng_tfct(kk)*(a_vp(k,i,j)-v0(k))
                   a_wt(k,i,j)=a_wt(k,i,j) - spng_wfct(kk)*(a_wp(k,i,j))
                end if
             end do
          end do
       end do
    end if

  end subroutine sponge

  !
  ! ---------------------------------------------------------------------
  ! SALSA_diagnostics: Update properties for the current timestep:
  !                    E.g. if enough water has evaporated from droplets,
  !                    deplete the cloud droplet bins and move CCN material
  !                    back to the aerosol regime.
  !                    In addition, update the diagnostic scalars for total grid-cell
  !                    liquid water contents.
  !
  ! Juha Tonttila, FMI, 2014
  ! Tomi Raatikainen, FMI, 2016

  SUBROUTINE SALSA_diagnostics(reset_stats)
    USE grid, ONLY : nxp,nyp,nzp,    &
                     a_naerop,a_maerop,a_ncloudp,a_mcloudp,a_nprecpp,a_mprecpp,      &
                     a_gaerop, prtcl, a_rh, a_temp, a_rhi, a_dn,                     &
                     a_nicep,a_micep,a_nsnowp,a_msnowp, diss, mws, dens, dens_ice, dens_snow, level, &
                     dtl, diag_ra, diag_na, diag_rc, diag_nc, diag_rr, diag_nr, diag_ri, diag_ni, diag_rs, diag_ns
    USE mo_submctl, ONLY : nbins,ncld,nprc,ica,fca,icb,fcb,ira,fra,in2a,fn2a,    &
                               nice,nsnw,iia,fia,iib,fib,isa,fsa,        &
                               msu,moc,mno,mnh,avog,pi6,                     &
                               surfw0, rg, nlim, prlim, pi, &
                               lscndgas, aerobins, calc_correlation
    USE class_ComponentIndex, ONLY : GetIndex, GetNcomp, IsUsed

    IMPLICIT NONE

    LOGICAL :: reset_stats

    INTEGER :: i,j,k,bc,ba,bb,s,sc,sa,nc,nn

    REAL :: zvol, ra, rb
    REAL :: ns, cd
    REAL, DIMENSION(nzp,nxp,nyp) :: tmp_ra, tmp_na, tmp_rc, tmp_nc, tmp_rr, tmp_nr, tmp_ri, tmp_ni, tmp_rs, tmp_ns

    nn = GetNcomp(prtcl)+1 ! total number of species

    ! Change in concentrations due to diagnostics (e.g. release of cloud/rain/ice/snow into aerosol,
    ! cleaning particles without mass or negligible number concentration)
    tmp_ra(:,:,:)=SUM(a_maerop(:,:,:,(nn-1)*nbins+1:nn*nbins),DIM=4)
    tmp_na(:,:,:)=SUM(a_naerop,DIM=4)
    tmp_rc(:,:,:)=SUM(a_mcloudp(:,:,:,(nn-1)*ncld+1:nn*ncld),DIM=4)
    tmp_nc(:,:,:)=SUM(a_ncloudp,DIM=4)
    tmp_rr(:,:,:)=SUM(a_mprecpp(:,:,:,(nn-1)*nprc+1:nn*nprc),DIM=4)
    tmp_nr(:,:,:)=SUM(a_nprecpp,DIM=4)
    tmp_ri(:,:,:)=SUM(a_micep(:,:,:,(nn-1)*nice+1:nn*nice),DIM=4)
    tmp_ni(:,:,:)=SUM(a_nicep,DIM=4)
    tmp_rs(:,:,:)=SUM(a_msnowp(:,:,:,(nn-1)*nsnw+1:nn*nsnw),DIM=4)
    tmp_ns(:,:,:)=SUM(a_nsnowp,DIM=4)
    IF (reset_stats) THEN ! Reset outputs (there are two SALSA_diagnostics calls during one time step)
        diag_ra=0.; diag_na=0.
        diag_rc=0.; diag_nc=0.
        diag_rr=0.; diag_nr=0.
        diag_ri=0.; diag_ni=0.
        diag_rs=0.; diag_ns=0.
    ENDIF


    ! Remove negative values
    a_naerop = MAX(0.,a_naerop)
    a_ncloudp = MAX(0.,a_ncloudp)
    a_nprecpp = MAX(0.,a_nprecpp)
    a_maerop = MAX(0.,a_maerop)
    a_mcloudp = MAX(0.,a_mcloudp)
    a_mprecpp = MAX(0.,a_mprecpp)

    a_nicep = MAX(0.,a_nicep)
    a_nsnowp = MAX(0.,a_nsnowp)
    a_micep = MAX(0.,a_micep)
    a_msnowp = MAX(0.,a_msnowp)

    ! Remove particles that have number but no mass. Also remove particles that have
    ! insignificant concentration indicated by nlim and prlim (note: #/m^3)
    DO j = 3,nyp-2
       DO i = 3,nxp-2
          DO k = 1,nzp
             ! Aerosols
             DO bc = 1,nbins
                IF (a_naerop(k,i,j,bc) > 0. .AND. SUM(a_maerop(k,i,j,bc:(nn-1)*nbins+bc:nbins)) <= 0. .OR. &
                        a_naerop(k,i,j,bc)*a_dn(k,i,j) < nlim) THEN
                   a_naerop(k,i,j,bc) = 0.
                   a_maerop(k,i,j,bc:(nn-1)*nbins+bc:nbins) = 0.
                END IF
             END DO

             ! Clouds
             DO bc = 1,ncld
                IF (a_ncloudp(k,i,j,bc) > 0. .AND. SUM(a_mcloudp(k,i,j,bc:(nn-1)*ncld+bc:ncld)) <= 0. .OR. &
                        a_ncloudp(k,i,j,bc)*a_dn(k,i,j) < nlim) THEN
                   a_ncloudp(k,i,j,bc) = 0.
                   a_mcloudp(k,i,j,bc:(nn-1)*ncld+bc:ncld) = 0.
                END IF
             END DO ! ncld

             ! Precipitation
             DO bc = 1,nprc
                IF (a_nprecpp(k,i,j,bc) > 0. .AND. a_mprecpp(k,i,j,(nn-1)*nprc+bc) <= 0. .OR. &
                        a_nprecpp(k,i,j,bc)*a_dn(k,i,j) < prlim) THEN
                   a_nprecpp(k,i,j,bc) = 0.
                   a_mprecpp(k,i,j,bc:(nn-1)*nprc+bc:nprc) = 0.
                END IF
             END DO ! nprc

             ! Ice
             IF (level<5) CYCLE
             DO bc = 1,nice
                IF (a_nicep(k,i,j,bc) > 0. .AND. SUM(a_micep(k,i,j,bc:(nn-1)*nice+bc:nice)) <= 0. .OR. &
                        a_nicep(k,i,j,bc)*a_dn(k,i,j) < prlim ) THEN
                   a_nicep(k,i,j,bc) = 0.
                   a_micep(k,i,j,bc:(nn-1)*nice+bc:nice) = 0.
                END IF
             END DO ! ncld

             ! Snow
             DO bc = 1,nsnw
                IF (a_nsnowp(k,i,j,bc) > 0. .AND. a_msnowp(k,i,j,(nn-1)*nsnw+bc) <= 0. .OR. &
                        a_nsnowp(k,i,j,bc)*a_dn(k,i,j) < prlim) THEN
                   a_nsnowp(k,i,j,bc) = 0.
                   a_msnowp(k,i,j,bc:(nn-1)*nsnw+bc:nsnw) = 0.
                END IF
             END DO ! nsnw

          END DO !k
       END DO !i
    END DO !j

    ! Ghost species
    DO j = 3,nyp-2
       DO i = 3,nxp-2
          DO k = 1,nzp
             ! Loop over cloud droplet bins
             DO bc = ica%cur,fcb%cur

                IF ( a_ncloudp(k,i,j,bc)*a_dn(k,i,j) > nlim .AND. a_rh(k,i,j)<0.999 .AND. &
                        a_mcloudp(k,i,j,(nn-1)*ncld+bc)<1e-5 ) THEN
                   ! Critical diameter (assuming soluble CCN)
                   ns = SUM( diss(1:nn-1)*a_mcloudp(k,i,j,bc:(nn-2)*ncld+bc:ncld)/mws(1:nn-1) )/a_ncloudp(k,i,j,bc)
                   cd = 3.*SQRT(ns*rg*a_temp(k,i,j)/(2.*pi*surfw0))

                   ! Wet diameter
                   zvol = (SUM( a_mcloudp(k,i,j,bc:(nn-1)*ncld+bc:ncld)/dens(1:nn) )/a_ncloudp(k,i,j,bc)/pi6)**(1./3.)

                   ! Lose the droplets if smaller than 0.2*critical diameter or 2 um or if there is no water
                   IF ( zvol < MAX(0.2*cd,2.e-6) .OR. a_mcloudp(k,i,j,(nn-1)*ncld+bc)<1e-25*a_ncloudp(k,i,j,bc) ) THEN
                      IF (bc<=fca%cur) THEN
                          ba = ica%par + (bc-ica%cur) ! Index for parallel aerosol bin
                      ELSE
                          ba = icb%par + (bc-icb%cur) ! Index for parallel aerosol bin
                      ENDIF
                      ! Move the number of particles from cloud to aerosol bins
                      a_naerop(k,i,j,ba) = a_naerop(k,i,j,ba) + a_ncloudp(k,i,j,bc)
                      a_ncloudp(k,i,j,bc) = 0.

                      ! Move ccn material back to aerosol regime (including water)
                      DO s = 1,nn
                         sc = (s-1)*ncld + bc
                         sa = (s-1)*nbins + ba
                         a_maerop(k,i,j,sa) = a_maerop(k,i,j,sa) + a_mcloudp(k,i,j,sc)
                         a_mcloudp(k,i,j,sc) = 0.
                      END DO
                   END IF ! critical diameter

                END IF  ! blim

             END DO ! bc

             ! Loop over precipitation bins
             DO bc = ira,fra

                IF ( a_nprecpp(k,i,j,bc)*a_dn(k,i,j) > prlim .AND. a_rh(k,i,j)<0.999 .AND. &
                        a_mprecpp(k,i,j,(nn-1)*nprc+bc)<1e-6 ) THEN
                   ! Critical diameter
                   ns = SUM( diss(1:nn-1)*a_mprecpp(k,i,j,bc:(nn-2)*nprc+bc:nprc)/mws(1:nn-1) )/a_nprecpp(k,i,j,bc)
                   cd = 3.*SQRT(ns*rg*a_temp(k,i,j)/(2.*pi*surfw0))

                   ! Wet diameter
                   zvol = (SUM( a_mprecpp(k,i,j,bc:(nn-1)*nprc+bc:nprc)/dens(1:nn) )/a_nprecpp(k,i,j,bc)/pi6)**(1./3.)

                   ! Lose the droplets if smaller than 0.02*critical diameter or 2 um or if there is no water
                   IF ( zvol < MAX(0.02*cd,2.e-6) .OR. a_mprecpp(k,i,j,(nn-1)*nprc+bc)<1e-25*a_nprecpp(k,i,j,bc) ) THEN

                      ! Move evaporating precipitation to aerosol bin based on dry radius and chemical composition

                      ! 1) Find the closest matching bin based on dry particle radius (a and b bins)
                      cd = 0.5*(SUM( a_mprecpp(k,i,j,bc:(nn-2)*nprc+bc:nprc)/dens(1:nn-1) )/a_nprecpp(k,i,j,bc)/pi6)**(1./3.) ! Dry radius
                      ba=fn2a ! Ignore 1a and note that aerobins contains the lower limit of bin dry radius
                      DO WHILE (cd<aerobins(ba) .AND. ba>in2a)
                         ba=ba-1
                      ENDDO
                      ! Corresponding b bin is ba+(fn2a-fn1a)=ba+fn2a-(in2a-1)=ba+fn2a-in2a+1
                      bb=ba+fn2a-in2a+1
                      ! 2) Select a or b bin
                      IF (a_naerop(k,i,j,bb)*a_dn(k,i,j)<=nlim) THEN
                         ! Empty b bin so select a
                         !ba = ba
                      ELSEIF (a_naerop(k,i,j,ba)*a_dn(k,i,j)<=nlim) THEN
                         ! Empty a bin so select b
                         ba = bb
                      ELSE
                         ! Both are present - find bin based on compositional similarity
                         ra=calc_correlation(a_maerop(k,i,j,ba:(nn-2)*nbins+ba:nbins),a_mprecpp(k,i,j,bc:(nn-2)*nprc+bc:nprc),nn-1)
                         rb=calc_correlation(a_maerop(k,i,j,bb:(nn-2)*nbins+bb:nbins),a_mprecpp(k,i,j,bc:(nn-2)*nprc+bc:nprc),nn-1)
                         IF (ra<rb) ba = bb
                      ENDIF

                      ! Move the number of particles from precipitation to aerosol bins
                      a_naerop(k,i,j,ba) = a_naerop(k,i,j,ba) + a_nprecpp(k,i,j,bc)
                      a_nprecpp(k,i,j,bc) = 0.

                      ! Move ccn material back to aerosol regime (including water)
                      DO s = 1,nn
                         sc = (s-1)*nprc + bc
                         sa = (s-1)*nbins + ba
                         a_maerop(k,i,j,sa) = a_maerop(k,i,j,sa) + a_mprecpp(k,i,j,sc)
                         a_mprecpp(k,i,j,sc) = 0.
                      END DO

                   END IF ! Critical diameter

                END IF ! prlim

             END DO ! bc

             ! Loop over ice bins
             DO bc = iia%cur,fib%cur

                IF ( a_nicep(k,i,j,bc)*a_dn(k,i,j) > prlim .AND. a_rhi(k,i,j)<0.999 .AND. &
                        a_micep(k,i,j,(nn-1)*nice+bc)<1e-8 ) THEN
                   ! Diameter (assuming constant ice density)
                   cd = (SUM( a_micep(k,i,j,bc:(nn-1)*nice+bc:nice)/dens_ice(1:nn) )/a_nicep(k,i,j,bc)/pi6)**(1./3.)

                   ! Dry to total mass ratio
                   zvol = SUM( a_micep(k,i,j,bc:(nn-2)*nice+bc:nice) )/SUM( a_micep(k,i,j,bc:(nn-1)*nice+bc:nice) )

                   ! Ice and snow don't have a critical size, but lose particles smaller than 2e-6 m and particles which dry to total mass ratio is more than 0.5
                   IF ( zvol>0.5 .OR. cd<2e-6 ) THEN
                      IF (bc<=fia%cur) THEN
                         ba = iia%par + (bc-iia%cur) ! Index for parallel aerosol bin
                      ELSE
                         ba = iib%par + (bc-iib%cur) ! Index for parallel aerosol bin
                      ENDIF
                      ! Move the number of particles from ice to aerosol bins
                      a_naerop(k,i,j,ba) = a_naerop(k,i,j,ba) + a_nicep(k,i,j,bc)
                      a_nicep(k,i,j,bc) = 0.

                      ! Move mass to aerosol (including water)
                      DO s = 1,nn
                         sc = (s-1)*nice + bc
                         sa = (s-1)*nbins + ba
                         a_maerop(k,i,j,sa) = a_maerop(k,i,j,sa) + a_micep(k,i,j,sc)
                         a_micep(k,i,j,sc) = 0.
                      END DO
                   END IF

                END IF  ! prlim

             END DO ! bc

             ! Loop over snow bins
             DO bc = isa,fsa

                IF ( a_nsnowp(k,i,j,bc)*a_dn(k,i,j) > prlim .AND. a_rhi(k,i,j)<0.999 .AND. &
                        a_msnowp(k,i,j,(nn-1)*nsnw+bc)<1e-8 ) THEN
                   ! Diameter (assuming constant snow density)
                   cd = (SUM( a_msnowp(k,i,j,bc:(nn-1)*nsnw+bc:nsnw)/dens_snow(1:nn) )/a_nsnowp(k,i,j,bc)/pi6)**(1./3.)

                   ! Dry to total mass ratio
                   zvol = SUM( a_msnowp(k,i,j,bc:(nn-2)*nsnw+bc:nsnw) )/SUM( a_msnowp(k,i,j,bc:(nn-1)*nsnw+bc:nsnw) )

                   ! Lose particles smaller than 2e-6 m and particles which dry to total mass ratio is more than 0.5
                   IF ( zvol>0.5  .OR. cd<2.e-6 ) THEN

                      ! Move evaporating snow to aerosol bin based on dry radius and chemical composition

                      ! 1) Find the closest matching bin based on dry particle radius (a and b bins)
                      cd = 0.5*(SUM( a_msnowp(k,i,j,bc:(nn-2)*nsnw+bc:nsnw)/dens(1:nn-1) )/a_nsnowp(k,i,j,bc)/pi6)**(1./3.) ! Dry radius
                      ba=fn2a ! Ignore 1a and note that aerobins contains the lower limit of bin dry radius
                      DO WHILE (cd<aerobins(ba) .AND. ba>in2a)
                         ba=ba-1
                      ENDDO
                      ! Corresponding b bin is ba+(fn2a-fn1a)=ba+fn2a-(in2a-1)=ba+fn2a-in2a+1
                      bb=ba+fn2a-in2a+1
                      ! 2) Select a or b bin
                      IF (a_naerop(k,i,j,bb)*a_dn(k,i,j)<=nlim) THEN
                         ! Empty b bin so select a
                         !ba = ba
                      ELSEIF (a_naerop(k,i,j,ba)*a_dn(k,i,j)<=nlim) THEN
                         ! Empty a bin so select b
                         ba = bb
                      ELSE
                         ! Both are present - find bin based on compositional similarity
                         ra=calc_correlation(a_maerop(k,i,j,ba:(nn-2)*nbins+ba:nbins),a_msnowp(k,i,j,bc:(nn-2)*nsnw+bc:nsnw),nn-1)
                         rb=calc_correlation(a_maerop(k,i,j,bb:(nn-2)*nbins+bb:nbins),a_msnowp(k,i,j,bc:(nn-2)*nsnw+bc:nsnw),nn-1)
                         IF (ra<rb) ba = bb
                      ENDIF

                      ! Move the number of particles from snow to aerosol bins
                      a_naerop(k,i,j,ba) = a_naerop(k,i,j,ba) + a_nsnowp(k,i,j,bc)
                      a_nsnowp(k,i,j,bc) = 0.

                      ! Move mass to aerosol (including water)
                      DO s = 1,nn
                         sc = (s-1)*nsnw + bc
                         sa = (s-1)*nbins + ba
                         a_maerop(k,i,j,sa) = a_maerop(k,i,j,sa) + a_msnowp(k,i,j,sc)
                         a_msnowp(k,i,j,sc) = 0.
                      END DO
                   END IF

                END IF ! prlim

             END DO ! bc

             ! Loop over aerosol bins
             DO ba = 1,nbins
                IF (a_naerop(k,i,j,ba)*a_dn(k,i,j) > nlim) THEN
                   zvol = SUM( a_maerop(k,i,j,ba:(nn-2)*nbins+ba:nbins)/dens(1:nn-1) )/a_naerop(k,i,j,ba) ! Dry volume

                   ! Particles smaller than 0.1 nm diameter are set to zero
                   IF ( zvol < pi6*1.e-10**3 ) THEN
                      ! Volatile species to the gas phase
                      IF (IsUsed(prtcl,'SO4') .AND. lscndgas) THEN
                         nc = GetIndex(prtcl,'SO4')
                         s = (nc-1)*nbins + ba
                         a_gaerop(k,i,j,1) = a_gaerop(k,i,j,1) + a_maerop(k,i,j,s) / msu * avog
                      END IF
                      IF (IsUsed(prtcl,'OC') .AND. lscndgas) THEN
                         nc = GetIndex(prtcl,'OC')
                         s = (nc-1)*nbins + ba
                         a_gaerop(k,i,j,5) = a_gaerop(k,i,j,5) + a_maerop(k,i,j,s) / moc * avog
                      END IF
                      IF (IsUsed(prtcl,'NO') .AND. lscndgas) THEN
                         nc = GetIndex(prtcl,'NO')
                         s = (nc-1)*nbins + ba
                         a_gaerop(k,i,j,2) = a_gaerop(k,i,j,2) + a_maerop(k,i,j,s) / mno * avog
                      END IF
                      IF (IsUsed(prtcl,'NH') .AND. lscndgas) THEN
                         nc = GetIndex(prtcl,'NH')
                         s = (nc-1)*nbins + ba
                         a_gaerop(k,i,j,3) = a_gaerop(k,i,j,3) + a_maerop(k,i,j,s) / mnh * avog
                      END IF

                      ! Mass and number to zero (insolube species and water are lost)
                      a_maerop(k,i,j,ba:(nn-1)*nbins+ba:nbins) = 0.
                      a_naerop(k,i,j,ba) = 0.
                   END IF
                END IF
             END DO

          END DO   ! k
       END DO   ! i
    END DO   ! j

    diag_ra(:,:,:)=diag_ra(:,:,:)+( SUM(a_maerop(:,:,:,(nn-1)*nbins+1:nn*nbins),DIM=4)-tmp_ra(:,:,:) )/dtl
    diag_na(:,:,:)=diag_na(:,:,:)+( SUM(a_naerop,DIM=4)-tmp_na(:,:,:) )/dtl
    diag_rc(:,:,:)=diag_rc(:,:,:)+( SUM(a_mcloudp(:,:,:,(nn-1)*ncld+1:nn*ncld),DIM=4)-tmp_rc(:,:,:) )/dtl
    diag_nc(:,:,:)=diag_nc(:,:,:)+( SUM(a_ncloudp,DIM=4)-tmp_nc(:,:,:) )/dtl
    diag_rr(:,:,:)=diag_rr(:,:,:)+( SUM(a_mprecpp(:,:,:,(nn-1)*nprc+1:nn*nprc),DIM=4)-tmp_rr(:,:,:) )/dtl
    diag_nr(:,:,:)=diag_nr(:,:,:)+( SUM(a_nprecpp,DIM=4)-tmp_nr(:,:,:) )/dtl
    diag_ri(:,:,:)=diag_ri(:,:,:)+( SUM(a_micep(:,:,:,(nn-1)*nice+1:nn*nice),DIM=4)-tmp_ri(:,:,:) )/dtl
    diag_ni(:,:,:)=diag_ni(:,:,:)+( SUM(a_nicep,DIM=4)-tmp_ni(:,:,:) )/dtl
    diag_rs(:,:,:)=diag_rs(:,:,:)+( SUM(a_msnowp(:,:,:,(nn-1)*nsnw+1:nn*nsnw),DIM=4)-tmp_rs(:,:,:) )/dtl
    diag_ns(:,:,:)=diag_ns(:,:,:)+( SUM(a_nsnowp,DIM=4)-tmp_ns(:,:,:) )/dtl

    !!!!!!!!!!!!!!!!!!!!!!!
    ! Update diagnostic tracers
    !!!!!!!!!!!!!!!!!!!!!!!

    CALL SALSA_diag_update

  END SUBROUTINE SALSA_diagnostics

  !
  ! ---------------------------------------------------------------------
  ! SALSA_diag_update: Update diagnostic concentration tracers

  SUBROUTINE SALSA_diag_update
    USE class_ComponentIndex, ONLY : GetIndex
    USE grid, ONLY : a_maerop,a_mcloudp,a_mprecpp,a_nprecpp, &
            a_micep,a_msnowp,a_nsnowp, a_rc, a_srp,a_snrp, a_ri,a_srs,a_snrs,prtcl
    USE mo_submctl, ONLY : nbins,in1a,fn2b,ncld,ica,fcb,nprc,ira,fra,nice,iia,fib,nsnw,isa,fsa
    INTEGER :: nc, str, end

    ! Liquid water content
    nc = GetIndex(prtcl,'H2O')
    ! Aerosols, regimes a and b
    str = (nc-1)*nbins + in1a
    end = (nc-1)*nbins + fn2b
    a_rc(:,:,:) = SUM(a_maerop(:,:,:,str:end),DIM=4)
    ! Clouds, regime a and b
    str = (nc-1)*ncld+ica%cur
    end = (nc-1)*ncld+fcb%cur
    a_rc(:,:,:) = a_rc(:,:,:) + SUM(a_mcloudp(:,:,:,str:end),DIM=4)
    ! Precipitation
    str = (nc-1)*nprc+ira
    end = (nc-1)*nprc+fra
    a_srp(:,:,:) = SUM(a_mprecpp(:,:,:,str:end),DIM=4)
    a_snrp(:,:,:) = SUM(a_nprecpp(:,:,:,ira:fra),DIM=4)

    ! ice, regimes a and b
    str = (nc-1)*nice+iia%cur
    end = (nc-1)*nice+fib%cur
    a_ri(:,:,:) = SUM(a_micep(:,:,:,str:end),DIM=4)
    ! Snow
    str = (nc-1)*nsnw+isa
    end = (nc-1)*nsnw+fsa
    a_srs(:,:,:) = SUM(a_msnowp(:,:,:,str:end),DIM=4)
    a_snrs(:,:,:) = SUM(a_nsnowp(:,:,:,isa:fsa),DIM=4)

  END SUBROUTINE SALSA_diag_update


end module step
