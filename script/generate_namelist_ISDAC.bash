#!/bin/bash

# Exit on error
set -e


##########################
###                    ###
### write the namelist ###
###                    ###
##########################

cat > ${dir}/NAMELIST <<EOF
 &version
  ver="v1.0.4"
 /

 &model
  nxp =   ${nxp:-68}
  nyp =   ${nyp:-68}
  nzp =   ${nzp:-140}
  deltax = ${deltax:-50.}
  deltay = ${deltay:-50.}
  deltaz = ${deltaz:-10.}
  nxpart = ${nxpart:-.false.}
  dzmax  = ${dzmax:-1200.}
  dzrat  = ${dzrat:-1.05}
  dtlong = ${dtlong:-1.}
  distim = ${distim:-100.}
  timmax = ${timmax:-28800.}
  Tspinup = ${Tspinup:-7200.}
${notJJA}  minispinup01 = ${minispinup01:-0.}
${notJJA}  minispinup02 = ${minispinup02:-0.}
${notJJA}  minispinupCase01 = ${minispinupCase01:-3}
${notJJA}  minispinupCase02 = ${minispinupCase02:-3}
  runtype = ${runtype:-'"INITIAL"'}
  level = ${level:-5}
  CCN = ${CCN:-30.e6}
  prndtl = ${prndtl:--0.3333333}
  filprf = ${filprf:-"'isdac'"}
  hfilin = ${hfilin:-"'isdac.rst'"}
  ssam_intvl = ${ssam_intvl:-120.}
  savg_intvl = ${savg_intvl:-120.}
  mcflg = ${mcflg:-.FALSE.}
  frqhis  = ${frqhis:-3600.}
  istpfl  = ${istpfl:-1}
  lbinanl = ${lbinanl:-.false.}
  frqanl = ${frqanl:-5400.}
  corflg = ${corflg:-.false.}
  ipsflg = ${ipsflg:-1}
  itsflg = ${itsflg:-1}
  strtim = ${strtim:-180.0}
  sed_aero = ${sed_aero:-.FALSE.}
  sed_cloud = ${sed_cloud:-.TRUE.}
  sed_precp = ${sed_precp:-.TRUE.}
  sed_ice = ${sed_ice:-.FALSE.}
  sed_snow = ${sed_snow:-.FALSE.}
  iradtyp = ${iradtyp:-3}                ! 1 = no radiation, only large-scale forcing, 3 = radiation + large-scale forcing 
  case_name = ${case_name:-"'ascos'"}            ! Case-specific large-scale forcing: none = not used, 
                                      ! default = simple divergence forcing with specified div 
  div = ${div:-1.5e-6}              ! Divergence for e.g. case_name = 'default'
  sfc_albedo = ${sfc_albedo:-0.7}
  radsounding = ${radsounding:-"'datafiles/ksaw.lay'"}
  
  cntlat = ${cntlat:-71.32}
  strtim = ${strtim:-117.75}
  

  isfctyp = ${isfctyp:-0} ! surface fluxes
  sst = ${sst:-267.}

  dthcon = ${dthcon:-0.} ! heat flux
  drtcon = ${drtcon:-0.}  ! latent

  ubmin  = ${ubmin:--0.25}
  zrough = ${zrough:-4e-4}
  th00 = ${th00:-267.}
  umean =  ${umean:--7.0}
  vmean = ${vmean:-2.45554452055}
 /

 &salsa	
   nlcoag = ${nlcoag:-.FALSE.}       ! Master coagulation switch
   
   !! selfcoagulation processes   
   nlcgcc = ${nlcgcc:-T}       ! Self-collection of cloud droplets
   nlcgpp = ${nlcgpp:-T}       ! Self-collection of rain drops
   nlcgaa = ${nlcgaa:-.FALSE.}      ! Aerosol coagulation
   nlcgii = ${nlcgii:-T}       ! Self-collection of ice
   nlcgss = ${nlcgss:-.FALSE.}       ! Self-collection of snow

   !! coagulation between different particles   
   nlcgpc = ${nlcgpc:-T}       ! Rain collection of cloud droplets
   nlcgca = ${nlcgca:-T}       ! Cloud collection of aerosols
   nlcgpa = ${nlcgpa:-T}       ! Rain collection of aerosols

   ! ice related
   nlcgia = ${nlcgia:-T}       ! Ice collection of aerosols
   nlcgic = ${nlcgic:-T}       ! Ice collection of cloud droplets
   nlcgip = ${nlcgip:-T}       ! Ice collection of rain drops
   
   ! snow related
   nlcgsa = ${nlcgsa:-.FALSE.}       ! Snow collection of aerosols
   nlcgsc = ${nlcgsc:-.FALSE.}       ! Snow collection of cloud droplets
   nlcgsi = ${nlcgsi:-.FALSE.}       ! Snow collection of ice particles
   nlcgsp = ${nlcgsp:-.FALSE.}       ! Snow collection of rain drops

   nlcnd       = ${nlcnd:-.TRUE.}  ! Master condensation switch
   nlcndgas    = ${nlcndgas:-.FALSE.}  ! --Aerosol precursor gas codensation
   nlcndh2oae  = ${nlcndh2oae:-.TRUE.}  ! --Condensation of water on aerosols (if FALSE, equilibrium assumed)
   nlcndh2ocl  = ${nlcndh2ocl:-.TRUE.}  ! --Condensation of water on cloud droplets (and drizzle)
   nlcndh2oic  = ${nlcndh2oic:-.TRUE.}  ! --Condensation of water on ice particles
   nlauto      = ${nlauto:-.TRUE.}  ! Master autoconversion switch
   nlautosnow  = ${nlautosnow:-.FALSE.} ! Master snow autoconversion switch
   nlactiv     = ${nlactiv:-.TRUE.}  ! Master cloud activation switch
   nlactbase   = ${nlactbase:-.FALSE.}  ! --Switch for parameterized cloud base activation
   nlactintst  = ${nlactintst:-.TRUE.}  ! --Switch for interstitial activation based on host model Smax

   nlichom     = ${nlichom:-.FALSE.}     ! Switch for homogeneous ice nucleation
   nlichet     = ${nlichet:-.FALSE.}     ! Switch for heterogeneous ice nucleation
   nlicimmers  = ${nlicimmers:-.FALSE.}   ! Switch for ice nucleation by immersion
   nlicmelt    = ${nlicmelt:-.FALSE.}    ! Switch for ice'n' snow melting
${notJJA}   nlicbasic   = ${nlicbasic:-.FALSE.}
   
${notJJA}   nlfixinc   = ${nlfixinc:-.TRUE.}      ! Fix ice number concentration to be over given limit fixINC
${notJJA}   fixINC     = ${fixINC:-1.0}         ! fixed ice number concentration #/kg, nlfixinc should be set to true inorder to have this working

   rhlim = ${rhlim:-1.2}          ! RH limit for SALSA during initialization and spinup

   isdtyp = ${isdtyp:-0}
   nspec = ${nspec:-1}
   listspec = ${listspec:-"'SO4','','','','','',''"}            !!!! "'SO4','DU','OC','','','',''"
   volDistA = ${volDistA:-1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}   
   volDistB = ${volDistB:-0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0}
   nf2a = ${nf2a:-1.0}

   sigmag = ${sigmag:- 1.5,    2.45, 2.0, 2.0, 2.0, 2.0, 2.0}  ! Stdev for initial aerosol size distribution for isdtyp == 0 (uniform)  
   dpg    = ${dpg:-    0.2,     0.7, 0.2, 0.2, 0.2, 0.2, 0.2}     ! Mode mean diameters in micrometers
   n      = ${n:-   156.42,    6.42,  0.,  0.,  0.,  0.,  0.}  ! Mode number concentrations in #/mg
 /

EOF
 
exit