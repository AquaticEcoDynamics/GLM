!###############################################################################
!#                                                                             #
!# glm_api_aed.F90                                                             #
!#                                                                             #
!# The interface between glm and libaed-xxx                                    #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Agriculture and Environment                                   #
!#     The University of Western Australia                                     #
!#                                                                             #
!#     http://aquatic.science.uwa.edu.au/                                      #
!#                                                                             #
!# Copyright 2024 - The University of Western Australia                        #
!#                                                                             #
!#  This file is part of GLM (General Lake Model)                              #
!#                                                                             #
!#  GLM is free software: you can redistribute it and/or modify                #
!#  it under the terms of the GNU General Public License as published by       #
!#  the Free Software Foundation, either version 3 of the License, or          #
!#  (at your option) any later version.                                        #
!#                                                                             #
!#  GLM is distributed in the hope that it will be useful,                     #
!#  but WITHOUT ANY WARRANTY; without even the implied warranty of             #
!#  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the              #
!#  GNU General Public License for more details.                               #
!#                                                                             #
!#  You should have received a copy of the GNU General Public License          #
!#  along with this program.  If not, see <http://www.gnu.org/licenses/>.      #
!#                                                                             #
!###############################################################################

#include "aed_api.h"

#undef MISVAL
#ifndef _FORTRAN_SOURCE_
#define _FORTRAN_SOURCE_ 1
#endif

#ifndef __GFORTRAN__
#  ifndef isnan
#    define isnan(x) ieee_is_nan(x)
#  endif
#endif

#include "glm.h"

!-------------------------------------------------------------------------------
MODULE glm_api_aed
!
   USE ISO_C_BINDING
   USE IEEE_ARITHMETIC

   USE aed_util
   USE aed_common
   USE glm_types
   USE aed_api
   USE aed_zones
   USE glm_api_zones


   IMPLICIT NONE

   PRIVATE ! By default, make everything private
!
#include "glm_globals.h"
#include "glm_plot.h"
#include "glm_ncdf.h"
#include "glm_csv.h"
#include "glm_mobl.h"
!
#if USE_DL_LOADER
# define _WQ_INIT_GLM        "wq_init_glm"
# define _WQ_SET_GLM_DATA    "wq_set_glm_data"
# define _WQ_DO_GLM          "wq_do_glm"
# define _WQ_CLEAN_GLM       "wq_clean_glm"
# define _WQ_INIT_GLM_OUTPUT "wq_init_glm_output"
# define _WQ_WRITE_GLM_      "wq_write_glm"
# define _WQ_VAR_INDEX_C     "wq_var_index_c"
# define _WQ_IS_VAR          "wq_is_var"
# define _WQ_INFLOW_UPDATE   "wq_inflow_update"
#else
# define _WQ_INIT_GLM        "api_init_glm"
# define _WQ_SET_GLM_DATA    "api_set_glm_data"
# define _WQ_DO_GLM          "api_do_glm"
# define _WQ_CLEAN_GLM       "api_clean_glm"
# define _WQ_INIT_GLM_OUTPUT "api_init_glm_output"
# define _WQ_WRITE_GLM_      "api_write_glm"
# define _WQ_VAR_INDEX_C     "api_var_index_c"
# define _WQ_IS_VAR          "api_is_var"
# define _WQ_INFLOW_UPDATE   "api_update_inflow_wq"
#endif
!
!-------------------------------------------------------------------------------
!MODULE DATA

   AED_REAL :: par_fraction =  0.450
   AED_REAL :: nir_fraction =  0.510
   AED_REAL :: uva_fraction =  0.035
   AED_REAL :: uvb_fraction =  0.005

   !# Arrays for state and diagnostic variables
   AED_REAL,DIMENSION(:,:),ALLOCATABLE,TARGET :: cc !# water quality array: nlayers, nvars
   AED_REAL,DIMENSION(:),  ALLOCATABLE,TARGET :: cc_hz !# water quality array - benthic: nvars
   AED_REAL,DIMENSION(:,:),ALLOCATABLE,TARGET :: cc_diag
   AED_REAL,DIMENSION(:),  ALLOCATABLE,TARGET :: cc_diag_hz

   !# Arrays for work, vertical movement, and cross-boundary fluxes
!  AED_REAL,DIMENSION(:,:),ALLOCATABLE :: ws
   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: dz
   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: tss
   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: sed_zones

   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: depth
   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: col_depth

   !# Arrays for environmental variables not supplied externally.
   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: par
   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: pres
   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: uva
   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: uvb
   AED_REAL,DIMENSION(:),ALLOCATABLE,TARGET :: nir

   !# External variables
   AED_REAL,DIMENSION(:),POINTER :: temp
   AED_REAL,DIMENSION(:),POINTER :: salt
   AED_REAL,DIMENSION(:),POINTER :: rad
   AED_REAL,DIMENSION(:),POINTER :: rho
   AED_REAL,DIMENSION(:),POINTER :: area
   AED_REAL,DIMENSION(:),POINTER :: extc
   AED_REAL,DIMENSION(:),POINTER :: layer_stress
   AED_REAL,DIMENSION(:),POINTER :: cvel

   AED_REAL,DIMENSION(:),POINTER :: rain
   AED_REAL,DIMENSION(:),POINTER :: evap
   AED_REAL,DIMENSION(:),POINTER :: air_temp
   AED_REAL,DIMENSION(:),POINTER :: humidity
   AED_REAL,DIMENSION(:),POINTER :: I_0
   AED_REAL,DIMENSION(:),POINTER :: wnd
   AED_REAL,DIMENSION(:),POINTER :: air_pres

!  CHARACTER(len=48),ALLOCATABLE :: names(:)
!  CHARACTER(len=48),ALLOCATABLE :: bennames(:)
!  CHARACTER(len=48),ALLOCATABLE :: diagnames(:)

!  AED_REAL,DIMENSION(:),ALLOCATABLE :: min_
!  AED_REAL,DIMENSION(:),ALLOCATABLE :: max_

   INTEGER,DIMENSION(:),ALLOCATABLE :: externalid
   INTEGER,DIMENSION(:),ALLOCATABLE :: zexternalid
   INTEGER,DIMENSION(:),ALLOCATABLE :: plot_id_v, plot_id_sv, plot_id_d, plot_id_sd

!  LOGICAL :: reinited = .FALSE.

   INTEGER :: n_aed_vars, n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet
!  INTEGER :: zone_var = 0

   CHARACTER(len=64) :: NULCSTR = ""

!===============================================================================
CONTAINS


!###############################################################################
SUBROUTINE api_init_glm(i_fname, len, NumWQ_Vars, NumWQ_Ben)                   &
                                                      BIND(C, name=_WQ_INIT_GLM)
!-------------------------------------------------------------------------------
! Initialize the GLM-AED driver by reading settings from aed.nml.
!-------------------------------------------------------------------------------
!ARGUMENTS
   CCHARACTER,INTENT(in) :: i_fname(*)
   CSIZET,INTENT(in)     :: len
   CINTEGER,INTENT(out)  :: NumWQ_Vars, NumWQ_Ben
!
!LOCALS
   CHARACTER(len=80) :: fname

   TYPE(api_config_t) :: conf
!
!-------------------------------------------------------------------------------
!BEGIN
   CALL make_string(fname, i_fname, len)

   conf%MaxLayers = MaxLayers
   conf%par_fraction =  0.450
   conf%nir_fraction =  0.510
   conf%uva_fraction =  0.035
   conf%uvb_fraction =  0.005

   conf%mobility_off = mobility_off
   conf%bioshade_feedback = bioshade_feedback
   conf%repair_state = repair_state
   conf%link_rain_loss = link_rain_loss
   conf%link_solar_shade = link_solar_shade
   conf%link_bottom_drag = link_bottom_drag

   conf%split_factor = split_factor
   conf%benthic_mode = benthic_mode

   conf%rain_factor = rain_factor
   conf%sw_factor = sw_factor
   conf%friction = friction

   conf%Kw = Kw

   CALL aed_config_model(conf)

   CALL api_set_glm_env()

   n_aed_vars = aed_init_model(fname, n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet)

   NumWQ_Vars  = n_vars
   NumWQ_Ben   = n_vars_ben

   ALLOCATE(plot_id_v(n_vars))
   ALLOCATE(plot_id_sv(n_vars_ben))
   ALLOCATE(plot_id_d(n_vars_diag))
   ALLOCATE(plot_id_sd(n_vars_diag_sheet))
   plot_id_v = -1; plot_id_sv = -1; plot_id_d = -1; plot_id_sd = -1

   ALLOCATE(externalid(n_aed_vars))
   ALLOCATE(zexternalid(n_aed_vars))
END SUBROUTINE api_init_glm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE api_set_glm_env()
!-------------------------------------------------------------------------------
!ARGUMENTS
!
!LOCALS
   INTEGER :: status

   TYPE(api_env_t) :: env
!
!-------------------------------------------------------------------------------
!BEGIN
   !# Save pointers to external dynamic variables that we need later (in do_glm_wq)
   lheights => theLake%Height
   temp     => theLake%Temp
   salt     => theLake%Salinity
   rho      => theLake%Density
   area     => theLake%LayerArea
   rad      => theLake%Light
   cvel     => theLake%Umean
   extc     => theLake%ExtcCoefSW
   layer_stress => theLake%LayerStress

   !# Provide pointers to arrays with environmental variables to aed.
   wnd      => aMetData%WindSpeed
   rain     => aMetData%Rain
   I_0      => aMetData%ShortWave
   air_temp => aMetData%AirTemp
   air_pres => aMetData%AirPres
   humidity => aMetData%RelHum

   evap => aSurfData%Evap

   !# Allocate array for photosynthetically active radiation (PAR).
   !# This will be calculated internally during each time step.
   ALLOCATE(par(MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (par)'
   par = zero_

   ALLOCATE(nir(MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (nir)'
   nir = zero_
   ALLOCATE(uva(MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (uva)'
   uva = zero_
   ALLOCATE(uvb(MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (uvb)'
   uvb = zero_

   !# Allocate array for local pressure.
   !# This will be calculated [approximated] from layer depths internally
   !# during each time step.
   ALLOCATE(pres(MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (pres)'
   pres = zero_

   ALLOCATE(dz(MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (dz)'
   dz = zero_
   ALLOCATE(depth(MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (depth)'
   depth = zero_
   ALLOCATE(sed_zones(MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (sed_zones)'
   sed_zones = zero_

   ALLOCATE(col_depth(1),stat=status)

   ALLOCATE(tss(MaxLayers),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (tss)'
   tss = zero_

   env%yearday   => yearday
   timestep = dt
   env%timestep  => timestep

   env%longitude => longitude
   env%latitude  => latitude

   env%temp          => temp
   env%salt          => salt
   env%rho           => rho
   env%dz            => dz
   env%height        => lheights
   env%area          => area
   env%depth         => depth
   env%col_depth     => col_depth
   env%extc          => extc
   env%tss           => tss
!  env%ss1           => ss1
!  env%ss2           => ss2
!  env%ss3           => ss3
!  env%ss4           => ss4
   env%cvel          => cvel
!  env%vvel          => vvel
!  env%bio_drag      => bio_drag
   env%rad           => rad
   env%I_0           => I_0
   env%wnd           => wnd
   env%air_temp      => air_temp
   env%air_pres      => air_pres
   env%rain          => rain
   env%evap          => evap
   env%humidity      => humidity
!  env%longwave      => longwave
!  env%bathy         => bathy
!  env%rainloss      => rainloss
!  env%ustar_bed     => ustar_bed
!  env%wv_uorb       => wv_uorb
!  env%wv_t          => wv_t
   env%layer_stress  => layer_stress
   env%sed_zones     => sed_zones
   env%par => par
   env%nir => nir
   env%uva => uva
   env%uvb => uvb

   env%pres => pres

   CALL aed_set_model_env(env)
END SUBROUTINE api_set_glm_env
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE api_set_glm_data()                     BIND(C, name=_WQ_SET_GLM_DATA)
!-------------------------------------------------------------------------------
!ARGUMENTS
!
!LOCALS
   INTEGER :: status

   TYPE(api_data_t) :: dat
!
!-------------------------------------------------------------------------------
!BEGIN
   IF (n_zones .GT. 0) &
      CALL api_set_glm_zones(n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet)

   !# Now that we know how many vars we need, we can allocate space for them
   ALLOCATE(cc(MaxLayers, (n_vars + n_vars_ben)),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (cc)'
   cc = zero_         !# initialise to zero

   CALL set_c_wqvars_ptr(cc)

!  ALLOCATE(cc_hz(n_vars_ben),stat=status)
!  IF (status /= 0) STOP 'allocate_memory(): Error allocating (cc_hz)'
!  cc_hz = zero_         !# initialise to zero

   !# Allocate diagnostic variable array and set all values to zero.
   !# (needed because time-integrated/averaged variables will increment rather than set the array)
   ALLOCATE(cc_diag(MaxLayers, n_vars_diag),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (cc_diag)'
   cc_diag = zero_

   !# Allocate diagnostic variable array and set all values to zero.
   !# (needed because time-integrated/averaged variables will increment rather than set the array)
   ALLOCATE(cc_diag_hz(n_vars_diag_sheet),stat=status)
   IF (status /= 0) STOP 'allocate_memory(): Error allocating (cc_diag_hz)'
   cc_diag_hz = zero_

   CALL set_c_wqdvars_ptr(cc_diag, cc_diag_hz, n_vars_diag, n_vars_diag_sheet)

   dat%cc => cc
!  dat%cc_hz => cc_hz
   dat%cc_diag => cc_diag
   dat%cc_diag_hz => cc_diag_hz

   CALL aed_set_model_data(dat)
END SUBROUTINE api_set_glm_data
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
CINTEGER FUNCTION api_is_var(id,i_vname,len)            BIND(C, name=_WQ_IS_VAR)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER,INTENT(in)   :: id
   CCHARACTER,INTENT(in) :: i_vname(*)
   CSIZET,INTENT(in)     :: len
!
!LOCALS
   CHARACTER(len=45) :: vname
   TYPE(aed_variable_t),POINTER :: tvar
   INTEGER :: i, v, sv, d, sd
!
!-------------------------------------------------------------------------------
!BEGIN
   CALL make_string(vname, i_vname, len)

   v = 0; sv = 0; d = 0; sd = 0
   DO i=1,n_aed_vars
      IF ( aed_get_var(i, tvar) ) THEN
         IF ( .NOT. (tvar%diag .OR. tvar%extern) ) THEN
            IF ( tvar%sheet ) THEN ; sv=sv+1; ELSE ; v=v+1 ; ENDIF
            IF ( TRIM(tvar%name) == vname ) THEN
               IF (tvar%sheet) THEN
                  api_is_var=-sv
                  plot_id_sv(sv) = id;
               ELSE
                  api_is_var=v
                  plot_id_v(v) = id;
               ENDIF
               RETURN
            ENDIF
         ELSEIF ( tvar%diag ) THEN
            IF ( tvar%sheet ) THEN ; sd=sd+1; ELSE ; d=d+1 ; ENDIF
            IF ( TRIM(tvar%name) == vname ) THEN
               IF (tvar%sheet) THEN
                  api_is_var=-sd
                  plot_id_sd(sd) = id;
               ELSE
                  api_is_var=d
                  plot_id_d(d) = id;
               ENDIF
               RETURN
            ENDIF
         ENDIF
      ENDIF
   ENDDO
   api_is_var = 0
END FUNCTION api_is_var
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE doMobilityF(N,dt,h,A,ww,min_C,cc)
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in)     :: N       !# number of vertical layers
   AED_REAL,INTENT(in)    :: dt      !# time step (s)
   AED_REAL,INTENT(in)    :: h(*)    !# layer thickness (m)
   AED_REAL,INTENT(in)    :: A(*)    !# layer areas (m2)
   AED_REAL,INTENT(in)    :: ww(*)   !# vertical speed (m/s)
   AED_REAL,INTENT(in)    :: min_C   !# minimum allowed cell concentration
   AED_REAL,INTENT(inout) :: cc(*)   !# cell concentration
!
!LOCALS
   CINTEGER :: NC
!
!-------------------------------------------------------------------------------
!BEGIN
   NC = N
   CALL doMobility(NC,dt,h,A,ww,min_C,cc)
END SUBROUTINE doMobilityF
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE api_do_glm(wlev, pIce)                       BIND(C, name=_WQ_DO_GLM)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER,INTENT(in) :: wlev
   CLOGICAL,INTENT(in) :: pIce
!
!LOCALS
   LOGICAL :: doSurface
   INTEGER :: lev
   AED_REAL :: surf

   PROCEDURE(aed_mobility_t),POINTER :: doMobilityP
!
!-------------------------------------------------------------------------------
!BEGIN
   col_depth = lheights(wlev)
   surf = lheights(wlev)
   !# re-compute the layer heights and depths
   dz(1) = lheights(1)
   depth(1) = surf - lheights(1)
   DO lev=2,wlev
      dz(lev) = lheights(lev) - lheights(lev-1)
      depth(lev) = surf - lheights(lev)
   ENDDO

   !# Calculate local pressure
   pres(1:wlev) = -lheights(1:wlev)

   doSurface = .NOT. pIce

   doMobilityP => doMobilityF
   CALL aed_set_mobility(doMobilityP)
   CALL aed_run_model(1, wlev, doSurface)
END SUBROUTINE api_do_glm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE api_clean_glm()                           BIND(C, name=_WQ_CLEAN_GLM)
!-------------------------------------------------------------------------------
! Finish biogeochemical model
!-------------------------------------------------------------------------------
!BEGIN
   CALL aed_clean_model()
   ! Deallocate internal arrays
!  IF (ALLOCATED(cc))         DEALLOCATE(cc)
!  IF (ALLOCATED(cc_hz))      DEALLOCATE(cc_hz)
!  IF (ALLOCATED(cc_diag))    DEALLOCATE(cc_diag)
!  IF (ALLOCATED(cc_diag_hz)) DEALLOCATE(cc_diag_hz)
   IF (ALLOCATED(par))        DEALLOCATE(par)
   IF (ALLOCATED(nir))        DEALLOCATE(nir)
   IF (ALLOCATED(uva))        DEALLOCATE(uva)
   IF (ALLOCATED(uvb))        DEALLOCATE(uvb)
   IF (ALLOCATED(pres))       DEALLOCATE(pres)
   IF (ALLOCATED(dz))         DEALLOCATE(dz)
   IF (ALLOCATED(tss))        DEALLOCATE(tss)
   IF (ALLOCATED(sed_zones))  DEALLOCATE(sed_zones)
   IF (ALLOCATED(depth))      DEALLOCATE(depth)
   IF (ALLOCATED(col_depth))  DEALLOCATE(col_depth)
END SUBROUTINE api_clean_glm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE api_init_glm_output(ncid,x_dim,y_dim,z_dim,zone_dim,time_dim)       &
                                               BIND(C, name=_WQ_INIT_GLM_OUTPUT)
!-------------------------------------------------------------------------------
!  Initialize the output by defining biogeochemical variables.
!-------------------------------------------------------------------------------
!
!ARGUMENTS
   CINTEGER,INTENT(in) :: ncid,x_dim,y_dim,z_dim,zone_dim,time_dim
!
!LOCALS
   INTEGER i !, v, d
   INTEGER dims(4)

   TYPE(aed_variable_t),POINTER :: tv
!
!-------------------------------------------------------------------------------
!BEGIN
   !# Put NetCDF library in define mode.
   CALL define_mode_on(ncid)

   !# Set up dimension indices for 3D (+ time) variables (longitude,latitude,depth,time).
   dims(1) = x_dim
   dims(2) = y_dim
   dims(3) = z_dim
   dims(4) = time_dim

!  v = 0; d = 0
   DO i=1,n_aed_vars
      IF ( aed_get_var(i, tv) ) THEN
         IF ( .NOT. (tv%sheet .OR. tv%extern) ) THEN
            !# only for state and diag vars that are not sheet
            externalid(i) = NEW_NC_VARIABLE(ncid, TRIM(tv%name), LEN_TRIM(tv%name), NF90_REALTYPE, 4, dims(1:4))
            CALL set_nc_attributes(ncid, externalid(i), MYTRIM(tv%units), MYTRIM(tv%longname) PARAM_FILLVALUE)
         ENDIF
      ENDIF
   ENDDO

   !# Set up dimension indices for 2D (+ time) variables (longitude,latitude,time).
   dims(1) = x_dim
   dims(2) = y_dim
   dims(3) = time_dim

!  v = 0; d = 0
   DO i=1,n_aed_vars
      IF ( aed_get_var(i, tv) ) THEN
         IF ( tv%sheet .AND. .NOT. tv%extern ) THEN
            !# only for state and diag sheet vars
            externalid(i) = NEW_NC_VARIABLE(ncid, TRIM(tv%name), LEN_TRIM(tv%name), NF90_REALTYPE, 3, dims(1:3))
            CALL set_nc_attributes(ncid, externalid(i), MYTRIM(tv%units), MYTRIM(tv%longname) PARAM_FILLVALUE)
         ENDIF
      ENDIF
   ENDDO

   !# Set up dimension indices for 2D+zone (+ time) variables (longitude,latitude,zone,time).
   IF ( n_zones .GT. 0 ) THEN
      dims(1) = x_dim
      dims(2) = y_dim
      dims(3) = zone_dim
      dims(4) = time_dim

!print*,x_dim,y_dim,zone_dim,time_dim
!     v = 0; d = 0
      DO i=1,n_aed_vars
         IF ( aed_get_var(i, tv) ) THEN
            IF ( tv%sheet .AND. .NOT. tv%extern ) THEN
               !# only for state and diag sheet vars
               zexternalid(i) = NEW_NC_VARIABLE(ncid, TRIM(tv%name)//"_Z", LEN_TRIM(tv%name)+2, NF90_REALTYPE, 4, dims(1:4))
               CALL set_nc_attributes(ncid, zexternalid(i), MYTRIM(tv%units), MYTRIM(tv%longname) PARAM_FILLVALUE)
            ENDIF
         ENDIF
      ENDDO
   ENDIF

   !# Take NetCDF library out of define mode (ready for storing data).
   CALL define_mode_off(ncid)
END SUBROUTINE api_init_glm_output
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
!  Save properties of biogeochemical model, including state variable
!  values, diagnostic variable values, and sums of conserved quantities.
!-------------------------------------------------------------------------------
SUBROUTINE api_write_glm(ncid,wlev,nlev,lvl,point_nlevs) BIND(C, name=_WQ_WRITE_GLM_)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER,INTENT(in) :: ncid, wlev, nlev
   CINTEGER,INTENT(in) :: lvl(*), point_nlevs
!
!LOCALS
   TYPE(aed_variable_t),POINTER :: tv

   INTEGER  :: i, j, v, d, sv, sd
   INTEGER  :: z
   AED_REAL :: val_out
   CLOGICAL :: last = .FALSE.
!
!-------------------------------------------------------------------------------
!BEGIN
   v = 0; d = 0; sv = 0; sd = 0
   DO i=1,n_aed_vars
      IF ( aed_get_var(i, tv) ) THEN
         IF ( tv%diag ) THEN
            !# Process and store diagnostic variables.
            IF ( tv%sheet ) THEN
               sd = sd + 1
               !# Process and store diagnostic variables defined on horizontal slices of the domain.
               IF ( n_zones .GT. 0 ) THEN
                  z_diag_hz(n_zones+1,sd) = cc_diag_hz(sd)
                  CALL store_nc_array(ncid, zexternalid(i), XYNT_SHAPE, n_zones, n_zones, array=z_diag_hz(1:n_zones+1,sd))
               ENDIF
               CALL store_nc_scalar(ncid, externalid(i), XYT_SHAPE, scalar=cc_diag_hz(sd))
               IF ( do_plots .AND. plot_id_sd(sd).GE.0 ) THEN
                  IF ( n_zones .GT. 0 ) THEN
                     DO z=1,n_zones
                        CALL put_glm_val_z(plot_id_sd(sd),z_diag_hz(z, sd), z)
                     ENDDO
                  ENDIF
                  CALL put_glm_val_s(plot_id_sd(sd),cc_diag_hz(sd))
               ENDIF
               DO j=1,point_nlevs
                  val_out = missing
                  IF ((lvl(j) .EQ. wlev) .AND. tv%top) val_out = cc_diag(1, v)
                  IF ((lvl(j) .EQ. 0)    .AND. tv%bot) val_out = cc_diag(1, v)
                  CALL write_csv_point(j, tv%name, len_trim(tv%name), val_out, NULCSTR, 0, last=last)
               ENDDO
            ELSE  !# not sheet
               d = d + 1
               !# Store diagnostic variable values defined on the full domain.
               CALL store_nc_array(ncid, externalid(i), XYZT_SHAPE, wlev, nlev, array=cc_diag(:, d))
               IF ( do_plots .AND. plot_id_d(d).GE.0 ) &
                  CALL put_glm_val(plot_id_d(d), cc_diag(1:wlev, d))
               DO j=1,point_nlevs
                  IF (lvl(j) .GE. 0) THEN ; val_out = cc_diag(lvl(j)+1, d)
                  ELSE                    ; val_out = missing     ; ENDIF
                  CALL write_csv_point(j, tv%name, len_trim(tv%name), val_out, NULCSTR, 0, last=last)
                  CALL write_csv_point_avg(j, tv%name, len_trim(tv%name), cc_diag(:, d), NULCSTR, 0, last=last)
               ENDDO
            ENDIF
         ELSE IF ( .NOT. tv%extern ) THEN  ! not diag
            IF ( tv%sheet ) THEN
               sv = sv + 1
               !# Store benthic biogeochemical state variables.
               IF ( n_zones .GT. 0 ) THEN
                  CALL store_nc_array(ncid, zexternalid(i), XYNT_SHAPE, n_zones, n_zones, array=z_cc(1:n_zones, 1, n_vars+sv))
               ENDIF
               CALL store_nc_scalar(ncid, externalid(i), XYT_SHAPE, scalar=cc(1, n_vars+sv))
               IF ( do_plots .AND. plot_id_sv(sv).GE.0 ) THEN
                  IF ( n_zones .GT. 0 ) THEN
                     DO z=1,n_zones
                        CALL put_glm_val_z(plot_id_sv(sv), z_cc(z, 1, n_vars+sv), z)
                     ENDDO
                  ENDIF
                  CALL put_glm_val_s(plot_id_sv(sv), cc(1, n_vars+sv))
               ENDIF
               DO j=1,point_nlevs
                  val_out = missing
                  IF ((lvl(j) .EQ. wlev) .AND. tv%top) val_out = cc(1, v)
                  IF ((lvl(j) .EQ. 0)    .AND. tv%bot) val_out = cc(1, v)
                  CALL write_csv_point(j, tv%name, len_trim(tv%name), val_out, NULCSTR, 0, last=last)
               ENDDO
            ELSE     !# not sheet
               v = v + 1
               !# Store pelagic biogeochemical state variables.
               CALL store_nc_array(ncid, externalid(i), XYZT_SHAPE, wlev, nlev, array=cc(:, v))
               IF ( do_plots .AND. plot_id_v(v).GE.0 ) CALL put_glm_val(plot_id_v(v), cc(1:wlev, v))
               DO j=1,point_nlevs
                  IF (lvl(j) .GE. 0) THEN ; val_out = cc(lvl(j)+1, v)
                  ELSE                    ; val_out = missing     ; ENDIF
                  CALL write_csv_point(j, tv%name, len_trim(tv%name), val_out, NULCSTR, 0, last=last)
                  CALL write_csv_point_avg(j, tv%name, len_trim(tv%name), cc(:, v), NULCSTR, 0, last=last)
               ENDDO
            ENDIF
         ENDIF
      ENDIF
   ENDDO
END SUBROUTINE api_write_glm
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE wq_inflow_update(wqinf, nwqVars, temp, salt)                        &
                                                 BIND(C, name=_WQ_INFLOW_UPDATE)
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(C_PTR),VALUE :: wqinf
   CINTEGER,INTENT(in) :: nwqVars
   AED_REAL,INTENT(inout) :: temp, salt
!
!LOCALS
   AED_REAL,DIMENSION(:),POINTER :: wqInfF
!
!-------------------------------------------------------------------------------
!BEGIN
   CALL C_F_POINTER(wqinf, wqInfF, [nwqVars])
   CALL aed_inflow_update(wqInfF, temp, salt)
END SUBROUTINE wq_inflow_update
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
CINTEGER FUNCTION aed_var_index_c(name, len)       BIND(C, name=_WQ_VAR_INDEX_C)
!-------------------------------------------------------------------------------
!ARGUMENTS
   CCHARACTER,INTENT(in) :: name(*)
   CSIZET,INTENT(in)     :: len
!
!LOCALS
   CHARACTER(len=len) :: tn
!
!-------------------------------------------------------------------------------
!BEGIN
   tn = trim(transfer(name(1:len),tn))
   aed_var_index_c = aed_var_index(tn) - 1
END FUNCTION aed_var_index_c
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE glm_api_aed
