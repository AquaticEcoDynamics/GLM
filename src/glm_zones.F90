!###############################################################################
!#                                                                             #
!# glm_zones.F90                                                               #
!#                                                                             #
!# The sediment zone processing bit for WQ                                     #
!#                                                                             #
!# Developed by :                                                              #
!#     AquaticEcoDynamics (AED) Group                                          #
!#     School of Agriculture and Environment                                   #
!#     The University of Western Australia                                     #
!#                                                                             #
!#     http://aquatic.science.uwa.edu.au/                                      #
!#                                                                             #
!# Copyright 2013-2026 : The University of Western Australia                   #
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

#undef MISVAL

#include "glm.h"

MODULE glm_zones

   USE ISO_C_BINDING

   USE glm_types
   USE aed_common

   IMPLICIT NONE

   PRIVATE ! By default, make everything private

   AED_REAL,DIMENSION(:,:,:),ALLOCATABLE,TARGET :: z_cc      ! (n_vars, n_levs, nsed_zones)
!  AED_REAL,DIMENSION(:,:),  ALLOCATABLE,TARGET :: z_cc_hz   ! (n_vars,         nsed_zones)
   AED_REAL,DIMENSION(:,:,:),ALLOCATABLE,TARGET :: z_diag    ! (n_vars, n_levs, nsed_zones)
   AED_REAL,DIMENSION(:,:),  ALLOCATABLE,TARGET :: z_diag_hz ! (n_vars,         nsed_zones)

   AED_REAL,DIMENSION(:),POINTER :: lheights

   !# Layer-zone contact shares for the current host step, computed once by
   !# set_zone_shares (from calc_zone_areas) and used by every mapping, including
   !# the flux delivery in glm_aed.F90:
   !#   zshare(lev,zon)  fraction of layer lev's bed lying in zone zon (zone_bed_frac)
   !#   zhshare(lev,zon) fraction of layer lev's height lying in zone zon (zone_frac)
   !#   zwbed(zon), zwvol(zon)  the zone's wet bed area and water volume from them
   AED_REAL,DIMENSION(:,:),ALLOCATABLE :: zshare, zhshare
   AED_REAL,DIMENSION(:),  ALLOCATABLE :: zwbed, zwvol
   LOGICAL :: shares_set = .FALSE.
   !# indices (into x_diag) of the zone-averaged (zavg) non-sheet diagnostics
   INTEGER,DIMENSION(:),ALLOCATABLE :: zavg_idx
   INTEGER,DIMENSION(:),ALLOCATABLE :: bflux_idx   !# bflux layer diagnostics (K88)

   AED_REAL,DIMENSION(:),POINTER :: zone_heights
   INTEGER :: nvars, nbenv, nvdiag, nvdiag_hz

   PUBLIC n_zones, zone_heights, lheights
   PUBLIC wq_set_glm_zones, copy_from_zone, copy_to_zone, calc_zone_areas, zone_frac, zone_bed_frac
   PUBLIC zshare

!  PUBLIC z_cc, z_cc_hz, z_diag, z_diag_hz, theZones, theLake
   PUBLIC z_cc, z_diag, z_diag_hz, theZones, theLake

CONTAINS

!###############################################################################
SUBROUTINE wq_set_glm_zones(numVars, numBenV, numDiagV, numDiagHzV)            &
                                                BIND(C, name="wq_set_glm_zones")
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER,INTENT(in) :: numVars, numBenV, numDiagV, numDiagHzV
!
!LOCALS
!
!-------------------------------------------------------------------------------
!BEGIN
   nvars = numVars
   nbenv = numBenV
   nvdiag = numDiagV
   nvdiag_hz = numDiagHzV

   lheights => theLake%Height
   zone_heights => theZones%zheight

   ALLOCATE(z_cc(numVars+numBenV, MaxLayers, n_zones)) ; z_cc = 0.
!  ALLOCATE(z_cc_hz(numVars+numBenV, n_zones))         ; z_cc_hz = 0.
   ALLOCATE(z_diag(numDiagV, MaxLayers, n_zones))      ; z_diag = 0.
   ALLOCATE(z_diag_hz(numDiagHzV, n_zones+1))          ; z_diag_hz = 0.
   theZones%zarea = 0.
END SUBROUTINE wq_set_glm_zones
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE set_zone_shares(wlev)
!-------------------------------------------------------------------------------
! The layer-zone contact shares for the current layer structure (see zshare).
! Layers move only between host steps, so this runs once per step, from
! calc_zone_areas; the routines that map between layers and zones read the
! stored shares instead of recomputing them.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: wlev
!LOCALS
   INTEGER  :: lev, zon
   AED_REAL :: bed
!-------------------------------------------------------------------------------
!BEGIN
   IF ( .NOT. ALLOCATED(zshare) ) THEN
      ALLOCATE(zshare(MaxLayers, n_zones), zhshare(MaxLayers, n_zones))
      ALLOCATE(zwbed(n_zones), zwvol(n_zones))
   ENDIF
   zshare = 0. ; zhshare = 0. ; zwbed = 0. ; zwvol = 0.
   DO lev=1,wlev
      bed = theLake(lev)%LayerArea ; IF ( lev > 1 ) bed = bed - theLake(lev-1)%LayerArea
      DO zon=1,n_zones
         zshare(lev,zon)  = zone_bed_frac(lev, zon)
         zhshare(lev,zon) = zone_frac(lev, zon)
         zwbed(zon) = zwbed(zon) + zshare(lev,zon) * bed
         zwvol(zon) = zwvol(zon) + zhshare(lev,zon) * theLake(lev)%LayerVol
      ENDDO
   ENDDO
   shares_set = .TRUE.
END SUBROUTINE set_zone_shares
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE calc_zone_areas(areas, wlev, surf)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,DIMENSION(:),INTENT(in) :: areas
   INTEGER,INTENT(in) :: wlev
   AED_REAL,INTENT(in) :: surf
!
!LOCALS
   INTEGER  :: lev, zon
   LOGICAL  :: w_zones
   AED_REAL :: bed
!
!-------------------------------------------------------------------------------
!BEGIN
   theZones%zarea = 0.  ; theZones%z_pc_wet = 0.
   w_zones = .FALSE.

   !# Zone areas: each layer's bed (areas() is CUMULATIVE here) is shared
   !# between zones by the part of that bed lying in each zone
   !# (zone_bed_frac), the same shares used to map zone fluxes onto layers. This
   !# handles a layer whose bottom sits exactly on a boundary, a bottom layer
   !# that straddles a boundary and a layer spanning several zones.
   CALL set_zone_shares(wlev)
   DO lev=1, wlev
      bed = areas(lev) ; IF ( lev > 1 ) bed = areas(lev) - areas(lev-1)
      DO zon=1, n_zones
         theZones(zon)%zarea = theZones(zon)%zarea + bed * zshare(lev,zon)
      ENDDO
   ENDDO

   !# Wet fraction of each zone (unchanged): zones below the surface are fully
   !# wet, the first zone reaching above it is partly wet, those above are dry.
   DO zon=1, n_zones
      IF ( zone_heights(zon) <= surf ) THEN
         theZones(zon)%z_pc_wet = 1.0
      ELSEIF ( .NOT. w_zones ) THEN
         w_zones = .TRUE.
         IF ( zon > 1 ) THEN
            theZones(zon)%z_pc_wet = surf / (zone_heights(zon) - zone_heights(zon-1))
         ELSE
            theZones(zon)%z_pc_wet = surf / zone_heights(zon)
         ENDIF
      ENDIF
   ENDDO

   theZones(1:n_zones)%zpres = -zone_heights(1:n_zones)
END SUBROUTINE calc_zone_areas
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE copy_to_zone(x_cc, x_diag, wlev)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,DIMENSION(:,:),INTENT(in) :: x_cc
   AED_REAL,DIMENSION(:,:),INTENT(in) :: x_diag
   INTEGER,INTENT(in) :: wlev
!
!LOCALS
   INTEGER  :: zon, lev, a_zones
   AED_REAL :: surf, f, bed, w
   AED_REAL :: wbed(n_zones), wvol(n_zones)
   LOGICAL  :: w_zones
!
!-------------------------------------------------------------------------------
!BEGIN
   ! Reset zone data structure to zero, for copying in information from main
   ! lake data structures. Note that z_cc(:,:,nvars+1:nvars_ben) is not zeroed
   ! as this is the benthic zone data we need to preserve
   z_cc(1:nvars,:,:) = 0.
   z_diag(:,:,:) = 0.
   z_diag_hz(:,:) = 0.

   ! Initialise the zone environment
   theZones%zrad = 0.
   theZones%zsalt = 0.
   theZones%ztemp = 0.
   theZones%zrho = 0.
   theZones%zextc = 0.
   theZones%zlayer_stress = 0.
   theZones%ztss = 0.
   theZones%zpar = 0.
   theZones%znir = 0.
   theZones%zuva = 0.
   theZones%zuvb = 0.
   theZones(1)%z_sed_zones = 1.
   theZones%zvel = 0.

   ! Populate the 1st layer in each zone structure with the zone-averaged
   ! quantity (what each zone's sediment sees). Each layer is weighted by the
   ! bed area it touches within the zone (zone_bed_frac * bed area), the same
   ! shares used for the zone areas and to map zone fluxes onto layers. A wet
   ! zone with no bed of its own (vertical wall) falls back to water-volume
   ! weights, shared by height (zone_frac). Zones entirely above the water are
   ! left at zero, as before.
   ! Benthic vars (nvars+1 ->) are not touched: their zone values are kept.
   w_zones = .FALSE.
   IF ( .NOT. shares_set ) CALL set_zone_shares(wlev)

   a_zones = 1
   DO lev=1,wlev
      bed = theLake(lev)%LayerArea ; IF ( lev > 1 ) bed = bed - theLake(lev-1)%LayerArea
      DO zon=1,n_zones
         IF ( zwbed(zon) > 0. ) THEN
            f = zshare(lev,zon)
            IF ( f <= 0. ) CYCLE
            w = f * bed / zwbed(zon)
         ELSE
            f = zhshare(lev,zon)
            IF ( f <= 0. ) CYCLE
            w = f * theLake(lev)%LayerVol / zwvol(zon)
         ENDIF
         IF ( w <= 0. ) CYCLE
         a_zones = MAX(a_zones, zon)

         z_cc(1:nvars,1,zon) = z_cc(1:nvars,1,zon) + w * x_cc(1:nvars,lev)  ! 1:nvars is the water column "cc", not benthic
         z_diag(:,zon,zon)   = z_diag(:,zon,zon)   + w * x_diag(:,lev)      ! store at level=zon so cell(zon) is correct with layer_idx=zon

         theZones(zon)%ztemp         = theZones(zon)%ztemp + w * theLake(lev)%Temp
         theZones(zon)%zsalt         = theZones(zon)%zsalt + w * theLake(lev)%Salinity
         theZones(zon)%zrho          = theZones(zon)%zrho  + w * theLake(lev)%Density
         theZones(zon)%zrad          = theZones(zon)%zrad  + w * theLake(lev)%Light
         theZones(zon)%zpar          = theZones(zon)%zpar  + w * theLake(lev)%Light*0.45
         theZones(zon)%zvel          = theZones(zon)%zvel  + w * theLake(lev)%Umean
         theZones(zon)%zextc         = theZones(zon)%zextc + w * theLake(lev)%ExtcCoefSW
         theZones(zon)%zlayer_stress = theZones(zon)%zlayer_stress + w * theLake(lev)%LayerStress
      ENDDO
   ENDDO

   theZones%zdz = 0.
   surf = lheights(wlev)
   IF ( surf > zone_heights(1) ) THEN
      theZones(1)%zdepth = zone_heights(1)
      theZones(1)%zdz = zone_heights(1)
   ELSE
      theZones(1)%zdepth = surf
      theZones(1)%zdz = surf
      w_zones = .TRUE.
   ENDIF

   DO zon=2,a_zones
      theZones(zon)%z_sed_zones = zon
      IF ( .NOT. w_zones ) THEN
          IF ( surf > zone_heights(zon) ) THEN
             theZones(zon)%zdepth = zone_heights(zon)
             theZones(zon)%zdz =  zone_heights(zon) - zone_heights(zon-1)
          ELSE
             theZones(zon)%zdepth = surf
             theZones(zon)%zdz = surf - zone_heights(zon-1)
             w_zones = .TRUE.
          ENDIF
      ELSE
         theZones(zon)%zdepth = surf
         theZones(zon)%zdz = surf - zone_heights(zon-1)
      ENDIF
   ENDDO

END SUBROUTINE copy_to_zone
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE copy_from_zone(n_aed_vars, x_cc, x_diag, x_diag_hz, wlev)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,DIMENSION(:,:),INTENT(inout) :: x_cc
   AED_REAL,DIMENSION(:,:),INTENT(inout) :: x_diag
   AED_REAL,DIMENSION(:),INTENT(inout) :: x_diag_hz
   INTEGER,INTENT(in) :: wlev, n_aed_vars
!
!LOCALS
   INTEGER  :: zon, lev, v_start, v_end
   AED_REAL :: scale, area, ring
!  LOGICAL  :: column_benthic_var_averaging = .false.
!  INTEGER  :: water_column_zone = 1
   INTEGER  :: i, j, k
   TYPE(aed_variable_t),POINTER :: tvar
!
!-------------------------------------------------------------------------------
!BEGIN
   v_start = nvars+1 ; v_end = nvars+nbenv

!print*,"size x_cc 1", size(x_cc, 1)
!print*,"size x_cc 2", size(x_cc, 2)
!print*,"size x_cc_hz 1", size(x_cc_hz, 1)
!
!print*,"size z_cc 1", size(z_cc, 1)
!print*,"size z_cc 2", size(z_cc, 2)
!print*,"size z_cc 3", size(z_cc, 3)
!print*,"size z_cc_hz 1", size(z_cc_hz, 1)
!print*,"size z_cc_hz 2", size(z_cc_hz, 2)
!
!print*,"Z1 cc(:,1) = ", x_cc(:,1), " xch ",x_cc_hz(:), " zcc ",z_cc(:,1,1), " zch ",z_cc_hz(:,1), nvars

   ! Set water layer zone-averaged diagnostics (zavg == true) and benthic vars
   ! from the zones, each zone weighted by the share of the layer's height
   ! bed in it (zone_bed_frac, the same shares used for zone areas and fluxes).
   IF ( .NOT. ALLOCATED(zavg_idx) ) THEN
      !# which non-sheet diagnostics are zone-averaged (fixed once the model is set up)
      ALLOCATE(zavg_idx(n_aed_vars)) ; k = 0 ; j = 0
      DO i=1,n_aed_vars
         IF ( aed_get_var(i, tvar) ) THEN
            IF ( tvar%var_type == V_DIAGNOSTIC .AND. .NOT. tvar%sheet ) THEN
               j = j + 1
               IF ( tvar%zavg ) THEN ; k = k + 1 ; zavg_idx(k) = j ; ENDIF
            ENDIF
         ENDIF
      ENDDO
      zavg_idx = zavg_idx(1:k)
   ENDIF
   IF ( .NOT. shares_set ) CALL set_zone_shares(wlev)

   DO lev=1,wlev
      DO k=1,SIZE(zavg_idx)
         j = zavg_idx(k)
         x_diag(j, lev) = 0.
         DO zon=1,n_zones
            scale = zshare(lev,zon)
            IF ( scale > 0. ) x_diag(j, lev) = x_diag(j, lev) + z_diag(j, zon, zon) * scale   ! level=zon matches copy_to_zone
         ENDDO
      ENDDO
      x_cc(v_start:v_end, lev) = 0.
      DO zon=1,n_zones
         scale = zshare(lev,zon)
         IF ( scale > 0. ) x_cc(v_start:v_end, lev) = x_cc(v_start:v_end, lev) + z_cc(v_start:v_end, 1, zon) * scale
      ENDDO
   ENDDO


   !# Benthic-flux diagnostics (bflux, e.g. OXY_oxy_dsfv): the module sets them
   !# per zone as flux / zone thickness; rebuild each layer's value as the rate the
   !# coupler applies, sum over zones of flux x share x bed ring / layer volume (K88).
   IF ( .NOT. ALLOCATED(bflux_idx) ) THEN
      ALLOCATE(bflux_idx(n_aed_vars)) ; k = 0 ; j = 0
      DO i=1,n_aed_vars
         IF ( aed_get_var(i, tvar) ) THEN
            IF ( tvar%var_type == V_DIAGNOSTIC .AND. .NOT. tvar%sheet ) THEN
               j = j + 1
               IF ( tvar%bflux ) THEN ; k = k + 1 ; bflux_idx(k) = j ; ENDIF
            ENDIF
         ENDIF
      ENDDO
      bflux_idx = bflux_idx(1:k)
   ENDIF
   IF ( SIZE(bflux_idx) > 0 ) THEN
      DO lev=1,wlev
         ring = theLake(lev)%LayerArea ; IF ( lev > 1 ) ring = ring - theLake(lev-1)%LayerArea
         DO k=1,SIZE(bflux_idx)
            j = bflux_idx(k)
            x_diag(j, lev) = 0.
            IF ( theLake(lev)%LayerVol <= 0. ) CYCLE
            DO zon=1,n_zones
               scale = zshare(lev,zon)
               IF ( scale > 0. ) x_diag(j, lev) = x_diag(j, lev) + z_diag(j, zon, zon) * theZones(zon)%zdz * scale
            ENDDO
            x_diag(j, lev) = x_diag(j, lev) * ring / theLake(lev)%LayerVol
         ENDDO
      ENDDO
   ENDIF

   ! Set single-value (lake-wide) sheet diagnostics to the mean of the zones,
   ! weighted by area. Re-zero only the rezero=.TRUE. sheet diagnostics here:
   ! aed_run_model/aed_run_column also calls aed_calculate_riparian once,
   ! unconditionally, for the whole (non-zoned) column before this per-zone
   ! aggregation runs; that call writes a correct single-pass value into the
   ! same sheet-diagnostic slots via a plain assignment. Without re-zeroing,
   ! this loop's "+=" accumulates on top of that value instead of superseding
   ! it, silently doubling every zone-averaged, rezero=.TRUE. sheet diagnostic
   ! (e.g. aed_environ's ENV_air_temp, ENV_wind_speed, ENV_humidity, ...)
   ! whenever benthic zones are active. rezero=.FALSE. diagnostics (e.g. CGM
   ! running averages, ch4_ebb_dsfv) are deliberately left untouched here -
   ! this loop's accumulation across timesteps is how they are meant to persist.
   area = SUM(theZones(1:n_zones)%zarea)
   j = 0
   DO i=1,n_aed_vars
      IF ( aed_get_var(i, tvar) ) THEN
         IF ( tvar%var_type == V_DIAGNOSTIC .AND. tvar%sheet ) THEN
            j = j + 1
            IF ( tvar%rezero ) x_diag_hz(j) = 0.
         ENDIF
      ENDIF
   ENDDO
   DO zon=1, n_zones
      x_diag_hz = x_diag_hz + (z_diag_hz(:, zon) * (theZones(zon)%zarea/area))
   ENDDO
END SUBROUTINE copy_from_zone
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE GInitialTemp(m,depth,wv,topTemp,botTemp,nSPinUpDays,tNew) BIND(C, name="InitialTemp")
!-------------------------------------------------------------------------------
   USE aed_util, ONLY : InitialTemp
!ARGUMENTS
   INTEGER,intent(in)   :: m
   AED_REAL,intent(in)  :: wv,depth(0:m+1)
   AED_REAL,intent(in)  :: topTemp,botTemp,nSPinUpDays
   AED_REAL,intent(out) :: tNew(0:m+1)
!-------------------------------------------------------------------------------
!BEGIN
   CALL InitialTemp(m, depth, wv, topTemp, botTemp, nSPinUpDays, tNew)
END SUBROUTINE GInitialTemp
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE GInitialTempFV(n,depth,vwc_init,surfTemp,deepTemp,spinup_days,dt,tNew) &
                                                  BIND(C, name="InitialTempFV")
!-------------------------------------------------------------------------------
!# C-callable spin-up for the cell-centred soil-temperature model.  Reads the
!# soil thermal properties from the (Fortran-owned) glm_types globals and
!# returns the n interior-cell temperatures in tNew(1:n).
!-------------------------------------------------------------------------------
   USE aed_util, ONLY : InitialTempFV, soil_props_t
!ARGUMENTS
   INTEGER,intent(in)   :: n
   AED_REAL,intent(in)  :: depth(0:n+1), vwc_init, surfTemp, deepTemp, spinup_days, dt
   AED_REAL,intent(out) :: tNew(n)
!LOCALS
   TYPE(soil_props_t) :: props
!-------------------------------------------------------------------------------
!BEGIN
   props%k_mineral = sed_k_mineral ; props%k_water = sed_k_water
   props%c_mineral = sed_c_mineral ; props%c_water = sed_c_water ; props%c_air = sed_c_air
   props%bulk_density = sed_bulk_density ; props%mineral_density = sed_mineral_density
   props%porosity = sed_porosity
   CALL InitialTempFV(n, depth, vwc_init, surfTemp, deepTemp, spinup_days, dt, props, tNew)
END SUBROUTINE GInitialTempFV
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE ZSoilTemp(izone) BIND(C, name="zZSoilTemp")
!-------------------------------------------------------------------------------
   USE aed_util, ONLY : SoilTempFV, soil_props_t
!ARGUMENTS
   !# izone is a C pointer passed BY VALUE from glm_surface.c (ZSoilTemp(&theZones[z]));
   !# it MUST carry the VALUE attribute or Fortran treats it as pass-by-reference and
   !# dereferences the struct's first word as the address.
   TYPE(C_PTR),VALUE :: izone
!LOCALS
   TYPE(ZoneType),POINTER :: zone
   TYPE(SedLayerType),DIMENSION(:),POINTER :: layers
   INTEGER :: nsl, n, kk
   AED_REAL :: hf, deepTemp
   AED_REAL,DIMENSION(:),ALLOCATABLE :: sdepth, cellvwc, cellt
   TYPE(soil_props_t) :: props
!-------------------------------------------------------------------------------
!BEGIN
   CALL C_F_POINTER(izone, zone);
   nsl = zone%n_sed_layers
   n   = nsl - 2
   CALL C_F_POINTER(zone%c_layers, layers, [nsl]);
   !# Cell-centred layout: node depths d(0:n+1) = layers(1:nsl)%depth; the n
   !# interior cells (the prognostic state) live in slots 2..nsl-1; layers(1) is
   !# the surface boundary node and layers(nsl) the deep boundary node.
   ALLOCATE(sdepth(0:n+1), cellvwc(n), cellt(n))
   DO kk = 1, nsl
      sdepth(kk-1) = layers(kk)%depth
   ENDDO
   DO kk = 1, n
      cellvwc(kk) = layers(kk+1)%vwc
      cellt(kk)   = layers(kk+1)%temp
   ENDDO
   deepTemp = layers(nsl)%temp
   props%k_mineral = sed_k_mineral ; props%k_water = sed_k_water
   props%c_mineral = sed_c_mineral ; props%c_water = sed_c_water ; props%c_air = sed_c_air
   props%bulk_density = sed_bulk_density ; props%mineral_density = sed_mineral_density
   props%porosity = sed_porosity
   CALL SoilTempFV(n, sdepth, cellvwc, zone%ztemp, deepTemp, soil_dt, props, cellt, hf)
   zone%heatflux = hf
   DO kk = 1, n
      layers(kk+1)%temp = cellt(kk)
   ENDDO
   layers(1)%temp = zone%ztemp
   DEALLOCATE(sdepth, cellvwc, cellt)
END SUBROUTINE ZSoilTemp
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE GSoilTemp(m,depth,wv,topTemp,temp,heatflux) BIND(C, name="SoilTemp")
!-------------------------------------------------------------------------------
   USE aed_util, ONLY : SoilTemp
!ARGUMENTS
   INTEGER,intent(in) :: m
   AED_REAL,intent(in) :: depth(0:m+1), wv(m), topTemp
   AED_REAL,intent(inout) :: temp(m+1)
   AED_REAL, intent(out) :: heatflux

!-------------------------------------------------------------------------------
!BEGIN
   CALL SoilTemp(m, depth, wv, topTemp, temp, heatflux)
END SUBROUTINE GSoilTemp
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

!###############################################################################
AED_REAL FUNCTION zone_frac(lev, zon)
!-------------------------------------------------------------------------------
! Fraction of layer lev (heights lheights(lev-1)..lheights(lev), from 0 for the
! bottom layer) lying in zone zon. The bottom zone extends down to the bed and
! the top zone upward without limit, so the fractions of every layer sum to 1.
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: lev, zon
!
!LOCALS
   AED_REAL :: lo, hi, zlo, zhi
!
!-------------------------------------------------------------------------------
!BEGIN
   lo = 0. ; IF ( lev > 1 ) lo = lheights(lev-1)
   hi = lheights(lev)
   zlo = -1.0e30 ; IF ( zon > 1 ) zlo = zone_heights(zon-1)
   zhi =  1.0e30 ; IF ( zon < n_zones ) zhi = zone_heights(zon)
   IF ( hi <= lo ) THEN
      zone_frac = 0. ; IF ( hi > zlo .AND. hi <= zhi ) zone_frac = 1.
   ELSE
      zone_frac = MAX(0., MIN(hi, zhi) - MAX(lo, zlo)) / (hi - lo)
   ENDIF
END FUNCTION zone_frac
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
AED_REAL FUNCTION zone_bed_frac(lev, zon)
!-------------------------------------------------------------------------------
! Fraction of layer lev's BED (A(top) - A(bottom), the full area for the bottom
! layer) lying in zone zon. The bed between two heights is the difference of
! the plan areas there: the layer's own areas at its bottom and top, and the
! zone's ztop_area (same table interpolation as the layer areas) where a zone
! boundary cuts the layer. Shares of every layer sum to 1 and summed over the
! layers give each zone its true bed area. A layer with no bed (zero thickness
! or a vertical wall) is shared by height (zone_frac).
!-------------------------------------------------------------------------------
!ARGUMENTS
   INTEGER,INTENT(in) :: lev, zon
!
!LOCALS
   AED_REAL :: lo, hi, zlo, zhi, alo, ahi, a1, a2
!
!-------------------------------------------------------------------------------
!BEGIN
   lo = 0. ; IF ( lev > 1 ) lo = lheights(lev-1)
   hi = lheights(lev)
   alo = 0. ; IF ( lev > 1 ) alo = theLake(lev-1)%LayerArea
   ahi = theLake(lev)%LayerArea
   zlo = -1.0e30 ; IF ( zon > 1 ) zlo = zone_heights(zon-1)
   zhi =  1.0e30 ; IF ( zon < n_zones ) zhi = zone_heights(zon)
   IF ( hi <= lo .OR. ahi <= alo ) THEN
      zone_bed_frac = zone_frac(lev, zon)
   ELSEIF ( hi <= zlo .OR. lo >= zhi ) THEN
      zone_bed_frac = 0.
   ELSE
      a1 = alo ; IF ( lo < zlo ) a1 = theZones(zon-1)%ztop_area
      a2 = ahi ; IF ( hi > zhi ) a2 = theZones(zon)%ztop_area
      zone_bed_frac = MIN(1., MAX(0., (a2 - a1) / (ahi - alo)))
   ENDIF
END FUNCTION zone_bed_frac
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE glm_zones
