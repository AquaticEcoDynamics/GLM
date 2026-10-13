!###############################################################################
!#                                                                             #
!# glm_api_zones.F90                                                           #
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
!# Copyright 2024-2026 : The University of Western Australia                   #
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

#define _UNUSED(x) if (.FALSE.) print*,shape(x)

MODULE glm_api_zones

   USE ISO_C_BINDING

   USE glm_types
   USE aed_common
   USE aed_zones

   IMPLICIT NONE

   PRIVATE ! By default, make everything private

   AED_REAL,DIMENSION(:,:,:),ALLOCATABLE,TARGET :: z_cc      ! (n_vars, n_levs, nsed_zones)
   AED_REAL,DIMENSION(:,:),  ALLOCATABLE,TARGET :: z_cc_hz   ! (n_vars, nsed_zones)
   AED_REAL,DIMENSION(:,:,:),ALLOCATABLE,TARGET :: z_diag    ! (n_vars, n_levs, nsed_zones)
   AED_REAL,DIMENSION(:,:),  ALLOCATABLE,TARGET :: z_diag_hz ! (n_vars, nsed_zones)

   AED_REAL,DIMENSION(:),POINTER :: lheights

   !# Layer-zone contact shares for the current host step, computed once by
   !# set_zone_shares (from api_calc_zone_areas) and used by every mapping:
   !#   zshare(lev,zon)  fraction of layer lev's bed lying in zone zon (zone_bed_frac)
   !#   zhshare(lev,zon) fraction of layer lev's height lying in zone zon (zone_frac),
   !#                    the volume weights for a wet zone with no bed of its own
   !#   zwbed(zon), zwvol(zon)  the zone's wet bed area and water volume from them
   !# aedZones(zon)%z_share points at zshare(:,zon), so libaed-api's glm_benthics
   !# delivers the zone fluxes with the same shares.
   AED_REAL,DIMENSION(:,:),ALLOCATABLE,TARGET :: zshare, zhshare
   AED_REAL,DIMENSION(:),  ALLOCATABLE :: zwbed, zwvol
   LOGICAL :: shares_set = .FALSE.
   !# indices (into x_diag) of the zone-averaged (zavg) non-sheet diagnostics
   INTEGER,DIMENSION(:),ALLOCATABLE :: zavg_idx

   AED_REAL,DIMENSION(:),POINTER :: zone_heights
   INTEGER :: n_aed_vars
   INTEGER :: nvars, nbenv, nvdiag, nvdiag_hz

   PUBLIC n_zones, zone_heights, lheights
   PUBLIC api_set_glm_zones, api_copy_from_zone, api_copy_to_zone, api_calc_zone_areas

   PUBLIC z_cc, z_cc_hz, z_diag, z_diag_hz, theZones, theLake

CONTAINS


!###############################################################################
SUBROUTINE api_set_glm_zones(numVars, numBenV, numDiagV, numDiagHzV, nAEDvars) &
                                               BIND(C, name="api_set_glm_zones")
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER,INTENT(in) :: numVars, numBenV, numDiagV, numDiagHzV, nAEDvars
!
!LOCALS
   INTEGER :: zon

   PROCEDURE(copy_to_zone_t),POINTER    :: copy_to
   PROCEDURE(copy_from_zone_t),POINTER  :: copy_from
   PROCEDURE(calc_zone_areas_t),POINTER :: calc_areas
!
!-------------------------------------------------------------------------------
!BEGIN
   nvars = numVars
   nbenv = numBenV
   nvdiag = numDiagV
   nvdiag_hz = numDiagHzV
   n_aed_vars = nAEDvars

   lheights => theLake%Height

   ALLOCATE(z_cc(numVars, MaxLayers, n_zones))    ; z_cc = 0.
   ALLOCATE(z_cc_hz(numBenV, n_zones))            ; z_cc_hz = 0.
   ALLOCATE(z_diag(numDiagV, MaxLayers, n_zones)) ; z_diag = 0.
   ALLOCATE(z_diag_hz(numDiagHzV, n_zones+1))     ; z_diag_hz = 0.

   CALL aed_init_zones(n_zones, 1, z_cc, z_cc_hz, z_diag, z_diag_hz)

   copy_to => api_copy_to_zone
   copy_from => api_copy_from_zone
   calc_areas => api_calc_zone_areas
   CALL api_set_zone_funcs(copy_to, copy_from, calc_areas)

   ALLOCATE(zshare(MaxLayers, n_zones))  ; zshare = 0.
   ALLOCATE(zhshare(MaxLayers, n_zones)) ; zhshare = 0.
   ALLOCATE(zwbed(n_zones), zwvol(n_zones)) ; zwbed = 0. ; zwvol = 0.

   DO zon=1,n_zones
      aedZones(zon)%z_env%z_area = 0.
      aedZones(zon)%z_env%z_height = theZones(zon)%zheight
      aedZones(zon)%z_env%z_top_area = theZones(zon)%ztop_area
      aedZones(zon)%z_share => zshare(:,zon)
   ENDDO
   theZones%zarea = 0.
   zone_heights => theZones%zheight
END SUBROUTINE api_set_glm_zones
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE set_zone_shares(aedZones, n_zones, wheights, wlev)
!-------------------------------------------------------------------------------
! The layer-zone contact shares for the current layer structure (see zshare).
! Layers move only between host steps, so this runs once per step, from
! api_calc_zone_areas; the routines that map between layers and zones read the
! stored shares instead of recomputing them.
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(api_zone_t),DIMENSION(:),INTENT(in) :: aedZones
   INTEGER,INTENT(in) :: n_zones, wlev
   AED_REAL,DIMENSION(:),POINTER,INTENT(in) :: wheights
!LOCALS
   INTEGER  :: lev, zon
   AED_REAL :: bed
!-------------------------------------------------------------------------------
!BEGIN
   zshare = 0. ; zhshare = 0.     !# rows above wlev stay 0 (no layer there)
   zwbed = 0. ; zwvol = 0.
   DO lev=1,wlev
      bed = theLake(lev)%LayerArea ; IF ( lev > 1 ) bed = bed - theLake(lev-1)%LayerArea
      DO zon=1,n_zones
         zshare(lev,zon)  = zone_bed_frac(aedZones, n_zones, wheights, lev, zon)
         zhshare(lev,zon) = zone_frac(aedZones, n_zones, wheights, lev, zon)
         zwbed(zon) = zwbed(zon) + zshare(lev,zon) * bed
         zwvol(zon) = zwvol(zon) + zhshare(lev,zon) * theLake(lev)%LayerVol
      ENDDO
   ENDDO
   shares_set = .TRUE.
END SUBROUTINE set_zone_shares
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE api_calc_zone_areas(aedZones, n_zones, areas, wheights, wlev)
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(api_zone_t),DIMENSION(:),INTENT(inout) :: aedZones
   INTEGER,INTENT(in) :: n_zones
!
   AED_REAL,DIMENSION(:),POINTER,INTENT(in) :: areas
   AED_REAL,DIMENSION(:),POINTER,INTENT(in) :: wheights
   INTEGER,INTENT(in) :: wlev
!
!LOCALS
   INTEGER  :: lev, zon
   LOGICAL  :: w_zones
   AED_REAL :: surf
!
!-------------------------------------------------------------------------------
!BEGIN
   DO zon=1, n_zones
      aedZones(zon)%z_env%z_area = 0.  ; aedZones(zon)%z_env%z_pc_wet = 0.
   ENDDO
   w_zones = .FALSE.
   surf = wheights(wlev)

   !# Zone areas: each layer's bed (areas() is the INCREMENTAL layer area) is
   !# shared between zones by the part of that bed lying in each zone
   !# (zone_bed_frac), the same shares glm_benthics uses to map zone fluxes onto
   !# layers. This handles a layer whose bottom sits exactly on a boundary, a
   !# bottom layer that straddles a boundary and a layer spanning several zones.
   CALL set_zone_shares(aedZones, n_zones, wheights, wlev)
   DO lev=1, wlev
      DO zon=1, n_zones
         aedZones(zon)%z_env%z_area = aedZones(zon)%z_env%z_area + areas(lev) * zshare(lev,zon)
      ENDDO
   ENDDO

   !# Wet fraction of each zone (unchanged): zones below the surface are fully
   !# wet, the first zone reaching above it is partly wet, those above are dry.
   DO zon=1, n_zones
      IF ( aedZones(zon)%z_env%z_height <= surf ) THEN
         aedZones(zon)%z_env%z_pc_wet = 1.0
      ELSEIF ( .NOT. w_zones ) THEN
         w_zones = .TRUE.
         IF ( zon > 1 ) THEN
            aedZones(zon)%z_env%z_pc_wet = surf / (aedZones(zon)%z_env%z_height - aedZones(zon-1)%z_env%z_height)
         ELSE
            aedZones(zon)%z_env%z_pc_wet = surf / aedZones(zon)%z_env%z_height
         ENDIF
      ENDIF
   ENDDO

   DO zon=1, n_zones
      aedZones(zon)%z_env%z_pres = -aedZones(zon)%z_env%z_height
      aedZones(zon)%z_env%z_col_area = aedZones(zon)%z_env%z_area
      theZones(zon)%zarea = aedZones(zon)%z_env%z_area     ! for netCDF 'zarea' output
   ENDDO
END SUBROUTINE api_calc_zone_areas
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE api_copy_to_zone(aedZones, n_zones, wheights, x_cc, x_cc_hz, x_diag, x_diag_hz, wlev)
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(api_zone_t),DIMENSION(:),INTENT(inout) :: aedZones
   INTEGER,INTENT(in) :: n_zones
!
   AED_REAL,DIMENSION(:),  POINTER,INTENT(in) :: wheights
   AED_REAL,DIMENSION(:,:),POINTER,INTENT(in) :: x_cc
   AED_REAL,DIMENSION(:),  POINTER,INTENT(in) :: x_cc_hz
   AED_REAL,DIMENSION(:,:),POINTER,INTENT(in) :: x_diag
   AED_REAL,DIMENSION(:),  POINTER,INTENT(in) :: x_diag_hz
!
   INTEGER,INTENT(in) :: wlev
!
!LOCALS
   INTEGER  :: zon, lev, a_zones
   AED_REAL :: surf, f, bed, w
   LOGICAL  :: w_zones
   INTEGER  :: i, sd
   TYPE(aed_variable_t),POINTER :: tv
!
!-------------------------------------------------------------------------------
!BEGIN
   _UNUSED(x_cc_hz)
   _UNUSED(x_diag)
   _UNUSED(x_diag_hz)

   DO zon=1,n_zones
      z_cc(1:nvars,:,zon) = 0.
      z_diag(:,:,zon) = 0.
      !# z_diag_hz is the HOME of the per-zone benthic diagnostics (bound to
      !# aedZones(zon)%z_cc_diag_hz), so clear it PER-VARIABLE: diagnostics
      !# flagged rezero=.FALSE. (e.g. the CGM running averages cgm_tavg/lavg/savg
      !# and slough counters) must persist across timesteps; the rest are cleared.
      sd = 0
      DO i=1,n_aed_vars
         IF ( aed_get_var(i, tv) ) THEN
            IF ( tv%var_type == V_DIAGNOSTIC .AND. tv%sheet ) THEN
               sd = sd + 1
               IF ( tv%rezero ) z_diag_hz(sd,zon) = 0.
            ENDIF
         ENDIF
      ENDDO

      aedZones(zon)%z_env%z_temp = 0.
      aedZones(zon)%z_env%z_salt = 0.
      aedZones(zon)%z_env%z_rho = 0.
      aedZones(zon)%z_env%z_rad = 0.
      aedZones(zon)%z_env%z_extc = 0.
      aedZones(zon)%z_env%z_layer_stress = 0.
      aedZones(zon)%z_env%z_tss = 0.
      aedZones(zon)%z_env%z_par = 0.
      aedZones(zon)%z_env%z_nir = 0.
      aedZones(zon)%z_env%z_uva = 0.
      aedZones(zon)%z_env%z_uvb = 0.
      aedZones(zon)%z_env%z_sed_zones = zon
      aedZones(zon)%z_env%z_vel = 0.
   ENDDO

   !# Zone-averaged water conditions (what each zone's sediment sees). Each
   !# layer is weighted by the bed area it touches within the zone
   !# (zone_bed_frac * bed area), the same shares used for the zone areas and
   !# to map zone fluxes onto layers. A wet zone with no bed of its own
   !# (vertical wall) falls back to water-volume weights, shared by height
   !# (zone_frac). Zones entirely above the water have no weight and are left
   !# at zero, as before.
   w_zones = .FALSE.
   IF ( .NOT. shares_set ) CALL set_zone_shares(aedZones, n_zones, wheights, wlev)

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

         z_cc(1:nvars,1,zon) = z_cc(1:nvars,1,zon) + w * x_cc(1:nvars,lev)  ! level=1 to match aedZones(zon)%z_cc => z_cc(:,1,:) pointer
         z_diag(:,1,zon)     = z_diag(:,1,zon)     + w * x_diag(:,lev)      ! level=1 matches library wiring z_cc_diag => z_diag(:,1,:)

         aedZones(zon)%z_env%z_temp         = aedZones(zon)%z_env%z_temp + w * theLake(lev)%Temp
         aedZones(zon)%z_env%z_salt         = aedZones(zon)%z_env%z_salt + w * theLake(lev)%Salinity
         aedZones(zon)%z_env%z_rho          = aedZones(zon)%z_env%z_rho  + w * theLake(lev)%Density
         aedZones(zon)%z_env%z_rad          = aedZones(zon)%z_env%z_rad  + w * theLake(lev)%Light
         aedZones(zon)%z_env%z_par          = aedZones(zon)%z_env%z_par  + w * theLake(lev)%Light*0.45
         aedZones(zon)%z_env%z_vel          = aedZones(zon)%z_env%z_vel  + w * theLake(lev)%Umean
         aedZones(zon)%z_env%z_extc         = aedZones(zon)%z_env%z_extc + w * theLake(lev)%ExtcCoefSW
         aedZones(zon)%z_env%z_layer_stress = aedZones(zon)%z_env%z_layer_stress + w * theLake(lev)%LayerStress
      ENDDO
   ENDDO

   DO zon=1,n_zones
      aedZones(zon)%z_env%z_dz = 0.
   ENDDO

   surf = wheights(wlev)

   ! Populate lake-wide sheet globals into each zone's env struct.
   ! In glm_aed.F90 these are wired as direct pointer-to-pointer; here the zone
   ! struct holds scalar fields so we copy the live values each timestep.
   DO zon=1,n_zones
      aedZones(zon)%z_env%z_col_depth = surf
      aedZones(zon)%z_env%z_air_pres  = MetData%AirPres
      aedZones(zon)%z_env%z_air_temp  = MetData%AirTemp
      aedZones(zon)%z_env%z_wind      = MetData%WindSpeed
      aedZones(zon)%z_env%z_humidity  = MetData%RelHum
      aedZones(zon)%z_env%z_rain      = MetData%Rain
      aedZones(zon)%z_env%z_evap      = SurfData%Evap
      aedZones(zon)%z_env%z_I_0       = MetData%ShortWave
      aedZones(zon)%z_env%z_longwave  = MetData%LongWave
   ENDDO

   IF ( surf > aedZones(1)%z_env%z_height ) THEN
      aedZones(1)%z_env%z_depth = aedZones(1)%z_env%z_height
      aedZones(1)%z_env%z_dz = aedZones(1)%z_env%z_height
   ELSE
      aedZones(1)%z_env%z_depth = surf
      aedZones(1)%z_env%z_dz = surf
      w_zones = .TRUE.
   ENDIF

   DO zon=2,a_zones
      aedZones(zon)%z_env%z_sed_zones = zon
      IF ( .NOT. w_zones ) THEN
          IF ( surf > aedZones(zon)%z_env%z_height ) THEN
             aedZones(zon)%z_env%z_depth = aedZones(zon)%z_env%z_height
             aedZones(zon)%z_env%z_dz =  aedZones(zon)%z_env%z_height - aedZones(zon-1)%z_env%z_height
          ELSE
             aedZones(zon)%z_env%z_depth = surf
             aedZones(zon)%z_env%z_dz = surf - aedZones(zon-1)%z_env%z_height
             w_zones = .TRUE.
          ENDIF
      ELSE
         aedZones(zon)%z_env%z_depth = surf
         aedZones(zon)%z_env%z_dz = surf - aedZones(zon-1)%z_env%z_height
      ENDIF
   ENDDO

END SUBROUTINE api_copy_to_zone
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE api_copy_from_zone(aedZones, n_zones, wheights, x_cc, x_cc_hz, x_diag, x_diag_hz, wlev)
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(api_zone_t),DIMENSION(:),INTENT(in) :: aedZones
   INTEGER,INTENT(in) :: n_zones
!
   AED_REAL,DIMENSION(:),  POINTER,INTENT(in) :: wheights
   AED_REAL,DIMENSION(:,:),POINTER,INTENT(inout) :: x_cc
   AED_REAL,DIMENSION(:),  POINTER,INTENT(inout) :: x_cc_hz
   AED_REAL,DIMENSION(:,:),POINTER,INTENT(inout) :: x_diag
   AED_REAL,DIMENSION(:),  POINTER,INTENT(inout) :: x_diag_hz
!
   INTEGER,INTENT(in) :: wlev
!
!LOCALS
   INTEGER  :: zon, lev, i, j, k
   AED_REAL :: scale, area
   TYPE(aed_variable_t),POINTER :: tvar
!
!-------------------------------------------------------------------------------
!BEGIN
   _UNUSED(x_cc)

   ! Set water layer zone-averaged diagnostics (zavg == true) from the zones,
   ! each zone weighted by its share of the layer's bed (zone_bed_frac, the
   ! same shares used for zone areas and flux mapping).
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
   IF ( .NOT. shares_set ) CALL set_zone_shares(aedZones, n_zones, wheights, wlev)

   DO lev=1,wlev
      DO k=1,SIZE(zavg_idx)
         j = zavg_idx(k)
         x_diag(j,lev) = 0.
         DO zon=1,n_zones
            scale = zshare(lev,zon)
            IF ( scale > 0. ) x_diag(j,lev) = x_diag(j,lev) + z_diag(j,1,zon) * scale   ! level=1 matches library z_cc_diag => z_diag(:,1,:)
         ENDDO
      ENDDO
   ENDDO

   ! Single-value (lake-wide) benthic state output: mean of the zones weighted
   ! by bed area (it previously ended up holding the bottom zone's values).
   area = 0.
   DO zon=1,n_zones
      area = area + aedZones(zon)%z_env%z_area
   ENDDO
   IF ( area > 0. ) THEN
      x_cc_hz(:) = 0.
      DO zon=1,n_zones
         x_cc_hz(:) = x_cc_hz(:) + z_cc_hz(:,zon) * (aedZones(zon)%z_env%z_area/area)
      ENDDO
   ENDIF

   ! Set the normal sheet diagnostics to the mean of the zone, weighted by area.
   ! Re-zero only the rezero=.TRUE. sheet diagnostics here: aed_run_model calls
   ! aed_calculate_riparian once, unconditionally, for the whole (non-zoned)
   ! column before this per-zone aggregation runs; that call writes a correct
   ! single-pass value into the same sheet-diagnostic slots via a plain
   ! assignment. Without re-zeroing, this loop's "+=" accumulates on top of
   ! that value instead of superseding it, silently doubling every
   ! zone-averaged, rezero=.TRUE. sheet diagnostic (e.g. aed_environ's
   ! ENV_air_temp, ENV_wind_speed, ENV_humidity, ...) whenever benthic zones
   ! are active. rezero=.FALSE. sheet diagnostics are deliberately left
   ! untouched here - this loop's accumulation across timesteps is how they
   ! are meant to persist. (ch4_ebb_dsfv is NOT an example of this: it is a
   ! per-layer diagnostic, not a sheet one, so tvar%sheet is .FALSE. for it
   ! and this loop never touches it at all.)
   area = 0.
   DO zon=1,n_zones
      area = area + aedZones(zon)%z_env%z_area
   ENDDO
   j = 0
   DO i=1,n_aed_vars
      IF ( aed_get_var(i, tvar) ) THEN
         IF ( tvar%var_type == V_DIAGNOSTIC .AND. tvar%sheet ) THEN
            j = j + 1
            IF ( tvar%rezero ) x_diag_hz(j) = 0.
         ENDIF
      ENDIF
   ENDDO
   DO zon=1,n_zones
      x_diag_hz = x_diag_hz + (z_diag_hz(:,zon) * (aedZones(zon)%z_env%z_area/area))
   ENDDO

END SUBROUTINE api_copy_from_zone
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
AED_REAL FUNCTION zone_frac(aedZones, n_zones, wheights, lev, zon)
!-------------------------------------------------------------------------------
! Fraction of layer lev (heights wheights(lev-1)..wheights(lev), from 0 for the
! bottom layer) lying in zone zon. The bottom zone extends down to the bed and
! the top zone upward without limit, so the fractions of every layer sum to 1.
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(api_zone_t),DIMENSION(:),INTENT(in) :: aedZones
   INTEGER,INTENT(in) :: n_zones, lev, zon
   AED_REAL,DIMENSION(:),POINTER,INTENT(in) :: wheights
!
!LOCALS
   AED_REAL :: lo, hi, zlo, zhi
!
!-------------------------------------------------------------------------------
!BEGIN
   lo = 0. ; IF ( lev > 1 ) lo = wheights(lev-1)
   hi = wheights(lev)
   zlo = -1.0e30 ; IF ( zon > 1 ) zlo = aedZones(zon-1)%z_env%z_height
   zhi =  1.0e30 ; IF ( zon < n_zones ) zhi = aedZones(zon)%z_env%z_height
   IF ( hi <= lo ) THEN
      zone_frac = 0. ; IF ( hi > zlo .AND. hi <= zhi ) zone_frac = 1.
   ELSE
      zone_frac = MAX(0., MIN(hi, zhi) - MAX(lo, zlo)) / (hi - lo)
   ENDIF
END FUNCTION zone_frac
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
AED_REAL FUNCTION zone_bed_frac(aedZones, n_zones, wheights, lev, zon)
!-------------------------------------------------------------------------------
! Fraction of layer lev's BED (A(top) - A(bottom), the full area for the bottom
! layer) lying in zone zon. The bed between two heights is the difference of
! the plan areas there: the layer's own areas at its bottom and top, and the
! zone's z_top_area (same table interpolation as the layer areas) where a zone
! boundary cuts the layer. Shares of every layer sum to 1 and summed over the
! layers give each zone its true bed area. A layer with no bed (zero thickness
! or a vertical wall) is shared by height (zone_frac).
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(api_zone_t),DIMENSION(:),INTENT(in) :: aedZones
   INTEGER,INTENT(in) :: n_zones, lev, zon
   AED_REAL,DIMENSION(:),POINTER,INTENT(in) :: wheights
!
!LOCALS
   AED_REAL :: lo, hi, zlo, zhi, alo, ahi, a1, a2
!
!-------------------------------------------------------------------------------
!BEGIN
   lo = 0. ; IF ( lev > 1 ) lo = wheights(lev-1)
   hi = wheights(lev)
   alo = 0. ; IF ( lev > 1 ) alo = theLake(lev-1)%LayerArea
   ahi = theLake(lev)%LayerArea
   zlo = -1.0e30 ; IF ( zon > 1 ) zlo = aedZones(zon-1)%z_env%z_height
   zhi =  1.0e30 ; IF ( zon < n_zones ) zhi = aedZones(zon)%z_env%z_height
   IF ( hi <= lo .OR. ahi <= alo ) THEN
      zone_bed_frac = zone_frac(aedZones, n_zones, wheights, lev, zon)
   ELSEIF ( hi <= zlo .OR. lo >= zhi ) THEN
      zone_bed_frac = 0.
   ELSE
      a1 = alo ; IF ( lo < zlo ) a1 = aedZones(zon-1)%z_env%z_top_area
      a2 = ahi ; IF ( hi > zhi ) a2 = aedZones(zon)%z_env%z_top_area
      zone_bed_frac = MIN(1., MAX(0., (a2 - a1) / (ahi - alo)))
   ENDIF
END FUNCTION zone_bed_frac


END MODULE glm_api_zones
