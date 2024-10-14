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

#undef MISVAL
#ifndef _FORTRAN_SOURCE_
#define _FORTRAN_SOURCE_ 1
#endif

#include "glm.h"

MODULE glm_api_zones

   USE ISO_C_BINDING

   USE glm_types
   USE aed_zones

   IMPLICIT NONE

   PRIVATE ! By default, make everything private

   AED_REAL,DIMENSION(:,:,:),ALLOCATABLE,TARGET :: z_cc      ! (nsed_zones, n_levs, n_vars)
   AED_REAL,DIMENSION(:,:),  ALLOCATABLE,TARGET :: z_cc_hz   ! (nsed_zones, n_vars)
   AED_REAL,DIMENSION(:,:,:),ALLOCATABLE,TARGET :: z_diag    ! (nsed_zones, n_levs, n_vars)
   AED_REAL,DIMENSION(:,:),  ALLOCATABLE,TARGET :: z_diag_hz ! (nsed_zones, n_vars)

   AED_REAL,DIMENSION(:),POINTER :: lheights

   INTEGER :: nvars, nbenv, nvdiag, nvdiag_hz

   PUBLIC n_zones, lheights
   PUBLIC api_set_glm_zones, api_copy_from_zone, api_copy_to_zone, api_calc_zone_areas

   PUBLIC z_cc, z_cc_hz, z_diag, z_diag_hz, theZones, theLake

CONTAINS


!###############################################################################
SUBROUTINE api_set_glm_zones(numVars, numBenV, numDiagV, numDiagHzV)           &
                                               BIND(C, name="api_set_glm_zones")
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER,INTENT(in) :: numVars, numBenV, numDiagV, numDiagHzV
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

   lheights => theLake%Height

   ALLOCATE(z_cc(n_zones, MaxLayers, numVars+numBenV))  ; z_cc = 0.
   ALLOCATE(z_cc_hz(n_zones, numVars+numBenV))          ; z_cc_hz = 0.
   ALLOCATE(z_diag(n_zones, MaxLayers, numDiagV))       ; z_diag = 0.
   ALLOCATE(z_diag_hz(n_zones+1, numDiagHzV))           ; z_diag_hz = 0.

   CALL aed_init_zones(n_zones, 1, z_cc, z_cc_hz, z_diag, z_diag_hz)

   copy_to => api_copy_to_zone
   copy_from => api_copy_from_zone
   calc_areas => api_calc_zone_areas
   CALL api_set_zone_funcs(copy_to, copy_from, calc_areas)

   DO zon=1,n_zones
      aedZones(zon)%z_env(1)%z_area = 0.
      aedZones(zon)%z_env(1)%z_height = theZones(zon)%zheight
   ENDDO
END SUBROUTINE api_set_glm_zones
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

   zon = 1
   aedZones(1)%z_env(1)%z_area = areas(1)
   DO lev=2, wlev
      IF ( wheights(lev) > aedZones(zon)%z_env(1)%z_height ) zon = zon + 1

      aedZones(zon)%z_env(1)%z_area = aedZones(zon)%z_env(1)%z_area + areas(lev) - areas(lev-1)

      IF ( aedZones(zon)%z_env(1)%z_height > surf ) THEN
         IF (.NOT. w_zones) THEN
            w_zones = .TRUE.
            IF ( zon > 1 ) THEN
               aedZones(zon)%z_env(1)%z_pc_wet = surf / (aedZones(zon)%z_env(1)%z_height - aedZones(zon-1)%z_env(1)%z_height)
            ELSE
               aedZones(zon)%z_env(1)%z_pc_wet = surf / aedZones(zon)%z_env(1)%z_height
            ENDIF
         ENDIF
      ELSE
         aedZones(zon)%z_env(1)%z_pc_wet = 1.0
      ENDIF
   ENDDO

   DO zon=1, n_zones
      aedZones(zon)%z_env(1)%z_pres = -aedZones(zon)%z_env(1)%z_height
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
   AED_REAL :: surf
   INTEGER  :: zcount(n_zones)
   LOGICAL  :: w_zones
!
!-------------------------------------------------------------------------------
!BEGIN
   DO zon=1,n_zones
      z_cc(zon,:,1:nvars) = 0.
      z_diag(zon,:,:) = 0.
      z_diag_hz(zon,:) = 0.

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
      aedZones(zon)%z_env%z_sed_zones = 1.
      aedZones(zon)%z_env%z_vel = 0.
   ENDDO

   a_zones = 1
   zcount = 0
   w_zones = .FALSE.
   zon = 1
   DO lev=1,wlev
      IF ( lev > 1 .AND. wheights(lev) > aedZones(zon)%z_env(1)%z_height ) THEN
         zon = zon + 1
         IF (zon > n_zones) STOP 'Water level height is higher than highest zone height'
         aedZones(zon)%z_env(1)%z_sed_zones = zon
      ENDIF

      ! Pelagic variables need zonifying, in case benthic people want them
      ! Note that this should not be done on benthic vars because it will
      ! introduce errors in z_cc (split layers)
      ! Ideally this average would be based on volume weighting

      z_cc(zon,lev,1:nvars) = z_cc(zon,lev,1:nvars) + x_cc(lev,1:nvars)
      z_diag(zon,lev,:)     = z_diag(zon,lev,:) + x_diag(lev,:)
      z_diag_hz(zon,:)      = z_diag_hz(zon,:) + x_diag_hz(:)

      aedZones(zon)%z_env(1)%z_temp         = aedZones(zon)%z_env(1)%z_temp + theLake(lev)%Temp
      aedZones(zon)%z_env(1)%z_salt         = aedZones(zon)%z_env(1)%z_salt + theLake(lev)%Salinity
      aedZones(zon)%z_env(1)%z_rho          = aedZones(zon)%z_env(1)%z_rho  + theLake(lev)%Density
      aedZones(zon)%z_env(1)%z_rad          = aedZones(zon)%z_env(1)%z_rad  + theLake(lev)%Light
      aedZones(zon)%z_env(1)%z_vel          = aedZones(zon)%z_env(1)%z_vel  + theLake(lev)%Umean
      aedZones(zon)%z_env(1)%z_extc         = aedZones(zon)%z_env(1)%z_extc + theLake(lev)%ExtcCoefSW
      aedZones(zon)%z_env(1)%z_layer_stress = aedZones(zon)%z_env(1)%z_layer_stress + theLake(lev)%LayerStress

      zcount(zon) = zcount(zon) + 1
   ENDDO
   a_zones = zon

   DO zon=1,a_zones
      z_cc(zon,:,1:nvars) = z_cc(zon,:,1:nvars)/zcount(zon)
      z_diag(zon,:,:)     = z_diag(zon,:,:)/zcount(zon)
      z_diag_hz(zon,:)  = z_diag_hz(zon,:)/zcount(zon)
   ENDDO

   DO zon=1,a_zones
      IF (zcount(zon) /= 0) THEN
         aedZones(zon)%z_env(1)%z_temp         = aedZones(zon)%z_env(1)%z_temp / zcount(zon)
         aedZones(zon)%z_env(1)%z_salt         = aedZones(zon)%z_env(1)%z_salt / zcount(zon)
         aedZones(zon)%z_env(1)%z_rho          = aedZones(zon)%z_env(1)%z_rho  / zcount(zon)
         aedZones(zon)%z_env(1)%z_rad          = aedZones(zon)%z_env(1)%z_rad  / zcount(zon)
         aedZones(zon)%z_env(1)%z_vel          = aedZones(zon)%z_env(1)%z_vel  / zcount(zon)
         aedZones(zon)%z_env(1)%z_extc         = aedZones(zon)%z_env(1)%z_extc / zcount(zon)
         aedZones(zon)%z_env(1)%z_layer_stress = aedZones(zon)%z_env(1)%z_layer_stress / zcount(zon)
      ELSE
         aedZones(zon)%z_env(1)%z_temp         = 0.
         aedZones(zon)%z_env(1)%z_salt         = 0.
         aedZones(zon)%z_env(1)%z_rho          = 0.
         aedZones(zon)%z_env(1)%z_rad          = 0.
         aedZones(zon)%z_env(1)%z_vel          = 0.
         aedZones(zon)%z_env(1)%z_extc         = 0.
         aedZones(zon)%z_env(1)%z_layer_stress = 0.
      ENDIF

      aedZones(zon)%z_env(1)%z_dz = 0.
   ENDDO

   surf = wheights(wlev)
   IF ( surf > aedZones(1)%z_env(1)%z_height ) THEN
      aedZones(1)%z_env(1)%z_depth = aedZones(1)%z_env(1)%z_height
      aedZones(1)%z_env(1)%z_dz = aedZones(1)%z_env(1)%z_height
   ELSE
      aedZones(1)%z_env(1)%z_depth = surf
      aedZones(1)%z_env(1)%z_dz = surf
      w_zones = .TRUE.
   ENDIF

   DO zon=2,a_zones
      aedZones(zon)%z_env(1)%z_sed_zones = zon
      IF ( .NOT. w_zones ) THEN
          IF ( surf > aedZones(zon)%z_env(1)%z_height ) THEN
             aedZones(zon)%z_env(1)%z_depth = aedZones(zon)%z_env(1)%z_height
             aedZones(zon)%z_env(1)%z_dz =  aedZones(zon)%z_env(1)%z_height - aedZones(zon-1)%z_env(1)%z_height
          ELSE
             aedZones(zon)%z_env(1)%z_depth = surf
             aedZones(zon)%z_env(1)%z_dz = surf - aedZones(zon-1)%z_env(1)%z_height
             w_zones = .TRUE.
          ENDIF
      ELSE
         aedZones(zon)%z_env(1)%z_depth = surf
         aedZones(zon)%z_env(1)%z_dz = surf - aedZones(zon-1)%z_env(1)%z_height
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
   INTEGER  :: zon, lev, v_start, v_end
   AED_REAL :: scale, area
   LOGICAL  :: splitZone
!
!-------------------------------------------------------------------------------
!BEGIN
   v_start = nvars+1 ; v_end = nvars+nbenv

   zon = n_zones
   DO lev=wlev,1,-1
      IF ( zon .GT. 1 ) THEN
         IF (lev .GT. 1) THEN
            splitZone = wheights(lev-1) < aedZones(zon-1)%z_env(1)%z_height
         ELSE
            splitZone = 0.0 < aedZones(zon-1)%z_env(1)%z_height
         ENDIF
      ELSE
         splitZone = .FALSE.
      ENDIF

      IF (splitZone) THEN
         IF (lev .GT. 1) THEN
            scale = (aedZones(zon-1)%z_env(1)%z_height - wheights(lev-1)) / (wheights(lev) - wheights(lev-1))
         ELSE
            scale = (aedZones(zon-1)%z_env(1)%z_height - 0.0) / (wheights(lev) - 0.0)
         ENDIF

         WHERE(z_diag(zon,lev,:) /= 0.) &
            x_diag(lev,:) = z_diag(zon,lev,:) * scale
         x_cc(lev,v_start:v_end) = z_cc(zon,lev,v_start:v_end) * scale

         zon = zon - 1

         WHERE(z_diag(zon,lev,:) /= 0.) &
            x_diag(lev,:) = x_diag(lev,:) + (z_diag(zon,lev,:) * (1.0 - scale))
         x_cc(lev,v_start:v_end) = x_cc(lev,v_start:v_end) + &
                                   z_cc(zon,lev,v_start:v_end) * (1.0 - scale)
      ELSE
         WHERE(z_diag(zon,lev,:) /= 0.) &
            x_diag(lev,:) = z_diag(zon,lev,:)
         x_cc(lev,v_start:v_end) = z_cc(zon,lev,v_start:v_end)
      ENDIF
   ENDDO

   ! Set the normal sheet diagnostics to the mean of the zone, weighted by area
   area = 0.
   DO zon=1,n_zones
      area = area + aedZones(zon)%z_env(1)%z_area
   ENDDO
   DO zon=1,n_zones
      x_diag_hz = x_diag_hz + (z_diag_hz(zon,:) * (aedZones(zon)%z_env(1)%z_area/area))
   ENDDO
END SUBROUTINE api_copy_from_zone
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE glm_api_zones
