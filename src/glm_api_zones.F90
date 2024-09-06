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
!# Copyright 2013 - 2024 -  The University of Western Australia                #
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

   USE aed_common

   USE glm_types
   USE glm_zones

   USE aed_api
   USE aed_zones

   IMPLICIT NONE

   PRIVATE ! By default, make everything private

   INTEGER :: nvars, nbenv

   PUBLIC api_set_glm_zones

CONTAINS

!###############################################################################
SUBROUTINE api_set_glm_zones(numVars, numBenV, numDiagV, numDiagHzV)           &
                                               BIND(C, name="api_set_glm_zones")
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER,INTENT(in) :: numVars, numBenV, numDiagV, numDiagHzV
!
!LOCALS
   INTEGER :: i

   PROCEDURE(copy_to_zone_t),POINTER    :: copy_to
   PROCEDURE(copy_from_zone_t),POINTER  :: copy_from
   PROCEDURE(calc_zone_areas_t),POINTER :: calc_areas
!
!-------------------------------------------------------------------------------
!BEGIN
   nvars = numVars
   nbenv = numBenV

   CALL wq_set_glm_zones(numVars, numBenV, numDiagV, numDiagHzV)

   CALL aed_init_zones(n_zones, 1, z_cc, z_cc_hz, z_diag, z_diag_hz)

   copy_to => api_copy_to_zone
   copy_from => api_copy_from_zone
   calc_areas => api_calc_zone_areas
   CALL api_set_zone_funcs(copy_to, copy_from, calc_areas)

   DO i=1,n_zones
      aedZones(i)%z_heights(1) = theZones(i)%zheight
   ENDDO
END SUBROUTINE api_set_glm_zones
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE api_calc_zone_areas(aedZones, n_zones, areas, heights, wlev)
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(api_zone_t),DIMENSION(:),INTENT(inout) :: aedZones
   INTEGER,INTENT(in) :: n_zones
!
   AED_REAL,DIMENSION(:),POINTER,INTENT(in) :: areas
   AED_REAL,DIMENSION(:),POINTER,INTENT(in) :: heights
   INTEGER,INTENT(in) :: wlev
!
!LOCALS
   INTEGER  :: lev, zon
   AED_REAL :: surf
   LOGICAL  :: w_zones
!
!-------------------------------------------------------------------------------
!BEGIN
   DO zon=1, n_zones
      aedZones(zon)%z_area = 0.  ; aedZones(zon)%z_pc_wet = 0.
   ENDDO
   w_zones = .FALSE.
   surf = heights(wlev)

   zon = 1
   aedZones(1)%z_area(1) = areas(1)
   DO lev=2, wlev
!print*,"heights(",lev,") = ",heights(lev)," zheight(",zon,") = ",aedZones(zon)%z_heights(1)
      IF ( heights(lev) > aedZones(zon)%z_heights(1) ) zon = zon + 1

      aedZones(zon)%z_area(1) = aedZones(zon)%z_area(1) + areas(lev) - areas(lev-1)

      IF ( aedZones(zon)%z_heights(1) > surf ) THEN
         IF (.NOT. w_zones) THEN
            w_zones = .TRUE.
            IF ( zon > 1 ) THEN
               aedZones(zon)%z_pc_wet(1) = surf / (aedZones(zon)%z_heights(1) - aedZones(zon-1)%z_heights(1))
            ELSE
               aedZones(zon)%z_pc_wet(1) = surf / aedZones(zon)%z_heights(1)
            ENDIF
         ENDIF
      ELSE
         aedZones(zon)%z_pc_wet(1) = 1.0
      ENDIF
   ENDDO

   DO zon=1, n_zones
      aedZones(zon)%z_pres(1) = -aedZones(zon)%z_heights(1)
   ENDDO
END SUBROUTINE api_calc_zone_areas
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE api_copy_to_zone(aedZones, n_zones, heights, x_cc, x_cc_hz, x_diag, x_diag_hz, wlev)
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(api_zone_t),DIMENSION(:),INTENT(inout) :: aedZones
   INTEGER,INTENT(in) :: n_zones
!
   AED_REAL,DIMENSION(:),  POINTER,INTENT(in) :: heights
   AED_REAL,DIMENSION(:,:),POINTER,INTENT(in) :: x_cc
   AED_REAL,DIMENSION(:),  POINTER,INTENT(in) :: x_cc_hz
   AED_REAL,DIMENSION(:,:),POINTER,INTENT(in) :: x_diag
   AED_REAL,DIMENSION(:),  POINTER,INTENT(in) :: x_diag_hz
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
! print*,'copy_to_zone', wlev
   DO zon=1,n_zones
      aedZones(zon)%z_cc(:,:) = 0.
      aedZones(zon)%z_cc_hz(:) = 0.
      aedZones(zon)%z_cc_diag(:,:) = 0.
      aedZones(zon)%z_cc_diag_hz(:) = 0.

      aedZones(zon)%z_rad(1) = 0.
      aedZones(zon)%z_salt(1) = 0.
      aedZones(zon)%z_temp(1) = 0.
      aedZones(zon)%z_rho(1) = 0.
      aedZones(zon)%z_extc(1) = 0.
      aedZones(zon)%z_layer_stress(1) = 0.
      aedZones(zon)%z_tss(1) = 0.
      aedZones(zon)%z_par(1) = 0.
      aedZones(zon)%z_nir(1) = 0.
      aedZones(zon)%z_uva(1) = 0.
      aedZones(zon)%z_uvb(1) = 0.
      aedZones(zon)%z_sed_zones(1) = 1.
      aedZones(zon)%z_vel(1) = 0.
   ENDDO

   a_zones = 1
   zcount = 0
   w_zones = .FALSE.
   zon = 1
   DO lev=1,wlev
      IF ( lev > 1 .AND. heights(lev) > aedZones(zon)%z_heights(1) ) THEN
         zon = zon + 1
         IF (zon > n_zones) STOP 'Water level height is higher than highest zone height'
         aedZones(zon)%z_sed_zones(1) = zon
      ENDIF

      ! Pelagic variables need zonifying, in case benthic people want them
      ! Note that this should not be done on benthic vars because it will
      ! introduce errors in z_cc (split layers)
      ! Ideally this average would be based on volume weighting

      aedZones(zon)%z_cc(1,1:nvars) = aedZones(zon)%z_cc(1,1:nvars) + x_cc(lev,1:nvars)
      aedZones(zon)%z_cc_hz(:)      = aedZones(zon)%z_cc_hz(:) + x_cc_hz(:)
      aedZones(zon)%z_cc_diag(1,:)  = aedZones(zon)%z_cc_diag(1,:) + x_diag(lev,:)
      aedZones(zon)%z_cc_diag_hz(:) = aedZones(zon)%z_cc_diag_hz(:) + x_diag_hz(:)

      aedZones(zon)%z_temp(1)         = aedZones(zon)%z_temp(1) + theLake(lev)%Temp
      aedZones(zon)%z_salt(1)         = aedZones(zon)%z_salt(1) + theLake(lev)%Salinity
      aedZones(zon)%z_rho(1)          = aedZones(zon)%z_rho(1)  + theLake(lev)%Density
      aedZones(zon)%z_rad(1)          = aedZones(zon)%z_rad(1)  + theLake(lev)%Light
      aedZones(zon)%z_vel(1)          = aedZones(zon)%z_vel(1)  + theLake(lev)%Umean
      aedZones(zon)%z_extc(1)         = aedZones(zon)%z_extc(1) + theLake(lev)%ExtcCoefSW
      aedZones(zon)%z_layer_stress(1) = aedZones(zon)%z_layer_stress(1) + theLake(lev)%LayerStress

! print*,'to   lev = ',lev, 'zon = ',zon,'x_diag_hz(1) = ',x_diag_hz(1),' z_diag_hz(zon,1) = ',z_diag_hz(zon,1)
      zcount(zon) = zcount(zon) + 1
   ENDDO
   a_zones = zon

   DO zon=1,a_zones
      aedZones(zon)%z_cc(1,1:nvars) = aedZones(zon)%z_cc(1,1:nvars)/zcount(zon)
      aedZones(zon)%z_cc_hz(:)      = aedZones(zon)%z_cc_hz(:)/zcount(zon)
      aedZones(zon)%z_cc_diag(1,:)  = aedZones(zon)%z_cc_diag(1,:)/zcount(zon)
      aedZones(zon)%z_cc_diag_hz(:) = aedZones(zon)%z_cc_diag_hz(:)/zcount(zon)
   ENDDO

   DO zon=1,n_zones
      IF (zcount(zon) /= 0) THEN
         aedZones(zon)%z_temp(1)         = aedZones(zon)%z_temp(1) / zcount(zon)
         aedZones(zon)%z_salt(1)         = aedZones(zon)%z_salt(1) / zcount(zon)
         aedZones(zon)%z_rho(1)          = aedZones(zon)%z_rho(1)  / zcount(zon)
         aedZones(zon)%z_rad(1)          = aedZones(zon)%z_rad(1)  / zcount(zon)
         aedZones(zon)%z_vel(1)          = aedZones(zon)%z_vel(1)  / zcount(zon)
         aedZones(zon)%z_extc(1)         = aedZones(zon)%z_extc(1) / zcount(zon)
         aedZones(zon)%z_layer_stress(1) = aedZones(zon)%z_layer_stress(1) / zcount(zon)
      ELSE
         aedZones(zon)%z_temp(1)         = 0.
         aedZones(zon)%z_salt(1)         = 0.
         aedZones(zon)%z_rho(1)          = 0.
         aedZones(zon)%z_rad(1)          = 0.
         aedZones(zon)%z_vel(1)          = 0.
         aedZones(zon)%z_extc(1)         = 0.
         aedZones(zon)%z_layer_stress(1) = 0.
      ENDIF

      aedZones(zon)%z_dz(1) = 0.
   ENDDO

   surf = heights(wlev)
   IF ( surf > aedZones(1)%z_heights(1) ) THEN
      aedZones(1)%z_depth = aedZones(1)%z_heights(1)
      aedZones(1)%z_dz = aedZones(1)%z_heights(1)
   ELSE
      aedZones(1)%z_depth(1) = surf
      aedZones(1)%z_dz(1) = surf
      w_zones = .TRUE.
   ENDIF

   DO zon=2,a_zones
      aedZones(zon)%z_sed_zones(1) = zon
      IF ( .NOT. w_zones ) THEN
          IF ( surf > aedZones(zon)%z_heights(1) ) THEN
             aedZones(zon)%z_depth(1) = aedZones(zon)%z_heights(1)
             aedZones(zon)%z_dz(1) =  aedZones(zon)%z_heights(1) - aedZones(zon-1)%z_heights(1)
          ELSE
             aedZones(zon)%z_depth(1) = surf
             aedZones(zon)%z_dz(1) = surf - aedZones(zon-1)%z_heights(1)
             w_zones = .TRUE.
          ENDIF
      ELSE
         aedZones(zon)%z_depth(1) = surf
         aedZones(zon)%z_dz(1) = surf - aedZones(zon-1)%z_heights(1)
      ENDIF
   ENDDO
END SUBROUTINE api_copy_to_zone
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE api_copy_from_zone(aedZones, n_zones, heights, x_cc, x_cc_hz, x_diag, x_diag_hz, wlev)
!-------------------------------------------------------------------------------
!ARGUMENTS
   TYPE(api_zone_t),DIMENSION(:),INTENT(in) :: aedZones
   INTEGER,INTENT(in) :: n_zones
!
   AED_REAL,DIMENSION(:),  POINTER,INTENT(in) :: heights
   AED_REAL,DIMENSION(:,:),POINTER,INTENT(inout) :: x_cc
   AED_REAL,DIMENSION(:),  POINTER,INTENT(inout) :: x_cc_hz
   AED_REAL,DIMENSION(:,:),POINTER,INTENT(inout) :: x_diag
   AED_REAL,DIMENSION(:),  POINTER,INTENT(inout) :: x_diag_hz
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
            splitZone = heights(lev-1) < aedZones(zon-1)%z_heights(1)
         ELSE
            splitZone = 0.0 < aedZones(zon-1)%z_heights(1)
         ENDIF
      ELSE
         splitZone = .FALSE.
      ENDIF

      IF (splitZone) THEN
         IF (lev .GT. 1) THEN
            scale = (aedZones(zon-1)%z_heights(1) - heights(lev-1)) / (heights(lev) - heights(lev-1))
         ELSE
            scale = (aedZones(zon-1)%z_heights(1) - 0.0) / (heights(lev) - 0.0)
         ENDIF

         WHERE(aedZones(zon)%z_cc_diag(1,:) /= 0.) &
            x_diag(lev,:) = aedZones(zon)%z_cc_diag(1,:) * scale
         x_cc(lev,v_start:v_end) = aedZones(zon)%z_cc(1,v_start:v_end) * scale

         zon = zon - 1

         WHERE(aedZones(zon)%z_cc_diag(1,:) /= 0.) &
            x_diag(lev,:) = x_diag(lev,:) + (aedZones(zon)%z_cc_diag(1,:) * (1.0 - scale))
         x_cc(lev,v_start:v_end) = x_cc(lev,v_start:v_end) + &
                                   aedZones(zon)%z_cc(1,v_start:v_end) * (1.0 - scale)
      ELSE
         WHERE(aedZones(zon)%z_cc_diag(1,:) /= 0.)
            x_diag(lev,:) = aedZones(zon)%z_cc_diag(1,:)
         ENDWHERE
         x_cc(lev,v_start:v_end) = aedZones(zon)%z_cc(1,v_start:v_end)
      ENDIF
   ENDDO

   ! Set the normal sheet diagnostics to the mean of the zone, weighted by area
   area = 0.
   DO zon=1,n_zones
      area = area + SUM(aedZones(zon)%z_area)
   ENDDO
   DO zon=1,n_zones
      x_diag_hz = x_diag_hz + (aedZones(zon)%z_cc_diag_hz(:) * (aedZones(zon)%z_area(1)/area))
   ENDDO
END SUBROUTINE api_copy_from_zone
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE glm_api_zones
