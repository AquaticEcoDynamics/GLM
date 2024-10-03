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

#define _VOLUME_SCALING_ 0

#undef MISVAL
#ifndef _FORTRAN_SOURCE_
#define _FORTRAN_SOURCE_ 1
#endif

#include "glm.h"

MODULE glm_zones

   USE ISO_C_BINDING

   USE glm_types

   IMPLICIT NONE

   PRIVATE ! By default, make everything private

   AED_REAL,DIMENSION(:,:,:),ALLOCATABLE,TARGET :: z_cc      ! (nsed_zones, n_levs, n_vars)
   AED_REAL,DIMENSION(:,:),  ALLOCATABLE,TARGET :: z_cc_hz   ! (nsed_zones, n_vars)
   AED_REAL,DIMENSION(:,:,:),ALLOCATABLE,TARGET :: z_diag    ! (nsed_zones, n_levs, n_vars)
   AED_REAL,DIMENSION(:,:),  ALLOCATABLE,TARGET :: z_diag_hz ! (nsed_zones, n_vars)

   AED_REAL,DIMENSION(:),POINTER :: zz

   AED_REAL,DIMENSION(:),POINTER :: zone_heights
   INTEGER :: nvars, nbenv, nvdiag, nvdiag_hz

   PUBLIC n_zones, zone_heights, zz
   PUBLIC wq_set_glm_zones, copy_from_zone, copy_to_zone, calc_zone_areas

   PUBLIC z_cc, z_cc_hz, z_diag, z_diag_hz, theZones, theLake
!  PUBLIC n_vars, n_vars_ben, n_vars_diag, n_vars_diag_sheet, zone_var

CONTAINS

!###############################################################################
SUBROUTINE wq_set_glm_zones(numVars, numBenV, numDiagV, numDiagHzV)            &
                                                BIND(C, name="wq_set_glm_zones")
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER,INTENT(in) :: numVars, numBenV, numDiagV, numDiagHzV
!
!LOCALS
   ! none
!
!-------------------------------------------------------------------------------
!BEGIN
   nvars = numVars
   nbenv = numBenV
   nvdiag = numDiagV
   nvdiag_hz = numDiagHzV

   zz => theLake%Height
   zone_heights => theZones%zheight

   ALLOCATE(z_cc(n_zones, MaxLayers, numVars+numBenV))  ; z_cc = 0.
   ALLOCATE(z_cc_hz(n_zones, numVars+numBenV))          ; z_cc_hz = 0.
   ALLOCATE(z_diag(n_zones, MaxLayers, numDiagV))       ; z_diag = 0.
   ALLOCATE(z_diag_hz(n_zones+1, numDiagHzV))           ; z_diag_hz = 0.
   theZones%zarea = 0.
END SUBROUTINE wq_set_glm_zones
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE calc_zone_areas(areas, wlev, surf)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,DIMENSION(:),INTENT(in) :: areas
   INTEGER,INTENT(in) :: wlev
   AED_REAL :: surf
!
!LOCALS
   INTEGER  :: lev, zon
   LOGICAL  :: w_zones
#if _VOLUME_SCALING_
   INTEGER  :: ij
   AED_REAL :: x, y
#endif
!
!-------------------------------------------------------------------------------
!BEGIN
   theZones%zarea = 0.  ; theZones%z_pc_wet = 0.
   w_zones = .FALSE.

#if _VOLUME_SCALING_

   DO zon=1, n_zones
      x = zone_heights(zon) * 10.0
      y = AMOD(x, 1.0)
      ij = INT(x - y)
      IF (ij .GT. NMorph) THEN
         y = y + FLOAT(ij - NMorph)
         ij = NMorph
      ENDIF

      zvol(zon) = MphLevelVol(ij) + y * dMphLevelVol(ij)
      theZones(zon)%zarea = MphLevelArea(ij) + y * dMphLevelArea(ij)
   ENDDO

   zon = 1
   DO lev=2, wlev
      IF ( zz(lev) > zone_heights(zon) ) zon = zon + 1

      IF ( zone_heights(zon) > surf ) THEN
         IF (.NOT. w_zones ) THEN
            w_zones = .TRUE.
            IF ( zon > 1 ) THEN
               theZones(zon)%z_pc_wet = surf / (zone_heights(zon) - zone_heights(zon-1))
            ELSE
               theZones(zon)%z_pc_wet = surf / zone_heights(zon)
            ENDIF
         ENDIF
      ELSE
         theZones(zon)%z_pc_wet = 1.0
      ENDIF
   ENDDO

#else

   zon = 1
   theZones(1)%zarea = areas(1)
   DO lev=2, wlev
      IF ( zz(lev) > zone_heights(zon) ) zon = zon + 1

      theZones(zon)%zarea = theZones(zon)%zarea + areas(lev) - areas(lev-1)

      IF ( zone_heights(zon) > surf ) THEN
         IF (.NOT. w_zones) THEN
            w_zones = .TRUE.
            IF ( zon > 1 ) THEN
               theZones(zon)%z_pc_wet = surf / (zone_heights(zon) - zone_heights(zon-1))
            ELSE
               theZones(zon)%z_pc_wet = surf / zone_heights(zon)
            ENDIF
         ENDIF
      ELSE
         TheZones(zon)%z_pc_wet = 1.0
      ENDIF
   ENDDO

#endif

   theZones(1:n_zones)%zpres = -zone_heights(1:n_zones)
END SUBROUTINE calc_zone_areas
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE copy_to_zone(x_cc, x_diag, x_diag_hz, wlev)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,DIMENSION(:,:),INTENT(in) :: x_cc
   AED_REAL,DIMENSION(:,:),INTENT(in) :: x_diag
   AED_REAL,DIMENSION(:),INTENT(in) :: x_diag_hz
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
   z_cc(:,:,1:nvars) = 0.
   z_diag(:,:,:) = 0.
   z_diag_hz(:,:) = 0.
   theZones%zrad = 0. ; theZones%zsalt = 0. ; theZones%ztemp = 0. ; theZones%zrho = 0.
   theZones%zextc = 0. ; theZones%zlayer_stress = 0. ; theZones%ztss = 0.
   theZones%zpar = 0. ; theZones%znir = 0. ; theZones%zuva = 0. ; theZones%zuvb = 0.
   theZones%z_sed_zones = 1. ; theZones%zvel = 0.

   a_zones = 1
   zcount = 0
   w_zones = .FALSE.
   zon = 1
   DO lev=1,wlev
      IF ( lev > 1 .AND. zz(lev) > zone_heights(zon) ) THEN
         zon = zon + 1
         IF (zon > n_zones) STOP 'Water level height is higher than highest zone height'
         theZones(zon)%z_sed_zones = zon
      ENDIF

      ! Pelagic variables need zonifying, in case benthic people want them
      ! Note that this should not be done on benthic vars because it will
      ! introduce errors in z_cc (split layers)
      ! Ideally this average would be based on volume weighting

      z_cc(zon,lev,1:nvars) = z_cc(zon,lev,1:nvars) + x_cc(lev,1:nvars)
      z_diag(zon,lev,:)     = z_diag(zon,lev,:) + x_diag(lev,:)
      z_diag_hz(zon,:)  = z_diag_hz(zon,:) + x_diag_hz(:)

      theZones(zon)%ztemp         = theZones(zon)%ztemp + theLake(lev)%Temp
      theZones(zon)%zsalt         = theZones(zon)%zsalt + theLake(lev)%Salinity
      theZones(zon)%zrho          = theZones(zon)%zrho  + theLake(lev)%Density
      theZones(zon)%zrad          = theZones(zon)%zrad  + theLake(lev)%Light
      theZones(zon)%zvel          = theZones(zon)%zvel  + theLake(lev)%Umean
      theZones(zon)%zextc         = theZones(zon)%zextc + theLake(lev)%ExtcCoefSW
      theZones(zon)%zlayer_stress = theZones(zon)%zlayer_stress + theLake(lev)%LayerStress

      zcount(zon) = zcount(zon) + 1
   ENDDO
   a_zones = zon

   DO zon=1,a_zones
      z_cc(zon,:,1:nvars) = z_cc(zon,:,1:nvars)/zcount(zon)
      z_diag(zon,:,:)     = z_diag(zon,:,:)/zcount(zon)
      z_diag_hz(zon,:)  = z_diag_hz(zon,:)/zcount(zon)
   ENDDO

   WHERE (zcount /= 0)
      theZones%ztemp         = theZones%ztemp / zcount
      theZones%zsalt         = theZones%zsalt / zcount
      theZones%zrho          = theZones%zrho  / zcount
      theZones%zrad          = theZones%zrad  / zcount
      theZones%zvel          = theZones%zvel  / zcount
      theZones%zextc         = theZones%zextc / zcount
      theZones%zlayer_stress = theZones%zlayer_stress / zcount
   ELSEWHERE
      theZones%ztemp         = 0.
      theZones%zsalt         = 0.
      theZones%zrho          = 0.
      theZones%zrad          = 0.
      theZones%zvel          = 0.
      theZones%zextc         = 0.
      theZones%zlayer_stress = 0.
   ENDWHERE

   theZones%zdz = 0.
   surf = zz(wlev)
   IF ( surf > zone_heights(1) ) THEN
      theZones(1)%zdepth = zone_heights(1)
      theZones(1)%zdz = zone_heights(1)
   ELSE
      theZones%zdepth = surf
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
SUBROUTINE copy_from_zone(x_cc, x_diag, x_diag_hz, wlev)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,DIMENSION(:,:),INTENT(inout) :: x_cc
   AED_REAL,DIMENSION(:,:),INTENT(inout) :: x_diag
   AED_REAL,DIMENSION(:),INTENT(inout) :: x_diag_hz
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
            splitZone = zz(lev-1) < zone_heights(zon-1)
         ELSE
            splitZone = 0.0 < zone_heights(zon-1)
         ENDIF
      ELSE
         splitZone = .FALSE.
      ENDIF

      IF (splitZone) THEN
         IF (lev .GT. 1) THEN
            scale = (zone_heights(zon-1) - zz(lev-1)) / (zz(lev) - zz(lev-1))
         ELSE
            scale = (zone_heights(zon-1) - 0.0) / (zz(lev) - 0.0)
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
   area = SUM(theZones(1:n_zones)%zarea)
   DO zon=1,n_zones
      x_diag_hz = x_diag_hz + (z_diag_hz(zon,:) * (theZones(zon)%zarea/area))
   ENDDO

END SUBROUTINE copy_from_zone
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE GInitialTemp(m,depth,wv,topTemp,botTemp,nSPinUpDays,tNew) BIND(C, name="InitialTemp")
!-------------------------------------------------------------------------------
#ifdef AED2
   USE aed2_util, ONLY : InitialTemp
#else
   USE aed_util, ONLY : InitialTemp
#endif
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
SUBROUTINE ZSoilTemp(izone) BIND(C, name="zZSoilTemp")
!-------------------------------------------------------------------------------
#ifdef AED2
   USE aed2_util, ONLY : SoilTemp
#else
   USE aed_util, ONLY : SoilTemp
#endif
!ARGUMENTS
   TYPE(C_PTR),INTENT(inout) :: izone
!LOCALS
   TYPE(ZoneType),POINTER :: zone
   TYPE(SedLayerType),DIMENSION(:),POINTER :: layers
!-------------------------------------------------------------------------------
!BEGIN
   CALL C_F_POINTER(izone, zone);
   CALL C_F_POINTER(zone%c_layers, layers, [zone%n_sed_layers]);
   CALL SoilTemp(zone%n_sed_layers, layers%depth, layers%vwc, zone%ztemp, layers%temp, zone%heatflux)
END SUBROUTINE ZSoilTemp
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE GSoilTemp(m,depth,wv,topTemp,temp,heatflux) BIND(C, name="SoilTemp")
!-------------------------------------------------------------------------------
#ifdef AED2
   USE aed2_util, ONLY : SoilTemp
#else
   USE aed_util, ONLY : SoilTemp
#endif
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

END MODULE glm_zones
