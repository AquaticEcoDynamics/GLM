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
!# Copyright 2013 - 2018 -  The University of Western Australia                #
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

   AED_REAL,ALLOCATABLE,DIMENSION(:,:),TARGET :: z_cc   !(nsed_zones, nsed_vars)
   AED_REAL,ALLOCATABLE,DIMENSION(:,:),TARGET :: z_diag    !(nsed_zones, n_vars)
   AED_REAL,ALLOCATABLE,DIMENSION(:,:),TARGET :: z_diag_hz !(nsed_zones, n_vars)

   AED_REAL,DIMENSION(:),POINTER :: zz

   INTEGER :: n_zones, w_zones

!  TYPE(ZoneType),ALLOCATABLE,DIMENSION(:),TARGET :: theZones
   TYPE(ZoneType),DIMENSION(:),POINTER :: theZones

   AED_REAL,DIMENSION(:),POINTER :: zone_heights
   INTEGER :: nvars, nbenv

   TYPE(LakeDataType),DIMENSION(:),POINTER :: theLake

   PUBLIC n_zones, zone_heights, zz, z_cc, theLake
   PUBLIC wq_set_glm_zones, copy_from_zone, copy_to_zone, calc_zone_areas

   PUBLIC z_diag, z_diag_hz, theZones

CONTAINS

!###############################################################################
SUBROUTINE wq_set_glm_zones(Zones, numZones, numVars, numBenV) BIND(C, name="wq_set_glm_zones")
!-------------------------------------------------------------------------------
!ARGUMENTS
!  AED_REAL,TARGET,INTENT(in) :: z_heights(1:numZones)
   TYPE(C_PTR),VALUE :: Zones
   CINTEGER,INTENT(in) :: numZones, numVars, numBenV
!
!LOCALS
   INTEGER :: i
   INTEGER :: n_sed_layers = 20
!  AED_REAL :: surf
   TYPE(SedLayerType),DIMENSION(:),POINTER :: layers
!
!-------------------------------------------------------------------------------
!BEGIN
   n_zones = numZones
   nvars = numVars
   nbenv = numBenV
   CALL C_F_POINTER(Zones, theZones, [numZones]);

   zone_heights => theZones%zheight

!  print *,'C_F_POINTER',zone_heights(:),theZones(1)%zheight,theZones(2)%zheight

!  CALL C_F_POINTER(z_heights, zone_heights, [numZones]);
   DO i=1,n_zones
!     ALLOCATE(theZones(i)%layers(n_sed_layers))
      CALL C_F_POINTER(theZones(i)%c_layers, layers, [n_sed_layers]);
   ENDDO
   ALLOCATE(z_cc(n_zones, numVars+numBenV))
   z_cc = 900.!   !MH if i initialise this in init then nothing happens so doing it here.

!  zdz(1) = z_dep(1)
!  DO i=2,n_zones
!     zdz(i) = z_dep(i) - z_dep(i-1)
!  ENDDO
!  DO i=1,n_zones ; z_sed_zones(i) = i; ENDDO
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
#if _VOLUME_SCALING_
   INTEGER :: ij
   AED_REAL :: x, y
#endif
!
!-------------------------------------------------------------------------------
!BEGIN
   theZones%zarea = 0.  ; theZones%z_pc_wet = 0. ; w_zones = 0

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
         IF (w_zones == 0) THEN
            w_zones = lev
            theZones(zon)%z_pc_wet = surf / (zone_heights(zon) - zone_heights(zon-1))
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
         IF (w_zones == 0) THEN
            w_zones = lev
            theZones(zon)%z_pc_wet = surf / (zone_heights(zon) - zone_heights(zon-1))
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
   AED_REAL :: scale
   LOGICAL  :: splitZone
!
!-------------------------------------------------------------------------------
!BEGIN
   v_start = nvars+1 ; v_end = nvars+nbenv

   zon = n_zones
   DO lev=wlev,1,-1
      IF ( zon .NE. 1 ) THEN
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

         WHERE(z_diag(zon,:) /= 0.) &
            x_diag(lev,:) = z_diag(zon,:) * scale
         x_cc(lev,v_start:v_end) = z_cc(zon,v_start:v_end) * scale

         zon = zon - 1

         WHERE(z_diag(zon,:) /= 0.) &
            x_diag(lev,:) = x_diag(lev,:) + (z_diag(zon,:) * (1.0 - scale))
         x_cc(lev,v_start:v_end) = x_cc(lev,v_start:v_end) + &
                                   z_cc(zon,v_start:v_end) * (1.0 - scale)
      ELSE
         WHERE(z_diag(zon,:) /= 0.) &
            x_diag(lev,:) = z_diag(zon,:)
         x_cc(lev,v_start:v_end) = z_cc(zon,v_start:v_end)
      ENDIF
   ENDDO
   DO lev=1,n_zones
!     x_diag_hz(v_start:v_end) = x_diag_hz(v_start:v_end) + z_diag_hz(lev,1:nbenv)
      x_diag_hz = x_diag_hz + z_diag_hz(lev,:)
   ENDDO
END SUBROUTINE copy_from_zone
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE copy_to_zone(x_cc, wlev)
!-------------------------------------------------------------------------------
!ARGUMENTS
   AED_REAL,DIMENSION(:,:),INTENT(in) :: x_cc
   INTEGER,INTENT(in) :: wlev
!
!LOCALS
   INTEGER  :: zon, lev
   AED_REAL :: surf
   INTEGER  :: zcount(n_zones)
!
!-------------------------------------------------------------------------------
!BEGIN
   zon = 1 ; z_cc(:,1:nvars) = 0.
   theZones%zrad = 0. ; theZones%zsalt = 0. ; theZones%ztemp = 0. ; theZones%zrho = 0.
   theZones%zextc_coef = 0. ; theZones%zlayer_stress = 0. ; theZones%ztss = 0. ; theZones%zpar = 0.
   theZones%znir = 0. ; theZones%zuva = 0. ; theZones%zuvb = 0. ; theZones%z_sed_zones = 1.

   zcount = 0
   DO lev=1,wlev
      IF ( zz(lev) > zone_heights(zon) ) THEN
         zon = zon + 1
         theZones(zon)%z_sed_zones = zon
      ENDIF

      ! Pelagic variables need zonifying, in case benthic people want them
      ! Note that this should not be done on benthic vars because it will
      ! introduce errors in z_cc (split layers)
      ! Ideally this average would be based on volume weighting

      z_cc(zon,1:nvars) = z_cc(zon,1:nvars) + x_cc(lev,1:nvars)

      theZones(zon)%ztemp         = theZones(zon)%ztemp + theLake(lev)%Temp
      theZones(zon)%zsalt         = theZones(zon)%zsalt + theLake(lev)%Salinity
      theZones(zon)%zrho          = theZones(zon)%zrho  + theLake(lev)%Density
      theZones(zon)%zrad          = theZones(zon)%zrad  + theLake(lev)%Light
      theZones(zon)%zextc_coef    = theZones(zon)%zextc_coef + theLake(lev)%ExtcCoefSW
      theZones(zon)%zlayer_stress = theZones(zon)%zlayer_stress + theLake(lev)%LayerStress

      zcount(zon) = zcount(zon) + 1
   ENDDO

   DO zon=1,n_zones
      z_cc(zon,1:nvars) = z_cc(zon,1:nvars)/zcount(zon)
   ENDDO

   theZones%ztemp         = theZones%ztemp / zcount
   theZones%zsalt         = theZones%zsalt / zcount
   theZones%zrho          = theZones%zrho  / zcount
   theZones%zrad          = theZones%zrad  / zcount
   theZones%zextc_coef    = theZones%zextc_coef / zcount
   theZones%zlayer_stress = theZones%zlayer_stress / zcount

   surf = zz(wlev)
   IF ( surf > zone_heights(1) ) THEN
      theZones(1)%zdepth = zone_heights(1)
   ELSE
      theZones%zdepth = surf
      w_zones = 1
   ENDIF

   DO zon=2,n_zones
      theZones(zon)%z_sed_zones = zon
      IF ( w_zones == 0 ) THEN
          IF ( surf > zone_heights(zon) ) THEN
             theZones(zon)%zdepth = zone_heights(zon)
          ELSE
             theZones(zon)%zdepth = surf
             w_zones = zon
          ENDIF
      ELSE
         theZones(zon)%zdepth = surf
      ENDIF
   ENDDO
END SUBROUTINE copy_to_zone
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE GInitialTemp(m,depth,wv,topTemp,botTemp,nSPinUpDays,tNew) BIND(C, name="InitialTemp")
!-------------------------------------------------------------------------------
   USE aed2_util
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
SUBROUTINE ZSoilTemp(izone) BIND(C, name="ZSoilTemp")
!-------------------------------------------------------------------------------
   USE aed2_util
!ARGUMENTS
   TYPE(C_PTR),INTENT(inout) :: izone
!LOCALS
   TYPE(ZoneType),POINTER :: zone
   TYPE(SedLayerType),DIMENSION(:),POINTER :: layers
!-------------------------------------------------------------------------------
!BEGIN
   CALL C_F_POINTER(izone, zone);
   CALL C_F_POINTER(zone%c_layers, layers, [zone%n_sedLayers]);
   CALL SoilTemp(zone%n_sedLayers, layers%depth, layers%vwc, zone%ztemp, layers%temp, zone%heatflux)
END SUBROUTINE ZSoilTemp
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
SUBROUTINE GSoilTemp(m,depth,wv,topTemp,temp,heatflux) BIND(C, name="SoilTemp")
!-------------------------------------------------------------------------------
   USE aed2_util
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
