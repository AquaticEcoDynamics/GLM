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

   AED_REAL,ALLOCATABLE,DIMENSION(:,:),TARGET :: z_cc  ! (nsed_zones, nsed_vars)
   AED_REAL,ALLOCATABLE,DIMENSION(:,:),TARGET :: z_diag    ! (nsed_zones, n_vars)
   AED_REAL,ALLOCATABLE,DIMENSION(:,:),TARGET :: z_diag_hz ! (nsed_zones, n_vars)

   AED_REAL,DIMENSION(:),POINTER :: zz
   AED_REAL,ALLOCATABLE,DIMENSION(:),TARGET :: z_sed_zones

   AED_REAL,ALLOCATABLE,DIMENSION(:) :: z_pc_wet
   AED_REAL,ALLOCATABLE,DIMENSION(:),TARGET :: zrad, zsalt, ztemp, zrho, zarea
   AED_REAL,ALLOCATABLE,DIMENSION(:),TARGET :: zextc_coef, zlayer_stress, ztss

   AED_REAL,ALLOCATABLE,DIMENSION(:),TARGET :: zdz, zpar, zdepth, zpres
   AED_REAL,ALLOCATABLE,DIMENSION(:),TARGET :: znir, zuva, zuvb

   INTEGER :: n_zones, w_zones
   AED_REAL,DIMENSION(:),POINTER :: zone_heights
   INTEGER :: nvars, nbenv


   TYPE(LakeDataType),DIMENSION(:),POINTER :: theLake

   PUBLIC n_zones, zone_heights, zz, z_cc, theLake
   PUBLIC wq_set_glm_zones, copy_from_zone, copy_to_zone, calc_zone_areas

   PUBLIC zrad, zsalt, ztemp, zrho, zarea, zextc_coef, zlayer_stress, ztss, zdz, zpar
   PUBLIC zdepth, zpres, z_pc_wet, z_sed_zones, z_diag, z_diag_hz
   PUBLIC znir, zuva, zuvb

CONTAINS

!###############################################################################
SUBROUTINE wq_set_glm_zones(z_heights, numZones, numVars, numBenV) BIND(C, name="wq_set_glm_zones")
!-------------------------------------------------------------------------------
!ARGUMENTS
   CINTEGER,INTENT(in) :: numZones, numVars, numBenV
   AED_REAL,TARGET,INTENT(in) :: z_heights(1:numZones)
!
!LOCALS
!  INTEGER :: i
!  AED_REAL :: surf
!
!-------------------------------------------------------------------------------
!BEGIN
   n_zones = numZones
   zone_heights => z_heights
   nvars = numVars
   nbenv = numBenV
   ALLOCATE(z_cc(n_zones, numVars+numBenV))
   z_cc = 900.!   !MH if i intiialise this in iinit then nothing happens so doing it here.
   ALLOCATE(zrad(n_zones))
   ALLOCATE(zsalt(n_zones))
   ALLOCATE(ztemp(n_zones))
   ALLOCATE(zrho(n_zones))
   ALLOCATE(zarea(n_zones))
   ALLOCATE(zextc_coef(n_zones))
   ALLOCATE(zlayer_stress(n_zones))
   ALLOCATE(ztss(n_zones))
   ALLOCATE(zdz(n_zones))
   ALLOCATE(zpar(n_zones))
   ALLOCATE(znir(n_zones))
   ALLOCATE(zuva(n_zones))
   ALLOCATE(zuvb(n_zones))
   ALLOCATE(zpres(n_zones))
   ALLOCATE(zdepth(n_zones))
   ALLOCATE(z_sed_zones(n_zones))
   ALLOCATE(z_pc_wet(n_zones))
!  zdz(1) = z_dep(1)
!  DO i=2,n_zones
!     zdz(i) = z_dep(i) - z_dep(i-1)
!  ENDDO
!  DO i=1,n_zones ; z_sed_zones(i) = i; ENDDO
   zarea = 0.
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
   AED_REAL :: x, y
#endif
!
!-------------------------------------------------------------------------------
!BEGIN
   zarea = 0.  ; z_pc_wet = 0. ; w_zones = 0

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
      zarea(zon) = MphLevelArea(ij) + y * dMphLevelArea(ij)
   ENDDO

   zon = 1
   DO lev=2, wlev
      IF ( zz(lev) > zone_heights(zon) ) zon = zon + 1

      IF ( zone_heights(zon) > surf ) THEN
         IF (w_zones == 0) THEN
            w_zones = lev
            z_pc_wet(zon) = surf / (zone_heights(zon) - zone_heights(zon-1))
         ENDIF
      ELSE
         z_pc_wet(zon) = 1.0
      ENDIF
   ENDDO

#else

   zon = 1
   zarea(1) = areas(1)
   DO lev=2, wlev
      IF ( zz(lev) > zone_heights(zon) ) zon = zon + 1

      zarea(zon) = zarea(zon) + areas(lev) - areas(lev-1)

      IF ( zone_heights(zon) > surf ) THEN
         IF (w_zones == 0) THEN
            w_zones = lev
            z_pc_wet(zon) = surf / (zone_heights(zon) - zone_heights(zon-1))
         ENDIF
      ELSE
         z_pc_wet(zon) = 1.0
      ENDIF
   ENDDO

#endif

   zpres(1:n_zones) = -zone_heights(1:n_zones)
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
         splitZone = zz(lev-1) < zone_heights(zon-1)
      ELSE
         splitZone = .FALSE.
      ENDIF

      IF (splitZone) THEN
         scale = (zone_heights(zon-1) - zz(lev-1)) / (zz(lev) - zz(lev-1))
         x_diag(lev,:) = z_diag(zon,:) * scale
         x_cc(lev,v_start:v_end) = z_cc(zon,v_start:v_end) * scale

         zon = zon - 1

         x_diag(lev,:) = x_diag(lev,:) + (z_diag(zon,:) * (1.0 - scale))
         x_cc(lev,v_start:v_end) = x_cc(lev,v_start:v_end) + &
                                   z_cc(zon,v_start:v_end) * (1.0 - scale)
      ELSE
         x_diag(lev,:) = z_diag(zon,:)
         x_cc(lev,v_start:v_end) = z_cc(zon,v_start:v_end)
      ENDIF
   ENDDO
   DO lev=1,n_zones
      x_diag_hz(v_start:v_end) = x_diag_hz(v_start:v_end) + z_diag_hz(lev,1:nbenv)
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
   zrad = 0. ; zsalt = 0. ; ztemp = 0. ; zrho = 0.
   zextc_coef = 0. ; zlayer_stress = 0. ; ztss = 0. ; zpar = 0.
   znir = 0. ; zuva = 0. ; zuvb = 0. ; z_sed_zones = 1.

   zcount = 0
   DO lev=1,wlev
      IF ( zz(lev) > zone_heights(zon) ) THEN
         zon = zon + 1
         z_sed_zones(zon) = zon
      ENDIF

      ! Pelagic variables need zonifying, in case benthic people want them
      ! Note that this should not be done on benthic vars becasue it will
      ! introduce errors in z_cc (split layers)
      ! Ideally this average would be based on volume weighting

      z_cc(zon,1:nvars) = z_cc(zon,1:nvars) + x_cc(lev,1:nvars)

      ztemp(zon)         = ztemp(zon) + theLake(lev)%Temp
      zsalt(zon)         = zsalt(zon) + theLake(lev)%Salinity
      zrho(zon)          = zrho(zon)  + theLake(lev)%Density
      zrad(zon)          = zrad(zon)  + theLake(lev)%Light
      zextc_coef(zon)    = zextc_coef(zon) + theLake(lev)%ExtcCoefSW
      zlayer_stress(zon) = zlayer_stress(zon) + theLake(lev)%LayerStress

      zcount(zon) = zcount(zon) + 1
   ENDDO

   DO zon=1,n_zones
      z_cc(zon,1:nvars) = z_cc(zon,1:nvars)/zcount(zon)
   ENDDO

   ztemp         = ztemp / zcount
   zsalt         = zsalt / zcount
   zrho          = zrho  / zcount
   zrad          = zrad  / zcount
   zextc_coef    = zextc_coef / zcount
   zlayer_stress = zlayer_stress / zcount

   surf = zz(wlev)
   IF ( surf > zone_heights(1) ) THEN
      zdepth(1) = zone_heights(1)
   ELSE
      zdepth = surf
      w_zones = 1
   ENDIF

   DO zon=2,n_zones
      z_sed_zones(zon) = zon
      IF ( w_zones == 0 ) THEN
          IF ( surf > zone_heights(zon) ) THEN
             zdepth(zon) = zone_heights(zon)
          ELSE
             zdepth(zon) = surf
             w_zones = zon
          ENDIF
      ELSE
         zdepth(zon) = surf
      ENDIF
   ENDDO
END SUBROUTINE copy_to_zone
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


END MODULE glm_zones
