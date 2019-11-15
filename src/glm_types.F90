!###############################################################################
!#                                                                             #
!#  glm_types.F90                                                              #
!#                                                                             #
!#  A module to define constants and types.                                    #
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
#include "glm.h"


!************************* Important Note **************************************
!* The order of entries in these structures MUST match those in glm_types.h    *
!************************* Important Note **************************************

!*******************************************************************************
MODULE glm_types

   USE ISO_C_BINDING

   IMPLICIT NONE

!===============================================================================
!GLOBAL CONSTANTS

   AED_REAL,PARAMETER :: missing = MISVAL

!===============================================================================
!TYPE DECLARATIONS

   !#===========================================================#!
   TYPE,BIND(C) :: StringT
      CINTEGER :: Len
      CCHARACTER :: S(40)
   END TYPE StringT

   !#===========================================================#!
   !# Structured type for key global lake environmental vars
   !# A Lake will be an allocated array of MaxLayers of these
   TYPE,BIND(C) :: LakeDataType
      AED_REAL :: Density          !# density kg/m3
      AED_REAL :: Temp             !# temperature
      AED_REAL :: Salinity         !# salinity
      AED_REAL :: Height           !# layer heights above the bottom
      AED_REAL :: MeanHeight       !# mean height of a layer
      AED_REAL :: LayerVol         !# volume of layer
      AED_REAL :: LayerArea        !# area of layer

      AED_REAL :: Light            !# PAR, photosynthetically active radiation
      AED_REAL :: ExtcCoefSW       !# Kd, light extinction coefficient

      AED_REAL :: Vol1             !# Cumulative volume to this layer top
      AED_REAL :: Epsilon          !# Diffusivity

      AED_REAL :: Umean            !# Mean velocity
      AED_REAL :: Uorb             !# Orbital velocity
      AED_REAL :: LayerStress      !# Layer Stress
   END TYPE LakeDataType

   !#===========================================================#!
   !# Structured type for Met vars
   TYPE,BIND(C) :: MetDataType
      AED_REAL :: Rain             !# rainfall
      AED_REAL :: RelHum           !# relative humidty
      AED_REAL :: SatVapDef        !# vapour pressure
      AED_REAL :: LongWave         !# longwave radiation
      AED_REAL :: ShortWave        !# shortwave radiation
      AED_REAL :: AirTemp          !# temperature
      AED_REAL :: WindSpeed        !# windspeed
      AED_REAL :: Snow             !# snowfall
      AED_REAL :: RainConcPO4      !# Concentration of PO4 in rain
      AED_REAL :: RainConcTP       !# Concentration of TP in rain
      AED_REAL :: RainConcNO3      !# Concentration of NO3 in rain
      AED_REAL :: RainConcNH4      !# Concentration of NH4 in rain
      AED_REAL :: RainConcTN       !# Concentration of TN in rain
      AED_REAL :: RainConcSi       !# Concentration of SI in rain
      AED_REAL :: WindDir          !# Wind direction
      AED_REAL :: As               !# Area of sheltering
   END TYPE MetDataType

   !#===========================================================#!
   !# Structured type for Surface Data vars
   TYPE,BIND(C) :: SurfaceDataType
      AED_REAL :: Evap             !# Evaporation
      AED_REAL :: delzBlueIce      !# Thickness of blue ice layer
      AED_REAL :: delzWhiteIce     !# Thickness of white ice layer
      AED_REAL :: delzSnow         !# Thickness of snow layer
      AED_REAL :: dHt              !# Change in thickness of snow / ice layer
      AED_REAL :: RhoSnow          !# Density of snow layer (kg/m^3)
      AED_REAL :: dailyEvap        !# Daily Evaporation (m3/day)
      AED_REAL :: dailyRain        !# Daily Rain (m3/day)
      AED_REAL :: dailyRunoff      !# Daily Rain (m3/day)
      AED_REAL :: dailySnow        !# Daily Snow (m3/day)
      AED_REAL :: dailyQsw         !# Daily Short Wave Radiation (J/day)
      AED_REAL :: dailyQe          !# Daily Latent Heat(J/day)
      AED_REAL :: dailyQh          !# Daily Sensible Heat (J/day)
      AED_REAL :: dailyQlw         !# Daily Long Wave Radiation (J/day)
      AED_REAL :: dailyInflow      !# Total Daily Inflow (m3/day)
      AED_REAL :: dailyOutflow     !# Total Daily Outflow (m3/day)
      AED_REAL :: dailyOverflow    !# Total Daily Overflow (m3/day)
      AED_REAL :: albedo           !# Daily surface albedo
      AED_REAL :: dailyzonL        !# Average z/L value, atmos stability
   END TYPE SurfaceDataType

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !#===========================================================#!
   !# Structured type for Sediment Layer
   TYPE,BIND(C) :: SedLayerType
      AED_REAL :: depth            !# Layer depth
      AED_REAL :: temp             !# Layer temperature
      AED_REAL :: vwc
      AED_REAL :: wq
   END TYPE SedLayerType

   !#===========================================================#!
   !# Structured type for iSediment Zones
   TYPE,BIND(C) :: ZoneType
      AED_REAL :: zheight
      AED_REAL :: zrad
      AED_REAL :: zsalt
      AED_REAL :: ztemp
      AED_REAL :: zrho
      AED_REAL :: zarea
      AED_REAL :: zextc_coef
      AED_REAL :: zlayer_stress
      AED_REAL :: ztss
      AED_REAL :: zdz
      AED_REAL :: zpar
      AED_REAL :: znir
      AED_REAL :: zuva
      AED_REAL :: zuvb
      AED_REAL :: zpres
      AED_REAL :: zdepth
      AED_REAL :: z_sed_zones
      AED_REAL :: z_pc_wet
      AED_REAL :: heatflux
      CINTEGER :: n_sedLayers;     !# number of sediment layers
      TYPE(C_PTR) :: c_layers      !# array of sed layers
   END TYPE ZoneType

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

CONTAINS

!#
!# These are 2 useful routines for converting between fortran and C strings
!# They are in here because, well, I guess they are sort of type conversions
!#

!###############################################################################
SUBROUTINE make_string(s1,s2,len)
   CHARACTER(len=*),INTENT(out) :: s1
   CHARACTER,INTENT(in) :: s2(*)
   CSIZET,INTENT(in)    :: len
!LOCALS
   CSIZET :: i
!
!-------------------------------------------------------------------------------
!BEGIN
   s1 = ''
   DO i=1,len
      s1 = s1 // " "
      s1(i:i) = s2(i)
   ENDDO
END SUBROUTINE make_string
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


!###############################################################################
FUNCTION make_c_string(s1,s2) RESULT(len)
   CCHARACTER,INTENT(out) :: s1(*)
   CHARACTER(len=*),INTENT(in) :: s2
!LOCALS
   INTEGER :: i
   INTEGER :: len
!
!-------------------------------------------------------------------------------
!BEGIN
   len = len_trim(s2)
   DO i=1,len
      s1(i) = s2(i:i)
   ENDDO
   s1(len+1) = ACHAR(0)
END FUNCTION make_c_string
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE glm_types
