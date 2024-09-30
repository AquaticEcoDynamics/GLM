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
      AED_REAL :: AirPres          !# air pressure
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
      AED_REAL :: dailyzonL        !# Average z/L value, daily atmospheric stability
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
      AED_REAL :: zextc
      AED_REAL :: zlayer_stress
      AED_REAL :: ztss
      AED_REAL :: zdz
      AED_REAL :: zvel
      AED_REAL :: zpar
      AED_REAL :: znir
      AED_REAL :: zuva
      AED_REAL :: zuvb
      AED_REAL :: zpres
      AED_REAL :: zdepth
      AED_REAL :: z_sed_zones
      AED_REAL :: z_pc_wet
      AED_REAL :: heatflux
      CINTEGER :: n_sed_layers;    !# number of sediment layers
      TYPE(C_PTR) :: c_layers      !# array of sed layers
   END TYPE ZoneType

!  !#===========================================================#!
!  !# Structured type for Particle Transport Model (PTM)
!  TYPE,BIND(C) :: ParticleDataType
!      INTEGER  :: Status         ! indivdual particle status
!      AED_REAL :: Height
!      AED_REAL :: Mass
!      AED_REAL :: Diam
!      AED_REAL :: Density
!      AED_REAL :: Velocity
!      AED_REAL :: vvel
!      CINTEGER :: Layer
!  END TYPE ParticleDataType

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

   !================================================================
   !# variables in C code of GLM
   !----------------------------------------------------------------
   TYPE(CINTEGER),BIND(C, name="MaxLayers") :: MaxLayers
   TYPE(C_PTR),BIND(C, name="Lake")         :: cLake
   TYPE(CINTEGER),BIND(C, name="n_zones")   :: n_zones
   TYPE(C_PTR),BIND(C, name="theZones")     :: cZones

   TYPE(C_PTR),BIND(C, name="pMetData")     :: cMetData
   TYPE(C_PTR),BIND(C, name="pSurfData")    :: cSurfData

   TYPE(LakeDataType),   DIMENSION(:),POINTER :: theLake
   TYPE(MetDataType),                 POINTER :: MetData   !# Meteorological data
   TYPE(MetDataType),    DIMENSION(:),POINTER :: aMetData  !# Meteorological data
   TYPE(SurfaceDataType),             POINTER :: SurfData  !# Surface Data
   TYPE(SurfaceDataType),DIMENSION(:),POINTER :: aSurfData !# Surface Data
   TYPE(ZoneType),       DIMENSION(:),POINTER :: theZones

   TYPE(CLOGICAL),BIND(C, name="mobility_off")     :: mobility_off
   TYPE(CLOGICAL),BIND(C, name="bioshade_feedback"):: bioshade_feedback
   TYPE(CLOGICAL),BIND(C, name="repair_state")     :: repair_state
   TYPE(CLOGICAL),BIND(C, name="do_plots")         :: do_plots
   TYPE(CLOGICAL),BIND(C, name="link_rain_loss")   :: link_rain_loss
   TYPE(CLOGICAL),BIND(C, name="link_solar_shade") :: link_solar_shade
   TYPE(CLOGICAL),BIND(C, name="link_bottom_drag") :: link_bottom_drag
   TYPE(CLOGICAL),BIND(C, name="ice")              :: ice

   TYPE(CINTEGER),BIND(C, name="split_factor")     :: split_factor
   TYPE(CINTEGER),BIND(C, name="ode_method")       :: ode_method
   TYPE(CINTEGER),BIND(C, name="benthic_mode")     :: benthic_mode

   TYPE(AED_REAL),BIND(C, name="rain_factor") :: rain_factor
   TYPE(AED_REAL),BIND(C, name="sw_factor")   :: sw_factor
   TYPE(AED_REAL),BIND(C, name="friction")    :: friction

   TYPE(AED_REAL),BIND(C, name="Kw") :: Kw
   TYPE(AED_REAL),BIND(C, name="dt") :: dt

   TYPE(AED_REAL),TARGET,BIND(C, name="yearday")   :: yearday
   TYPE(AED_REAL),TARGET,BIND(C, name="timestep")  :: timestep
   TYPE(AED_REAL),TARGET,BIND(C, name="Longitude") :: longitude
   TYPE(AED_REAL),TARGET,BIND(C, name="Latitude")  :: latitude

CONTAINS

!#
!# These are 2 useful routines for converting between fortran and C strings
!# They are in here because, well, I guess they are sort of type conversions
!#

!###############################################################################
SUBROUTINE make_string(s1,s2,len)
   CHARACTER(len=*),INTENT(out) :: s1
   CCHARACTER,INTENT(in) :: s2(*)
   CSIZET,INTENT(in)    :: len
!LOCALS
   CHARACTER(len=len) :: s3
!
!-------------------------------------------------------------------------------
!BEGIN
   s1 = trim(transfer(s2(1:len),s3))
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


!###############################################################################
SUBROUTINE glm_init_fortran_support() BIND(C, name="glm_init_fortran_support")
!
!-------------------------------------------------------------------------------
!BEGIN
   CALL C_F_POINTER(cLake, theLake, [MaxLayers]);
   CALL C_F_POINTER(cMetData, aMetData, [1])
   CALL C_F_POINTER(cMetData, MetData)
   CALL C_F_POINTER(cSurfData, aSurfData, [1])
   CALL C_F_POINTER(cSurfData, SurfData)
   CALL C_F_POINTER(cZones, theZones, [n_zones]);
END SUBROUTINE glm_init_fortran_support
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

END MODULE glm_types
