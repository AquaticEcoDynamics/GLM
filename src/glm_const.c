/******************************************************************************
 *                                                                            *
 * glm_const.c                                                                *
 *                                                                            *
 * Declaration of constants                                                   *
 *                                                                            *
 * Developed by :                                                             *
 *     AquaticEcoDynamics (AED) Group                                         *
 *     School of Agriculture and Environment                                  *
 *     The University of Western Australia                                    *
 *                                                                            *
 *     http://aquatic.science.uwa.edu.au/                                     *
 *                                                                            *
 * Copyright 2013 - 2018 -  The University of Western Australia               *
 *                                                                            *
 *  This file is part of GLM (General Lake Model)                             *
 *                                                                            *
 *  GLM is free software: you can redistribute it and/or modify               *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 *  GLM is distributed in the hope that it will be useful,                    *
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of            *
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the             *
 *  GNU General Public License for more details.                              *
 *                                                                            *
 *  You should have received a copy of the GNU General Public License         *
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.     *
 *                                                                            *
 ******************************************************************************/

#include "glm.h"
#include "glm_const.h"

// General constants
const AED_REAL g               = 9.81; // 9.780310.*(1+(0.00530239.*(sin(abs(lat)).^2) - 0.00000587.*(sin(abs(2*lat)).^2) - (31.55e-8)*alt)); % gravitational acceleration, m/s2
const AED_REAL Pi              = PI;
const AED_REAL two_Pi          = (2.*PI);
const AED_REAL halfPi          = (PI/2.);
const AED_REAL Stefan_Boltzman = 5.67E-8;  //# Stefan-Boltzman constant
const AED_REAL Kelvin          = 273.15;     // Used for Celsius to Kelvin conversion

// Water
const AED_REAL rho0            = 1.0e3; //* Density of water standard to convert specific density to density in kg/m3
const AED_REAL Visc            = 0.00000114;
const AED_REAL Latent_Heat_Evap= 2.453E+6;  // Latent heat of evaporation J/kg   // 2.501e6-2370*Ts
const AED_REAL SPHEAT          = 4185.5;      // Specific heat of water  J/(kg·K) (15 °C, 101.325 kPa)

// Air
/******************************************************************************
 * -- Ratio of the molecular (or molar) weight of water to dry air [-]:       *
 *    mwrw2a    =    18.016    /    28.966;      // = 0.62197                 *
 * -- The universal gas constant  [J mol^-1 K^-1] = 8.31436;                  *
 * -- Gas constant for dry air in terms of mass of gas rather than moles      *
 * -- [J kg^-1 K^-1]:                                                         *
 *    c_gas   = 1.0E3 *    8.31436     /    28.966;                           *
 ******************************************************************************/
const AED_REAL mwrw2a          = 18.016 / 28.966;
const AED_REAL c_gas           = 1.0E3 * 8.31436 / 28.966;
const AED_REAL cp_air          = 1005.0;  // Specific heat of air
const AED_REAL atm_pressure_sl = 1013.25;  //# Atmospheric pressure in hectopascals @ sea level ==101300 Pa
//const AED_REAL Rspecific       = 287.058; // Gas constant  J/kg/K

// Factors and conversions (space and time)
const AED_REAL AreaFactor      = 1.0e6; //* Multiplicative factor to get area to m**2
const AED_REAL MLday2m3sec     = 1.0e3/86400.0; //* Multiplicative factor to get ML/day to m**3/s
const AED_REAL ML2m3           = 1.0e3; //* Multiplicative factor to get ML to m**3
const AED_REAL PiDeg           = 180.0;
const AED_REAL deg2rad         = PI/180.;
const AED_REAL rad2deg         = 180./PI;
const AED_REAL SecsPerDay      = 86400.0;
const int      iSecsPerDay     = 86400;
const AED_REAL SecsPerHr       = 3600.0;
const int      iSecsPerHr      = 3600;

// Numeric / other
const AED_REAL missing         = MISVAL;
const AED_REAL zero            = 0.0;
