/******************************************************************************
 *                                                                            *
 * glm_const.h                                                                *
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
#ifndef _GLM_CONST_H_
#define _GLM_CONST_H_

#include "glm.h"

extern const AED_REAL Visc;

extern const AED_REAL AreaFactor;       //# Multiplicative factor to get area to m**2
extern const AED_REAL MLday2m3sec;      //# Multiplicative factor to get ML/day to m**3/s
extern const AED_REAL ML2m3;            //# Multiplicative factor to get ML to m**3
extern const AED_REAL rho0;             //# Density of water standard to convert specific density to density in kg/m3

extern const AED_REAL g;

#define PI 3.14159265358979323846

extern const AED_REAL Pi;
extern const AED_REAL two_Pi;
extern const AED_REAL halfPi;
extern const AED_REAL PiDeg;

extern const AED_REAL deg2rad;
extern const AED_REAL rad2deg;

extern const AED_REAL SecsPerDay;
extern const int      iSecsPerDay;
extern const AED_REAL SecsPerHr;
extern const int      iSecsPerHr;

extern const AED_REAL missing;

extern const AED_REAL Kelvin;           //# kelvin to celsius conversion
extern const AED_REAL Latent_Heat_Evap; //# Latent heat of evaporation J/kg
extern const AED_REAL SPHEAT;           //# Specific heat of water  J/(kg·K) (15 °C, 101.325 kPa)
extern const AED_REAL c_gas;            //# Gas constant  J/kg/K
extern const AED_REAL Stefan_Boltzman;  //# Stefan-Boltzman constant
extern const AED_REAL mwrw2a;
extern const AED_REAL zero;

extern const AED_REAL atm_pressure_sl;
extern const AED_REAL cp_air;

#endif
