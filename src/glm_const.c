/******************************************************************************
 *                                                                            *
 * glm_const.c                                                                *
 *                                                                            *
 * Declaration of constants                                                   *
 *                                                                            *
 * Developed by :                                                             *
 *     AquaticEcoDynamics (AED) Group                                         *
 *     School of Earth & Environment                                          *
 *     The University of Western Australia                                    *
 *                                                                            *
 *     http://aed.see.uwa.edu.au/                                             *
 *                                                                            *
 * Copyright 2013 - 2016 -  The University of Western Australia               *
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

const AED_REAL Visc    =  0.00000114;

const AED_REAL AreaFactor = 1.0e6; //* Multiplicative factor to get area to m**2
const AED_REAL MLday2m3sec = 1.0e3/86400.0; //* Multiplicative factor to get ML/day to m**3/s
const AED_REAL ML2m3       = 1.0e3; //* Multiplicative factor to get ML to m**3
const AED_REAL rho0		   = 1.0e3; //* Density of water standard to convert specific density to density in kg/m3

const AED_REAL g          = 9.81;

const AED_REAL Pi         = PI;
const AED_REAL two_Pi     = (2.*PI);
const AED_REAL halfPi     = (PI/2.);
const AED_REAL PiDeg      = 180.0;

const AED_REAL deg2rad = PI/180.;
const AED_REAL rad2deg = 180./PI;

const AED_REAL SecsPerDay = 86400.0;

const AED_REAL missing    = MISVAL;

const AED_REAL Kelvin = 273.15;     // kelvin to celsius conversion
const AED_REAL Latent_Heat_Evap = 2.453E+6;  // Latent heat of evaporation J/kg
const AED_REAL SPHEAT = 4185.5;      // Specific heat of water  J/(kg·K) (15 °C, 101.325 kPa)
const AED_REAL Rspecific = 287.058; // Gas constant  J/kg/K
const AED_REAL Stefan_Boltzman = 5.67E-8;  //# Stefan-Boltzman constant

const AED_REAL zero    =  0.0;
