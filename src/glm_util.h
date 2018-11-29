/******************************************************************************
 *                                                                            *
 * glm_util.h                                                                 *
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
#ifndef _GLM_UTIL_H_
#define _GLM_UTIL_H_

#include "glm.h"

AED_REAL zenith_angle(AED_REAL lon, AED_REAL lat, int day, int iclock, AED_REAL TZ);

// #ifndef sqr
//    AED_REAL sqr(AED_REAL x);
// #endif

// #ifndef gprime
//    AED_REAL gprime(AED_REAL d1, AED_REAL d2);
// #endif

AED_REAL combine(AED_REAL c1, AED_REAL v1, AED_REAL d1,
                    AED_REAL c2, AED_REAL v2, AED_REAL d2);

// #ifndef combine_vol
//    AED_REAL combine_vol(AED_REAL c1, AED_REAL v1, AED_REAL c2, AED_REAL v2);
// #endif

AED_REAL calculate_density(AED_REAL temp, AED_REAL salt);
AED_REAL saturated_vapour(AED_REAL temp);

int internal_var(const char *name);

#endif
