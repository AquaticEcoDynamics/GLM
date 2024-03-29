/******************************************************************************
 *                                                                            *
 * glm_sizes.c                                                                *
 *                                                                            *
 *  routines for getting layers sizes (area/volume)                           *
 *                                                                            *
 * Developed by :                                                             *
 *     AquaticEcoDynamics (AED) Group                                         *
 *     School of Agriculture and Environment                                  *
 *     The University of Western Australia                                    *
 *                                                                            *
 *     http://aquatic.science.uwa.edu.au/                                     *
 *                                                                            *
 * Copyright 2013 - 2024 -  The University of Western Australia               *
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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "glm.h"

#include "glm_types.h"
#include "glm_const.h"
#include "glm_globals.h"

#include "glm_util.h"
#include "glm_mixu.h"

//#define dbgprt(...) fprintf(stderr, __VA_ARGS__)
#define dbgprt(...) /* __VA_ARGS__ */


int _storage_index(AED_REAL height, AED_REAL *y)
{
    AED_REAL _y = AMOD((height * MphInc), 1.0);
    int ij = height - _y;
    if (ij >= Nmorph) {
        _y += (ij - Nmorph);
        ij = Nmorph;
    }
    ij--;
    *y = _y;
    return ij;
}

AED_REAL _storage_volume(AED_REAL height)
{
    AED_REAL y = 0.;
    int ij = _storage_index(height, &y);

    if (ij >= 0) return MphLevelVol[ij] + y * dMphLevelVol[ij];
    return MphLevelVol[0];
}

AED_REAL _storage_area(AED_REAL height)
{
    AED_REAL y = 0.;
    int ij = _storage_index(height, &y);

    if (ij >= 0) return MphLevelArea[ij] + y * dMphLevelArea[ij];
    return MphLevelArea[0];
}
