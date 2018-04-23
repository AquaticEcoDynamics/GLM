/******************************************************************************
 *                                                                            *
 * libglm.h                                                                   *
 *                                                                            *
 * Developed by :                                                             *
 *     AquaticEcoDynamics (AED) Group                                         *
 *     School of Agriculture and Environment                                  *
 *     The University of Western Australia                                    *
 *                                                                            *
 *     http://aquatic.science.uwa.edu.au/                                     *
 *                                                                            *
 * Copyright 201 -  The University of Western Australia                       *
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
#ifndef _LIB_GLM_H_
#define _LIB_GLM_H_

#include "glm.h"

/*----------------------------------------------------------------------------*/
typedef void (*init_model_t)(int *jstart, int *nsave);
typedef void (*do_model_coupled_t)(int step_start, int step_end,
	AED_REAL *FlowNew, AED_REAL *DrawNew, AED_REAL *elevation, int nsave);
typedef void (*end_model_t)(void);

#endif
