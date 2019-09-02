/******************************************************************************
 *                                                                            *
 * glm_output.h                                                               *
 *                                                                            *
 * Developed by :                                                             *
 *     AquaticEcoDynamics (AED) Group                                         *
 *     School of Agriculture and Environment                                  *
 *     The University of Western Australia                                    *
 *                                                                            *
 *     http://aquatic.science.uwa.edu.au/                                     *
 *                                                                            *
 * Copyright 2013 - 2019 -  The University of Western Australia               *
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
#ifndef _GLM_OUTPUT_H_
#define _GLM_OUTPUT_H_

void init_output(int jstart, const char *out_dir, const char *out_fn,
                   int MaxLayers, AED_REAL Longitude, AED_REAL Latitude);
void write_output(int jday, int iclock, int nsave, int stepnum);
void write_diags(int jday, AED_REAL LakeNum);
//void write_outflow(int of_idx, int jday, AED_REAL DrawHeight, AED_REAL vol)
void write_outflow(int of_idx, int jday, AED_REAL DrawHeight,
                 AED_REAL vol, AED_REAL vol_bc, AED_REAL hwBot, AED_REAL hwTop);
void close_output(void);

int intern_is_var(int id, const char *v);

#endif
