/******************************************************************************
 *                                                                            *
 * glm_mixu.h                                                                 *
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
#ifndef _GLM_MIXU_H_
#define _GLM_MIXU_H_

void add_this_layer(AED_REAL *VMsum, AED_REAL *Tsum,
                    AED_REAL *Ssum, AED_REAL *VMLOC, AED_REAL *TMLOC,
                    AED_REAL *SMLOC, AED_REAL *DFLOC, int idx);
void average_layer(int *j1, int *k1,
                           AED_REAL MeanTemp, AED_REAL MeanSalt, AED_REAL Dens);
void resize_internals(int icode, int lnu);

#endif
