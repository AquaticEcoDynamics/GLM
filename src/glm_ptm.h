/******************************************************************************
 *                                                                            *
 * glm_ptm.h                                                                  *
 *                                                                            *
 * Contains Particle Tracking Model                                           *
 *                                                                            *
 * Developed by :                                                             *
 *                                                                            *
 * Copyright 2024 - The University of Western Australia                       *
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
#ifndef _GLM_PTM_H_
#define _GLM_PTM_H_

#include "glm.h"

//void ptm_init_glm(char *fname, size_t *len, int *kk, int *NumWQVars, int *NumWQBen, AED_REAL *pKw);
void ptm_init_glm(void);
AED_REAL do_particle_tracking(void);
void ptm_init_glm_output(int ncid, int time_dim);
//void ptm_write_glm(int *ncid, int *wlev, int *nlev, int *lvl, int *point_nlevs);
void ptm_write_glm(int ncid, int ptm_write_glm);

void do_ptm_update(void);

void ptm_redistribute(AED_REAL upper_height, AED_REAL lower_height);
void ptm_addparticles(int new_particles, AED_REAL upper_height, AED_REAL lower_height);
void ptm_layershift(AED_REAL shift_height, AED_REAL shift_amount);
void ptm_update_layerid(void);

// void ptm_destroyparticles()
// void ptm_restart()
extern int ptm_sw;
extern AED_REAL settling_velocity;
extern int init_particle_num;
extern int sed_deactivation;

#endif
