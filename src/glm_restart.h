/******************************************************************************
 *                                                                            *
 * glm_restart.h                                                              *
 *                                                                            *
 * Developed by :                                                             *
 *     AquaticEcoDynamics (AED) Group                                         *
 *     School of Agriculture and Environment                                  *
 *     The University of Western Australia                                    *
 *                                                                            *
 *     http://aquatic.science.uwa.edu.au/                                     *
 *                                                                            *
 * Copyright 2013-2026 : The University of Western Australia                  *
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
#ifndef _GLM_RESTART_H_
#define _GLM_RESTART_H_

/* Name of the restart file (NULL => no NetCDF restart) */
extern char *restart_fname;
/* Write restart every restart_nsave steps (<=0 => only at end) */
extern int   restart_nsave;

/*
 * write_glm_restart: write all model state to a NetCDF restart file.
 */
void write_glm_restart(const char *fn);

/*
 * read_glm_restart: read model state from a NetCDF restart file.
 * Returns 1 on success, 0 if the file does not exist.
 * Aborts on format/dimension mismatch errors.
 *
 * Does NOT restore PTM particle state - see read_glm_restart_ptm() below.
 * PTM_Stat is not yet allocated (ptm_init_glm() has not run) at the point in
 * initialise_lake() where this is called, so any attempt to load it here
 * would silently find PTM_Stat == NULL and skip particle state entirely.
 */
int  read_glm_restart(const char *fn);

/*
 * read_glm_restart_ptm: read ONLY the PTM particle state (ptm_stat, ptm_vars)
 * from a NetCDF restart file, re-opening it independently of
 * read_glm_restart(). Call this after ptm_init_glm() has allocated PTM_Stat/
 * PTM_Vars (i.e. from init_glm(), not from initialise_lake()) - see the
 * comment on read_glm_restart() above for why the two cannot share one call.
 * Returns 1 on success, 0 if the file does not exist or has no PTM state.
 */
int  read_glm_restart_ptm(const char *fn);

#endif
