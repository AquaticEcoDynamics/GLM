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
 */
int  read_glm_restart(const char *fn);

#endif
