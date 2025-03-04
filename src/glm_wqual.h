/******************************************************************************
 *                                                                            *
 * glm_wqual.h                                                                *
 *                                                                            *
 * The interface between glm and water quality code                           *
 *                                                                            *
 * Developed by :                                                             *
 *     AquaticEcoDynamics (AED) Group                                         *
 *     School of Agriculture and Environment                                  *
 *     The University of Western Australia                                    *
 *                                                                            *
 *     http://aquatic.science.uwa.edu.au/                                     *
 *                                                                            *
 * Copyright 2013 - 2025 -  The University of Western Australia               *
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
#ifndef _GLM_WQUAL_H_
#define _GLM_WQUAL_H_

#ifndef _FORTRAN_SOURCE_

typedef void (*wq_init_glm_t)(char *fname, size_t *len, int *NumWQVars, int *NumWQBen);
typedef void (*wq_set_glm_data_t)(void);
typedef void (*wq_do_glm_t)(int *wlev, CLOGICAL *pIce);
typedef void (*wq_clean_glm_t)(void);
typedef void (*wq_init_glm_output_t)(int *ncid, int *x_dim, int *y_dim, int *z_dim,
                                                          int *zone_dim, int *time_dim);
typedef void (*wq_write_glm_t)(int *ncid, int *wlev, int *nlev, int *lvl, int *point_nlevs);
typedef int  (*wq_var_index_c_t)(const char*name, size_t *len);
typedef int (*wq_is_var_t)(int *id, const char *v, size_t *len);
typedef void (*wq_set_glm_zones_t)(int *numVars, int *numBenV, int *numDiagV, int *numDiagHzV);
typedef void (*wq_ZSoilTemp_t)(ZoneType *zone);

typedef void (*wq_inflow_update_t)(AED_REAL *wqinf, int *nwqVars, AED_REAL *temp, AED_REAL *salt);


extern wq_init_glm_t        p_wq_init_glm;
extern wq_set_glm_data_t    p_wq_set_glm_data;
extern wq_do_glm_t          p_wq_do_glm;
extern wq_clean_glm_t       p_wq_clean_glm;
extern wq_init_glm_output_t p_wq_init_glm_output;
extern wq_write_glm_t       p_wq_write_glm;
extern wq_var_index_c_t     p_wq_var_index_c;
extern wq_is_var_t          p_wq_is_var;
extern wq_ZSoilTemp_t       p_wq_ZSoilTemp;
extern wq_inflow_update_t   p_wq_inflow_update;

#define wq_init_glm        (*p_wq_init_glm)
#define wq_set_glm_data    (*p_wq_set_glm_data)
#define wq_do_glm          (*p_wq_do_glm)
#define wq_clean_glm       (*p_wq_clean_glm)
#define wq_init_glm_output (*p_wq_init_glm_output)
#define wq_write_glm_      (*p_wq_write_glm)
#define wq_var_index_c     (*p_wq_var_index_c)
#define wq_is_var          (*p_wq_is_var)
#define ZSoilTemp          (*p_wq_ZSoilTemp)
#define wq_inflow_update   (*p_wq_inflow_update)

int prime_wq(const char *which);

#if USE_DL_LOADER
extern wq_set_glm_zones_t p_wq_set_glm_zones;
#define wq_set_glm_zones  (*p_wq_set_glm_zones)

void wq_init_glm(char *fname, size_t *len, int *NumWQVars, int *NumWQBen);
void wq_set_glm_data(void);
void wq_do_glm(int *wlev, CLOGICAL *pIce);
void wq_clean_glm(void);
void wq_init_glm_output(int *ncid, int *x_dim, int *y_dim, int *z_dim, int *zone_dim, int *time_dim);
void wq_write_glm_(int *ncid, int *wlev, int *nlev, int *lvl, int *point_nlevs);
int  wq_var_index_c(const char*name, size_t *len);
int wq_is_var(int *id, const char *v, size_t *len);
void wq_inflow_update(AED_REAL *wqinf, int *nwqVars, AED_REAL *temp, AED_REAL *salt);

#else

#if FABM
void fabm_init_glm(char *fname, size_t *len, int *NumWQVars, int *NumWQBen);
void fabm_set_glm_data(void);
void fabm_do_glm(int *wlev, CLOGICAL *pIce);
void fabm_clean_glm(void);
void fabm_init_glm_output(int *ncid, int *x_dim, int *y_dim, int *z_dim, int *zone_dim, int *time_dim);
void fabm_write_glm(int *ncid, int *wlev, int *nlev, int *lvl, int *point_nlevs);
int  fabm_var_index_c(const char*name, size_t *len);
int  fabm_is_var(int *id, const char *v, size_t *len);
#endif

#if API
void api_init_glm(char *fname, size_t len, int *NumWQVars, int *NumWQBen);
void api_set_glm_data(void);
void api_set_glm_where(AED_REAL *Longitude, AED_REAL *Latitude, AED_REAL *yearday, AED_REAL *timestep);
void api_do_glm(int *wlev, CLOGICAL *pIce);
void api_clean_glm(void);
void api_init_glm_output(int *ncid, int *x_dim, int *y_dim, int *z_dim, int *zone_dim, int *time_dim);
void api_write_glm(int *ncid, int *wlev, int *nlev, int *lvl, int *point_nlevs);
int  api_var_index_c(const char*name, size_t *len);
int  api_is_var(int *id, const char *v, size_t *len);
void api_update_inflow_wq(AED_REAL *wqinf, int *nwqVars, AED_REAL *temp, AED_REAL *salt);
void api_set_glm_ptm(int *num_particle_groups, int *max_particle_num);
#endif

#if AED
void aed_init_glm(char *fname, size_t *len, int *NumWQVars, int *NumWQBen);
void aed_set_glm_data(void);
void aed_do_glm(int *wlev, CLOGICAL *pIce);
void aed_clean_glm(void);
void aed_init_glm_output(int *ncid, int *x_dim, int *y_dim, int *z_dim, int *zone_dim, int *time_dim);
void aed_write_glm(int *ncid, int *wlev, int *nlev, int *lvl, int *point_nlevs);
int  aed_var_index_c(const char*name, size_t *len);
int  aed_is_var(int *id, const char *v, size_t *len);
void aed_update_inflow_wq(AED_REAL *wqinf, int *nwqVars, AED_REAL *temp, AED_REAL *salt);
#endif

#if AED2
void aed2_init_glm(char *fname, size_t *len, int *NumWQVars, int *NumWQBen);
void aed2_set_glm_data(void);
void aed2_do_glm(int *wlev, CLOGICAL *pIce);
void aed2_clean_glm(void);
void aed2_init_glm_output(int *ncid, int *x_dim, int *y_dim, int *z_dim, int *zone_dim, int *time_dim);
void aed2_write_glm(int *ncid, int *wlev, int *nlev, int *lvl, int *point_nlevs);
int  aed2_var_index_c(const char*name, size_t *len);
int  aed2_is_var(int *id, const char *v, size_t *len);
#endif

void InitialTemp(int *m, const AED_REAL *depth, const AED_REAL *wv,
                         const AED_REAL *topTemp, const AED_REAL *botTemp,
                         const AED_REAL *nSPinUpDays, AED_REAL *tNew);
void zZSoilTemp(ZoneType *zone);
void SoilTemp(int *m, const AED_REAL *depth, const AED_REAL *wv,
                      const AED_REAL *topTemp, AED_REAL *temp, const AED_REAL *heatflux);

#endif

//void wq_set_glm_zones(int *numVars, int *numBenV, int *numDiagV, int *numDiagHzV);

extern int ode_method;
extern int split_factor;
extern CLOGICAL bioshade_feedback;
extern CLOGICAL repair_state;
extern int benthic_mode;
extern CLOGICAL no_zones;
extern CLOGICAL do_plots;

#endif

#endif
