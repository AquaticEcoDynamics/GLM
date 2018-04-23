/******************************************************************************
 *                                                                            *
 * glm_plugin.h                                                               *
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
#ifndef _GLM_PLUGIN_H_
#define _GLM_PLUGIN_H_

#include "glm.h"

typedef void (*set_c_wqvars_ptr_t)(AED_REAL *iwqvars);
typedef void (*Mobility_t)(int *N, AED_REAL *dt, AED_REAL *h, AED_REAL *A,
                         AED_REAL *ww, AED_REAL *min_C, AED_REAL *cc);

typedef void (*define_mode_on_t)(int *ncid);
typedef void (*define_mode_off_t)(int *ncid);
typedef int (*new_nc_variable_t)(int ncid, const char *name, int data_type,
                                                      int ndim, const int *dims);
typedef void (*set_nc_attributes_t)(int ncid, int id, const char *units,
                                      const char *long_name, AED_REAL FillValue);
typedef void (*store_nc_array_t)(int ncid, int id, int var_shape, int nvals,
                                                   int maxvals, AED_REAL *array);
typedef void (*store_nc_scalar_t)(int ncid, int id, int var_shape, AED_REAL scalar);

typedef void (*write_csv_point_t)(int f, const char *name, AED_REAL val, const char *cval, int last);
typedef void (*put_glm_val_t)(int plot_id, AED_REAL *val);
typedef void (*put_glm_val_s_t)(int plot_id, AED_REAL *val);

typedef struct _plugin_funcs_ {
    set_c_wqvars_ptr_t  set_c_wqvars_ptr;
    Mobility_t          Mobility;

    define_mode_on_t    define_mode_on;
    define_mode_off_t   define_mode_off;
    new_nc_variable_t   new_nc_variable;
    set_nc_attributes_t set_nc_attributes;
    store_nc_array_t    store_nc_array;
    store_nc_scalar_t   store_nc_scalar;

    write_csv_point_t   write_csv_point;
    put_glm_val_t       put_glm_val;
    put_glm_val_s_t     put_glm_val_s;
} plugin_funcs;

typedef void (*set_funcs_t)(
    set_c_wqvars_ptr_t  set_c_wqvars_ptr,
    Mobility_t          Mobility,
    define_mode_on_t    define_mode_on,
    define_mode_off_t   define_mode_off,
    new_nc_variable_t   new_nc_variable,
    set_nc_attributes_t set_nc_attributes,
    store_nc_array_t    store_nc_array,
    store_nc_scalar_t   store_nc_scalar,
    write_csv_point_t   write_csv_point,
    put_glm_val_t       put_glm_val,
    put_glm_val_s_t     put_glm_val_s);

#ifdef _WIN32
  __declspec(dllexport)
#endif
void set_funcs(
    set_c_wqvars_ptr_t  set_c_wqvars_ptr,
    Mobility_t          Mobility,
    define_mode_on_t    define_mode_on,
    define_mode_off_t   define_mode_off,
    new_nc_variable_t   new_nc_variable,
    set_nc_attributes_t set_nc_attributes,
    store_nc_array_t    store_nc_array,
    store_nc_scalar_t   store_nc_scalar,
    write_csv_point_t   write_csv_point,
    put_glm_val_t       put_glm_val,
    put_glm_val_s_t     put_glm_val_s);

#endif
