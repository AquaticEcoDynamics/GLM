/******************************************************************************
 *                                                                            *
 * glm_plugin.c                                                               *
 *                                                                            *
 * This file is only needed for the WIN32 code and only in the DLL water      *
 *  quality modules.  Since they currently dont work anyway, it's not used.   *
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

#ifdef _WIN32

#include <string.h>
#include <stdlib.h>

#include "glm_types.h"

#include "glm_plugin.h"


static plugin_funcs funcs;

char *strndup(const char *s, size_t len);

/******************************************************************************/
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
    put_glm_val_s_t     put_glm_val_s)
{
    funcs.set_c_wqvars_ptr  = set_c_wqvars_ptr;
    funcs.Mobility          = Mobility;

    funcs.define_mode_on    = define_mode_on;
    funcs.define_mode_off   = define_mode_off;
    funcs.new_nc_variable   = new_nc_variable;
    funcs.set_nc_attributes = set_nc_attributes;
    funcs.store_nc_array    = store_nc_array;
    funcs.store_nc_scalar   = store_nc_scalar;

    funcs.write_csv_point   = write_csv_point;
    funcs.put_glm_val       = put_glm_val;
    funcs.put_glm_val_s     = put_glm_val_s;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
void set_c_wqvars_ptr(AED_REAL *iwqvars) { (*funcs.set_c_wqvars_ptr)(iwqvars); }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void Mobility(int *N, AED_REAL *dt, AED_REAL *h, AED_REAL *A,
                         AED_REAL *ww, AED_REAL *min_C, AED_REAL *cc)
{ (*funcs.Mobility)(N, dt, h, A, ww, min_C, cc); }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void define_mode_on(int *ncid) { (*funcs.define_mode_on)(ncid); }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void define_mode_off(int *ncid) { (*funcs.define_mode_off)(ncid); }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
int new_nc_variable_(int *ncid, const char *name, int *len, int *data_type,
                                                     int *ndim, const int *dims)
{
   int i, n = *ndim, ret;
   char *s = strndup(name,*len);
   int dim[6];

   for (i = 0; i < n; i++) dim[n-i-1] = dims[i];
   ret = (*funcs.new_nc_variable)(*ncid, s, *data_type, n, dim);
   free(s);

   return ret;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
void set_nc_attributes_(int *ncid, int *id, const char *units,
                                     const char *long_name, AED_REAL *FillValue)
{ (*funcs.set_nc_attributes)(*ncid, *id, units, long_name, *FillValue); }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void store_nc_scalar_(int *ncid, int *id, int *var_shape, AED_REAL *array)
{ (*funcs.store_nc_scalar)(*ncid, *id, *var_shape, *array); }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void store_nc_array_(int *ncid, int *id, int *var_shape, int *nvals,
                                                  int *maxvals, AED_REAL *array)
{ (*funcs.store_nc_array)(*ncid, *id, *var_shape, *nvals, *maxvals, array); }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
static char *make_c_string(const char *in, int len)
{
    char *t = malloc(len + 1);
    strncpy(t, in, len); t[len] = 0;
    return t;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
void write_csv_point_(int *f, const char *name, int *len, AED_REAL *val,
                              const char *cval, int *vlen, int *last)
{
    char *n = make_c_string(name, *len);
    char *v = make_c_string(cval, *vlen);
    (*funcs.write_csv_point)(*f-1, n, *val, v, *last);
    free(n); free(v);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
void put_glm_val_s_(int *plot_id, AED_REAL *val)
{   (*funcs.put_glm_val_s)(*plot_id, val); }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
void put_glm_val_(int *plot_id, AED_REAL *val)
{   (*funcs.put_glm_val)(*plot_id, val); }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
char *strndup(const char *s, size_t len)
{
    size_t l = strlen(s);
    char *t = malloc(min(l,len)+1);
    if ( t != NULL ) { strncpy(t,s,min(l,len)+1); t[min(l,len)] = 0; }
    return t;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#endif
