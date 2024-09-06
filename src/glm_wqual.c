/******************************************************************************
 *                                                                            *
 * glm_wqual.c                                                                *
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
#include <string.h>
#include <stdlib.h>
#include <ctype.h>

#include "glm.h"
#include "glm_types.h"
#include "glm_wqual.h"

#include "glm_globals.h"
#include "glm_ncdf.h"
#include "glm_csv.h"
#include "glm_plot.h"
#include "glm_mobl.h"

#if USE_DL_LOADER
#include <dlfcn.h>
#ifdef _WIN32
#include "glm_plugin.h"
set_funcs_t          p_set_funcs          = NULL;
#endif

static void *glm_wq_handle = NULL;
#endif
static void dummy_inflow_update(AED_REAL *wqinf, int *nwqVars, AED_REAL *temp, AED_REAL *salt);

wq_init_glm_t        p_wq_init_glm        = NULL;
wq_set_glm_data_t    p_wq_set_glm_data    = NULL;
wq_do_glm_t          p_wq_do_glm          = NULL;
wq_clean_glm_t       p_wq_clean_glm       = NULL;
wq_init_glm_output_t p_wq_init_glm_output = NULL;
wq_write_glm_t       p_wq_write_glm       = NULL;
wq_var_index_c_t     p_wq_var_index_c     = NULL;
wq_is_var_t          p_wq_is_var          = NULL;
wq_set_glm_zones_t   p_wq_set_glm_zones   = NULL;
wq_ZSoilTemp_t       p_wq_ZSoilTemp       = NULL;
wq_inflow_update_t   p_wq_inflow_update   = (wq_inflow_update_t)dummy_inflow_update;

void glm_init_fortran_support(void);

#if USE_DL_LOADER
/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void *find_entry(void *glm_wq_handle, const char *entry)
{
    void *ret = dlsym(glm_wq_handle, entry);
    if ( ret == NULL ) {
        fputs(dlerror(), stderr); fputc('\n', stderr);
        exit(1);
    }
    return ret;
}

void dummyZSoilTemp(ZoneType *zone)
{
}
#endif
static void dummy_inflow_update(AED_REAL *wqinf, int *nwqVars, AED_REAL *temp, AED_REAL *salt)
{
}

/******************************************************************************
 *                                                                            *
 ******************************************************************************/
int prime_wq(const char *which)
{
#ifndef PLOTS
    CLOGICAL do_plots = FALSE;
#endif
    glm_init_fortran_support();
#if USE_DL_LOADER
    char dirname[1024];
    char *wq_name = NULL;
    char *ts = NULL;

    if ( glm_wq_handle != NULL )
        return 0;

    dirname[0] = 0;
    wq_name = strncat(dirname, "libglm_wq_", 1024);
    ts = &wq_name[strlen(wq_name)];
    strcat(wq_name, which);
    while (*ts) { *ts = tolower(*ts); ts++; }
#ifdef _WIN32
    strcat(wq_name, ".dll");
#else
    strcat(wq_name, ".so");
#endif

    fprintf(stderr, "Attempting to load WQ code from '%s'\n", wq_name);
    if ( (glm_wq_handle = dlopen(wq_name, RTLD_NOW)) == NULL ) {
        fputs(dlerror(), stderr); fputc('\n', stderr);
        exit(1);
    }

    p_wq_init_glm        =        (wq_init_glm_t) find_entry(glm_wq_handle, "wq_init_glm");
    p_wq_set_glm_data    =    (wq_set_glm_data_t) find_entry(glm_wq_handle, "wq_set_glm_data");
    p_wq_do_glm          =          (wq_do_glm_t) find_entry(glm_wq_handle, "wq_do_glm");
    p_wq_clean_glm       =       (wq_clean_glm_t) find_entry(glm_wq_handle, "wq_clean_glm");
    p_wq_init_glm_output = (wq_init_glm_output_t) find_entry(glm_wq_handle, "wq_init_glm_output");
    p_wq_write_glm       =       (wq_write_glm_t) find_entry(glm_wq_handle, "wq_write_glm");
    p_wq_var_index_c     =     (wq_var_index_c_t) find_entry(glm_wq_handle, "wq_var_index_c");
    p_wq_is_var          =          (wq_is_var_t) find_entry(glm_wq_handle, "wq_is_var");
    p_wq_set_glm_zones   =   (wq_set_glm_zones_t) find_entry(glm_wq_handle, "wq_set_glm_zones");
    p_wq_ZSoilTemp       =       (wq_ZSoilTemp_t) find_entry(glm_wq_handle, "zZSoilTemp");
    if ( p_wq_ZSoilTemp == NULL ) p_wq_ZSoilTemp = (wq_ZSoilTemp_t) dummyZSoilTemp;
    p_wq_inflow_update   =   (wq_inflow_update_t) find_entry(glm_wq_handle, "wq_inflow_update");
#ifdef _WIN32
    p_set_funcs          =          (set_funcs_t) find_entry(glm_wq_handle, "set_funcs");

    (*p_set_funcs)(set_c_wqvars_ptr, Mobility, define_mode_on,
                   define_mode_off, new_nc_variable, set_nc_attributes,
                   store_nc_array, store_nc_scalar, write_csv_point,
                   put_glm_val, put_glm_val_s);
#endif

#else
//  # !USE_DL_LOADER

    if ( strcmp(which, "fabm") == 0 ) {
#ifdef FABM
        p_wq_init_glm        =        (wq_init_glm_t) fabm_init_glm;
        p_wq_set_glm_data    =    (wq_set_glm_data_t) fabm_set_glm_data;
        p_wq_do_glm          =          (wq_do_glm_t) fabm_do_glm;
        p_wq_clean_glm       =       (wq_clean_glm_t) fabm_clean_glm;
        p_wq_init_glm_output = (wq_init_glm_output_t) fabm_init_glm_output;
        p_wq_write_glm       =       (wq_write_glm_t) fabm_write_glm;
        p_wq_var_index_c     =     (wq_var_index_c_t) fabm_var_index_c;
        p_wq_is_var          =          (wq_is_var_t) fabm_is_var;
#else
        fprintf(stderr, "FABM not supported in this build\n");
        exit(1);
#endif
    } else if ( strcmp(which, "aed2") == 0 ) {
#ifdef AED2
        p_wq_init_glm        =        (wq_init_glm_t) aed2_init_glm;
        p_wq_set_glm_data    =    (wq_set_glm_data_t) aed2_set_glm_data;
        p_wq_do_glm          =          (wq_do_glm_t) aed2_do_glm;
        p_wq_clean_glm       =       (wq_clean_glm_t) aed2_clean_glm;
        p_wq_init_glm_output = (wq_init_glm_output_t) aed2_init_glm_output;
        p_wq_write_glm       =       (wq_write_glm_t) aed2_write_glm;
        p_wq_var_index_c     =     (wq_var_index_c_t) aed2_var_index_c;
        p_wq_is_var          =          (wq_is_var_t) aed2_is_var;
#else
        fprintf(stderr, "AED2 not supported in this build\n");
        exit(1);
#endif

        p_wq_ZSoilTemp       =       (wq_ZSoilTemp_t) ZSoilTemp;
    } else if ( strcmp(which, "aed") == 0 ) {
#ifdef AED
        p_wq_init_glm        =        (wq_init_glm_t) aed_init_glm;
        p_wq_set_glm_data    =    (wq_set_glm_data_t) aed_set_glm_data;
        p_wq_do_glm          =          (wq_do_glm_t) aed_do_glm;
        p_wq_clean_glm       =       (wq_clean_glm_t) aed_clean_glm;
        p_wq_init_glm_output = (wq_init_glm_output_t) aed_init_glm_output;
        p_wq_write_glm       =       (wq_write_glm_t) aed_write_glm;
        p_wq_var_index_c     =     (wq_var_index_c_t) aed_var_index_c;
        p_wq_is_var          =          (wq_is_var_t) aed_is_var;
        p_wq_inflow_update   =   (wq_inflow_update_t) aed_update_inflow_wq;
#else
        fprintf(stderr, "AED not supported in this build\n");
        exit(1);
#endif

        p_wq_ZSoilTemp       =       (wq_ZSoilTemp_t) zZSoilTemp;
    } else if ( strcmp(which, "api") == 0 ) {
#ifdef API
        p_wq_init_glm        =        (wq_init_glm_t) api_init_glm;
        p_wq_set_glm_data    =    (wq_set_glm_data_t) api_set_glm_data;
        p_wq_do_glm          =          (wq_do_glm_t) api_do_glm;
        p_wq_clean_glm       =       (wq_clean_glm_t) api_clean_glm;
        p_wq_init_glm_output = (wq_init_glm_output_t) api_init_glm_output;
        p_wq_write_glm       =       (wq_write_glm_t) api_write_glm;
        p_wq_var_index_c     =     (wq_var_index_c_t) api_var_index_c;
        p_wq_is_var          =          (wq_is_var_t) api_is_var;
        p_wq_inflow_update   =   (wq_inflow_update_t) api_update_inflow_wq;
#else
        fprintf(stderr, "API not supported in this build\n");
        exit(1);
#endif
        p_wq_ZSoilTemp       =       (wq_ZSoilTemp_t) zZSoilTemp;
    } else if ( *which != 0 ) {
        fprintf(stderr, "\"%s\" not a water quality module supported in this build\n", which);
        exit(1);
    }
#endif
    if ( p_wq_inflow_update == NULL ) p_wq_inflow_update = (wq_inflow_update_t) dummy_inflow_update;

    return 0;
}
