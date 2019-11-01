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

wq_init_glm_t        p_wq_init_glm        = NULL;
wq_set_glm_data_t    p_wq_set_glm_data    = NULL;
wq_do_glm_t          p_wq_do_glm          = NULL;
wq_clean_glm_t       p_wq_clean_glm       = NULL;
wq_init_glm_output_t p_wq_init_glm_output = NULL;
wq_write_glm_t       p_wq_write_glm       = NULL;
wq_var_index_c_t     p_wq_var_index_c     = NULL;
wq_set_flags_t       p_wq_set_flags       = NULL;
wq_is_var_t          p_wq_is_var          = NULL;
wq_set_glm_zones_t   p_wq_set_glm_zones   = NULL;


int ode_method = 1, split_factor = 1;
int benthic_mode = 0;
CLOGICAL bioshade_feedback = TRUE, repair_state = TRUE;
CLOGICAL mobility_off = FALSE;     //  !# flag to turn mobility off


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
#endif

/******************************************************************************
 *                                                                            *
 ******************************************************************************/
int prime_wq(const char *which)
{
#ifndef PLOTS
    CLOGICAL do_plots = FALSE;
#endif
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
    p_wq_set_flags       =       (wq_set_flags_t) find_entry(glm_wq_handle, "wq_set_flags");
    p_wq_is_var          =          (wq_is_var_t) find_entry(glm_wq_handle, "wq_is_var");
    p_wq_set_glm_zones   =   (wq_set_glm_zones_t) find_entry(glm_wq_handle, "wq_set_glm_zones");
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
        p_wq_set_flags       =       (wq_set_flags_t) fabm_set_flags;
        p_wq_is_var          =          (wq_is_var_t) fabm_is_var;
#else
        fprintf(stderr, "FABM not supported in this build\n");
        exit(1);
#endif
    } else {
#ifdef AED2
        p_wq_init_glm        =        (wq_init_glm_t) aed2_init_glm;
        p_wq_set_glm_data    =    (wq_set_glm_data_t) aed2_set_glm_data;
        p_wq_do_glm          =          (wq_do_glm_t) aed2_do_glm;
        p_wq_clean_glm       =       (wq_clean_glm_t) aed2_clean_glm;
        p_wq_init_glm_output = (wq_init_glm_output_t) aed2_init_glm_output;
        p_wq_write_glm       =       (wq_write_glm_t) aed2_write_glm;
        p_wq_var_index_c     =     (wq_var_index_c_t) aed2_var_index_c;
        p_wq_set_flags       =       (wq_set_flags_t) aed2_set_flags;
        p_wq_is_var          =          (wq_is_var_t) aed2_is_var;
#else
        fprintf(stderr, "AED2 not supported in this build\n");
        exit(1);
#endif
    }
#endif

    // This is weird. Comment out the debug fprintf below and the flags come out wrong, leave the debug in and
    // they are OK ....
    fprintf(stderr,
        "     'wq_setup': split_factor %d mobility_off %d bioshade_feedback %d repair_state %d ode_method %d benthic_mode %d do_plots %d\n",
                    split_factor, mobility_off, bioshade_feedback,repair_state, ode_method, benthic_mode, do_plots);

    (*p_wq_set_flags)(&split_factor, &mobility_off, &bioshade_feedback,
                                     &repair_state, &ode_method, &benthic_mode, &do_plots,
                                     &link_rain_loss, &link_solar_shade, &link_bottom_drag);

    return 0;
}
