/******************************************************************************
 *                                                                            *
 * glm_csv.c                                                                  *
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

#include "glm.h"
#include "glm_globals.h"
#include "glm_const.h"
#include "glm_csv.h"
#include "aed_csv.h"
#include "aed_time.h"
#include "glm_util.h"
#include "glm_wqual.h"

typedef char VARNAME[40];
typedef char FILNAME[80];


/*----------------------------------------------------------------------------*/

int csv_point_nlevs = 0;
int csv_point_nvars = 0;
static int csv_points[MaxPointCSV];
AED_REAL csv_point_at[MaxPointCSV+1];
static char * csv_point_fname = NULL;
CLOGICAL csv_point_frombot[MaxPointCSV+1];
CLOGICAL csv_point_depth_avg[MaxPointCSV+1];
AED_REAL csv_point_zone_upper[MaxPointCSV+1];
AED_REAL csv_point_zone_lower[MaxPointCSV+1];
static VARNAME csv_point_vars[MaxCSVOutVars];

static CLOGICAL csv_point_depth_run[MaxPointCSV+1];
static int csv_point_depth_top[MaxPointCSV+1];
static int csv_point_depth_bot[MaxPointCSV+1];
static int csv_point_depth_cur[MaxPointCSV+1];
static AED_REAL csv_point_depth_tvol[MaxPointCSV+1];
static AED_REAL csv_point_top_scaler[MaxPointCSV+1];
static AED_REAL csv_point_bot_scaler[MaxPointCSV+1];

static char * csv_lake_fname = NULL;
int csv_lake_file = -1;

static int csv_outfls[MaxOut+1];
static char * csv_outfl_fname = NULL;
static char * csv_ovrfl_fname = NULL;
int csv_outlet_allinone;
int csv_outfl_nvars = 0;
static VARNAME csv_outfl_vars[MaxCSVOutVars];
int ofl_wq_idx[MaxCSVOutVars];

#define BUFCHUNK    10240


/*============================================================================*/

/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void configure_csv(int point_nlevs, AED_REAL *point_at, const char *point_fname,
                    int *point_frombot, int point_nvars, int *point_depth_avg,
                    AED_REAL *point_zone_upper, AED_REAL *point_zone_lower,
                    const char *lake_fname)
{
    int i;
    csv_point_nlevs = point_nlevs;
    for (i = 0; i < csv_point_nlevs; i++) {
        csv_point_at[i] =         point_at[i];
        csv_point_frombot[i] =    point_frombot[i];
        if ( point_depth_avg == NULL )
            csv_point_depth_avg[i] =  FALSE;
        else
            csv_point_depth_avg[i] =  point_depth_avg[i];

        if ( point_zone_upper == NULL )
            csv_point_zone_upper[i] = NaN;
        else
            csv_point_zone_upper[i] = point_zone_upper[i];

        if ( point_zone_lower == NULL )
            csv_point_zone_lower[i] = NaN;
        else
            csv_point_zone_lower[i] = point_zone_lower[i];

        csv_point_depth_run[i] = FALSE;
    }
    if ( point_fname != NULL ) csv_point_fname = strdup(point_fname);
    csv_point_nvars = point_nvars;
    if ( lake_fname != NULL ) csv_lake_fname = strdup(lake_fname);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void configure_outfl_csv(int outlet_allinone, const char *outfl_fname,
                                       int outfl_nvars, const char *ovrfl_fname)
{
    csv_outlet_allinone = outlet_allinone;
    if ( outfl_fname != NULL ) csv_outfl_fname = strdup(outfl_fname);
    csv_outfl_nvars = outfl_nvars;
    if ( ovrfl_fname != NULL ) csv_ovrfl_fname = strdup(ovrfl_fname);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void set_csv_point_varname(int which, const char *point_varnam)
{ strncpy(csv_point_vars[which], point_varnam, 40); }
/*----------------------------------------------------------------------------*/
void set_csv_outfl_varname(int which, const char *outfl_varnam)
{ strncpy(csv_outfl_vars[which], outfl_varnam, 40); }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
static void create_outflow_csv(int ofidx, const char *out_dir, const char *fname)
{
    int i;

    csv_outfls[ofidx] = open_csv_output(out_dir, fname);

    csv_header_start(csv_outfls[ofidx]);
    for (i = 0; i < csv_outfl_nvars; i++)
        csv_header_var(csv_outfls[ofidx], csv_outfl_vars[i]);
    csv_header_end(csv_outfls[ofidx]);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void init_csv_output(const char *out_dir)
{
    int i,j;
    char fname[256];

    for (j = 0; j < MaxPointCSV; j++) csv_points[j] = -1;

    if ( csv_point_fname != NULL ) {
        for (j = 0; j < csv_point_nlevs; j++) {
            i = csv_point_at[j];

            if ( i == 0 ) {
                snprintf(fname,20,"%s%s", csv_point_fname, (csv_point_frombot[j])?"ben":"surf");
            } else
                snprintf(fname,20,"%s%d", csv_point_fname, i);
            csv_points[j] = open_csv_output(out_dir, fname);

            csv_header_start(csv_points[j]);
            for (i = 0; i < csv_point_nvars; i++)
                csv_header_var(csv_points[j], csv_point_vars[i]);

            csv_header_end(csv_points[j]);
        }
    }

    if ( csv_lake_fname != NULL ) {
        csv_lake_file = open_csv_output(out_dir, csv_lake_fname);

        csv_header_start(csv_lake_file);
        csv_header_var(csv_lake_file, "Volume");   //, "m3");
        csv_header_var(csv_lake_file, "Vol Snow"); //, "m3");
        csv_header_var(csv_lake_file, "Vol Blue Ice"); //, "m3");
        csv_header_var(csv_lake_file, "Vol White Ice"); //, "m3");
        csv_header_var(csv_lake_file, "Tot Inflow Vol"); //, "m3");
        csv_header_var(csv_lake_file, "Tot Outflow Vol"); //, "m3");
        csv_header_var(csv_lake_file, "Overflow Vol"); //, "m3");
        csv_header_var(csv_lake_file, "Evaporation"); //, "m3");
        csv_header_var(csv_lake_file, "Rain"); //, "m3");
        csv_header_var(csv_lake_file, "Local Runoff"); //, "m3");
        csv_header_var(csv_lake_file, "Snowfall"); //, "m3");
        csv_header_var(csv_lake_file, "Lake Level"); //, "m");
        csv_header_var(csv_lake_file, "Surface Area"); //, "m2");
        csv_header_var(csv_lake_file, "Blue Ice Thickness");
        csv_header_var(csv_lake_file, "Snow Thickness");
        csv_header_var(csv_lake_file, "Snow Density");
        csv_header_var(csv_lake_file, "White Ice Thickness");
        csv_header_var(csv_lake_file, "Albedo");
        csv_header_var(csv_lake_file, "Max Temp");
        csv_header_var(csv_lake_file, "Min Temp");
        csv_header_var(csv_lake_file, "Surface Temp");
        csv_header_var(csv_lake_file, "Daily Qsw");
        csv_header_var(csv_lake_file, "Daily Qe");
        csv_header_var(csv_lake_file, "Daily Qh");
        csv_header_var(csv_lake_file, "Daily Qlw");
        csv_header_var(csv_lake_file, "Light");
        csv_header_var(csv_lake_file, "Benthic Light");
        csv_header_var(csv_lake_file, "Surface Wave Height"); //, "m");
        csv_header_var(csv_lake_file, "Surface Wave Length");
        csv_header_var(csv_lake_file, "Surface Wave Period");
        csv_header_var(csv_lake_file, "LakeNumber");
        csv_header_var(csv_lake_file, "Max dT/dz");
        csv_header_var(csv_lake_file, "CD");
        csv_header_var(csv_lake_file, "CHE");
        csv_header_var(csv_lake_file, "z/L");
        csv_header_end(csv_lake_file);
    } else
        csv_lake_file = -1;

    for (j = 0; j <= MaxOut; j++) csv_outfls[j] = -1;

    if ( csv_outfl_fname != NULL ) {
        if ( csv_outlet_allinone ) {
            create_outflow_csv(0, out_dir, csv_outfl_fname);
        } else {
            for (j = 0; j < NumOut; j++) {
                snprintf(fname,20,"%s%02d", csv_outfl_fname, j);
                create_outflow_csv(j, out_dir, fname);
            }

            if ( csv_ovrfl_fname != NULL )
                create_outflow_csv(MaxOut, out_dir, csv_ovrfl_fname);
        }

        // Now map the wq state var names to their indices
        for (i = 0; i < csv_outfl_nvars; i++) {
            size_t k = strlen(csv_outfl_vars[i]);
            if ( internal_var(csv_outfl_vars[i]) )
                ofl_wq_idx[i] = -1;
            else if ((ofl_wq_idx[i] = wq_var_index_c(csv_outfl_vars[i], &k)) < 0)
                fprintf(stderr, "Cannot find \"%s\" for outflow value\n", csv_outfl_vars[i]);
        }
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void write_csv_lake(const char *name, AED_REAL val, const char *cval, int last)
{ write_csv_var(csv_lake_file, name, val, cval, last); }

/*----------------------------------------------------------------------------*/
void write_csv_point(int p, const char *name, AED_REAL val, const char *cval, int last)
{
    if ( ! csv_point_depth_avg[p] )
        write_csv_var(csv_points[p], name, val, cval, last);
}

/*----------------------------------------------------------------------------*/
void write_csv_point_avg(int p, const char *name, AED_REAL *vals,
                                                    const char *cval, int last)
{
    int i;
    AED_REAL totl, val, tvol, h1, h2, vh = 0., vl = 0.;

    if ( csv_point_depth_avg[p] ) {
        // if no scales, init first...
        if ( !csv_point_depth_run[p] ) {
            csv_point_depth_run[p] = TRUE;

            // Find which layer "point_at" is in
            csv_point_depth_cur[p] = -1;
            if ( csv_point_at[p] <= Lake[0].Height ) {
                csv_point_depth_cur[p] = 0;
            } else {
                for (i = 1; i < NumLayers; i++) {
                    if ( csv_point_at[p] <= Lake[i].Height &&
                         csv_point_at[p] > Lake[i-1].Height ) {
                        csv_point_depth_cur[p] = i;
                    }
                }
            }

            // Find which layer "upper" is in
            csv_point_depth_top[p] = -1;
            if ( csv_point_zone_upper[p] <= Lake[0].Height ) {
                csv_point_depth_top[p] = 0;
            } else {
                for (i = 1; i < NumLayers; i++) {
                    if ( csv_point_zone_upper[p] <= Lake[i].Height &&
                         csv_point_zone_upper[p] > Lake[i-1].Height ) {
                        csv_point_depth_top[p] = i;
                    }
                }
            }

            // Find which layer "lower" is in
            csv_point_depth_bot[p] = -1;
            if ( csv_point_zone_lower[p] <= Lake[0].Height ) {
                csv_point_depth_bot[p] = 0;
            } else {
                for (i = 1; i < NumLayers; i++) {
                    if ( csv_point_zone_lower[p] <= Lake[i].Height &&
                         csv_point_zone_lower[p] > Lake[i-1].Height ) {
                        csv_point_depth_bot[p] = i;
                    }
                }
            }

            // now change val to be the average over the range for this var
            // run through levels from first and last included
            tvol = 0.0;

            // top and bottom layers will only be a proportion of the volume
            i = csv_point_depth_bot[p];
            h1 = Lake[i].Height;
            h2 = csv_point_zone_lower[p];
            if ( i > 0 ) {
                h1 -= Lake[i-1].Height;
                h2 -= Lake[i-1].Height;
            }
            csv_point_bot_scaler[p] = Lake[i].LayerVol * ((h1 - h2) / h1);
            tvol += csv_point_bot_scaler[p];

            // other layers will be full volume contributions
            for (i = csv_point_depth_bot[p]+1; i < csv_point_depth_top[p]; i++) {
                tvol +=  Lake[i].LayerVol;
            }

            // top and bottom layers will only be a proportion of the volume
            i = csv_point_depth_top[p];
            h1 = Lake[i].Height;
            h2 = csv_point_zone_upper[p];
            if ( i > 0 ) {
                h1 -= Lake[i-1].Height;
                h2 -= Lake[i-1].Height;
            }
            csv_point_top_scaler[p] = Lake[i].LayerVol * (h2 / h1);
            tvol += csv_point_top_scaler[p];

            csv_point_depth_tvol[p] = tvol;
        }

        totl = 0.0;
        if ( vals != NULL ) {
            for (i = csv_point_depth_bot[p]+1; i < csv_point_depth_top[p]; i++) {
                totl += (Lake[i].LayerVol * vals[i]);
            }
        } else {
            for (i = csv_point_depth_bot[p]+1; i < csv_point_depth_top[p]; i++) {
                if ( strcmp(name, "temp") == 0 ) {
                    totl += (Lake[i].LayerVol * Lake[i].Temp);
                } else if ( strcmp(name, "salt") == 0 ) {
                    totl += (Lake[i].LayerVol * Lake[i].Salinity);
                } else if ( strcmp(name, "dens") == 0 ) {
                    totl += (Lake[i].LayerVol * Lake[i].Density);
                }
            }
        }

        if ( vals != NULL ) {
            vh = vals[csv_point_depth_top[p]];
            vl = vals[csv_point_depth_bot[p]];
        } else {
            if ( strcmp(name, "temp") == 0 ) {
                vh = Lake[csv_point_depth_top[p]].Temp;
                vl = Lake[csv_point_depth_bot[p]].Temp;
            } else if ( strcmp(name, "salt") == 0 ) {
                vh = Lake[csv_point_depth_top[p]].Salinity;
                vl = Lake[csv_point_depth_bot[p]].Salinity;
            } else if ( strcmp(name, "dens") == 0 ) {
                vh = Lake[csv_point_depth_top[p]].Density;
                vl = Lake[csv_point_depth_bot[p]].Density;
            }
        }

        totl += vh * csv_point_top_scaler[p];
        totl += vl * csv_point_bot_scaler[p];

        val = (totl / csv_point_depth_tvol[p]) ;

        write_csv_var(csv_points[p], name, val, cval, last);

        if ( last ) csv_point_depth_run[p] = FALSE; // reset scales
    }
}

/*----------------------------------------------------------------------------*/
void write_csv_outfl(int ofl, const char *name, AED_REAL val, const char *cval, int last)
{
    if (csv_outfls[ofl] >= 0)
        write_csv_var(csv_outfls[ofl], name, val, cval, last);
}
/*----------------------------------------------------------------------------*/
void write_csv_outfl_idx(int ofl, int var, AED_REAL val, const char *cval, int last)
{
    if (csv_outfls[ofl] >= 0)
        write_csv_var(csv_outfls[ofl], csv_outfl_vars[var], val, cval, last);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void glm_close_csv_output()
{
    int i;

    for (i = 0; i < csv_point_nlevs; i++)
        if ( csv_points[i] >= 0 ) close_csv_output(csv_points[i]);

    if ( csv_lake_file >= 0 ) close_csv_output(csv_lake_file);

    for (i = 0; i < NumOut; i++)
        if ( csv_outfls[i] >= 0 ) close_csv_output(csv_outfls[i]);

    if ( csv_outfls[MaxOut] >= 0 ) close_csv_output(csv_outfls[MaxOut]);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/



/******************************************************************************
 *                                                                            *
 * for Fortran                                                                *
 *                                                                            *
 ******************************************************************************/

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
    write_csv_var(csv_points[*f-1], n, *val, v, *last);
    free(n); free(v);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
void write_csv_point_avg_(int *f, const char *name, int *len, AED_REAL *vals,
                                         const char *cval, int *vlen, int *last)
{
    char *n = make_c_string(name, *len);
    char *v = make_c_string(cval, *vlen);
    write_csv_point_avg(csv_points[*f-1], n, vals, v, *last);
    free(n); free(v);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
