/******************************************************************************
 *                                                                            *
 * glm_output.c                                                               *
 *                                                                            *
 *  output for glm                                                            *
 *                                                                            *
 * Developed by :                                                             *
 *     AquaticEcoDynamics (AED) Group                                         *
 *     School of Agriculture and Environment                                  *
 *     The University of Western Australia                                    *
 *                                                                            *
 *     http://aquatic.science.uwa.edu.au/                                     *
 *                                                                            *
 * Copyright 2013 - 2019 -  The University of Western Australia               *
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
#include <stdlib.h>
#include <string.h>
#include <sys/types.h>
#include <sys/stat.h>
#ifndef _WIN32
  #include <unistd.h>
#else
  #include <direct.h>
  #define S_ISDIR(mode) (mode & _S_IFDIR)
  #define mkdir(path, mode) _mkdir(path)
#endif


#include "glm.h"

#include "glm_types.h"
#include "glm_const.h"
#include "glm_globals.h"

#include "aed_time.h"
#include "aed_csv.h"
#include "glm_csv.h"
#include "glm_ncdf.h"
#include "glm_wqual.h"
#include "glm_plot.h"

#ifdef PLOTS
#include "libplot.h"
#endif

/* for WQ interface */
void wq_write_glm(int ncid, int wlev, int nlev, int *lvl, int point_nlevs)
{ wq_write_glm_(&ncid, &wlev, &nlev, lvl, &point_nlevs); }


extern AED_REAL XLW, XCO, XEV, QSW;

static int plot_id[10];

/******************************************************************************
 * Initialise output streams                                                  *
 ******************************************************************************/
void init_output(int jstart, const char *out_dir, const char *out_fn,
                   int oMaxLayers, AED_REAL Longitude, AED_REAL Latitude)
{
    char ts[20];
    char path[1024];
    struct stat sb;
    extern int startTOD;

    if ( out_dir != NULL && stat(out_dir, &sb) ) {
        fprintf(stderr, "Directory \"%s\" does not exist - attempting to create it\n", out_dir);
        if ( mkdir(out_dir, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH) ) {
            fprintf(stderr, "mkdir failed\n");
            exit(1);
        }
    } else {
        if ( ! S_ISDIR(sb.st_mode) ) {
            fprintf(stderr, "Name given in out_dir (%s) is not a directory\n", out_dir);
            exit(1);
        }
    }

    MaxLayers = oMaxLayers;
    write_time_string(ts,jstart,startTOD);
    snprintf(path, 1024, "%s/%s.nc", out_dir, out_fn);
    ncid = init_glm_ncdf(path, "glm run", Latitude, Longitude, MaxLayers, ts);

    init_csv_output(out_dir);

    //# Initialize WQ output (creates NetCDF variables)
    if (wq_calc) wq_init_glm_output(&ncid, &x_dim, &y_dim, &z_dim, &zone_dim, &time_dim);

#ifdef PLOTS
    if ( do_plots ) {
        int i;
        for (i = 0; i < 10; i++) plot_id[i] = -1;
        init_plots(jstart,nDays,MaxHeight);
    }
#endif
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
static AED_REAL min_temp(LakeDataType *Lake, int count)
{
    int i;
    AED_REAL min;

    min = Lake[0].Temp;
    for (i = 1; i < count; i++)
        if ( min > Lake[i].Temp ) min = Lake[i].Temp;
    return min;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
static AED_REAL max_temp(LakeDataType *Lake, int count)
{
    int i;
    AED_REAL max;

    max = Lake[0].Temp;
    for (i = 1; i < count; i++)
        if ( max < Lake[i].Temp ) max = Lake[i].Temp;
    return max;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
static AED_REAL sum_lake_layervol()
{
    AED_REAL sum = 0.;
    int i;
    for (i = 0; i < NumLayers; i++)
        sum += Lake[i].LayerVol;
    return sum;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
static AED_REAL max_dtdz_at(LakeDataType *Lake, int count)
{
    int  i, max_at;
    AED_REAL dtdz, max_dtdz;

    if (count < 2)
        return Lake[0].Height / 2;

    max_at = 1;
    max_dtdz = (Lake[1].Temp - Lake[0].Temp) / (Lake[1].Height - Lake[0].Height);
    for (i = 2; i < count; i++) {
        dtdz = (Lake[i].Temp - Lake[i-1].Temp) / (Lake[i].Height - Lake[i-1].Height);
        if ( max_dtdz < dtdz ) {
            max_dtdz = dtdz;
            max_at = i;
        }
    }
    return (Lake[max_at].Height - Lake[max_at - 1].Height) / 2;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#ifdef PLOTS
/******************************************************************************
 * Return a number > 0 if regular var, < 0 if sheet var, 0 if not a var       *
 ******************************************************************************/
int _intern_is_var(const char *v)
{
    if ( do_plots ) {
        if (strcasecmp("temp", v) == 0) return 1;
        if (strcasecmp("salt", v) == 0) return 2;
        if (strcasecmp("radn", v) == 0) return 3;
        if (strcasecmp("extc", v) == 0) return 4;
        if (strcasecmp("dens", v) == 0) return 5;
        if (strcasecmp("uorb", v) == 0) return 6;
        if (strcasecmp("taub", v) == 0) return 7;
    }
    return 0;
}
int intern_is_var(int id, const char *v)
{
    int vid = _intern_is_var(v);
    if ( !vid ) return 0;
    plot_id[abs(vid)-1] = id;
    return vid;
}
#endif
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/******************************************************************************
 * Write the various output files and/or plots.                               *
 ******************************************************************************/
void write_output(int jday, int iclock, int nsave, int stepnum)
{
    int i, n, lvl[MaxPointCSV];
    char ts[20];

    if ( csv_point_nlevs > 0 ) {
        for (i = 0; i < MaxPointCSV; i++) lvl[i] = -1;

        //# Output at end of time step so add noSecs to time string
        write_time_string(ts, jday, iclock + noSecs);
        for (i = 0; i < csv_point_nlevs; i++) {
            write_csv_point(i, "time", 0.0, ts, FALSE);

            //# find which level csv_at is in
            if (csv_point_frombot[i]) { //Measure as height from bottom layer
                for (n = 0; n < NumLayers; n++) {
                    if ( Lake[n].Height >= csv_point_at[i] ) {
                        lvl[i] = n;
                        break;
                    }
                }
            } else { //Measure as depth from surface
                for (n = NumLayers-1; n >= botmLayer; n--) {
                    if ( (Lake[surfLayer].Height - Lake[n].Height) >= csv_point_at[i] ) {
                        lvl[i] = n + 1;
                        break;
                    }
                    //If reached bottom layer
                    if (lvl[i] == -1) lvl[i] = botmLayer;
                }
            }

            write_csv_point(i, "temp", Lake[lvl[i]].Temp,     NULL, FALSE);
            write_csv_point(i, "salt", Lake[lvl[i]].Salinity, NULL, FALSE);
        }
    }

#ifdef PLOTS
    if ( do_plots ) do_internal_plots(plot_id);
#endif

    write_glm_ncdf(ncid, NumLayers, MaxLayers, stepnum, timestep);

    //# outputs WQ vars to NetCDF
    if (wq_calc)
        wq_write_glm(ncid, NumLayers, MaxLayers, lvl, csv_point_nlevs);

    if (csv_point_nlevs > 0) {
        for (i = 0; i < csv_point_nlevs; i++)
            write_csv_point(i, "", 0.0, NULL, TRUE);
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 * Write the various output files and/or plots.                               *
 ******************************************************************************/
void write_diags(int jday, AED_REAL LakeNum)
{
    char ts[20];
    extern AED_REAL Hs, L, T;

    if ( csv_lake_file < 0 ) return;

    AED_REAL lake_level = Lake[surfLayer].Height;
    if (lake_level<0.011) lake_level=zero;

    //# Output at end of day
    write_time_string(ts, jday, SecsPerDay);

    write_csv_lake("Time",            0.0,                       ts,   FALSE);
    write_csv_lake("Volume",          sum_lake_layervol(),       NULL, FALSE);
    write_csv_lake("Vol Snow",        SurfData.delzSnow*Lake[surfLayer].LayerArea*SurfData.RhoSnow/1e3, NULL, FALSE);
    //# Magic numbers for ice density are from glm_surface.c
    write_csv_lake("Vol Blue Ice",    SurfData.delzBlueIce*Lake[surfLayer].LayerArea*917.0/1e3, NULL, FALSE);
    write_csv_lake("Vol White Ice",   SurfData.delzWhiteIce*Lake[surfLayer].LayerArea*890.0/1e3, NULL, FALSE);
    write_csv_lake("Tot Inflow Vol",  SurfData.dailyInflow,      NULL, FALSE);
    write_csv_lake("Tot Outflow Vol", SurfData.dailyOutflow,     NULL, FALSE);
    write_csv_lake("Overflow Vol",    SurfData.dailyOverflow,    NULL, FALSE);
    write_csv_lake("Evaporation",     SurfData.dailyEvap,        NULL, FALSE);
    write_csv_lake("Rain",            SurfData.dailyRain,        NULL, FALSE);
    write_csv_lake("Local Runoff",    SurfData.dailyRunoff,      NULL, FALSE);
    write_csv_lake("Snowfall",        SurfData.dailySnow,        NULL, FALSE);
    write_csv_lake("Lake Level",      lake_level,                NULL, FALSE);
    write_csv_lake("Surface Area",    Lake[surfLayer].LayerArea, NULL, FALSE);
    write_csv_lake("Blue Ice Thickness", SurfData.delzBlueIce,   NULL, FALSE);
    write_csv_lake("White Ice Thickness",SurfData.delzWhiteIce,  NULL, FALSE);
    write_csv_lake("Snow Thickness",  SurfData.delzSnow,         NULL, FALSE);
    write_csv_lake("Snow Density",    SurfData.RhoSnow,          NULL, FALSE);
    write_csv_lake("Albedo",          SurfData.albedo,           NULL, FALSE);
    write_csv_lake("Max Temp",        max_temp(Lake, NumLayers), NULL, FALSE);
    write_csv_lake("Min Temp",        min_temp(Lake, NumLayers), NULL, FALSE);
    write_csv_lake("Surface Temp",    Lake[surfLayer].Temp,      NULL, FALSE);
    write_csv_lake("Daily Qsw",       SurfData.dailyQsw / Lake[surfLayer].LayerArea/SecsPerDay, NULL, FALSE);
    write_csv_lake("Daily Qe",        SurfData.dailyQe / Lake[surfLayer].LayerArea/SecsPerDay, NULL, FALSE);
    write_csv_lake("Daily Qh",        SurfData.dailyQh / Lake[surfLayer].LayerArea/SecsPerDay, NULL, FALSE);
    write_csv_lake("Daily Qlw",       SurfData.dailyQlw / Lake[surfLayer].LayerArea/SecsPerDay, NULL, FALSE);
    write_csv_lake("Light",           Lake[surfLayer].Light,     NULL, FALSE);
    write_csv_lake("Benthic Light",   Benthic_Light_pcArea,      NULL, FALSE);

    write_csv_lake("Surface Wave Height",Hs,                     NULL, FALSE);
    write_csv_lake("Surface Wave Length",L,                      NULL, FALSE);
    write_csv_lake("Surface Wave Period",T,                      NULL, FALSE);

    write_csv_lake("LakeNumber",      LakeNum,                   NULL, FALSE);
    write_csv_lake("Max dT/dz",    max_dtdz_at(Lake, NumLayers), NULL, FALSE);
    write_csv_lake("CD",              coef_wind_drag,            NULL, FALSE);
    write_csv_lake("CHE",             coef_wind_chwn,            NULL, FALSE);
    write_csv_lake("z/L",             SurfData.dailyzonL*(noSecs/SecsPerDay), NULL, TRUE);
//    write_csv_lake("coef_wind_drag",  coef_wind_drag,            NULL, TRUE);

}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 * Write the outflow data file with WQ variables.                             *
 ******************************************************************************/
void write_outflow(int of_idx, int jday, AED_REAL DrawHeight, AED_REAL vol,
                                AED_REAL vol_bc, AED_REAL hwBot, AED_REAL hwTop)
{
    char ts[64];
    int i, lvl;
    extern int csv_outlet_allinone, csv_outfl_nvars;
    extern int ofl_wq_idx[];
    static AED_REAL vol_tot, state_of_v[MaxCSVOutVars];

    //# work out which level the outflow is at now.
    for (lvl = botmLayer; lvl <= surfLayer; lvl++)
        if (Lake[lvl].Height >= DrawHeight) break;

    if ( lvl > surfLayer ) vol = 0.; // DrawHeight above the lake top

    else {
        if ( csv_outlet_allinone ) {
            if ( of_idx == 0 ) // initialize
                vol_tot = 0;

            for (i = 0; i < csv_outfl_nvars; i++) {
                state_of_v[i] *= vol_tot;
                state_of_v[i] += (vol * _WQ_Vars(ofl_wq_idx[i], lvl));
            }
            vol_tot += vol;
            for (i = 0; i < csv_outfl_nvars; i++)
                state_of_v[i] /= vol_tot;

            if ( of_idx != MaxOut ) return;
            of_idx = 0;
        } else if (wq_calc) {
            for (i = 0; i < csv_outfl_nvars; i++)
                if ( ofl_wq_idx[i] >= 0 ) state_of_v[i] = _WQ_Vars(ofl_wq_idx[i], lvl);
        }
    }

    write_time_string(ts, jday, 0);

    write_csv_outfl(of_idx, "time",       0.0,                     ts,   FALSE);
    write_csv_outfl(of_idx, "flow",       vol,                     NULL, FALSE);

    if (vol > 0. ) {
        write_csv_outfl(of_idx, "Temp",   Lake[lvl].Temp,          NULL, FALSE);
        write_csv_outfl(of_idx, "Salt",   Lake[lvl].Salinity,      NULL, FALSE);

        if (wq_calc) {   //# must do each of the WQ vars
            // # the first 3 vars are flow, temp and salt
            for (i = 3; i < csv_outfl_nvars; i++)
                write_csv_outfl_idx(of_idx, i,  state_of_v[i],     NULL, FALSE);
        }

        write_csv_outfl(of_idx, "hbot",      hwBot,                NULL, FALSE);
        write_csv_outfl(of_idx, "htop",      hwTop,                NULL, FALSE);
        write_csv_outfl(of_idx, "flbc",      vol_bc,               NULL, FALSE);
    }

    //# force a newline
    write_csv_outfl(of_idx, "",           0.0,                     NULL, TRUE);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 * Close files and plots                                                      *
 ******************************************************************************/
void close_output()
{
    close_glm_ncdf(ncid);
    glm_close_csv_output();

#ifdef PLOTS
    if ( do_plots ) {
        extern char *all_plots_name;
        if ( saveall > 1 && all_plots_name ) {
            save_all_plots_named(all_plots_name);
            saveall = 0;
        }
        do_cleanup(saveall);
    }
#endif
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
