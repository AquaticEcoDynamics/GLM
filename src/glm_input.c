/******************************************************************************
 *                                                                            *
 * glm_input.c                                                                *
 *                                                                            *
 * Routines to read some input files.  These will be replaced                 *
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
#include <math.h>
#include <stdlib.h>
#include <string.h>

#include "glm.h"
#include "glm_types.h"
#include "glm_globals.h"
#include "aed_time.h"
#include "glm_util.h"
#include "aed_csv.h"
#include "glm_bird.h"
#include "glm_input.h"


//#define dbgprt(...) fprintf(stderr, __VA_ARGS__)
#define dbgprt(...) /* __VA_ARGS__ */

typedef struct _inf_data_ {
    int inf;  /* Inflow Number */
    int n_vars;   /* number of vars in this file */
    int flow_idx;
    int temp_idx;
    int salt_idx;
    int in_vars[MaxVars];
} InflowDataT;

typedef struct _out_data_ {
    int outf;
    int draw_idx;
} OutflowDataT;

typedef struct _withdrTemp_data_ {
    int withdrTempf;
    int wtemp_idx;
} WithdrawalTempDataT;

static InflowDataT inf[MaxInf];
static OutflowDataT outf[MaxOut];
static WithdrawalTempDataT withdrTempf = { -1, -1 };

static int metf = -1, kwf = -1;
static int rain_idx = -1, hum_idx  = -1, lwav_idx = -1, sw_idx   = -1,
           atmp_idx = -1, wind_idx = -1, snow_idx = -1, rpo4_idx = -1,
           rtp_idx  = -1, rno3_idx = -1, rnh4_idx = -1, rtn_idx  = -1,
           rsi_idx  = -1, wdir_idx = -1, kw_idx   = -1, time_idx = -1;

int lw_ind = 0;
static int have_snow = FALSE, have_rain_conc = FALSE;
static int have_fetch = FALSE;

static int n_steps;

static MetDataType *submet = NULL;
static AED_REAL loaded_day;

AED_REAL wind_factor = 1.0;   //# Windspeed scaling factor
AED_REAL sw_factor   = 1.0;
//AED_REAL lw_factor   = 1.0;
//AED_REAL lw_offset   = 0.0;
AED_REAL at_factor   = 1.0;
AED_REAL at_offset   = 0.0;   //# add this to airtemp in met
AED_REAL rh_factor   = 1.0;
AED_REAL rain_factor = 1.0;


int *WQ_VarsIdx = NULL;

/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void read_daily_inflow(int julian, int NumInf, AED_REAL *flow, AED_REAL *temp,
                                               AED_REAL *salt, AED_REAL *wq)
{
    int csv;
    int i,j,k;

    for (i = 0; i < NumInf; i++) {
        int n_invars = inf[i].n_vars;
        csv = inf[i].inf;
        find_day(csv, time_idx, julian);

        flow[i] = get_csv_val_r(csv,inf[i].flow_idx);
        temp[i] = get_csv_val_r(csv,inf[i].temp_idx);
        salt[i] = get_csv_val_r(csv,inf[i].salt_idx);

        for (j = 0; j < n_invars; j++) {
            if (WQ_VarsIdx[j] < 0) k = j; else k = WQ_VarsIdx[j];
            if (inf[i].in_vars[k] == -1 )
                WQ_INF_(wq, i, k) = 0.;
            else
                WQ_INF_(wq, i, k) = get_csv_val_r(csv,inf[i].in_vars[j]);
        }
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void read_daily_outflow(int julian, int NumOut, AED_REAL *draw)
{
    int csv, i;

    for (i = 0; i < NumOut; i++) {
        csv = outf[i].outf;
        find_day(csv, time_idx, julian);
        draw[i] = get_csv_val_r(csv,outf[i].draw_idx);
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void read_daily_withdraw_temp(int julian, AED_REAL *withdrTemp)
{
    int csv;

    if ( (csv = withdrTempf.withdrTempf) > -1 ) {
        find_day(csv, time_idx, julian);
        *withdrTemp = get_csv_val_r(csv, withdrTempf.wtemp_idx);
    } else
        *withdrTemp = 0.;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void read_daily_met(int julian, MetDataType *met)
{
    int csv, i, idx ;//, err = 0;
    AED_REAL now, tomorrow, t_val, sol;
    AED_REAL eff_area, surf_area, ld, x_ws;

    surf_area = Lake[surfLayer].LayerArea;

//  fprintf(stderr, "READ DAILY MET for %d\n", julian);
    now = julian;
    tomorrow = now + 1.0;
    loaded_day = now;

    dbgprt("read_daily_met (SUBDAY_MET) in\n");

    csv = metf;
    find_day(csv, time_idx, julian);

    for (i = 0; i < n_steps; i++)
        memset(&submet[i], 0, sizeof(MetDataType));

    i = 0;
    while ( (t_val = get_csv_val_r(csv, time_idx)) < tomorrow) {
        if ( i >= n_steps ) {
            int dd,mm,yy;
            calendar_date(now,&yy,&mm,&dd);
            fprintf(stderr, "Warning! Too many steps in met for %4d-%02d-%02d\n", yy,mm,dd);
            break;
        }

        idx = floor((t_val-floor(t_val))*24+1.e-8); // add 1.e-8 to compensate for rounding error
        // fprintf(stderr, "Read met for %16.8f ; %15.12f (%2d)\n", t_val, (t_val-floor(t_val))*24., idx);
//      if ( idx != i ) {
//          if ( !err ) {
//             int dd,mm,yy;
//             calendar_date(now,&yy,&mm,&dd);
//             fprintf(stderr, "Possible sequence issue in met for day %4d-%02d-%02d\n", yy,mm,dd);
//          }
//          idx = i;
//          err = 1;
//      }
        if (idx >= n_steps) {
            int dd,mm,yy;
            calendar_date(now,&yy,&mm,&dd);
            fprintf(stderr, "Step error for %4d-%02d-%02d!\n", yy,mm,dd);
            break;
        }

        // Rain is the exception - goes as is
        submet[idx].Rain        = get_csv_val_r(csv, rain_idx) * rain_factor;
        submet[idx].RelHum      = get_csv_val_r(csv, hum_idx)  * rh_factor;

        if ( submet[idx].RelHum > 100. ) submet[idx].RelHum = 100.;

        if ( lwav_idx != -1 )
            submet[idx].LongWave  = get_csv_val_r(csv, lwav_idx); // * lw_factor;
        else
            submet[idx].LongWave  = 0.;
        if ( sw_idx != -1 )
            submet[idx].ShortWave = get_csv_val_r(csv, sw_idx) * sw_factor;
        else
            submet[idx].ShortWave = 0.;

        switch ( rad_mode ) {
            case 0 : // use the value already read.
            case 1 :
            case 2 :
                break;
            case 3 : // Solar radiation supplied, calculate cloud cover.
            case 4 :
            case 5 :
                sol = calc_bird(Longitude, Latitude, julian, idx*3600, timezone_m);
                sol = sol * sw_factor;
                if ( rad_mode == 4 )
                    sol = clouded_bird(sol, submet[idx].LongWave);
                if ( rad_mode == 3 )
                    submet[idx].LongWave = cloud_from_bird(sol, submet[idx].ShortWave);
                else
                    submet[idx].ShortWave = sol;
                break;
        }

    //  submet[idx].AirTemp     = get_csv_val_r(csv, atmp_idx) * at_factor;
        submet[idx].AirTemp     = get_csv_val_r(csv, atmp_idx) * at_factor + at_offset;
        submet[idx].WindSpeed   = get_csv_val_r(csv, wind_idx) * wind_factor;

        // Read in rel humidity into svd (%), and convert to satvap
        submet[idx].SatVapDef   =  (submet[idx].RelHum/100.) * saturated_vapour(submet[idx].AirTemp);

        if ( have_snow )
             submet[idx].Snow = get_csv_val_r(csv, snow_idx);
        else submet[idx].Snow = 0. ;

        if ( have_rain_conc ) {
            submet[idx].RainConcPO4 = get_csv_val_r(csv, rpo4_idx);
            submet[idx].RainConcTp  = get_csv_val_r(csv, rtp_idx);
            submet[idx].RainConcNO3 = get_csv_val_r(csv, rno3_idx);
            submet[idx].RainConcNH4 = get_csv_val_r(csv, rnh4_idx);
            submet[idx].RainConcTn  = get_csv_val_r(csv, rtn_idx);
            submet[idx].RainConcSi  = get_csv_val_r(csv, rsi_idx);
        } else {
            submet[idx].RainConcPO4 = 0.;
            submet[idx].RainConcTp  = 0.;
            submet[idx].RainConcNO3 = 0.;
            submet[idx].RainConcNH4 = 0.;
            submet[idx].RainConcTn  = 0.;
            submet[idx].RainConcSi  = 0.;
        }


        if ( fetch_mode == 2 || fetch_mode == 3 )
             submet[idx].WindDir = get_csv_val_r(csv, wdir_idx);
        else submet[idx].WindDir = 0. ;

        // wind-sheltering
        switch ( fetch_mode ) {
            case 0 : // no change to the wind-speed
            case 1 : // effective area based on general relationship
                eff_area = surf_area* tanh(surf_area/fetch_aws);
                submet[idx].WindSpeed = submet[idx].WindSpeed * eff_area/surf_area;
                break;
            case 2 : // Markfort et al 2009 model, adapted for wind direction
                ld = LenAtCrest+WidAtCrest/2;
                x_ws = get_fetch(submet[idx].WindDir);
                eff_area = fmax((ld*ld)/2 *acos(x_ws/ld) - (x_ws/2)*pow(ld*ld-x_ws*x_ws,0.5),0.0);
                submet[idx].WindSpeed = submet[idx].WindSpeed * eff_area/surf_area;
                break;
            case 3 : // just scale based on the wind direction directly
                eff_area = get_fetch(submet[idx].WindDir);
                submet[idx].WindSpeed = submet[idx].WindSpeed * eff_area;
                break;
        }

        i++;

        if (!load_csv_line(csv) ) break;
    }

    *met = submet[0];
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 * get_fetch                                                                  *
 ******************************************************************************/
AED_REAL get_fetch(AED_REAL windDir)
{
    AED_REAL fetch_s;
    int i = 0;

    while (windDir > 360.) windDir -= 360.;

    while (i < fetch_ndirs && windDir < fetch_dirs[i])
        i++;

    if ( i > 0 && i < fetch_ndirs ) {
        // the easy case.
        fetch_s = fetch_scale[i-1] +
                   (fetch_scale[i] - fetch_scale[i-1]) *
                    (windDir - fetch_dirs[i-1]) / (fetch_dirs[i] - fetch_dirs[i-1]);
    } else if ( i == 0 ) fetch_s = fetch_scale[0];
    else {
        fetch_s = fetch_scale[fetch_ndirs-1] +
                   (fetch_scale[0] - fetch_scale[fetch_ndirs-1]) *
                    (windDir - fetch_dirs[fetch_ndirs-1]) /
                       (fetch_dirs[0] - fetch_dirs[fetch_ndirs-1]);
    }

    return fetch_s;
}

/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void read_sub_daily_met(int julian, int iclock, MetDataType *met)
{
    AED_REAL now ;
    int  idx = 0;

    now = julian;
    if ( now != loaded_day ) {
        fprintf(stderr, "Loaded day %12.6f Not equal to %12.4f\n", loaded_day, now);
        exit(1);
    }

    idx = iclock/3600;

    *met = submet[idx];
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
static void locate_time_column(int csv, const char *which, const char *fname)
{
    int lt_idx = -1;

    if ( (lt_idx = find_csv_var(csv, "time")) < 0 )
        lt_idx = find_csv_var(csv, "date");

    if (lt_idx != 0) {
        fprintf(stderr,
                 "Error in %s file '%s': 'Time (Date)' is not first column!\n",
                                                                  which, fname);
        exit(1);
    }
    if (time_idx < 0) time_idx = lt_idx;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void open_met_file(const char *fname, int snow_sw, int rain_sw,
                                                            const char *timefmt)
{
    char *LWTypeName;

    have_snow = snow_sw;
    have_rain_conc = rain_sw;
    have_fetch = fetch_sw;
    if ( (metf = open_csv_input(fname, timefmt)) < 0 ) {
        fprintf(stderr, "Failed to open '%s'\n", fname);
        exit(1);
    }

    locate_time_column(metf, "met", fname);

    if ( (rain_idx = find_csv_var(metf,"Rain")) < 0 ) {
        fprintf(stderr,"Error in met file, Rain not found!\n");
        exit(1);
    }
    if ( (hum_idx = find_csv_var(metf,"RelHum")) < 0 ) {
        fprintf(stderr,"Error in met file, RelHum not found!\n");
        exit(1);
    }
    if ( lw_ind == LW_CC )
        LWTypeName = "Cloud";
    else
        LWTypeName = "LongWave";

    if ( rad_mode != 3 && rad_mode != 5 ) {
        if ((lwav_idx = find_csv_var(metf, LWTypeName)) < 0 ) {
            fprintf(stderr,"Error in met file, '%s' not found!\n", LWTypeName);
            exit(1);
        }
    }
    sw_idx = find_csv_var(metf,"ShortWave");
    atmp_idx = find_csv_var(metf,"AirTemp");
    wind_idx = find_csv_var(metf,"WindSpeed");
    if ( have_snow ) {
        snow_idx = find_csv_var(metf,"Snow");
        if ( snow_idx < 0 ) {
            have_snow = FALSE;
            fprintf(stderr,"Warning in met file, snowice is enabled but Snow column not found!\n");
        }

    }
    if (have_rain_conc) {
        rpo4_idx = find_csv_var(metf,"rainPO4");
        rtp_idx  = find_csv_var(metf,"rainTP");
        rno3_idx = find_csv_var(metf,"rainNO3");
        rnh4_idx = find_csv_var(metf,"rainNH4");
        rtn_idx  = find_csv_var(metf,"rainTN");
        rsi_idx  = find_csv_var(metf,"rainSi");
    }

    wdir_idx = find_csv_var(metf,"WindDir");
    if (wdir_idx != -1 ) {
        if ( !have_fetch ) {
            fprintf(stderr, "Met file has wind direction - but fetch_sw is off\n");
        }
    } else {
        if ( have_fetch ) {
            fprintf(stderr, "fetch_sw is on but there is no wind direction in the met file\n");
            have_fetch = FALSE;
        }
    }

    n_steps = 86400.0 / timestep;
    // Allocate sub daily met array with an element for each timestep
    submet = calloc(n_steps, sizeof(MetDataType));

    if (subdaily) {
        if (rad_mode == 0) { //Then need to determine rad_mode from longwave type
            if ( sw_idx != -1 )  { // we have solar data
                if ( lwav_idx == -1 ) rad_mode = 3;
                else {
                    if ( lw_ind == LW_CC ) rad_mode = 1;
                    else                   rad_mode = 2;
                }
            } else { // no solar data
                if ( lwav_idx == -1 ) rad_mode = 5;
                else {
                    if ( lw_ind == LW_CC ) rad_mode = 4;
                //  else                   rad_mode = X;
                }
            }
        }
    }

}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void open_kw_file(const char *fname, const char *timefmt)
{
    if ( (kwf = open_csv_input(fname, timefmt)) < 0 ) {
        fprintf(stderr, "Failed to open '%s'\n", fname);
        exit(1);
    }
    locate_time_column(kwf, "Kd", fname);
    if ( (kw_idx = find_csv_var(kwf, "Kd")) < 0 ) {
        fprintf(stderr,"Error in Kd file, Kd not found!\n");
        exit(1);
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void read_daily_kw(int julian, AED_REAL *kwout)
{
    int csv;
    if ( (csv = kwf) > -1 ) {
        find_day(csv, time_idx, julian);
        *kwout = get_csv_val_r(csv, kw_idx);
    } else
        *kwout = Kw; //just use the Kw supplied in the file
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void open_inflow_file(int idx, const char *fname,
                             int nvars, const char *vars[], const char *timefmt)
{
    int j,k,l;

    if ( (inf[idx].inf = open_csv_input(fname, timefmt)) < 0 ) {
        fprintf(stderr, "Failed to open '%s'\n", fname);
        exit(1) ;
    }
    locate_time_column(inf[idx].inf, "inflow", fname);

    inf[idx].flow_idx = find_csv_var(inf[idx].inf,"flow");
    inf[idx].temp_idx = find_csv_var(inf[idx].inf,"temp");
    inf[idx].salt_idx = find_csv_var(inf[idx].inf,"salt");
    l = 0;
    for (j = 0; j < nvars; j++) {
        k = find_csv_var(inf[idx].inf, vars[j]);
        if (k == -1)
            fprintf(stderr, "No match for '%s' in file '%s'\n", vars[j], fname);
        else {
            if ( k != inf[idx].flow_idx && k != inf[idx].temp_idx &&
                                           k != inf[idx].salt_idx )
                inf[idx].in_vars[l++] = k;
        }
    }
    inf[idx].n_vars = l;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void open_outflow_file(int idx, const char *fname, const char *timefmt)
{
    if ( (outf[idx].outf = open_csv_input(fname, timefmt)) < 0 ) {
        fprintf(stderr, "Failed to open '%s'\n", fname);
        exit(1) ;
    }

    locate_time_column(outf[idx].outf, "outflow", fname);

    outf[idx].draw_idx = find_csv_var(outf[idx].outf,"flow");
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void open_withdrtemp_file(const char *fname, const char *timefmt)
{
    if ( (withdrTempf.withdrTempf = open_csv_input(fname, timefmt)) < 0 ) {
        fprintf(stderr, "Failed to open '%s'\n", fname);
        exit(1) ;
    }

    locate_time_column(withdrTempf.withdrTempf, "withdrTemp", fname);

    withdrTempf.wtemp_idx = find_csv_var(withdrTempf.withdrTempf,"temp");
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void close_kw_files()
{ if ( kwf >= 0 ) close_csv_input(kwf); }
/*----------------------------------------------------------------------------*/
void close_met_files()
{ if ( metf >= 0 ) close_csv_input(metf); }
/*----------------------------------------------------------------------------*/
void close_inflow_files()
{ int i; for (i = 0; i < NumInf; i++) close_csv_input(inf[i].inf); }
/*----------------------------------------------------------------------------*/
void close_outflow_files()
{ int i; for (i = 0; i < NumOut; i++) close_csv_input(outf[i].outf); }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void close_withdrtemp_files()
{ if ( withdrTempf.withdrTempf > -1 ) close_csv_input(withdrTempf.withdrTempf); }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
