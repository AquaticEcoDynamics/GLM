/******************************************************************************
 *                                                                            *
 * glm_debug.c                                                                *
 *                                                                            *
 * Developed by :                                                             *
 *     AquaticEcoDynamics (AED) Group                                         *
 *     School of Agriculture and Environment                                  *
 *     The University of Western Australia                                    *
 *                                                                            *
 *     http://aquatic.science.uwa.edu.au/                                     *
 *                                                                            *
 * Copyright 2018 -  The University of Western Australia                      *
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
#include <stdarg.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "glm.h"
#include "glm_globals.h"

#include "glm_debug.h"

#include <aed_time.h>

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


static int debug_on = FALSE;
static int ldbg_mix  = FALSE;

#if DEBUG
/******************************************************************************/
void crash_()
{
    int *x = (int*)1;
    *x = 0;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
static FILE *dbgfile = NULL;
void _glm_dbg(const char *fmt, ...)
{
    va_list nap;

    if ( ! debug_on ) return;
    if ( dbgfile == NULL ) dbgfile = fopen("glm_debug.log", "w");

    if ( dbgfile == NULL ) dbgfile = stderr;

    va_start(nap, fmt);
    vfprintf(dbgfile, fmt, nap);
    va_end(nap);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void LakeCheck(const char*str)
{
    if ( isnan(Lake[surfLayer].Height)   ||
         isnan(Lake[surfLayer].Density)  ||
         isnan(Lake[surfLayer].Salinity) ||
         isnan(Lake[surfLayer].Temp)     ||
         isnan(Lake[surfLayer].LayerVol) ||
         Lake[surfLayer].Salinity < 0 || Lake[surfLayer].Salinity > 200 ||
           Lake[surfLayer].Temp < -10 || Lake[surfLayer].Temp > 80 ) {
       fflush(stdout);
       fprintf(stderr, "\n\n LakeCheck temp/salt failed %s\n", str);
       fprintf(stderr, "  Lake[surfLayer].Height   = %f\n", Lake[surfLayer].Height);
       fprintf(stderr, "  Lake[surfLayer].Density  = %f\n", Lake[surfLayer].Density);
       fprintf(stderr, "  Lake[surfLayer].Salinity = %f\n", Lake[surfLayer].Salinity);
       fprintf(stderr, "  Lake[surfLayer].Temp     = %f\n", Lake[surfLayer].Temp);
       fprintf(stderr, "  Lake[surfLayer].LayerVol = %f\n", Lake[surfLayer].LayerVol);

       exit(1);
    }
//  dbgprt("%s temp = %f vol = %f dens = %f hgt = %f\n",str,
//        Lake[surfLayer].Temp, Lake[surfLayer].LayerVol, Lake[surfLayer].Density, Lake[surfLayer].Height);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#endif


/******************************************************************************/
void _glm_dbg_on()  { debug_on = 1; }
void _glm_dbg_off() { debug_on = 0; }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/

static FILE *mixcsv = NULL;
static char *More_Fields = NULL;
static int fields_inited = 0;


/******************************************************************************/
void _dbg_mix_init_fields()
{
    if (! dbg_mix) return;
    if (fields_inited) return;
    _dbg_mix_add_field(",Epi_dz");
    _dbg_mix_add_field(",MeanSalt");
    _dbg_mix_add_field(",MeanTemp");
    _dbg_mix_add_field(",gPrimeTwoLayer");
    _dbg_mix_add_field(",Vol_Epi");
    _dbg_mix_add_field(",Mass_Epi");
    _dbg_mix_add_field(",Time_end_shear");
    _dbg_mix_add_field(",Time_start_shear");
    _dbg_mix_add_field(",Half_Seiche_Period");
    _dbg_mix_add_field(",Thermocline_Height");
    _dbg_mix_add_field(",u_avg");
    _dbg_mix_add_field(",WindSpeedX");
    _dbg_mix_add_field(",Dens_Epil");
    _dbg_mix_add_field(",Epi_Thick");
    _dbg_mix_add_field(",dMdz");
    _dbg_mix_add_field(",q_cub");
    _dbg_mix_add_field(",LengthAtThermo");
    _dbg_mix_add_field(",Hypl_Thick");
    _dbg_mix_add_field(",Dens_Hypl");
    _dbg_mix_add_field(",Energy_Conv");
    _dbg_mix_add_field(",Energy_WindStir");
    _dbg_mix_add_field(",Energy_TotStir");
    _dbg_mix_add_field(",Energy_Deepen");
    _dbg_mix_add_field(",Energy_Shear");
    _dbg_mix_add_field(",del_u");
    _dbg_mix_add_field(",u_avgSQ");
    _dbg_mix_add_field(",u_eff");
    _dbg_mix_add_field(",u0_old");
    _dbg_mix_add_field(",u_avg_old");
    _dbg_mix_add_field(",deltaKH");
    _dbg_mix_add_field(",del_deltaKH");
    _dbg_mix_add_field(",GPEFFC");
    _dbg_mix_add_field(",accn");
    _dbg_mix_add_field(",zsml_tilda");
    _dbg_mix_add_field(",Slope");
    _dbg_mix_add_field(",VMsum");
    _dbg_mix_add_field(",Tsum");
    _dbg_mix_add_field(",Ssum");
    _dbg_mix_add_field(",DepMX");
    _dbg_mix_add_field(",PrevThick");
    _dbg_mix_add_field(",OldSlope");
    _dbg_mix_add_field(",Time_count_end_shear");
    _dbg_mix_add_field(",Time_count_sim");
    _dbg_mix_add_field(",FO");
    _dbg_mix_add_field(",FSUM");
    _dbg_mix_add_field(",u_f");
    _dbg_mix_add_field(",u0");
    _dbg_mix_add_field(",coef_mix_KH");
    _dbg_mix_add_field(",coef_mix_conv");
    _dbg_mix_add_field(",coef_wind_stir");
    _dbg_mix_add_field(",coef_mix_shear");
    _dbg_mix_add_field(",coef_mix_turb");
    _dbg_mix_add_field(",U_star");
    _dbg_mix_add_field(",U_star_sqr");
    _dbg_mix_add_field(",U_star_cub");
    _dbg_mix_add_field(",Epilimnion_Mid_Ht");
    _dbg_mix_add_field(",q_sqr");
    _dbg_mix_add_field(",IntWaveSpeed");
    _dbg_mix_add_field(",ZeroMom");
    _dbg_mix_add_field(",FirstMom");
    _dbg_mix_add_field(",delzkm1");
    _dbg_mix_add_field(",Vol_Hypl");
    _dbg_mix_add_field(",Hypl_Mass");
    fields_inited = 1;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
void _dbg_mix_str(const char *fmt, ...)
{
    va_list nap;

    if ( ! ldbg_mix ) return;
    if ( mixcsv == NULL ) {
        mixcsv = fopen("glm_mixer.csv", "w");
        if ( mixcsv == NULL ) ldbg_mix = 0;
        else fprintf(mixcsv,
         "time,step,where,loop,num_layers,epi_bot,meta_top,SurfTemp,Energy_AvailableMix,Energy_RequiredMix,redg%s\n",
            (More_Fields==NULL)?"":More_Fields);
    }

    if ( mixcsv == NULL ) return;

    va_start(nap, fmt);
    vfprintf(mixcsv, fmt, nap);
    va_end(nap);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
void _dbg_mix_add_field(const char *f)
{   int olen = (More_Fields==NULL)?0:strlen(More_Fields);
    char *tstr = realloc(More_Fields, olen+strlen(f)+2);
    tstr[olen] = 0;
    strcat(tstr, f);
    More_Fields = tstr;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
static char _ts[20];
void _dbg_time(int jday, int iclock)
{
    if ( ! ldbg_mix ) return;
    write_time_string(_ts, jday, iclock + noSecs);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
void _dbg_mixer_s(int d1, int d2, int ebl, int mtl)
{
    static int loop_count;
    if (d2 != _DBG_LOOP_) loop_count = 0;
    if ( ! ldbg_mix ) return;
    _dbg_mix_str("%s,%d,%d,%d,%d,%d,%d,%e",
       _ts, d1, d2, loop_count++, NumLayers, ebl, mtl, Lake[surfLayer].Temp);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
void _dbg_mixer_a(AED_REAL e1) { _dbg_mix_str(",%e", e1); }
void _dbg_mixer_e() { _dbg_mix_str("\n"); }

void _mix_dbg_on()  { ldbg_mix = dbg_mix; } // turn on only if global flag is set
void _mix_dbg_off() { ldbg_mix = 0; }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


#if DEBUG
/******************************************************************************/
static FILE *lakefl = NULL;
static void _lak_dbg(const char *fmt, ...)
{
    va_list nap;

    if ( lakefl == NULL )
        lakefl = fopen("glm_lake_dbg.csv", "w");

    if ( lakefl == NULL ) return;

    va_start(nap, fmt);
    vfprintf(lakefl, fmt, nap);
    va_end(nap);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
void _DumpLake(int where, int extra)
{
    int i;
    _lak_dbg("%s,%c,%d,%d,%3d", _ts,'H',where,extra,NumLayers);
    for (i = surfLayer; i >= botmLayer; i--) {
        _lak_dbg(",%.*e", 16, Lake[i].Height);
    }
    _lak_dbg("\n");
    _lak_dbg("%s,%c,%d,%d,%3d", _ts,'T',where,extra,NumLayers);
    for (i = surfLayer; i >= botmLayer; i--) {
        _lak_dbg(",%.*e", 16, Lake[i].Temp);
    }
    _lak_dbg("\n");
    _lak_dbg("%s,%c,%d,%d,%3d", _ts,'D',where,extra,NumLayers);
    for (i = surfLayer; i >= botmLayer; i--) {
        _lak_dbg(",%.*e", 16, Lake[i].Density);
    }
    _lak_dbg("\n");
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#endif
