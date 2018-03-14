/******************************************************************************
 *                                                                            *
 * glm_main.c                                                                 *
 *                                                                            *
 * Developed by :                                                             *
 *     AquaticEcoDynamics (AED) Group                                         *
 *     School of Earth & Environment                                          *
 *     The University of Western Australia                                    *
 *                                                                            *
 *     http://aed.see.uwa.edu.au/                                             *
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
#include <stdlib.h>
#include <string.h>

#include "glm.h"

#include "glm_types.h"
#include "glm_globals.h"

#ifdef PLOTS
   #include <libplot.h>
   extern CLOGICAL do_plots, saveall;
   extern char *plots_nml_name;
#ifdef XPLOTS
   extern int xdisp;
#endif
#endif

char *all_plots_name = NULL;
static int debug_on = FALSE;
static int dbg_mix  = FALSE;

extern char glm_nml_file[];
extern void run_model(void);

/******************************************************************************/
int main(int argc, char *argv[])
{
    char *nmlfile = NULL;
    int all_ok = 1, show_options = 0;

#ifdef PLOTS
    saveall = 0;
#ifdef XPLOTS
    xdisp = 0;
#endif
#endif

    argv++; argc--;
    while (argc > 0) {
        if ( strcmp(*argv, "--nml") == 0) {
            argv++; argc--;
            nmlfile = *argv;
        }
        else if (strcmp(*argv, "--debug") == 0) {
            debug_on = TRUE;
        }
        else if (strcmp(*argv, "--dbgmix") == 0) {
            dbg_mix = TRUE;
        }
#ifdef PLOTS
#ifdef XPLOTS
        else if (strcmp(*argv, "--xdisp") == 0) {
            xdisp = 1;
            if ( argc > 1 && strncmp(argv[1], "--", 2) != 0 ) {
                argv++; argc--;
                plots_nml_name = *argv;
            }
        }
#endif
        else if (strcmp(*argv, "--saveall") == 0) {
            if ( saveall == 0) saveall = 1;
        }
        else if (strcmp(*argv, "--save-all-in-one") == 0) {
            saveall = 2;
            if ( argc > 1 && strncmp(argv[1], "--", 2) != 0 ) {
                argv++; argc--;
                all_plots_name = *argv;
            }
        }
#endif
        else {
            if (strcmp(*argv, "--help") != 0)
                printf("Unknown flag %s\n", *argv);
            show_options = 1;
            all_ok = 0;
        }
        argc--; argv++;
    }

#ifdef PLOTS
# ifdef XPLOTS
    do_plots = xdisp || saveall;
# else
    do_plots = saveall;
# endif
#endif

    printf("       ----------------------------------------------------\n");
    printf("       |  General Lake Model (GLM)   Version %-12s |\n", GLM_VERSION);
    printf("       ----------------------------------------------------\n");

#ifdef __GNUC__
    printf("glm built using gcc version %d.%d.%d\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#elif defined(_MSC_VER) && !defined(__INTEL_COMPILER)
    printf("glm built using MSC version %ld\n", _MSC_VER);
#endif

    if ( show_options ) {
       printf("--help  : show this blurb\n");
       printf("--nml <nmlfile> : get parameters from nmlfile\n");
#ifdef PLOTS
#ifdef XPLOTS
       printf("--xdisp : display temp/salt and selected others in x-window\n");
       printf("--xdisp <plotsfile> : like --xdisp, but use <plotsfile> instead of plots.nml\n");
#endif
       printf("--saveall : save plots to png files\n");
       printf("--save-all-in-one : save all plots to png file\n");
       printf("--save-all-in-one <destfile> : save all plots to png file <destfile>\n");
#endif
    }
    else if ( all_ok ) {
        if ( nmlfile != NULL ) strncpy(glm_nml_file, nmlfile, 256);
        run_model();

        printf("------------------------------------------------\n");
        printf("              Run Complete\n");
    }

    exit(0);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#if DEBUG
void crash_()
{
    int *x = (int*)1;
    *x = 0;
}
#endif

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#include <stdarg.h>
#include <aed_time.h>
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

void _glm_dbg_on()  { debug_on = 1; }
void _glm_dbg_off() { debug_on = 0; }

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
static FILE *mixcsv = NULL;
void _dbg_mix_str(const char *fmt, ...)
{
    va_list nap;

    if ( ! dbg_mix ) return;
    if ( mixcsv == NULL ) {
        mixcsv = fopen("glm_mixer.csv", "w");
        if ( mixcsv == NULL ) dbg_mix = 0;
        else fprintf(mixcsv,
         "time,step,where,loop,bottom,surface,epi_bot,meta_top,Energy_AvailableMix,Energy_RequiredMix,redg\n");
    }

    if ( mixcsv == NULL ) return;

    va_start(nap, fmt);
    vfprintf(mixcsv, fmt, nap);
    va_end(nap);
}
static char _ts[20];
void _dbg_time(int jday, int iclock)
{
    if ( ! dbg_mix ) return;
    write_time_string(_ts, jday, iclock + noSecs);
}
void _dbg_mixer(int d1, int d2, int d3,
       int bl, int sl, int ebl, int mtl,
       AED_REAL e1, AED_REAL e2, AED_REAL e3)
{
    if ( ! dbg_mix ) return;
    _dbg_mix_str("%s,%d,%d,%d,%d,%d,%d,%e,%e,%e\n",
       _ts,d1,d2,d3, bl,sl,ebl,mtl,e1,e2,e3);
}

void _mix_dbg_on()  { dbg_mix = 1; }
void _mix_dbg_off() { dbg_mix = 0; }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
