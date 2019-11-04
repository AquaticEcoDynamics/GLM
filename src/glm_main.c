/******************************************************************************
 *                                                                            *
 * glm_main.c                                                                 *
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

#include "glm_debug.h"

char *all_plots_name = NULL;

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
            _glm_dbg_on();
        }
        else if (strcmp(*argv, "--dbgmix") == 0) {
            _mix_dbg_on();
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
        else if (strcmp(*argv, "--quiet") == 0) {
            quiet = 5;
            if ( argc > 1 && strncmp(argv[1], "--", 2) != 0 ) {
                argv++; argc--;
                quiet = atoi(*argv);
            }
        }
        else {
            if (strcmp(*argv, "--help") != 0)
                fprintf(stderr, "Unknown flag %s\n", *argv);
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

    if (quiet < 10) {
        printf("     \n");
        printf("    -------------------------------------------------------\n");
        printf("    |  General Lake Model (GLM)   Version %-12s    |\n", GLM_VERSION);
        printf("    -------------------------------------------------------\n");
        printf("     \n");

#ifdef __GNUC__
        printf("     glm built using gcc version %d.%d.%d\n", __GNUC__, __GNUC_MINOR__, __GNUC_PATCHLEVEL__);
#elif defined(_MSC_VER) && !defined(__INTEL_COMPILER)
        printf("     glm built using MSC version %ld\n", _MSC_VER);
#endif
    }

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

       printf("--quiet   : less messages\n");
       printf("--quiet <level> : set quiet level (1-10)\n");
#endif
    }
    else if ( all_ok ) {
        if ( nmlfile != NULL ) strncpy(glm_nml_file, nmlfile, 256);

        run_model();

        if (quiet < 10) {
            printf("    Model Run Complete\n");
            printf("    -------------------------------------------------------\n\n");
        }
    }

    exit(0);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
