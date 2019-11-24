/******************************************************************************
 *                                                                            *
 * glm_plot.c                                                                 *
 *                                                                            *
 * plotting for glm                                                           *
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
#ifdef PLOTS

#include <stdio.h>
#include <string.h>
#include <stdlib.h>

#include "glm.h"

#include "glm_types.h"
#include "glm_globals.h"
#include "glm_plot.h"

#include "glm_output.h"
#include "glm_wqual.h"

#include "namelist.h"
#include "libplot.h"

#ifdef XPLOTS
    int xdisp = 0;
#endif

CLOGICAL do_plots, saveall = 0;
static int nplots = 9;
static int max_plots = 16, *theplots = NULL;
static char **vars = NULL;
int today = 0;
int plotstep = 0;
AED_REAL psubday = 1;
char * plots_nml_name = "plots.nml";


/******************************************************************************/
void init_plots(int jstart, int ndays, AED_REAL crest)
{
    int        maxx, maxy, width, height,i,w,h,acrs;
    int        plot_width, plot_height;
    int        namlst;
    char     **title;
    AED_REAL  *min_z;
    AED_REAL  *max_z;
    AED_REAL   min_x, max_x;
    AED_REAL   min_y = -1, max_y = -1;
    char      *glm_vers = NULL;
    char      *title_font = NULL;
    char      *label_font = NULL;
    int        tsz, lsz;
    int        default_t = FALSE;

    NAMELIST plots_window[] = {
          { "plots_window",   TYPE_START,            NULL               },
          { "width",          TYPE_INT,              &width             },
          { "height",         TYPE_INT,              &height            },
          { NULL,             TYPE_END,              NULL               }
    };
    NAMELIST plots[] = {
          { "plots",          TYPE_START,            NULL               },
          { "nplots",         TYPE_INT,              &nplots            },
          { "plot_width",     TYPE_INT,              &plot_width        },
          { "plot_height",    TYPE_INT,              &plot_height       },
          { "title",          TYPE_STR|MASK_LIST,    &title             },
          { "title_font",     TYPE_STR,              &title_font        },
          { "title_size",     TYPE_INT,              &tsz               },
          { "label_font",     TYPE_STR,              &label_font        },
          { "label_size",     TYPE_INT,              &lsz               },
          { "vars",           TYPE_STR|MASK_LIST,    &vars              },
          { "min_z",          TYPE_DOUBLE|MASK_LIST, &min_z             },
          { "max_z",          TYPE_DOUBLE|MASK_LIST, &max_z             },
          { "min_y",          TYPE_DOUBLE|MASK_LIST, &min_y             },
          { "max_y",          TYPE_DOUBLE|MASK_LIST, &max_y             },
          { NULL,             TYPE_END,              NULL               }
    };

/*----------------------------------------------------------------------------*/

    if ( (namlst = open_namelist(plots_nml_name)) < 0 ) {
        fprintf(stderr, "Error opening namelist file plots.nml\n");
        fprintf(stderr, "Using defaults\n");

        width = 1000; height=  300;

        nplots = 2; plot_width = 400; plot_height = 200;

        // title = { "Temperature", "Salinity" };
        title = malloc(3*sizeof(char*));
        title[0] = "Temperature"; title[1] = "Salinity"; title[2] = NULL;
        // vars  = { "temp","salt" };
        vars = malloc(3*sizeof(char*));
        vars[0] = "temp"; vars[1] = "salt"; vars[2] = NULL;
        // min_z = { -1.0, 0.0 };
        min_z = malloc(3 * sizeof(double));
        min_z[0] = -1.0; min_z[1] = 0.0; min_z[2] = 0.0;
        // max_z = { 30.0, 0.6 };
        max_z = malloc(3 * sizeof(double));
        max_z[0] = 30.0; max_z[1] = 0.6; max_z[2] = 0.0;

        default_t = TRUE;
    } else {
        /*-----------------------------------------------*/
        if ( get_namelist(namlst, plots_window) != 0 ) {
            fprintf(stderr, "Error reading the 'plots_window' namelist from plots.nml\n");
            fprintf(stderr, "Using defaults\n");

            width = 1000; height=  300;
        }

        if ( get_namelist(namlst, plots) != 0 ) {
            fprintf(stderr,"Error reading the 'plots' namelist from plots.nml\n");
            fprintf(stderr, "Using defaults\n");

            nplots = 2; plot_width = 400; plot_height = 200;

            // title = { "Temperature", "Salinity" };
            title = malloc(3*sizeof(char*));
            title[0] = "Temperature"; title[1] = "Salinity"; title[2] = NULL;
            // vars  = { "temp","salt" };
            vars = malloc(3*sizeof(char*));
            vars[0] = "temp"; vars[1] = "salt"; vars[2] = NULL;
            // min_z = { -1.0, 0.0 };
            min_z = malloc(3 * sizeof(double));
            min_z[0] = -1.0; min_z[1] = 0.0; min_z[2] = 0.0;
            // max_z = { 30.0, 0.6 };
            max_z = malloc(3 * sizeof(double));
            max_z[0] = 30.0; max_z[1] = 0.6; max_z[2] = 0.0;

            default_t = TRUE;
        }
    }

    max_plots = nplots + 12;
    theplots = malloc(max_plots*sizeof(int));

    for (i = 0; i < max_plots; i++) theplots[i] = -1;

    glm_vers = malloc(strlen(GLM_VERSION) + 10);
    sprintf(glm_vers, "GLM-%s", GLM_VERSION);

    maxx = width;
    maxy = height;

    set_progname(glm_vers);
    set_shortprogname("GLM");
    set_aboutmessage("General Lake Model\n(C) The University of Western Australia\nhttp://aquatic.science.uwa.edu.au/");
#ifdef PLOTS
#ifdef XPLOTS
    if ( xdisp ) {
        if ( init_plotter_max(max_plots, &maxx, &maxy) < 0 ) exit(1);
    } else
#endif
        if (do_plots) init_plotter_no_gui();
#endif

    acrs = (maxx + 90) / (100 + plot_width);

    if (title_font != NULL) set_plot_font(PF_TITLE, tsz, title_font);
    if (label_font != NULL) set_plot_font(PF_LABEL, lsz, label_font);

    min_x = jstart;
    max_x = jstart + ndays + (1-psubday);
    if ( min_y == -1 ) min_y = 0;
    if ( max_y == -1 ) max_y = crest;
    w = 10;
    h = 10;
    for (i = 0; i < nplots; i++) {
        int vn = 0;
        if ( !(vn = intern_is_var(i, vars[i])) ) {
            if (WQ_Vars != NULL) {
                size_t l = strlen(vars[i]);
                if ( ! (vn = wq_is_var(&i, vars[i], &l)) ) {
                    fprintf(stderr, "No plottable var \"%s\"\n", vars[i]);
                    continue;
                }
//              else if ( vn < 0 ) fprintf(stderr, "WQ sheet var \"%s\"\n", vars[i]);
//              else  fprintf(stderr, "WQ var \"%s\"\n", vars[i]);
            }
        }
//      else if ( vn < 0 ) fprintf(stderr, "Internal sheet var \"%s\"\n", vars[i]);
//      else  fprintf(stderr, "Internal var \"%s\"\n", vars[i]);

        theplots[i] = create_plot(w, h, plot_width, plot_height, title[i]);
        w += plot_width + 100;
        if ( ((i+1) % acrs) == 0 ) {
            w = 10;
            h += plot_height + 100;
        }
        set_plot_x_label(theplots[i], "Time");
        set_plot_x_limits(theplots[i], min_x, max_x);
        if ( vn > 0 ) {
            set_plot_y_label(theplots[i], "Depth");
            set_plot_y_limits(theplots[i], min_y, max_y);
            set_plot_z_limits(theplots[i], min_z[i], max_z[i]);
        } else {
            set_plot_y_label(theplots[i], "mmol/m3");
            set_plot_y_limits(theplots[i], min_z[i], max_z[i]);
        }
        set_plot_version(theplots[i], glm_vers);
        set_plot_varname(theplots[i], vars[i]);
        if ( n_zones > 0 ) {
            int j;
            for (j = 0; j < n_zones; j++) show_h_line(i, theZones[j].zheight);
        } else {
            show_h_line(i, CrestHeight);
        }
    }
    free(glm_vers);
    if (default_t) {
        free(title); free(vars); free(min_z); free(max_z);
        title=NULL; vars = NULL; min_z = NULL; max_z = NULL;
    }
    if (namlst >= 0) close_namelist(namlst);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
void do_internal_plots(const int plot_id[])
{
    AED_REAL todayish;
    int i, j;

    todayish = psubday;
    todayish *= plotstep;
    todayish += today;
    for (i = 0; i < NumLayers; i++) {
        if ( (j=plot_id[0]) >= 0 ) plot_value(theplots[j], todayish, Lake[i].Height, Lake[i].Temp);
        if ( (j=plot_id[1]) >= 0 ) plot_value(theplots[j], todayish, Lake[i].Height, Lake[i].Salinity);
        if ( (j=plot_id[2]) >= 0 ) plot_value(theplots[j], todayish, Lake[i].Height, Lake[i].Light);
        if ( (j=plot_id[3]) >= 0 ) plot_value(theplots[j], todayish, Lake[i].Height, Lake[i].ExtcCoefSW);
        if ( (j=plot_id[4]) >= 0 ) plot_value(theplots[j], todayish, Lake[i].Height, Lake[i].Density);
        if ( (j=plot_id[5]) >= 0 ) plot_value(theplots[j], todayish, Lake[i].Height, Lake[i].Uorb);
        if ( (j=plot_id[6]) >= 0 ) plot_value(theplots[j], todayish, Lake[i].Height, Lake[i].LayerStress);
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
void put_glm_val_s_(int *plot_id, AED_REAL *val)
{ put_glm_val_s(*plot_id, val); }
/******************************************************************************/
void put_glm_val_s(int plot_id, AED_REAL *val)
{
    AED_REAL todayish;

/*----------------------------------------------------------------------------*/
    if ( !do_plots || plot_id >= nplots || today <= 0 ) return;

    todayish = psubday;
    todayish *= plotstep;
    todayish += today;

    plot_value(theplots[plot_id], todayish, val[0], 0.);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
void put_glm_val_(int *plot_id, AED_REAL *val)
{ put_glm_val(*plot_id, val); }
/******************************************************************************/
void put_glm_val(int plot_id, AED_REAL *val)
{
    int i;
    AED_REAL todayish;
    AED_REAL height_for_plot, val_for_plot;

/*----------------------------------------------------------------------------*/
    if ( !do_plots || plot_id >= nplots || today <= 0 ) return;

    todayish = psubday;
    todayish *= plotstep;
    todayish += today;

    for (i = 0; i < NumLayers; i++) {
        height_for_plot = Lake[i].Height;
        val_for_plot = val[i];
        if (Lake[surfLayer].Height < 0.01005) {
            height_for_plot = 0.0;
            val_for_plot = 0.0;
        }
        plot_value(theplots[plot_id], todayish, height_for_plot, val_for_plot);
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#endif
