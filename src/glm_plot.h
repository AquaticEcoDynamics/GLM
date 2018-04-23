/******************************************************************************
 *                                                                            *
 * glm_plot.h                                                                 *
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
#ifndef _GLM_PLOT_H_
#define _GLM_PLOT_H_

#ifdef PLOTS

#ifdef _FORTRAN_SOURCE_

  INTERFACE

     !##########################################################################
     SUBROUTINE init_plots(jstart,ndays,crest) BIND(C,name="init_plots_")
        USE ISO_C_BINDING
        CINTEGER,INTENT(in) :: jstart,ndays
        AED_REAL,INTENT(in) :: crest
     END SUBROUTINE init_plots

     SUBROUTINE put_glm_val_s(plot_id,val) BIND(C,name="put_glm_val_s_")
#       if defined( _WIN32 ) && USE_DL_LOADER
        !DEC$ ATTRIBUTES DLLIMPORT :: put_glm_val_s_
#       endif
        USE ISO_C_BINDING
        CINTEGER,INTENT(in) :: plot_id
        AED_REAL,INTENT(in) :: val(*)
     END SUBROUTINE put_glm_val_s

     SUBROUTINE put_glm_val(plot_id,val) BIND(C,name="put_glm_val_")
#       if defined( _WIN32 ) && USE_DL_LOADER
        !DEC$ ATTRIBUTES DLLIMPORT :: put_glm_val_
#       endif
        USE ISO_C_BINDING
        CINTEGER,INTENT(in) :: plot_id
        AED_REAL,INTENT(in) :: val(*)
     END SUBROUTINE put_glm_val
     !++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

  END INTERFACE

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#else
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#ifdef XPLOTS
    extern int xdisp;
#endif
    extern CLOGICAL do_plots, saveall;
    extern int today, plotstep;
    extern AED_REAL psubday;

/******************************************************************************/
void init_plots(int jstart, int ndays, AED_REAL crest);
void put_glm_val_s(int plot_id, AED_REAL *val);
void put_glm_val(int plot_id, AED_REAL *val);
void do_internal_plots(const int plot_id[]);

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void init_plots_(int *jstart, int *ndays, AED_REAL *crest);
void put_glm_val_s_(int *plot_id, AED_REAL *val);
void put_glm_val_(int *plot_id, AED_REAL *val);

#endif

#endif

#endif
