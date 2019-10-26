/******************************************************************************
 *                                                                            *
 * glm_balance.c                                                              *
 *                                                                            *
 *  mass balance file for glm                                                 *
 *                                                                            *
 * Developed by :                                                             *
 *     AquaticEcoDynamics (AED) Group                                         *
 *     School of Agriculture and Environment                                  *
 *     The University of Western Australia                                    *
 *                                                                            *
 *     http://aquatic.science.uwa.edu.au/                                     *
 *                                                                            *
 * Copyright 2019 -  The University of Western Australia                      *
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
#include "glm_balance.h"

static int mbf = -1;
static int mbnv = 0;

//------------------------------------------------------------------------------
// These for mass balance
//------------------------------------------------------------------------------
AED_REAL *mb_ifvar = NULL;
AED_REAL *mb_ofvar = NULL;
AED_REAL *mb_lkvar = NULL;
int      *mb_idx = NULL;


/******************************************************************************
 ******************************************************************************/
void mb_add_inflows(AED_REAL vol, AED_REAL *wq_vars)
{
    int i;
    if ( mbf < 0 ) return;

    for (i = 0; i < mbnv; i++) {
        mb_ifvar[i] += (wq_vars[mb_idx[i]] * vol);
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 ******************************************************************************/
void mb_sub_outflows(int layer, AED_REAL subvol)
{
    int i;
    if ( mbf < 0 ) return;

    for (i = 0; i < mbnv; i++) {
        mb_ifvar[i] += (_WQ_Vars(mb_idx[i], layer) * subvol);
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 * open the mass balance debug file                                           *
 ******************************************************************************/
void open_balance(const char *out_dir, const char *balance_fname,
              int balance_varnum, const char**balance_vars, const char *timefmt)
{
    int i;
    size_t l;
    VARNAME mbs;

    if ( (mbf = open_csv_output(out_dir, balance_fname)) < 0 ) {
        fprintf(stderr, "Failed to create '%s'\n", balance_fname);
        exit(1) ;
    }

    mb_ifvar = calloc(balance_varnum, sizeof(AED_REAL));
    mb_ofvar = calloc(balance_varnum, sizeof(AED_REAL));
    mb_lkvar = calloc(balance_varnum, sizeof(AED_REAL));
    mb_idx   = calloc(balance_varnum, sizeof(int));

    csv_header_start(mbf);
    for (i = 0; i < balance_varnum; i++) {
        sprintf(mbs, "Inf %s", balance_vars[i]);
        csv_header_var(mbf, mbs);
        sprintf(mbs, "Lake %s", balance_vars[i]);
        csv_header_var(mbf, mbs);
        sprintf(mbs, "Out %s", balance_vars[i]);
        csv_header_var(mbf, mbs);

        mb_ifvar[i] = 0.0;
        mb_ofvar[i] = 0.0;
        mb_lkvar[i] = 0.0;

        l = strlen(balance_vars[i]);
        mb_idx[i] = wq_var_index_c(balance_vars[i], &l);
    }
    csv_header_end(mbf);

    mbnv = balance_varnum;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 * Write the mass balance file.                                               *
 ******************************************************************************/
void write_balance(int jday)
{
    char ts[64];
    int i, j;

    if ( mbf < 0 ) return;

    write_time_string(ts, jday, 0);

    write_csv_start(mbf, ts);

    for (i = 0; i < mbnv; i++) {
        write_csv_val(mbf, mb_ifvar[i]);
        write_csv_val(mbf, mb_ofvar[i]);

        for (j = 0; j < surfLayer; j++) {
            mb_lkvar[i] +=  (_WQ_Vars(mb_idx[i], j) * Lake[j].LayerVol);
        }
        write_csv_val(mbf, mb_lkvar[i]);

        mb_ifvar[i] = 0.0;
        mb_ofvar[i] = 0.0;
        mb_lkvar[i] = 0.0;
    }

    write_csv_end(mbf);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 * Close files and plots                                                      *
 ******************************************************************************/
void close_balance()
{
    if ( mbf >= 0 ) close_csv_output(mbf);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
