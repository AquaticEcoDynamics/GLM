/******************************************************************************
 *                                                                            *
 * glm_mobl.c                                                                 *
 *                                                                            *
 * code for mobility                                                          *
 *                                                                            *
 * Assumptions :                                                              *
 *                                                                            *
 *  1) movement direction has at most one change down the layers              *
 *  2) sides of the lake slope inward (ie bottom is narrower than top)        *
 *                                                                            *
 * -------------------------------------------------------------------------- *
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
#include <math.h>
#include "glm.h"
#include "glm_types.h"
#include "glm_mobl.h"



/******************************************************************************/
static AED_REAL Rising(AED_REAL *Y, AED_REAL *cc, AED_REAL *ww, AED_REAL *vols,
                              AED_REAL *A,  AED_REAL *mins, AED_REAL dt,
                                                           int start, int end);
static AED_REAL Sinking(AED_REAL *Y, AED_REAL *cc, AED_REAL *ww, AED_REAL *vols,
                              AED_REAL *A,  AED_REAL *mins, AED_REAL dt,
                                                           int start, int end);




/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void Mobility(int *N_in,          // number of vertical layers
              AED_REAL *dt_in,    // time step (s)
              AED_REAL *h,        // array of layer thicknesses (m)
              AED_REAL *A,        // array of layer areas (m^2)
              AED_REAL *ww,       // array of vertical mobility speeds (m/s)
              AED_REAL *min_C_in, // minimum concentration (mmol/m^3)
              AED_REAL *cc)       // array of cell concentrations (mmol/m^3)
{
    /**************************************************************************
     * Since this routine is called only from the fortran the arguments are   *
     * in the "by reference" form used by fortran, and assigned to local      *
     * dereferenced variables                                                 *
     **************************************************************************/
    int      N = *N_in;
    AED_REAL dt = *dt_in;
    AED_REAL min_C = *min_C_in;

    int      i;

    AED_REAL *vols;  // layer volume (m^3)
    AED_REAL *mins;  // minimum layer mass of variable (mmol)
    AED_REAL *Y;     // total mass of variable (mmol) in layer

    AED_REAL dtMax = dt, tdt, tmp;
    int dirChng, sign;

/*----------------------------------------------------------------------------*/
    vols = calloc(N, sizeof(AED_REAL));
    mins = calloc(N, sizeof(AED_REAL));
    Y = calloc(N, sizeof(AED_REAL));

    /**************************************************************************
     * determine mobility timestep i.e. maximum time step that particles      *
     * will not pass through more than one layer in a time step               *
     **************************************************************************/

    dtMax = dt;
    dirChng = 0;   //this represents the layer at which direction switches from sinking to rising or visa versa
    sign = ( ww[0] > 0. );  // positive for rising, negative for sinking

    for (i = 0; i < N; i++) {
        // for convenience
        vols[i] = (h[i] * A[i]);
        mins[i] = (min_C * vols[i]);
        Y[i] = cc[i] * vols[i];

        // look for the change of direction
        if ( sign != ( ww[i] > 0. ) ) { sign = !sign; dirChng = i-1; }

        /**********************************************************************
         * check if all movement can be from within this cell                 *
         **********************************************************************/
        if ( (fabs(ww[i]) * dt) > h[i] ) {
            tdt = h[i]/fabs(ww[i]);
            if ( tdt < dtMax ) dtMax = tdt;
        }

        /**********************************************************************
         * check if movement can all be into the next cell.                   *
         * if movement is settling, next is below, otherwise next is above.   *
         * check also for top or bottom in case of oopsies                    *
         **********************************************************************/
        if ( ww[i] > 0. ) {
            if ( (i < N-1) && (fabs(ww[i]) * dt) > h[i+1] ) {
                tdt = h[i+1]/fabs(ww[i]);
                if ( tdt < dtMax ) dtMax = tdt;
            }
        } else {
            if ( (i > 0) && (fabs(ww[i]) * dt) > h[i-1] ) {
                tdt = h[i-1]/fabs(ww[i]);
                if ( tdt < dtMax ) dtMax = tdt;
            }
        }
    } // end find maximum time step dtMax
    if ( dirChng == 0 && ww[0] > 0. ) dirChng = N; //all rising
    if ( dirChng == 0 && ww[N-1] < 0. ) dirChng = N; //all sinking

    /**************************************************************************/

    tdt = dtMax;
    do  {   // do this in steps of dtMax, but at least once
        // each time tdt is dtMax, except, possibly, the last which is whatever
        // was left.
        if ( (dt -= dtMax) < 0. ) tdt = dtMax + dt;

        /**********************************************************************
         * 2 possibilities,                                                   *
         *   1) lower levels rising, upper levels sinking                     *
         *   2) lower levels sinking, upper levels rising                     *
         **********************************************************************/
        if ( ww[0] > 0. ) { // lower levels rising
            if ( ww[N-1] < 0. ) { // top levels are sinking
                tmp = Sinking(Y, cc, ww, vols, A, mins, tdt, N-1, dirChng);
                Y[dirChng-1] += tmp;
                cc[dirChng-1] = Y[dirChng-1] / vols[dirChng-1];
            }
            Rising( Y, cc, ww, vols, A, mins, tdt, 0,   dirChng-1);
        } else { // lower levels sinking
            Sinking(Y, cc, ww, vols, A, mins, tdt, dirChng-1, 0);
            if ( ww[N-1] > 0. )   // top levels are rising
                Rising( Y, cc, ww, vols, A, mins, tdt, dirChng,   N-1);
        }
    } while ( dt > 0. );

    /**************************************************************************/
    free(vols); free(mins); free(Y);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 * Rising is the easier on the two since the slope means we dont need to look *
 * at relative areas (the cell above will always be >= to the current cell)   *
 * all matter is moved to the next cell                                       *
 *  for each cell except the top :                                            *
 *    1) calculate how much is going to move                                  *
 *    2) subtract amount that must now move                                   *
 *    3) add the amount moved from previous cell to current cell              *
 *    4) fix concentration                                                    *
 *  for the top cell :                                                        *
 *    1) add the amount moved from previous cell to current cell              *
 *    2) fix concentration                                                    *
 ******************************************************************************/
AED_REAL Rising(AED_REAL *Y, AED_REAL *cc, AED_REAL *ww, AED_REAL *vols,
                             AED_REAL *A,  AED_REAL *mins, AED_REAL dt,
                                                           int start, int end)
{
    int i;
    AED_REAL mov = 0., moved = 0.;

    // go from start to one below end
    for (i = start; i < end; i++) {
        // speed times time (=h) time area * concen = mass to move
        mov = (ww[i] * dt) * A[i] * cc[i];
        //  if removing that much would bring it below min conc
        if ( (Y[i] + moved - mov) < mins[i] ) mov = Y[i] + moved - mins[i];

        Y[i] = Y[i] + moved - mov;
        cc[i] = Y[i] / vols[i]; // return it to a concentration
        moved = mov; // for the next step
    }
    // nothing rises out of the end cell, but we still add that which came from below
    Y[end] += moved;
    cc[end] = Y[end] / vols[end];

    return 0.; // nothing has moved out of this cell
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *  for each cell :                                                           *
 *    1) calculate how much is going to move                                  *
 *    2) subtract amount that must now move                                   *
 *    3) add the amount moved from previous cell to current cell              *
 *    4) fix concentration                                                    *
 *    5) compute the amount that will go to the next cell                     *
 ******************************************************************************/
AED_REAL Sinking(AED_REAL *Y, AED_REAL *cc, AED_REAL *ww, AED_REAL *vols,
                              AED_REAL *A,  AED_REAL *mins, AED_REAL dt,
                                                           int start, int end)
{
    int i;
    AED_REAL mov = 0., moved = 0.;

    for (i = start; i > end; i--) {
        // speed times time (=h) time area * concen = mass to move
        mov = (fabs(ww[i]) * dt) * A[i] * cc[i];
        //  if removing that much would bring it below min conc
        if ( (Y[i] + moved - mov) < mins[i] ) mov = Y[i] + moved - mins[i];

        Y[i] = Y[i] + moved - mov;
        cc[i] = Y[i] / vols[i]; // return it to a concentration

        // so now mov has how much has moved out of the cell, but not all
        // of that will go into the next cell
        if ( i > 0 ) moved = mov * (A[i-1] / A[i]); // for the next step
        else moved = mov; // we are about to exit anyway.
    }

    return moved; // this is how much left this cell
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
