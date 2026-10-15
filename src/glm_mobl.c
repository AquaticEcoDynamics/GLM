/******************************************************************************
 *                                                                            *
 * glm_mobl.c                                                                 *
 *                                                                            *
 * code for mobility                                                          *
 *                                                                            *
 * Assumptions :                                                              *
 *                                                                            *
 *  1) sides of the lake slope inward (ie bottom is narrower than top)        *
 *  2) layer 0 is the bottom layer, layer N-1 the surface layer               *
 *                                                                            *
 * Deposition is an EXPORT: material that settles onto the lake bed (the bed  *
 * ring of a layer, or the flat bed under the bottom layer) leaves the water  *
 * column and is not credited to any sediment pool or zone here. AED's        *
 * benthic sedimentation terms (e.g. Psed_poc / swi_poc in                    *
 * aed_organic_matter) are computed separately as w*C of the bottom water and *
 * are reporting estimates, not a flux balanced against this removal.         *
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
 * Copyright 2013-2026 : The University of Western Australia                  *
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



/******************************************************************************
 *                                                                            *
 * Vertical settling/rising of one variable through the layers.               *
 *                                                                            *
 * Each layer moves at its own velocity ww[i] (+ve up, -ve down), so any      *
 * pattern of rising and sinking layers is handled. In each sub-step, from    *
 * the concentrations at the start of the sub-step:                           *
 *   sinking layer i : mov = |w| dt A[i] C[i] leaves the layer; the share     *
 *                     A[i-1]/A[i] enters layer i-1 (the part over the        *
 *                     layer's bottom face) and the rest deposits on the bed  *
 *                     ring between A[i-1] and A[i]. For the bottom layer     *
 *                     all of it deposits on the flat bed.                    *
 *   rising layer i  : mov = w dt A[i] C[i] enters layer i+1; nothing rises   *
 *                     out of the surface layer.                              *
 * A[i] is the area at the top of layer i and vols are the TRUE layer         *
 * volumes, so the mass moved is exactly the change in C * V.                 *
 *                                                                            *
 * (The previous version handled at most one change of direction and had two *
 * index errors: with all layers sinking the bottom layer never settled and   *
 * never received the material falling from the layer above (it was          *
 * discarded), and where the direction changed, layers were assigned to the  *
 * wrong direction - sinking layers were moved upward, giving negative        *
 * concentrations - or material was delivered two layers down.)               *
 *                                                                            *
 ******************************************************************************/
void doMobility(const CINTEGER *N_in,   // number of vertical layers
              const AED_REAL *dt_in,    // time step (s)
              const AED_REAL *h,        // array of layer thicknesses (m)
              const AED_REAL *A,        // array of layer areas (m^2)
              const AED_REAL *V,        // array of layer volumes (m^3)
              const AED_REAL *ww,       // array of vertical mobility speeds (m/s)
              const AED_REAL *min_C_in, // minimum concentration (mmol/m^3)
              AED_REAL *mcc)            // array of cell concentrations (mmol/m^3)
{
    doMobilityDep(N_in, dt_in, h, A, V, ww, min_C_in, mcc, NULL);
}

/******************************************************************************
 * As doMobility, and if dep is not NULL also adds to dep[i] the mass (mmol)  *
 * that layer i deposited on its bed strip during the call: what leaves a     *
 * sinking layer and does not pass into the layer below (K19, K20). The host  *
 * hands it to the sediment zones in contact with each strip.                 *
 ******************************************************************************/
void doMobilityDep(const CINTEGER *N_in,   // number of vertical layers
              const AED_REAL *dt_in,    // time step (s)
              const AED_REAL *h,        // array of layer thicknesses (m)
              const AED_REAL *A,        // array of layer areas (m^2)
              const AED_REAL *V,        // array of layer volumes (m^3)
              const AED_REAL *ww,       // array of vertical mobility speeds (m/s)
              const AED_REAL *min_C_in, // minimum concentration (mmol/m^3)
              AED_REAL *mcc,            // array of cell concentrations (mmol/m^3)
              AED_REAL *dep)            // deposited mass per layer (mmol), added to; may be NULL
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

    AED_REAL *Y;     // total mass of variable (mmol) in layer
    AED_REAL *mov;   // mass leaving each layer in this sub-step (mmol)

    AED_REAL dtMax, tdt, tmp, w;

/*----------------------------------------------------------------------------*/
    if ( N < 1 ) return;

    Y   = calloc(N, sizeof(AED_REAL));
    mov = calloc(N, sizeof(AED_REAL));

    /**************************************************************************
     * determine mobility timestep, the maximum time step for which:          *
     *  - no layer loses more than it holds  (|w| dt A[i] <= V[i]), and       *
     *  - particles do not pass through more than one layer (|w| dt <= h of   *
     *    this layer and of the layer they move into)                         *
     **************************************************************************/
    dtMax = dt;
    for (i = 0; i < N; i++) {
        Y[i] = mcc[i] * V[i];
        w = fabs(ww[i]);
        if ( w <= 0. ) continue;

        if ( A[i] > 0. && (w * dt * A[i]) > V[i] ) {
            tdt = V[i] / (w * A[i]);
            if ( tdt < dtMax ) dtMax = tdt;
        }
        if ( (w * dt) > h[i] ) {
            tdt = h[i] / w;
            if ( tdt < dtMax ) dtMax = tdt;
        }
        if ( ww[i] > 0. ) {
            if ( (i < N-1) && (w * dt) > h[i+1] ) {
                tdt = h[i+1] / w;
                if ( tdt < dtMax ) dtMax = tdt;
            }
        } else {
            if ( (i > 0) && (w * dt) > h[i-1] ) {
                tdt = h[i-1] / w;
                if ( tdt < dtMax ) dtMax = tdt;
            }
        }
    }

    /**************************************************************************/

    tdt = dtMax;
    if ( dtMax > 0.0001 ) { // in case a cell height was very small or the ww value very large
        do  {   // do this in steps of dtMax, but at least once
            // each time tdt is dtMax, except, possibly, the last which is
            // whatever was left.
            if ( (dt -= dtMax) < 0. ) tdt = dtMax + dt;

            // 1) how much leaves each layer, from the start-of-step state,
            //    never taking a layer below its minimum
            for (i = 0; i < N; i++) {
                mov[i] = 0.;
                if ( ww[i] > 0. && i == N-1 ) continue;  // nothing rises out of the surface
                if ( ww[i] == 0. ) continue;
                mov[i] = fabs(ww[i]) * tdt * A[i] * mcc[i];
                tmp = Y[i] - min_C * V[i];
                if ( mov[i] > tmp ) mov[i] = ( tmp > 0. ) ? tmp : 0.;
            }
            // 2) move it: up into the layer above, or down into the layer below
            //    (the share over its bottom face) with the rest deposited
            for (i = 0; i < N; i++) {
                if ( mov[i] == 0. ) continue;
                Y[i] -= mov[i];
                if ( ww[i] > 0. )
                    Y[i+1] += mov[i];
                else if ( i > 0 ) {
                    Y[i-1] += mov[i] * (A[i-1] / A[i]);
                    if ( dep != NULL ) dep[i] += mov[i] * (1. - A[i-1] / A[i]);
                } else if ( dep != NULL )
                    dep[i] += mov[i];             // the bottom layer: all of it
            }
            for (i = 0; i < N; i++) mcc[i] = Y[i] / V[i];
        } while ( dt > 0. );
    }

    /**************************************************************************/
    free(Y); free(mov);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
