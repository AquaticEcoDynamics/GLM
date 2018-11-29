/******************************************************************************
 *                                                                            *
 * glm_mixu.c                                                                 *
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
 #                                                                            *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include "glm.h"

#include "glm_types.h"
#include "glm_const.h"
#include "glm_globals.h"

#include "aed_time.h"
#include "glm_util.h"
#include "glm_mixu.h"


/******************************************************************************
 *                                                                            *
 * VMsum                              [UPDATED]                               *
 * Tsum                               [UPDATED]                               *
 * Ssum                               [UPDATED]                               *
 * VMLOC  :- combined volumetric mass [RETURNED]                              *
 * TMLOC  :- Combined temperature     [RETURNED]                              *
 * SMLOC  :- Combined salinity        [RETURNED]                              *
 * DFLOC  :- Density                  [RETURNED]                              *
 * idx    :- layer index              [USED]                                  *
 *                                                                            *
 ******************************************************************************/
void add_this_layer(AED_REAL *VMsum, AED_REAL *Tsum,
                    AED_REAL *Ssum, AED_REAL *VMLOC, AED_REAL *TMLOC,
                    AED_REAL *SMLOC, AED_REAL *DFLOC, int idx)
{
    AED_REAL Layer_Mass;

    Layer_Mass = Lake[idx].Density * Lake[idx].LayerVol;

    *VMsum += Layer_Mass;
    *Tsum  += (Lake[idx].Temp * Layer_Mass);
    *Ssum  += (Lake[idx].Salinity * Layer_Mass);

    *VMLOC = *VMsum;            //# Combined volumetic mass
    *TMLOC = *Tsum / (*VMLOC);  //# Combined temperature
    *SMLOC = *Ssum / (*VMLOC);  //# Combined salinity

    *DFLOC = calculate_density(*TMLOC,*SMLOC);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 * This routine assigns the mean mixed layer properties to all of the         *
 * layers in the epilimnion and increments the layer pointers j, k            *
 ******************************************************************************/
void average_layer(int *j, int *k,
                            AED_REAL MeanTemp, AED_REAL MeanSalt, AED_REAL Dens)
{
    int i, jl = *j;

    for (i = jl; i < NumLayers; i++) {
        Lake[i].Temp = MeanTemp;
        Lake[i].Salinity = MeanSalt;
        Lake[i].Density = Dens;
    }

    (*j) -= 1;
    (*k) -= 1;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
static void one_layer(int i, AED_REAL *xx, AED_REAL *dxx)
{
    AED_REAL x, y;

    int ij;

    x = Lake[i].Height * MphInc;
    y = AMOD(x, 1.0);
    ij = (x - y);
    if (ij > Nmorph) {
        y += (ij - Nmorph);
        ij = Nmorph;
    }
    if (ij > 0) {  // CAB limit to 0 min.
        ij--; // offset for 0 based index
        Lake[i].Vol1 = xx[ij] + y * dxx[ij];
        Lake[i].LayerArea = MphLevelArea[ij] + y * dMphLevelArea[ij];
    } else {
        Lake[i].Vol1 = xx[0];
        Lake[i].LayerArea = MphLevelArea[0];
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 * From the given physical data, this function evaluates arrays of depths     *
 * and areas corresponding to an array of volumes (icode=2) or arrays of      *
 * volume and areas from depths (icode=1); starting at layer LNU              *
 ******************************************************************************/
void resize_internals(int icode, int lnu)
{
   AED_REAL VolSum, x, y;

   int i, j, ij, k, l, ln;

//-------------------------------------------------------------------------------

    ln = lnu;
    if (icode == 1) {
        /**********************************************************************
         * Find volumes and areas given depths; ij points to storage table    *
         * entry closest to but less than the given depth, Y is the relative  *
         * location of DEPTH between storage table entries ij and ij+1.       *
         * Begin by calculating the total volume in the layer structure.      *
         **********************************************************************/
        VolSum = 0.0;
        for (k = 0; k < NumInf; k++)
            VolSum += Inflows[k].TotIn;

        while(NumLayers >= 1) { // stop at 1
            one_layer(surfLayer, MphLevelVol, dMphLevelVol);

            Lake[surfLayer].Vol1 -= VolSum;

            /* compute to one below current surface layer */
            for (i = ln; i < surfLayer; i++)
                one_layer(i, MphLevelVoldash, dMphLevelVolda);

            if (surfLayer <= botmLayer)
                break;
            if (Lake[surfLayer].Vol1 > Lake[surfLayer-1].Vol1)
                break;

            Lake[surfLayer-1].Vol1 = Lake[surfLayer].Vol1 + VolSum;
            Lake[surfLayer-1].LayerArea = Lake[surfLayer].LayerArea;
            Lake[surfLayer-1].Height = Lake[surfLayer].Height;
            Lake[surfLayer-1].Temp =
                    combine(Lake[surfLayer-1].Temp, Lake[surfLayer-1].LayerVol, Lake[surfLayer-1].Density,
                            Lake[surfLayer].Temp,   Lake[surfLayer].LayerVol,   Lake[surfLayer].Density);
            Lake[surfLayer-1].Salinity =
                    combine(Lake[surfLayer-1].Salinity, Lake[surfLayer-1].LayerVol, Lake[surfLayer-1].Density,
                            Lake[surfLayer].Salinity,   Lake[surfLayer].LayerVol,   Lake[surfLayer].Density);
            Lake[surfLayer-1].Density = calculate_density(Lake[surfLayer-1].Temp, Lake[surfLayer-1].Salinity);

            NumLayers--;
            if (surfLayer < ln) {
                if (--ln < botmLayer) {
                    fprintf(stderr,"surface layer less than bottom layer.\n");
                    exit(1);
                }
            }
        }
    } else {
        /**********************************************************************
         * calculate depths given volumes; J points to storage table volume   *
         * closest to but not greater than the given layer volume             *
         **********************************************************************/
        VolSum = Lake[surfLayer].Vol1;
        for (i = 0; i < NumInf; i++)
            VolSum += Inflows[i].TotIn;

        j = 0;
        while (j < Nmorph) {
            if (VolSum <= MphLevelVol[j]) {
                j--;
                break;
            }
            j++;
        }
        if (j >= Nmorph) j = Nmorph - 1;
        Lake[surfLayer].Height = ((j+1) + ((VolSum - MphLevelVol[j]) / dMphLevelVol[j])) / MphInc;

        if (lnu <= botmLayer) {
            l = 0;
            j = 0;

            /* find lowest layer (l) whose volume exceeds the first table entry */
            while (l <= surfLayer && Lake[l].Vol1 <= MphLevelVoldash[0]) {
                Lake[l].Height = (Lake[l].Vol1 / MphLevelVoldash[0]) / MphInc;
                l++;
            }
        } else {
            x = Lake[lnu-1].Height * MphInc;
            y = AMOD(x, 1.0);
            j = (x - y) - 1;
            l = lnu;
        }

        /* compute to one below current surface */
        for (k = l; k < surfLayer; k++) {
            while (j < Nmorph) {
                if (Lake[k].Vol1 <= MphLevelVoldash[j]) {
                    j--;
                    break;
                }
                j++;
            }
            if (j >= Nmorph) j = Nmorph - 1;
            if (j < 0) j = 0;

            Lake[k].Height = ((j+1) + ( Lake[k].Vol1 - MphLevelVoldash[j] ) / dMphLevelVolda[j]) / MphInc;
        }

        //# determine areas
        for (i = lnu; i <= surfLayer; i++) {
            x = Lake[i].Height * MphInc;
            y = AMOD(x, 1.0);
            ij = (x - y) - 1;
            if (ij >= Nmorph) ij = Nmorph - 1;

            if (ij >= 0) Lake[i].LayerArea = MphLevelArea[ij] + y * dMphLevelArea[ij];
            else         Lake[i].LayerArea = MphLevelArea[0] * Lake[i].Height * MphInc;
        }
    }

    //# calculate layer volumes and mean depths
    ln = lnu;
    if (lnu == botmLayer) {
        Lake[botmLayer].LayerVol  = Lake[botmLayer].Vol1;
        Lake[botmLayer].MeanHeight = Lake[botmLayer].Height / 2.;
        ln = lnu+1;
    }

    for (i=ln; i <= surfLayer; i++) {
        Lake[i].LayerVol  = Lake[i].Vol1 - Lake[i-1].Vol1;
        Lake[i].MeanHeight = (Lake[i].Height + Lake[i-1].Height) / 2.;
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
