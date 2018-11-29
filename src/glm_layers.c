/******************************************************************************
 *                                                                            *
 * glm_layers.c                                                               *
 *                                                                            *
 * Contains layer utility routines                                            *
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
 *  Algorithms implmented in the below functions are adapted from those       *
 *  originally developed by the Centre for Water Research, University of      *
 *  Western Australia. Readers are referred to Imberger and Patterson (1981)  *
 *  for a comprehensive overview of the original mixing model.                *
 *                                                                            *
 *  Imberger, J. and Patterson, J.C., 1981. A dynamic reservoir simulation    *
 *    model-DYRESM:5. In: H.B. Fisher (ed.), Transport Models for Inland      *
 *    and Coastal Waters. Academic Press, New York: 310-361.                  *
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
#include <math.h>

#include "glm.h"

#include "glm_types.h"
#include "glm_const.h"
#include "glm_globals.h"

#include "glm_util.h"
#include "glm_mixu.h"

//#define dbgprt(...) fprintf(stderr, __VA_ARGS__)
#define dbgprt(...) /* __VA_ARGS__ */

/******************************************************************************
 * This subroutine checks the reservoir layer structure for compliance        *
 * with the specified volume and depth limits.  Adjustments are made          *
 * as required. Layer structure is checked for minimum limits first,          *
 * combining the  checked layer with the smallest adjacent layer if           *
 * necessary. This may result in the formation of a layer exceeding           *
 * maximum limits.Layer structure is checked for maximum limits,              *
 * splitting the checked layer if necessary.  After splitting, the            *
 * resulting layers will all be greater than their minimum limits             *
 * provided VMax >= 2 VMin.  The default value is VMax = 2 VMin.              *
 * split depths ! volumes                                                     *
 *                                                                            *
 * DMax                maximum allowable layer thickness                      *
 * DMin                minimum allowable layer thickness                      *
 * KB                  first layer above split or mixed layers                *
 * KT                  new top layer number                                   *
 * VMax                max allowed volume                                     *
 * VMin                min allowed volume                                     *
 ******************************************************************************/
void check_layer_thickness(void)
{
//LOCALS
    AED_REAL D;
    AED_REAL DELDP;
    AED_REAL V;         // Split volume
    AED_REAL Vdown;     // Volume of layer below amalgamation
    AED_REAL Vup;       // Volume of layer above amalgamation

    int i;
    int iadd;
    int j;
    int wqidx;
    int jold;
    int k;
    int KB;
    int KLAST;
    int KT;
    int M;        // Number of layers after splitting
    int VSUMCHK;

/*----------------------------------------------------------------------------*/
    dbgprt(" CHKLAY 01 lake[44].depth = %20.15f\n", Lake[44].Height);

    //# Check against vmin
    KLAST=botmLayer;
    // while (1) { //
    while (NumLayers > 1) { // stop at 1 layer
        for (i = KLAST; i <= surfLayer; i++) {
            if (i == botmLayer)
                 DELDP = Lake[i].Height;
            else
                 DELDP = Lake[i].Height - Lake[i-1].Height;
            if ((Lake[i].LayerVol < VMin) && (DELDP < DMin)) break;
        }

        if (i > surfLayer) break;

        // Layer i is amalgamated with its smallest neighbour
        if (i == botmLayer) {
            Vup = zero;
            Vdown = 1.0;
        } else {
            if (i == surfLayer) {
                Vup = 1.0;
                Vdown = zero;
            } else {
                Vup = Lake[i+1].LayerVol;
                Vdown = Lake[i-1].LayerVol;
            }
        }

        j = i;
        if (Vup > Vdown) j = i-1;

        Lake[j].Salinity = combine(Lake[j].Salinity,   Lake[j].LayerVol,   Lake[j].Density,
                                   Lake[j+1].Salinity, Lake[j+1].LayerVol, Lake[j+1].Density);
        Lake[j].Temp = combine(Lake[j].Temp,   Lake[j].LayerVol,   Lake[j].Density,
                               Lake[j+1].Temp, Lake[j+1].LayerVol, Lake[j+1].Density);

        for (wqidx = 0; wqidx < Num_WQ_Vars; wqidx++)
            _WQ_Vars(wqidx,j) = combine_vol(_WQ_Vars(wqidx,j), Lake[j].LayerVol, _WQ_Vars(wqidx,j+1), Lake[j+1].LayerVol);

        Lake[j].Density = calculate_density(Lake[j].Temp, Lake[j].Salinity);
        Lake[j].LayerVol = Lake[j].LayerVol + Lake[j+1].LayerVol;
        Lake[j].Height = Lake[j+1].Height;
        Lake[j].Vol1 = Lake[j+1].Vol1;
        Lake[j].LayerArea = Lake[j+1].LayerArea;
        Lake[j].Epsilon = Lake[j+1].Epsilon;
        KLAST=j;

        // Renumber layers j+2,j+3,---,surfLayer
        if (j != (surfLayer-1)) {
            KT = surfLayer-1;
            KB = j + 1;
            for (k = KB; k <= KT; k++) {
                Lake[k].Height = Lake[k+1].Height;
                Lake[k].Density = Lake[k+1].Density;
                Lake[k].Temp = Lake[k+1].Temp;
                Lake[k].Salinity = Lake[k+1].Salinity;

                for (wqidx=0; wqidx < Num_WQ_Vars; wqidx++)
                    _WQ_Vars(wqidx, k) = _WQ_Vars(wqidx, k+1);

                Lake[k].LayerVol = Lake[k+1].LayerVol;
                Lake[k].Vol1 = Lake[k+1].Vol1;
                Lake[k].LayerArea = Lake[k+1].LayerArea;
                Lake[k].Epsilon = Lake[k+1].Epsilon;
            }
        }
        NumLayers--;
    }

    // here when all layers have been checked for VMin, DMin
    if (surfLayer != botmLayer) {
        for (i = botmLayer+1; i <= surfLayer; i++)
            Lake[i].MeanHeight=(Lake[i].Height+Lake[i-1].Height)/ 2.0;
    }
    Lake[botmLayer].MeanHeight = Lake[botmLayer].Height/ 2.0;

    // check layers for VMax
    //sgs Flag to prevent top layer splitting more than once
    VSUMCHK = FALSE;
    KLAST=botmLayer;
    while(1) {
        if (VSUMCHK) return;

        for (i = KLAST; i <= surfLayer; i++) {
            if (i == botmLayer)
                DELDP=Lake[i].Height;
            else
                DELDP=Lake[i].Height-Lake[i-1].Height;

            if (i == surfLayer) VSUMCHK = TRUE;

            if (Lake[i].LayerVol > VMax || DELDP > DMax) break;
        }

        // return to calling program when all layers have been checked
        if (i > surfLayer) return;

        // layer i is split into M layers
        M = 2;
        while (1) {
            V = Lake[i].LayerVol/M;
            D = DELDP/M;
            if (V <= VMax && D <= DMax) break;
            if (Lake[surfLayer].Height<0.3) break;
            M++;

            // if M+surfLayer is greater than the max no. of layers, a mistake will occur
            //  - an array bound error
            if (M + NumLayers > MaxLayers) {
                fprintf(stderr, "Array bounds error - too many layers. NumLayers = %d, M = %d\n", NumLayers, M);
                fprintf(stderr, "i = %d V = %20.15f VMax = %20.15f D = %20.15f DMax = %20.15f\n",
                                      i, V, VMax, D, DMax);
                exit(1);
            }
        }

        // renumber layers above split. iadd is the number of added layers
        // j is the new layer number, jold is the old layer number
        // include water quality
        iadd = M - 1;
        KLAST = i;
        if (i != surfLayer) {
            KT = surfLayer+iadd;
            KB = i+M;
            for (j = KT; j >= KB; j--) {
                jold = j-iadd;
                Lake[j].Vol1 = Lake[jold].Vol1;
                Lake[j].Density = Lake[jold].Density;
                Lake[j].Temp = Lake[jold].Temp;
                Lake[j].Salinity = Lake[jold].Salinity;

                for (wqidx = 0; wqidx < Num_WQ_Vars; wqidx++)
                    _WQ_Vars(wqidx,j) = _WQ_Vars(wqidx,jold);

                Lake[j].LayerVol = Lake[jold].LayerVol;
                Lake[j].Epsilon = Lake[jold].Epsilon;
            }
        }

        // process the added layers, include water quality
        for (k=i; k <= i+iadd; k++) {
            Lake[k].LayerVol = V;
            if (k == botmLayer)
                Lake[k].Vol1=Lake[k].LayerVol;
            else
                Lake[k].Vol1=Lake[k-1].Vol1+Lake[k].LayerVol;

            Lake[k].Density = Lake[i].Density;
            Lake[k].Temp = Lake[i].Temp;
            Lake[k].Salinity = Lake[i].Salinity;

            for (wqidx = 0; wqidx < Num_WQ_Vars; wqidx++)
                _WQ_Vars(wqidx,k)=_WQ_Vars(wqidx,i);

            Lake[k].Epsilon = Lake[i].Epsilon;
        }
        NumLayers += iadd;

        // get new depths for layers i thru surfLayer
        resize_internals(2,i);

    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 #CAB#  These comments are clearly WRONG!                                     *
 * This subroutine finds the level of Neutral Bouyancy for a given inflow     *
 * and returns the layer number (i), the half-thickness (B0), basin length    *
 * at the intrusion midpoint (AL), basin width at the intrusion               *
 * midpoint, and the mean intrusion velocity (IntrusionVelocity) in m/s.      *
 ******************************************************************************/
void insert(AED_REAL q, AED_REAL di, AED_REAL bsl, AED_REAL temp, AED_REAL salt,
                             AED_REAL *wqx, int ntims, AED_REAL *width, int *ll)
{
    AED_REAL AHLE;
    AED_REAL LenAtInsert;
    AED_REAL ALSQ;
    AED_REAL BL;
    AED_REAL B0;
    AED_REAL DBB,DELT,DELB;
    AED_REAL DHLE;
    AED_REAL DT;

    AED_REAL *DVR;

    AED_REAL DZ;
    AED_REAL F;
    AED_REAL GD;
    AED_REAL GR;
    AED_REAL R;
    AED_REAL TDASH;
    AED_REAL IntrusionVelocity;
    AED_REAL Viscosity;
    AED_REAL XN;
    AED_REAL XNSQ;
    AED_REAL ZP,ZT,ZB;

    int wqidx;
    int i,j,k;
    int iz;
    int JB,JT;
    int KX,KY;
    int BotLayerIndex;  //
    int TopLayerIndex;

/*----------------------------------------------------------------------------*/

    DVR = calloc(MaxLayers, sizeof(AED_REAL));

    DELT=0.;
    DELB=0.;

    for (iz = botmLayer; iz <= surfLayer; iz++)
       if (Lake[iz].Height > 0.) break;
    if (iz > surfLayer) iz = surfLayer;

    DHLE = Lake[iz].Height;
    AHLE = Lake[iz].LayerArea;

    // FInd the appropriate layer depth, based on the inflow density
    if (di <= Lake[surfLayer].Density) {
        // Here for surface overflow
        i = surfLayer;
        *ll = i;
        B0 = (Lake[surfLayer].Height-Lake[surfLayer-1].Height)/ 2.0;
        LenAtInsert = LenAtCrest;
        if (*width <= 1E-7) *width = Lake[surfLayer].LayerArea / LenAtInsert;
        IntrusionVelocity = q / ( 2.0 * B0 * (*width));
    } else {
        // Find level of neutral buoyancy
        for (i = surfLayer; i > botmLayer; i--)
            if (di <= Lake[i-1].Density) break;

        if (i <= botmLayer) {
            //  Here for underflow
            i = botmLayer;
            *ll = i;
            LenAtInsert = DHLE/sin(bsl);
            if ((*width) <= 1E-7) (*width) = AHLE / LenAtInsert;
            B0 = (DHLE/ 2.0)/ 2.0;
            IntrusionVelocity = q / ( 2.0 * B0 * (*width));
        } else {
            //  Here for intrusion
            JT = i;
            *ll = i;
            JB = i-1;
            LenAtInsert = Lake[i].Height/sin(bsl);
            ALSQ = sqr(LenAtInsert);
            if ((*width) <= 1E-7) (*width) = Lake[i].LayerArea/LenAtInsert;

            while(1) {
                DT = Lake[JT].MeanHeight;
                DBB = Lake[JB].MeanHeight;
                DZ = DT - DBB;
                XNSQ = g*(Lake[JB].Density-Lake[JT].Density)/(di*DZ);
                if (XNSQ  <=  zero)
                    //  Here for unstable stratification
                    BL=LenAtInsert;
                else {
                    //  here for stable stratification
                    XN = sqrt(XNSQ);
                    F = q/((*width)*ntims*XN*ALSQ);
                    Viscosity = Lake[i].Epsilon * 20.0;
                    if (Viscosity <= 0.) Viscosity=Visc;
                    GR = XNSQ*sqr(ALSQ)/sqr(Viscosity);

                    R = F*pow(GR,(1.0/3.0));
                    TDASH=ntims*XN/(pow(GR,(1.0/6.0)));
                    R /= TDASH;
                    if (R > 1.)
                        BL = 0.44*TDASH*LenAtInsert*sqrt(R);
                    else
                        BL = 0.57*LenAtInsert*pow(R, (3.0/ 2.0))*pow((TDASH/R), (5.0/6.0));

                    BL = MIN(BL,LenAtInsert);
                    if (BL < 1.0) BL = 1.0;
                }

                // B0 is 1/2 the intrusion thickness
                B0 = q/((*width)*BL);
                if (B0 > DZ) {
                    if ( !((JT == surfLayer) && (JB == botmLayer)) ) {
                        if (JT != surfLayer) JT++;
                        if (JB != botmLayer) JB--;
                        continue;
                    }
                }
                break;
            }

            if (Lake[i].Height < (DHLE + 1.0)) {
                LenAtInsert = DHLE/sin(bsl);
                if ((*width) <= 1E-7) (*width) = AHLE / LenAtInsert;
                B0 = (DHLE+ 2.0)/ 2.0;
            }
            IntrusionVelocity = q/( 2.0*B0*(*width));
        }
    }

    //  Mix the inflow with the appropriate layers.
    TopLayerIndex=i;
    ZP=Lake[i].MeanHeight;
    if ( ! (i > botmLayer && i < surfLayer) ) {
        //  Here for underflow, overflow, or fully mixed
        BotLayerIndex=i;
        DVR[i]=q;
    } else {
        // Here for intrusion
        while(1) {
            TopLayerIndex++;
            DELT = 0.0;
            if (TopLayerIndex  !=  surfLayer) {
                GD=gprime(Lake[TopLayerIndex].Density,Lake[i].Density);
                if (GD > zero) DELT = 0.15 * pow((IntrusionVelocity/(ntims)),2) / GD;
                if (DELT > (Lake[TopLayerIndex].MeanHeight-ZP) || GD <= zero) continue;
                if (Lake[TopLayerIndex-1].Height > (ZP+DELT)) TopLayerIndex--;
            }
            break;
        }

        ZT=Lake[TopLayerIndex].Height;
        BotLayerIndex=i;
        while(1) {
            BotLayerIndex--;
            if (BotLayerIndex != botmLayer) {
                GD = gprime(Lake[i].Density,Lake[BotLayerIndex].Density);
                if (GD > zero) DELB = 0.15 * sqr((IntrusionVelocity/(ntims))) / GD;
                if (DELB > (ZP-Lake[BotLayerIndex].MeanHeight) || GD <= zero) continue;
                if (Lake[BotLayerIndex].Height < (ZP-DELB)) BotLayerIndex++;
            }
            break;
        }

        ZB = zero;
        if (BotLayerIndex > botmLayer) ZB = Lake[BotLayerIndex-1].Height;
        if (BotLayerIndex == TopLayerIndex) {
            // Here if intrusion is entirely within layer i
            DVR[BotLayerIndex]=q;
        } else {
            // Aportion inflow amongst layers BotLayerIndex,BotLayerIndex+1,---,TopLayerIndex
            DELT = ZT - ZP;
            DELB = ZP - ZB;
            if (BotLayerIndex != i) {
                DVR[BotLayerIndex]=q*(Lake[BotLayerIndex].Height-ZB+DELB*sin(Pi*(Lake[BotLayerIndex].Height-ZP)/DELB)/Pi)/(ZT-ZB);
                if (BotLayerIndex != (i-1)) {
                    KX = BotLayerIndex+1;
                    KY = i-1;
                    for (k = KX; k <= KY; k++)
                        DVR[k]=q*(Lake[k].Height-Lake[k-1].Height+DELB *
                            (sin(Pi*(ZP-Lake[k-1].Height)/DELB)-sin(Pi*(ZP-Lake[k].Height)/DELB))/Pi)/(ZT-ZB);
                }
            }
            DVR[i]=q*(ZP-Lake[i-1].Height+DELB*sin(Pi*(ZP-Lake[i-1].Height)/DELB)/Pi)/(ZT-ZB);
            DVR[i]=DVR[i]+q*(Lake[i].Height-ZP+DELT*sin(Pi*(Lake[i].Height-ZP)/DELT)/Pi)/(ZT-ZB);
            if (TopLayerIndex != i) {
                if (TopLayerIndex != (i+1)) {
                    KX=i+1;
                    KY=TopLayerIndex-1;
                    for (k = KX; k<=KY; k++)
                        DVR[k]=q*(Lake[k].Height-Lake[k-1].Height+DELT*(sin(Pi*(Lake[k].Height-ZP)/
                                  DELT)-sin(Pi*(Lake[k-1].Height-ZP)/DELT))/Pi)/(ZT-ZB);
                }
                DVR[TopLayerIndex]=q*(ZT-Lake[TopLayerIndex-1].Height+DELT*sin(Pi*(ZP-Lake[TopLayerIndex-1].Height)/DELT)/Pi)/(ZT-ZB);
            }
        }
    }

    // Insert new intriusion layer into water column layers and adjust layer properties
    // (i.e., "combine"). Include water quality attributes when combining.
    for (k = BotLayerIndex; k <= TopLayerIndex; k++) {
        // Temp and Salinity
        Lake[k].Temp = combine(Lake[k].Temp,Lake[k].LayerVol,Lake[k].Density,temp,DVR[k],di);
        Lake[k].Salinity = combine(Lake[k].Salinity,Lake[k].LayerVol,Lake[k].Density,salt,DVR[k],di);
        // Water QUality
        for (wqidx = 0; wqidx < Num_WQ_Vars; wqidx++)
            _WQ_Vars(wqidx,k) = combine_vol(_WQ_Vars(wqidx,k),Lake[k].LayerVol,wqx[wqidx],DVR[k]);
        // Density and Volume
        Lake[k].Density=calculate_density(Lake[k].Temp,Lake[k].Salinity);
        Lake[k].LayerVol=Lake[k].LayerVol+DVR[k];
    }

    Lake[botmLayer].Vol1 = Lake[botmLayer].LayerVol;
    if (surfLayer != botmLayer) {
        for (j = (botmLayer+1); j <= surfLayer; j++)
            Lake[j].Vol1 = Lake[j-1].Vol1 + Lake[j].LayerVol;
    }

    free(DVR);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
