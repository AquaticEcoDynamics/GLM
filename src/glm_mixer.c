/******************************************************************************
 *                                                                            *
 * glm_mixer.c                                                                *
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
#include "glm_mixu.h"
#include "glm_util.h"

#include "aed_time.h"

//#define dbgprt(...) fprintf(stderr, __VA_ARGS__)
#define dbgprt(...) /* __VA_ARGS__ */

#include "glm_debug.h"

/******************************************************************************/

#define DEEPENED_BOT     1
#define MOMENTUM_CUT     2
#define IS_MIXED         3

/*============================================================================*/

AED_REAL DepMX   = 0.;

static   AED_REAL Epi_dz;     //# Thickness of epilimnion [m]
static   AED_REAL MeanSalt;   //# MeanSalt ... mean salinity
static   AED_REAL MeanTemp;   //# MeanTemp ... mass averaged mean temperature of epilimnion

static   AED_REAL PrevThick = 0.;   //# mixed layer thickness from previous time step

static   AED_REAL gPrimeTwoLayer = 0.;  //# Reduced gravity for int wave estimate

static   AED_REAL Energy_RequiredMix = 0.;   //# Energy required to entrain next layer into the epilimnion
static   AED_REAL Energy_AvailableMix = 0.;  //# Total available energy to mix (carries over from previous timesteps)

static   AED_REAL Vol_Epi = 0.;  //# Volume of Epilimnion (surface layer after Kelvin-Helmholtz) m3
static   AED_REAL Mass_Epi = 0.; //# Sigma mass of Epilimnion (surface layer after Kelvin-Helmholtz) kg

static   AED_REAL OldSlope   = 0.;
static   AED_REAL Time_end_shear   = 0.;  //# Time left before shear cut off [hours]
static   AED_REAL Time_start_shear = 0.;  //# Time count since start of sim for shear period start [hours]
static   AED_REAL Time_count_end_shear = 0.;  //# Time count since start of sim for shear period end [hours]
static   AED_REAL Time_count_sim     = 0.;  //# Time count since start of simulation [hours]

static   AED_REAL Half_Seiche_Period = 0.; //# One half the seiche period
static   AED_REAL Thermocline_Height = 0.; //# Height at the top of the metalimnion [m]
static   AED_REAL FO      = 0.;
static   AED_REAL FSUM    = 0.;
static   AED_REAL u_f     = 0.;
static   AED_REAL u0      = 0.;
static   AED_REAL u_avg   = 0.;


extern AED_REAL coef_mix_KH;
extern AED_REAL coef_mix_conv,coef_wind_stir,coef_mix_shear,coef_mix_turb;

/*============================================================================*/

static AED_REAL kelvin_helmholtz(int *Meta_topLayer,
                         int *Epi_botmLayer, AED_REAL Dens, AED_REAL *WQ_VarsM);


/******************************************************************************
 * Perform mixing by calculating potential energy                             *
 * released by bouyancy flux and surface wind stress                          *
 ******************************************************************************/
int mixed_layer_deepening(AED_REAL *WQ_VarsM, int Mixer_Count, int *_Meta_topLayer, AED_REAL *_Dens_Epil)
{
    const AED_REAL twelve = 12.0,
                   twfour = 24.0,
                   half = 0.500,
                   tdfac = 8.33E-4;

    //* Wind parameters
    AED_REAL WindSpeedX;         //# Actual wind speed, accounting for wind factor or ice [m s-1]
    AED_REAL U_star;             //# U*, wind induced surface water shear speed [m s-1]
    AED_REAL U_star_sqr;         //# U*^2 [m2 s-2]
    AED_REAL U_star_cub;         //# U*^3 [m3 s-3]

    //* Epilimnion variables
    AED_REAL Dens_Epil;          //# Mean epilimnion density [kg/m3]
    AED_REAL Epi_Thick;          //# Effective thickness of the epilimnion (Volume/Area)
    AED_REAL Epilimnion_Mid_Ht;  //# Epilimnion height measured to middle of epilimnion
    AED_REAL dMdz;               //# Delta mass vertical gradient (kg/m), for w*
    AED_REAL q_cub;              //# q^3 ...  q*3 = w*3 + C_w u*3
    AED_REAL q_sqr;              //# q^2 ...  (q*3)^(2/3)

    //* Metalimnion variables
    AED_REAL LengthAtThermo;     //# Effective length of lake at thermocline
    AED_REAL IntWaveSpeed;       //# Wave speed along the thermocline for a two layer fluid
    AED_REAL ZeroMom;            //# Metalimnion 0th momement of density above the bottom (dz * density)
    AED_REAL FirstMom;           //# Metalimnion 1st momement of density above the bottom (dz * height * density)
    AED_REAL delzkm1;            //# Delta height (layer thickness) of layer getting entrained

    //* Hypolimnion variables
    AED_REAL Hypl_Thick;         //# Effective thickness of the hypolimnion (Volume/Area)
    AED_REAL Vol_Hypl;           //# Volume of hypolimnion [m^3]
    AED_REAL Dens_Hypl;          //# Mean hypolimnion density [kg/m3]
    AED_REAL Hypl_Mass;          //* Mass of layers contributing to hypolimnion

    //* Energy vars
    AED_REAL Energy_Conv;        //# Energy released by convective overturn (Kraus Turner)
    AED_REAL Energy_WindStir;    //# Energy available from wind stirring
    AED_REAL Energy_TotStir;     //# Total energy available for stirring
    AED_REAL Energy_Deepen;      //# Energy available from the rate of deepening due to shear?
    AED_REAL Energy_Shear;       //# Energy available from shear production

    //* Velocity vars (for shear calculation)
    AED_REAL redg;               //# Reduced gravity; g' or g prime
    AED_REAL u_avgSQ;            //# Average shear velocity (u_b)
    AED_REAL u_eff;              //# Effective velocity
    AED_REAL u0_old;             //# Previous base velocity
    AED_REAL u_avg_old;          //# Previous avg velocity computed
    AED_REAL deltaKH;            //# deltaKH : K-H Billow length scale
    AED_REAL del_deltaKH;        //# deldeltaKH : Change in deltaKH
    AED_REAL GPEFFC;             //# Effective gPrime
    AED_REAL del_u;              //# Change in velocity over time
    AED_REAL accn;               //# Acceleration
    AED_REAL zsml_tilda;         //# SML depth tilda
    AED_REAL Slope;              //# Slope

    //# Layer properties for combining
    AED_REAL VMsum;              //# Cumulative volumetric mass of mixed layer
    AED_REAL Tsum;               //# Mean Temperature of cumulative volume of mixed layer
    AED_REAL Ssum;               //# Mean Salinity of cumulative volume of mixed layer

    int i, wqvidx;               //# water quality variable index
    int Meta_topLayer;           //# Index for top layer of hypolimnion
    int Epi_botmLayer;           //# Index for bottom layer of epilimnion

    int loop_count = 0; // for debugging only

    /**************************************************************************
     * Initialise                                                             *
     **************************************************************************/
    Tsum  = zero;
    Ssum  = zero;
    VMsum = zero;
    ZeroMom = zero;
    FirstMom = zero;

//#if DEBUG
    _dbg_mix_init_fields();
    /* set to missing for debugger */
    redg = missing; Epi_Thick = missing; dMdz = missing; q_cub = missing;
    LengthAtThermo = missing; Hypl_Thick = missing; Dens_Hypl = missing;
    Energy_Conv = missing; Energy_WindStir = missing; Energy_TotStir = missing;
    Energy_Deepen = missing; Energy_Shear = missing; del_u = missing;
    u_avgSQ = missing; u_eff = missing; u0_old = missing; u_avg_old = missing;
    deltaKH = missing; del_deltaKH = missing; GPEFFC = missing; accn = missing;
    zsml_tilda = missing; Slope = missing; Epilimnion_Mid_Ht = missing;
    q_sqr = missing; IntWaveSpeed = missing; delzkm1 = missing;
    Vol_Hypl = missing; Hypl_Mass = missing; Dens_Epil = *_Dens_Epil;
    Epi_botmLayer = -9; Meta_topLayer = -9;
//#endif

    /**************************************************************************
     * Determine actual surface wind accounting for wind factor or ice        *
     **************************************************************************/
    if (ice) WindSpeedX = 1e-5;
    else     WindSpeedX = MetData.WindSpeed;

    /**************************************************************************
     * Calculate shear velocity U*, U*^2 and U*^3                             *
     **************************************************************************/
    // CAB - need to debug this.
    U_star = coef_wind_drag * sqrt(WindSpeedX*WindSpeedX);
    // U_star = sqrt( coef_wind_drag * WindSpeedX * WindSpeedX );
    U_star_sqr = U_star*U_star;        //# U*^2 handy in mixing calcs
    U_star_cub = U_star*U_star*U_star; //# U*^3 handy in mixing calcs

    _DBG_MIXER_(1, 0, 0,     // step 1, before loop, loop_step
           NumLayers, Epi_botmLayer, Meta_topLayer,
           Energy_AvailableMix, Energy_RequiredMix, redg);

    /**************************************************************************
     *                                                                        *
     * STEP 1 - CONVECTIVE OVERTURN                                           *
     *                                                                        *
     * Loop from bottom to surface layer and determine the density of all layer
     * volume from layer Epi_botmLayer and above and check if greater than    *
     * density of next layer down if it is then have reached bottom of        *
     * epilimnion and break from loop.                                        *
     *                                                                        *
     * This is done following algorithm in Imberger & Patterson (1981) p326   *
     *                                                                        *
     **************************************************************************/
    for (i = botmLayer; i <= surfLayer; i++) {
        Epi_botmLayer = surfLayer - i;

        //# Add new layers to the upper mixed layer, and return updated density
        add_this_layer(&VMsum, &Tsum, &Ssum, &Mass_Epi,
                       &MeanTemp, &MeanSalt, &Dens_Epil, Epi_botmLayer);

        //# Compute the 0th & 1st moments of density (about the bottom) (Eq 32)
        if (Epi_botmLayer != botmLayer) {
            AED_REAL tRho = (Lake[Epi_botmLayer].Density - rho0) *
                            (Lake[Epi_botmLayer].Height - Lake[Epi_botmLayer-1].Height);
            ZeroMom  = ZeroMom  + tRho;
            FirstMom = FirstMom + tRho * Lake[Epi_botmLayer].MeanHeight;

            if (Dens_Epil < Lake[Epi_botmLayer-1].Density+1e-6) break;
        } else {
            AED_REAL tRho = (Lake[botmLayer].Density - rho0) * Lake[botmLayer].Height;
            ZeroMom  = ZeroMom  + tRho;
            FirstMom = FirstMom + tRho * Lake[Epi_botmLayer].MeanHeight;
        }
    }

    _DBG_MIXER_(1, 2, 0,     // step 1, after loop, loop_step
           NumLayers, Epi_botmLayer, Meta_topLayer,
           Energy_AvailableMix, Energy_RequiredMix, redg);

    /**************************************************************************
     * Epi_botmLayer is now the bottom layer of the epilimnion / mixed layer  *
     * Dens_Epil is the density of the new epilimnion / mixed layer.          *
     **************************************************************************/

    /**************************************************************************
     * Meta_topLayer is now the top layer of the metalimnion                  *
     * Epi_botmLayer is the bottom layer of the epilimnion                    *
     * Therefore, Meta_topLayer+1 == Epi_botmLayer                            *
     *                                                                        *
     * Test if mixing hit the bottom, and update time counters                *
     **************************************************************************/

//  if ( Epi_botmLayer == botmLayer ) {
        //# This means lake is fully mixed so set all layers to mean properties
        for (i = Epi_botmLayer; i < surfLayer; i++) {
            Lake[i].Temp = MeanTemp;
            Lake[i].Salinity = MeanSalt;
            Lake[i].Density = Dens_Epil;
        }
//  }

    if (Epi_botmLayer != botmLayer) //# stratified lake
        //# Epilimnion mean height measured to middle of epilimnion
        Epilimnion_Mid_Ht = (Lake[surfLayer].Height + Lake[Epi_botmLayer-1].Height) / 2.0;
    else //# fully mixed epilimnion height == mid lake height
        //# Epilimnion mean height is the middle of the single layer
        Epilimnion_Mid_Ht = (Lake[surfLayer].Height) / 2.0;

    Meta_topLayer = Epi_botmLayer - 1;
    if (Epi_botmLayer == botmLayer) {
        //# This means that lake fully mixed: exit now with DEEPENED_BOT
        //# Add num_hours to sim time counter (in hours since sim)
        Time_count_sim += noSecs / SecsPerHr;
        *_Dens_Epil = Dens_Epil; *_Meta_topLayer = Meta_topLayer;
        return DEEPENED_BOT;
    }
    //# Update the delta mass vertical gradient (Eq 32 in Imberger and Patterson)
    dMdz = (FirstMom-Epilimnion_Mid_Ht*ZeroMom);

    /**************************************************************************
     * Energy_Conv measures energy released by convective overturn            *
     * (ie. cooled dense water falling; Kraus Turner deepening)               *
     **************************************************************************/
    Energy_Conv = half * coef_mix_conv * g * dMdz/(Dens_Epil*noSecs)*noSecs;
    if (Energy_Conv < zero) Energy_Conv = zero;

    /**************************************************************************
     *                                                                        *
     * STEP 2 - STIRRING                                                      *
     *                                                                        *
     * This algorithm computes the total energy available for stirring and    *
     * adds the amount to any previous amount that is stored up               *
     * This includes energy from w_star (above) plus wind stirring energy     *
     * Note coef_wind_stir is equivalent to eta                               *
     **************************************************************************/
    Energy_WindStir = half * coef_wind_stir * U_star_cub * noSecs;

    Energy_TotStir = Energy_Conv + Energy_WindStir;

    q_cub = 2.0 * Energy_TotStir / ( MAX(coef_mix_conv,0.05) * noSecs) ;

    if (q_cub <= zero) q_cub = 1e-10;
    q_sqr = pow(q_cub, (2.0/3.0));

    //# Add stirring energy to available mixing energy
    Energy_AvailableMix += Energy_TotStir;

    loop_count = 0;
    _DBG_MIXER_(2, 0, 0,     // step 2, before loop, loop_step
           NumLayers, Epi_botmLayer, Meta_topLayer,
           Energy_AvailableMix, Energy_RequiredMix, redg);

    /**************************************************************************
     * Now loop through layers using the stirring energy to mix. Computes     *
     * the energy required to mix next layer and compares with available      *
     * energy. If it hits the bottom it will return.                          *
     **************************************************************************/
    while (TRUE) {
        //# Compute energy required to mix k-1 layer
        Epi_dz = Lake[surfLayer].Height - Lake[Meta_topLayer].Height;

        if (Meta_topLayer > botmLayer)
            delzkm1 = Lake[Meta_topLayer].Height - Lake[Meta_topLayer-1].Height;
        else
            delzkm1 = Lake[Meta_topLayer].Height;

        redg = gprime(Dens_Epil, Lake[Meta_topLayer].Density);
        Energy_RequiredMix = half * (redg * Epi_dz + coef_mix_turb * q_sqr) * delzkm1 ;

        _DBG_MIXER_(2, 1, ++loop_count,     // step 2, in loop, loop_step
           NumLayers, Epi_botmLayer, Meta_topLayer,
           Energy_AvailableMix, Energy_RequiredMix, redg);

        //# Not enough energy to entrain any more layers into the mixed layer
        if (Energy_AvailableMix < Energy_RequiredMix) break;

        //# Entrain the layer Meta_topLayer (k-1) into the SML
        add_this_layer(&VMsum,&Tsum,&Ssum,&Mass_Epi,&MeanTemp,&MeanSalt,&Dens_Epil,Meta_topLayer);

        //# Now remove energy used to entrain k-1 layer into the mixed layer
        Energy_AvailableMix -= Energy_RequiredMix;
        average_layer(&Meta_topLayer,&Epi_botmLayer,MeanTemp,MeanSalt,Dens_Epil);
        if (Meta_topLayer < botmLayer) {
            //# Here if mixed layer deepening gets to bottom
            Time_count_sim += noSecs/SecsPerHr;  //#Add num_hours to sim time counter (in hours since sim)
            *_Dens_Epil = Dens_Epil;
            *_Meta_topLayer = Meta_topLayer;
            return DEEPENED_BOT;
        }
    }

    _DBG_MIXER_(2, 2, 0,     // step 2, after loop, loop_step
           NumLayers, Epi_botmLayer, Meta_topLayer,
           Energy_AvailableMix, Energy_RequiredMix, redg);

    /**************************************************************************
     *                                                                        *
     * STEP 3 - SHEAR PRODUCTION                                              *
     *                                                                        *
     * Calculate parameters needed for momentum/shear computation             *
     * IntWaveSpeed is wave speed along the thermocline for a two layer fluid *
     * Epi_Thick, Hypl_Thick are the thicknesses of the epilimnion and        *
     *     hypolimnion respectively calculated as volume/mean area            *
     * PrevThick is the mixed layer thickness from previous time step         *
     * Dens_Hypl is the mean hypolimnion density                              *
     * Dens_Epi is the epilimnion density                                     *
     **************************************************************************/

    //# Cutoff shear production if thermocline occurs below bottom
    if (Lake[Meta_topLayer].Height <= 0.) {
        *_Dens_Epil = Dens_Epil;
        *_Meta_topLayer = Meta_topLayer;
        return MOMENTUM_CUT;
    }

    //# Calculate Kraus-Turner depth (change in SML depth over dt)
    zsml_tilda = MAX( DepMX - Lake[Meta_topLayer].Height, zero );

    Hypl_Mass = zero;
    for (i = botmLayer; i <= Meta_topLayer; i++)
         Hypl_Mass = Lake[i].LayerVol * Lake[i].Density + Hypl_Mass;

    Dens_Hypl = Hypl_Mass / Lake[Meta_topLayer].Vol1;
    gPrimeTwoLayer = gprime(Dens_Epil, Dens_Hypl);  //# gprime_effective
    if (gPrimeTwoLayer <= 1E-7) gPrimeTwoLayer = 1E-7;
    Vol_Epi = Lake[surfLayer].Vol1 - Lake[Meta_topLayer].Vol1;
    Vol_Hypl = Lake[Meta_topLayer].Vol1;
    Epi_Thick = (Vol_Epi/(Lake[Meta_topLayer].LayerArea+Lake[surfLayer].LayerArea))*2.0;
    Hypl_Thick = (Vol_Hypl/Lake[Meta_topLayer].LayerArea)*2.0;
    IntWaveSpeed = sqrt((fabs(gPrimeTwoLayer)*Epi_Thick*Hypl_Thick)/(Epi_Thick+Hypl_Thick));

    /**************************************************************************
     * Adjust momentum for water entrained by the deepening process           *
     * if Epi_Thick > PrevThick                                               *
     **************************************************************************/
    u_avg_old = u_avg;
    if (u_avg > zero && PrevThick < Epi_Thick) {
         u_f = u_f * PrevThick / Epi_Thick;
         u0 = u_f;
    }
    u0_old = u0;

    /**************************************************************************
     * Compute effective length at the level of the thermocline               *
     * Assume:                                                                *
     *      1) that lake approximates as an ellipse                           *
     *      2) area = pi/4 * Length * Width                                   *
     *      3) ratio Length:Width at thermocline is same as at crest          *
     **************************************************************************/
    LengthAtThermo = sqrt(Lake[Meta_topLayer].LayerArea*4.0/Pi*(LenAtCrest/WidAtCrest));

    /**************************************************************************
     * Check momentum time counters                                           *
     * Half_Seiche_Period is one half the seiche period                       *
     * Half_Seiche_Period > 0 indicates a current shear event                 *
     * EffectiveForceTime is the effective forcing time                       *
     * Time_start_shear is the start of shear forcing (hours from sim start)  *
     * Time_count_end_shear is the end of shear forcing (hours from sim start)*
     * dt_damp is the damping time                                            *
     **************************************************************************/
    if (Half_Seiche_Period <= zero) {
         AED_REAL EffectiveForceTime;
         Half_Seiche_Period = LengthAtThermo / (2.0 * IntWaveSpeed * SecsPerHr);
         EffectiveForceTime = Half_Seiche_Period;
         if (U_star > zero) {
            AED_REAL Lake_Depth = Epi_Thick+Hypl_Thick;
            AED_REAL dt_damp = 2.0*(Lake_Depth/tdfac)*(Lake_Depth/U_star_sqr) *
                        (Hypl_Thick/Epi_Thick) *
                        pow((gPrimeTwoLayer*Epi_Thick*Hypl_Thick/Lake_Depth), 0.25) /
                        (sqrt(2.0*LengthAtThermo)*SecsPerHr)*sqrt(Visc);
            EffectiveForceTime = 1.59 * Half_Seiche_Period;
            if (dt_damp/Half_Seiche_Period < 10.0)
                // Add damping factor to effective time
                EffectiveForceTime = (1.0+ 0.59 *(1.0-(1/cosh(dt_damp/Half_Seiche_Period-1.0))))*Half_Seiche_Period;
         }
         Time_start_shear = Time_count_sim;
         Time_count_end_shear = Time_start_shear + EffectiveForceTime;
    }

    /**************************************************************************
     * Calculate momentum forcing parameters for current time step            *
     * accn is the acceleration (m/s**2) of the mixed layer by wind stress    *
     **************************************************************************/
    accn = U_star_sqr / Epi_Thick;
    FSUM = FSUM + accn;
    Slope = (accn-FO) + OldSlope;
    if (accn <= zero) Slope = zero;
    else if (fabs(Slope/accn) <= 1e-5) Slope = zero;
    if (Slope < zero) Slope = zero;

    /**************************************************************************
     * Check for momentum cutoff within current time step. Calculate time step*
     * for forcing, Time_end_shear, and reset parameters for next time step   *
     **************************************************************************/
    //# Add num_hours to sim time counter (in hours since sim start)
    Time_count_sim += noSecs/SecsPerHr;
    if (Time_count_sim >= Time_count_end_shear) {
         //# Here if cutoff within current time step
         Time_end_shear = Time_count_end_shear - Time_count_sim + noSecs/SecsPerHr;
         OldSlope = accn - (FSUM / (Mixer_Count));
         if (OldSlope < zero) OldSlope = zero;
    } else {
         Time_end_shear = noSecs/SecsPerHr;
         OldSlope = Slope;
    }
    FO = accn;

    //# Momentum update
    if (u0 < 1E-7) u0 = zero;
    if (Slope < 1E-7) Slope = zero;
    u_f = u0 + Slope * Time_end_shear * SecsPerHr;
    u_avgSQ = (u_f*u_f + u_f*u0 + u0*u0) / 3.0;
    if (u_avgSQ < 1E-7) u_avgSQ = 1E-7;
    u_avg = sqrt(u_avgSQ);
    u0 = u_f;
    del_u = u_avg - u_avg_old;

    //# K-H length scales (deltaKH and del_deltaKH) for shear production
    del_deltaKH = 2.0 * coef_mix_KH * u_avg * del_u / gPrimeTwoLayer;
    deltaKH = coef_mix_KH * u_avg * u_avg / gPrimeTwoLayer;
    if (deltaKH < 1.0E-10) deltaKH = 0.0;
    if (del_deltaKH < 1.0E-10) del_deltaKH = 0.0;

    //# Energy from shear-mixing
    Energy_Shear = half * coef_mix_shear * (u_avg * u_avg * (zsml_tilda + del_deltaKH / 6.0) + u_avg * deltaKH * del_u / 3.0) +
           redg * deltaKH * (deltaKH * zsml_tilda / (twfour * (Lake[surfLayer].Height - Lake[Meta_topLayer].Height)) -
           del_deltaKH / twelve);

    if (Energy_Shear < zero) Energy_Shear = zero;

    //# Add available kinetic energy
    Energy_AvailableMix += Energy_Shear;
    GPEFFC = gPrimeTwoLayer * Epi_Thick;

    loop_count = 0;
    _DBG_MIXER_(3, 0, 0,     // step 2, before loop, loop_step
           NumLayers, Epi_botmLayer, Meta_topLayer,
           Energy_AvailableMix, Energy_RequiredMix, redg);

    //# Now undergo deepening loop for shear production
    del_u = zero;
    while (TRUE) {
        //# Save current values of Epi_Thick and u_avg in case of mixing
        PrevThick = Epi_Thick;
        u_eff = u_avg;

        //# Thickness of surface layer and next layer down
        if (Meta_topLayer > botmLayer)
            delzkm1 = Lake[Meta_topLayer].Height - Lake[Meta_topLayer-1].Height;
        else
            delzkm1 = Lake[Meta_topLayer].Height;
        Epi_dz = Lake[surfLayer].Height - Lake[Meta_topLayer].Height;

        //# Compute energy available for mixing next layer, Ea
        Energy_Deepen = half * coef_mix_shear * (u_eff * u_eff * (delzkm1 + del_deltaKH / 6.0) + u_eff * deltaKH * del_u / 3.0)
                      + redg * deltaKH * (deltaKH * delzkm1 / (twfour * Epi_dz) - del_deltaKH / twelve);

        if (Energy_Deepen < zero) Energy_Deepen = zero;

        //# Add available kinetic energy
        Energy_AvailableMix += Energy_Deepen;

        //# Compute energy required to entrain next layer, Er
        Energy_RequiredMix = (redg * Epi_dz + coef_mix_turb * q_sqr) * delzkm1 / 2.0;

        //# Compare energy available with energy required
        dbgprt("Energy_AvailableMix = %10.5f\n",Energy_AvailableMix);
        dbgprt("Energy_RequiredMix = %10.5f\n",Energy_RequiredMix);
        dbgprt( (Energy_AvailableMix < Energy_RequiredMix) ? "About to break\n":"");

        _DBG_MIXER_(3, 1, ++loop_count,     // step 2, in loop, loop_step
           NumLayers, Epi_botmLayer, Meta_topLayer,
           Energy_AvailableMix, Energy_RequiredMix, redg);

        //# If there is insufficient energy to deepen then break
        if (Energy_AvailableMix < Energy_RequiredMix) break;

        //# Entrain layer Meta_topLayer
        dbgprt("NumLayers, surfLayer, botmLayer, Meta_topLayer, Epi_botmLayer = %d,%d,%d,%d,%d\n",
                                      NumLayers, surfLayer, botmLayer, Meta_topLayer, Epi_botmLayer);
        dbgprt("Epilimnion Mass, Temp, Salt, Dens = %10.5f,%10.5f,%10.5f,%10.5f\n",
                                      Vol_Epi,MeanTemp,MeanSalt,Dens_Epil);
        dbgprt("Meta_topLayer Vol, Temp, Salt, Dens = %10.5f,%10.5f,%10.5f,%10.5f\n",
                                      Lake[Meta_topLayer].LayerVol,
                                      Lake[Meta_topLayer].Temp,Lake[Meta_topLayer].Salinity,
                                      Lake[Meta_topLayer].Density);

        add_this_layer(&VMsum,&Tsum,&Ssum,&Mass_Epi,&MeanTemp,&MeanSalt,&Dens_Epil,Meta_topLayer);
        average_layer(&Meta_topLayer,&Epi_botmLayer,MeanTemp,MeanSalt,Dens_Epil);

        dbgprt("New Epilimnion Mass, Temp, Salt, Dens = %10.5f,%10.5f,%10.5f,%10.5f\n",
                                      Vol_Epi,MeanTemp,MeanSalt,Dens_Epil);
        dbgprt("NumLayers, surfLayer, botmLayer, Meta_topLayer, Epi_botmLayer = %d,%d,%d,%d,%d\n",
                                      NumLayers,surfLayer,botmLayer,Meta_topLayer,Epi_botmLayer);

        //# Just used energy to entrain another layer so no longer available
        Energy_AvailableMix -= Energy_RequiredMix;

        if (Meta_topLayer < botmLayer) {
            *_Dens_Epil = Dens_Epil; *_Meta_topLayer = Meta_topLayer;
            return DEEPENED_BOT;
        }

        //# Adjust u_f and u_avg for entrained mass
        Hypl_Mass = zero;

        for (i = botmLayer; i <= Meta_topLayer; i++)
            Hypl_Mass = Lake[i].LayerVol*Lake[i].Density+Hypl_Mass;

        Dens_Hypl = Hypl_Mass/Lake[Meta_topLayer].Vol1;
        redg = gprime(Dens_Epil,Lake[Meta_topLayer].Density);

        Vol_Hypl = Lake[Meta_topLayer].Vol1;
        Vol_Epi = Vol_Epi + Lake[Meta_topLayer+1].LayerVol;
        Epi_Thick = (Vol_Epi / (Lake[Meta_topLayer].LayerArea + Lake[surfLayer].LayerArea)) * 2.0;
        u_f = u_f * PrevThick / Epi_Thick;
        u_avg = sqrt((u0_old*u0_old + u0_old*u_f + u_f*u_f)/3.0);
        del_u = u_avg-u_eff;
        deltaKH = coef_mix_KH*u_avg*u_avg/(redg);
        del_deltaKH = 2.0*coef_mix_KH*u_avg*del_u/(redg);
        u0 = u_f;
    }

    _DBG_MIXER_(3, 2, 0,     // step 2, in loop, loop_step
           NumLayers, Epi_botmLayer, Meta_topLayer,
           Energy_AvailableMix, Energy_RequiredMix, redg);

    //# Here if insufficient energy to entrain next layer

    //# Check momentum time counters for cutoff
    Thermocline_Height = Lake[Meta_topLayer].Height;
    gPrimeTwoLayer = GPEFFC / Epi_Thick;

    //# Now update and average water quality
    for (wqvidx=0; wqvidx < Num_WQ_Vars; wqvidx++) {
         WQ_VarsM[wqvidx] = 0.0;
         // sum over the layers
         for (i = Meta_topLayer+1; i <= surfLayer; i++)
            WQ_VarsM[wqvidx] = _WQ_Vars(wqvidx, i) * Lake[i].LayerVol + WQ_VarsM[wqvidx];
         // divide by the volume
         WQ_VarsM[wqvidx] = WQ_VarsM[wqvidx] / (Lake[surfLayer].Vol1-Lake[Meta_topLayer].Vol1);
    }

    /**************************************************************************
     *                                                                        *
     * STEP 4 - KELVIN-HELMHOLTZ MIXING                                       *
     *                                                                        *
     * Check if the interface is unstable, such that shear would induce       *
     * billows. This is done in the separate routine kelvin_helmholtz()       *
     **************************************************************************/
    if (coef_mix_KH > zero && surface_mixing==1)
         Dens_Epil = kelvin_helmholtz(&Meta_topLayer, &Epi_botmLayer,
                                      Dens_Epil, WQ_VarsM);

    /**************************************************************************
     *                                                                        *
     * FINAL JOBS                                                             *
     *                                                                        *
     **************************************************************************/
    DepMX = Lake[Meta_topLayer].Height;

    *_Dens_Epil = Dens_Epil;
    *_Meta_topLayer = Meta_topLayer;

    //# Check if insufficient energy to mix this time step,
    //# and keep count of mixing model time count
    if (Time_count_sim < Time_count_end_shear)  return IS_MIXED;

    //# All available energy used, reset mixing model step count
    return MOMENTUM_CUT;
}



/******************************************************************************
 * Performs the surface mixing due to wind forcing                            *
 ******************************************************************************/
void do_mixing()
{
#ifndef _VISUAL_C_
    // The visual c compiler doesn't like this so must malloc manually
    AED_REAL WQ_VarsM[Num_WQ_Vars];
#else
    AED_REAL *WQ_VarsM;
#endif

    int i, wqvidx;                //# Counters
    int Meta_topLayer;            //# Index for top layer of hypolimnion
    AED_REAL Dens_Epil;           //# Mean epilimnion density [kg/m3]
    AED_REAL Vol_Hypl;            //# Volume of hypolimnion [m^3]
    static int Mixer_Count = 0;   //# Mixer model step counter
    int res = -1;

/*----------------------------------------------------------------------------*/

#ifdef _VISUAL_C_
    WQ_VarsM = malloc(sizeof(AED_REAL) * Num_WQ_Vars);
#endif

    Mixer_Count++;  //# Increment mixing step counter

    switch ( (res = mixed_layer_deepening(WQ_VarsM, Mixer_Count, &Meta_topLayer, &Dens_Epil)) ) {
        case DEEPENED_BOT:
            //# Here if deepened to bottom
            OldSlope = zero; //# Old slope = zero as fully mixed
            Energy_AvailableMix = zero;   //# Total available energy to mix reset to zero as lake fully mixed
            FO = zero;

            //# mix all water quality variables to average value
            for (wqvidx=0; wqvidx < Num_WQ_Vars; wqvidx++) {
                 WQ_VarsM[wqvidx] = 0.0;
                 for (i = botmLayer; i <= surfLayer; i++)
                     WQ_VarsM[wqvidx] = _WQ_Vars(wqvidx, i) * Lake[i].LayerVol + WQ_VarsM[wqvidx];
                 WQ_VarsM[wqvidx] = WQ_VarsM[wqvidx] / Lake[surfLayer].Vol1;
            }
            /***** fall through ******/

         case MOMENTUM_CUT:
            //# Here if momentum cutoff
            Mixer_Count = 0;  //# Reset mixing model step count (note: not reset if case IS_MIXED)
            FSUM = zero;
            Half_Seiche_Period = zero;  //# Reset one half the seiche period
            u_f = zero;
            u0 = zero;
            u_avg = zero;
            Time_start_shear = Time_count_sim;     //# Reset start time of shear forcing to hours from sim start
            Time_count_end_shear = Time_count_sim; //# Reset finish time of shear forcing to hours from sim start

            //# Mark 2 ends here. At this stage layers Meta_topLayer+1,Meta_topLayer+2,---Epi_botmLayer-1 are also mixed
            //# So make epilimnion into one big layer adding Meta_topLayer+1 ... surfLayer

            //# Renumber mixed layers
            /***** fall through ******/

        case IS_MIXED:
            //# Meta_topLayer+1 becomes the surface layer == mixed epilimnion layers
            Lake[Meta_topLayer+1].Height = Lake[surfLayer].Height;
            Lake[Meta_topLayer+1].Temp = MeanTemp;
            Lake[Meta_topLayer+1].Salinity = MeanSalt;

            //# water quality and particles
            for (wqvidx=0; wqvidx < Num_WQ_Vars; wqvidx++)
                _WQ_Vars(wqvidx,Meta_topLayer+1) = WQ_VarsM[wqvidx];

            //# reset the layer volume, density and area for the surface layer
            Lake[Meta_topLayer+1].Vol1 = Lake[surfLayer].Vol1;

            //# calculate cumulative volume of layers bottom .. Meta_topLayer
            //# technically volume of hypolimnion + metalimnion
            Vol_Hypl = zero; //# The case when fully mixed Meta_topLayer<bottom layer
            if (Meta_topLayer >= botmLayer) Vol_Hypl = Lake[Meta_topLayer].Vol1;

            //# volume of new mixed mega epilimnion layer
            Lake[Meta_topLayer+1].LayerVol = Lake[surfLayer].Vol1 - Vol_Hypl;
            Lake[Meta_topLayer+1].Density = Dens_Epil;
            Lake[Meta_topLayer+1].LayerArea = Lake[surfLayer].LayerArea;

            NumLayers = Meta_topLayer + 2; //# add 2 as count from 0 (ie bottom layer == 0)

            dbgprt("Time_count_sim = %10.5f\n",Time_count_sim);
            dbgprt("Time_start_shear = %10.5f\n",Time_start_shear);
            dbgprt("Time_count_end_shear = %10.5f\n",Time_count_end_shear);
            dbgprt("End mix KH NumLayers, surfLayer, botmLayer, Meta_topLayer, Epi_botmLayer = %d,%d,%d,%d,%d\n",
                                               NumLayers, surfLayer, botmLayer, Meta_topLayer, Epi_botmLayer);

            /***** fall through ******/
        default :
            break;
    }

#ifdef _VISUAL_C_
    free(WQ_VarsM);
#endif
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
static AED_REAL kelvin_helmholtz(int *Meta_topLayer, int *Epi_botmLayer, AED_REAL Dens, AED_REAL *WQ_VarsM)
{
    AED_REAL Surface_Height; //# Height of lake surface
    AED_REAL Delta_Mix;      //# Thickness of mixing layer [m]
    AED_REAL eps = 0.02;     //# Minimum tolerance
    AED_REAL eps6;           //# Six times minimum tolerance (0.12)
    AED_REAL t_billow;       //# Time period for billowing effects
    AED_REAL bilshear;       //# Ratio time period for billowing to shear
    AED_REAL T;
    AED_REAL THB;
    AED_REAL THT;
    AED_REAL DepthB;
    AED_REAL DNL;
    AED_REAL D3;
    AED_REAL eBot;
    AED_REAL eTop;
    AED_REAL HMIN;

    int Meta_botmLayer;      //# Layer index for the layer intersecting the bottom of the metalimnion
    int j1, k1;              //# j1 = top layer of the metalimnion, k1 = bottom layer of the epilimnion
    int i;
    int ir;
    int kl1;
    int lbi;
    int n;
    int nl;
    int top, up;
    int wqvidx;

/*----------------------------------------------------------------------------*/

    j1 = *Meta_topLayer; k1 = *Epi_botmLayer;
    dbgprt("Start KH NumLayers, surfLayer, botmLayer, j1, k1 = %d,%d,%d,%d,%d\n",NumLayers, surfLayer, botmLayer, j1, k1);

    //# Set Tolerances
    eps6 = 6. * eps;

    //# Compute Delta_Mix
    t_billow = u_avg / (gPrimeTwoLayer * SecsPerHr);
    bilshear = t_billow / Time_end_shear;
    Delta_Mix = 0.;

    if (surfLayer <= botmLayer || bilshear > 10.0) return Dens;

    Delta_Mix = (coef_mix_KH*u_avg*u_avg)/(gPrimeTwoLayer*2.0*cosh(bilshear));

    //# Limit the thickness of the mixing layer to less than either the hypolimnion or epilimnion
    HMIN = MIN(Epi_dz, Thermocline_Height);
    if (Delta_Mix > HMIN) Delta_Mix = HMIN;
    eTop = Thermocline_Height + Delta_Mix;
    eBot = Thermocline_Height - Delta_Mix;
    top = FALSE;
    if (eTop > Lake[surfLayer-1].Height) eTop = Lake[surfLayer-1].Height;
    if ((eTop-Thermocline_Height) < eps6/2.0) return Dens;
    if (eBot < Lake[botmLayer].Height) eBot = Lake[botmLayer].Height;
    if ((Thermocline_Height-eBot) < eps6/2.0) return Dens;
    Surface_Height = Lake[surfLayer].Height;

    //# Find layer intersecting ebot
    for (Meta_botmLayer = botmLayer; Meta_botmLayer <= (k1-1); Meta_botmLayer++)
        if (Lake[Meta_botmLayer].Height > eBot) break;

    DepthB = zero;
    if (Meta_botmLayer > botmLayer) DepthB = Lake[Meta_botmLayer-1].Height;

    //# Check to see if ebot coincides with existing depth value
    T = fabs(eBot-DepthB);
    if (T <= eps) {
        eBot = DepthB;
        Meta_botmLayer--;
    } else {
        T = fabs(Lake[Meta_botmLayer].Height - eBot);
        if (T <= eps) {
            eBot = Lake[Meta_botmLayer].Height;
            if (Meta_botmLayer == k1-1) return Dens;
        } else {
            //# Here if new layer must be added below mixed region
            *Epi_botmLayer = (++k1);
            kl1 = k1-Meta_botmLayer-1; /* = number of layers */

            for (i = botmLayer; i < kl1; i++) {
                *Meta_topLayer = (j1 = k1-i-1);
                Lake[j1].Height = Lake[j1-1].Height;
                Lake[j1].Density = Lake[j1-1].Density;
                Lake[j1].Temp = Lake[j1-1].Temp;
                Lake[j1].Salinity = Lake[j1-1].Salinity;

                for (wqvidx=0; wqvidx < Num_WQ_Vars; wqvidx++)
                    _WQ_Vars(wqvidx,j1) = _WQ_Vars(wqvidx,j1-1);
            }
            Lake[Meta_botmLayer].Height = eBot;
        }
    }

    //# Here after position (ebot) of bottom of shear zone has been
    //# determined and extra layer added, if necessary
    T = fabs(eTop-eBot);
    if (T < eps6) return Dens;

    //# Check number of layers in bottom half of shear layer - there must
    //# be at least three
    nl = k1-Meta_botmLayer-1;
    if (nl < 3) {
        if (nl != 2) {
            //# Here if nl=1
            nl = 3;
            D3 = (Lake[k1-1].Height - eBot)/3.0;
            Lake[Meta_botmLayer+3].Height = Lake[k1-1].Height;
            Lake[Meta_botmLayer+2].Height = eBot + D3 + D3;
            Lake[Meta_botmLayer+1].Height = eBot + D3;
            for (i = 2; i <= 3; i++) {
                lbi = Meta_botmLayer + i;
                Lake[lbi].Density = Lake[k1-1].Density;
                Lake[lbi].Temp = Lake[k1-1].Temp;
                Lake[lbi].Salinity = Lake[k1-1].Salinity;
                for (wqvidx=0; wqvidx < Num_WQ_Vars; wqvidx++)
                    _WQ_Vars(wqvidx,lbi)=_WQ_Vars(wqvidx,k1-1);
            }
            *Epi_botmLayer = (k1 = Meta_botmLayer+4);
        } else {
            //# Here if nl=2.0
            nl = 3;

            //# Find thicker layer
            THT = Lake[k1-1].Height - Lake[k1-2].Height;
            THB = Lake[k1-2].Height - eBot;
            if (THT <= THB) {
                //# Here if tht .le. thb - divide layer k1-2
                Lake[Meta_botmLayer+1].Height = eBot + THB/2.0;
                Lake[Meta_botmLayer+2].Height = eBot + THB;
                Lake[Meta_botmLayer+3].Height = eBot + THT + THB;
                *Epi_botmLayer = (++k1);
                Lake[k1-1].Density = Lake[k1-2].Density;
                Lake[k1-1].Temp = Lake[k1-2].Temp;
                Lake[k1-1].Salinity = Lake[k1-2].Salinity;

                for (wqvidx=0; wqvidx < Num_WQ_Vars; wqvidx++)
                    _WQ_Vars(wqvidx,k1-1)=_WQ_Vars(wqvidx,k1-2);

                Lake[k1-2].Density  = Lake[k1-3].Density;
                Lake[k1-2].Temp     = Lake[k1-3].Temp;
                Lake[k1-2].Salinity = Lake[k1-3].Salinity;

                for (wqvidx=0; wqvidx < Num_WQ_Vars; wqvidx++)
                    _WQ_Vars(wqvidx,k1-2) = _WQ_Vars(wqvidx,k1-3);
            } else {
                //# Here if tht > thb - split layer k1-1
                *Epi_botmLayer = (++k1);
                Lake[k1-1].Height = Lake[k1-2].Height;
                Lake[k1-1].Temp = Lake[k1-2].Temp;
                Lake[k1-1].Salinity = Lake[k1-2].Salinity;

                for (wqvidx=0; wqvidx < Num_WQ_Vars; wqvidx++)
                    _WQ_Vars(wqvidx,k1-1) = _WQ_Vars(wqvidx,k1-2);

                Lake[k1-1].Density = Lake[k1-2].Density;
                Lake[k1-2].Height = Lake[k1-3].Height + THT/2.0;
            }
        }
    }

    //# Here after bottom half of shear zone has at least three layers,
    //# divide top half of shear zone into nl layers
    DNL = (eTop - Lake[k1-1].Height)/(nl);
    for (i = 1; i <= nl; i++) {
        *Meta_topLayer = (j1 = i + k1 - 1);
        Lake[j1].Height = Lake[k1-1].Height + (i)*DNL;
        Lake[j1].Density = Dens;
        Lake[j1].Temp = MeanTemp;
        Lake[j1].Salinity = MeanSalt;

        for (wqvidx=0; wqvidx < Num_WQ_Vars; wqvidx++)
            _WQ_Vars(wqvidx,j1) = WQ_VarsM[wqvidx];
    }

    //# The number of the layer just below the mixed region is now j1=k1+nl-1
    //# unless top = T ... if top=T, then layer j1 is the mixed region
    if (!top) {
        NumLayers = k1 + nl + 1;
        Lake[surfLayer].Temp = MeanTemp;
        Lake[surfLayer].Density = Dens;
        Lake[surfLayer].Salinity = MeanSalt;

        for (wqvidx=0; wqvidx < Num_WQ_Vars; wqvidx++)
            _WQ_Vars(wqvidx,surfLayer) = WQ_VarsM[wqvidx];
    } else {
        *Meta_topLayer = (--j1);
        NumLayers = k1+nl;
    }
    Lake[surfLayer].Height = Surface_Height;
    Epi_dz = Lake[surfLayer].Height - Lake[surfLayer-1].Height;

    //# Calculate volumes
    resize_internals(1,botmLayer);

    //# Relax density structure within shear zone
    up = TRUE;
    ir = nl-2;
    while (TRUE) {
        //# Mix middle two layers k1, k1-1
        Lake[k1].Temp = combine(Lake[k1].Temp,   Lake[k1].LayerVol,   Lake[k1].Density,
                                Lake[k1-1].Temp, Lake[k1-1].LayerVol, Lake[k1-1].Density);
        Lake[k1].Salinity = combine(Lake[k1].Salinity,   Lake[k1].LayerVol,   Lake[k1].Density,
                                    Lake[k1-1].Salinity, Lake[k1-1].LayerVol, Lake[k1-1].Density);
        Lake[k1].Density = calculate_density(Lake[k1].Temp,Lake[k1].Salinity);

        for (wqvidx=0; wqvidx < Num_WQ_Vars; wqvidx++)
            _WQ_Vars(wqvidx,k1) = combine_vol(_WQ_Vars(wqvidx,k1),   Lake[k1].LayerVol,
                                              _WQ_Vars(wqvidx,k1-1), Lake[k1-1].LayerVol);

        Lake[k1-1].Density = Lake[k1].Density;
        Lake[k1-1].Temp = Lake[k1].Temp;
        Lake[k1-1].Salinity = Lake[k1].Salinity;

        for (wqvidx=0; wqvidx < Num_WQ_Vars; wqvidx++)
            _WQ_Vars(wqvidx,k1-1)=_WQ_Vars(wqvidx,k1);

        if (ir == botmLayer) break;

        while (TRUE) {
            //# do loop mixes up (up = .true.) or down (up = .false.)
            for (i = botmLayer; i <= ir; i++) {
                if (up) n = k1 + i - 1;
                else    n = k1 - i - 1;
                Lake[n].Temp = combine(Lake[n].Temp,Lake[n].LayerVol,Lake[n].Density,
                                       Lake[n+1].Temp,Lake[n+1].LayerVol,Lake[n+1].Density);
                Lake[n].Salinity = combine(Lake[n].Salinity,Lake[n].LayerVol,
                                           Lake[n].Density,Lake[n+1].Salinity,
                                           Lake[n+1].LayerVol,Lake[n+1].Density);
                Lake[n].Density = calculate_density(Lake[n].Temp,Lake[n].Salinity);

                for (wqvidx=0; wqvidx < Num_WQ_Vars; wqvidx++)
                     _WQ_Vars(wqvidx,n) = combine_vol(_WQ_Vars(wqvidx,n),   Lake[n].LayerVol,
                                             _WQ_Vars(wqvidx,n+1), Lake[n+1].LayerVol);
                Lake[n+1].Density = Lake[n].Density;
                Lake[n+1].Temp = Lake[n].Temp;
                Lake[n+1].Salinity = Lake[n].Salinity;

                for (wqvidx=0; wqvidx < Num_WQ_Vars; wqvidx++)
                     _WQ_Vars(wqvidx,n+1) = _WQ_Vars(wqvidx,n);
            }

            up = !up;
            if (up) {
                ir--;
                break;
            }
        }
    }

    //# Here after relaxation complete, reset the mixed region variables
    Dens = Lake[surfLayer].Density;
    MeanTemp = Lake[surfLayer].Temp;
    MeanSalt = Lake[surfLayer].Salinity;
    for (wqvidx = 0; wqvidx < Num_WQ_Vars; wqvidx++)
        WQ_VarsM[wqvidx] = _WQ_Vars(wqvidx, surfLayer);

//  These 2 values are not used
//  Vol_Epi = Lake[surfLayer].LayerVol;
//  Mass_Epi = Vol_Epi*Dens;

    dbgprt("End KH NumLayers, surfLayer, botmLayer, j1, k1 = %d,%d,%d,%d,%d\n",NumLayers, surfLayer, botmLayer, j1, k1);

    return Dens;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
