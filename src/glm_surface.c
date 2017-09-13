/******************************************************************************
 *                                                                            *
 * glm_surface.c                                                              *
 *                                                                            *
 * Developed by :                                                             *
 *     AquaticEcoDynamics (AED) Group                                         *
 *     School of Earth & Environment                                          *
 *     The University of Western Australia                                    *
 *                                                                            *
 *     http://aed.see.uwa.edu.au/                                             *
 *                                                                            *
 * Copyright 2013 - 2016 -  The University of Western Australia               *
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
#include <math.h>

#include "glm.h"

#include "glm_types.h"
#include "glm_const.h"
#include "glm_globals.h"

#include "glm_util.h"
#include "aed_time.h"
#include "glm_mixu.h"

//#define dbgprt(...) fprintf(stderr, __VA_ARGS__)
#define dbgprt(...) /* __VA_ARGS__ */


AED_REAL  atmos_stability(AED_REAL *Q_latentheat,
                          AED_REAL *Q_sensible,
                          AED_REAL  WindSp,
                          AED_REAL  WaterTemp,
                          AED_REAL  AirTemp,
                          AED_REAL  p_atm,
                          AED_REAL  humidity_surface,
                          AED_REAL  humidity_altitude);

/******************************************************************************
 * Variables for the ice cover components                                     *
 * density_s = density of snow                                                *
 * Tf = temperature at which the water freezes                                *
 * K_snow = thermal conductivity of snow                                      *
 * K_I = thermal conductivity of ice                                          *
 * K_w = molecular thermal conductivity of water                              *
 * Tmelt = temperature of ice melt                                            *
 * L_i = the latent heat of fusion for ice                                    *
 * L_s = the latent heat of fusion for snow                                   *
 * L_fusion = a switch between the snow and ice latent heats of fusion        *
 ******************************************************************************/

int ice = FALSE;          // flag that tells if there is ice formation

// These are made available for the lake.csv output
AED_REAL Q_shortwave;     // Solar radiation at water surface
AED_REAL Q_sensibleheat;  // Sensible heat flux
AED_REAL Q_latentheat;    // Evaporative heat flux
AED_REAL Q_longwave;      // Net longwave heat flux

static AED_REAL  Q_iceout;     // heat flux out of snow/ice surface
static AED_REAL  Q_icewater;   // heat flux across water/ice interface
static AED_REAL  Q_surflayer;  // heat flux through the water surface
static AED_REAL  Q_watermol;   // heat flux through water due to molecular conductivity
static AED_REAL  Q_underflow;  // heat flux through water due to flow under the ice
static AED_REAL  Q_icemet;     // the heat flux at the surface due to meteorological forcing
//static AED_REAL  U_FLOW;     // the velocity of the underflow



void recalc_surface_salt(void);
AED_REAL calculate_qsw(int kDays, int mDays, int iclock,
                        AED_REAL Latitude, AED_REAL SWOld, AED_REAL ShortWave, AED_REAL WindSp);



/******************************************************************************
 *                                                                            *
 *  atm_density  [kg m^3]                                                     *
 *     = 0.348 * (1 + r)/(1 + 1.61*r) * p/T   of TVA, p 5.12                  *
 *                                                                            *
 *  Returns the density of the atmosphere given total atmospheric pressure    *
 *  and water vapour pressure.                                                *
 *                                                                            *
 *  Ref: TVA, [1972,p5.12]                                                    *
 *                                                                            *
 * -- Ratio of the molecular (or molar) weight of water to dry air [-]:       *
 *    mwrw2a    =    18.016    /    28.966;      // = 0.62197                 *
 *                                                                            *
 * -- The universal gas constant  [J mol^-1 K^-1] = 8.31436;                  *
 * -- Gas constant for dry air in terms of mass of gas rather than moles      *
 * -- [J kg^-1 K^-1]:                                                         *
 *    c_gas   = 1.0E3 *    8.31436     /    28.966;                           *
 *                                                                            *
 ******************************************************************************/
const AED_REAL mwrw2a  = 18.016 / 28.966;               // = 0.62197
const AED_REAL c_gas   = 1.0E3 * 8.31436 / 28.966;

/******************************************************************************/
AED_REAL atm_density(AED_REAL atmosPressure, // (total) atmospheric pressure     [Pa]
                     AED_REAL vapPressure,   // water vapour (partial) pressure  [Pa]
                     AED_REAL AirTemp)       // dry bulb air temperature         [Cel]
{
/*
    // Dry air
    return (p_atm)/(287.058 * (AirTemp + Kelvin));
*/

    // Moist air
    // mixing ratio: r = Mwater/(Mwater+Mdry_air)
    AED_REAL r = mwrw2a * vapPressure/(atmosPressure - vapPressure);

    return 1.0/c_gas * (1 + r)/(1 + r/mwrw2a) * atmosPressure/(AirTemp + Kelvin);
}


/******************************************************************************
 * Performs the thermal transfers at the surface                              *
 * ice cover version, Brett Wallace 1996, David Hamilton 1998                 *
 * time step fixed to one hour, David Hamilton, 1998                          *
 * Light modified to sinusoidal distribution, David Hamilton 1993             *
 * The day commences at midday, variable daylength (Time of year and latitude)*
 ******************************************************************************/
void do_surface_thermodynamics(int jday, int iclock, int LWModel,
                              AED_REAL Latitude, AED_REAL SWOld, AED_REAL ShortWave)
// LWModel   // type of longwave radiation
// SWOld     // Total solar radiation at the surface for yesterday
// ShortWave // Total solar radiation at the surface for today
{
/*----------------------------------------------------------------------------*/
    const AED_REAL  K_ice_white = 2.3;  //# thermal conductivity of white ice
    const AED_REAL  K_ice_blue = 2.0;   //# thermal conductivity of blue ice
    const AED_REAL  K_water = 0.57;     //# molecular thermal conductivity of water
//  const AED_REAL  CSEN = 0.0014;      //# sensible heat?  only used in underflow => redundant?
    const AED_REAL  Latent_Heat_Fusion = 334000.0;  // latent heat of fusion J/kg both ice & snow
    const AED_REAL  Temp_melt = 0.0;    //# temperature at which ice melts ~= temperature at which water freezes = 0.0oC
    const AED_REAL  f_sw_wl1 = 0.7;     //# fraction of short wave radiation in first wavelength band
    const AED_REAL  f_sw_wl2 = 0.3;     //# fraction of short wave radiation in second wavelength band
    const AED_REAL  attn_ice_blue_wl1 = 1.5;   //# attenuation coefficient of the ice in the first spectral band
    const AED_REAL  attn_ice_blue_wl2 = 20.;   //# attenuation coefficient of the ice in the second spectral band
    const AED_REAL  attn_ice_white_wl1 = 6.0;  //# attenuation coefficient of the white ice in the first spectral band
    const AED_REAL  attn_ice_white_wl2 = 20.;  //# attenuation coefficient of the white ice in the second spectral band
    const AED_REAL  attn_snow_wl1 = 6.0;  //# attenuation coefficient of the snow in the first spectral band
    const AED_REAL  attn_snow_wl2 = 20.;  //# attenuation coefficient of the snow in the second spectral band
    const AED_REAL  rho_ice_blue = 917.0; //# density of blue ice
    const AED_REAL  rho_ice_white = 890.0; //# density of white ice
    const AED_REAL  eps_water = 0.985;     //# emissivity of the water surface
    const AED_REAL  p_atm = 1013.0;        //# Atmospheric pressure in hectopascals ==101300 Pa

/*----------------------------------------------------------------------------*/

    AED_REAL CloudCover;  //# Cloud cover as a fraction of 1 used to calculate incident long wave radiation
    AED_REAL dTemp;       //# Change in temperature due to atmospheric heat exchange light only when ice on

#ifndef _VISUAL_C_
    // The visual c compiler on doesn't like this so must malloc manually
    AED_REAL LayerThickness[MaxLayers],     //# Layer thickness (m)
             heat[MaxLayers];
#else
    AED_REAL *LayerThickness;
    AED_REAL *heat;
#endif
    AED_REAL dHt_WhiteIce;   //# Change in height of white ice when snow weight exceeds buoyancy potential
    AED_REAL Q_lw_in;        //# Long wave radiation from atmosphere W/m2
    AED_REAL Q_whiteice;     //# heat flux due to snow turning into white ice
    AED_REAL Q_snowice;      //# heat flux due to snow turning into white ice
    AED_REAL Q_Rain;         //# heat flux due to rain either melting ice or freezing into ice
    AED_REAL Q_lw_out;       //# Long wave radiation emission W/m2
    AED_REAL SatVap_surface; //# Saturated vapour pressure at surface layer or top ice layer
    AED_REAL WindSp;         //# Wind speed corrected by multiplicative factor

    AED_REAL AAA, BBB, CCC, DDD, EEE, FFF, GGG;
    AED_REAL Temp_ice;        //# Temperaure of the top layer of the ice cover model oC
    AED_REAL T01_OLD,T01_NEW; //# Not sure how these are used as Temp_ice always == 0 >= melting
    AED_REAL rho_snow;        //# Density of snow layer kg/m3
    AED_REAL rho_snow_old;    //# Density of snow layer compacted by rain kg/m3
    AED_REAL compact_snow;    //# Fraction of snow compacted by rain []
    AED_REAL K_snow;          //# thermal conductivity of snow [W/m2/oC]
    AED_REAL SUMPO4, SUMTP, SUMNO3, SUMNH4, SUMTN, SUMSI;
    AED_REAL Q_latent_ice;    //# Latent heat of melting or freezing of ice [W/m2]
    AED_REAL rho_air;         //# atm_density

   //# New parameters for heat flux estimate
    AED_REAL KSED;
    AED_REAL TYEAR;
    AED_REAL ZSED;

    int i;
//  int underFlow;
    int kDays;

/*----------------------------------------------------------------------------*/

#ifdef _VISUAL_C_
    LayerThickness = malloc(sizeof(AED_REAL) * MaxLayers);
    heat = malloc(sizeof(AED_REAL) * MaxLayers);
#endif

    Q_longwave = 0.;
    SurfData.Evap = 0.;
    memset(heat, 0, sizeof(AED_REAL)*MaxLayers);

    T01_NEW = 0.;
    T01_OLD = 0.;
    Q_latent_ice = 0.0;
    SUMPO4 = 0.0;
    SUMTP = 0.0;
    SUMNO3 = 0.0;
    SUMNH4 = 0.0;
    SUMTN = 0.0;
    SUMSI = 0.0;
    Q_whiteice = 0.0;
    Q_snowice = 0.0;
    Q_Rain = 0.0;
    rho_snow = SurfData.RhoSnow;
//  underFlow = FALSE;

    //# Snow conductivity (eqn. 12 Rogers et al. 1995)
    //# To avoid division by zero set to 0.31 when no snow
    if (SurfData.HeightSnow  ==  0.)
        K_snow = 0.31;
    else
        K_snow = 0.021+(0.0042*rho_snow)+(2.2E-09*pow(rho_snow, 3));

    // Initialize Temp_ice to the Air Temperature
    Temp_ice = MetData.AirTemp;

    if (ice) WindSp = 0.00001;
    else     WindSp = MetData.WindSpeed;

    //# Convert julian days to days since the start of the year

    //# Get shortwave - either from sub-daily data, or approximate
    Q_shortwave = calculate_qsw(day_of_year(jday - 1), day_of_year(jday), iclock, Latitude, SWOld, ShortWave, WindSp);

    // Into layer surfLayer goes qsw(surfLayer)-qsw(surfLayer-1) over the
    // area common to layers surfLayer and surfLayer-1, and qsw(surfLayer)
    // over the rest of area(surfLayer) units of heat[surfLayer] are
    // joules/sec; units of area(surfLayer) are 10**6 m**2
    for (i = (botmLayer+1); i <= surfLayer; i++)
        LayerThickness[i] = Lake[i].Height-Lake[i-1].Height;
    LayerThickness[botmLayer] = Lake[botmLayer].Height;


    if (!ice)
        // Assume PAR only (45% of the incident short wave radiation penetrates beyond the surface layer)
        Lake[surfLayer].Light = 0.45 * Q_shortwave;
    else
        // If there is ice cover the heating due to short wave radiation is attenuated across all three ice layers
        Lake[surfLayer].Light =
             (0.45 * Q_shortwave * ((f_sw_wl1 * exp(-attn_snow_wl1 * SurfData.HeightSnow)) +
                     (f_sw_wl2 * exp(-attn_snow_wl2 * SurfData.HeightSnow)))) *
                    ((f_sw_wl1 * exp(-attn_ice_white_wl1 * SurfData.HeightWhiteIce)) +
                     (f_sw_wl2 * exp(-attn_ice_white_wl2 * SurfData.HeightWhiteIce))) *
                    ((f_sw_wl1 * exp(-attn_ice_blue_wl1  * SurfData.HeightBlackIce)) +
                     (f_sw_wl2 * exp(-attn_ice_blue_wl2  * SurfData.HeightBlackIce)));

    Lake[surfLayer-1].Light = Lake[surfLayer].Light * exp(-Lake[surfLayer].ExtcCoefSW*LayerThickness[surfLayer]);

    // Heating due to PAR
    heat[surfLayer] = ( Lake[surfLayer-1].LayerArea * (Lake[surfLayer].Light - Lake[surfLayer-1].Light) +
                       (Lake[surfLayer].LayerArea - Lake[surfLayer-1].LayerArea) * Lake[surfLayer].Light * 0.1);
    // Heating due to NIR/UV (all absorbed in surf layer only)
    //heat[surfLayer] = heat[surfLayer] + 0.55 * Q_shortwave * Lake[surfLayer].LayerArea;
    heat[surfLayer] = heat[surfLayer] + (1 - 0.45) * Lake[surfLayer].Light / 0.45 * Lake[surfLayer].LayerArea;

    //Total daily short wave radiation (J/day), stored for lake.csv
    SurfData.dailyQsw += Lake[surfLayer].LayerArea * Q_shortwave * noSecs;

    // Surface heat exchange: Units are Joules/m**2/s or W/m**2

    // Get atmospheric pressure and density
    //        atm_density = (p_atm*100.0)/(287.058 * (MetData.AirTemp+Kelvin));
    rho_air = atm_density(p_atm*100.0,MetData.SatVapDef,MetData.AirTemp);

        if (iclock == 0 && ice) {
        AED_REAL BuoyantPotential ;

        // Snow cover as well as ice cover
        if (SurfData.HeightSnow > 0.0) {
            if (MetData.Snow > 0.0 && MetData.Rain >= 0.0) {
//fprintf(stdout, "Snow on snow\n");
                // SNOWFALL ON Snow
                if (MetData.Rain == 0.0) MetData.Rain = MetData.Snow*0.10; //Use 10:1 snow volume to water equivalent
                if (MetData.AirTemp > 0.0)
                    compact_snow = 0.166+0.834*(1.-exp(-1.*MetData.Rain));
                else
                    compact_snow = 0.088+0.912*(1.-exp(-1.*MetData.Rain));

                //Compact snow first
                rho_snow_old = rho_snow+(snow_rho_max-rho_snow)*compact_snow;
                SurfData.HeightSnow = SurfData.HeightSnow*rho_snow/rho_snow_old;
                //Determine avg density of combined snow
                rho_snow = 1000.0*((rho_snow_old*SurfData.HeightSnow/1000.0)+MetData.Rain)/(SurfData.HeightSnow+(MetData.Snow));

                SurfData.HeightSnow = SurfData.HeightSnow+(MetData.Snow);

                //This should be a good estimate of water equivalent in from snow
                // Use rain because it should be volumetric equivalent of water (converted above, or supplied by user)
                SurfData.dailySnow += MetData.Rain * Lake[surfLayer].LayerArea;

                if (rho_snow > snow_rho_max) rho_snow = snow_rho_max;
                if (rho_snow < snow_rho_min) rho_snow = snow_rho_min;
                Q_Rain = 0.0;
            } else if (MetData.Snow == 0.0 && MetData.Rain == 0.0) {
                // No snowfall or rainfall, compaction only of snow
                if (MetData.AirTemp > 0.0) compact_snow = 0.166;
                else                        compact_snow = 0.088;

                rho_snow_old = rho_snow+(snow_rho_max-rho_snow)*compact_snow;
                SurfData.HeightSnow = SurfData.HeightSnow*rho_snow/rho_snow_old;
                rho_snow = rho_snow_old;
                Q_Rain = 0.0;
            } else if (MetData.Snow == 0.0 && MetData.Rain > 0.0) {
                /********************************************************************
                 * Rainfall on Snow. Check the air temperature. If AirTemp > 0 then *
                 * add the Rain. If AirTemp < 0 then add the rainfall to Snow.      *
                 ********************************************************************/
                if (MetData.AirTemp > 0.0) {
                    compact_snow = 0.166+0.834*(1.-exp(-1.*MetData.Rain));
                    rho_snow_old = rho_snow+(snow_rho_max-rho_snow)*compact_snow;
                    SurfData.HeightSnow = SurfData.HeightSnow*rho_snow/rho_snow_old;
                    rho_snow = rho_snow_old;
                    Lake[surfLayer].Height = Lake[surfLayer].Height+MetData.Rain;
                    SurfData.dailyRain += MetData.Rain * Lake[surfLayer].LayerArea;

                    recalc_surface_salt();

                    if (Temp_ice == Temp_melt)
                        Q_Rain = SPHEAT*(MetData.AirTemp-Temp_ice)*(MetData.Rain)/noSecs;
                    else
                        Q_Rain = Latent_Heat_Fusion*MetData.Rain/noSecs;
                } else {
                    MetData.Snow = MetData.Rain/snow_rho_max;
                    rho_snow_old = snow_rho_max;
                    SurfData.HeightSnow = SurfData.HeightSnow*rho_snow/rho_snow_old;
                    rho_snow = snow_rho_max;
                    SurfData.HeightSnow = SurfData.HeightSnow+MetData.Snow;
                    Q_Rain = 0.0;
                    SurfData.dailyRain += MetData.Rain * Lake[surfLayer].LayerArea;
                }
            }
        } else {  // no snow
            if (MetData.Snow > 0.0 && MetData.Rain >= 0.0) {
//fprintf(stdout,"Snow add on ice %f\n", SurfData.HeightSnow);
                // Snowfall on ice
                if (MetData.Rain == 0.0)MetData.Rain = MetData.Snow*0.10;
                SurfData.HeightSnow = MetData.Snow;


                rho_snow = rho0*MetData.Rain/MetData.Snow;
                SurfData.dailySnow += MetData.Rain * Lake[surfLayer].LayerArea;

                if (rho_snow > snow_rho_max)rho_snow = snow_rho_max;
                if (rho_snow < snow_rho_min)rho_snow = snow_rho_min;
                Q_Rain = 0.0;
            } else if (MetData.Snow == 0.0 && MetData.Rain == 0.0) {
                // No snowfall or rainfall
                Q_Rain = 0.0;
            } else if (MetData.Snow == 0.0 && MetData.Rain > 0.0) {
                // Rainfall on ice - need to know whether it will contribute to water (T
                if (MetData.AirTemp > 0.0) {
                    Lake[surfLayer].Height = Lake[surfLayer].Height+MetData.Rain;
                    SurfData.dailyRain += MetData.Rain * Lake[surfLayer].LayerArea;
                    recalc_surface_salt();

                    if (Temp_ice == Temp_melt)
                        Q_Rain = SPHEAT*(MetData.AirTemp-Temp_ice)*(MetData.Rain)/noSecs;
                    else
                        Q_Rain = Latent_Heat_Fusion*MetData.Rain/noSecs;
                } else {
                    MetData.Snow = MetData.Rain/snow_rho_max;
                    rho_snow = snow_rho_max;
                    SurfData.HeightSnow = MetData.Snow;
                    SurfData.dailySnow += MetData.Snow * Lake[surfLayer].LayerArea * rho_snow/1000.0;
                    Q_Rain = 0.0;
                }
            }
        }

        // Archmides principle for weight of Snow on Ice:
        // When the weight of ice and snow exceeds the buoyancy of the ice cover,
        // the ice will crack, surface water will seep into snow layer leading to
        // formation of White Ice
        if (rho_snow == 0.) rho_snow = 0.00000000000000001; //CAB# fudge factor!

        BuoyantPotential = ((SurfData.HeightBlackIce*(rho0-rho_ice_blue) +
                                    SurfData.HeightWhiteIce*(rho0-rho_ice_white))/rho_snow);

        if (SurfData.HeightSnow > BuoyantPotential) {

            dHt_WhiteIce = SurfData.HeightSnow-BuoyantPotential;

            Q_whiteice = ((Lake[surfLayer].Temp*SPHEAT+Latent_Heat_Fusion)*Lake[surfLayer].Density
                 *dHt_WhiteIce*(1.-(rho_snow/Lake[surfLayer].Density)))/SecsPerDay;

            SurfData.HeightWhiteIce = SurfData.HeightWhiteIce + SurfData.HeightSnow - BuoyantPotential;
            //Adjust surface layer height down based on water moving into white ice
            Lake[surfLayer].Height -= dHt_WhiteIce * (rho_ice_white - rho_snow) / Lake[surfLayer].Density;

//fprintf(stdout,"Reduced surface layer thickness by %f to %f\n", dHt_WhiteIce * (rho_ice_white - rho_snow) / Lake[surfLayer].Density, Lake[surfLayer].Height);
            Q_snowice = (Q_whiteice*SurfData.HeightSnow)/(2.0*K_snow);
            SurfData.HeightSnow = BuoyantPotential;

            recalc_surface_salt();
            //Missing: Q_whiteice and Q_snowice go unused, probably missing that part of energy budget
            //LCB: Actually used above in heat budget for ice model above only calculated daily
            //      need to make sure Q_whiteice not reset

        } else {
            dHt_WhiteIce = 0.0;
            Q_whiteice = 0.0;
            Q_snowice = 0.0;
        }

        // Rainfall/snowfall and dry deposition effects on nutrients
        SUMPO4 = SUMPO4 + 0.145 + MetData.RainConcPO4 * MetData.Rain ;
        SUMTP  = SUMTP  + 0.368 + MetData.RainConcTp  * MetData.Rain ;
        SUMNO3 = SUMNO3 + 0.507 + MetData.RainConcNO3 * MetData.Rain ;
        SUMNH4 = SUMNH4 + 0.507 + MetData.RainConcNH4 * MetData.Rain ;
        SUMTN  = SUMTN  + 2.576 + MetData.RainConcTn  * MetData.Rain ;
        SUMSI  = SUMSI  + 0.0   + MetData.RainConcSi  * MetData.Rain ;
    } // end iclock == 0 && ice

    if (!ice) {
        //#  Evaporative heat flux affects top layer only
        SatVap_surface = saturated_vapour(Lake[surfLayer].Temp);
        //# Q_latentheat [W/m2] = CE * rho_air * latent heat * psychro const / air_presssure * windspeed * VPD
        Q_latentheat = -CE * rho_air * Latent_Heat_Evap * (0.622/p_atm) * WindSp * (SatVap_surface - MetData.SatVapDef);

        if (Q_latentheat > 0.0) Q_latentheat = 0.0;

        // Evaporative flux in m/s
        if ( no_evap )
            SurfData.Evap = 0.0;
        else
            SurfData.Evap = Q_latentheat / Latent_Heat_Evap / rho0;

        // Conductive Heat Gain only affects top layer.
        //# Q_sensibleheat [W/m2] = CH * rho_air * specific heat * windspeed * temp diff
        Q_sensibleheat = -CH * (rho_air * 1005.) * WindSp * (Lake[surfLayer].Temp - MetData.AirTemp);

        // If chosen by user, do atmospheric stability correction routine
//fprintf(stderr, "coef_wind_drag = %e Q_l %e Q_s %e\n", coef_wind_drag, Q_latentheat, Q_sensibleheat);
        if (atm_stab)
            coef_wind_drag =
               atmos_stability(&Q_latentheat,
                               &Q_sensibleheat,
                                WindSp,
                                Lake[surfLayer].Temp,
                                MetData.AirTemp,
                                p_atm*100.,
                                SatVap_surface,
                                MetData.SatVapDef);
//fprintf(stderr, " and now = %e Q_l %e Q_s %e\n", coef_wind_drag, Q_latentheat, Q_sensibleheat);

        // Long Wave emission (ie. longwave out) affects only top layer.
        Q_lw_out = -Stefan_Boltzman * eps_water * pow((Kelvin+Lake[surfLayer].Temp), 4.0);

        // Long Wave absorption (ie. longwave in) affects only top layer
        // See Henderson-Sellers 1986 for a good summary
        if (LWModel  ==  LW_CC){
            // Cloud data is available
            AED_REAL eps_star = 0.8;  // default in case of a duff cloudmode value
            CloudCover = MetData.LongWave;
            switch (cloud_mode) {
                case 1:
                    // Idso and Jackson (1969)
                    // eps_star = (1.0 + 0.275*CloudCover)*(1.0 - 0.261 * exp(-0.000777 * pow(-MetData.AirTemp, 2.0))); //
                    eps_star = (1.0 + 0.17 * CloudCover * CloudCover)*(1.0 - 0.261 * exp(-0.000777 * pow(-MetData.AirTemp, 2.0)));
                    break;
                case  2:
                    // Swinbank (1963)
                    eps_star = (1.0 + 0.17 * CloudCover * CloudCover) * (9.365e-6 * pow(MetData.AirTemp + Kelvin, 2.0));
                    break;
                case 3:
                    // Brutsaert (1975)
                    eps_star = (1.0 + 0.275*CloudCover) * 0.642 * pow(MetData.SatVapDef/(MetData.AirTemp+Kelvin), 1/7) ;
                    break;
                case 4:
                    // Yajima 2014 - Tono Dam
                    eps_star = (1.0 - pow(CloudCover, 2.796) ) * 0.642 * pow(MetData.SatVapDef/(MetData.AirTemp+Kelvin), 1/7) +
                                                           0.955 * pow(CloudCover, 2.796) ;
                    break;
            }

            Q_lw_in = (1 - 0.03) * eps_star * Stefan_Boltzman * pow((Kelvin+MetData.AirTemp), 4.0);
            Q_longwave = Q_lw_out + Q_lw_in;
        } else if (LWModel  ==  LW_IN) {
            Q_lw_in = (1 - 0.03) * MetData.LongWave;
            Q_longwave = Q_lw_out+Q_lw_in;
        } else if (LWModel  ==  LW_NET)
            Q_longwave = MetData.LongWave;

//fprintf(stderr, "%e %e \n%e %e \n%e\n", Q_latentheat, Q_sensibleheat, Q_longwave, Lake[surfLayer].LayerArea);
        heat[surfLayer] = heat[surfLayer]+(Q_latentheat+Q_sensibleheat+Q_longwave)*Lake[surfLayer].LayerArea;

        // Daily heat budget (J/day)
        SurfData.dailyQe += Q_latentheat * Lake[surfLayer].LayerArea * noSecs;
        SurfData.dailyQh += Q_sensibleheat * Lake[surfLayer].LayerArea * noSecs;
        SurfData.dailyQlw += Q_longwave * Lake[surfLayer].LayerArea * noSecs;
    } else {
        // The various atmospheric fluxes - evaporative, sensible heat, longwave
        // for ice cover all depend on the ice surface temperature.  Note that wind
        // here is not set to zero (USE MetData.WindSpeed, NOT WindSp].
        // Emissivity for LW OUT is 0.985 (eps_water)
        T01_NEW = 50.;
        T01_OLD = -50.;
        Temp_ice = 0.;
        while (1) {
            //Saturated vapor pressure above snow and ice is different than above water
            //From Jeong S (2009), ultimately from Mellor (1964) "Properties of Snow"
            SatVap_surface = (1+(0.00972*Temp_ice)+(0.000042*pow(Temp_ice, 2)))*saturated_vapour(Temp_ice);

            //I think this might be wrong, resulting value seems way too small, even for ice
            //Q_latentheat = -3.9 * MetData.WindSpeed * (SatVap_surface - MetData.SatVapDef);
            //# Q_latentheat [W/m2] = CE * rho_air * latent heat * psychro const / air_presssure * windspeed * VPD
            Q_latentheat = -CE * rho_air * Latent_Heat_Evap * (0.622/p_atm) * MetData.WindSpeed * (SatVap_surface - MetData.SatVapDef);
            if (Q_latentheat > 0.0) Q_latentheat = 0.0;

            //LCB: changed sign of Q_sensible heat to match that of no ice
            //Q_sensibleheat = -CH * MetData.WindSpeed * (Temp_ice - MetData.AirTemp);
            //# Q_sensibleheat [W/m2] = CH * rho_air * specific heat * windspeed * temp diff
            Q_sensibleheat = -CH * (rho_air * 1005.) * WindSp * (Temp_ice - MetData.AirTemp);

            Q_lw_out = -Stefan_Boltzman * eps_water * pow((Kelvin+Temp_ice), 4.0);

            if (LWModel  ==  LW_CC) {
                CloudCover = MetData.LongWave;
                Q_lw_in = (1.0 + 0.275 * CloudCover) * Stefan_Boltzman * pow((Kelvin+MetData.AirTemp), 4.0) *
                      (1.0 - 0.261 * exp(-0.000777E-4 * pow(MetData.AirTemp, 2.0)));
                Q_longwave = Q_lw_out+Q_lw_in;
            } else if (LWModel  ==  LW_IN) {
                Q_lw_in = MetData.LongWave;
                Q_longwave = Q_lw_out+Q_lw_in;
            } else if (LWModel  ==  LW_NET)
                Q_longwave = MetData.LongWave;

            // This is the net meteorological flux that drives the ice algorithm
            // flux due to precipitation on ice or Snow
            Q_icemet = (Q_latentheat+Q_sensibleheat+Q_longwave+Q_Rain);

            // Now determine the new ice/Snow surface temperature based on the balance
            // fluxes h_ice and the expression for the upward heat flux as given by p

            AAA = (1. - exp(-attn_snow_wl1 * SurfData.HeightSnow)) / (K_snow * attn_snow_wl1);
            BBB = exp(-attn_snow_wl1 * SurfData.HeightSnow) * (1.-exp(-attn_ice_white_wl1 * SurfData.HeightWhiteIce))/
                     (attn_ice_white_wl1 * K_ice_blue);
            FFF = exp((-attn_snow_wl1*SurfData.HeightSnow)-(attn_ice_white_wl1 * SurfData.HeightWhiteIce))*
                     (1.-exp(-attn_ice_blue_wl1*SurfData.HeightBlackIce))/(attn_ice_blue_wl1*K_ice_white);
            CCC = (1.-exp(-attn_snow_wl2*SurfData.HeightSnow))/(K_snow*attn_snow_wl2);
            DDD = exp(-attn_snow_wl2*SurfData.HeightSnow)*(1.-exp(-attn_ice_white_wl2*SurfData.HeightWhiteIce))/
                     (attn_ice_white_wl2*K_ice_blue);
            GGG = exp((-attn_snow_wl2*SurfData.HeightSnow)-(attn_ice_white_wl2*SurfData.HeightWhiteIce))*
                     (1.-exp(-attn_ice_blue_wl2*SurfData.HeightBlackIce))/(attn_ice_blue_wl2*K_ice_white);
            EEE = (K_snow*K_ice_white*K_ice_blue)/((SurfData.HeightSnow*K_ice_white*K_ice_blue)+(SurfData.HeightBlackIce*K_snow*K_ice_blue)+
                     (SurfData.HeightWhiteIce*K_snow*K_ice_white));

            Q_iceout = ((Temp_melt-Temp_ice-(Q_shortwave*(f_sw_wl1*(AAA+BBB+FFF)+f_sw_wl2*(CCC+DDD+GGG)))+Q_snowice)*EEE)+Q_shortwave+Q_whiteice;

            // Now compare the balance between the surface fluxes - and iterate
            // by bisection - NOTE the loop back to 200 in the iteration
            //if (fabs(Q_iceout+Q_icemet) > 1.0 && fabs(T01_NEW-T01_OLD) > 0.001) {
            //    if ((Q_iceout+Q_icemet) < 0.0) T01_NEW = Temp_ice;
            //    if ((Q_iceout+Q_icemet) > 0.0) T01_OLD = Temp_ice;
            //    Temp_ice = (T01_NEW+T01_OLD)/2.0;
            //} else {
            //    break;
            //}
            if (fabs(Q_iceout+Q_icemet) > 1.0 && fabs(T01_NEW-T01_OLD) > 0.001) {
                if ((Q_iceout+Q_icemet) < 0.0) T01_NEW = Temp_ice;
                if ((Q_iceout+Q_icemet) > 0.0) T01_OLD = Temp_ice;
                Temp_ice = (T01_NEW+T01_OLD)/2.0;
            } else {
                break; //Start melting
            }
        }// end while

        T01_NEW = 50.0;
        T01_OLD = -50.0;

        //Evaporation in water equivalent (so use rho0)
        SurfData.Evap = Q_latentheat / Latent_Heat_Evap / rho0;

        SurfData.dailyEvap += (SurfData.Evap * noSecs * Lake[surfLayer].LayerArea);

        if (SurfData.HeightSnow > 0.){
            SurfData.HeightSnow += Q_latentheat/Latent_Heat_Evap/rho_snow*noSecs;
        } else if (SurfData.HeightWhiteIce > 0.) {
            SurfData.HeightWhiteIce += Q_latentheat/Latent_Heat_Evap/rho_ice_white*noSecs;
        } else {
            SurfData.HeightBlackIce += Q_latentheat/Latent_Heat_Evap/rho_ice_blue*noSecs;
        }

        //--------------------------------------------------------------------+
        // Now compare the ice/Snow surface temperature with the melting      |
        // temperature - if it is above the melting temp. then adjust         |
        // the thickness of the surface ice/Snow layer accordingly            |
        //--------------------------------------------------------------------+
        if (Temp_ice >= Temp_melt) {
            Temp_ice = Temp_melt;
            SatVap_surface = (1+(0.00972*Temp_ice)+(0.000042*pow((Temp_ice), 2)))*saturated_vapour(Temp_ice);
            //# Q_latentheat [W/m2] = CE * rho_air * latent heat * psychro const / air_presssure * windspeed * VPD
            Q_latentheat = -CE * rho_air * Latent_Heat_Evap * (0.622/p_atm) * WindSp * (SatVap_surface - MetData.SatVapDef);
            if (Q_latentheat > 0.0) Q_latentheat = 0.0;

            //LCB: changed sign of Q_sensible heat to match that of no ice
            //Q_sensibleheat = -CH * MetData.WindSpeed * (Temp_ice - MetData.AirTemp);
            //# Q_sensibleheat [W/m2] = CH * rho_air * specific heat * windspeed * temp diff
            Q_sensibleheat = -CH * (rho_air * 1005.) * WindSp * (Temp_ice - MetData.AirTemp);

            Q_lw_out = -Stefan_Boltzman * eps_water * pow((Kelvin+Temp_ice), 4.0);


            if (LWModel  ==  LW_CC) {
                CloudCover = MetData.LongWave;
                Q_lw_in = (1.0 + 0.275 * CloudCover) * Stefan_Boltzman * pow((Kelvin+MetData.AirTemp), 4.0) *
                      (1.0 - 0.261 * exp(-0.000777E-4 * pow(MetData.AirTemp, 2.0)));
                Q_longwave = Q_lw_out + Q_lw_in;
            } else if (LWModel  ==  LW_IN) {
                Q_lw_in = MetData.LongWave;
                Q_longwave = Q_lw_out + Q_lw_in;
            } else if (LWModel  ==  LW_NET)
                Q_longwave = MetData.LongWave;

            Q_icemet = (Q_latentheat+Q_sensibleheat+Q_longwave+Q_Rain);
            Q_iceout = ((Temp_melt-Temp_ice-(Q_shortwave*(f_sw_wl1*(AAA+BBB+FFF)+f_sw_wl2*(CCC+DDD+GGG))) +Q_snowice)*EEE)+Q_shortwave+Q_whiteice;

            //------------------------------------------------------------------
            // Now determine the new ice/Snow/water surface temperature and height ad
            // based on the balance between the surface fluxes h_ice and the expression
            // the upward heat flux as given by Patterson and Hamblin (1988]. note on
            // melting can occur at the surface and therefore water level will increase
            // note assumption that ice thickness won't change more than 1 cm in a ti
            //------------------------------------------------------------------
            // LAW: reformatted this to read easier. Also fixed usage of densities for
            //      for melting ice

            // If there is snow, melt that first
            if (SurfData.HeightSnow > 0.0) {
            //if (SurfData.HeightSnow != 0.0) {   //!MH Getting -ve Snow Height in lake.csv

                if (rho_snow == 0.0) rho_snow = snow_rho_max;
                SurfData.dHt = (1/(Latent_Heat_Fusion*rho_snow))*(Q_icemet+Q_iceout)*noSecs;
                if ((SurfData.HeightSnow-SurfData.dHt) < 0.0)SurfData.dHt = SurfData.HeightSnow;
                SurfData.HeightSnow = SurfData.HeightSnow-SurfData.dHt;
                Lake[surfLayer].Height = Lake[surfLayer].Height+SurfData.dHt*(rho_snow/Lake[surfLayer].Density);
                recalc_surface_salt();

            // Otherwise melt the white ice
            } else if (SurfData.HeightWhiteIce > 0.){

                SurfData.dHt = (1/(Latent_Heat_Fusion*rho_ice_white))*(Q_icemet+Q_iceout)*noSecs;

                if ((SurfData.HeightWhiteIce-SurfData.dHt) < 0.0)
                    SurfData.dHt = SurfData.HeightWhiteIce;

                SurfData.HeightWhiteIce -= SurfData.dHt;

                Lake[surfLayer].Height += SurfData.dHt*(rho_ice_white/Lake[surfLayer].Density);
                recalc_surface_salt();
            // Lastly, melt the blue ice
            } else {

                SurfData.dHt = (1/(Latent_Heat_Fusion*rho_ice_blue))*(Q_icemet+Q_iceout)*noSecs;
                if ((SurfData.HeightBlackIce-SurfData.dHt) < 0.) {
                    SurfData.dHt = SurfData.HeightBlackIce;
                }
                SurfData.HeightBlackIce = SurfData.HeightBlackIce-SurfData.dHt;

                Lake[surfLayer].Height += SurfData.dHt*(rho_ice_blue/Lake[surfLayer].Density);
                recalc_surface_salt();
            }// End melt snow/white/blueice if

        }// End Melting If

       // Daily heat budget (MJ/day)
       SurfData.dailyQe += Q_latentheat * Lake[surfLayer].LayerArea * noSecs;
       SurfData.dailyQh += Q_sensibleheat * Lake[surfLayer].LayerArea * noSecs;
       SurfData.dailyQlw += Q_longwave * Lake[surfLayer].LayerArea * noSecs;

    }

    // Now look at the ice or water interface
    if (surfLayer > botmLayer) {
        for (i = surfLayer-1; i >= botmLayer; i-- )
            Lake[i].Light = Lake[i+1].Light * exp(-Lake[i+1].ExtcCoefSW*LayerThickness[i+1]);

        /*--------------------------------------------------------------------*
         * Into layer i goes QSW[i]-QSW(i-1) over the area common to layers   *
         * i and i-1 and QSW[i] over the rest of AREA[i]                      *
         * units of heat[i) are joules/sec; units of AREA[i] are 10**6 m**2   *
         *--------------------------------------------------------------------*/
        for (i = surfLayer-1; i >= (botmLayer+1); i--)
            heat[i] = (Lake[i-1].LayerArea * (Lake[i].Light - Lake[i-1].Light) +
                       Lake[i].Light * (Lake[i].LayerArea - Lake[i-1].LayerArea) * 0.1) ;
        heat[botmLayer] = Lake[botmLayer].Light * Lake[botmLayer].LayerArea;
    }

    // Compute the temperature increase in non-surface layers over noSecs
    for (i = botmLayer; i <= surfLayer; i++) {
        if (fabs(heat[i]) >= 1E-20 && Lake[i].Density != 0.0 && Lake[i].LayerVol != 0.0)
            dTemp = heat[i]*noSecs/(SPHEAT*Lake[i].Density*Lake[i].LayerVol);
        else
            dTemp = 0.;

        Lake[i].Temp += dTemp;
    }

//  if (ice) { printf("Surf Temp post light = %10.5f\n",Lake[surfLayer].Temp);}

    // The change in ice thickness at the bottom can now be determined
    // with the temperature of the lake water readjusted for the surface
    // heat exchange influence of the Ice/Snow
    if (ice) {
        // Both ablation and accretion of the ice can occur at the ice-water interface
        // change dht will be governed by the flux coming down from the surface
        // The water below the ice , See Patterson + Hamblin (1988)
        // first calculate the flux through the ice
        AAA = (1. - exp(-attn_snow_wl1 * SurfData.HeightSnow)) / (K_snow * attn_snow_wl1);
        BBB = exp(-attn_snow_wl1*SurfData.HeightSnow)*(1.-exp(-attn_ice_white_wl1*SurfData.HeightWhiteIce))/ (attn_ice_white_wl1*K_ice_blue);
        FFF = exp((-attn_snow_wl1*SurfData.HeightSnow)-(attn_ice_white_wl1*SurfData.HeightWhiteIce))*
                                      (1.-exp(-attn_ice_blue_wl1*SurfData.HeightBlackIce))/(attn_ice_blue_wl1*K_ice_white);
        CCC = (1.-exp(-attn_snow_wl2*SurfData.HeightSnow))/(K_snow*attn_snow_wl2);
        DDD = exp(-attn_snow_wl2*SurfData.HeightSnow)*(1.-exp(-attn_ice_white_wl2*SurfData.HeightWhiteIce))/ (attn_ice_white_wl2*K_ice_blue);
        GGG = exp((-attn_snow_wl2*SurfData.HeightSnow)-(attn_ice_white_wl2*SurfData.HeightWhiteIce))*
                                      (1.-exp(-attn_ice_blue_wl2*SurfData.HeightBlackIce))/(attn_ice_blue_wl2*K_ice_white);
        EEE = (K_snow*K_ice_white*K_ice_blue)/((SurfData.HeightSnow*K_ice_white*K_ice_blue)+(SurfData.HeightBlackIce*K_snow*K_ice_blue)+
                                       (SurfData.HeightWhiteIce*K_snow*K_ice_white));
        Q_iceout = ((Temp_melt-Temp_ice-(Q_shortwave*(f_sw_wl1*(AAA+BBB+FFF)+f_sw_wl2*(CCC+DDD+GGG))) +Q_snowice)*EEE)+Q_shortwave+Q_whiteice;
        Q_icewater = Q_iceout-(Q_shortwave*f_sw_wl1*(1.-exp(-(attn_snow_wl1*SurfData.HeightSnow+attn_ice_blue_wl1*SurfData.HeightBlackIce
                            +attn_ice_white_wl1*SurfData.HeightWhiteIce))))-(Q_shortwave*f_sw_wl2*(1.-exp(-(attn_snow_wl2*
                             SurfData.HeightSnow+attn_ice_blue_wl2*SurfData.HeightBlackIce+attn_ice_white_wl2*SurfData.HeightWhiteIce))));

        // Now determine the flux through the lake water below the ice
        // for temperature flux at the ice water interface use an exponential dec
        // fickian diffusion so use a gaussian distribution and adjust the diffuse
        // Accordingly - see farmer (1978) from p and h (1988)
        // see eq 22 of rogers et al and note that DZ has been adjusted to ...
        Q_watermol = -K_water*(Temp_melt-Lake[surfLayer].Temp)/0.039;

        // Heat transfer due to underflow seems not to be used as underFlow == .FALSE.
        //if (underFlow)
        //    Q_underflow = CSEN*(Lake[surfLayer].Density - rho0)*SPHEAT*U_FLOW*(Lake[surfLayer].Temp-Temp_melt);
        //else
            Q_underflow = 0.0;

        //# LCB: Need to check as Q_latent_ice == 0.0
        Q_surflayer = Q_watermol + Q_underflow + Q_latent_ice;

        // Now we can determine the amount of ablation or accretion of ice at
        // the given by qf-qw.  once the ice has been melted or formed assume
        // fluxes are in equilibrium - thus determining the heat loss/gain from
        // the surf calorimetry].  correction for thermal contraction given by
        // ratios of d is constant and area can be considered as constant as
        // change in depth again because time step is small, error for complete
        // melting of ice wi

        SurfData.dHt = (Q_icewater-Q_surflayer)*noSecs/(Latent_Heat_Fusion*rho_ice_blue);
        if (SurfData.dHt < -1.*SurfData.HeightBlackIce) {
            if (SurfData.dHt < -1. * (SurfData.HeightBlackIce+SurfData.HeightWhiteIce) )
                SurfData.dHt = -1.*(SurfData.HeightBlackIce+SurfData.HeightWhiteIce);
            SurfData.HeightWhiteIce = SurfData.HeightWhiteIce+SurfData.HeightBlackIce+SurfData.dHt;
            SurfData.HeightBlackIce = 0.0;
        } else
            SurfData.HeightBlackIce = SurfData.HeightBlackIce+SurfData.dHt;

        // Adjust water temperature for ice/water exchange, determine latent heat
        // LCB: moved calculation of latent heat of ice above Temperature calculation
        Q_latent_ice = Latent_Heat_Fusion*rho_ice_blue*SurfData.dHt/noSecs;
        //Lake[surfLayer].Temp = Lake[surfLayer].Temp+(((K_water*((Temp_melt-Lake[surfLayer].Temp)/0.039))*Lake[surfLayer].LayerArea*
        //                         noSecs)+Q_latent_ice)/(SPHEAT*Lake[surfLayer].Density*Lake[surfLayer].LayerVol);
        Lake[surfLayer].Temp = Lake[surfLayer].Temp+((-Q_watermol + Q_latent_ice)*Lake[surfLayer].LayerArea*noSecs)
                                 /(SPHEAT*Lake[surfLayer].Density*Lake[surfLayer].LayerVol);
//      printf("Surf Temp = %10.5f\n",Lake[surfLayer].Temp);
        Lake[surfLayer].Height = Lake[surfLayer].Height-SurfData.dHt*(rho_ice_blue/Lake[surfLayer].Density);

        recalc_surface_salt();
    }

    /******************************************************************************
     * Default sediment "heating" factor in glm_surface.c but since it is based   *
     * on Mendota ends up cooling Kinneret! Beware when using default values.     *
     * LAW: Added a switch and parameterization so we can experiment with this.   *
     ******************************************************************************/

    if(sed_heat_sw){
        // Input of heat from the sediments - based on rogers et al heat flux
        // sediment temperature measurements for lake mendota (birge et al 1927)
        // and basically invariant at 5m at deep station, LayerThickness is required
        // LAW: Modified to use cosine so user and specify peak day coefficient
        kDays = day_of_year(jday);
        TYEAR = sed_temp_mean + sed_temp_amplitude * cos(((kDays-sed_temp_peak_doy)*2.*Pi)/365.);
        ZSED = 6.;
        KSED = 1.2;
        for (i = botmLayer+1; i <= surfLayer; i++) {
            Lake[i].Temp += ((KSED*(TYEAR-Lake[i].Temp)/ZSED)*
                      (Lake[i].LayerArea-Lake[i-1].LayerArea)*
                       LayerThickness[i]*noSecs)/(SPHEAT*Lake[i].Density*Lake[i].LayerVol);
        }
        Lake[botmLayer].Temp += ((KSED*(TYEAR-Lake[botmLayer].Temp)/ZSED)*
                                   Lake[botmLayer].LayerArea*LayerThickness[botmLayer] *
                               noSecs)/(SPHEAT*Lake[botmLayer].Density*Lake[botmLayer].LayerVol);
    }

    // precipitation, evaporation in the absence of ice
    if (! ice) {
    AED_REAL catch_runoff = 0.;

        if ( catchrain ) {
            // compute runoff in m3 for this time step
            if ( MetData.Rain > rain_threshold )
                catch_runoff = (MaxArea - Lake[surfLayer].LayerArea) *
                    ( MetData.Rain - rain_threshold )  *  (noSecs/SecsPerDay) * runoff_coef;
        }

        SurfData.dailyEvap += SurfData.Evap * noSecs * Lake[surfLayer].LayerArea ;
        SurfData.dailyRain += MetData.Rain * (noSecs / SecsPerDay) * Lake[surfLayer].LayerArea;
        if ( catchrain ) SurfData.dailyRain += catch_runoff;

        // Rainfall composition.  NOTE that evaporation leaves salts (nutrients)
        // deposits them at a rate dependent on the input composition of rainfall
        // with changes in depth, not area, for the surface layer. therefore just
        // depths to get new composition. firstly evaporation
        Lake[surfLayer].Height += (SurfData.Evap * noSecs + (MetData.Rain) * (noSecs / SecsPerDay));
        if ( catchrain ) Lake[surfLayer].Height += catch_runoff / Lake[surfLayer].LayerArea;

        // Add snow directly to surface layer height if there is no ice.
        // If there is ice, snow will be handled in the next block
        // Use 1:10 rule for snow water equivalent (Any better out there??)
        Lake[surfLayer].Height += (MetData.Snow * Lake[surfLayer].LayerArea * (1/10) * (noSecs / SecsPerDay));
        //SurfData.dailySnow += (MetData.Snow * Lake[surfLayer].LayerArea * (1.0/10.0) * (noSecs / SecsPerDay));

// CAB  // MH VOLUME HACK
// CAB  if (Lake[surfLayer].Height < 0.1)  Lake[surfLayer].Height = 0.1001 ;
        recalc_surface_salt();
    }

    // Recalculate densities
    for (i = botmLayer; i <= surfLayer; i++)
        Lake[i].Density = calculate_density(Lake[i].Temp,Lake[i].Salinity);

    // Set ice cover flag
    if (Lake[surfLayer].Temp <= 0.0 && SurfData.HeightBlackIce == 0.) {
        ice = TRUE;
        SurfData.HeightBlackIce = 0.05;
        SurfData.HeightWhiteIce = 0.0;
        Lake[surfLayer].Height -= 0.05*(rho_ice_blue/Lake[surfLayer].Density);

        recalc_surface_salt();

        SurfData.HeightSnow = 0.0;
    }

    if ((SurfData.HeightBlackIce+SurfData.HeightWhiteIce) < 0.05  &&  ice) {
        Lake[surfLayer].Height = Lake[surfLayer].Height+SurfData.HeightBlackIce*
              (rho_ice_blue/Lake[surfLayer].Density)+SurfData.HeightWhiteIce*
              (rho_ice_white/Lake[surfLayer].Density)+SurfData.HeightSnow*
              (rho_snow/Lake[surfLayer].Density);

        recalc_surface_salt();

        ice = FALSE;
        SurfData.HeightBlackIce = 0.0;
        SurfData.HeightWhiteIce = 0.0;
        SurfData.HeightSnow = 0.0;
    }
    SurfData.RhoSnow = rho_snow;

    // for (i = botmLayer; i <= surfLayer; i++)
    //     printf("Light = %10.5f\n",Lake[surfLayer].Light);
    //     printf("surfLayer = %d\n",surfLayer);

#ifdef _VISUAL_C_
    free(LayerThickness);  free(heat);
#endif

}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
AED_REAL calculate_qsw(int kDays,     // Days since start of year for yesterday
                       int mDays,     // Days since start of year for today
                       int iclock,    // #seconds into the day
                       AED_REAL Latitude,
                       AED_REAL SWOld,     // Total solar radiation at the surface for yesterday
                       AED_REAL ShortWave, // Total solar radiation at the surface for today
                       AED_REAL WindSp)    // current wind speed (adjusted for possible ice)
{
    AED_REAL RADTIM0;  //# Start of time step as fraction of day in radians
    AED_REAL RADTIM1;  //# End of time step as fraction of day in radians
    AED_REAL SOLAR0;   //# Solar angle yesterday
    AED_REAL SOLAR1;   //# Solar angle today
    AED_REAL SOLARW;   //# Sun's arc from sunrise to sunset for yesterday
    AED_REAL SOLARY;   //# Sun's arc from sunrise to sunset for today
    AED_REAL SW0;      //# Total solar radiation at the surface for yesterday after ALBEDO
    AED_REAL SW1;      //# Total solar radiation at the surface for today after ALBEDO
    AED_REAL TD0;      //# Day length for yesterday
    AED_REAL TD1;      //# Day length for yesterday
    int jClock;        //# End of the current time step

    AED_REAL Albedo0;  //# Photosynthetically active radiation scattered previous day
    AED_REAL Albedo1;  //# Photosynthetically active radiation scattered present day
    AED_REAL Temp_ice; //# Temperature of the top layer of the ice cover model [0C]
    AED_REAL lQSW;

/*----------------------------------------------------------------------------*/

    lQSW = 0.;

    RADTIM1 = 0.;
    RADTIM0 = 0.;

    Albedo0 = 0.;
    Albedo1 = 0.;

    //# Initialize Temp_ice to the Air Temperature
    Temp_ice = MetData.AirTemp;

    // Put in wind factor, set ice flag, determine resultant albedo
    if (ice) {
        if ((SurfData.HeightBlackIce+SurfData.HeightWhiteIce) > 0.55) {
            if (Temp_ice <= -5.)                   Albedo0 = 0.6;
            else if (Temp_ice > -5.0 && Temp_ice < 0.) Albedo0 = 0.44 - 0.032 * Temp_ice;
            else if (Temp_ice >= 0.)               Albedo0 = 0.44;
        } else{
            Albedo0 = 0.08 + 0.44 * pow((SurfData.HeightBlackIce+SurfData.HeightWhiteIce-0.05), 0.28);

        }

        if (SurfData.HeightSnow > 0.0) {
            if (Temp_ice <= -5.)                   Albedo1 = 0.7;
            else if (Temp_ice > -5.0 && Temp_ice < 0.) Albedo1 = 0.5 - 0.04 * Temp_ice;
            else if (Temp_ice >= 0.)               Albedo1 = 0.5;

            if (SurfData.HeightSnow < 0.1)
                Albedo0 = Albedo1-(((0.1-SurfData.HeightSnow)/0.1)*(Albedo1-Albedo0));
            else
                Albedo0 = Albedo1;

            //Adjust albedo based on multiplicative factor and back down to 1 if too large
            Albedo0 = snow_albedo_factor * Albedo0;
            if(Albedo0 > 1.0){
                Albedo0 = 1.0;
            }

        }

        Albedo1 = Albedo0;
    } else {
        Albedo0 = 0.08;
        Albedo1 = 0.08;

        switch (albedo_mode) {
            case 1: // simple daily albedo calc as in Hamilton and Schladow
                if (Latitude > two_Pi) { //# Lake is in the northern hemisphere
                    Albedo0 = 0.08 - 0.02 * sin(two_Pi * kDays/365 - (halfPi));
                    Albedo1 = 0.08 - 0.02 * sin(two_Pi * mDays/365 - (halfPi));
                } else { //# Lake is in the southern hemisphere
                    Albedo0 = 0.08 - 0.02 * sin(two_Pi * kDays/365 + (halfPi));
                    Albedo1 = 0.08 - 0.02 * sin(two_Pi * mDays/365 + (halfPi));
                }
                break;
            case 2: { // Briegleb et al. (1986), B scheme in Scinocca et al 2006.
                    AED_REAL sza = zenith_angle(Longitude, Latitude, mDays, iclock, timezone_r);  //degrees
                    if (sza < 80) {
                        AED_REAL csza = cos(Pi * sza/180);
                        Albedo1 = ((2.6 / (1.1 * pow(csza, 1.7) - 0.065)) +
                                  (15 * (csza-0.1) * (csza-0.5) * (csza-1))) /100.;
                    } else
                        Albedo1 = 0.3;
                }
                break;
             case 3: { // Yajima and Yamamoto (2014).
                    AED_REAL sza = zenith_angle(Longitude, Latitude, mDays, iclock, timezone_r);  //degrees
                    AED_REAL csza = cos(Pi * sza/180);

                    // from excel, need fixing
                    Albedo1 = MAX(0.01, 0.001 * MetData.RelHum * pow(1-csza, 0.33) -
                                        0.001 * WindSp * pow(1-csza, -0.57) -
                                        0.001 * 6 * pow(1-csza, 0.829));
                }
                break;
            default : break;
        }
    }

    //Add albedo tracking to help debug ice/snow duration
    SurfData.albedo = Albedo1;

    if ( subdaily )
        lQSW = ShortWave * (1.0-Albedo1);
    else {
        //# Determine the daylength (hours) from the latitude (entered as -ve for
        //# and the number of days since the start of the year. firstly determine
        SOLAR0 = -23.45 * sin((284 + kDays) * 2.0 * Pi / 365.0);
        SOLAR1 = -23.45 * sin((284 + mDays) * 2.0 * Pi / 365.0);

        // Determine the angle of the sun's arc from sunrise to sunset after
        // converting to radians, first determine solar declination
        SOLAR0 = SOLAR0 * Pi / 180;
        SOLAR1 = SOLAR1 * Pi / 180;
        SOLARW = 2.0 * acos(-tan(Latitude) * tan(SOLAR0));
        SOLARY = 2.0 * acos(-tan(Latitude) * tan(SOLAR1));

        // Calculate the daylength (Secs) using angular velocity 15 Degrees/Hour
        TD0 = 3600 * SOLARW / (15 * Pi/180);
        TD1 = 3600 * SOLARY / (15 * Pi/180);

        // Set the end of the glm time step
        jClock = iclock + noSecs;

        // Convert the start and end of the day to radians
        if (iclock < 0.5*TD0)
            RADTIM0 = iclock*Pi/TD0+Pi/2;
        else if (iclock >= 0.5*TD0 && iclock <= SecsPerDay-0.5*TD0)
            RADTIM0 = Pi;
        else if (iclock > SecsPerDay-0.5*TD0)
            RADTIM0 = (0.5*TD0+iclock-SecsPerDay)*Pi/TD0;

        if (jClock < 0.5*TD1)
            RADTIM1 = jClock*Pi/TD1+Pi/2;
        else if (jClock >= 0.5*TD1 && jClock <= SecsPerDay-0.5*TD1)
            RADTIM1 = Pi;
        else if (jClock > SecsPerDay-0.5*TD1)
            RADTIM1 = (0.5*TD1+jClock-SecsPerDay)*Pi/TD1;

        SW0 = SWOld*(1.0-Albedo0) * 86.4;
        SW1 = ShortWave*(1.0-Albedo1) * 86.4;

        // Determine the area under a sinusoidal curve and convert to a rate
        if (iclock <= (SecsPerDay - 0.5 * TD0) && jClock > (SecsPerDay - 0.5 * TD1) )
            lQSW = (0.5*SW0*1000.0*(cos(-1.0*RADTIM0)+1.0) + 0.5*SW1*1000.0*(1.0-cos(-1.0*RADTIM1)))/noSecs;
        else if (iclock <= (SecsPerDay - 0.5 * TD0) && jClock <= (SecsPerDay - 0.5 * TD1))
            lQSW = (0.5*SW0*1000.0*(cos(-1.0*RADTIM0)-cos(-1.0*RADTIM1)))/noSecs;
        else if (iclock > (SecsPerDay - 0.5 * TD0) && jClock > (SecsPerDay - 0.5 * TD1))
            lQSW = (0.5*SW1*1000.0*(cos(-1.0*RADTIM0)-cos(-1.0*RADTIM1)))/noSecs;
    }

    if (lQSW < 0.1) lQSW = 0.0;

    return lQSW;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
void recalc_surface_salt()
{
    AED_REAL OldVol, AddDensity, WaterMass;

    OldVol = Lake[surfLayer].LayerVol;

    resize_internals(1, surfLayer);

    AddDensity = calculate_density(Lake[surfLayer].Temp, zero);

    WaterMass = OldVol * (Lake[surfLayer].Density - AddDensity) +
                Lake[surfLayer].LayerVol * AddDensity;

    Lake[surfLayer].Salinity = (Lake[surfLayer].Density / WaterMass) *
                                                   OldVol * Lake[surfLayer].Salinity;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/



/******************************************************************************
 ******************************************************************************/

//#define gamma_air 0.0
//#define k_air     0.0
//#define mu_air    0.0
#define gamma_air 0.1
#define k_air     0.0548
#define mu_air    0.077

#define WIND_HEIGHT 10.0
#define HUMIDITY_HEIGHT 10.0

#define SIGN(a,b)  ( ((b) >= 0.) ? fabs(a) : -fabs(a) )

static AED_REAL psi_m(AED_REAL zL);
static AED_REAL psi_hw(AED_REAL zL);


/******************************************************************************/
AED_REAL  atmos_stability(AED_REAL *Q_latentheat,
                          AED_REAL *Q_sensible,
                          AED_REAL  WindSp,
                          AED_REAL  WaterTemp,
                          AED_REAL  AirTemp,
                          AED_REAL  p_atm,
                          AED_REAL  humidity_surface,
                          AED_REAL  humidity_altitude)
{
    AED_REAL U10, U_sensM, CD4;
    AED_REAL zL, L, zL0, z0, zS, G1, G2, G3, G5, G6;
    AED_REAL U_sensH, CDN10, CHWN10, CDN4, CDN3;
    AED_REAL CHW, CHWN, rCDN, P1, P2, P4;
    AED_REAL rho_air, T_virt, Ux, dT, dq;
    AED_REAL SH, LH, momn_flux;
    AED_REAL r_o, r_a, rho_a, rho_o, alpha_e;
    AED_REAL Q_latentheat_still, Q_sensible_still;

    int NCOUNT;

    AED_REAL vonK = 0.41;      // Von Karmans constant
    AED_REAL c_z0 = 0.0001;    // Default roughness
    AED_REAL cp_air = 1005.0;  // Specific heat of air

    AED_REAL zL_MAX;

/*----------------------------------------------------------------------------*/

    // Bound the iteration (e.g. 15 for 10m, 3 for 2m)
    zL_MAX = -15.0;

    if (fabs(WIND_HEIGHT-10.0) > 0.5)
        U10 = WindSp * (log(10.0/c_z0)/log(WIND_HEIGHT/c_z0));
    else
        U10 = WindSp;

    U_sensM = WindSp;

    CHWN10 = CH;

    // Calculate still air approximations and use this as a minimum
    // Fluxes to sill air see TVA Section 5.311 and 5.314

    // mixing ratios
    r_o = humidity_surface/(1-humidity_surface/c_gas);
    r_a = humidity_altitude/(1-humidity_altitude/c_gas);

    // density
    // 0.01 is for conversion from pascal to millibars
    rho_a = 0.348*((1+r_a)/(1+1.61*r_a))*(p_atm*0.01/(AirTemp+Kelvin));  // should use atm_density
    rho_o = 0.348*((1+r_o)/(1+1.61*r_o))*(p_atm*0.01/(WaterTemp+Kelvin));

    dT = WaterTemp - AirTemp;
    dq = humidity_surface - humidity_altitude;

    if (rho_a - rho_o > zero) {
        alpha_e = 0.137*0.5*(gamma_air/cp_air)* pow((9.81*(rho_a-rho_o)/(rho_a*mu_air*k_air)), (1/3.0));
        Q_sensible_still = -alpha_e * dT;
        Q_latentheat_still = -alpha_e * dq * Latent_Heat_Evap;
    } else {
        Q_sensible_still = zero;
        Q_latentheat_still = zero;
    }

    CDN10 = 0.001;
    if (U_sensM<0.01) {
        *Q_sensible = Q_sensible_still;
        *Q_latentheat = Q_latentheat_still;
    } else {
        // Neutral Drag Coefficient is a function of windspeed @ 10m
        //Option1
          //  if (U10 > 5.0)
          //      CDN10 = (1.0 + 0.07*(U10-5.0))/1000.0;
        //Option2
        CDN10 = 1.92E-7 * U10*U10*U10 + 0.00096;

        //Check
        if (CDN10>0.0025) CDN10 = 0.025;
    }

    // Atmospheric Density
    //?rho_air = RHO_AIR*Kelvin/(AirTemp+Kelvin);
    //rho_air = (p_atm)/(287.058 * (AirTemp+Kelvin));
    rho_air = rho_a;

    // Charnock computation of roughness, from Ux estimate.
    // 0.00001568 is kin viscosity
    // 0.012 is Charnock constant (alpha)
    Ux = sqrt(CDN10  * U_sensM * U_sensM);
    z0 = (0.012*Ux*Ux/9.81) + 0.11*0.00001568/Ux;
    CDN10 = pow(vonK/log(10./z0),2.0);

    // Estimate Surface Roughness Lengths
    z0 = 10.0/(exp(vonK/sqrt(CDN10)));
    zS = 10.0/(exp(vonK*vonK/(CHWN10*log(10.0/z0))));

    // Height Correction Factors
    G1 = log(10.0/z0);
    G2 = log(10.0/zS);
    G3 = log(HUMIDITY_HEIGHT/zS);
    G5 = log(HUMIDITY_HEIGHT/z0);
    G6 = log(WIND_HEIGHT/z0);

    CDN4 = CDN10*(G1*G1)/(G6*G6);    // Scale down to sensor heights
    CDN3 = CDN10*(G1*G1)/(G5*G5);
    CHWN = CHWN10*(G1*G2)/(G5*G3);
    CD4  = CDN4;                   // Initialize
    CHW  = CHWN;

    // Windspeed at the humidity sensor height
    U_sensH = U_sensM*(G5/G6);

    // Virtual Air Temperature
    T_virt = (AirTemp + Kelvin) * (1.0 + 0.61*humidity_altitude);

    // Heat Fluxes
    dT = WaterTemp - AirTemp;
    dq = humidity_surface - humidity_altitude;
    SH = CHW * rho_air * cp_air * U_sensH  * dT;
    LH = CHW * rho_air * U_sensH * dq;

    // Friction Velocity
    momn_flux = CD4 * rho_air * U_sensM*U_sensM;
    Ux = sqrt(momn_flux/rho_air);

    // Monin - Obukhov Length
    L = -rho_air *Ux*Ux*Ux * T_virt / (vonK * 9.81
                            * ((SH/cp_air) + 0.61*(AirTemp+Kelvin)*LH));
    if (fabs(L) < 0.5)
        L = SIGN(1.0e-20, dT);
    zL = HUMIDITY_HEIGHT/L;

    // Start Iterative Sequence for Heat Flux Calculations
    NCOUNT = 1;
    zL0 = zero;
    while ((fabs(zL - zL0) >= 0.0001*fabs(zL)) && (fabs(zL) <= fabs(zL_MAX))) {
        zL0 = zL;
        zL = WIND_HEIGHT/L;

        if (++NCOUNT>=15)
            break;

        // Calculate Drag Coefficient, CD
        P4 = psi_m(zL);
        rCDN = sqrt(CDN4);
        CD4 = CDN4/(1.0+CDN4*(P4*P4 - 2.0*vonK*P4/rCDN)/(vonK*vonK));

        // Calculate Humdity/Temp Coefficient, CHW
        zL = HUMIDITY_HEIGHT/L;

        P1 = psi_m(zL);
        P2 = psi_hw(zL);
        rCDN = sqrt(CDN3);
        CHW = CHWN/(1.0 + CHWN*(P1*P2 - (vonK*P2/rCDN)
                            - vonK*P1*rCDN/CHWN)/(vonK*vonK));

        // Recalculate Heat Fluxes
        SH = CHW * rho_air * cp_air * U_sensH * dT;
        LH = CHW * rho_air * U_sensH  * dq;
        momn_flux = CD4 * rho_air * U_sensM*U_sensM;

        // Recalculate Friction Velocity
        Ux = sqrt(momn_flux/rho_air);

        // Recalculate Monin - Obukhov Length
        L = -rho_air *Ux*Ux*Ux * T_virt / (vonK * 9.81
                                * ((SH/cp_air) + 0.61*(AirTemp+Kelvin)*LH));

        if (fabs(L) < 0.5)
            L = SIGN(1.0e-20,dT);
        zL = HUMIDITY_HEIGHT/L;
    } // enddo

    // Last Calculation - But 1st, check for high values
    if (fabs(zL)>fabs(zL_MAX))
        zL = SIGN(fabs(zL_MAX),zL);
    else
        zL = WIND_HEIGHT/L;

    P4 = psi_m(zL);
    rCDN = sqrt(CDN4);
    CD4 = CDN4/(1.0+CDN4*(P4*P4 - 2.0*vonK*P4/rCDN)/(vonK*vonK));
    zL = zL*HUMIDITY_HEIGHT/WIND_HEIGHT;
    P1 = psi_m(zL);
    P2 = psi_hw(zL);
    rCDN = sqrt(CDN3);
    CHW = CHWN/(1.0 + CHWN*(P1*P2 - (vonK*P2/rCDN)
            - vonK*P1*rCDN/CHWN)/(vonK*vonK));

    *Q_sensible = -CHW * rho_air * cp_air * U_sensH * dT;
    *Q_latentheat = -CHW * rho_air * U_sensH * dq * Latent_Heat_Evap;

    // Limit minimum to still air value
    if (Q_sensible_still < *Q_sensible)
        *Q_sensible = Q_sensible_still;
    if (Q_latentheat_still < *Q_latentheat)
        *Q_latentheat = Q_latentheat_still;

    // Link stability corrected drag to main code
    return CD4;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
static AED_REAL psi_m(AED_REAL zL)
{
    AED_REAL X;

    AED_REAL AA = 5.0;
    AED_REAL PIE = 3.14159;

/*----------------------------------------------------------------------------*/

    if (zL < 0.0) {
        X = pow(1.0 - 16.0*zL, 0.25);
        return 2.0*log((1.0+X)/2.0)+log((1.0+X*X)/2.0) - 2.0*atan(X) + PIE/2.0;
    } else if (zL > 0.0) {
        if (zL > 0.5) {
            if (zL > 10.0)
                return log(1.0*zL) - 0.76*zL - 12.093;
            return (0.5/(zL*zL)) -4.25/zL -7.0*log(zL)-0.852;
        }
        return -AA * zL;
    }
    return 0.0;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
static AED_REAL psi_hw(AED_REAL zL)
{
    AED_REAL X;

    AED_REAL AA = 5.0;

/*----------------------------------------------------------------------------*/

    if (zL < 0.0) {
        X = pow(1.0 - 16.0*zL, 0.25);
        return 2.0*log((1.0+X*X)/2.0);
    } else if (zL > 0.0) {
        if (zL > 0.5) {
            if (zL > 10.0)
                return log(zL) - 0.76*zL - 12.093;
            return (0.5/(zL*zL)) -4.25/zL -7.0*log(zL) -0.852;
        }
        return -AA * zL;
    }

    return 0.0;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
