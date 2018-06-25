/******************************************************************************
 *                                                                            *
 * glm_surface.c                                                              *
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
#include <math.h>

#include "glm.h"

#include "glm_types.h"
#include "glm_const.h"
#include "glm_globals.h"

#include "glm_util.h"
#include "aed_time.h"
#include "glm_mixu.h"
#include "glm_layers.h"

#include "glm_wqual.h"  // for sediment heating - move?


//#define dbgprt(...) fprintf(stderr, __VA_ARGS__)
#define dbgprt(...) /* __VA_ARGS__ */


void solpond(int nband, int npoint,
             double depth, double *delz, double rb, double hdir, double anglei, double hdif,
             double *energy, double *absorb, double *gx);

static double expint(double ri, double cr, double h, double cmu, int n);
static double fct(double ri, double cr, double h, double cmu, int n);
static double fdif(double rindex, double critw, double h, double cmu, int n);
static double ref(double ri, double cr, double cmu);
static void gaus10(double ri, double cr, double c, double d, double (*f)(), int n, double h, double *v);

int  atmos_stability(     AED_REAL *Q_latentheat,
                          AED_REAL *Q_sensible,
                          AED_REAL  wind_speed,
                          AED_REAL  temp_water,
                          AED_REAL  temp_air,
                          AED_REAL  humidity_surface,
                          AED_REAL  humidity_altitude,
                          AED_REAL  rho_air,
                          AED_REAL  rho_o,
                          AED_REAL  p_atm,
                          AED_REAL  latent_heat_vap,
                          AED_REAL *coef_wind_drag,
                          AED_REAL *coef_wind_chwn,
                          AED_REAL *zonL );

/******************************************************************************
 * Module variables                                                           *
 ******************************************************************************/

int ice = FALSE;          // flag that tells if there is ice cover

// These are made available for the lake.csv output
AED_REAL zonL;                 // average z/L - MO atmospheric stability
AED_REAL Q_shortwave;          // Solar radiation at water surface
AED_REAL Q_sensibleheat;       // Sensible heat flux
AED_REAL Q_latentheat;         // Evaporative heat flux
AED_REAL Q_longwave;           // Net longwave heat flux
static AED_REAL  Q_iceout;     // Heat flux out of snow/ice surface
static AED_REAL  Q_icemet;     // Heat flux at the surface due to meteorological forcing
static AED_REAL  Q_watermol;   // Heat flux through water due to molecular conductivity
static AED_REAL  Q_icewater;   // Heat flux across water/ice interface
static AED_REAL  Q_surflayer;  // Heat flux through the water surface
static AED_REAL  Q_underflow;  // Heat flux through water due to flow under the ice
//static AED_REAL  U_flow;     // the velocity of the underflow

static AED_REAL  snow_rain_compact = 1. ; //update based on timestep and scaling

void recalc_surface_salt(void);

AED_REAL calculate_qsw(int kDays, int mDays, int iclock,
                       AED_REAL Latitude, AED_REAL SWOld,
                       AED_REAL ShortWave, AED_REAL WindSp);


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
 ******************************************************************************/
AED_REAL atm_density(AED_REAL atmosPressure, // (total) atmospheric pressure     [Pa]
                     AED_REAL vapPressure,   // water vapour (partial) pressure  [Pa]
                     AED_REAL AirTemp)       // dry bulb air temperature         [Cel]
{
/*
    // Dry air
    return (p_atm)/(287.058 * (AirTemp + Kelvin));
*/
  //r_o = humidity_surface/(1-humidity_surface/c_gas);

    // Moist air
    // mixing ratio: r = Mwater/(Mwater+Mdry_air)
    // (the grams of vapor per kg of dry air)
    AED_REAL r = mwrw2a * vapPressure/(atmosPressure - vapPressure);

    //rho_a = 100.*p_atm / (c_gas*(AirTemp+Kelvin));

    return 1.0/c_gas * (1 + r)/(1 + r/mwrw2a) * atmosPressure/(AirTemp+Kelvin);
}


/******************************************************************************
 * Performs thermal transfers across the lake surface (water and ice)         *
 * Ice cover extensions based on extension made by Brett Wallace &            *
 * David Hamilton 1996-1998                                                   *
 ******************************************************************************/
void do_surface_thermodynamics(int jday, int iclock, int LWModel,
                              AED_REAL Latitude, AED_REAL SWOld, AED_REAL ShortWave)

   // LWModel   // type of longwave radiation
   // SWOld     // Total solar radiation at the surface for yesterday
   // ShortWave // Total solar radiation at the surface for today
{
/*----------------------------------------------------------------------------*/
    const AED_REAL  K_ice_white = 2.3;     //# thermal conductivity of white ice
    const AED_REAL  K_ice_blue = 2.0;      //# thermal conductivity of blue ice
    const AED_REAL  K_water = 0.57;        //# molecular thermal conductivity of water
    const AED_REAL  Latent_Heat_Fusion = 334000.0;  // latent heat of fusion J/kg both ice & snow
    const AED_REAL  Temp_melt = 0.0;       //# temperature at which ice melts ~= temperature at which water freezes = 0.0oC
    const AED_REAL  f_sw_wl1 = 0.7;        //# fraction of short wave radiation in first wavelength band
    const AED_REAL  f_sw_wl2 = 0.3;        //# fraction of short wave radiation in second wavelength band
    const AED_REAL  attn_ice_blue_wl1 = 1.5;   //# attenuation coefficient of the ice in the first spectral band
    const AED_REAL  attn_ice_blue_wl2 = 20.;   //# attenuation coefficient of the ice in the second spectral band
    const AED_REAL  attn_ice_white_wl1 = 6.0;  //# attenuation coefficient of the white ice in the first spectral band
    const AED_REAL  attn_ice_white_wl2 = 20.;  //# attenuation coefficient of the white ice in the second spectral band
    const AED_REAL  attn_snow_wl1 = 6.0;   //# attenuation coefficient of the snow in the first spectral band
    const AED_REAL  attn_snow_wl2 = 20.;   //# attenuation coefficient of the snow in the second spectral band
    const AED_REAL  rho_ice_blue = 917.0;  //# density of blue ice
    const AED_REAL  rho_ice_white = 890.0; //# density of white ice
    const AED_REAL  eps_water = 0.985;     //# emissivity of the water surface
/*----------------------------------------------------------------------------*/

//  int nband = 4;
    int npoint = 10;

#ifndef _VISUAL_C_
    // The visual c compiler on doesn't like this so must malloc manually
    AED_REAL LayerThickness[MaxLayers],     //# Layer thickness (m)
             heat[MaxLayers];
    int layer_zone[MaxLayers];

//  AED_REAL energy[nband];
//  AED_REAL absorb[nband];
    AED_REAL gx[MaxLayers];
//  AED_REAL depths[MaxLayers];

#else
    AED_REAL *LayerThickness;
    AED_REAL *heat;
    int *layer_zone;

    //AED_REAL **solar;
//  AED_REAL *energy;
//  AED_REAL *absorb;
    AED_REAL *gx;
//  AED_REAL *depths;
#endif

    AED_REAL p_atm;          //# Atmospheric pressure in hectopascals, eg. 101300 Pa
    AED_REAL rho_air, rho_o; //# atm_density
    AED_REAL SatVap_surface; //# Saturated vapour pressure at surface layer or top ice layer
    AED_REAL altitude;
//  AED_REAL c_gas_m;
    AED_REAL latent_heat_vap;

    AED_REAL WindSp;         //# Wind speed corrected by multiplicative factor
    AED_REAL CloudCover;     //# Cloud cover as a fraction of 1 used to calculate incident long wave radiation

    AED_REAL dTemp;          //# Change in temperature due to atmospheric heat exchange
    AED_REAL dHt_WhiteIce;   //# Change in height of white ice when snow weight exceeds buoyancy potential

    AED_REAL Q_lw_in;        //# Long wave radiation from atmosphere W/m2
    AED_REAL Q_lw_out;       //# Long wave radiation emission W/m2
    AED_REAL Q_whiteice;     //# heat flux due to snow turning into white ice
    AED_REAL Q_snowice;      //# heat flux due to snow turning into white ice
    AED_REAL Q_rain;         //# heat flux due to rain either melting ice or freezing into ice

    AED_REAL AAA, BBB, CCC, DDD, EEE, FFF, GGG;
    AED_REAL Temp_ice;        //# temperature of the top layer of the ice cover model oC
    AED_REAL T01_OLD,T01_NEW; //# Not sure how these are used as Temp_ice always == 0 >= melting
    AED_REAL rho_snow;        //# Density of snow layer kg/m3
    AED_REAL rho_snow_old;    //# Density of snow layer compacted by rain kg/m3
    AED_REAL compact_snow;    //# Fraction of snow compacted by rain []
    AED_REAL K_snow;          //# thermal conductivity of snow [W/m2/oC]
    AED_REAL Q_latent_ice;    //# Latent heat of melting or freezing of ice [W/m2]

    AED_REAL SUMPO4, SUMTP, SUMNO3, SUMNH4, SUMTN, SUMSI;

    AED_REAL flankArea, NotPARLight_s, NotPARLight_sm1;

    //# NB These are set to 1.0 to remove a "possibly uninitialised" warning
    //#     They are initialised if litoral_sw is true, and only used if litoral_sw
    //#     so it should be OK, but should check from time to time.
    AED_REAL onshoreDensity = 1.0, offshoreDensity = 1.0;
    AED_REAL onshoreVol = 1.0, offshoreVol = 1.0;

    //# New parameters for sediment heat flux estimate
    AED_REAL KSED, TYEAR, ZSED;

    int i, z, j, wqidx;
    int kDays;
//  int non_nuetral_converged;

    //int nband, npoint;
    AED_REAL depth, rb, anglei, hdir, hdif;

/*----------------------------------------------------------------------------*/


    /**********************************************************************
    * Prior to surface heating, check for local runoff inputs from rain
    * and merge with the surface layer properties
    ***********************************************************************/

    AED_REAL catch_runoff = zero;

    if ( catchrain  && MetData.Rain > rain_threshold) {

//      AED_REAL VolSum = Lake[surfLayer].Vol1;  //# Total lake volume before this inflow

        // Compute runoff in m3 for this time step
        catch_runoff = ( MaxArea - Lake[surfLayer].LayerArea ) *
                         ( MetData.Rain - rain_threshold )  *  (noSecs/SecsPerDay) * runoff_coef;

        // Catchment runoff as a special inflow to the surface layer
        AED_REAL RunoffFlowRate = catch_runoff;
        Lake[surfLayer].Temp = combine(Lake[surfLayer].Temp, Lake[surfLayer].LayerVol, Lake[surfLayer].Density,
                                   MetData.AirTemp, RunoffFlowRate, calculate_density(MetData.AirTemp, zero+0.01));
        Lake[surfLayer].Salinity = combine(Lake[surfLayer].Salinity, Lake[surfLayer].LayerVol, Lake[surfLayer].Density,
                                           zero+0.01, RunoffFlowRate, calculate_density(MetData.AirTemp, zero+0.01));

        for (wqidx = 0; wqidx < Num_WQ_Vars; wqidx++)
            _WQ_Vars(wqidx, surfLayer) = combine_vol(_WQ_Vars(wqidx, surfLayer), Lake[surfLayer].LayerVol,
                                           zero, RunoffFlowRate);

        Lake[surfLayer].Density = calculate_density(Lake[surfLayer].Temp, Lake[surfLayer].Salinity);
        Lake[surfLayer].LayerVol = Lake[surfLayer].LayerVol+RunoffFlowRate;

        Lake[botmLayer].Vol1 = Lake[botmLayer].LayerVol;
        if (surfLayer != botmLayer) {
            for (j = (botmLayer+1); j <= surfLayer; j++)
                    Lake[j].Vol1 = Lake[j-1].Vol1 + Lake[j].LayerVol;
        }

        //# Make adjustments to correct layer heights following vol addition
        resize_internals(2,botmLayer);
        check_layer_thickness();

        SurfData.dailyRunoff += catch_runoff;

    }

    /**********************************************************************
    * Now get ready for surface heating
    ***********************************************************************/

#ifdef _VISUAL_C_
    LayerThickness = malloc(sizeof(AED_REAL) * MaxLayers);
    heat = malloc(sizeof(AED_REAL) * MaxLayers);
    layer_zone = malloc(sizeof(int) * MaxLayers);
//  energy = malloc(sizeof(AED_REAL) * nband);
//  absorb = malloc(sizeof(AED_REAL) * nband);
    gx = malloc(sizeof(AED_REAL) * MaxLayers);
#endif

    Q_longwave = 0.;
    SurfData.Evap = 0.;
    memset(heat, 0, sizeof(AED_REAL)*MaxLayers);
    memset(layer_zone, 0, sizeof(int)*MaxLayers);

    SUMPO4 = 0.0;
    SUMTP = 0.0;
    SUMNO3 = 0.0;
    SUMNH4 = 0.0;
    SUMTN = 0.0;
    SUMSI = 0.0;
    T01_NEW = 0.;
    T01_OLD = 0.;
    Q_latent_ice = 0.0;
    Q_whiteice = 0.0;
    Q_snowice = 0.0;
    Q_rain = 0.0;
    rho_snow = SurfData.RhoSnow;
//  underFlow = FALSE;

    //# Snow conductivity (eqn. 12 Rogers et al. 1995)
    //# To avoid division by zero set to 0.31 when no snow
    if (SurfData.HeightSnow  ==  0.)
        K_snow = 0.31;
    else
        K_snow = 0.021+(0.0042*rho_snow)+(2.2E-09*pow(rho_snow, 3));

    // Initialize the ice temperature to the air temperature
    Temp_ice = MetData.AirTemp;

    // Modify wind speed experienced by the surafce layer, if ice is present
    if (ice) WindSp = 0.00001;
    else     WindSp = MetData.WindSpeed;


    /**********************************************************************
     * ATMOSPHERIC CONDITIONS
     * Set atmospheric properties
     *********************************************************************/

    // Get atmospheric pressure and density

    // get altitude
    altitude = Base + Lake[surfLayer].Height + 5.;
    // update air pressure at altitude
    p_atm = ((100.*atm_pressure_sl) * pow((1 - 2.25577e-5*altitude),5.25588))/100.; //> hPa
    // gte the saturation vapour pressure at the water surface
    SatVap_surface = saturated_vapour(Lake[surfLayer].Temp); //> hPa
    // calculate gas constant for moist air
// This is not used anywhere
//  c_gas_m = c_gas*(1 + 0.608*(mwrw2a*MetData.SatVapDef/p_atm)); //> J/kg/K
    // update the latent heat of vaporisation, based on temperature
    latent_heat_vap = 2.501e6 - 2370*Lake[surfLayer].Temp; //> J/kg
    // update air density, right above the lake surface and at the ref height
    rho_air = atm_density(p_atm*100.0,MetData.SatVapDef*100.0,MetData.AirTemp); //>kg/m3
    rho_o = atm_density(p_atm*100.0,SatVap_surface*100.0,Lake[surfLayer].Temp); //>kg/m3

/*      printf(">>>>");
        printf(">altitude = %10.5f\n",altitude);
        printf(">p_atm = %10.5f\n",p_atm);
        printf(">SatVap_surface = %10.5f\n",SatVap_surface);
        printf(">c_gas_m = %10.5f\n",c_gas_m);
        printf(">latent_heat_vap = %10.5f\n",latent_heat_vap);
        printf(">rho_air = %10.5f\n",rho_air);
        printf(">rho_o = %10.5f\n",rho_o);
        printf(">MetData.SatVapDef = %10.5f\n",MetData.SatVapDef);
*/


    /**********************************************************************
     * SHORTWAVE INPUT
     * Assess surface shortwave flux hitting the water or ice
     * and add to the surface layer, plus work out heating
     *********************************************************************/

    //# Get shortwave intensity (W/m2) - either from provided sub-daily data,
    //  or approximate from daily total
    //  note, convert "julian days" to the "days since the start of the year"
    Q_shortwave = calculate_qsw(day_of_year(jday - 1), day_of_year(jday),
                                    iclock, Latitude, SWOld, ShortWave, WindSp);


    //# Into the surfLayer goes Qsw(surfLayer)-Qsw(surfLayer-1) over the
    //  area common to layers surfLayer and surfLayer-1, and Qsw(surfLayer)
    //  over the rest of the area. The units of heat[surfLayer] are
    // joules/sec; units of area(surfLayer) are m**2
    for (i = (botmLayer+1); i <= surfLayer; i++)
        LayerThickness[i] = Lake[i].Height-Lake[i-1].Height;

    LayerThickness[botmLayer] = Lake[botmLayer].Height;

    //# PAR entering the upper most (surface) layer
    if (!ice)
        // Assume PAR only (45% of the incident shortwave radiation
        // penetrates beyond the surface layer)
        Lake[surfLayer].Light = 0.45 * Q_shortwave;
    else
        // If there is ice cover the heating due to shortwave radiation
        // is attenuated across all three ice layers
        Lake[surfLayer].Light =
            (0.45 * Q_shortwave * ((f_sw_wl1 * exp(-attn_snow_wl1 * SurfData.HeightSnow)) +
                                  (f_sw_wl2 * exp(-attn_snow_wl2 * SurfData.HeightSnow)))) *
                                 ((f_sw_wl1 * exp(-attn_ice_white_wl1 * SurfData.HeightWhiteIce)) +
                                  (f_sw_wl2 * exp(-attn_ice_white_wl2 * SurfData.HeightWhiteIce))) *
                                 ((f_sw_wl1 * exp(-attn_ice_blue_wl1  * SurfData.HeightBlueIce)) +
                                  (f_sw_wl2 * exp(-attn_ice_blue_wl2  * SurfData.HeightBlueIce)));


    //# PAR entering the layer below the surface layer
    Lake[surfLayer-1].Light = Lake[surfLayer].Light * exp(-Lake[surfLayer].ExtcCoefSW*LayerThickness[surfLayer]);

    //# Not using Q_shortwave here as it hasnt been attenuated through ice.
    NotPARLight_s = (1-0.45)*(Lake[surfLayer].Light/0.45) ;
    NotPARLight_sm1 = NotPARLight_s * exp(-2.*Lake[surfLayer].ExtcCoefSW*LayerThickness[surfLayer]);

    //# Heating due to PAR
    flankArea = Lake[surfLayer].LayerArea - Lake[surfLayer-1].LayerArea;
    heat[surfLayer] = Lake[surfLayer-1].LayerArea * (Lake[surfLayer].Light - (1.-MAX(1.,sed_reflectivity[0]))*Lake[surfLayer-1].Light) ;
    //heat[surfLayer] = heat[surfLayer] + flankArea * (Lake[surfLayer].Light - Lake[surfLayer-1].Light*0.9);
    heat[surfLayer] = heat[surfLayer] + flankArea * (Lake[surfLayer].Light);
    //heat[surfLayer] = heat[surfLayer] + (Lake[surfLayer].LayerArea - Lake[surfLayer-1].LayerArea) * Lake[surfLayer].Light * 0.1);

    //# Heating due to NIR/UV
    //heat[surfLayer] = heat[surfLayer] + 0.55 * Q_shortwave * Lake[surfLayer].LayerArea;
    heat[surfLayer] = heat[surfLayer] + Lake[surfLayer-1].LayerArea * (NotPARLight_s - NotPARLight_sm1);
    heat[surfLayer] = heat[surfLayer] + flankArea * (NotPARLight_s - (1.-MAX(1.,sed_reflectivity[0]))*NotPARLight_sm1);

    if (littoral_sw) {
      onshoreDensity = calculate_density(Lake[onshoreLayer].Temp,Lake[surfLayer].Salinity);
      offshoreDensity = calculate_density(Lake[offshoreLayer].Temp,Lake[surfLayer].Salinity);
      offshoreVol = Lake[surfLayer-1].LayerArea * LayerThickness[surfLayer];
      onshoreVol = Lake[surfLayer].LayerVol - offshoreVol;

      heat[offshoreLayer] = Lake[surfLayer-1].LayerArea * (Lake[surfLayer].Light - Lake[surfLayer-1].Light)
                          + Lake[surfLayer-1].LayerArea * (NotPARLight_s - NotPARLight_sm1);
      heat[onshoreLayer]  = flankArea * (Lake[surfLayer].Light)
                          + flankArea * (NotPARLight_s - NotPARLight_sm1*0.95);
    }

    //# Summarise total daily short wave radiation (J/day), stored for lake.csv
    SurfData.dailyQsw += Lake[surfLayer].LayerArea * Q_shortwave * noSecs;



    // ---- MH TEST SOLPOND IN PROGRESS ---- //
    if(light_mode == 2){
    //# Advanced option - compute the light penetration suing the integral of light adsorption
    depth = Lake[surfLayer].Height;
    rb = 0.3;
    anglei = 10;
    hdir = Q_shortwave * 0.9;
    hdif = Q_shortwave * 0.1;


    npoint = NumLayers+1;
  //  nband = n_bands;

  //  printf(">solpond = \n");
  //  memset(energy, 0, sizeof(AED_REAL)*n_bands);
  //  memset(absorb, 0, sizeof(AED_REAL)*n_bands);
    memset(gx, 0, sizeof(AED_REAL)*MaxLayers);
  //  memset(depths, 0, sizeof(AED_REAL)*(surfLayer+1));
  //  for (i = botmLayer; i <= surfLayer; i++)
  //      Lake[i].Density = calculate_density(Lake[i].Temp,Lake[i].Salinity);
  //  printf(">solpond 2 = \n");
//    energy[0] = 0.51; energy[1] = 0.45; energy[2] = 0.035; energy[3] = 0.005;
//    absorb[0] = 1.; absorb[1] = 0.5; absorb[2] = 4.; absorb[3] = 4.;
//    for (i = 0; i < nband; i++) {
//      energy[i] = energy_frac[i];
//      absorb[i] = light_extc[i];
//    }

//    printf(">solpond 3 = \n");

//    printf(">solpond = %10.1f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",energy[0],energy[1],energy[2],energy[3],absorb[0],absorb[1],absorb[2],absorb[3]);

    //solpond(n_bands, npoint, depth, rb, hdir, anglei, hdif, energy, absorb, gx);

    solpond(n_bands, npoint, depth, LayerThickness, rb, hdir, anglei, hdif, energy_frac, light_extc, gx);

    //printf(">solpond = %10.1f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f %6.2f\n",
    //                       gx[0],gx[1],gx[2],gx[3],gx[4],gx[5],gx[6],gx[7],gx[8]);

  }
    // MH


    /**********************************************************************
     * SURFACE MASS FLUXES ON ICE
     * Assess surface mass fluxes of snow and rain on the ice
     * and check ice buoyancy for crackign and white ice formation
     *********************************************************************/
    if (iclock == 0 && ice) {
        AED_REAL BuoyantPotential ;

        if (SurfData.HeightSnow > 0.0) {

          // Snow cover as well as ice cover
          if (MetData.Snow > 0.0 && MetData.Rain >= 0.0) {
                // Snowfall on snow
                if (MetData.Rain == 0.0) MetData.Rain = MetData.Snow*0.10;  //Use 10:1 snow volume to water equivalent
                if (MetData.AirTemp > 0.0)
                    compact_snow = 0.166+0.834*(1.-exp(-1.*MetData.Rain/snow_rain_compact));
                else
                    compact_snow = 0.088+0.912*(1.-exp(-1.*MetData.Rain/snow_rain_compact));

                // Compact snow first
                rho_snow_old = rho_snow+(snow_rho_max-rho_snow)*compact_snow;
                SurfData.HeightSnow = SurfData.HeightSnow*rho_snow/rho_snow_old;

                // Determine the avg density of combined snow
                rho_snow = 1000.0*((rho_snow_old*SurfData.HeightSnow/1000.0)+MetData.Rain)/(SurfData.HeightSnow+(MetData.Snow));

                SurfData.HeightSnow = SurfData.HeightSnow+(MetData.Snow);

                // This should be a good estimate of water equivalent in from snow
                // Use rain because it should be volumetric equivalent of water (converted above, or supplied by user)
                SurfData.dailySnow += MetData.Rain * Lake[surfLayer].LayerArea;

                if (rho_snow > snow_rho_max) rho_snow = snow_rho_max;
                if (rho_snow < snow_rho_min) rho_snow = snow_rho_min;
                Q_rain = 0.0;

            } else if (MetData.Snow == 0.0 && MetData.Rain == 0.0) {
                // No snowfall or rainfall, compaction only of snow
                if (MetData.AirTemp > 0.0) compact_snow = 0.166;
                else                       compact_snow = 0.088;

                rho_snow_old = rho_snow+(snow_rho_max-rho_snow)*compact_snow;
                SurfData.HeightSnow = SurfData.HeightSnow*rho_snow/rho_snow_old;
                rho_snow = rho_snow_old;
                Q_rain = 0.0;

            } else if (MetData.Snow == 0.0 && MetData.Rain > 0.0) {
                // Rainfall on Snow.
                if (MetData.AirTemp > 0.0) {
                    //Check the air temperature and if AirTemp > 0 then
                    //  add the rain into the lake;
                    compact_snow = 0.166+0.834*(1.-exp(-1.*MetData.Rain/snow_rain_compact));
                    rho_snow_old = rho_snow+(snow_rho_max-rho_snow)*compact_snow;
                    SurfData.HeightSnow = SurfData.HeightSnow*rho_snow/rho_snow_old;
                    rho_snow = rho_snow_old;

                    Lake[surfLayer].Height = Lake[surfLayer].Height+MetData.Rain;
                    SurfData.dailyRain += MetData.Rain * Lake[surfLayer].LayerArea;
                    recalc_surface_salt();

                    if (Temp_ice == Temp_melt)
                        Q_rain = SPHEAT*(MetData.AirTemp-Temp_ice)*(MetData.Rain)/noSecs;
                    else
                        Q_rain = Latent_Heat_Fusion*MetData.Rain/noSecs;

                } else {
                    // If AirTemp < 0 then add the rainfall to snow
                    //MetData.Snow = MetData.Rain/snow_rho_max;
                    MetData.Snow = MetData.Rain/(snow_rho_max/rho0);
                    rho_snow_old = snow_rho_max;
                    SurfData.HeightSnow = SurfData.HeightSnow*rho_snow/rho_snow_old;
                    rho_snow = snow_rho_max;
                    SurfData.HeightSnow = SurfData.HeightSnow+MetData.Snow;
                    Q_rain = 0.0;
                    SurfData.dailyRain += MetData.Rain * Lake[surfLayer].LayerArea;
                }
            }
        } else {

            // No snow cover on the ice
            if (MetData.Snow > 0.0 && MetData.Rain >= 0.0) {

                // Snowfall on ice
                if (MetData.Rain == 0.0)MetData.Rain = MetData.Snow*0.10;
                SurfData.HeightSnow = MetData.Snow;

                rho_snow = rho0*MetData.Rain/MetData.Snow;
                SurfData.dailySnow += MetData.Rain * Lake[surfLayer].LayerArea;

                if (rho_snow > snow_rho_max)rho_snow = snow_rho_max;
                if (rho_snow < snow_rho_min)rho_snow = snow_rho_min;
                Q_rain = 0.0;

            } else if (MetData.Snow == 0.0 && MetData.Rain == 0.0) {

                // No snowfall or rainfall
                Q_rain = 0.0;

            } else if (MetData.Snow == 0.0 && MetData.Rain > 0.0) {

                // Rainfall on ice - need to know whether it will contribute to water
                if (MetData.AirTemp > 0.0) {
                    Lake[surfLayer].Height = Lake[surfLayer].Height+MetData.Rain;
                    SurfData.dailyRain += MetData.Rain * Lake[surfLayer].LayerArea;
                    recalc_surface_salt();

                    if (Temp_ice == Temp_melt)
                        Q_rain = SPHEAT*(MetData.AirTemp-Temp_ice)*(MetData.Rain)/noSecs;
                    else
                        Q_rain = Latent_Heat_Fusion*MetData.Rain/noSecs;

                } else {
                    // If AirTemp < 0
                    MetData.Snow = MetData.Rain/snow_rho_max;
                    rho_snow = snow_rho_max;
                    SurfData.HeightSnow = MetData.Snow;
                    SurfData.dailySnow += MetData.Snow * Lake[surfLayer].LayerArea * rho_snow/1000.0;
                    Q_rain = 0.0;
                }
            }
        }


        //# Archmides principle for the weight of Snow sitting on Ice -
        //  when the weight of snow exceeds the buoyancy of the ice
        //  cover, the ice will crack, and surface water will seep into snow
        //  layer leading to formation of White Ice
        if (rho_snow == 0.) rho_snow = snow_rho_min;

        BuoyantPotential = ((SurfData.HeightBlueIce*(rho0-rho_ice_blue) +
                             SurfData.HeightWhiteIce*(rho0-rho_ice_white))/rho_snow);

        if (SurfData.HeightSnow > BuoyantPotential) {

            dHt_WhiteIce = SurfData.HeightSnow-BuoyantPotential;

            Q_whiteice = ((Lake[surfLayer].Temp*SPHEAT+Latent_Heat_Fusion)*Lake[surfLayer].Density
                 *dHt_WhiteIce*(1.-(rho_snow/Lake[surfLayer].Density)))/SecsPerDay;

            SurfData.HeightWhiteIce = SurfData.HeightWhiteIce + SurfData.HeightSnow - BuoyantPotential;

            //Adjust surface layer height down based on water moving into white ice
            Lake[surfLayer].Height -= dHt_WhiteIce * (rho_ice_white - rho_snow) / Lake[surfLayer].Density;

            Q_snowice = (Q_whiteice*SurfData.HeightSnow)/(2.0*K_snow);
            SurfData.HeightSnow = BuoyantPotential;

            recalc_surface_salt();

        } else {
            dHt_WhiteIce = 0.0;
            Q_whiteice = 0.0;
            Q_snowice = 0.0;
        }

        // rainfall/snowfall and dry deposition effects on nutrients
        SUMPO4 = SUMPO4 + 0.145 + MetData.RainConcPO4 * MetData.Rain ;
        SUMTP  = SUMTP  + 0.368 + MetData.RainConcTp  * MetData.Rain ;
        SUMNO3 = SUMNO3 + 0.507 + MetData.RainConcNO3 * MetData.Rain ;
        SUMNH4 = SUMNH4 + 0.507 + MetData.RainConcNH4 * MetData.Rain ;
        SUMTN  = SUMTN  + 2.576 + MetData.RainConcTn  * MetData.Rain ;
        SUMSI  = SUMSI  + 0.0   + MetData.RainConcSi  * MetData.Rain ;

    }  // end iclock == 0 && ice


    /**********************************************************************
     * NON-PENETRATIVE HEAT FLUXES
     * Now do non-pentrative heat fluxes, depending on ice presence
     *********************************************************************/
    if (!ice) {
        //#  Evaporative heat flux affects top layer only
        SatVap_surface = saturated_vapour(Lake[surfLayer].Temp);

        //# Q_latentheat [W/m2] = CE * rho_air * latent heat * psychro const / air_presssure * windspeed * VPD
        Q_latentheat = -CE * rho_air * latent_heat_vap * (mwrw2a/p_atm) * WindSp * (SatVap_surface - MetData.SatVapDef);

        if (Q_latentheat > 0.0) Q_latentheat = 0.0; // no condensation

        //# Conductive/sensible heat gain only affects the top layer.
        //  Q_sensibleheat [W/m2] = CH * rho_air * specific heat * windspeed * temp diff
        Q_sensibleheat = -CH * rho_air * cp_air * WindSp * (Lake[surfLayer].Temp - MetData.AirTemp);

        //# If chosen by user, do atmospheric stability correction routine
        coef_wind_chwn = CH;
        coef_wind_drag = CD;

        if (atm_stab>0) {
//           non_nuetral_converged =
               atmos_stability(&Q_latentheat,
                               &Q_sensibleheat,
                                WindSp,
                                Lake[surfLayer].Temp,
                                MetData.AirTemp,
                                mwrw2a*SatVap_surface/p_atm,
                                mwrw2a*MetData.SatVapDef/p_atm,
                                rho_air,
                                rho_o,
                                p_atm,
                                latent_heat_vap,
                                &coef_wind_drag,
                                &coef_wind_chwn,
                                &zonL );

             SurfData.dailyzonL += zonL;
             //fprintf(stderr, " and now z/L %e = ...%e Q_l %e Q_s %e\n", zonL, coef_wind_drag, Q_latentheat, Q_sensibleheat);
        }

        //# Compute evaporative mass flux, in m/s, based on heat flux
        if ( no_evap )
            SurfData.Evap = 0.0;
        else
            SurfData.Evap = Q_latentheat / (latent_heat_vap * Lake[surfLayer].Density);

        //# Longwave emission (ie. longwave out), affecting only the top layer
        Q_lw_out = -Stefan_Boltzman * eps_water * pow((Kelvin+Lake[surfLayer].Temp), 4.0);

        //# Longwave absorption (ie. longwave in) also affects only top layer
        //  see Henderson-Sellers (1986) or Flerchinger (2009) for a good summary
        if (LWModel  ==  LW_CC) {
            // Cloud data is available
            AED_REAL eps_star = 0.8;  // default in case of a duff cloudmode value
            CloudCover = MetData.LongWave;
            switch (cloud_mode) {
                case 1:
                    // Idso and Jackson (1969)
                    // eps_star = (1.0 + 0.275*CloudCover)*(1.0 - 0.261 * exp(-0.000777 * pow(-MetData.AirTemp, 2.0))); //
                    eps_star = (1.0 + 0.275*CloudCover)*(1.0 - 0.261 * exp(-0.000777 * pow(MetData.AirTemp, 2.0)));
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
                    eps_star = (1.0 - pow(CloudCover, 2.796) ) * 1.24 * pow(MetData.SatVapDef/(MetData.AirTemp+Kelvin), 1/7) +
                                                           0.955 * pow(CloudCover, 2.796) ;
                    break;
            }

            Q_lw_in = (1 - 0.03) * eps_star * Stefan_Boltzman * pow((Kelvin+MetData.AirTemp), 4.0);
            Q_longwave = Q_lw_out + Q_lw_in;

        } else if (LWModel  ==  LW_IN) {
            // Incoming long-wave data is provided (assumes not corrected for lw albedo)
            Q_lw_in = (1 - 0.03) * MetData.LongWave;
            Q_longwave = Q_lw_out+Q_lw_in;

        } else if (LWModel  ==  LW_NET)
            // Net long-wave data is provided (ignore internally computed Q_lw_out)
            Q_longwave = MetData.LongWave;

        //# Update surface layer heat (J/s)
        heat[surfLayer] = heat[surfLayer]+(Q_latentheat+Q_sensibleheat+Q_longwave)*Lake[surfLayer].LayerArea;

        if (littoral_sw) {
          heat[offshoreLayer] = heat[offshoreLayer] +
                                (Q_latentheat+Q_sensibleheat+Q_longwave)*Lake[surfLayer-1].LayerArea;
          heat[onshoreLayer] = heat[onshoreLayer] +
                                (Q_latentheat+Q_sensibleheat+Q_longwave)*flankArea;
        }

        //# Daily heat budget (J/day)
        SurfData.dailyQe += Q_latentheat * Lake[surfLayer].LayerArea * noSecs;
        SurfData.dailyQh += Q_sensibleheat * Lake[surfLayer].LayerArea * noSecs;
        SurfData.dailyQlw += Q_longwave * Lake[surfLayer].LayerArea * noSecs;

    } else {

        //# Ice cover specific computation:
        //  The various atmospheric fluxes - evaporative, sensible heat, longwave
        //  for ice cover all depend on the ice surface temperature.  Note that
        //  wind here is not set to zero (USE MetData.WindSpeed, NOT WindSp].
        //  Emissivity for longwave OUT is 0.985 (eps_water)
        T01_NEW =  50.;
        T01_OLD = -50.;
        Temp_ice = zero;
        while (1) {
            //# Saturated vapor pressure above snow and ice is different than above
            //  water. Algoirthm taken from Jeong (2009), ultimately from
            //  Mellor (1964) "Properties of Snow"
            SatVap_surface = (1.+(0.00972*Temp_ice)+(0.000042*pow(Temp_ice, 2.)))
                             *saturated_vapour(Temp_ice);

            // I think this might be wrong, resulting value seems way too small, even for ice
            // Q_latentheat = -3.9 * MetData.WindSpeed * (SatVap_surface - MetData.SatVapDef);
            //# Q_latentheat [W/m2] = CE * rho_air * latent heat * psychro const / air_presssure * windspeed * VPD
            Q_latentheat = -CE * rho_air * latent_heat_vap * (mwrw2a/p_atm) * MetData.WindSpeed
                           * (SatVap_surface - MetData.SatVapDef);
            if (Q_latentheat > zero) Q_latentheat = zero;

            //Q_sensibleheat = -CH * MetData.WindSpeed * (Temp_ice - MetData.AirTemp);
            //# Q_sensibleheat [W/m2] = CH * rho_air * specific heat * windspeed * temp diff
            Q_sensibleheat = -CH * rho_air * cp_air * WindSp * (Temp_ice - MetData.AirTemp);

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

            //# This is the net meteorological flux that drives the ice algorithm
            //  flux due to precipitation on ice or snow
            Q_icemet = Q_latentheat+Q_sensibleheat+Q_longwave+Q_rain;

            //# Now determine the new ice/snow surface temperature based on the balance
            //  fluxes h_ice and the expression for the upward heat flux
            AAA = (1. - exp(-attn_snow_wl1 * SurfData.HeightSnow)) / (K_snow * attn_snow_wl1);
            BBB = exp(-attn_snow_wl1 * SurfData.HeightSnow)
                * (1.-exp(-attn_ice_white_wl1 * SurfData.HeightWhiteIce))/(attn_ice_white_wl1 * K_ice_white);
            FFF = exp((-attn_snow_wl1*SurfData.HeightSnow)-(attn_ice_white_wl1*SurfData.HeightWhiteIce))
                * (1.-exp(-attn_ice_blue_wl1*SurfData.HeightBlueIce))/(attn_ice_blue_wl1*K_ice_blue);

            CCC = (1.-exp(-attn_snow_wl2*SurfData.HeightSnow))/(K_snow*attn_snow_wl2);
            DDD = exp(-attn_snow_wl2*SurfData.HeightSnow)
                *(1.-exp(-attn_ice_white_wl2*SurfData.HeightWhiteIce))/(attn_ice_white_wl2*K_ice_white);
            GGG = exp((-attn_snow_wl2*SurfData.HeightSnow)-(attn_ice_white_wl2*SurfData.HeightWhiteIce))
                *(1.-exp(-attn_ice_blue_wl2*SurfData.HeightBlueIce))/(attn_ice_blue_wl2*K_ice_blue);
            EEE = (K_snow*K_ice_white*K_ice_blue)
                / (SurfData.HeightSnow*K_ice_white*K_ice_blue)+(SurfData.HeightBlueIce*K_snow*K_ice_blue)+(SurfData.HeightWhiteIce*K_snow*K_ice_white);

            Q_iceout = ((Temp_melt-Temp_ice-(Q_shortwave*(f_sw_wl1*(AAA+BBB+FFF)
                           + f_sw_wl2*(CCC+DDD+GGG)))+Q_snowice)*EEE)+Q_shortwave+Q_whiteice;

            //# Now compare the balance between the surface fluxes - and iterate
            //  by bisection (note the loop back in the iteration)

            if (fabs(Q_iceout+Q_icemet) > 1.0 && fabs(T01_NEW-T01_OLD) > 0.001) {
                if ((Q_iceout+Q_icemet) < 0.0) T01_NEW = Temp_ice;
                if ((Q_iceout+Q_icemet) > 0.0) T01_OLD = Temp_ice;
                Temp_ice = (T01_NEW+T01_OLD)/2.0;
            } else {
                break;  // start melting
            }
        } // end while

        //# Reset
        T01_NEW =  50.0;
        T01_OLD = -50.0;

        //# Evaporation in water equivalents (so use rho0)
        SurfData.Evap = Q_latentheat / Latent_Heat_Evap / rho0;

        //# Increment daily summaries
        SurfData.dailyEvap += (SurfData.Evap * noSecs * Lake[surfLayer].LayerArea);
        if (SurfData.HeightSnow > 0.){
            SurfData.HeightSnow += Q_latentheat/Latent_Heat_Evap/rho_snow*noSecs;
        } else if (SurfData.HeightWhiteIce > 0.) {
            SurfData.HeightWhiteIce += Q_latentheat/Latent_Heat_Evap/rho_ice_white*noSecs;
        } else {
            SurfData.HeightBlueIce += Q_latentheat/Latent_Heat_Evap/rho_ice_blue*noSecs;
        }

        //--------------------------------------------------------------------+
        // ICE MELTING & FREEZING @ SURFACE
        // Now compare the ice/snow surface temperature with the melting
        // temperature - if it is above the melting temp, then adjust
        // the thickness of the surface ice/snow layer accordingly
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

            Q_icemet = (Q_latentheat+Q_sensibleheat+Q_longwave+Q_rain);
            Q_iceout = ((Temp_melt-Temp_ice-(Q_shortwave*(f_sw_wl1*(AAA+BBB+FFF)+f_sw_wl2*(CCC+DDD+GGG))) +Q_snowice)*EEE)+Q_shortwave+Q_whiteice;

            //------------------------------------------------------------------
            // Now determine the new ice/snow/water surface temperature and height,
            // based on the balance between the surface fluxes h_ice and the expression
            // the upward heat flux as given by Patterson and Hamblin (1988]. Note
            // melting can occur at the surface and therefore water level will increase.
            // Note assumption that ice thickness won't change more than 1 cm in a timestep
            //------------------------------------------------------------------
            // LAW: reformatted this to read easier. Also fixed usage of densities for
            //      for melting ice

            if (SurfData.HeightSnow > 0.0) {
            //if (SurfData.HeightSnow != 0.0) {   //!MH Getting -ve Snow Height in lake.csv

                // If there is snow, melt that first
                if (rho_snow == 0.0) rho_snow = snow_rho_max;
                SurfData.dHt = (1/(Latent_Heat_Fusion*rho_snow))*(Q_icemet+Q_iceout)*noSecs;
                if ((SurfData.HeightSnow-SurfData.dHt) < 0.0)SurfData.dHt = SurfData.HeightSnow;
                SurfData.HeightSnow = SurfData.HeightSnow-SurfData.dHt;
                Lake[surfLayer].Height = Lake[surfLayer].Height+SurfData.dHt*(rho_snow/Lake[surfLayer].Density);
                recalc_surface_salt();

            } else if (SurfData.HeightWhiteIce > 0.){
                // Otherwise melt the white ice
                SurfData.dHt = (1/(Latent_Heat_Fusion*rho_ice_white))*(Q_icemet+Q_iceout)*noSecs;

                if ((SurfData.HeightWhiteIce-SurfData.dHt) < 0.0)
                    SurfData.dHt = SurfData.HeightWhiteIce;

                SurfData.HeightWhiteIce -= SurfData.dHt;

                Lake[surfLayer].Height += SurfData.dHt*(rho_ice_white/Lake[surfLayer].Density);
                recalc_surface_salt();

            } else {
                // Lastly, melt the blue ice
                SurfData.dHt = (1/(Latent_Heat_Fusion*rho_ice_blue))*(Q_icemet+Q_iceout)*noSecs;
                if ((SurfData.HeightBlueIce-SurfData.dHt) < 0.) {
                    SurfData.dHt = SurfData.HeightBlueIce;
                }
                SurfData.HeightBlueIce = SurfData.HeightBlueIce-SurfData.dHt;

                Lake[surfLayer].Height += SurfData.dHt*(rho_ice_blue/Lake[surfLayer].Density);
                recalc_surface_salt();
            } // end melting snow/white/blueice

        } // End melting if

       // Increment daily heat budget (MJ/day)
       SurfData.dailyQe += Q_latentheat * Lake[surfLayer].LayerArea * noSecs;
       SurfData.dailyQh += Q_sensibleheat * Lake[surfLayer].LayerArea * noSecs;
       SurfData.dailyQlw += Q_longwave * Lake[surfLayer].LayerArea * noSecs;
    }


    /**********************************************************************
     * APPLY HEATING TO WATER LAYERS
     * Now look at the ice or water interface
     *********************************************************************/
    if (surfLayer > botmLayer) {
        AED_REAL npl_p1, npli;

        for (i = surfLayer-1; i >= botmLayer; i-- )
            Lake[i].Light = Lake[i+1].Light * exp(-Lake[i+1].ExtcCoefSW*LayerThickness[i+1]);

        /*--------------------------------------------------------------------*
         * Into layer i goes QSW[i]-QSW(i-1) over the area common to layers   *
         * i and i-1 and QSW[i] over the rest of AREA[i]                      *
         * units of heat[i) are joules/sec; units of AREA[i] are 10**6 m**2   *
         *--------------------------------------------------------------------*/
        npl_p1 = NotPARLight_sm1;
        for (i = surfLayer-1; i >= (botmLayer+1); i--){
            flankArea = Lake[i].LayerArea - Lake[i-1].LayerArea;

            heat[i] = Lake[i-1].LayerArea * (Lake[i].Light - Lake[i-1].Light);
            heat[i] = heat[i] + flankArea * ( Lake[i].Light - Lake[i-1].Light*(1.-sed_reflectivity[0]));

            //MH non PAR temporary hack:
            npli = npl_p1 * exp(-Lake[i+1].ExtcCoefSW*2*LayerThickness[i+1]) ;
            heat[i] = heat[i] + (npl_p1-npli)*Lake[i].LayerArea;
            npl_p1 = npli;
        }
        heat[botmLayer] = Lake[botmLayer].Light * (1-(1.-sed_reflectivity[0])*exp(-Lake[botmLayer].ExtcCoefSW*LayerThickness[botmLayer])) * Lake[botmLayer].LayerArea;
        //heat[botmLayer] = heat[botmLayer] + npl_p1 * Lake[botmLayer].LayerArea;
    }

    //# Compute the temperature increase in non-surface layers over noSecs
    for (i = botmLayer; i <= surfLayer; i++) {
        if (fabs(heat[i]) >= 1E-20 && Lake[i].Density != 0.0 && Lake[i].LayerVol != 0.0)
            dTemp = heat[i]*noSecs/(SPHEAT*Lake[i].Density*Lake[i].LayerVol);
        else
            dTemp = 0.;

        Lake[i].Temp += dTemp;
    }

    if (littoral_sw) {
      dTemp = heat[offshoreLayer]*noSecs/(SPHEAT*offshoreDensity*offshoreVol);
      Lake[offshoreLayer].Temp += dTemp;

      dTemp = heat[onshoreLayer]*noSecs/(SPHEAT*onshoreDensity*onshoreVol);
      Lake[onshoreLayer].Temp += dTemp;
    }


    /**********************************************************************
     * ICE MELTING & FREEZING @ BOTTOM
     * The change in ice thickness at the bottom can now be determined
     * since the temperature of the lake water is readjusted for the surface
     * heat exchange, including the influence of the Ice/Snow
     *********************************************************************/
    if (ice) {
        //# Both ablation and accretion of the ice can occur at the ice-water
        //  interface change dht will be governed by the flux coming down from
        //  the surface to the water below the ice; see Patterson & Hamblin (1988)

        //# First calculate the flux through the ice
        AAA = (1. - exp(-attn_snow_wl1 * SurfData.HeightSnow)) / (K_snow * attn_snow_wl1);
        BBB = exp(-attn_snow_wl1*SurfData.HeightSnow)*(1.-exp(-attn_ice_white_wl1*SurfData.HeightWhiteIce))/ (attn_ice_white_wl1*K_ice_blue);
        FFF = exp((-attn_snow_wl1*SurfData.HeightSnow)-(attn_ice_white_wl1*SurfData.HeightWhiteIce))*
                                      (1.-exp(-attn_ice_blue_wl1*SurfData.HeightBlueIce))/(attn_ice_blue_wl1*K_ice_white);
        CCC = (1.-exp(-attn_snow_wl2*SurfData.HeightSnow))/(K_snow*attn_snow_wl2);
        DDD = exp(-attn_snow_wl2*SurfData.HeightSnow)*(1.-exp(-attn_ice_white_wl2*SurfData.HeightWhiteIce))/ (attn_ice_white_wl2*K_ice_blue);
        GGG = exp((-attn_snow_wl2*SurfData.HeightSnow)-(attn_ice_white_wl2*SurfData.HeightWhiteIce))*
                                      (1.-exp(-attn_ice_blue_wl2*SurfData.HeightBlueIce))/(attn_ice_blue_wl2*K_ice_white);
        EEE = (K_snow*K_ice_white*K_ice_blue)/((SurfData.HeightSnow*K_ice_white*K_ice_blue)+(SurfData.HeightBlueIce*K_snow*K_ice_blue)+
                                       (SurfData.HeightWhiteIce*K_snow*K_ice_white));

        Q_iceout = ((Temp_melt-Temp_ice-(Q_shortwave*(f_sw_wl1*(AAA+BBB+FFF)+f_sw_wl2*(CCC+DDD+GGG))) +Q_snowice)*EEE)+Q_shortwave+Q_whiteice;

        Q_icewater = Q_iceout-(Q_shortwave*f_sw_wl1*(1.-exp(-(attn_snow_wl1*SurfData.HeightSnow+attn_ice_blue_wl1*SurfData.HeightBlueIce+attn_ice_white_wl1*SurfData.HeightWhiteIce))))
                             -(Q_shortwave*f_sw_wl2*(1.-exp(-(attn_snow_wl2*SurfData.HeightSnow+attn_ice_blue_wl2*SurfData.HeightBlueIce+attn_ice_white_wl2*SurfData.HeightWhiteIce))));

        //# Now determine the flux through the lake water below the ice
        //  for temperature flux at the ice water interface use an exponential decrease
        //  fickian diffusion so use a gaussian distribution and adjust the diffuse
        //  Accordingly - see Farmer (1978) from p and h (1988)
        //  see eq 22 of Rogers et al and note that DZ has been adjusted
        Q_watermol = -K_water*(Temp_melt-Lake[surfLayer].Temp) / 0.039;

        //# Heat transfer due to ice underflow
        //  if (underFlow) Q_underflow = CSEN*(Lake[surfLayer].Density - rho0)*SPHEAT*U_FLOW*(Lake[surfLayer].Temp-Temp_melt);
        //  else
        Q_underflow = 0.0; //  Currently not being used as underFlow == .FALSE.

        //# Heat transfer to the water
        Q_surflayer = Q_watermol + Q_underflow + Q_latent_ice; //# LCB: Need to check as Q_latent_ice == 0.0

        //# Now we determine the amount of ablation or accretion of ice as
        //  given by qf-qw.  Once the ice has melted or formed, assume fluxes
        //  are in equilibrium. Correction for thermal contraction given
        //  by ratios of d is constant and area can be considered as constant
        //  as change in depth and time step is small
        SurfData.dHt = (Q_icewater-Q_surflayer)*noSecs/(Latent_Heat_Fusion*rho_ice_blue);
        if (SurfData.dHt < -1.*SurfData.HeightBlueIce) {
            if (SurfData.dHt < -1. * (SurfData.HeightBlueIce+SurfData.HeightWhiteIce) )
                SurfData.dHt = -1.*(SurfData.HeightBlueIce+SurfData.HeightWhiteIce);
            SurfData.HeightWhiteIce = SurfData.HeightWhiteIce+SurfData.HeightBlueIce+SurfData.dHt;
            SurfData.HeightBlueIce = 0.0;
        } else
            SurfData.HeightBlueIce = SurfData.HeightBlueIce+SurfData.dHt;

        //# Adjust water temperature for ice/water exchange, determine latent heat
        // LCB: moved calculation of latent heat of ice above temperature calculation ???
        Q_latent_ice = Latent_Heat_Fusion*rho_ice_blue*SurfData.dHt/noSecs;

        //Lake[surfLayer].Temp = Lake[surfLayer].Temp+(((K_water*((Temp_melt-Lake[surfLayer].Temp)/0.039))*Lake[surfLayer].LayerArea*
        //                         noSecs)+Q_latent_ice)/(SPHEAT*Lake[surfLayer].Density*Lake[surfLayer].LayerVol);
        Lake[surfLayer].Temp = Lake[surfLayer].Temp+((-Q_watermol + Q_latent_ice)*Lake[surfLayer].LayerArea*noSecs)
                                 /(SPHEAT*Lake[surfLayer].Density*Lake[surfLayer].LayerVol);

        Lake[surfLayer].Height = Lake[surfLayer].Height-SurfData.dHt*(rho_ice_blue/Lake[surfLayer].Density);

        recalc_surface_salt();
    }

    /**************************************************************************
     * SEDIMENT HEATING
     * Sediment "heating" factor now applied to any layer based on which zone
     * it sits above.
     **************************************************************************/
    if(sed_heat_sw){
        //# Input of heat from the sediments. Originally based on Rogers et al
        //  heat flux and sediment temperature measurements for Lake Mendota
        //  (Birge et al 1927), but made generic.  LayerThickness is required
        //  LAW: Modified to use cosine so user and specify peak day coefficient
        kDays = day_of_year(jday);
        ZSED = sed_temp_depth;
        KSED = sed_heat_Ksoil;

        if(benthic_mode == 1){
          //# Apply the same sediment heating parameters across all layers
          kDays = day_of_year(jday);
          TYEAR = sed_temp_mean[0] + sed_temp_amplitude[0] * cos(((kDays-sed_temp_peak_doy[0])*2.*Pi)/365.);
          for (i = botmLayer+1; i <= surfLayer; i++) {
             Lake[i].Temp += ((KSED*(TYEAR-Lake[i].Temp)/ZSED)*
                             (Lake[i].LayerArea-Lake[i-1].LayerArea)*
                             noSecs)/(SPHEAT*Lake[i].Density*Lake[i].LayerVol);
                             //LayerThickness[i]*noSecs)/(SPHEAT*Lake[i].Density*Lake[i].LayerVol);
          }
          Lake[botmLayer].Temp += ((KSED*(TYEAR-Lake[botmLayer].Temp)/ZSED)*
                                   Lake[botmLayer].LayerArea *
                                   //Lake[botmLayer].LayerArea*LayerThickness[botmLayer] *
                                   noSecs)/(SPHEAT*Lake[botmLayer].Density*Lake[botmLayer].LayerVol);
        }
        else  if(benthic_mode == 2){
          //# Apply the sediment zone specific heating parameters to overlying
          //  layers. First find which layers correspond to which zone
          for (i = botmLayer; i <= surfLayer; i++) {
             layer_zone[i] = 0;
             for (z = 0; z < n_zones; z++) {
                 if (Lake[i].Height<zone_heights[z] && Lake[i].Height>zone_heights[z-1])
                    layer_zone[i] = z;
             }
          }
          //# Now compute layer-specifc sed heating and increment temperature
          for (i = botmLayer+1; i <= surfLayer; i++) {
             TYEAR = sed_temp_mean[layer_zone[i]] + sed_temp_amplitude[layer_zone[i]] * cos(((kDays-sed_temp_peak_doy[layer_zone[i]])*2.*Pi)/365.);
             Lake[i].Temp += ((KSED*(TYEAR-Lake[i].Temp)/ZSED)*
                      (Lake[i].LayerArea-Lake[i-1].LayerArea)*
                      //LayerThickness[i]*noSecs)/(SPHEAT*Lake[i].Density*Lake[i].LayerVol);
                       noSecs)/(SPHEAT*Lake[i].Density*Lake[i].LayerVol);
          }
          TYEAR = sed_temp_mean[0] + sed_temp_amplitude[0] * cos(((kDays-sed_temp_peak_doy[0])*2.*Pi)/365.);
          Lake[botmLayer].Temp += ((KSED*(TYEAR-Lake[botmLayer].Temp)/ZSED)*
                                   //Lake[botmLayer].LayerArea*LayerThickness[botmLayer] *
                                   Lake[botmLayer].LayerArea *
                               noSecs)/(SPHEAT*Lake[botmLayer].Density*Lake[botmLayer].LayerVol);
        }
        if (littoral_sw) {
          TYEAR = sed_temp_mean[2] + sed_temp_amplitude[2] * cos(((kDays-sed_temp_peak_doy[2])*2.*Pi)/365.);
          Lake[onshoreLayer].Temp += ((KSED*(TYEAR-Lake[onshoreLayer].Temp)/ZSED)*onshoreVol * noSecs)
                                     /(SPHEAT*onshoreDensity*onshoreVol);
       }
    }


    /**************************************************************************
     * SURFACE MASS FLUXES (NO ICE)
     * Precipitation, evaporation in the absence of ice
     **************************************************************************/
    if (!ice) {

        AED_REAL rainvol = MAX( MetData.Rain,zero ) * (noSecs / SecsPerDay) * Lake[surfLayer].LayerArea;

        SurfData.dailyRain += rainvol;
        SurfData.dailyEvap += SurfData.Evap * noSecs * Lake[surfLayer].LayerArea ;

        //# Rainfall composition.  Note that evaporation leaves consituents &
        //  deposits them at a rate dependent on the composition of rainfall
        //  with changes in depth, not area, for the surface layer. Therefore just
        //  depths to get new composition. firstly evaporation
        Lake[surfLayer].Height += MAX( SurfData.Evap*noSecs,-0.9*Lake[surfLayer].Height )
                                  + rainvol / Lake[surfLayer].LayerArea;

        Lake[surfLayer].Temp = combine(Lake[surfLayer].Temp, Lake[surfLayer].LayerVol, Lake[surfLayer].Density,
                                       MetData.AirTemp, rainvol, calculate_density(MetData.AirTemp, zero+0.001));
        Lake[surfLayer].Salinity = combine(Lake[surfLayer].Salinity, Lake[surfLayer].LayerVol, Lake[surfLayer].Density,
                                       zero+0.001, rainvol, calculate_density(MetData.AirTemp, zero+0.001));
        for (wqidx = 0; wqidx < Num_WQ_Vars; wqidx++)
            _WQ_Vars(wqidx, surfLayer) = combine_vol(_WQ_Vars(wqidx, surfLayer), Lake[surfLayer].LayerVol,
                                                                   zero, rainvol);


        //# Add snow directly to surface layer height if there is no ice.
        //  If there is ice, snow will be handled in the next block
        //  Use 1:10 rule for snow water equivalent (Any better out there??)
        Lake[surfLayer].Height += MAX( MetData.Snow, zero) * Lake[surfLayer].LayerArea * (1./10.) * (noSecs / SecsPerDay);

        recalc_surface_salt();
    }

    //# Recalculate densities
    for (i = botmLayer; i <= surfLayer; i++)
        Lake[i].Density = calculate_density(Lake[i].Temp,Lake[i].Salinity);

    if (littoral_sw) {
        onshoreDensity = calculate_density(Lake[onshoreLayer].Temp,Lake[surfLayer].Salinity);
        offshoreDensity = calculate_density(Lake[onshoreLayer].Temp,Lake[surfLayer].Salinity);
    }

    //# Check and set ice cover flag
    if (Lake[surfLayer].Temp <= 0.0 && SurfData.HeightBlueIce == 0. && Lake[surfLayer].Height>0.1) {
        // Start a new blue ice layer
        ice = TRUE;
        SurfData.HeightBlueIce = 0.05;
        SurfData.HeightWhiteIce = 0.0;
        SurfData.HeightSnow = 0.0;
        Lake[surfLayer].Height -= 0.05*(rho_ice_blue/Lake[surfLayer].Density);

        recalc_surface_salt();
    }
    if ((SurfData.HeightBlueIce+SurfData.HeightWhiteIce) < 0.05  &&  ice) {
        Lake[surfLayer].Height = Lake[surfLayer].Height+SurfData.HeightBlueIce*
              (rho_ice_blue/Lake[surfLayer].Density)+SurfData.HeightWhiteIce*
              (rho_ice_white/Lake[surfLayer].Density)+SurfData.HeightSnow*
              (rho_snow/Lake[surfLayer].Density);

        recalc_surface_salt();

        ice = FALSE;
        SurfData.HeightBlueIce = 0.0;
        SurfData.HeightWhiteIce = 0.0;
        SurfData.HeightSnow = 0.0;
    }
    SurfData.RhoSnow = rho_snow;

#ifdef _VISUAL_C_
    free(LayerThickness);  free(heat);  free(layer_zone);

//  free(energy); free(absorb); free(gx);
    free(gx);
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
    AED_REAL TD1;      //# Day length for today
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

    //# Get current Solar Zenith Angle
    AED_REAL sza = MIN(90.0,zenith_angle(Longitude, Latitude, mDays, iclock, timezone_r));  //degrees
    AED_REAL csza = cos(Pi * sza/180);

    //# Put in wind factor, set ice flag, determine resultant albedo
    if (ice) {
        //# Albedo of ice cover depends on ice thickness, cover type & temperature
        //  This algotihm is based on Table 1 of Vavrus et al (1996).
        if ((SurfData.HeightBlueIce+SurfData.HeightWhiteIce) > 0.55) {
            if (Temp_ice <= -5.)                       Albedo0 = 0.6;
            else if (Temp_ice > -5.0 && Temp_ice < 0.) Albedo0 = 0.44 - 0.032*Temp_ice;
            else if (Temp_ice >= 0.)                   Albedo0 = 0.44;
        } else{
            Albedo0 = 0.08 + 0.44 * pow((SurfData.HeightBlueIce+SurfData.HeightWhiteIce-0.05), 0.28);
        }

        if (SurfData.HeightSnow > 0.0) {
            if (Temp_ice <= -5.)                   Albedo1 = 0.7;
            else if (Temp_ice > -5.0 && Temp_ice < 0.) Albedo1 = 0.5 - 0.04 * Temp_ice;
            else if (Temp_ice >= 0.)               Albedo1 = 0.5;

            if (SurfData.HeightSnow < 0.1)
                Albedo0 = Albedo1-(((0.1-SurfData.HeightSnow)/0.1)*(Albedo1-Albedo0));
            else
                Albedo0 = Albedo1;

            // Adjust albedo based on multiplicative factor and back down to 1 if too large
            Albedo0 = snow_albedo_factor * Albedo0;
            if(Albedo0 > 1.0){
                Albedo0 = 1.0;
            }

        }
        Albedo1 = Albedo0;

    } else {
        // Open water conditions
        Albedo0 = 0.08;
        Albedo1 = 0.08;

        switch (albedo_mode) {
            case 1: // simple daily albedo calc as in Hamilton and Schladow (1997)
                if (Latitude > two_Pi) { //# Lake is in the northern hemisphere
                    Albedo0 = 0.08 - 0.02 * sin(two_Pi * kDays/365 - (halfPi));
                    Albedo1 = 0.08 - 0.02 * sin(two_Pi * mDays/365 - (halfPi));
                } else { //# Lake is in the southern hemisphere
                    Albedo0 = 0.08 - 0.02 * sin(two_Pi * kDays/365 + (halfPi));
                    Albedo1 = 0.08 - 0.02 * sin(two_Pi * mDays/365 + (halfPi));
                }
                break;
            case 2: { // Briegleb et al. (1986), B scheme in Scinocca et al 2006
                    Albedo1 = ((2.6/(1*pow(csza, 1.7)+0.065)) + (15*(csza-0.1)*(csza-0.5)*(csza-1)))/100;
                }
                break;
             case 3: { // Yajima and Yamamoto (2014)
                    Albedo1 = MAX(0.02, 0.001 * MetData.RelHum * pow(1-csza, 0.33) -
                                        0.001 * WindSp * pow(1-csza, -0.57) -
                                        0.001 * 6 * pow(1-csza, 0.829));
                }
                break;
            default : break;
        }
    }
    //# Increment the daily albedo for this time-step, if within daylight hours
    if (sza<89.) SurfData.albedo += Albedo1;

    //# Now either assign BC value if sub-daily, or compute if BC is daily only
    if ( subdaily )
        lQSW = ShortWave * (1.0-Albedo1);
    else {
        //# Determine the daylength (hours) from the latitude (entered as -ve for
        //  and the number of days since the start of the year. firstly determine
        SOLAR0 = -23.45 * sin((284 + kDays) * 2.0 * Pi / 365.0);
        SOLAR1 = -23.45 * sin((284 + mDays) * 2.0 * Pi / 365.0);

        //# Determine the angle of the sun's arc from sunrise to sunset after
        //  converting to radians, first determine solar declination
        SOLAR0 = SOLAR0 * Pi / 180;
        SOLAR1 = SOLAR1 * Pi / 180;
        SOLARW = 2.0 * acos(-tan(Latitude) * tan(SOLAR0));
        SOLARY = 2.0 * acos(-tan(Latitude) * tan(SOLAR1));

        //# Calculate the daylength (Secs) using angular velocity 15 Degrees/Hour
        TD0 = 3600 * SOLARW / (15 * Pi/180);
        TD1 = 3600 * SOLARY / (15 * Pi/180);

        //# Set the end of the GLM time step
        jClock = iclock + noSecs;

        //# Convert the start and end of the day to radians
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

        //# Determine the area under a sinusoidal curve and convert to a rate
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

    Lake[surfLayer].Salinity = Lake[surfLayer].Salinity *
                               (Lake[surfLayer].Density / WaterMass) * OldVol ;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/




/******************************************************************************
 * ATMOSPHERIC STABILITY
 ******************************************************************************/

#define K_air    0.0280
#define D_air    0.0000214

#define WIND_HEIGHT 10.0
#define HUMIDITY_HEIGHT 10.0

#define SIGN(a,b)  ( ((b) >= 0.) ? fabs(a) : -fabs(a) )

static AED_REAL psi_m(AED_REAL zL);
static AED_REAL psi_hw(AED_REAL zL);

/******************************************************************************/
int atmos_stability(      AED_REAL *Q_latentheat,
                          AED_REAL *Q_sensible,
                          AED_REAL  wind_speed,
                          AED_REAL  temp_water,
                          AED_REAL  temp_air,
                          AED_REAL  humidity_surface,
                          AED_REAL  humidity_altitude,
                          AED_REAL  rho_air,
                          AED_REAL  rho_o,
                          AED_REAL  p_atm,
                          AED_REAL  latent_heat_vap,
                          AED_REAL *coef_wind_drag,
                          AED_REAL *coef_wind_chwn,
                          AED_REAL *zonL                )
{

    AED_REAL U10, U_sensM, U_sensH, Ux;
    AED_REAL zL, L, zL0, z0, zS, G1, G2, G3, G5, G6, Ldenom;
    AED_REAL CDN10, CHWN10, CDN4, CDN3, CHWN, rCDN, CD4, CHW;
    AED_REAL P1, P2, P4;
    AED_REAL T_virt, dT, dq;
    AED_REAL SH, LH, momn_flux;
    AED_REAL alpha_e, alpha_h, visc_k_air;
    AED_REAL Q_latentheat_still, Q_sensible_still;

    int atmos_count, atmos_status;

    AED_REAL vonK = 0.4100;    // von Karman's constant
    AED_REAL c_z0 = 0.0001;    // Default surafce roughness
    AED_REAL zL_MAX = -15.0;   // Bound the iteration (eg. 15 for 10m, 3 for 2m)

/*----------------------------------------------------------------------------*/
    atmos_status = 0;

    //# Some initial windspeed checks
    U_sensM = wind_speed;
    if (fabs(WIND_HEIGHT-10.0) > 0.5)
        U10 = wind_speed * (log(10.0/c_z0)/log(WIND_HEIGHT/c_z0));
    else
        U10 = wind_speed;

    CHWN10 = CH;

    //# Surface temperature and humidity gradients
    dT = temp_water - temp_air;
    dq = humidity_surface - humidity_altitude;

    visc_k_air = (1./rho_air)*(4.94e-8*temp_air + 1.7184e-5); //>m2/s

    //# First calculate still air approximations (free convenction) and use
    //  this as a minimum limit to compare with forced convection value
    //  For the fluxes to still air, see TVA Section 5.311 and 5.314

    //# Still air flux computation
    if (rho_air - rho_o > zero) {
        alpha_h = 0.137*0.5*K_air * pow( (g* fabs(rho_air-rho_o)/(rho_air*visc_k_air*D_air)), (1/3.0));
        alpha_e = alpha_h/cp_air;
        Q_sensible_still = -alpha_h * dT;
        Q_latentheat_still = -alpha_e * latent_heat_vap * dq ;
        // printf(">alpha_e = %10.5f\n",alpha_e);
    } else {
        Q_sensible_still = zero;
        Q_latentheat_still = zero;
    }

    if(atm_stab==2){
      // Assign free vs forced
      if (Q_sensible_still < *Q_sensible)
         *Q_sensible = Q_sensible_still;
      if (Q_latentheat_still < *Q_latentheat)
         *Q_latentheat = Q_latentheat_still;
      *zonL = zero;
      atmos_status = -2;

      return atmos_status;
    }

    /////////////

    printf("top = %10.5f\n",0.137*0.5*K_air);
    printf("bit = %10.5f\n",pow( (g* fabs(rho_air-rho_o)/(rho_air*visc_k_air*D_air)), (1/3.0)));
    printf("visc_k_air = %10.5f\n",visc_k_air);
    printf("D_air = %10.5f\n",D_air);
    printf("K_air = %10.5f\n",K_air);
    printf("dq = %10.5f\n",dq);
    printf("*Q_sensible_still = %10.5f\n",Q_sensible_still);
    printf("*Q_latentheat_still = %10.5f\n",Q_latentheat_still);

    /////////////



    //# Now check windspeed
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
        if (CDN10>0.0025) CDN10 = 0.0025;
    }


    //# Charnock computation of roughness, from Ux estimate.
    //  0.00001568 is kinematic viscosity ?
    //  0.012 is Charnock constant (alpha)
    Ux = sqrt(CDN10  * U_sensM * U_sensM);
    z0 = (0.012*Ux*Ux/g) + 0.11*visc_k_air/Ux;
    CDN10 = pow(vonK/log(10./z0),2.0);

    //# Estimate surface roughness lengths
    z0 = 10.0/(exp(vonK/sqrt(CDN10)));
    zS = 10.0/(exp(vonK*vonK/(CHWN10*log(10.0/z0))));

    //# Height correction factors
    G1 = log(10.0/z0);
    G2 = log(10.0/zS);
    G3 = log(HUMIDITY_HEIGHT/zS);
    G5 = log(HUMIDITY_HEIGHT/z0);
    G6 = log(WIND_HEIGHT/z0);

    CDN4 = CDN10*(G1*G1)/(G6*G6);    // Scale down to sensor heights
    CDN3 = CDN10*(G1*G1)/(G5*G5);
    CHWN = CHWN10*(G1*G2)/(G5*G3);
    CD4  = CDN4;                     // Initialize
    CHW  = CHWN;

    //# Windspeed at the humidity sensor height
    U_sensH = U_sensM*(G5/G6);

    //# Virtual air temperature
    T_virt = (temp_air+Kelvin) * (1.0 + 0.61*humidity_altitude);

    //# Heat fluxes based on bulk transfer (forced convection)
    SH = CHW * rho_air * cp_air * U_sensH  * dT;  //> W/m2
    LH = CHW * rho_air * U_sensH * dq;            //>

    //# Friction velocity
    momn_flux = CD4 * rho_air * U_sensM*U_sensM;
    Ux = sqrt(momn_flux/rho_air);

    //# Monin-Obukhov Length
    Ldenom = (vonK * g * ((SH/cp_air) + 0.61*(temp_air+Kelvin)*LH));
    if (fabs(Ldenom) < 1e-5) {
        zL = SIGN(zL_MAX,dT);
        L = HUMIDITY_HEIGHT/zL;
      }
    else {
        L = -(rho_air*Ux*Ux*Ux*T_virt) / Ldenom;
        zL = HUMIDITY_HEIGHT/L;
      }

      printf("U_sensM = %10.5f\n",U_sensM);
      printf("L = %10.5f\n",L);
      printf("zL = %10.5f\n",zL);

    //# Start iterative sequence for heat flux calculations
    atmos_count = 1;
    atmos_status = 1;
    zL0 = zero;
    while ((fabs(zL - zL0) >= 0.0001) ){ // && (fabs(zL) <= fabs(zL_MAX)+1)) {
        zL0 = zL;
        zL = WIND_HEIGHT/L;

        if (++atmos_count>=100){
            atmos_status = -1;
            break;
        }

        // Calculate drag coefficient, CD
        P4 = psi_m(zL);
        rCDN = sqrt(CDN4);
        CD4 = CDN4/(1.0+CDN4*(P4*P4 - 2.0*vonK*P4/rCDN)/(vonK*vonK));

        // Calculate Humdity/Temp coefficient, CHW
        zL = HUMIDITY_HEIGHT/L;

        P1 = psi_m(zL);
        P2 = psi_hw(zL);
        rCDN = sqrt(CDN3);
        CHW = CHWN/(1.0 + CHWN*(P1*P2 - (vonK*P2/rCDN)
                            - vonK*P1*rCDN/CHWN)/(vonK*vonK));

        // Recalculate heat and momn fluxes
        SH = CHW * rho_air * cp_air * U_sensH * dT;
        LH = CHW * rho_air * U_sensH  * dq;
        momn_flux = CD4 * rho_air * U_sensM*U_sensM;

        // Recalculate friction velocity
        Ux = sqrt(momn_flux/rho_air);

        // Recalculate Monin - Obukhov length
        L = -rho_air *Ux*Ux*Ux * T_virt / (vonK * g
                                * ((SH/cp_air) + 0.61*(temp_air+Kelvin)*LH));

        //printf("L = %10.5f\n",L);
        //printf("CHW = %10.5f\n",CHW);
        //printf("CD4 = %10.5f\n",CD4);

      //  if (fabs(L) < 0.5)
      //      L = SIGN(1.0e-20,dT);
        zL = HUMIDITY_HEIGHT/L;
    } // enddo

    if (atmos_status==-1)
       return atmos_status;


    //# Last calculation ... but 1st, check for high values
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

    *coef_wind_drag = CD4;
    *coef_wind_chwn = CHW;

    *Q_sensible = -CHW * rho_air * cp_air * U_sensH * dT;
    *Q_latentheat = -CHW * rho_air * U_sensH * dq * latent_heat_vap;

    printf("*atmos_count = %10d\n",atmos_count);
    printf("*CHW = %10.6f\n",CHW);
    printf("*CD4 = %10.6f\n",CD4);
    printf("zL = %10.5f\n",zL);
    printf("dT = %10.5f\n",dT);
    printf("dq = %10.5f\n",dq);
    printf("*Q_sensible = %10.5f\n",*Q_sensible);
    printf("*Q_latentheat = %10.5f\n",*Q_latentheat);

    *zonL = zL;
    if (atm_stab==3)
       return atmos_status;

    //# Limit minimum to still air value
    if (Q_sensible_still < *Q_sensible)
        *Q_sensible = Q_sensible_still;
    if (Q_latentheat_still < *Q_latentheat)
        *Q_latentheat = Q_latentheat_still;

    return atmos_status;

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



/******************************************************************************
 *                                                                            *
 * solpond watfiv ncs. ez. mae,cengel, time = (, 10), p = 20                  *
 *                                                                            *
 *       This program is to determine the local rate of solar energy          *
 *       absorption by water in a water column over specified number of       *
 *       positions.                                                           *
 *                                                                            *
 *       energy(k) = fraction of energy in wavelength band k                  *
 *       absorb(k) = volume absorption coefficient of water in band k (l/m)   *
 *       hdir   = direct solar rad. incident on a harizontal surface (w/m2)   *
 *       hdif   = diffuse solar rad. incident on a horizontal surface (w/m2)  *
 *       nband  = number of wavelength bands                                  *
 *       depth  = solar pond depth (m)                                        *
 *       gx(j)  = rate of solar energy absorption in water at level d (w/m3)  *
 *       npoint = number of positions at which gx i s to be calculated        *
 *       rindex = relative refractive index of water                          *
 *       rinver = i / rindex                                                  *
 *       critw  = cosine of critical angle of incidence for water side        *
 *       crita  = cosine of critical angle of incidence for air side (=0)     *
 *       anqlei = angle of incidence of direct radiation                      *
 *       canqle = cos(anqlei)                                                 *
 *       rmu    = cosine of the angle of refraction of direct solar rad.      *
 *       rb     = reflectivity of bottom surface (diffuse)                    *
 *                                                                            *
 ******************************************************************************/




/******************************************************************************
 *                                                                            *
 * This is the main routine with the inputs passed as parameters and the      *
 *  result returned in the array "gx"                                         *
 *                                                                            *
 ******************************************************************************/
void solpond(int nband, int npoint,
    double depth, double *delz, double rb, double hdir, double anglei, double hdif,
    double *energy, double *absorb, double *gx)
{
// implicit real*8(a-h,o-z),integer*4(i,n)

    const double crita = 0.0;
    const double rindex = 1.333;
    const double pi = 3.14159265358979;

    /* local variables */
    int i, j, l;
//  double rinver, critw, cangle, rmu, del, hkdif, hkdir;
    double rinver, critw, cangle, rmu, hkdif, hkdir;

    double a, a2, vdif, vint1, vint2, em2, vdifx, vdirx, vdir1, vdir2;
    double alpha, gbk, x, xpa, xm, gxk;

/*
    // These were originally read in:
    read, nband, depth, rb, hdir, anglei, hdif, npoint
    read, (energy(i), i=1, nband), (absorb(j), j=1, nband)
*/

    for (i = 0; i < npoint; i++)
        gx[i] = zero;

    rinver = 1.0 / rindex;
    critw = sqrt(1.0 - pow(rinver, 2));
    cangle = cos(anglei * pi/180.0);
    rmu = sqrt(pow(rindex, 2) + pow(cangle, 2) - 1.0) / rindex;



    for (l = 0; l < nband; l++) {
        hkdir = energy[l] * hdir;
        hkdif = energy[l] * hdif;
        a = depth * absorb[l];

        gaus10(rindex, critw, critw, 1.0, fdif, 1, a, &vdif);
        a2 = a * 2.0;
        gaus10(rinver, critw, 0.0, critw, fct, 1, a2, &vint1);
        gaus10(rinver, critw, critw, 1.0, fct, 1, a2, &vint2);
        alpha = vint1 + vint2;

        gbk = ((1.0-ref(rindex, crita, cangle)) * hkdir * exp(-a/rmu) +
           2.0 * pow(rindex, 2) * hkdif*vdif) / (1.0-2.0*rb*alpha);

        x = 0.0;
        //del = a/(npoint-1);
        //del = zero;


        gx[0] += hkdir + hkdif;

        for (j = 0; j < npoint; j++) {
            xpa = a+x;
            xm = a-x;

            gaus10(rindex, critw, 0.0, 1.0, expint, 2, xm, &em2);
            gaus10(rindex, critw, critw, 1.0, fdif, 0, x, &vdifx);
            gaus10(rinver, critw, 0.0, critw, fct, 1, xpa, &vdir1);
            gaus10(rinver, critw, critw, 1.0, fct, 1, xpa, &vdir2);
            vdirx = vdir1 + vdir2;

            gxk = absorb[l]*((1.0-ref(rindex, crita, cangle))*exp(-x/rmu)*
                   hkdir/rmu+2.0 * pow(rindex, 2) * hkdif*vdifx+2.0*rb*gbk*(vdirx+em2));

            gx[j+1] += gxk;

            //del += delz[((npoint-1)-j)-1];
            //x = (j+1) * del;
            x += delz[((npoint-1)-j)-1];
        }
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *    This subroutine performs numerical integration                          *
 *    using a 10 - point gausian quadrature                                   *
 *                                                                            *
 *    c = lower bound on integration                                          *
 *    d = upper bound on integration                                          *
 *    v = value of integral                                                   *
 *    f = function to be integrated                                           *
 *                                                                            *
 ******************************************************************************/
void gaus10(double ri, double cr, double c, double d, double (*f)(), int n, double h, double *v)
{
//  implicit real*8(a-h,o-z),integer*4(i,n)

    double x[10], w[10];
    double a, b;
    int i;

    *v = 0.0;

    if (c == d) return;

    x[0] = 0.148874338981631;
    x[1] = 0.433395394129247;
    x[2] = 0.679409568299024;
    x[3] = 0.865063366688985;
    x[4] = 0.973906528517172;
    w[0] = 0.295524224714753;
    w[1] = 0.269266719309996;
    w[2] = 0.219086362515982;
    w[3] = 0.149451349150581;
    w[4] = 0.066671344308688;

    for (i = 0; i < 5; i++) {
        x[i+5] = -x[i];
        w[i+5] =  w[i];
    }
    for (i = 0; i < 10; i++) {
        a = ((d-c) * x[i] + d + c) / 2.0;
        b = (d-c) * w[i] / 2.0;
        *v += (*f)(ri, cr, h, a, n) * b;
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/******************************************************************************
 *                                                                            *
 ******************************************************************************/
static double fct(double ri, double cr, double h, double cmu, int n)
{
//  implicit real*8(a-h,o-z), integer*4(i,n)
    return ref(ri, cr, cmu) * pow(cmu, n) * exp(-h/cmu);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/******************************************************************************
 *                                                                            *
 ******************************************************************************/
static double expint(double ri, double cr, double h, double cmu, int n)
//    this function is to determine exponential
{
//  implicit real*8(a-h,o-z),integer*4(i,n)
    return pow(cmu, (n-2)) * exp(-h/cmu);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/******************************************************************************
 *                                                                            *
 *    This function calculates the fresnel reflectivity                       *
 *    cmu = cosine of the angle of incidence                                  *
 *    ri = relative refractive index of the refracting medium                 *
 *    cr = critical cmu value below which reflectivity is unity               *
 *                                                                            *
 ******************************************************************************/
//    integral functions
static double ref(double ri, double cr, double cmu)
{
//  implicit real*8(a-h,o-z),integer*4(i,n)
    double ssine, rcos;

    if (cmu < cr) return 1.0;

    ssine = sqrt(pow(ri, 2) - 1.0+ pow(cmu, 2));
    rcos = pow(ri, 2) * cmu;
    return 0.50 * pow(((ssine-cmu)/(ssine+cmu)), 2) +
           0.50 * pow(((rcos-ssine)/(rcos+ssine)), 2);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/******************************************************************************
 *                                                                            *
 *    This function is to determine the refracted diffuse radiation           *
 *                                                                            *
 ******************************************************************************/
static double fdif(double rindex, double critw, double h, double cmu, int n)
{
//  implicit real*8(a-h,o-z),integer*4(i,n)
    double cmup;

    if (cmu < critw) return 0.0;

    cmup = sqrt(1.0 - pow(rindex, 2) * (1.0 - pow(cmu, 2)));
    return (1.0-ref(rindex, 0.0, cmup)) * pow(cmu, n) * exp(-h/cmu);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
#if 0
// Not currently used ?
static void soiltemp(int m, double wv)
{
#ifndef _VISUAL_C_
    // The visual c compiler doesn't like this so must malloc manually
    AED_REAL  w[m+1],t[m+1],tn[m+1],k[m+1],cp[m],a[m+1],b[m],c[m],d[m],z[m+1];
#else
    AED_REAL *w, *t, *tn, *k, *cp, *a, *b, *c, *d, *z;
#endif
    int ta = 20, am = 15, bd = 1.3, tb = 20;
    int i;

    AED_REAL ti, dt, da, f, g, mc, c1, c2, c3, c4;

#ifdef _VISUAL_C_
    w =  malloc(sizeof(AED_REAL) * (m+1));
    t =  malloc(sizeof(AED_REAL) * (m+1));
    tn = malloc(sizeof(AED_REAL) * (m+1));
    k =  malloc(sizeof(AED_REAL) * (m+1));
    cp = malloc(sizeof(AED_REAL) * m);
    a =  malloc(sizeof(AED_REAL) * (m+1));
    b =  malloc(sizeof(AED_REAL) * m);
    c =  malloc(sizeof(AED_REAL) * m);
    d =  malloc(sizeof(AED_REAL) * m);
    z =  malloc(sizeof(AED_REAL) * (m+1));
#endif

    k[0] = 20;    // boundary layer conductance in w/(m^2 k)

    for (i = 1; i <= m; i++) {
        z[i+1] = z[i]+.005* pow(1.5, (i-1));
        t[i] = tb;
    }

    t[m+1] = tb ; tn[m+1] = t[m+1] ; t[0] = tb;

    // ti is time of day;dt is time step (sec);da is day number
    ti = 0 ; dt = 3600 ; da = 0;
    f = .6 ; g = 1-f;

    mc = .12;   // clay fraction
    c1 = .65-.78*bd+.6*bd*bd ; c2 = 1.06*bd ; c3 = 1+2.6/sqrt(mc) ; c4 = .3+.1*bd*bd;

    for (i = 1; i <= m; i++) {
        cp[i] = (2400000*bd/2.65+4180000*wv)*(z[i+1]-z[i-1])/(2*dt);
        k[i] = (c1+c2*wv-(c1-c4)*exp(-pow((c3*wv), 4)))/(z[i+1]-z[i]);
    }

    do {
        ti = ti+dt/3600 ; if ( ti > 24 ) { ti = ti-24 ; da = da+1; }
        tn[0] = ta+am*sin(.261799*(ti-6));
        for (i = 1; i <= m; i++) {
            c[i] = -k[i]*f  ;  a[i+1] = c[i];
            b[i] = f*(k[i]+k[i-1])+cp[i];
            d[i] = g*k[i-1]*t[i-1]+(cp[i]-g*(k[i]+k[i-1]))*t[i]+g*k[i]*t[i+1];
        }
        d[1] = d[1]+k[0]*tn[0]*f;
        d[m] = d[m]+k[m]*f*tn[m+1];
        for (i = 1; i <= m-1; i++) {
            c[i] = c[i]/b[i];
            d[i] = d[i]/b[i];
            b[i+1] = b[i+1]-a[i+1]*c[i];
            d[i+1] = d[i+1]-a[i+1]*d[i];
        }
        tn[m] = d[m]/b[m];
        for (i = m-1; i>= 1; i--) {
            tn[i] = d[i]-c[i]*tn[i+1];
        }
//      print "day =";da,"hour =";ti
//      print"heat flux =";k[0]*(g*(t[0]-t[1])+f*(tn[0]-tn[1]));"w/m2"
//      print "depth","temperature","k[i]"
        for (i = 0; i <= m+1; i++) {
//          print z[i],tn[i],k[i];
            t[i] = tn[i];
        }
    } while (da < 5);
#ifdef _VISUAL_C_
    free(w); free(t); free(tn); free(k); free(cp);
    free(a); free(b); free(c);  free(d); free(z);
#endif
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#endif


#if 0
SUBROUTINE CalcZeroWindHeatFluxes(airTemp_C, surfTemp_C, vapPressure,      &
                                  atmosPressure, surfWaterDensity,         &
                                  zeroWindLatentHeatFlux,                  &
                                  zeroWindSensHeatFlux)
  IMPLICIT NONE
  REAL, INTENT(IN)  :: airTemp_C               !air temperature              [Cel]
  REAL, INTENT(IN)  :: surfTemp_C              !surface water temperature    [Cel]
  REAL, INTENT(IN)  :: atmosPressure           !(total) atmospheric pressure [hPa]
  REAL, INTENT(IN)  :: vapPressure             !water vapour pressure        [hPa]
  REAL, INTENT(IN)  :: surfWaterDensity        !density of surface water     [kg m^-3]
  REAL, INTENT(OUT) :: zeroWindLatentHeatFlux  !sensible heat flux density   [W m^-2]
  REAL, INTENT(OUT) :: zeroWindSensHeatFlux    !latent heat flux density     [W m^-2]

  ! ------ Local Constants ----------------
  ! -- For evaporation into still air: equation at top of p5.14, TVA Report:
  REAL, PARAMETER :: AIR_KINEMATIC_VISC  = 5.48E-02    ![m^2 hr^-1]
  REAL, PARAMETER :: AIR_MOL_HEAT_COND   = 0.1         ![kJ m^-1 hr^-1 K^-1]
  REAL, PARAMETER :: AIR_MOL_HEAT_DIFF   = 7.7E-02     ![m^2 hr^-1]
  REAL, PARAMETER :: AIR_SP_HEAT_CONST_P = 1.0         ![kJ kg^-1]
  REAL, PARAMETER :: CORR_COEFF          = 0.5         ![-]
  REAL, PARAMETER :: G_ACCN              = 1.28E+08    ![m hr^-2]

  ! ------ Local Variables ----------------
  REAL :: alpha_e          //coeff. for evap. into still air (TVA,p5.14)
  REAL :: alpha_h          //coeff. for sens. heat into still air (TVA,p5.17)
  REAL :: C_0              //water vapour conc. at water surface      [-]
  REAL :: C_a              //water vapour conc. of the ambient air    [-]
  REAL :: delta_conc       //difference in water vapour concs.        [-]
  REAL :: delta_rho        //difference in the air densities    [kg m^-3]
  REAL :: delta_T          //difference in the air temperatures     [Cel]
  REAL :: factor_1e        //variable for intermediate calculation step
  REAL :: factor_1h        //   "      "       "            "       "
  REAL :: factor_2         //   "      "       "            "       "
  REAL :: htEvap           //height of water evap. in this time step  [m]
  REAL :: rho_0            //air density at water surface       [kg m^-3]
  REAL :: rho_a            //air density of the ambient air     [kg m^-3]
  REAL :: satVapPressure   //saturated vapour pressure          [hPa]

  !--------------------------------------------------

  !---------
  ! 1. Calculate air densities:

  ! -- Saturated atmospheric density - i.e., atmospheric density at the
  ! -- water surface:
  satVapPressure = SatnVapourPressure(surfTemp_C)
  rho_0 = AirDensity_General(atmosPressure,satVapPressure,surfTemp_C)

  ! -- At some height above the water surface:
  rho_a = AirDensity_General(atmosPressure,vapPressure,airTemp_C)


  !---------
  ! 2. Calculate evaporation and sensible heat transfer coefficients:

  delta_rho = rho_a - rho_0

  factor_1e =   &
    0.137 * CORR_COEFF * AIR_MOL_HEAT_COND/(AIR_SP_HEAT_CONST_P*surfWaterDensity)
  factor_1h =   &
    0.137 * CORR_COEFF * AIR_MOL_HEAT_COND
  factor_2 =    &
    ( G_ACCN*ABS(delta_rho)/(rho_a*AIR_KINEMATIC_VISC*AIR_MOL_HEAT_DIFF) )**ONE_THIRD

  alpha_e = factor_1e * factor_2     !TVA, section 5.311, p5.14    [m hr^-1]
  alpha_h = factor_1h * factor_2     !TVA, section 5.314, p5.17    [m hr^-1]


  !---------
  ! 3. Calculate evaporation and sensible heat power flux densities:

  ! -- i) Evaporative heat:
  C_0 = MOL_WT_RATIO_WATER_AIR * satVapPressure/atmosPressure
  C_a = MOL_WT_RATIO_WATER_AIR *    vapPressure/atmosPressure
  delta_conc = C_0 - C_a       !TVA, section 5.311, p5.13           [-]

  htEvap = alpha_e * delta_conc * timeStep/NUM_SECS_PER_HR

  ! -- Denominator is to convert from [m hr^-1] to a power flux density
  ! -- in SI units, [W m^-2]:
  zeroWindLatentHeatFlux =   &
         htEvap * surfWaterDensity * LATENT_HEAT_VAP_WATER/timeStep


  ! -- ii) Sensible heat:
  delta_T = surfTemp_C - airTemp_C

  IF (delta_rho > 0) THEN
    ! -- Denominator ("3.6") is to convert from [kJ m^-2 hr^-1] to a
    ! -- power flux density in SI units [W m^-2]:
    zeroWindSensHeatFlux = - alpha_h*delta_T/3.6
  ELSE
    zeroWindSensHeatFlux = 0.0
  END IF


  //$Diagnostic output:
  //$WRITE(LOG_FILE,*)
  //$WRITE(LOG_FILE,'(1X,"====lws==============================")')
  //$WRITE(LOG_FILE,'(1X,"CalcZeroWindHeatFluxes:   at end of this S/R")')
  //$WRITE(LOG_FILE,'(1X,3X,"airTemp_C = ",F7.3,3X,"surfTemp_C = ",F7.3,3X, &
  //$          &"delta_T = ",F7.3)')  airTemp_C, surfTemp_C, delta_T
  //$WRITE(LOG_FILE,'(1X,3X,"alpha_e = ",ES13.6,3X,"alpha_h = ",ES13.6)')     &
  //$                                                    alpha_e, alpha_h
  //$WRITE(LOG_FILE,'(1X,"Uwind = 0",3X,"C_0 = ",ES13.6,3X,"C_a = ",ES13.6, &
  //$        &3X,"delta_conc = ",ES13.6)') C_0, C_a, delta_conc
  //$WRITE(LOG_FILE,'(1X,3X,"htEvap = ",ES13.6)')  htEvap
  //$WRITE(LOG_FILE,'(1X,3X,"zeroWindLatentHeatFlux = ",ES13.6,"  [W m^-2]",3X, &
  //$                   &"zeroWindSensHeatFlux = ",ES13.6,"  [W m^-2]")')    &
  //$                     zeroWindLatentHeatFlux, zeroWindSensHeatFlux
  //$WRITE(LOG_FILE,'(1X,"====/lws=============================",/)')


END SUBROUTINE CalcZeroWindHeatFluxes
#endif
