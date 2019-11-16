/******************************************************************************
 *                                                                            *
 * glm_globals.c                                                              *
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

#include "glm.h"
#include "glm_types.h"
#include "glm_globals.h"
#include "aed_csv.h"


int MaxLayers;   //# Maximum number of layers in this sim
int NumLayers;   //# current number of layers
LakeDataType *Lake = NULL;


AED_REAL Latitude, Longitude;

AED_REAL DMin;    //# minimum layer thickness
AED_REAL DMax;    //# maximum layer thickness
AED_REAL VMin;    //# minimum layer volume
AED_REAL VMax;    //# maximum layer volume

int wq_calc = FALSE;

AED_REAL Kw;      //# background light attenuation (m**-1)

int Num_WQ_Vars;  //# number of water quality variables
int Num_WQ_Ben;   //# number of benthic water quality variables
int atm_stab = 0; //# Account for non-neutral atmospheric stability

//------------------------------------------------------------------------------

AED_REAL Base;        //# bottom elevation of reservoir
AED_REAL MaxHeight;   //# maxmimum height of reservoir
AED_REAL CrestHeight; //# crest height of reservoir
AED_REAL LenAtCrest;  //# length of reservoir at crest
AED_REAL WidAtCrest;  //# width of reservoir at crest
AED_REAL VolAtCrest;  //# volume at crest level
AED_REAL MaxVol;      //# volume at maximum level
AED_REAL MaxArea;
AED_REAL Benthic_Light_pcArea;
AED_REAL Benthic_Imin = 0.;

//------------------------------------------------------------------------------

int NumInf = 0;                    //# number of inflows
InflowDataType Inflows[MaxInf];    //# Array of Inflows

int NumOut = 0;                    //# Number of outflows
OutflowDataType Outflows[MaxOut];  //# Array of Outflows
int O2crit;
int O2critdep;
int O2critdays;
CLOGICAL MIXwithdraw = FALSE;
CLOGICAL COUPLoxy = FALSE;
AED_REAL WithdrawalTemp;
AED_REAL fac_range_upper = -1, fac_range_lower = -1;
AED_REAL MINlaketemp;

AED_REAL crest_width = 6.0;
AED_REAL crest_factor = 0.61;

CLOGICAL single_layer_draw = FALSE;
AED_REAL outflow_thick_limit = 100.0;

//------------------------------------------------------------------------------

int NumDif;
AED_REAL mol_diffusivity[MaxDif];
//# CD is the coef wind drag specified in the config, coef_wind_drag gets set to
//# CD every time around the daily loop, coef_wind_drag is used in the loop
AED_REAL coef_wind_drag = 0.0013;
AED_REAL coef_wind_chwn = 0.0013;
AED_REAL CD = 0.0013;
AED_REAL CE = 0.0013;
AED_REAL CH = 0.0013;

//------------------------------------------------------------------------------

MetDataType MetData;         //# Meteorological data
SurfaceDataType SurfData;    //# Surface Data

int subdaily = FALSE;

//------------------------------------------------------------------------------

int Nmorph = 0;             //# Number of data points in internal morphometry vector

AED_REAL  MphInc = 10.0;
AED_REAL *MphLevelArea    = NULL; //# area at each internal levels determined by linear interpolation
AED_REAL *dMphLevelArea   = NULL; //# gradients of area between 0.1m levels
AED_REAL *dMphLevelVol    = NULL; //# gradients of volume between 0.1m levels
AED_REAL *dMphLevelVolda  = NULL; //#
AED_REAL *MphLevelVol     = NULL; //# volume at each level determined by linear interpolation
AED_REAL *MphLevelVoldash = NULL; //#

//------------------------------------------------------------------------------
AED_REAL vel;
AED_REAL WaveNumSquared;
AED_REAL XMoment1;

//------------------------------------------------------------------------------

AED_REAL einff;   //# change in potential energy (see do_inflows)
AED_REAL coef_mix_KH = 0.3;     //# Kelvin-Helmholtz billowing effects
AED_REAL coef_mix_conv = 0.125; //# convective overturn
AED_REAL coef_mix_shear = 0.2;  //# shear efficiency
AED_REAL coef_mix_turb = 0.51;  //# unsteady effects
AED_REAL coef_wind_stir = 0.23; //# wind stirring
AED_REAL coef_mix_hyp = 0.5;    //# efficiency of hypolimnetic mixing
AED_REAL coef_mix_shreq = 1.0;  //# unsteady effects

CLOGICAL non_avg = FALSE;
int deep_mixing = 2;
int surface_mixing = 1;

//
CLOGICAL catchrain = FALSE;
AED_REAL rain_threshold = 0.04;
AED_REAL runoff_coef = 0.3;

int       rad_mode = 0;
int       albedo_mode = 1;
int       cloud_mode = 1;
int       light_mode = 1;
int       n_bands = 2;
AED_REAL  *light_extc = NULL;
AED_REAL  *energy_frac = NULL;

LOGICAL link_solar_shade = FALSE;
LOGICAL link_rain_loss   = FALSE;
LOGICAL link_bottom_drag = FALSE;
AED_REAL biodrag = 0.0;

AED_REAL salt_fall = 0.0;

int      density_model = 0;

AED_REAL albedo_mean = 0.08;         //# mean albedo
AED_REAL albedo_amplitude = 0.02 ;   //#  albedo seasonal amplitude
AED_REAL lw_factor   = 1.0;
AED_REAL lw_offset   = 0.0;

//------------------------------------------------------------------------------
// SNOWICE
AED_REAL snow_albedo_factor = 1.0;    //# scaling multiplier for computed albedo
AED_REAL snow_rho_max       = 300.;   //# maximum snow density allowed
AED_REAL snow_rho_min       = 50.;    //# minimum snow density allowed
AED_REAL snow_water_equivalent = 0.1; //# snow volume to water equivalent, 10:1
AED_REAL snow_rain_compact = 1.;      //# set at module level
AED_REAL K_ice_white = 2.3;           //# thermal conductivity of white ice
AED_REAL K_ice_blue = 2.0;            //# thermal conductivity of blue ice
AED_REAL K_water = 0.57;              //# molecular thermal conductivity of water
AED_REAL f_sw_wl1 = 0.7;              //# fraction of short wave radiation in first wavelength band
AED_REAL f_sw_wl2 = 0.3;              //# fraction of short wave radiation in second wavelength band
AED_REAL attn_ice_blue_wl1 = 1.5;     //# attenuation coefficient of the ice in the first spectral band
AED_REAL attn_ice_blue_wl2 = 20.;     //# attenuation coefficient of the ice in the second spectral band
AED_REAL attn_ice_white_wl1 = 6.0;    //# attenuation coefficient of the white ice in the first spectral band
AED_REAL attn_ice_white_wl2 = 20.;    //# attenuation coefficient of the white ice in the second spectral band
AED_REAL attn_snow_wl1 = 6.0;         //# attenuation coefficient of the snow in the first spectral band
AED_REAL attn_snow_wl2 = 20.;         //# attenuation coefficient of the snow in the second spectral band
AED_REAL rho_ice_blue = 917.0;        //# density of blue ice
AED_REAL rho_ice_white = 890.0;       //# density of white ice
AED_REAL min_ice_thickness = 0.05;    //# threshold thickness for new ice-on, or ice-off
AED_REAL dt_iceon_avg = 0.5;          //# moving average time-scale of water temp to identify ice-on transition

//------------------------------------------------------------------------------
// SEDIMENT
CLOGICAL sed_heat_sw        = FALSE;
int      sed_heat_model     = 0;
//AED_REAL sed_temp_mean        = 9.7;
//AED_REAL sed_temp_amplitude   = 2.7;
//AED_REAL sed_temp_peak_doy    = 151.;
AED_REAL  sed_heat_Ksoil    = 5.0;
AED_REAL  sed_temp_depth    = 0.1;
AED_REAL *sed_temp_mean     = NULL;
AED_REAL *sed_temp_amplitude = NULL;
AED_REAL *sed_temp_peak_doy = NULL;
AED_REAL *sed_reflectivity  = NULL;
AED_REAL *sed_roughness     = NULL;

//------------------------------------------------------------------------------
// FETCH
LOGICAL     fetch_sw = FALSE;
int         fetch_ndirs = 0;
AED_REAL   *fetch_dirs = NULL;
AED_REAL   *fetch_scale = NULL;
AED_REAL    fetch_height = 0.;
AED_REAL    fetch_porosity = 1.;

int         fetch_mode = 0;
AED_REAL    fetch_aws = 0.;
AED_REAL    fetch_xws = 0.;
char *      fetch_fws = NULL;

//------------------------------------------------------------------------------
// LITTORAL
CLOGICAL littoral_sw        = FALSE;

//------------------------------------------------------------------------------

AED_REAL timezone_r = 0.0, timezone_m = 0.0, timezone_i = 0.0, timezone_o = 0.0;

//------------------------------------------------------------------------------

int nDays;          //# number of days to simulate
AED_REAL timestep;
int noSecs;

//------------------------------------------------------------------------------

AED_REAL *WQ_Vars = NULL;  //# water quality array, nlayers, nvars

int       n_zones = 0;
//AED_REAL *zone_heights = NULL;
ZoneType *theZones = NULL;


//------------------------------------------------------------------------------
//  These for debugging
//------------------------------------------------------------------------------
CLOGICAL dbg_mix = FALSE;   //# debug output from mixer
CLOGICAL no_evap = FALSE;   //# turn off evaporation
int      quiet   = 0;       //# turn down output messages

void set_c_wqvars_ptr(AED_REAL *iwqv) { WQ_Vars = iwqv; }

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
void allocate_storage()
{
    MphLevelArea =    calloc(Nmorph, sizeof(AED_REAL));
    dMphLevelArea =   calloc(Nmorph, sizeof(AED_REAL));
    dMphLevelVol =    calloc(Nmorph, sizeof(AED_REAL));
    dMphLevelVolda =  calloc(Nmorph, sizeof(AED_REAL));
    MphLevelVol =     calloc(Nmorph, sizeof(AED_REAL));
    MphLevelVoldash = calloc(Nmorph, sizeof(AED_REAL));
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


#if DEBUG
/******************************************************************************/
void _debug_print_lake(FILE *of) {
    int i;

/*
    fprintf(of, "----------DEPTH----------------TEMP-----------------SALT-----------------DENS-----------------LVol------\n");
    for (i = 0; i < NumLayers; i++)
        fprintf(of, "%3d %16.11f %20.11f %20.11f %20.11f %20.11f %16.10f\n",
                    i, Lake[i].Height, Lake[i].Temp, Lake[i].Salinity, Lake[i].Density, Lake[i].LayerVol,Lake[i].Vol1);
*/
    fprintf(of, "MaxLayers %d NumLayers %d\n", MaxLayers, NumLayers);
    fprintf(of, "----------DEPTH----------------TEMP-----------------SALT-----------------DENS-----------------LVol--------------LArea----\n");
    for (i = 0; i < NumLayers; i++)
        fprintf(of, "%3d %16.11f %20.11f %20.11f %20.11f %20.11f %16.10f\n",
                    i, Lake[i].Height, Lake[i].Temp, Lake[i].Salinity, Lake[i].Density, Lake[i].LayerVol, Lake[i].LayerArea);
    fprintf(of, "-------------------------------------------------------------------------------------------------------------------------\n\n");
}
void debug_print_lake() { _debug_print_lake(stderr); }

/******************************************************************************/
void debug_initialisation(int which) {
    int i ;
    FILE *of = stderr;

    if (which) fprintf(of, "FORTRAN\n"); else fprintf(of, "C-----\n");

    _debug_print_lake(of);

    fprintf(of, "crest = %20.15f base = %20.15f VolAtCrest = %20.15f\n", CrestHeight, Base, VolAtCrest);

    fprintf(of, " Nmorph = %d\n", Nmorph);
    fprintf(of, "IDX----------StoLA----------------MphLevelVol--------------------dMphLevelVol------------------dMphLevelArea--------------\n");
    for (i = 0; i < Nmorph; i++)
        fprintf(of, "%3d, %20.15f %20.15f %20.15f %20.15f\n", i, MphLevelArea[i], MphLevelVol[i], dMphLevelVol[i], dMphLevelArea[i]);
    fputc('\n', of); fputc('\n', of);

    fprintf(of, "DMin %20.15f DMax %20.15f\nVMin %20.15f VMax %20.15f\n", DMin, DMax, VMin, VMax);
    fprintf(of, "EinFF %20.15f coef_mix_KH %20.15f coef_mix_conv %20.15f\nCS %f coef_mix_turb %f coef_wind_stir %f coef_mix_hyp %f\n",
                                       einff, coef_mix_KH, coef_mix_conv, coef_mix_shear, coef_mix_turb, coef_wind_stir, coef_mix_hyp);

    fprintf(of, "NumInf %d NumOut %d NumDif %d\n", NumInf, NumOut, NumDif);

    if (which) fprintf(of, "FORTRAN\n"); else fprintf(of, "C-----\n");

//exit(0);
}
void debug_initialisation_(int *which) { debug_initialisation(*which); }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#endif
