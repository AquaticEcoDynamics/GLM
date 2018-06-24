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

AED_REAL Kw;             //# background light attenuation (m**-1)

int Num_WQ_Vars;         //# number of water quality variables
int Num_WQ_Ben;          //# number of benthic water quality variables
CLOGICAL atm_stab = FALSE;   // Account for non-neutral atmospheric stability

//------------------------------------------------------------------------------

AED_REAL CrestLevel; //# crest elevation of reservoir
AED_REAL LenAtCrest; //# length of reservoir at crest
AED_REAL WidAtCrest; //# width of reservoir at crest
AED_REAL VolAtCrest; //# volume at crest level
AED_REAL Base;       //# bottom elevation of reservoir
AED_REAL Benthic_Light_pcArea;
AED_REAL Benthic_Imin = 0.;
AED_REAL MaxArea;

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

//------------------------------------------------------------------------------

int NumDif;
AED_REAL mol_diffusivity[MaxDif];
// CD is the coef wind drag specified in the config, coef_wind_drag gets set to
// CD every time around the daily loop, coef_wind_drag is used in the loop
AED_REAL coef_wind_drag = 0.0013;
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

CLOGICAL non_avg = FALSE;
int deep_mixing = 2;

//
CLOGICAL catchrain = FALSE;
AED_REAL rain_threshold = 0.04;
AED_REAL runoff_coef = 0.3;

int      rad_mode = 0;
int      albedo_mode = 1;
int      cloud_mode = 1;


//------------------------------------------------------------------------------
// SNOWICE
AED_REAL snow_albedo_factor = 1.0;
AED_REAL snow_rho_max       = 300.;
AED_REAL snow_rho_min       = 50.;

//------------------------------------------------------------------------------
// SED_HEAT
CLOGICAL sed_heat_sw        = FALSE;
AED_REAL sed_temp_mean      = 9.7;
AED_REAL sed_temp_amplitude = 2.7;
AED_REAL sed_temp_peak_doy  = 151.;

//------------------------------------------------------------------------------
// FETCH
LOGICAL     fetch_sw = FALSE;
int         fetch_ndirs = 0;
AED_REAL   *fetch_dirs = NULL;
AED_REAL   *fetch_scale = NULL;
AED_REAL    fetch_height = 0.;
AED_REAL    fetch_porosity = 1.;

//------------------------------------------------------------------------------
// LITTORAL
CLOGICAL littoral_sw        = TRUE;

//------------------------------------------------------------------------------

AED_REAL timezone_r = 0.0, timezone_m = 0.0, timezone_i = 0.0, timezone_o = 0.0;

//------------------------------------------------------------------------------

int nDays;          //# number of days to simulate
AED_REAL timestep;
int noSecs;

//------------------------------------------------------------------------------

AED_REAL *WQ_Vars = NULL;  //# water quality array, nlayers, nvars

int       n_zones;
AED_REAL *zone_heights = NULL;
AED_REAL *zone_area = NULL;


//------------------------------------------------------------------------------
//  These for debugging
//------------------------------------------------------------------------------
CLOGICAL no_evap = FALSE;

void set_c_wqvars_ptr(AED_REAL *iwqv) { WQ_Vars = iwqv; }

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
void allocate_storage()
{
    MphLevelArea =    malloc(sizeof(AED_REAL) * Nmorph);
    dMphLevelArea =   malloc(sizeof(AED_REAL) * Nmorph);
    dMphLevelVol =    malloc(sizeof(AED_REAL) * Nmorph);
    dMphLevelVolda =  malloc(sizeof(AED_REAL) * Nmorph);
    MphLevelVol =     malloc(sizeof(AED_REAL) * Nmorph);
    MphLevelVoldash = malloc(sizeof(AED_REAL) * Nmorph);

    memset(MphLevelArea,    0, sizeof(AED_REAL) * Nmorph);
    memset(dMphLevelArea,   0, sizeof(AED_REAL) * Nmorph);
    memset(dMphLevelVol,    0, sizeof(AED_REAL) * Nmorph);
    memset(dMphLevelVolda,  0, sizeof(AED_REAL) * Nmorph);
    memset(MphLevelVol,     0, sizeof(AED_REAL) * Nmorph);
    memset(MphLevelVoldash, 0, sizeof(AED_REAL) * Nmorph);
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

    fprintf(of, "crest = %20.15f base = %20.15f VolAtCrest = %20.15f\n", CrestLevel, Base, VolAtCrest);

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
