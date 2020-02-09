/******************************************************************************
 *                                                                            *
 * glm_globals.h                                                              *
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
#ifndef _GLM_GLOBALS_H_
#define _GLM_GLOBALS_H_

#ifdef _FORTRAN_SOURCE_

INTERFACE

     SUBROUTINE set_c_wqvars_ptr(iwqvars) BIND(C, name="set_c_wqvars_ptr")
        USE ISO_C_BINDING
#       if defined( _WIN32 ) && USE_DL_LOADER
        !DEC$ ATTRIBUTES DLLIMPORT :: set_c_wqvars_ptr
#       endif
        AED_REAL,INTENT(in) :: iwqvars(*)
     END SUBROUTINE set_c_wqvars_ptr

# if DEBUG
     SUBROUTINE debug_print_lake() BIND(C, name="debug_print_lake")
     END SUBROUTINE debug_print_lake

     SUBROUTINE debug_initialisation(which) BIND(C, name="debug_initialisation_")
        USE ISO_C_BINDING
        CINTEGER,INTENT(in) :: which
     END SUBROUTINE debug_initialisation
# endif


END INTERFACE

!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#else
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#include "glm.h"
#include "glm_types.h"

/* from glm_ncdf */
extern int ncid;

/* from glm_surf.F90 */
extern int ice;
extern AED_REAL AvgSurfTemp;
/*----------------------------------------------------------------------------*/
extern int MaxLayers;   //# Maximum number of layers in this sim
extern int NumLayers;   //# current number of layers
extern LakeDataType *Lake;


extern AED_REAL Latitude, Longitude;

extern AED_REAL DMin;    //# minimum layer thickness
extern AED_REAL DMax;    //# maximum layer thickness
extern AED_REAL VMin;    //# minimum layer volume
extern AED_REAL VMax;    //# maximum layer volume

extern AED_REAL Kw;            //# background light attenuation (m**-1)

extern int wq_calc;            //# are we doing water quality calcs
extern int Num_WQ_Vars;        //# number of water quality variables
extern int Num_WQ_Ben;         //# number of benthic water quality variables
extern AED_REAL *WQ_Vars;      //# water quality array : nlayers, nvars

extern int       n_zones;      //# number of sediment zones
//extern AED_REAL *zone_heights; //# heights for sed_zones
extern ZoneType *theZones;

/*----------------------------------------------------------------------------*/

extern AED_REAL MaxHeight;   //# maxmimum height of reservoir
extern AED_REAL CrestHeight; //# crest height of reservoir
extern AED_REAL LenAtCrest;  //# length of reservoir at crest
extern AED_REAL WidAtCrest;  //# width of reservoir at crest
extern AED_REAL VolAtCrest;  //# volume at crest level
extern AED_REAL MaxVol;      //# volume at max level
extern AED_REAL Base;        //# bottom elevation of reservoir

extern AED_REAL Benthic_Light_pcArea;
extern AED_REAL Benthic_Imin;
extern AED_REAL MaxArea;

/*----------------------------------------------------------------------------*/
// INFLOWS & OUTFLOWS
extern int NumInf;                 //# number of inflows
extern InflowDataType Inflows[];   //# Array of Inflows

extern int NumOut;                 //# Number of outflows
extern OutflowDataType Outflows[]; //# Array of Outflows
extern int O2crit;
extern int O2critdep;
extern int O2critdays;
extern CLOGICAL MIXwithdraw;
extern CLOGICAL COUPLoxy;
extern AED_REAL WithdrawalTemp;
extern AED_REAL fac_range_upper, fac_range_lower;
extern AED_REAL MINlaketemp;

extern AED_REAL crest_width;
extern AED_REAL crest_factor;

extern CLOGICAL single_layer_draw;
extern AED_REAL outflow_thick_limit;

/*----------------------------------------------------------------------------*/
//
extern int NumDif;
extern AED_REAL mol_diffusivity[];

/*----------------------------------------------------------------------------*/
// SURFACE
extern SurfaceDataType SurfData; //# Surface Data
extern MetDataType MetData;      //# Meteorological data
extern AED_REAL runoff_coef;
extern AED_REAL rain_threshold;
extern CLOGICAL catchrain;
extern int atm_stab;         //* Account for non-neutral atmospheric stability
extern AED_REAL coef_wind_drag;   //# = 0.0013;
extern AED_REAL coef_wind_chwn;   //# = 0.0013;
extern AED_REAL CD;               //# = 0.0013;
extern AED_REAL CE;               //# = 0.0013;
extern AED_REAL CH;               //# = 0.0013;
extern int      subdaily;         //# = FALSE;
extern int      rad_mode;
extern int      albedo_mode;
extern int      cloud_mode;
extern int      light_mode;
extern int      n_bands;
extern AED_REAL *light_extc;
extern AED_REAL *energy_frac;

extern LOGICAL  link_solar_shade;
extern LOGICAL  link_rain_loss;
extern LOGICAL  link_bottom_drag;
extern AED_REAL biodrag;

extern AED_REAL salt_fall;
/*----------------------------------------------------------------------------*/
// MORPHOMETRY
extern int Nmorph;                //# Number of data points

extern AED_REAL  MphInc;
extern AED_REAL *MphLevelArea;    //# area of each layer determined by linear interpolation
extern AED_REAL *dMphLevelArea;   //# gradients of area between 0.1m layers
extern AED_REAL *dMphLevelVol;    //# gradients of volume between 0.1m layers
extern AED_REAL *dMphLevelVolda;  //#
extern AED_REAL *MphLevelVol;     //# volume of each layer determined by linear interpolation
extern AED_REAL *MphLevelVoldash; //#

/*----------------------------------------------------------------------------*/
// MIXING
extern AED_REAL vel;
extern AED_REAL WaveNumSquared;
extern AED_REAL XMoment1;
extern AED_REAL einff;          //# change in potential energy (see do_inflows)
extern AED_REAL coef_mix_KH;    //# Kelvin-Helmholtz billows
extern AED_REAL coef_mix_conv;  //# convective overturn
extern AED_REAL coef_mix_shear; //# shear efficiency
extern AED_REAL coef_mix_turb;  //# unsteady effects
extern AED_REAL coef_wind_stir; //# wind stirring
extern AED_REAL coef_mix_hyp;   //# efficiency of hypolimnetic mixing
extern AED_REAL coef_mix_shreq;  //# unsteady effects

extern CLOGICAL mobility_off;
extern CLOGICAL non_avg;
extern int      deep_mixing;    //# = 0 => off > 0 => on
extern int      surface_mixing;

extern int      density_model;

extern AED_REAL albedo_mean;        //# mean albedo
extern AED_REAL albedo_amplitude;   //#  albedo seasonal amplitude
extern AED_REAL lw_factor ;
extern AED_REAL lw_offset ;

//# DepMX is the layer height of the meta top on the previous timestep.
extern AED_REAL DepMX;

extern AED_REAL PrevThick;   //# mixed layer thickness from previous time step

extern AED_REAL gPrimeTwoLayer;  //# Reduced gravity for int wave estimate

extern AED_REAL Energy_AvailableMix;  //# Total available energy to mix (carries over from previous timesteps)

extern AED_REAL Mass_Epi; //# Sigma mass of Epilimnion (surface layer after Kelvin-Helmholtz) kg

extern AED_REAL OldSlope;
extern AED_REAL Time_end_shear;  //# Time left before shear cut off [hours]
extern AED_REAL Time_start_shear; //# Time count since start of sim for shear period start [hours]
extern AED_REAL Time_count_end_shear;  //# Time count since start of sim for shear period end [hours]
extern AED_REAL Time_count_sim;  //# Time count since start of simulation [hours]

extern AED_REAL Half_Seiche_Period; //# One half the seiche period
extern AED_REAL Thermocline_Height; //# Height at the top of the metalimnion [m]
extern AED_REAL FO;
extern AED_REAL FSUM;
extern AED_REAL u_f;
extern AED_REAL u0;
extern AED_REAL u_avg;

/*----------------------------------------------------------------------------*/
// SNOWICE
extern AED_REAL snow_albedo_factor;
extern AED_REAL snow_rho_max;
extern AED_REAL snow_rho_min;
extern AED_REAL snow_water_equivalent;
extern AED_REAL snow_rain_compact;
extern AED_REAL  K_ice_white;
extern AED_REAL  K_ice_blue;
extern AED_REAL  K_water;
extern AED_REAL  f_sw_wl1;
extern AED_REAL  f_sw_wl2;
extern AED_REAL  attn_ice_blue_wl1;
extern AED_REAL  attn_ice_blue_wl2;
extern AED_REAL  attn_ice_white_wl1;
extern AED_REAL  attn_ice_white_wl2;
extern AED_REAL  attn_snow_wl1;
extern AED_REAL  attn_snow_wl2;
extern AED_REAL  rho_ice_blue;
extern AED_REAL  rho_ice_white;
extern AED_REAL  min_ice_thickness;
extern AED_REAL  dt_iceon_avg;

/*----------------------------------------------------------------------------*/
// SEDIMENT
extern CLOGICAL sed_heat_sw;
extern int      sed_heat_model;
extern AED_REAL sed_heat_Ksoil;
extern AED_REAL sed_temp_depth;
extern AED_REAL *sed_temp_mean;
extern AED_REAL *sed_temp_amplitude;
extern AED_REAL *sed_temp_peak_doy;
extern AED_REAL *sed_reflectivity;
extern AED_REAL *sed_roughness;

/*----------------------------------------------------------------------------*/
// FETCH
extern LOGICAL     fetch_sw;
extern int         fetch_ndirs;
extern AED_REAL   *fetch_dirs;
extern AED_REAL   *fetch_scale;
extern AED_REAL    fetch_height;
extern AED_REAL    fetch_porosity;

extern int         fetch_mode;
extern AED_REAL    fetch_aws;
extern AED_REAL    fetch_xws;
extern char       *fetch_fws;

/*----------------------------------------------------------------------------*/
// LITTORAL
extern CLOGICAL littoral_sw;

/*----------------------------------------------------------------------------*/
// TIME
extern AED_REAL timezone_r, timezone_m, timezone_i, timezone_o;
extern int nDays;          //# number of days to simulate
extern AED_REAL timestep;
extern int noSecs;

/*----------------------------------------------------------------------------*/
// DEBUGGING
extern CLOGICAL dbg_mix;   //# debug output from mixer
extern CLOGICAL no_evap;   //# turn off evaporation
extern int      quiet;     //# turn down output messages


/******************************************************************************/
void allocate_storage(void);
void set_c_wqvars_ptr(AED_REAL *iwqvars);
void debug_print_lake(void);
void debug_initialisation(int which);
void debug_initialisation_(int *which);

//#define _WQ_Vars(var,lyr) WQ_Vars[_IDX_2d(Num_WQ_Vars,MaxLayers,var,lyr)]
#define _WQ_Vars(var,lyr) WQ_Vars[_IDX_2d(MaxLayers,Num_WQ_Vars,lyr,var)]
#endif

/*============================================================================*/
#endif
