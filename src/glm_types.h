/******************************************************************************
 *                                                                            *
 * glm_types.h                                                                *
 *                                                                            *
 * Type declaration                                                           *
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
#ifndef _GLM_TYPES_H_
#define _GLM_TYPES_H_

#include "glm.h"

/************************* Important Note *************************************
 * The order of entries in these structures MUST match those in glm_types.F90 *
 ************************* Important Note *************************************/

#define MaxPar        37
#define MaxOut        20     /* Maximum number of outflows */
#define MaxInf        20     /* Maximum number of inflows */
#define MaxVars       60     /* Maximum number of variables */
#define MaxDif   (MaxVars+2) /* Maximum number of diffusing substances */

typedef int  LOGICAL;
typedef char varname[40];
typedef char filname[80];

/******************************************************************************
 *                                                                            *
 * Indexing macros                                                            *
 *                                                                            *
 * for an array declared in fortran as C(MaxI,MaxJ) and accessed as C(i,j)    *
 * we use this in C as C[_IDX_2d(MaxI,MaxJ,i,j)                               *
 *                                                                            *
 * a bit messy, but until we are all C we must live with the cards we are...  *
 *                                                                            *
 ******************************************************************************/
#define _IDX_2d(di,dj,i,j) (((di) * (j)) + (i))

#define _IDX_3d(di,dj,dk,   i,j,k)                         ((dk * dj * i) + (dk * j) + k)
#define _IDX_4d(di,dj,dk,dl,i,j,k,l)  ((dl * dk * dj * i) + (dl * dk * j) + (dl * k) + l)

/******************************************************************************
 * Macros for max and min - the first are the "standard" macros but suffer    *
 * from side effects - the latter are safer but have a greater overhead       *
 ******************************************************************************/
#ifdef _WIN32
  #define MAX(a,b) (((a) > (b)) ? (a):(b))
  #define MIN(a,b) (((a) < (b)) ? (a):(b))
#else
  #define MAX(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a > _b ? _a : _b; })
  #define MIN(a,b) \
   ({ __typeof__ (a) _a = (a); \
       __typeof__ (b) _b = (b); \
     _a < _b ? _a : _b; })
#endif

/******************************************************************************/

/*============================================================================*/
//TYPE DECLARATIONS

   /*===========================================================*/
   typedef struct StringT {
       int  Len;
       char S[40];
   } StringT;

   /*===========================================================*/
   // Structured type for inflow vars
   // An inflow will be an allocated array of MaxInf of these
   typedef struct InflowDataType {
       AED_REAL Alpha;           // half angle of stream
       AED_REAL DragCoeff;       // streambed drag coefficient
       AED_REAL Phi;             // streambed slope
       AED_REAL FlowRate;        // inflow flow rate
       AED_REAL Factor;          // scaling factor for inflow
       AED_REAL TemInf;          // inflow temperature
       AED_REAL SalInf;          // inflow salinity
       AED_REAL Dlwst;
       AED_REAL HFlow;
       AED_REAL TotIn;
       AED_REAL DIIns[MaxPar];   // inflow density
       AED_REAL DDown[MaxPar];   // downflow insertion depth
       AED_REAL QIns[MaxPar];    // inflow volume
       AED_REAL QDown[MaxPar];   // downflow volume
       AED_REAL TIns[MaxPar];    // inflow temperature
       AED_REAL TDown[MaxPar];   // downflow temperature
       AED_REAL SIns[MaxPar];    // inflow salinity
       AED_REAL SDown[MaxPar];   // downflow salinity
       AED_REAL DOld[MaxPar];

       AED_REAL WQIns[MaxPar][MaxVars];  // inflow water quality
       AED_REAL WQDown[MaxPar][MaxVars]; // downflow water quality
       AED_REAL WQInf[MaxVars];

       int  iCnt;
       int  NoIns;
       int  InPar[MaxPar];

       LOGICAL  SubmFlag;        // Is this a submerged inflow
   } InflowDataType;

   /*===========================================================*/
   // Structured type for outflow vars
   // An outflow will be an allocated array of MaxOut of these
   typedef struct OutflowDataType {
       int Type;                 // outflow type
       AED_REAL Hcrit;           // outlet height when crit O2
       int O2idx;                // O2 parameter idx in AED2/FABM
       char O2name;              // O2 parameter name in AED2/FABM
       AED_REAL TARGETtemp;      // Isotherm for withdrawal switch 4
       AED_REAL OLev;            // distance below surface level
       AED_REAL OLen;            // basin length at the outlet
       AED_REAL OWid;            // basin width at the outlet
       AED_REAL Draw;            // outflow volumes
       AED_REAL Factor;          // scaling factor for outflow
       LOGICAL  FloatOff;        // Is this a floating offtake
   } OutflowDataType;

   /*===========================================================*/
   // Structured type for key global lake environmental vars
   // A Lake will be an allocated array of MaxLayers of these
   typedef struct LakeDataType {
       AED_REAL Density;         // Density kg/m3
       AED_REAL Temp;            // temperature
       AED_REAL Salinity;        // salinity
       AED_REAL Height;          // 1-D depth array
       AED_REAL MeanHeight;      // Mean depth of a layer
       AED_REAL LayerVol;        // volume of layer
       AED_REAL LayerArea;       // area of layer

       AED_REAL Light;           // solar radiation over water layer depths
       AED_REAL ExtcCoefSW;      // light extinction coefficient

       AED_REAL Vol1;            // Cumulative volume to this layer top
       AED_REAL Epsilon;         // Diffusivity

       AED_REAL Umean;           // Mean velocity
       AED_REAL Uorb;            // Maximum orbital velocity
       AED_REAL LayerStress;     // Layer Stress
   } LakeDataType;

   /*===========================================================*/
   // Structured type for Met vars
   typedef struct MetDataType {
       AED_REAL Rain;            // raindfall
       AED_REAL RelHum;          // relative humidty
       AED_REAL SatVapDef;       // vapour pressure
       AED_REAL LongWave;        // longwave radiation
       AED_REAL ShortWave;       // shortwave radiation
       AED_REAL AirTemp;         // temperature
       AED_REAL WindSpeed;       // windspeed
       AED_REAL Snow;            // snowdfall
       AED_REAL RainConcPO4;     // Concentration of PO4 in rain
       AED_REAL RainConcTp;      // Concentration of TP in rain
       AED_REAL RainConcNO3;     // Concentration of NO3 in rain
       AED_REAL RainConcNH4;     // Concentration of NH4 in rain
       AED_REAL RainConcTn;      // Concentration of TN in rain
       AED_REAL RainConcSi;      // Concentration of SI in rain
       AED_REAL WindDir;         // Wind direction
       AED_REAL As;              // Area of sheltering
   } MetDataType;

   /*===========================================================*/
   // Structured type for Surface Data vars
   typedef struct SurfaceDataType {
       AED_REAL Evap;            // Evaporation
       AED_REAL delzBlueIce;     // Thickness of blue ice layer
       AED_REAL delzWhiteIce;    // Thickness of white ice layer
       AED_REAL delzSnow;        // Thickness of snow layer
       AED_REAL dHt;             // Change in thickness of the snow / ice layer
       AED_REAL RhoSnow;         // Density of snow layer (kg/m^3)
       AED_REAL dailyEvap;       // Daily Evaporation (m3/day)
       AED_REAL dailyRain;       // Daily Rain (m3/day)
       AED_REAL dailyRunoff;     // Daily Runoff (m3/day)
       AED_REAL dailySnow;       // Daily Snow (m3/day)
       AED_REAL dailyQsw;        // Daily Heat Flux (J/day)
       AED_REAL dailyQe;         // Daily Latent Heat(J/day)
       AED_REAL dailyQh;         // Daily Sensible Heat (J/day)
       AED_REAL dailyQlw;        // Daily Long Wave Radiation (J/day)
       AED_REAL dailyInflow;     // Total Daily Inflow (m3/day)
       AED_REAL dailyOutflow;    // Total Daily Outflow (m3/day)
       AED_REAL dailyOverflow;   // Total Daily Overflow (m3/day)
       AED_REAL albedo;          // Daily surface albedo
       AED_REAL dailyzonL;       // Daily atmospheric stability
   } SurfaceDataType;

   /*===========================================================*/
   // Structured type for Sediment Layer
   typedef struct  SedLayerType {
      AED_REAL depth;            // Layer height
      AED_REAL temp;             // Layer temperature
      AED_REAL vwc;
      AED_REAL wq;
   } SedLayerType;

   /*===========================================================*/
   // Structured type for iSediment Zones
   typedef struct ZoneType {
       AED_REAL zheight;
       AED_REAL zrad;
       AED_REAL zsalt;
       AED_REAL ztemp;
       AED_REAL zrho;
       AED_REAL zarea;
       AED_REAL zextc_coef;
       AED_REAL zlayer_stress;
       AED_REAL ztss;
       AED_REAL zdz;
       AED_REAL zpar;
       AED_REAL znir;
       AED_REAL zuva;
       AED_REAL zuvb;
       AED_REAL zpres;
       AED_REAL zdepth;
       AED_REAL z_sed_zones;
       AED_REAL z_pc_wet;
       AED_REAL heatflux;
       int n_sedLayers;      // number of sediment layers
       SedLayerType *layers;
   } ZoneType;

/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#endif
