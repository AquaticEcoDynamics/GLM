/******************************************************************************
 *                                                                            *
 * glm_util.c                                                                 *
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

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <math.h>

#include "glm.h"
#include "glm_types.h"
#include "glm_const.h"
#include "glm_globals.h"

/******************************************************************************
 * Calculate the Zenith Angle in degrees.                                     *
 * NB lat is in radians, lon in degrees [thats how they are in GLM]           *
 ******************************************************************************/
AED_REAL zenith_angle(AED_REAL lon, AED_REAL lat, int day, int iclock, AED_REAL TZ)
{
    AED_REAL Phi_day;          // Day Angle :- Position of the earth in sun's orbit
    AED_REAL SolarDeclination; // Solar Declination
    AED_REAL EquationOfTime;   // Equation of Time
    AED_REAL Phi_hr;           // Hour Angle

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    AED_REAL Hour = iclock / 3600.;

    /*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

    // Day Angle :- Position of the earth in sun's orbit
    Phi_day = (two_Pi*(day-1)/365);

    // Solar Declination (in radians)
    SolarDeclination = ( 0.006918 - 0.399912*cos(  Phi_day) +
                                    0.070257*sin(  Phi_day) -
                                    0.006758*cos(2*Phi_day) +
                                    0.000907*sin(2*Phi_day) -
                                    0.002697*cos(3*Phi_day) +
                                    0.00148 *sin(3*Phi_day) );

    // Equation of Time
    EquationOfTime = ( 0.0000075 + 0.001868*cos(  Phi_day) -
                                   0.032077*sin(  Phi_day) -
                                   0.014615*cos(2*Phi_day) -
                                   0.040849*sin(2*Phi_day) ) * 229.18;

    // Hour Angle (in degrees)
    Phi_hr = 15 * (Hour-12.5) + lon - TZ * 15 + (EquationOfTime/4);

    // Zenith Angle
    return acos( cos(SolarDeclination) * cos(lat) * cos(Phi_hr*deg2rad) +
                 sin(SolarDeclination) * sin(lat) ) * rad2deg;
}

#ifndef sqr
/******************************************************************************/
// This function computes the square of the argument
AED_REAL sqr(AED_REAL x) { return x*x; }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#endif


#ifndef gprime
/******************************************************************************
 * This function calculates reduced gravity (gprime) given two densities      *
 ******************************************************************************/
AED_REAL gprime(AED_REAL rho1, AED_REAL rho2)
//{ return g * (rho2 - rho1) / ((rho1 + rho2) / 2.0); }
{
#define GTOLERANCE  1E-7
    AED_REAL tmp = g * (rho2 - rho1) / ((rho1 + rho2) / 2.0);
    if (tmp < 0)
        fprintf(stderr, " gprime returns < 0 ( = %.12e) rho1 %.12e rho2 %.12e\n", tmp, rho1, rho2);
    if (fabs(tmp) < GTOLERANCE) {
        if (tmp < 0) return -GTOLERANCE;
        return GTOLERANCE;
    }
    return tmp;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#endif



/******************************************************************************
 * This function combines two layers and return the mean concentration        *
 ******************************************************************************/
AED_REAL combine(AED_REAL c1, AED_REAL v1, AED_REAL d1,
                 AED_REAL c2, AED_REAL v2, AED_REAL d2)
{
    AED_REAL MTotal;

/*----------------------------------------------------------------------------*/
    if (fabs(c1-c2) < 1e-5 && fabs(d1-d2) < 1e-5) return c1;

    MTotal = v1 * d1 + v2 * d2;

    return (c1 * v1 * d1 + c2 * v2 * d2) / MTotal;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


#ifndef combine_vol
/******************************************************************************
 * Function to combine two layers and return the mean                         *
 * concentration taking into account volumes                                  *
 ******************************************************************************/
AED_REAL combine_vol(AED_REAL c1,AED_REAL v1,AED_REAL c2,AED_REAL v2)
{ return (c1 * v1 + c2 * v2) / (v1 + v2); }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#endif


/******************************************************************************
 *                                                                            *
 * Function to calculate the density (rho) of water at a given temperature    *
 * (deg C) and salinity (ppm) based on UNESCO (1981) polynomial               *
 *                                                                            *
 ******************************************************************************/
static AED_REAL calculate_density_UNESCO(AED_REAL temp, AED_REAL salt)
{
    AED_REAL dpure;
    AED_REAL t1,t2,t3,t4,t5,tm;
    AED_REAL s1,s32,s2;
    AED_REAL csal1;
    AED_REAL csal32;
    AED_REAL csal2;

    AED_REAL term[15];

    const AED_REAL c1=999.842594,  c2=6.793952E-2, c3=9.095290E-3,
                   c4=1.001685E-4, c5=1.120083E-6, c6=6.536332E-9,
                   d1=8.24493E-1,  d2=4.0899E-3,   d3=7.6438E-5,
                   d4=8.2467E-7,   d5=5.3875E-9,   d6=5.72466E-3,
                   d7=1.0227E-4,   d8=1.6546E-6,   d9=4.8314E-4;

    t1 = (temp);
    s1 = (salt);
    tm = round(t1*10000.);
    t1 = (tm)/10000.;

    t2 = sqr(t1);
    t3 = t2*t1;
    t4 = t3*t1;
    t5 = t4*t1;
    s2 = s1*s1;
    s32 = pow(s1,1.5);

    term[0]  =  c1;
    term[1]  =  c2 * t1;
    term[2]  = -c3 * t2;
    term[3]  =  c4 * t3;
    term[4]  = -c5 * t4;
    term[5]  =  c6 * t5;
    term[6]  =  d1;
    term[7]  = -d2 * t1;
    term[8]  =  d3 * t2;
    term[9]  = -d4 * t3;
    term[10] =  d5 * t4;
    term[11] = -d6;
    term[12] =  d7 * t1;
    term[13] = -d8 * t2;
    term[14] =  d9;

    dpure  =  term[5]  + term[4]  + term[3]  + term[1] + term[2] + term[0];
    csal1  = (term[10] + term[9]  + term[8]  + term[7] + term[6]) * s1;
    csal32 = (term[13] + term[12] + term[11]) * s32;
    csal2  =  term[14] * s2;

    return dpure + csal1 + csal32 + csal2;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 * Calculates specific volume from Absolute Salinity, Conservative            *
 * Temperature and pressure, using the computationally-efficient              *
 * polynomial expression for specific volume (Roquet et al., 2014).           *
 *                                                                            *
 * sa     : Absolute Salinity                               [g/kg]            *
 * ct     : Conservative Temperature (ITS-90)               [deg C]           *
 * p      : sea pressure                                    [dbar]            *
 *          ( i.e. absolute pressure - 10.1325 dbar )                         *
 *                                                                            *
 * specvol: specific volume                                 [m^3/kg]          *
 *                                                                            *
 * This code was taken from the "Gibbs-SeaWater (GSW) Oceanographic Toolbox"  *
 *                                                                            *
 *   ---------------------------------------------------------------------    *
 *                                                                            *
 * Licence for the use of the Gibbs SeaWater (GSW) Oceanographic Toolbox      *
 *                                                                            *
 * Copyright (c) 2011, SCOR/IAPSO WG127 (Scientific Committee on Oceanic      *
 * Research/ International Association for the Physical Sciences of the       *
 * Oceans, Working Group 127).                                                *
 *                                                                            *
 * All rights reserved.                                                       *
 *                                                                            *
 * Redistribution and use, in source and binary forms, without modification,  *
 * is permitted provided that the following conditions are met:               *
 *                                                                            *
�*� * Redistributions of source code must retain the above copyright notice,  *
 *    this list of conditions and the following disclaimer.                   *
 *                                                                            *
�*� * Redistributions in binary form must reproduce the above copyright       *
 *    notice, this list of conditions and the following disclaimer in the     *
 *    documentation and/or other materials provided with the distribution.    *
 *                                                                            *
�*  * Neither the name of SCOR/IAPSO WG127 nor the names of its contributors  *
 *    may be used to endorse or promote products derived from this software   *
 *    without specific prior written permission.                              *
 *                                                                            *
 *                                                                            *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"*
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE  *
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE *
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE  *
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR        *
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF       *
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS   *
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN    *
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE *
 * POSSIBILITY OF SUCH DAMAGE.                                                *
 *                                                                            *
 * The software is available from http://www.TEOS-10.org                      *
 *                                                                            *
 ******************************************************************************/
#include "gsw_const.h"

static AED_REAL gsw_specvol(AED_REAL sa, AED_REAL ct, AED_REAL p)
{
    GSW_TEOS10_CONSTANTS;
    GSW_SPECVOL_COEFFICIENTS;
    AED_REAL  xs, ys, z, value;

    xs      = sqrt(gsw_sfac*sa + offset);
    ys      = ct*0.025;
    z       = p*1e-4;

    value = v000
        + xs*(v010 + xs*(v020 + xs*(v030 + xs*(v040 + xs*(v050
        + v060*xs))))) + ys*(v100 + xs*(v110 + xs*(v120 + xs*(v130 + xs*(v140
        + v150*xs)))) + ys*(v200 + xs*(v210 + xs*(v220 + xs*(v230 + v240*xs)))
        + ys*(v300 + xs*(v310 + xs*(v320 + v330*xs)) + ys*(v400 + xs*(v410
        + v420*xs) + ys*(v500 + v510*xs + v600*ys))))) + z*(v001 + xs*(v011
        + xs*(v021 + xs*(v031 + xs*(v041 + v051*xs)))) + ys*(v101 + xs*(v111
        + xs*(v121 + xs*(v131 + v141*xs))) + ys*(v201 + xs*(v211 + xs*(v221
        + v231*xs)) + ys*(v301 + xs*(v311 + v321*xs) + ys*(v401 + v411*xs
        + v501*ys)))) + z*(v002 + xs*(v012 + xs*(v022 + xs*(v032 + v042*xs)))
        + ys*(v102 + xs*(v112 + xs*(v122 + v132*xs)) + ys*(v202 + xs*(v212
        + v222*xs) + ys*(v302 + v312*xs + v402*ys))) + z*(v003 + xs*(v013
        + v023*xs) + ys*(v103 + v113*xs + v203*ys) + z*(v004 + v014*xs + v104*ys
        + z*(v005 + v006*z)))));

    return (value);
}
/******************************************************************************
 ******************************************************************************/
static AED_REAL calculate_density_TEOS(AED_REAL temp, AED_REAL salt)
{
    return (1.0/gsw_specvol(salt,temp,atm_pressure_sl/100.));
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 * Function to calculate the density (rho) of water at a given temperature    *
 * (deg C) and salinity (ppm)                                                 *
 *                                                                            *
 ******************************************************************************/
AED_REAL calculate_density(AED_REAL temp, AED_REAL salt)
{
    switch (density_model) {
        case 0 : return calculate_density_UNESCO(temp, salt);
        case 1 : return calculate_density_TEOS(temp, salt);
//      case 2 => custom etc .
    }
    return MISVAL;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 * Calculates the saturated vapour pressure (SatVap) corresponding to         *
 * temperature (temp, deg C)                                                  *
 ******************************************************************************/
AED_REAL saturated_vapour(AED_REAL temp)
{ return pow(10.0, 9.28603523 - (2322.37885/(temp+Kelvin))); }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************/
int internal_var(const char *name)
{
    if ( strncasecmp(name, "flow", 4) == 0 ||
         strncasecmp(name, "temp", 4) == 0 ||
         strncasecmp(name, "salt", 4) == 0 ||
         strncasecmp(name, "flbc", 4) == 0 ||
         strncasecmp(name, "htop", 4) == 0 ||
         strncasecmp(name, "hbot", 4) == 0 )
        return 1;
    return 0;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

#ifndef PLOTS
#ifdef _WIN32
/******************************************************************************/
char *strndup(const char *s, size_t len)
{
    size_t l = strlen(s);
    char *t = malloc(min(l, len) + 1);
    if (t != NULL) { strncpy(t, s, min(l, len) + 1); t[min(l, len)] = 0; }
    return t;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
char *ctime_r(const time_t *timep, char *buf)
{
    char *s = ctime(timep);
    strcpy(buf, s);
    return buf;
}

#endif
#endif
