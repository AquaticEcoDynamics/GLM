/******************************************************************************
 *                                                                            *
 * glm_ptm.c                                                                  *
 *                                                                            *
 * Contains Particle Tracking Model                                           *
 *                                                                            *
 * Developed by :                                                             *
 *                                                                            *
 * Copyright 2024 - The University of Western Australia                       *
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
#include <time.h>       /* time */

#include "glm.h"

#include "glm_types.h"      // check if these are needed
#include "glm_const.h"
#include "glm_globals.h"

#include "glm_ptm.h"

#include "glm_util.h"
#include "glm_ncdf.h"

            // Cell suspended in water column
#define WATER  0
            // Maybe forming a bed layer; set when hits bottom; turn off when resuspended
#define BED    1
            // Maybe forming a scum layer
#define SCUM   2
            // Maybe exiting the lake
#define EXIT   3

AED_REAL get_settling_velocity(AED_REAL settling_velocity);
AED_REAL random_walk(AED_REAL dt, AED_REAL Height, AED_REAL K_z, AED_REAL K_prime_z, AED_REAL vvel);

/*============================================================================*/

//CONSTANTS
//int max_particle_num;  // replace these from namelist
int num_particle_grp=1;
//int init_particle_num;
AED_REAL init_depth_min=0.0;
AED_REAL init_depth_max=2.0;
AED_REAL ptm_time_step=1/60;
AED_REAL ptm_diffusivity=1e-6;
//AED_REAL settling_velocity;
//AED_REAL settling_efficiency;

// VARIABLES
LOGICAL sed_deactivation = FALSE;


/*============================================================================*/


/******************************************************************************
 *                                                                            *
 *    This routine initialises the PTM data structure, and sets the initial   *
 *    particle properties and position                                        *
 *                                                                            *
 *                                                                            *
 *                                                                            *
 *    NOTES:                                                                  *
 *      all arrays begin with subscript of zero                               *
 *                                                                            *
 ******************************************************************************/
void ptm_init_glm()
{
//LOCALS
    int p;

    AED_REAL upper_height;
    AED_REAL lower_height;

/*----------------------------------------------------------------------------*/
//BEGIN

    // Allocate maximum number of particles
    Particle = calloc(max_particle_num, sizeof(ParticleDataType));

    // Set initial inactive particle status/properties (initial inactive particles)
    for (p = 0; p < max_particle_num; p++) {
        //partgroup[grp].istat[id_stat,p] id_stat is equivalent to Status, p is particle
        //partgroup[grp].istat[id_flag,p] id_flag is equivalent to Flag, p is particle
        Particle[p].Status = 0;
        Particle[p].Flag = 3;
        Particle[p].Mass = 0.0;
        Particle[p].Diam = 0.0;
        Particle[p].Density = 0.0;
        Particle[p].vvel = 0.0;
    }

    // Set initial active particle height within the water column
    upper_height = Lake[surfLayer].Height - init_depth_min;
    lower_height = Lake[surfLayer].Height - init_depth_max;

    ptm_addparticles(init_particle_num, max_particle_num, upper_height,lower_height);

    ptm_update_layerid();     // assign layers to active particles
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *    This routine updates the particle positions, based on random motions    *
 *    and settling/migration                                                  *
 *                                                                            *
 ******************************************************************************/
void do_ptm_update()
{
//LOCALS
    int p, tt, ij1, ij2, sub_steps;
    AED_REAL dt, K_z, K_above, K_prime_z;
    float rand_float, prob, prev_height, x1, x2, y1, y2, a1, a2;

/*----------------------------------------------------------------------------*/
//BEGIN

    // Update settling/migration velocity  ! Will overwrite AED
    for (p = 0; p < max_particle_num; p++) {
        if (Particle[p].Status>0) {
           Particle[p].vvel = get_settling_velocity(settling_velocity);
        }
    }

    // Loop through sub-timesteps, incrementing position
    sub_steps = 60;
    dt = 1;
    for (tt = 1; tt < sub_steps; tt++) {
        for (p = 0; p < max_particle_num; p++) {
          if (Particle[p].Status>0) {
            // Capture current height of particle to calculate probability of settling below
            prev_height = Particle[p].Height;

            // Update particle position based on diffusivity and vert velocity
            Particle[p].Flag= WATER;
            K_z = Lake[Particle[p].Layer].Epsilon;
            K_above = Lake[Particle[p+1].Layer].Epsilon;

            // determine whether to assume molecular diffusion K
            if(K_z < 1E-6){
                K_z = 1E-6;
            }

            if(K_above < 1E-6){
                K_above = 1E-6;
            }

            if(Particle[p].Layer == surfLayer){
                K_prime_z = 0;
                continue;
            }
            else {
                K_prime_z = fabs(K_z - K_above);
            }

            Particle[p].Height = random_walk(dt,Particle[p].Height, K_z, K_prime_z, Particle[p].vvel);

            if(prev_height > Particle[p].Height){

                // Get area at previous particle height
                x1 = prev_height * 10.0;
                y1 = x1 - (int)(x1 / 1.0) * 1.0;
                ij1 = (int)(x1 - y1) - 1;
                if(ij1 > Nmorph){
                    y1 = y1 + (float)(ij1 - Nmorph);
                    ij1 = Nmorph - 1;
                }
                a1 = MphLevelArea[ij1] + y1 * dMphLevelArea[ij1];

                // Get area at depth of current particle height
                x2 = Particle[p].Height * 10.0;
                y2 = x2 - (int)(x2 / 1.0) * 1.0;
                ij2 = (int)(x2 - y2) - 1;
                if(ij2 > Nmorph){
                    y2 = y2 + (float)(ij2 - Nmorph);
                    ij2 = Nmorph - 1;
                }
                a2 = MphLevelArea[ij2] + y2 * dMphLevelArea[ij2];

                // Calculate proportional difference between two areas
                prob = (a1 - a2) / a1 * settling_efficiency;

                // Bernoulli draw to determine if particle should be assigned as BED
                rand_float = ((float)rand())/RAND_MAX;
                if(rand_float < prob || Particle[p].Height < 0.02){
                    Particle[p].Flag= BED;
                    if(Particle[p].Height<0.0){
                        Particle[p].Height=0.0;
                    }
                    if(sed_deactivation){
                        Particle[p].Status = 0;
                    }
                }
            }

            // Determine if particle should be assigned as SCUM
            if(Particle[p].Height>Lake[surfLayer].Height){
                Particle[p].Flag= SCUM;                      // Maybe forming a scum layer
                Particle[p].Height=Lake[surfLayer].Height;
            }
          }
        }
    }

    ptm_update_layerid();     // assign layers to active particles
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *    This routine redistributes any particles within the provided depth      *
 *    range, used to capture the effect of layer mixing                       *
 *                                                                            *
 ******************************************************************************/
void ptm_redistribute(AED_REAL upper_height, AED_REAL lower_height)
{
//LOCALS
    int p;

    int rand_int;
    AED_REAL height_range;

/*----------------------------------------------------------------------------*/
//BEGIN
    // Get vertical range in the water column that mixed
    height_range = upper_height - lower_height;

    // Check for active particles in the height range
    for (p = 0; p < max_particle_num; p++) {
        if (Particle[p].Status>0) {
            if (Particle[p].Height>=lower_height && Particle[p].Height<=upper_height ) {
                // Particle is in the mixing zone, so re-position
                rand_int = rand() % 100 + 1;                            // random draw from unit distribution
                double random_double = (double)rand_int / 100;
                random_double = random_double * height_range;                 // scale unit random to requested range
                Particle[p].Height = lower_height + random_double;
            }
        }
    }
    ptm_update_layerid();     // assign layers to active particles
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *    This routine adds new particles within the provided depth               *
 *    range, used to capture the effect of layer mixing                       *
 *                                                                            *
 *                                                                            *
 ******************************************************************************/
void ptm_addparticles(int new_particles, int max_particle_num, AED_REAL upper_height,
                      AED_REAL lower_height)
{
//LOCALS
    int p, n;

    int rand_int;
    AED_REAL height_range;

/*----------------------------------------------------------------------------*/
//BEGIN
    // Get vertical range in the water column that mixed
    height_range = upper_height - lower_height;
    n = 0;

    // For each new particle, initialise their properties and height
    for (p = 0 ; p < max_particle_num; p++) {
        if(n == new_particles){
            break;
        }
        if(Particle[p].Status == 0 && Particle[p].Flag == 3){ // find the first inactive particles with EXIT flag
            Particle[p].Status = 1;
            Particle[p].Flag = 0;
            Particle[p].Mass = 1.0;
            Particle[p].Diam = 1e-6;
            Particle[p].Density = 1000.0;
            Particle[p].vvel = 0.0;

            // Assign particles initial height
            rand_int = rand() % 100 + 1;                            // random draw from unit distribution
            double random_double = (double)rand_int / 100;
            random_double = random_double * height_range;                 // scale unit random to requested range
            Particle[p].Height = lower_height + random_double;   // set particle height

            // Adjust counter
            n++;
        }
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

/******************************************************************************
 *                                                                            *
 *    This routine removes new particles from a specified layer, based        *
 *    on the proportion of layer volume that is removed through outflow       *
 *                                                                            *
 *                                                                            *
 ******************************************************************************/
void ptm_removeparticles(int layer_id, AED_REAL delta_vol, AED_REAL layer_vol, int max_particle_num)
{
//LOCALS
    AED_REAL layer_prop, rand_float;
    int p;

/*----------------------------------------------------------------------------*/
//BEGIN
    // For each particle, draw from Bernoulli distribution to see whether removed from layer
    layer_prop = delta_vol / layer_vol;
    for (p = 0; p < max_particle_num; p++) {
        if(Particle[p].Status == 1 && Particle[p].Layer == layer_id){
            rand_float = ((float)rand())/RAND_MAX;
            if(rand_float <= layer_prop){
                // If particle leaves through outflow, reset completely
                Particle[p].Status = 0;
                Particle[p].Flag = 3;
                Particle[p].Mass = 0.0;
                Particle[p].Diam = 0.0;
                Particle[p].Density = 0.0;
                Particle[p].vvel = 0.0;
            }
        }
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *    This routine raises or lowers particle vertical positions due to        *
 *    changes in layering associated with inflows/outflows                    *
 *                                                                            *
 *                                                                            *
 ******************************************************************************/
void ptm_layershift(AED_REAL shift_height, AED_REAL shift_amount)
{
//LOCALS
    int p;

    AED_REAL upper_height;
    AED_REAL lower_height;
//  AED_REAL height_range;

/*----------------------------------------------------------------------------*/
//BEGIN

    // Get vertical range in the water column that is impacted by the shift
    lower_height = shift_height;
    upper_height = shift_height + shift_amount;
//  height_range = upper_height - lower_height;

    // Check for active particles in the impacted height range
    for (p = 0; p < max_particle_num; p++) {
        if (Particle[p].Status>0) {
            if (Particle[p].Height>lower_height && Particle[p].Height<upper_height ) {
                // Particle is in the impacted zone, so re-position (lift or drop)
                Particle[p].Height = Particle[p].Height + shift_amount;   // adjust particle height
            }
        }
    }

    ptm_update_layerid();     // assign layers to active particles
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *    This routine sets the layer id for each particle based on its height    *
 *                                                                            *
 ******************************************************************************/
void ptm_update_layerid()
{
//LOCALS
    int p,i;

/*----------------------------------------------------------------------------*/
//BEGIN

    for (p = 0; p < max_particle_num; p++) {
        if (Particle[p].Status>0) {
            for (i = botmLayer; i < NumLayers; i++) {
                if (Particle[p].Height<Lake[i].Height) {
                    Particle[p].Layer = i;
                    break; // get out of layer loop
                }
            }
        }
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *      This routine redistributes particles using a random walk function     *
 *                                                                            *
 ******************************************************************************/
AED_REAL random_walk(AED_REAL dt, AED_REAL Height, AED_REAL K_z, AED_REAL K_prime_z, AED_REAL vvel)
{
//LOCALS

    AED_REAL updated_height;
    AED_REAL del_t;
    float random_float;

/*----------------------------------------------------------------------------*/
//BEGIN

    del_t = dt*60;

    random_float = -1+2*((float)rand())/RAND_MAX;            // random draw from uniform distribution [-1,1]

    updated_height = Height + K_prime_z * Height * del_t + random_float *
    sqrt((2 * K_z * (Height + 0.5 * K_prime_z * Height * del_t) * del_t) / (1.0/3)); // random walk

    updated_height = updated_height + vvel;                   // account for sinking/floating

    return updated_height;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


static int h_id, m_id, d_id, dn_id, vv_id, stat_id, flag_id;
static int set_no_p = -1;
static size_t start[2],edges[2];


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
//void ptm_write_glm(int *ncid, int *wlev, int *nlev, int *lvl, int *point_nlevs)
void ptm_write_glm(int ncid, int max_particle_num)
{
//LOCALS
    int p;
    AED_REAL *p_height, *mass, *diam, *density, *vvel;
    int *status, *flag;

/*----------------------------------------------------------------------------*/
//BEGIN


    set_no_p++;

    start[1] = 0;             edges[1] = max_particle_num;
    start[0] = set_no_p;      edges[0] = 1;

    p_height  = malloc(max_particle_num*sizeof(AED_REAL));
    mass  = malloc(max_particle_num*sizeof(AED_REAL));
    diam  = malloc(max_particle_num*sizeof(AED_REAL));
    density  = malloc(max_particle_num*sizeof(AED_REAL));
    vvel  = malloc(max_particle_num*sizeof(AED_REAL));

    status  = malloc(max_particle_num*sizeof(int));
    flag  = malloc(max_particle_num*sizeof(int));

    for (p = 0; p < max_particle_num; p++) {
		p_height[p] = Particle[p].Height;
        mass[p] = Particle[p].Mass;
        diam[p] = Particle[p].Diam;
        density[p] = Particle[p].Density;
        vvel[p] = Particle[p].vvel;
        status[p] = Particle[p].Status;
        flag[p] = Particle[p].Flag;
    }

     nc_put_vara(ncid, h_id, start, edges, p_height);
     nc_put_vara(ncid, m_id, start, edges, mass);
     nc_put_vara(ncid, d_id, start, edges, diam);
     nc_put_vara(ncid, dn_id, start, edges, density);
     nc_put_vara(ncid, vv_id, start, edges, vvel);
     nc_put_vara(ncid, stat_id, start, edges, status);
     nc_put_vara(ncid, flag_id, start, edges, flag);

    free(p_height);
    free(mass);
    free(diam);
    free(density);
    free(vvel);
    free(status);
    free(flag);

    check_nc_error(nc_sync(ncid));
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
//will also need to handle groups in write step; append group name onto
void ptm_init_glm_output(int ncid, int time_dim)
{
   int dims[2];

//
//------------------------------------------------------------------------------
//BEGIN
   define_mode_on(&ncid);   // Put NetCDF library in define mode.

   check_nc_error(nc_def_dim(ncid, "particles", max_particle_num, &ptm_dim));

   dims[1] = ptm_dim;
   dims[0] = time_dim;

   check_nc_error(nc_def_var(ncid, "particle_height", NC_REALTYPE, 2, dims, &h_id));
   set_nc_attributes(ncid, h_id, "meters", "Height of Particle" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_mass", NC_REALTYPE, 2, dims, &m_id));
   set_nc_attributes(ncid, m_id, "grams", "Mass of Particle" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_diameter", NC_REALTYPE, 2, dims, &d_id));
   set_nc_attributes(ncid, d_id, "meters", "Diameter of Particle" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_density", NC_REALTYPE, 2, dims, &dn_id));
   set_nc_attributes(ncid, dn_id, "g/m3", "Density of Particle" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_vvel", NC_REALTYPE, 2, dims, &vv_id));
   set_nc_attributes(ncid, vv_id, "m/s", "Settling Velocity of Particle" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_status", NC_INT, 2, dims, &stat_id));
   nc_put_att(ncid, stat_id, "long_name", NC_CHAR, 18, "Status of Particle");

   check_nc_error(nc_def_var(ncid, "particle_flag", NC_INT, 2, dims, &flag_id));
   nc_put_att(ncid, flag_id, "long_name", NC_CHAR, 18, "Location Flag of Particle");

   define_mode_off(&ncid);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *        This routine returns the settling velocity for a particle           *
 *                                                                            *
 ******************************************************************************/
AED_REAL get_settling_velocity(AED_REAL settling_velocity)
{
    return settling_velocity;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

