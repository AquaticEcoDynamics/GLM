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


AED_REAL get_settling_velocity(AED_REAL settling_velocity);
AED_REAL random_walk(AED_REAL dt, AED_REAL Height, AED_REAL Epsilon, AED_REAL vvel);

/*============================================================================*/

//CONSTANTS
int max_particle_num=1000000;  // replace these from namelist
int num_particle_grp=1;
int init_particle_num;
AED_REAL init_depth_min=0.0;
AED_REAL init_depth_max=2.0;
AED_REAL ptm_time_step=1/60;
AED_REAL ptm_diffusivity=1e-6;
<<<<<<< HEAD
AED_REAL settling_velocity=0.;
=======
AED_REAL settling_velocity;
AED_REAL settling_efficiency;
>>>>>>> settling on flanks and epsilon

// VARIABLES
int ptm_sw = FALSE;
int sed_deactivation = FALSE;

int last_particle = 0;

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

    last_particle = 0;

    // Allocate maximum number of particles
    Particle = calloc(max_particle_num, sizeof(ParticleDataType));

    // Pre-initialise all particle status index to 0 (all inactive)
    for (p = 0; p < max_particle_num; p++) Particle[p].Status = 0;

    // Set initial active particle status/properties (initial active particles)
    for (p = 0; p < init_particle_num; p++) {
        Particle[p].Status = 1;
        Particle[p].Flag = 0;
        Particle[p].Mass = 1.0;
        Particle[p].Diam = 1e-6;
        Particle[p].Density = 1000.0;
        Particle[p].vvel = 0.0;
    }

    num_particles = init_particle_num;

    // Set initial active particle height within the water column
    upper_height = Lake[surfLayer].Height - init_depth_min;
    lower_height = Lake[surfLayer].Height - init_depth_max;

    ptm_addparticles(init_particle_num,upper_height,lower_height);

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
    int p, tt, z, ij1, ij2, sub_steps;
    AED_REAL dt;
    float rand_float, prob, prev_height, x1, x2, y1, y2, a1, a2;

/*----------------------------------------------------------------------------*/
//BEGIN

    // Update settling/migration velocity  ! Will overwrite AED
    for (p = 0; p < num_particles; p++) {
        if (Particle[p].Status>0) {
           Particle[p].vvel = get_settling_velocity(settling_velocity);
        }
    }

    // Loop through sub-timesteps, incrementing position
    sub_steps = 60;
    dt = 1;
    for (tt = 1; tt < sub_steps; tt++) {
        for (p = 0; p < num_particles; p++) {
          if (Particle[p].Status>0) {
            // Capture current height of particle to calculate probability of settling below
            prev_height = Particle[p].Height;

            // Update particle position based on diffusivity and vert velocity
            Particle[p].Flag= WATER;  
            Particle[p].Height = random_walk(dt,Particle[p].Height, Lake[Particle[p].Layer].Epsilon,Particle[p].vvel);
            //fprintf(stderr, "  p   = %i\n", p);
            //fprintf(stderr, "  Particle[p].Layer   = %i\n", Particle[p].Layer);
            //fprintf(stderr, "  Lake[Particle[p].Layer].Epsilon   = %f\n", Lake[Particle[p].Layer].Epsilon);


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

            if(Particle[p].Flag == BED){
                fprintf(stderr, "  previous height of bed-flagged particle   = %f\n", prev_height);
                fprintf(stderr, "  height of bed-flagged particle   = %f\n", Particle[p].Height);
                fprintf(stderr, "  area at previous height   = %f\n", a1);
                fprintf(stderr, "  area at new height   = %f\n", a2);
                fprintf(stderr, "  prob   = %f\n", prob);
                fprintf(stderr, "  rand_float   = %f\n", rand_float);
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
    for (p = 0; p < last_particle; p++) {
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
void ptm_addparticles(int new_particles, AED_REAL upper_height, AED_REAL lower_height)
{
//LOCALS
    int p;

    int rand_int;
    AED_REAL height_range;

/*----------------------------------------------------------------------------*/
//BEGIN
    // Get vertical range in the water column that mixed
    height_range = upper_height - lower_height;

    // For each new particle, initialise their properties and height
    for (p = last_particle ; p < last_particle+new_particles; p++) {
        Particle[p].Status = 1;
        Particle[p].Mass = 1.0;
        Particle[p].Diam = 1e-6;
        Particle[p].Density = 1000.0;
        Particle[p].vvel = 0.0;

        // Assign particles initial height
        rand_int = rand() % 100 + 1;                            // random draw from unit distribution
        double random_double = (double)rand_int / 100;
        random_double = random_double * height_range;                 // scale unit random to requested range
        Particle[p].Height = lower_height + random_double;   // set particle height
    }

    last_particle = last_particle + new_particles;  // updates the last index number of active particle set
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
    for (p = 0; p < last_particle; p++) {
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

    for (p = 0; p < num_particles; p++) {
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
 ******************************************************************************/
AED_REAL random_walk(AED_REAL dt, AED_REAL Height, AED_REAL Epsilon, AED_REAL vvel)
{
//LOCALS

    AED_REAL updated_height;
    AED_REAL del_t;
    AED_REAL K;
    AED_REAL K_prime_z;
    float random_float;

/*----------------------------------------------------------------------------*/
//BEGIN

    del_t = dt*60;
    K = 1E-6;
    K_prime_z = 0.0;

    random_float = -1+2*((float)rand())/RAND_MAX;            // random draw from uniform distribution [-1,1]

<<<<<<< HEAD
    updated_height = Height + K_prime_z * Height * del_t + random_float *
=======
    // determine whether to use epsilon or K
    if(Epsilon > 1E-6){
        K = Epsilon;
    }

    updated_height = Height + K_prime_z * Height * del_t + random_float * 
>>>>>>> settling on flanks and epsilon
    sqrt((2 * K * (Height + 0.5 * K_prime_z * Height * del_t) * del_t) / (1.0/3)); // random walk

    updated_height = updated_height + vvel;                   // account for sinking/floating

    return updated_height;

    // Mary random walk equation in R:
            // following Visser 1997 https://www.int-res.com/articles/meps/158/m158p275.pdf
            // z_t1 <- z + K_prime_z * z * del_t + runif(1, min = -1,max = 1) * sqrt((2 * K_z * (z + 0.5 * K_prime_z * z * del_t) * del_t) / (1/3))
            // z is depth (or elevation in this case)
            // del_t is timestep of one minute in SECONDS
            // K_z is vertical diffusivity at depth/elevation z; we are assuming constant K of 1x10e-6 for now I think?
            // K_prime is delta K / delta z
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


static int h_id, m_id, d_id, dn_id, vv_id, stat_id, flag_id;
static int set_no_p = -1;
static size_t start[2],edges[2];


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
//void ptm_write_glm(int *ncid, int *wlev, int *nlev, int *lvl, int *point_nlevs)
void ptm_write_glm(int ncid, int num_particles)
{
//LOCALS
    int p;
    AED_REAL *p_height, *mass, *diam, *density, *vvel;
    int *status, *flag;

/*----------------------------------------------------------------------------*/
//BEGIN


    set_no_p++;

    start[1] = 0;             edges[1] = num_particles;
    start[0] = set_no_p;      edges[0] = 1;

    p_height  = malloc(num_particles*sizeof(AED_REAL));
    mass  = malloc(num_particles*sizeof(AED_REAL));
    diam  = malloc(num_particles*sizeof(AED_REAL));
    density  = malloc(num_particles*sizeof(AED_REAL));
    vvel  = malloc(num_particles*sizeof(AED_REAL));

    status  = malloc(num_particles*sizeof(int));
    
    for (p = 0; p < num_particles; p++) { 
		p_height[p] = Particle[p].Height;
        mass[p] = Particle[p].Mass;
        diam[p] = Particle[p].Diam;
        density[p] = Particle[p].Density;
        vvel[p] = Particle[p].vvel;
        status[p] = Particle[p].Status;
    } 
    
     nc_put_vara(ncid, h_id, start, edges, p_height);
     nc_put_vara(ncid, m_id, start, edges, mass);
     nc_put_vara(ncid, d_id, start, edges, diam);
     nc_put_vara(ncid, dn_id, start, edges, density);
     nc_put_vara(ncid, vv_id, start, edges, vvel);
     nc_put_vara(ncid, stat_id, start, edges, status);
     
    free(p_height); 
    free(mass);
    free(diam);
    free(density);
    free(vvel);
    free(status);
    
    check_nc_error(nc_sync(ncid));
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void ptm_init_glm_output(int ncid, int time_dim)
{
   int dims[2];

//
//------------------------------------------------------------------------------
//BEGIN
   define_mode_on(&ncid);   // Put NetCDF library in define mode.

   check_nc_error(nc_def_dim(ncid, "particles", num_particles, &ptm_dim));

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

