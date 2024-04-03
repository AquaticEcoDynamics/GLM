/******************************************************************************
 *                                                                            *
 * glm_ptm.h                                                                  *
 *                                                                            *
 * Contains Particle Tracking Model                                           *
 *                                                                            *
 * Developed by :                                                             *
 *                                                                            *
 * Copyright 2024 - 2024 -  The University of Western Australia               *
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

#include "glm_types.h"      // check if these are needed
#include "glm_const.h"
#include "glm_globals.h"

#include "glm_ptm.h"

#include "glm_util.h"


            // Maybe forming a bed layer
#define BED    1
            // Maybe forming a scum layer
#define SCUM   2


AED_REAL settling_velocity(void);
void random_walk(AED_REAL dt, AED_REAL Height, AED_REAL Epsilon, AED_REAL vvel);

/*============================================================================*/

//CONSTANTS
int max_particle_num=1000000;  // replace these from namelist
int num_particle_grp=1;
int init_particle_num=10;
AED_REAL init_depth_min=0.0;
AED_REAL init_depth_max=2.0;
AED_REAL ptm_time_step=1/60;
AED_REAL ptm_diffusivity=1e-6;

// VARIABLES
int ptm_sw = FALSE;

int last_particle = 0;
int num_particles = 0;

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
//  int i, j;
    int p;
//  int flag;

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
        Particle[p].Mass = 1.0;
        Particle[p].Diam = 1e-6;
        Particle[p].Density = 1000.0;
        Particle[p].vvel = 0.0;
    }

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
    int p, tt;
    int sub_steps;

//  AED_REAL upper_height;
//  AED_REAL lower_height;
    AED_REAL dt;

/*----------------------------------------------------------------------------*/
//BEGIN

   
    // Update settling/migration velocity  ! Will overwrite AED
    for (p = 0; p < num_particles; p++) { 
      Particle[p].vvel = settling_velocity();
    }
    
    // Loop through sub-timesteps, incrementing position
    sub_steps = 60;
    dt = 1;
    for (tt = 1; tt < sub_steps; tt++) { 
        for (p = 0; p < num_particles; p++) { 
          if (Particle[p].Status>0) {
            // Update particle position based on diffusivity and vert velocity
            random_walk(dt,Particle[p].Height, Lake[Particle[p].Layer].Epsilon,Particle[p].vvel);

            // Mary random walk equation in R:
            // following Visser 1997 https://www.int-res.com/articles/meps/158/m158p275.pdf
            // z_t1 <- z + K_prime_z * z * del_t + runif(1, min = -1,max = 1) * sqrt((2 * K_z * (z + 0.5 * K_prime_z * z * del_t) * del_t) / (1/3))
            // z is depth (or elevation in this case)
            // del_t is timestep of one minute in SECONDS
            // K_z is vertical diffusivity at depth/elevation z; we are assuming constant K of 1x10e-6 for now I think?
            // K_prime is delta K / delta z
            
            // Check if status change due to hitting bed, or surface
            if(Particle[p].Height<0){
                Particle[p].Status=BED ;                        // Maybe forming a bed layer
                Particle[p].Height=0.0 ; 
            }
            if(Particle[p].Height>Lake[surfLayer].Height){
                Particle[p].Status=SCUM  ;                      // Maybe forming a scum layer
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
 *                                                                            *
 ******************************************************************************/
void ptm_redistribute(AED_REAL upper_height, AED_REAL lower_height)
{
//LOCALS
    int p;

    AED_REAL rand;
    AED_REAL height_range;

/*----------------------------------------------------------------------------*/
//BEGIN

    // Get vertical range in the water column that mixed
    height_range = upper_height - lower_height;

    // Check for active particles in the height range
    for (p = 0; p < last_particle; p++) { 
        if (Particle[p].Status>0) {
          if (Particle[p].Height>lower_height && Particle[p].Height<upper_height ) {
            // Particle is in the mixing zone, so re-position
            rand = random();                            // random draw from unit distribution
            rand = rand * height_range;                 // scale unit random to requested range
            Particle[p].Height = lower_height + rand;   // reset particle height
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

    AED_REAL rand;
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
        rand = random();                            // random draw from unit distribution
        rand = rand * height_range;                 // scale unit random to requested range
        Particle[p].Height = lower_height + rand;   // set particle height
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
void random_walk(AED_REAL dt, AED_REAL Height, AED_REAL Epsilon, AED_REAL vvel)
{
}

//void ptm_write_glm(int *ncid, int *wlev, int *nlev, int *lvl, int *point_nlevs)
void ptm_write_glm()
{
}

void ptm_init_glm_output(int *ncid, int *time_dim)
{
}

AED_REAL settling_velocity()
{
    return 0.01;
}
