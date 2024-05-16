/******************************************************************************
 *                                                                            *
 * glm_ptm.h                                                                  *
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
        Particle[p].Mass = 1.0;
        Particle[p].Diam = 1e-6;
        Particle[p].Density = 1000.0;
        Particle[p].vvel = 0.0;
    }

    // Set initial active particle height within the water column
    upper_height = Lake[surfLayer].Height - init_depth_min;
    lower_height = Lake[surfLayer].Height - init_depth_max;
    fprintf(stderr, "  ptm_init_glm(): Particle[1].vvel Lake[surfLayer].Height   = %f %f\n", Particle[1].vvel, Lake[surfLayer].Height);


    
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
    srand (time(NULL));

    fprintf(stderr, "  ptm_redistribute(): upper %f lower %f range %f\n", upper_height, lower_height, height_range);


    // Check for active particles in the height range
    for (p = 0; p < last_particle; p++) { 
            fprintf(stderr, "  ptm_redistribute(): Particle[p].Status   = %d\n", Particle[p].Status);

        if (Particle[p].Status>0) {
          if (Particle[p].Height>lower_height && Particle[p].Height<upper_height ) {
            // Particle is in the mixing zone, so re-position
            rand_int = rand() % 100 + 1;                            // random draw from unit distribution
            double random_double = (double)rand_int / 100;
            random_double = random_double * height_range;                 // scale unit random to requested range
            Particle[p].Height = lower_height + random_double;
            fprintf(stderr, "  ptm_redistribute(): rand_int random_double   = %d %f\n", rand_int, random_double);
            fprintf(stderr, "  ptm_redistribute(): Particle[p].Height   = %f\n", Particle[p].Height);


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
    srand (time(NULL));

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
        fprintf(stderr, "  ptm_addparticles(): rand_int and random_double   = %d %f\n", rand_int, random_double);
        random_double = random_double * height_range;                 // scale unit random to requested range
        Particle[p].Height = lower_height + random_double;   // set particle height
        fprintf(stderr, "  ptm_addparticles(): Particle[p].Height   = %f %f %f \n", Particle[p].Height, lower_height, height_range);

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
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


static int h_id, m_id, d_id, dn_id, v_id, vv_id;
static size_t start[4],edges[4];
static int set_no = 1;

#if 1
/******************************************************************************
 *                                                                            *
 ******************************************************************************/
//void ptm_write_glm(int *ncid, int *wlev, int *nlev, int *lvl, int *point_nlevs)
void ptm_write_glm(int ncid, int nlev)
{
//LOCALS
    int p;
    AED_REAL *heights;

/*----------------------------------------------------------------------------*/
//BEGIN

    start[0] = set_no; edges[0] = 1;
    start[1] = 0;      edges[1] = max_particle_num;
    heights   = malloc(max_particle_num*sizeof(AED_REAL));

    for (p = 0; p < max_particle_num; p++) { 
        heights[p] = -1.0;
        if (Particle[p].Status>0) {
            heights[p] = Particle[p].Height;
            nc_put_vara(ncid, h_id, start, edges, heights);
        }
    }

    free(heights);
    check_nc_error(nc_sync(ncid));
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void ptm_init_glm_output(int ncid, int time_dim)
{
   int dims[4];
//
//------------------------------------------------------------------------------
//BEGIN
   define_mode_on(&ncid);   // Put NetCDF library in define mode.

   dims[1] = ptm_dim;
   dims[0] = time_dim;

   check_nc_error(nc_def_var(ncid, "Particle_Height", NC_REALTYPE, 2, dims, &h_id));
   set_nc_attributes(ncid, h_id, "meters", "Height of Particle" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "Particle_Mass", NC_REALTYPE, 2, dims, &m_id));
   set_nc_attributes(ncid, m_id, "grams", "Mass of Particle" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "Particle_Diameter", NC_REALTYPE, 2, dims, &d_id));
   set_nc_attributes(ncid, d_id, "meters", "Diameter of Particle" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "Particle_Density", NC_REALTYPE, 2, dims, &dn_id));
   set_nc_attributes(ncid, dn_id, "g/m3", "Density of Particle" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "Particle_Velocity", NC_REALTYPE, 2, dims, &v_id));
   set_nc_attributes(ncid, v_id, "m/s", "Height of Particle" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "Particle_vvel", NC_REALTYPE, 2, dims, &vv_id));
   set_nc_attributes(ncid, vv_id, "m/s", "Settling Velocity of Particle" PARAM_FILLVALUE);

   define_mode_off(&ncid);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
#endif


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
AED_REAL settling_velocity()
{
    return 0.01;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
