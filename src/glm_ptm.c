/******************************************************************************
 *                                                                            *
 * glm_ptm.c                                                                  *
 *                                                                            *
 * Contains Particle Tracking Model                                           *
 *                                                                            *
 * Developed by :                                                             *
 *                                                                            *
 * Copyright 2024-2026 : The University of Western Australia                  *
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
#include "glm_wqual.h"

            // Cell suspended in water column
#define WATER  0
            // Maybe forming a bed layer; set when hits bottom; turn off when resuspended
#define BED    1
            // Maybe forming a scum layer
#define SCUM   2
            // Maybe exiting the lake
#define EXIT   3

#define STAT   0
#define IDX2   1
#define IDX3   2
#define LAYR   3
#define FLAG   4
#define PTID   5

#define MASS   0
#define DIAM   1
#define DENS   2
#define VVEL   3
#define HGHT   4



#define _PTM_Stat(grp,part,var) PTM_Stat[_IDX_3d(1,max_particle_num,4,grp,part,var)]
#define _PTM_Vars(grp,part,var) PTM_Vars[_IDX_3d(1,max_particle_num,Num_WQ_Vars,grp,part,var)]

AED_REAL get_particle_density(AED_REAL particle_density);
AED_REAL get_particle_diameter(AED_REAL particle_diameter);
AED_REAL get_settling_velocity(AED_REAL settling_velocity);
AED_REAL random_walk(AED_REAL dt, AED_REAL Height, AED_REAL K_z, AED_REAL K_prime_z, AED_REAL vvel);

/*============================================================================*/

//CONSTANTS
int num_particle_grp=1;
AED_REAL init_depth_min=0.0;
AED_REAL init_depth_max=2.0;
AED_REAL ptm_time_step=1/60;
AED_REAL ptm_diffusivity=1e-6;

// VARIABLES
LOGICAL sed_deactivation = FALSE;
CINTEGER num_particle_groups = 1;

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

//  int bla;
    int pg;

/*----------------------------------------------------------------------------*/
//BEGIN

    //NEW allocate AED ptm data structures, and GLM pointers
    api_set_glm_ptm(&num_particle_groups,&max_particle_num);   //_WQ_SET_GLM_PTM

    //printf("_PTM_Stat(0,5000-1,0)  %d \n" ,_PTM_Stat(0,5000-1,0));
    //printf("_PTM_Stat(0,5000-1,1)  %d \n"  ,_PTM_Stat(0,5000-1,1));
    //printf("_PTM_Stat(0,5000-1,2)  %d \n"  ,_PTM_Stat(0,5000-1,2));
    //printf("_PTM_Stat(0,5000-1,3)  %d \n"  ,_PTM_Stat(0,5000-1,3));
    //printf("_PTM_Stat(0,5000,0)  %d \n"  ,_PTM_Stat(0,5000,0));
    //printf("_PTM_Stat(0,5000,1)  %d \n"  ,_PTM_Stat(0,5000,1));
    //printf("_PTM_Stat(0,5000,2)  %d \n"  ,_PTM_Stat(0,5000,2));
    //printf("_PTM_Stat(0,5000,3)  %d \n"  ,_PTM_Stat(0,5000,3));
    //printf("_PTM_Stat(0,5000,4)  %d \n" ,_PTM_Stat(0,5000,4));
    //printf("_PTM_Stat(0,0,0)  %d \n"  ,_PTM_Stat(0,0,0));
    //printf("_PTM_Stat(0,1,1)  %d \n"  ,_PTM_Stat(0,1,1));
    //printf("_PTM_Stat(0,2,2)  %d \n"  ,_PTM_Stat(0,2,2));
    //printf("_PTM_Stat(0,3,3)  %d \n"  ,_PTM_Stat(0,3,3));

    //printf("_PTM_Vars(0,0,0)  %f \n"  ,_PTM_Vars(0,0,0));
    //printf("_PTM_Vars(0,0,0)  %f \n"  ,_PTM_Vars(0,0,1));
    //printf("_PTM_Vars(0,0,0)  %f \n"  ,_PTM_Vars(0,0,2));
    //printf("_PTM_Vars(0,2,2)  %f \n"  ,_PTM_Vars(0,2,2));


    //NEW initialise (integer) status array for all particles to 0
    for (pg = 0; pg < num_particle_groups; pg++) {
      for (p = 0; p < max_particle_num; p++) {
         _PTM_Stat(pg,p,STAT) = 0;     // 1 idx_stat
         _PTM_Stat(pg,p,IDX2) = 0;     // 2 idx_2
         _PTM_Stat(pg,p,IDX3) = 0;     // 3 idx_3
         _PTM_Stat(pg,p,LAYR) = 0;     // 4 idx_layer
         _PTM_Stat(pg,p,FLAG) = 3;     // 5 flag
         _PTM_Vars(pg,p,MASS) = 0.0;   //
         _PTM_Vars(pg,p,DIAM) = 0.0;
         _PTM_Vars(pg,p,DENS)  = 0.0;
         _PTM_Vars(pg,p,VVEL)  = 0.0;
      }
    }

    // Set initial active particle height within the water column
    upper_height = Lake[surfLayer].Height - init_depth_min;
     if(upper_height > Lake[surfLayer].Height){
            fprintf(stderr, "     ERROR: Particles cannot be initialized above the surface of the lake");
            exit(1);
        }
    lower_height = Lake[surfLayer].Height - init_depth_max;
     if(lower_height < 0){
            fprintf(stderr, "     ERROR: Particles cannot be initialized at a depth greater than the maximum lake depth");
            exit(1);
        }

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
    int p, tt, ij1, ij2, sub_steps, pg;
    AED_REAL dt, K_z, K_above, K_prime_z;
    float rand_float, prob, prev_height, x1, x2, y1, y2, a1, a2;

/*----------------------------------------------------------------------------*/
//BEGIN
    pg = 0;

    // Loop through sub-timesteps, incrementing position
    sub_steps = 60;
    dt = 1;
    for (tt = 1; tt < sub_steps; tt++) {
        for (p = 0; p < max_particle_num; p++) {
          if (_PTM_Stat(pg,p,STAT)>0) {

            // printf("void do_ptm_update() %d %f \n"  , _PTM_Stat(pg,p,STAT),_PTM_Vars(pg,p,HGHT));

            // Capture current height of particle to calculate probability of settling below
            prev_height = _PTM_Vars(pg,p,HGHT);

            // Update particle position based on diffusivity and vert velocity
            _PTM_Stat(pg,p,FLAG)= WATER;
            K_z = Lake[_PTM_Stat(pg,p,LAYR)].Epsilon;
            K_above = Lake[_PTM_Stat(pg,p+1,LAYR)].Epsilon;

            // determine whether to assume molecular diffusion K
            if (K_z < ptm_diffusivity) K_z = ptm_diffusivity;

            if (K_above < ptm_diffusivity) K_above = ptm_diffusivity;

            if(_PTM_Stat(pg,p,LAYR) == surfLayer){
                K_prime_z = 0;
                continue;
            } else {
                K_prime_z = fabs(K_z - K_above);
            }

            _PTM_Vars(pg,p,HGHT) = random_walk(dt,_PTM_Vars(pg,p,HGHT), K_z, K_prime_z, _PTM_Vars(pg,p,VVEL));

            if (prev_height > _PTM_Vars(pg,p,HGHT)){

                // Get area at previous particle height
                x1 = prev_height * 10.0;
                y1 = x1 - (int)(x1 / 1.0) * 1.0;
                ij1 = (int)(x1 - y1) - 1;
                if(ij1 >= Nmorph){
                    y1 = y1 + (float)(ij1 - Nmorph);
                    ij1 = Nmorph - 1;
                } else if (ij1 < 0 ) ij1 = 0;
                a1 = MphLevelArea[ij1] + y1 * dMphLevelArea[ij1];

                // Get area at depth of current particle height
                x2 = _PTM_Vars(pg,p,HGHT) * 10.0;
                y2 = x2 - (int)(x2 / 1.0) * 1.0;
                ij2 = (int)(x2 - y2) - 1;
                if(ij2 >= Nmorph){
                    y2 = y2 + (float)(ij2 - Nmorph);
                    ij2 = Nmorph - 1;
                } else if (ij2 < 0 ) ij2 = 0;
                a2 = MphLevelArea[ij2] + y2 * dMphLevelArea[ij2];

                // Calculate proportional difference between two areas
                prob = (a1 - a2) / a1 * settling_efficiency;

                // Bernoulli draw to determine if particle should be assigned as BED
                rand_float = ((float)rand())/RAND_MAX;
                if(rand_float < prob || _PTM_Vars(pg,p,HGHT) < 0.02){
                    _PTM_Stat(pg,p,FLAG)= BED;
                    if(_PTM_Vars(pg,p,HGHT)<0.0){
                        _PTM_Vars(pg,p,HGHT)=0.0;
                    }
                    if(sed_deactivation){
                        _PTM_Stat(pg,p,STAT) = 0;
                    }
                }
            }

            // Determine if particle should be assigned as SCUM
            if(_PTM_Vars(pg,p,HGHT)>Lake[surfLayer].Height){
                _PTM_Stat(pg,p,FLAG)= SCUM;                      // Maybe forming a scum layer
                _PTM_Vars(pg,p,HGHT)=Lake[surfLayer].Height;
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

    int rand_int, pg;
    AED_REAL height_range;

/*----------------------------------------------------------------------------*/
//BEGIN
    // Get vertical range in the water column that mixed
    height_range = upper_height - lower_height;

    pg = 0;
    // Check for active particles in the height range
    for (p = 0; p < max_particle_num; p++) {
        if (_PTM_Stat(pg,p,STAT)>0) {
            if (_PTM_Vars(pg,p,HGHT)>=lower_height && _PTM_Vars(pg,p,HGHT)<=upper_height ) {
                // Particle is in the mixing zone, so re-position
                rand_int = rand() % 100 + 1;                            // random draw from unit distribution
                double random_double = (double)rand_int / 100;
                random_double = random_double * height_range;                 // scale unit random to requested range
                _PTM_Vars(pg,p,HGHT) = lower_height + random_double;
            }
        }
    }
    ptm_update_layerid();     // assign layers to active particles
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *    This routine adds new particles within the provided depth               *
 *    range                                                                   *                                                                        *
 *                                                                            *
 ******************************************************************************/
void ptm_addparticles(int new_particles, int max_particle_num, AED_REAL upper_height,
                      AED_REAL lower_height)
{
//LOCALS
    int p, n, pg, pid;

    int rand_int;
    AED_REAL height_range;

/*----------------------------------------------------------------------------*/
//BEGIN
    // Get vertical range in the water column that mixed
    height_range = upper_height - lower_height;
    n = 0;

    pg = 0;
    // For each new particle, initialise their properties and height
    for (p = 0 ; ; p++) {
        if(n == new_particles){
            break;
        }

        if(p >= max_particle_num){
            printf("ptm_addparticles(): WARNING no more available particles; skipping particle initialization");
            break;
        } else {
            if( _PTM_Stat(pg,p,STAT) == 0 && _PTM_Stat(pg,p,FLAG) == 3){ // find the first inactive particles with EXIT flag
                printf("ptm_addparticles() %d %d %d \n"  , p,n,new_particles);
                _PTM_Stat(pg,p,STAT) = 1;
                _PTM_Stat(pg,p,FLAG) = 0;
                _PTM_Vars(pg,p,MASS) = 1.0;
                _PTM_Vars(pg,p,DIAM) = get_particle_diameter(particle_diameter);
                _PTM_Vars(pg,p,DENS) = get_particle_density(particle_density);
                _PTM_Vars(pg,p,VVEL) = get_settling_velocity(settling_velocity);

                // Assign PTID
                if(_PTM_Stat(pg,p,PTID) < 0){
                   _PTM_Stat(pg,p,PTID) = p+1;
                } else {
                   pid = (int) floor(_PTM_Stat(pg,p,PTID) / max_particle_num);
                   _PTM_Stat(pg,p,PTID) = max_particle_num + pid*max_particle_num + (p+1 - pid*max_particle_num);
                }

                // Assign particles initial height
                rand_int = rand() % 100 + 1;                            // random draw from unit distribution
                double random_double = (double)rand_int / 100;
                random_double = random_double * height_range;           // scale unit random to requested range
                _PTM_Vars(pg,p,HGHT) = lower_height + random_double;     // set particle height

                // Adjust counter
                n++;
            }
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
    int p, pg;

/*----------------------------------------------------------------------------*/
//BEGIN
    pg = 0;
    // For each particle, draw from Bernoulli distribution to see whether removed from layer
    layer_prop = delta_vol / layer_vol;
    for (p = 0; p < max_particle_num; p++) {
        if(_PTM_Stat(pg,p,STAT) == 1 && _PTM_Stat(pg,p,LAYR) == layer_id){
            rand_float = ((float)rand())/RAND_MAX;
            if(rand_float <= layer_prop){
                // If particle leaves through outflow, reset completely
                _PTM_Stat(pg,p,STAT) = 0;
                _PTM_Stat(pg,p,FLAG) = 3;
                //_PTM_Vars(pg,p,MASS) = 0.0;
                //_PTM_Vars(pg,p,DIAM) = 0.0;
                //_PTM_Vars(pg,p,DENS) = 0.0;
                //_PTM_Vars(pg,p,VVEL) = 0.0;
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
    int p,pg;

    AED_REAL upper_height;
    AED_REAL lower_height;

/*----------------------------------------------------------------------------*/
//BEGIN

    // Get vertical range in the water column that is impacted by the shift
    lower_height = shift_height;
    upper_height = shift_height + shift_amount;

    pg = 0;
    // Check for active particles in the impacted height range
    for (p = 0; p < max_particle_num; p++) {
        if (_PTM_Stat(pg,p,STAT)>0) {
            if (_PTM_Vars(pg,p,HGHT)>lower_height && _PTM_Vars(pg,p,HGHT)<upper_height ) {
                // Particle is in the impacted zone, so re-position (lift or drop)
                _PTM_Vars(pg,p,HGHT) = _PTM_Vars(pg,p,HGHT) + shift_amount;   // adjust particle height
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
    int p,i,pg;

/*----------------------------------------------------------------------------*/
//BEGIN
    pg = 0;
    for (p = 0; p < max_particle_num; p++) {
        if (_PTM_Stat(pg,p,STAT)>0) {
            for (i = botmLayer; i < NumLayers; i++) {
                if (_PTM_Vars(pg,p,HGHT)<Lake[i].Height) {
                    _PTM_Stat(pg,p,LAYR) = i;
                    _PTM_Stat(pg,p,IDX3) = i;
                    _PTM_Stat(pg,p,IDX2) = 1;
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

    del_t = dt*60; // number of minutes (dt) times seconds per minute

    random_float = -1+2*((float)rand())/RAND_MAX;            // random draw from uniform distribution [-1,1]

    updated_height = Height + K_prime_z * Height * del_t + random_float *
    sqrt((2 * K_z * (Height + 0.5 * K_prime_z * Height * del_t) * del_t) / (1.0/3)); // random walk

    updated_height = updated_height + vvel * del_t;   // account for sinking/floating; 
                                                      // vvel is per second so needs to
                                                      // be multiplied by del_t which
                                                      // is number of seconds of random
                                                      // walk timestep

    return updated_height;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *        This routine returns the settling velocity for a particle in m/s    *
 *                                                                            *
 ******************************************************************************/
AED_REAL get_settling_velocity(AED_REAL settling_velocity)
{
    return settling_velocity / 86400;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *        This routine returns the density for a particle                     *
 *                                                                            *
 ******************************************************************************/
AED_REAL get_particle_density(AED_REAL particle_density)
{
    return particle_density;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *        This routine returns the diameter for a particle                    *
 *                                                                            *
 ******************************************************************************/
AED_REAL get_particle_diameter(AED_REAL particle_diameter)
{
    return particle_diameter;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


static int h_id, m_id, d_id, dn_id, vv_id, par_id, tem_id, no3_id, nh4_id, frp_id, c_id, n_id, pho_id, chl_id, num_id, cdiv_id, topt_id, lnalphachl_id, stat_id, flag_id, ptid_id;
static int set_no_p = -1;
static size_t start[2],edges[2];


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
//void ptm_write_glm(int *ncid, int *wlev, int *nlev, int *lvl, int *point_nlevs)
void ptm_write_glm(int ncid, int max_particle_num)
{
//LOCALS
    int p,pg;
    AED_REAL *p_height, *mass, *diam, *density, *vvel, *par, *tem, *no3, *nh4, *frp, *c, *n, *pho, *chl, *num, *cdiv, *topt, *lnalphachl;
    int *status, *flag, *ptid;

/*----------------------------------------------------------------------------*/
//BEGIN

    pg = 0;
    set_no_p++;

    start[1] = 0;             edges[1] = max_particle_num;
    start[0] = set_no_p;      edges[0] = 1;

    p_height  = malloc(max_particle_num*sizeof(AED_REAL));
    mass  = malloc(max_particle_num*sizeof(AED_REAL));
    diam  = malloc(max_particle_num*sizeof(AED_REAL));
    density  = malloc(max_particle_num*sizeof(AED_REAL));
    vvel  = malloc(max_particle_num*sizeof(AED_REAL));
    par  = malloc(max_particle_num*sizeof(AED_REAL));
    tem  = malloc(max_particle_num*sizeof(AED_REAL));
    no3  = malloc(max_particle_num*sizeof(AED_REAL));
    nh4  = malloc(max_particle_num*sizeof(AED_REAL));
    frp  = malloc(max_particle_num*sizeof(AED_REAL));
    c  = malloc(max_particle_num*sizeof(AED_REAL));
    n  = malloc(max_particle_num*sizeof(AED_REAL));
    pho  = malloc(max_particle_num*sizeof(AED_REAL));
    chl  = malloc(max_particle_num*sizeof(AED_REAL));
    num  = malloc(max_particle_num*sizeof(AED_REAL));
    cdiv  = malloc(max_particle_num*sizeof(AED_REAL));
    topt  = malloc(max_particle_num*sizeof(AED_REAL));
    lnalphachl  = malloc(max_particle_num*sizeof(AED_REAL));
    status  = malloc(max_particle_num*sizeof(int));
    flag  = malloc(max_particle_num*sizeof(int));
    ptid  = malloc(max_particle_num*sizeof(int));

    for (p = 0; p < max_particle_num; p++) {
        p_height[p]         = _PTM_Vars(pg,p,HGHT);    //Particle[p].Height;                REAL
        mass[p]             = _PTM_Vars(pg,p,MASS);    //Particle[p].Mass;                  REAL
        diam[p]             = _PTM_Vars(pg,p,DIAM);    // Particle[p].Diam;                 REAL
        density[p]          = _PTM_Vars(pg,p,DENS);    //Particle[p].Density;               REAL
        vvel[p]             = _PTM_Vars(pg,p,VVEL)*86400;  //Particle[p].vvel;              REAL
                                         //  VVEL+1 = HGHT
        par[p]              = _PTM_Vars(pg,p,VVEL+2);  //PAR experienced by particle;       REAL //ML why is it +2?
        tem[p]              = _PTM_Vars(pg,p,VVEL+3);  //temp experienced by particle;      REAL
        no3[p]              = _PTM_Vars(pg,p,VVEL+4);  //NO3 experienced by particle;       REAL
        nh4[p]              = _PTM_Vars(pg,p,VVEL+5);  //NH4 experienced by particle;       REAL
        frp[p]              = _PTM_Vars(pg,p,VVEL+6);  //FRP experienced by particle;       REAL
        c[p]                = _PTM_Vars(pg,p,VVEL+7);  //internal particle C;               REAL
        n[p]                = _PTM_Vars(pg,p,VVEL+8);  //internal particle N;               REAL
        pho[p]              = _PTM_Vars(pg,p,VVEL+9);  //internal particle P;               REAL
        chl[p]              = _PTM_Vars(pg,p,VVEL+10);  //internal particle chl;             REAL
        num[p]              = _PTM_Vars(pg,p,VVEL+11);  //number of cells in particle;       REAL
        cdiv[p]             = _PTM_Vars(pg,p,VVEL+12);  //internal C threshold for division; REAL
        topt[p]             = _PTM_Vars(pg,p,VVEL+13);  //particle temperature optimum;      REAL
        lnalphachl[p]       = _PTM_Vars(pg,p,VVEL+14); //ln alpha chl of particle;          REAL
        status[p]           = _PTM_Stat(pg,p,STAT);    //Particle[p].Status;                INT
        flag[p]             = _PTM_Stat(pg,p,FLAG);    //Particle[p].Flag;                  INT
        ptid[p]             = _PTM_Stat(pg,p,PTID);    //Particle[p].PTID;                  INT
    }

    nc_put_vara(ncid, h_id, start, edges, p_height);
    nc_put_vara(ncid, m_id, start, edges, mass);
    nc_put_vara(ncid, d_id, start, edges, diam);
    nc_put_vara(ncid, dn_id, start, edges, density);
    nc_put_vara(ncid, vv_id, start, edges, vvel);
    nc_put_vara(ncid, par_id, start, edges, par);
    nc_put_vara(ncid, tem_id, start, edges, tem);
    nc_put_vara(ncid, no3_id, start, edges, no3);
    nc_put_vara(ncid, nh4_id, start, edges, nh4);
    nc_put_vara(ncid, frp_id, start, edges, frp);
    nc_put_vara(ncid, c_id, start, edges, c);
    nc_put_vara(ncid, n_id, start, edges, n);
    nc_put_vara(ncid, pho_id, start, edges, pho);
    nc_put_vara(ncid, chl_id, start, edges, chl);
    nc_put_vara(ncid, num_id, start, edges, num);
    nc_put_vara(ncid, cdiv_id, start, edges, cdiv);
    nc_put_vara(ncid, topt_id, start, edges, topt);
    nc_put_vara(ncid, lnalphachl_id, start, edges, lnalphachl);
    nc_put_vara(ncid, stat_id, start, edges, status);
    nc_put_vara(ncid, flag_id, start, edges, flag);
    nc_put_vara(ncid, ptid_id, start, edges, ptid);

    free(p_height);
    free(mass);
    free(diam);
    free(density);
    free(vvel);
    free(par);
    free(tem);
    free(no3);
    free(nh4);
    free(frp);
    free(c);
    free(n);
    free(pho);
    free(chl);
    free(num);
    free(cdiv);
    free(topt);
    free(lnalphachl);
    free(status);
    free(flag);
    free(ptid);

    check_nc_error(nc_sync(ncid));
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
// will also need to handle groups in write step; append group name onto
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

   check_nc_error(nc_def_var(ncid, "particle_par", NC_REALTYPE, 2, dims, &par_id));
   set_nc_attributes(ncid, par_id, "ummol m2 sec", "particle layer PAR" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_tem", NC_REALTYPE, 2, dims, &tem_id));
   set_nc_attributes(ncid, tem_id, "degC", "particle layer temperature" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_no3", NC_REALTYPE, 2, dims, &no3_id));
   set_nc_attributes(ncid, no3_id, "mmolN", "particle layer NO3" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_nh4", NC_REALTYPE, 2, dims, &nh4_id));
   set_nc_attributes(ncid, nh4_id, "mmolN", "particle layer NH4" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_frp", NC_REALTYPE, 2, dims, &frp_id));
   set_nc_attributes(ncid, frp_id, "mmolP", "particle layer FRP" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_c", NC_REALTYPE, 2, dims, &c_id));
   set_nc_attributes(ncid, c_id, "pmolC/cell", "cell C concentration" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_n", NC_REALTYPE, 2, dims, &n_id));
   set_nc_attributes(ncid, n_id, "pmolN/cell", "cell N concentration" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_pho", NC_REALTYPE, 2, dims, &pho_id));
   set_nc_attributes(ncid, pho_id, "pmolP/cell", "cell P concentration" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_chl", NC_REALTYPE, 2, dims, &chl_id));
   set_nc_attributes(ncid, chl_id, "pgChl", "cell Chl concentration" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_num", NC_REALTYPE, 2, dims, &num_id));
   set_nc_attributes(ncid, num_id, "number", "number of cells/particle" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_cdiv", NC_REALTYPE, 2, dims, &cdiv_id));
   set_nc_attributes(ncid, cdiv_id, "pmol/cell", "cellular carbon content threshold for division" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_topt", NC_REALTYPE, 2, dims, &topt_id));
   set_nc_attributes(ncid, topt_id, "degC", "optimal temperature" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_lnalphachl", NC_REALTYPE, 2, dims, &lnalphachl_id));
   set_nc_attributes(ncid, lnalphachl_id, "(W m-2)-1(gChl molC)-1d-1", "slope of the P-I curve" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_status", NC_INT, 2, dims, &stat_id));
   nc_put_att(ncid, stat_id, "long_name", NC_CHAR, 18, "Status of Particle");

   check_nc_error(nc_def_var(ncid, "particle_flag", NC_INT, 2, dims, &flag_id));
   nc_put_att(ncid, flag_id, "long_name", NC_CHAR, 18, "Location Flag of Particle");

   check_nc_error(nc_def_var(ncid, "particle_ptid", NC_INT, 2, dims, &ptid_id));
   nc_put_att(ncid, ptid_id, "long_name", NC_CHAR, 18, "ID of Particle");

   define_mode_off(&ncid);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
