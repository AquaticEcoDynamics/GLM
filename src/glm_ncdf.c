/******************************************************************************
 *                                                                            *
 * glm_ncdf.c                                                                 *
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

#include <netcdf.h>

#include "glm.h"
#include "glm_types.h"
#include "glm_globals.h"
#include "glm_ncdf.h"

int ncid=-1;

static int set_no = -1;

//# dimension sizes
int x_dim, y_dim, z_dim, zone_dim, time_dim, restart_dim;

//# dimension lengths
static int lon_len=1;
static int lat_len=1;
static int height_len;
static int restart_len = 17;

static size_t start[4],edges[4];

static size_t start_r[1],edges_r[1];

//# variable ids
static int lon_id,lat_id,z_id,V_id,TV_id,Taub_id,NS_id,time_id;
static int HICE_id,HSNOW_id,HWICE_id, AvgSurfTemp_id;
static int precip_id,evap_id,rho_id,rad_id,extc_id,i0_id,wnd_id;
static int temp_id, salt_id, umean_id, uorb_id, restart_id;

#ifdef _WIN32
    char *strndup(const char *s, size_t len);
#endif

/*============================================================================*/
static void check_nc_error(int err);

/******************************************************************************
 *                                                                            *
 ******************************************************************************/
int init_glm_ncdf(const char *fn, const char *title, AED_REAL lat,
                                AED_REAL lon, int nlev, const char *start_time)
{
    int ncid;
    char time_str[128], history[128];
    extern char wq_lib[];

    int dims[4];
/*----------------------------------------------------------------------------*/

    if ( fn == NULL || strcmp(fn, "")==0 ) return -1;

    check_nc_error(nc_create(fn, NC_CLOBBER, &ncid));

    snprintf(time_str, 128, "hours since %s", start_time);
    snprintf(history,  128, "Created by glm3/%s v. %s", wq_lib, GLM_VERSION);

    height_len = nlev;

    //# define dimensions
    check_nc_error(nc_def_dim(ncid, "lon", 1, &x_dim));
    check_nc_error(nc_def_dim(ncid, "lat", 1, &y_dim));
    check_nc_error(nc_def_dim(ncid, "z", nlev, &z_dim));
    check_nc_error(nc_def_dim(ncid, "restart", 17, &restart_dim));
    if ( n_zones > 0 )
        check_nc_error(nc_def_dim(ncid, "nzones", n_zones, &zone_dim));
    check_nc_error(nc_def_dim(ncid, "time", NC_UNLIMITED, &time_dim));

    //# define coordinates
    dims[0] = x_dim;
    check_nc_error(nc_def_var(ncid, "lon", NC_REALTYPE, 1, dims, &lon_id));

    dims[0] = y_dim;
    check_nc_error(nc_def_var(ncid, "lat", NC_REALTYPE, 1, dims, &lat_id));

    dims[0] = time_dim;
    check_nc_error(nc_def_var(ncid, "NS",    NC_INT,      1, dims, &NS_id));
    check_nc_error(nc_def_var(ncid, "time",  NC_REALTYPE, 1, dims, &time_id));
    check_nc_error(nc_def_var(ncid, "hice",  NC_REALTYPE, 1, dims, &HICE_id));
    check_nc_error(nc_def_var(ncid, "hsnow", NC_REALTYPE, 1, dims, &HSNOW_id));
    check_nc_error(nc_def_var(ncid, "hwice", NC_REALTYPE, 1, dims, &HWICE_id));
    check_nc_error(nc_def_var(ncid, "avg_surf_temp", NC_REALTYPE, 1, dims, &AvgSurfTemp_id));
    
    dims[0] = restart_dim;
    check_nc_error(nc_def_var(ncid, "restart_variables", NC_REALTYPE, 1, dims, &restart_id));

    /**************************************************************************
     * define variables                                                       *
     **************************************************************************/
     
    //# x,y,t
    dims[2] = x_dim;
    dims[1] = y_dim;
    dims[0] = time_dim;

    check_nc_error(nc_def_var(ncid, "precip", NC_REALTYPE, 3, dims, &precip_id));
    check_nc_error(nc_def_var(ncid, "evap",   NC_REALTYPE, 3, dims, &evap_id));
    check_nc_error(nc_def_var(ncid, "I_0",    NC_REALTYPE, 3, dims, &i0_id));
    check_nc_error(nc_def_var(ncid, "wind",   NC_REALTYPE, 3, dims, &wnd_id));
    check_nc_error(nc_def_var(ncid, "Tot_V",  NC_REALTYPE, 3, dims, &TV_id));

    //# x,y,z,t
    dims[3] = x_dim;
    dims[2] = y_dim;
    dims[1] = z_dim;
    dims[0] = time_dim;

    check_nc_error(nc_def_var(ncid, "z",         NC_REALTYPE, 4, dims, &z_id));
    check_nc_error(nc_def_var(ncid, "V",         NC_REALTYPE, 4, dims, &V_id));
    check_nc_error(nc_def_var(ncid, "salt",      NC_REALTYPE, 4, dims, &salt_id));
    check_nc_error(nc_def_var(ncid, "temp",      NC_REALTYPE, 4, dims, &temp_id));

    check_nc_error(nc_def_var(ncid, "rho",       NC_REALTYPE, 4, dims, &rho_id));
    check_nc_error(nc_def_var(ncid, "rad",       NC_REALTYPE, 4, dims, &rad_id));
    check_nc_error(nc_def_var(ncid, "extc_coef", NC_REALTYPE, 4, dims, &extc_id));

    check_nc_error(nc_def_var(ncid, "u_mean",    NC_REALTYPE, 4, dims, &umean_id));
    check_nc_error(nc_def_var(ncid, "u_orb",     NC_REALTYPE, 4, dims, &uorb_id));
    check_nc_error(nc_def_var(ncid, "Taub",      NC_REALTYPE, 4, dims, &Taub_id));

    /**************************************************************************
     * assign attributes                                                      *
     **************************************************************************/

    //# coordinates
    set_nc_attributes(ncid, lon_id,    "degrees_east",  NULL        PARAM_FILLVALUE);
    set_nc_attributes(ncid, lat_id,    "degrees_north", NULL        PARAM_FILLVALUE);

    set_nc_attributes(ncid, time_id,   time_str,        NULL        PARAM_FILLVALUE);

    nc_put_att(ncid, NS_id, "long_name", NC_CHAR, 16, "Number of Layers");

    set_nc_attributes(ncid, HICE_id,   "meters",  "Height of Ice"   PARAM_FILLVALUE);
    set_nc_attributes(ncid, HSNOW_id,  "meters",  "Height of Snow"  PARAM_FILLVALUE);
    set_nc_attributes(ncid, HWICE_id,  "meters",  "Height of WhiteIce" PARAM_FILLVALUE);
    set_nc_attributes(ncid, AvgSurfTemp_id,  "celsius",  "Running average surface temperature" PARAM_FILLVALUE);
     
    set_nc_attributes(ncid, restart_id,  "various",  "dep_mx,prev_thick,g_prime_two_layer,\
energy_avail_max,mass_epi,old_slope,time_end_shear,time_start_shear,\
time_count_end,time_count_sim,half_seiche_period,thermocline_height,\
f0, fsum,u_f,u0,u_avg" PARAM_FILLVALUE);

    //# x,y,t
    set_nc_attributes(ncid, precip_id, "m/s",     "precipitation"   PARAM_FILLVALUE);
    set_nc_attributes(ncid, evap_id,   "m/s",     "evaporation"     PARAM_FILLVALUE);
    set_nc_attributes(ncid, i0_id,     "10E-6m",  "Shortwave"       PARAM_FILLVALUE);
    set_nc_attributes(ncid, wnd_id,    "m/s",     "wind"            PARAM_FILLVALUE);
    set_nc_attributes(ncid, TV_id,     "m3",      "lake volume"     PARAM_FILLVALUE);

    //# x,y,z,t
    set_nc_attributes(ncid, z_id,      "meters",  "layer heights"   PARAM_FILLVALUE);
    set_nc_attributes(ncid, V_id,      "m3",      "layer volume"    PARAM_FILLVALUE);
    set_nc_attributes(ncid, salt_id,   "g/kg",    "salinity"        PARAM_FILLVALUE);
    set_nc_attributes(ncid, temp_id,   "celsius", "temperature"     PARAM_FILLVALUE);

    set_nc_attributes(ncid, extc_id,   "unknown", "extc_coef"       PARAM_FILLVALUE);
    set_nc_attributes(ncid, rho_id,    "unknown", "density"         PARAM_FILLVALUE);
    set_nc_attributes(ncid, rad_id,    "unknown", "solar radiation" PARAM_FILLVALUE);

    set_nc_attributes(ncid, umean_id,  "m/s",     "mean velocity"    PARAM_FILLVALUE);
    set_nc_attributes(ncid, uorb_id,   "m/s",     "orbital velocity" PARAM_FILLVALUE);

    set_nc_attributes(ncid, Taub_id,   "N/m2",    "layer stress"    PARAM_FILLVALUE);

    //# global attributes
    nc_put_att(ncid, NC_GLOBAL, "Title", NC_CHAR, strlen(title), title);

    nc_put_att(ncid, NC_GLOBAL, "history", NC_CHAR, strlen(history), history);
    check_nc_error(nc_put_att(ncid, NC_GLOBAL, "Conventions", NC_CHAR, 6, "COARDS"));

    /**************************************************************************
     * Save some static data                                                  *
     **************************************************************************/

    //# leave define mode
    check_nc_error(nc_enddef(ncid));

    //# save latitude and logitude
    store_nc_scalar(ncid, lon_id, POINT, lon);
    store_nc_scalar(ncid, lat_id, POINT, lat);

    check_nc_error(nc_sync(ncid));

    return ncid;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void write_glm_ncdf(int ncid, int wlev, int nlev, int stepnum, AED_REAL timestep)
{
    AED_REAL temp_time, LakeVolume;
    AED_REAL *heights, *vols, *salts, *temps, *dens, *qsw, *extc_coef;
    AED_REAL *u_mean, *u_orb, *taub, *restart_variables;
    int i, littoralLayer = 0;

    if (ncid == -1) return;

    set_no++;

    store_nc_integer(ncid, NS_id, T_SHAPE, wlev);

    temp_time = (stepnum * timestep) / 3600.;
    LakeVolume = 0.0;
    for (i = 0; i < wlev; i++)
        LakeVolume += Lake[i].LayerVol;

    //# Time varying data : t
    /*------------------------------------------------------------------------*/
    store_nc_scalar(ncid,  time_id, T_SHAPE, temp_time);
    store_nc_scalar(ncid,  HICE_id, T_SHAPE, SurfData.delzBlueIce);
    store_nc_scalar(ncid, HWICE_id, T_SHAPE, SurfData.delzWhiteIce);
    store_nc_scalar(ncid, HSNOW_id, T_SHAPE, SurfData.delzSnow);
    store_nc_scalar(ncid, AvgSurfTemp_id, T_SHAPE, AvgSurfTemp);

    store_nc_scalar(ncid, precip_id, XYT_SHAPE, MetData.Rain);
    store_nc_scalar(ncid,   evap_id, XYT_SHAPE, SurfData.Evap);

    store_nc_scalar(ncid,  i0_id, XYT_SHAPE, MetData.ShortWave);
    store_nc_scalar(ncid, wnd_id, XYT_SHAPE, MetData.WindSpeed);

    store_nc_scalar(ncid,  TV_id, XYT_SHAPE, LakeVolume);

    //# Restart variables
    /*------------------------------------------------------------------------*/
    restart_variables   = malloc(17*sizeof(AED_REAL));
    restart_variables[0] = DepMX;
    restart_variables[1] = PrevThick;
    restart_variables[2] = gPrimeTwoLayer;
    restart_variables[3] = Energy_AvailableMix;
    restart_variables[4] = Mass_Epi;
    restart_variables[5] = OldSlope;
    restart_variables[6] = Time_end_shear;
    restart_variables[7] = Time_start_shear;
    restart_variables[8] = Time_count_end_shear;
    restart_variables[9] = Time_count_sim;
    restart_variables[10] = Half_Seiche_Period;
    restart_variables[11] = Thermocline_Height;
    restart_variables[12] = FO;
    restart_variables[13] = FSUM;
    restart_variables[14] = u_f;
    restart_variables[15] = u0;
    restart_variables[16] = u_avg;

    start_r[0] = 0; edges_r[0] = restart_len;

    //store_nc_scalar(ncid,  restart_id, R_SHAPE, restart_variables);
    check_nc_error(nc_put_vara(ncid,  restart_id, start_r, edges_r, restart_variables));

    //# Time varying profile data : z,t
    /*------------------------------------------------------------------------*/
    start[3] = 0;      edges[3] = lon_len;
    start[2] = 0;      edges[2] = lat_len;
    start[1] = 0;      edges[1] = height_len;
    start[0] = set_no; edges[0] = 1;

    heights   = malloc(nlev*sizeof(AED_REAL));
    vols      = malloc(nlev*sizeof(AED_REAL));
    salts     = malloc(nlev*sizeof(AED_REAL));
    temps     = malloc(nlev*sizeof(AED_REAL));
    dens      = malloc(nlev*sizeof(AED_REAL));
    qsw       = malloc(nlev*sizeof(AED_REAL));
    extc_coef = malloc(nlev*sizeof(AED_REAL));
    u_mean    = malloc(nlev*sizeof(AED_REAL));
    u_orb     = malloc(nlev*sizeof(AED_REAL));
    taub      = malloc(nlev*sizeof(AED_REAL));

    for (i = 0; i < wlev; i++) {
        heights[i] = Lake[i].Height;
        vols[i] = Lake[i].LayerVol;
        salts[i] = Lake[i].Salinity;
        temps[i] = Lake[i].Temp;
        dens[i] = Lake[i].Density;
        qsw[i] = Lake[i].Light;
        extc_coef[i] = Lake[i].ExtcCoefSW;
        u_mean[i] = Lake[i].Umean;
        u_orb[i] = Lake[i].Uorb;
        taub[i] = Lake[i].LayerStress;
    }
    for (i = wlev; i < nlev; i++) {
        heights[i] = NC_FILLER;
        vols[i] = NC_FILLER;
        salts[i] = NC_FILLER;
        temps[i] = NC_FILLER;
        dens[i] = NC_FILLER;
        qsw[i] = NC_FILLER;
        extc_coef[i] = NC_FILLER;
        u_mean[i] = NC_FILLER;
        u_orb[i] = NC_FILLER;
        taub[i] = NC_FILLER;
    }
    if ( littoral_sw ) {
        littoralLayer = wlev;
        heights[littoralLayer] = Lake[littoralLayer].Height;
        vols[littoralLayer] = Lake[littoralLayer].LayerVol;
        salts[littoralLayer] = Lake[littoralLayer].Salinity;
        temps[littoralLayer] = Lake[littoralLayer].Temp;
        dens[littoralLayer] = NC_FILLER;
        qsw[littoralLayer] = NC_FILLER;
        extc_coef[littoralLayer] = NC_FILLER;
        u_mean[littoralLayer] = NC_FILLER;
        u_orb[littoralLayer] = NC_FILLER;
        taub[littoralLayer] = NC_FILLER;
    }
    
    check_nc_error(nc_put_vara(ncid,     z_id, start, edges, heights));
    check_nc_error(nc_put_vara(ncid,     V_id, start, edges, vols));
    check_nc_error(nc_put_vara(ncid,  salt_id, start, edges, salts));
    check_nc_error(nc_put_vara(ncid,  temp_id, start, edges, temps));
    check_nc_error(nc_put_vara(ncid,   rho_id, start, edges, dens));
    check_nc_error(nc_put_vara(ncid,   rad_id, start, edges, qsw));
    check_nc_error(nc_put_vara(ncid,  extc_id, start, edges, extc_coef));
    check_nc_error(nc_put_vara(ncid, umean_id, start, edges, u_mean));
    check_nc_error(nc_put_vara(ncid,  uorb_id, start, edges, u_orb));
    check_nc_error(nc_put_vara(ncid,  Taub_id, start, edges, taub));

    free(heights); free(vols); free(salts);
    free(temps);  free(dens); free(qsw);
    free(extc_coef); free(u_mean); free(u_orb);
    free(taub); free(restart_variables);

    check_nc_error(nc_sync(ncid));
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void close_glm_ncdf(int ncid)
{
   if (ncid != -1)
      check_nc_error(nc_close(ncid));

   set_no = -1;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
int new_nc_variable(int ncid, const char *name, int data_type,
                                                      int ndim, const int *dims)
{
   int id;

    if (ncid == -1) return -1;
    check_nc_error(nc_def_var(ncid,name,data_type,ndim,dims,&id));
    return id;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void set_nc_attributes(int ncid, int id, const char *units,
                                      const char *long_name, AED_REAL FillValue)
{
    if (ncid == -1) return;
    nc_put_att(ncid, id, "units", NC_CHAR, strlen(units), units);
    if ( long_name != NULL ) {
        nc_put_att(ncid, id, "long_name", NC_CHAR, strlen(long_name), long_name);
        nc_put_att(ncid, id, "_FillValue", NC_REALTYPE, 1, &FillValue);
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void store_nc_integer(int ncid, int id, int var_shape, int iscalar)
{
    int iret;

    if (ncid == -1) return;
    /*------------------------------------------------------------------------*/
    if (var_shape == POINT)
        iret = nc_put_var(ncid, id, &iscalar);
    else if (var_shape == T_SHAPE) {
        start[0] = set_no; edges[0] = 1;
        iret = nc_put_vara(ncid, id, start, edges, &iscalar);
    } else {
        fprintf(stderr, "store_nc_integer : non valid shape %d\n", var_shape);
        exit(1);
    }
    check_nc_error(iret);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void store_nc_scalar(int ncid, int id, int var_shape, AED_REAL scalar)
{
    int iret;

    if (ncid == -1) return;
    /*------------------------------------------------------------------------*/
    if (var_shape == POINT)
        iret = nc_put_var(ncid, id, &scalar);
    else {
        if (var_shape == T_SHAPE) {
            start[0] = set_no; edges[0] = 1;
        } else if (var_shape == XYT_SHAPE) {
            start[2] = 0;      edges[2] = lon_len;
            start[1] = 0;      edges[1] = lat_len;
            start[0] = set_no; edges[0] = 1;
        } else {
            fprintf(stderr, "store_nc_scalar : non valid shape %d\n", var_shape);
            exit(1);
        }
        iret = nc_put_vara(ncid, id, start, edges, &scalar);
    }
    check_nc_error(iret);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void store_nc_array(int ncid, int id, int var_shape, int nvals,
                                                   int maxvals, AED_REAL *array)
{
    int iret, i;
    AED_REAL *tarr;

    if (ncid == -1) return;
    /*------------------------------------------------------------------------*/
    if (var_shape == Z_SHAPE) {
        start[0] = 0;      edges[0] = height_len;
    } else if (var_shape == XYZT_SHAPE || var_shape == XYNT_SHAPE) {
        start[3] = 0;      edges[3] = lon_len;
        start[2] = 0;      edges[2] = lat_len;
        if ( var_shape == XYZT_SHAPE ) {
            start[1] = 0;      edges[1] = height_len;
        } else {
            start[1] = 0;      edges[1] = n_zones;
        }
        start[0] = set_no; edges[0] = 1;
    } else {
        fprintf(stderr, "store_nc_array : non valid shape %d\n", var_shape);
        exit(1);
    }
    tarr = malloc(maxvals*sizeof(AED_REAL));
    for (i = 0; i < nvals; i++) tarr[i] = array[i];
    for (i = nvals; i < maxvals; i++) tarr[i] = NC_FILLER;
    iret = nc_put_vara(ncid, id, start, edges, tarr);
    free(tarr);

    check_nc_error(iret);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


#undef DEBUG
/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void check_nc_error(int err)
{
   if (err != NC_NOERR) {
      fprintf(stderr, "Error : %s (%d)\n", nc_strerror(err), err);
#if DEBUG
      CRASH("check_nc_error");
#endif
   }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 * for wq                                                                     *
 *                                                                            *
 ******************************************************************************/
void define_mode_on(int *ncid) { nc_redef(*ncid); }
void define_mode_off(int *ncid) { nc_enddef(*ncid); }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
int new_nc_variable_(int *ncid, const char *name, int *len, int *data_type,
                                                     int *ndim, const int *dims)
{
   int i, n = *ndim, ret;
   char *s = strndup(name,*len);
   int dim[6];

   for (i = 0; i < n; i++) dim[n-i-1] = dims[i];
   ret = new_nc_variable(*ncid, s, *data_type, n, dim);
   free(s);

   return ret;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void set_nc_attributes_(int *ncid, int *id, const char *units,
                                     const char *long_name, AED_REAL *FillValue)
{ set_nc_attributes(*ncid, *id, units, long_name, *FillValue); }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void store_nc_scalar_(int *ncid, int *id, int *var_shape, AED_REAL *array)
{ store_nc_scalar(*ncid, *id, *var_shape, *array); }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
void store_nc_array_(int *ncid, int *id, int *var_shape, int *nvals,
                                                  int *maxvals, AED_REAL *array)
{ store_nc_array(*ncid, *id, *var_shape, *nvals, *maxvals, array); }
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
