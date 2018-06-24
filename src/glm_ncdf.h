/******************************************************************************
 *                                                                            *
 * glm_ncdf.h                                                                 *
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
#ifndef _GLM_NCDF_H_
#define _GLM_NCDF_H_

#define NF90_FLOAT  5
#define NF90_DOUBLE 6

#ifdef _FORTRAN_SOURCE_
!-------------------------------------------------------------------------------

  INTERFACE

     CINTEGER FUNCTION init_glm_ncdf(fn,title,lat,lon,nlev,start_time) BIND(C, name="init_glm_ncdf_")
        USE ISO_C_BINDING
        CCHARACTER,INTENT(in) :: fn(*),title(*),start_time(*)
        AED_REAL,INTENT(in)   :: lat,lon
        CINTEGER,INTENT(in)   :: nlev
     END FUNCTION init_glm_ncdf

     SUBROUTINE glm_nc_definitions(x, y, z, tm) BIND(C, name="glm_nc_definitions_")
        USE ISO_C_BINDING
        CINTEGER,INTENT(out) :: x, y, z, tm
     END SUBROUTINE glm_nc_definitions

     SUBROUTINE define_mode_on(ncid) BIND(C, name="define_mode_on")
        USE ISO_C_BINDING
        CINTEGER,INTENT(in) :: ncid
     END SUBROUTINE define_mode_on

     SUBROUTINE define_mode_off(ncid) BIND(C, name="define_mode_off")
        USE ISO_C_BINDING
        CINTEGER,INTENT(in) :: ncid
     END SUBROUTINE define_mode_off

     CINTEGER FUNCTION new_nc_variable(ncid,name,len,data_type,ndim,dims) BIND(C, name="new_nc_variable_")
        USE ISO_C_BINDING
#       if defined( _WIN32 ) && USE_DL_LOADER
        !DEC$ ATTRIBUTES DLLIMPORT :: new_nc_variable_
#       endif
        CINTEGER,INTENT(in)   :: ncid
        CCHARACTER,INTENT(in) :: name(*)
        CINTEGER,INTENT(in)   :: len
        CINTEGER,INTENT(in)   :: data_type, ndim
        CINTEGER,INTENT(in)   :: dims(*)
     END FUNCTION new_nc_variable

     SUBROUTINE set_nc_attributes(ncid, id, units, long_name, FillValue) BIND(C, name="set_nc_attributes_")
        USE ISO_C_BINDING
        CINTEGER,INTENT(in)   :: ncid, id
        CCHARACTER,INTENT(in) :: units(*), long_name(*)
        AED_REAL,INTENT(in)   :: FillValue
     END SUBROUTINE set_nc_attributes

     SUBROUTINE store_nc_integer(ncid, id, var_shape, iscalar) BIND(C, name="store_nc_integer_")
        USE ISO_C_BINDING
        CINTEGER,INTENT(in)  :: ncid, id, var_shape, iscalar
     END SUBROUTINE store_nc_integer

     SUBROUTINE store_nc_scalar(ncid, id, var_shape, scalar) BIND(C, name="store_nc_scalar_")
        USE ISO_C_BINDING
#       if defined( _WIN32 ) && USE_DL_LOADER
        !DEC$ ATTRIBUTES DLLIMPORT :: store_nc_scalar_
#       endif
        CINTEGER,INTENT(in) :: ncid, id, var_shape
        AED_REAL,INTENT(in) :: scalar
     END SUBROUTINE store_nc_scalar

     SUBROUTINE store_nc_array(ncid, id, var_shape, nvals, maxvals, array) BIND(C, name="store_nc_array_")
        USE ISO_C_BINDING
#       if defined( _WIN32 ) && USE_DL_LOADER
        !DEC$ ATTRIBUTES DLLIMPORT :: store_nc_array_
#       endif
        CINTEGER,INTENT(in)    :: ncid, id, var_shape, nvals, maxvals
        AED_REAL,INTENT(inout) :: array(*)
     END SUBROUTINE store_nc_array

  END INTERFACE

#else

  int init_glm_ncdf(const char *fn, const char *title, AED_REAL lat,
                                  AED_REAL lon, int nlev, const char *start_time);
  void write_glm_ncdf(int ncid, int wlev, int nlev, int stepnum, AED_REAL timestep);
  void close_glm_ncdf(int ncid);
  void define_mode_on(int *ncid);
  void define_mode_off(int *ncid);
  int new_nc_variable(int ncid, const char *name, int data_type, int ndim, const int *dims);
  void set_nc_attributes(int ncid, int id, const char *units,
                                      const char *long_name, AED_REAL FillValue);
  void store_nc_integer(int ncid, int id, int var_shape, int iscalar);
  void store_nc_scalar(int ncid, int id, int var_shape, AED_REAL scalar);
  void store_nc_array(int ncid, int id, int var_shape, int nvals, int maxvals, AED_REAL *array);

  extern int ncid, x_dim, y_dim, z_dim, zone_dim, time_dim;

#endif

#endif
