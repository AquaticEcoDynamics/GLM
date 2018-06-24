/******************************************************************************
 *                                                                            *
 * glm_csv.h                                                                  *
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
#ifndef _GLM_CSV_H_
#define _GLM_CSV_H_

#define MaxPointCSV    10
#define MaxCSVOutVars   20

#ifdef _FORTRAN_SOURCE_

  INTERFACE

     SUBROUTINE write_csv_point(f, name, len, val, cval, vlen, last) BIND(C, name="write_csv_point_")
        USE ISO_C_BINDING
#       if defined( _WIN32 ) && USE_DL_LOADER
        !DEC$ ATTRIBUTES DLLIMPORT :: write_csv_point_
#       endif
        CINTEGER,INTENT(in)   :: f
        CCHARACTER,INTENT(in) :: name(*)
        CINTEGER,INTENT(in)   :: len
        AED_REAL,INTENT(in)   :: val
        CCHARACTER,INTENT(in) :: cval(*)
        CINTEGER,INTENT(in)   :: vlen
        CLOGICAL,INTENT(in)   :: last
     END SUBROUTINE write_csv_point

     SUBROUTINE write_csv_lake(name,len,val,cval,vlen,last) BIND(C,name="write_csv_lake_")
        USE ISO_C_BINDING
        CCHARACTER,INTENT(in) :: name(*)
        CINTEGER,INTENT(in)   :: len
        AED_REAL,INTENT(in)   :: val
        CCHARACTER,INTENT(in) :: cval(*)
        CINTEGER,INTENT(in)   :: vlen
        CLOGICAL,INTENT(in)   :: last
     END SUBROUTINE write_csv_lake

     SUBROUTINE close_csv_point_output() BIND(C, name="close_csv_point_output")
     END SUBROUTINE close_csv_point_output

     SUBROUTINE close_csv_lake_output() BIND(C, name="close_csv_lake_output")
     END SUBROUTINE close_csv_lake_output

  END INTERFACE
#else

/*############################################################################*/

extern AED_REAL csv_point_at[];
extern int csv_point_frombot[];
extern int csv_point_nlevs;
extern int csv_lake_file;

void init_csv_output(const char *out_dir);
void write_csv_point(int f, const char *name, AED_REAL val, const char *cval, int last);
void write_csv_point_(int *f, const char *name, int *len, AED_REAL *val, const char *cval, int *vlen, int *last);
void write_csv_lake(const char *name, AED_REAL val, const char *cval, int last);
void write_csv_lake_(const char *name, int *len, AED_REAL *val, const char *cval, int *vlen, int *last);
void glm_close_csv_output(void);
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/

void configure_csv(int point_nlevs, AED_REAL *point_at, const char *point_fname,
                             int *point_frombot, int point_nvars, const char *lake_fname);

void set_csv_point_varname(int which, const char *point_varnam);

void configure_outfl_csv(int outlet_allinone, const char *outfl_fname,
                                        int outfl_nvars, const char *ovrfl_fname);
void set_csv_outfl_varname(int which, const char *outfl_varnam);
void write_csv_outfl(int f, const char *name, AED_REAL val, const char *cval, int last);
void write_csv_outfl_idx(int ofl, int var, AED_REAL val, const char *cval, int last);

#endif

#endif
