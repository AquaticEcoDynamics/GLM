/******************************************************************************
 *                                                                            *
 * glm_mobl.h                                                                 *
 *                                                                            *
 * code for mobility                                                          *
 *                                                                            *
 * Developed by :                                                             *
 *     AquaticEcoDynamics (AED) Group                                         *
 *     School of Agriculture and Environment                                  *
 *     The University of Western Australia                                    *
 *                                                                            *
 *     http://aquatic.science.uwa.edu.au/                                     *
 *                                                                            *
 * Copyright 2013 - 2024 -  The University of Western Australia               *
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
#ifndef _GLM_MOBL_H_
#define _GLM_MOBL_H_

#include "glm.h"

#ifdef _FORTRAN_SOURCE_
!###############################################################################

  INTERFACE

    SUBROUTINE doMobility(N,dt,h,A,ww,min_C,cc) BIND(C, name="doMobility")
       USE ISO_C_BINDING
#      if defined( _WIN32 ) && USE_DL_LOADER
       !DEC$ ATTRIBUTES DLLIMPORT :: doMobility
#      endif
       CINTEGER,INTENT(in)     :: N       !# number of vertical layers
       AED_REAL,INTENT(in)     :: dt      !# time step (s)
       AED_REAL,INTENT(in)     :: h(*)    !# layer thickness (m)
       AED_REAL,INTENT(in)     :: A(*)    !# layer areas (m2)
       AED_REAL,INTENT(in)     :: ww(*)   !# vertical speed (m/s)
       AED_REAL,INTENT(in)     :: min_C   !# minimum allowed cell concentration
       AED_REAL,INTENT(inout)  :: cc(*)   !# cell concentration
    END SUBROUTINE doMobility

  END INTERFACE
!+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#else

   void doMobility(const int *N, const AED_REAL *dt,
                         const AED_REAL *h,  const AED_REAL *A,
                         const AED_REAL *ww, const AED_REAL *min_C, AED_REAL *cc);

#endif

#endif
