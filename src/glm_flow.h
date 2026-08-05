/******************************************************************************
 *                                                                            *
 * glm_flow.h                                                                 *
 *                                                                            *
 * Developed by :                                                             *
 *     AquaticEcoDynamics (AED) Group                                         *
 *     School of Agriculture and Environment                                  *
 *     The University of Western Australia                                    *
 *                                                                            *
 *     http://aquatic.science.uwa.edu.au/                                     *
 *                                                                            *
 * Copyright 2013-2026 : The University of Western Australia                  *
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
#ifndef _GLM_FLOW_H_
#define _GLM_FLOW_H_

#include "glm.h"

AED_REAL do_outflows(int jday, AED_REAL day_fraction);
AED_REAL do_overflow(int jday, AED_REAL day_fraction);
AED_REAL do_inflows(void);

//# Withdraw a volume at a given height (selective-withdrawal physics), used by
//# the oxygenation recirculation. Samples withdrawal-layer T/S/WQ into the
//# output args and returns the volume actually removed.
AED_REAL oxy_recirc_extract(AED_REAL height, AED_REAL want_vol,
                            AED_REAL *temp, AED_REAL *salt, AED_REAL *wq);

#endif
