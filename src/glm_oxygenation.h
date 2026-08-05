/******************************************************************************
 *                                                                            *
 * glm_oxygenation.h                                                          *
 *                                                                            *
 * Artificial oxygenation / aeration system for GLM.                          *
 *                                                                            *
 * Three approaches:                                                          *
 *   1. Direct addition  - add O2 at a user-defined rate to the layer at a    *
 *                         specified height. No water is added; only the O2   *
 *                         mass in that layer increases.                      *
 *   2. CSV-driven       - as (1) but the daily O2 input (and height) is read *
 *                         from a time-series CSV, like an inflow file.       *
 *   3. Withdraw+reinject- pull water from a depth (carrying its T, S and all *
 *                         WQ), add O2, and return it at another depth. Reuses*
 *                         the heat-pump capture/inject pattern.              *
 *                                                                            *
 *     AquaticEcoDynamics (AED) Group                                         *
 *     The University of Western Australia                                    *
 *                                                                            *
 *  Copyright 2013-2026 : The University of Western Australia                 *
 *                                                                            *
 *  This file is part of GLM (General Lake Model)                             *
 *                                                                            *
 ******************************************************************************/
#ifndef _GLM_OXYGENATION_H_
#define _GLM_OXYGENATION_H_

#include "glm_types.h"

//# Resolve O2 variable index and validate config (call from glm_init.c)
void init_oxygenation(void);

//# Open the CSV time-series for a direct-addition device (Approach 2)
void oxy_open_file(int idx, const char *fname, const char *timefmt);

//# Open the CSV time-series for the recirculation O2 loading (Approach 3)
void oxy_open_recirc_file(const char *fname, const char *timefmt);

//# Print a summary of the configured oxygenation system (call from glm_init.c)
void check_oxygenation_config(void);

//# Read today's CSV values for CSV-driven devices and the recirculation stream
void read_daily_oxygenation(int julian);

//# Approaches 1 & 2: add O2 to the target layer(s). Returns total O2 mass added.
AED_REAL do_oxygenation(AED_REAL day_fraction);

//# Approach 3 (self-contained): withdraw water at oxy_recirc_withdraw_height,
//# add O2, and return it at oxy_recirc_return_height. Call from glm_model.c
//# after do_outflows. Returns the O2 mass added.
AED_REAL oxy_do_recirculation(AED_REAL day_fraction);

#endif
