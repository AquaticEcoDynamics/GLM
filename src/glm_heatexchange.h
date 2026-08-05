/******************************************************************************
 *                                                                            *
 * glm_heatexchange.h                                                         *
 *                                                                            *
 * Developed by:                                                              *
 * Taynara Fernandes                                                          *
 * Matt Hipsey                                                                *
 * Casper Boon                                                                *
 *                                                                            *
 * Helmholtz Centre for Environmental Research (UFZ)                          *
 * Department of Lake Research (SEEFO)                                        *
 *                                                                            *
 ******************************************************************************/
#ifndef _GLM_HEATEXCHANGE_H_
#define _GLM_HEATEXCHANGE_H_

#include "glm_types.h"

void heat_pump_capture_outflow(int jday, AED_REAL DrawHeight, AED_REAL vol, AED_REAL temp, AED_REAL salt, AED_REAL *wq_vars);
void heat_pump_insert_inflow(void);
void init_heat_pump(void);
void check_heat_pump_config(void);

#endif
