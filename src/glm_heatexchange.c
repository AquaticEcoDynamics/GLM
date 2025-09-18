/******************************************************************************
*                                                                             *
* glm_heatexchange.c                                                          *
*                                                                             *
* Developed by:                                                               *
* Taynara Fernandes                                                           *
* Matt Hipsey                                                                 *
* Casper Boon                                                                 *
*                                                                             *
* Helmholtz Centre for Environmental Research (UFZ)                           *
* Department of Lake Research (SEEFO)                                         *
*                                                                             *
******************************************************************************/

#include <stdio.h>
#include <math.h>
#include "glm_types.h" 
#include "glm_const.h"
#include "glm_globals.h"
#include "glm_util.h"
#include "glm_input.h"

// Storage for intercepted outflow
static AED_REAL stored_flow_rate  = 0.0;
static AED_REAL stored_temp       = 0.0;
static AED_REAL stored_salt       = 0.0;
static AED_REAL stored_Drawheight = 0.0;
static int stored_jday            = -1;
static int print_counter          = 0;

static AED_REAL stored_WQ[MaxVars]; // WQ variables from outflow

/*****************************************************************************
 * Capture and store live flow data from GLM's output (called in glm_flow.c) *
 *****************************************************************************/
void heat_pump_capture_outflow(int jday, AED_REAL DrawHeight, AED_REAL vol, AED_REAL temp, AED_REAL salt, AED_REAL *wq_vars)
{
    // Only capture if heat pump is enabled
    if (heat_pump_switch <= 0) return;

    // Clear previous values when capturing new flow
    if (stored_jday != -1 && stored_jday != jday) {
        stored_flow_rate  = 0.0;
        stored_temp       = 0.0;
        stored_salt       = 0.0;
        stored_Drawheight = 0.0;
        // Clear WQ variables
        if (Num_WQ_Vars > 0 && wq_vars != NULL) {
            for (int wqidx = 0; wqidx < Num_WQ_Vars; wqidx++) {
                stored_WQ[wqidx] = 0.0;
            }
        }
    }

    // Store the captured flow data to modify heat pump backflow
    stored_flow_rate  = vol;
    stored_temp       = temp;
    stored_salt       = salt;
    stored_Drawheight = DrawHeight;
    stored_jday       = jday;
    
    // Store WQ variables
    if (Num_WQ_Vars > 0 && wq_vars != NULL) {
        for (int wqidx = 0; wqidx < Num_WQ_Vars; wqidx++) {
            stored_WQ[wqidx] = wq_vars[wqidx];
        }
    }
    print_counter++;  // Track capture events for printing
}
/*************************************************************************************
 * Insert previously stored flow (heat pump backflow) (called in glm_model.c l. 350) *
 *************************************************************************************/
void heat_pump_insert_inflow() 
{
    // Only proceed if heat pump is enabled
    if (heat_pump_switch <= 0) return;

    // Only proceed if we have captured flow data
    if (stored_jday == -1 || stored_flow_rate <= 0.0) return;
    
    // Check if the specified inflow index is valid
    if (heat_pump_inflow_idx < 0 || heat_pump_inflow_idx >= NumInf) {
        printf("ERROR: heat_pump_inflow_idx (%d) is out of range [0, %d]\n", 
               heat_pump_inflow_idx, NumInf-1);

        return;
    }

    // Calculate temperature change caused by the heat pump
    AED_REAL heated_temp;
    
    switch (heat_pump_switch) {
        case 1: {
            // Mode 1: Fixed temperature increase (defined in .nml file)
            heated_temp = stored_temp - heat_pump_temp_change;
            break;
        }
        case 2: {
            // Mode 2:  heat flux-based dT calculation
            AED_REAL flow_rate_m3s = stored_flow_rate / SecsPerDay; // m³/day to m³/s
            // Use dynamic heat flux if available, otherwise use static value
            AED_REAL current_heat_flux = (heat_pump_dynamic_heat_flux != 0.0) ? 
                                        heat_pump_dynamic_heat_flux : heat_pump_heat_flux;
            // ΔT = Q_heat / (ρ × Q_flow × c) Units: W / (kg/m³ × m³/s × J/(kg·K)) = J/s / (kg/s × J/(kg·K)) = K
            AED_REAL temp_change = current_heat_flux / (rho0 * flow_rate_m3s * SPHEAT);
            heated_temp = stored_temp - temp_change;
            break;
        }
        default: {
            // Default to mode 1 behavior for backward compatibility
            heated_temp = stored_temp - heat_pump_temp_change;
            break;
        }
    }
   
    // Modify the specified inflow to inject heated water at specified depth
    if (NumInf > 0) {
        
        // Assign flow, temp, salt, and submerged elevation
        Inflows[heat_pump_inflow_idx].FlowRate = stored_flow_rate;
        Inflows[heat_pump_inflow_idx].TemInf   = heated_temp;
        Inflows[heat_pump_inflow_idx].SalInf   = stored_salt + 0.5;                
        Inflows[heat_pump_inflow_idx].Factor   = 1.0;

        // Assign WQ variables
        if (Num_WQ_Vars > 0) {
            for (int wqidx = 0; wqidx < Num_WQ_Vars; wqidx++) {

                Inflows[heat_pump_inflow_idx].WQInf[wqidx] = stored_WQ[wqidx];
            }
        }

        // Minimal status output every 200 days to show dynamic heat flux changes
        if (print_counter % 200 == 0) {
            // First show extraction info
            if (heat_pump_outflow_idx >= 0 && heat_pump_outflow_idx < NumOut) {
                const char* outflow_type = (Outflows[heat_pump_outflow_idx].SubmElevDynamic) ? "dynamic" : "static";
                printf("Heat pump extracting at jday %d: Q=%8.4f m³/d, elev=%5.1f m[%s], T=%5.1f°C\n", 
                       stored_jday, stored_flow_rate, stored_Drawheight, outflow_type, stored_temp);
            } else {
                printf("Heat pump extracting at jday %d: Q=%8.4f m³/d, elev=%5.1f m[UNKNOWN], T=%5.1f°C\n", 
                       stored_jday, stored_flow_rate, stored_Drawheight, stored_temp);
            }
            
            // Then show injection info
            if (heat_pump_inflow_idx >= 0 && heat_pump_inflow_idx < NumInf) {
                const char* inflow_type = (Inflows[heat_pump_inflow_idx].SubmElevDynamic) ? "dynamic" : "static";
                // Show different output based on heat pump mode
                if (heat_pump_switch == 1) {
                    printf("Heat pump injecting  at jday %d: Q=%8.4f m³/d, elev=%5.1f m[%s], T=%5.1f°C (ΔT=%+6.2f°C)\n",
                           stored_jday, stored_flow_rate, Inflows[heat_pump_inflow_idx].SubmElev, inflow_type, heated_temp, heat_pump_temp_change);
                } else if (heat_pump_switch == 2) {
                    AED_REAL temp_change = heated_temp - stored_temp;
                    AED_REAL current_heat_flux = (heat_pump_dynamic_heat_flux != 0.0) ? 
                                                heat_pump_dynamic_heat_flux : heat_pump_heat_flux;
                    printf("Heat pump injecting  at jday %d: Q=%8.4f m³/d, elev=%5.1f m[%s], T=%5.1f°C (ΔT=%+6.3f°C, %5.0fW)\n",
                           stored_jday, stored_flow_rate, Inflows[heat_pump_inflow_idx].SubmElev, inflow_type, heated_temp, temp_change, current_heat_flux);
                } else {
                    printf("Heat pump (unknown) injecting at jday %d: Q=%8.4f m³/d, elev=%5.1f m[%s], T=%5.1f°C\n",
                           stored_jday, stored_flow_rate, Inflows[heat_pump_inflow_idx].SubmElev, inflow_type, heated_temp);
                }
            } else {
                if (heat_pump_switch == 1) {
                    printf("Heat pump injecting  at jday %d: Q=%8.4f m³/d, T=%5.1f°C (ΔT=%+6.2f°C)                [inflow=INVALID]\n",
                           stored_jday, stored_flow_rate, heated_temp, heat_pump_temp_change);
                } else if (heat_pump_switch == 2) {
                    AED_REAL temp_change = heated_temp - stored_temp;
                    AED_REAL current_heat_flux = (heat_pump_dynamic_heat_flux != 0.0) ? 
                                                heat_pump_dynamic_heat_flux : heat_pump_heat_flux;
                    printf("Heat pump injecting  at jday %d: Q=%8.4f m³/d, T=%5.1f°C (ΔT=%+6.3f°C, %5.0fW)        [inflow=INVALID]\n",
                           stored_jday, stored_flow_rate, heated_temp, temp_change, current_heat_flux);
                } else {
                    printf("Heat pump (unknown) injecting at jday %d: Q=%8.4f m³/d, T=%5.1f°C                     [inflow=INVALID]\n",
                           stored_jday, stored_flow_rate, heated_temp);
                }
            }
            
            // Show all WQ variables if available
            if (Num_WQ_Vars > 0) {
                printf("  Injected WQ: ");
                for (int i = 0; i < Num_WQ_Vars; i++) {
                    printf("WQ[%d]=%.6f ", i, Inflows[heat_pump_inflow_idx].WQInf[i]);
                }
            }
        }
    }

}

/***********************************************************
 * Initialize heat pump system (called from glm_init.c)   *
 ***********************************************************/
void init_heat_pump() 
{
    // Heat pump initialization
}

/****************************************************************
 * Check heat pump configuration (called from glm_init.c)      *
 ****************************************************************/
void check_heat_pump_config() 
{
    if (heat_pump_switch > 0) {
        // Display mode-specific configuration
        if (heat_pump_switch == 1) {
            printf("Heat pump enabled: ΔT = %.2f°C, inflow %d ↔ outflow %d\n", 
                   heat_pump_temp_change, heat_pump_inflow_idx, heat_pump_outflow_idx);
            
            // Validate mode 1 requirements
            if (heat_pump_temp_change == 0.0) {
                printf("WARNING: heat_pump_temp_change is zero - no heating will occur in Mode 1\n");
            }
        } else if (heat_pump_switch == 2) {
            printf("Heat pump enabled: Q = %.0f W, inflow %d ↔ outflow %d\n", 
                   heat_pump_heat_flux, heat_pump_inflow_idx, heat_pump_outflow_idx);
            
            // Validate mode 2 requirements  
            if (heat_pump_heat_flux == 0.0) {
                printf("WARNING: heat_pump_heat_flux is zero - no heating will occur in Mode 2\n");
            }
            if (heat_pump_heat_flux < 0.0) {
                printf("WARNING: heat_pump_heat_flux is negative (%.0f W) - this will cool the water\n", heat_pump_heat_flux);
            }
        } else {
            printf("Heat pump enabled (Unknown Mode %d): inflow %d ↔ outflow %d\n", 
                   heat_pump_switch, heat_pump_inflow_idx, heat_pump_outflow_idx);
            printf("WARNING: Unsupported heat pump mode - using Mode 1 (fixed ΔT) as fallback\n");
        }
               
        // Validate indices 
        if (heat_pump_inflow_idx < 0 || heat_pump_inflow_idx >= NumInf) {
            printf("ERROR: heat_pump_inflow_idx (%d) is out of range [0, %d]\n", 
                   heat_pump_inflow_idx, NumInf-1);
        }
        if (heat_pump_outflow_idx < 0 || heat_pump_outflow_idx >= NumOut) {
            printf("ERROR: heat_pump_outflow_idx (%d) is out of range [0, %d]\n", 
                   heat_pump_outflow_idx, NumOut-1);
        }
        
        // Ensure heat pump inflow is configured as submerged with dynamic elevation 
        if (heat_pump_inflow_idx >= 0 && heat_pump_inflow_idx < NumInf) {
            Inflows[heat_pump_inflow_idx].SubmFlag        = TRUE;
            Inflows[heat_pump_inflow_idx].SubmElevDynamic = TRUE;
            
        }
    } else {
        printf("Heat pump disabled (heat_pump_switch = %d)\n", heat_pump_switch);
    }
}