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
#include "glm.h"
#include "glm_types.h"
#include "glm_const.h"
#include "glm_globals.h"
#include "glm_util.h"
#include "glm_input.h"
#include "glm_mixu.h"

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
 * Insert heat pump water directly                *
 * Called AFTER do_outflows() in glm_model.c for mass conservation                   *
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
    AED_REAL temp_change_value;

    switch (heat_pump_switch) {
        case 1: {
            // Mode 1: Fixed temperature increase (defined in .nml file)
            temp_change_value = heat_pump_temp_change;
            heated_temp = stored_temp - temp_change_value;
            break;
        }
        case 2: {
            // Mode 2: heat flux-based dT calculation
            AED_REAL flow_rate_m3s = stored_flow_rate / SecsPerDay; // m³/day to m³/s
            // Use dynamic heat flux if available, otherwise use static value
            AED_REAL current_heat_flux = (heat_pump_dynamic_heat_flux != 0.0) ?
                                        heat_pump_dynamic_heat_flux : heat_pump_heat_flux;
            // ΔT = Q_heat / (ρ × Q_flow × c) Units: W / (kg/m³ × m³/s × J/(kg·K)) = J/s / (kg/s × J/(kg·K)) = K
            temp_change_value = current_heat_flux / (rho0 * flow_rate_m3s * SPHEAT);
            heated_temp = stored_temp - temp_change_value;
            break;
        }
        default: {
            // Default to mode 1 behavior for backward compatibility
            temp_change_value = heat_pump_temp_change;
            heated_temp = stored_temp - temp_change_value;
            break;
        }
    }

    // Get the injection height above lake bottom from the inflow configuration
    AED_REAL inject_height = Inflows[heat_pump_inflow_idx].SubmHeight;

    // Find the layer at injection height above lake bottom
    int Layer_inject;
    for (Layer_inject = botmLayer; Layer_inject <= surfLayer; Layer_inject++) {
        if (Lake[Layer_inject].Height >= inject_height) break;
    }
    if (Layer_inject > surfLayer) Layer_inject = surfLayer;

    // Calculate density of injected water
    AED_REAL inject_density = calculate_density(heated_temp, stored_salt);

    // Directly inject into the lake layer
    // This combines the injected water properties with the existing layer
    Lake[Layer_inject].Temp = combine(Lake[Layer_inject].Temp, Lake[Layer_inject].LayerVol, Lake[Layer_inject].Density,
                                      heated_temp, stored_flow_rate, inject_density);
    Lake[Layer_inject].Salinity = combine(Lake[Layer_inject].Salinity, Lake[Layer_inject].LayerVol, Lake[Layer_inject].Density,
                                          stored_salt, stored_flow_rate, inject_density);

    // Inject WQ variables
    if (Num_WQ_Vars > 0 && WQ_Vars != NULL) {
        for (int wqidx = 0; wqidx < Num_WQ_Vars; wqidx++) {
            _WQ_Vars(wqidx, Layer_inject) = combine_vol(_WQ_Vars(wqidx, Layer_inject), Lake[Layer_inject].LayerVol,
                                                        stored_WQ[wqidx], stored_flow_rate);
        }
    }

    // Update layer density after mixing
    Lake[Layer_inject].Density = calculate_density(Lake[Layer_inject].Temp, Lake[Layer_inject].Salinity);

    // Add the volume back to the layer (mass conservation: outflow removed it, now we add it back)
    Lake[Layer_inject].LayerVol = Lake[Layer_inject].LayerVol + stored_flow_rate;

    // Update cumulative volumes
    Lake[botmLayer].Vol1 = Lake[botmLayer].LayerVol;
    if (surfLayer != botmLayer) {
        for (int j = (botmLayer + 1); j <= surfLayer; j++) {
            Lake[j].Vol1 = Lake[j-1].Vol1 + Lake[j].LayerVol;
        }
    }

    // Recalculate layer heights from volumes (required for consistency)
    resize_internals(2, botmLayer);

    // Status output for monitoring
    print_counter++;
    if (print_counter % 200 == 0) {
        // Show extraction info
        if (heat_pump_outflow_idx >= 0 && heat_pump_outflow_idx < NumOut) {
            const char* outflow_type = (Outflows[heat_pump_outflow_idx].SubmElevDynamic) ? "dynamic" : "static";
            printf("Heat pump extracting at jday %d: Q=%8.4f m³/d, elev=%5.1f m[%s], T=%5.1f°C\n",
                   stored_jday, stored_flow_rate, stored_Drawheight, outflow_type, stored_temp);
        } else {
            printf("Heat pump extracting at jday %d: Q=%8.4f m³/d, elev=%5.1f m[UNKNOWN], T=%5.1f°C\n",
                   stored_jday, stored_flow_rate, stored_Drawheight, stored_temp);
        }

        // Show injection info
        if (heat_pump_switch == 1) {
            printf("Heat pump injecting  at jday %d: Q=%8.4f m³/d, elev=%5.1f m, T=%5.1f°C (ΔT=%+6.2f°C)\n",
                   stored_jday, stored_flow_rate, inject_height, heated_temp, -temp_change_value);
        } else if (heat_pump_switch == 2) {
            AED_REAL current_heat_flux = (heat_pump_dynamic_heat_flux != 0.0) ?
                                        heat_pump_dynamic_heat_flux : heat_pump_heat_flux;
            printf("Heat pump injecting  at jday %d: Q=%8.4f m³/d, elev=%5.1f m, T=%5.1f°C (ΔT=%+6.3f°C, %5.0fW)\n",
                   stored_jday, stored_flow_rate, inject_height, heated_temp, -temp_change_value, current_heat_flux);
        }
    }

    // Clear stored data after injection to prevent double injection
    stored_flow_rate = 0.0;
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

        // Configure heat pump inflow to use plunge dynamics (NOT submerged)
        // This allows cooler/denser water to find its neutral buoyancy level
        if (heat_pump_inflow_idx >= 0 && heat_pump_inflow_idx < NumInf) {
            Inflows[heat_pump_inflow_idx].SubmFlag        = FALSE;  // Use plunge dynamics
            Inflows[heat_pump_inflow_idx].SubmHeightDynamic = FALSE;
        }

        // For outflow: respect the nml configuration (don't force dynamic)
        // If you want static extraction, ensure elev_idx_outflow is NOT set in .nml
        if (heat_pump_outflow_idx >= 0 && heat_pump_outflow_idx < NumOut) {
            // Debug: show what was read from nml
            printf("Heat pump outflow %d: SubmElevDynamic = %s (from .nml)\n",
                   heat_pump_outflow_idx,
                   Outflows[heat_pump_outflow_idx].SubmElevDynamic ? "TRUE (dynamic)" : "FALSE (static)");
        }
    } else {
        printf("Heat pump disabled (heat_pump_switch = %d)\n", heat_pump_switch);
    }
}
