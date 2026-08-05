/******************************************************************************
 *                                                                            *
 * glm_oxygenation.c                                                          *
 *                                                                            *
 * Artificial oxygenation / aeration system for GLM.  See glm_oxygenation.h   *
 * for an overview of the three approaches.                                   *
 *                                                                            *
 *     AquaticEcoDynamics (AED) Group                                         *
 *     The University of Western Australia                                    *
 *                                                                            *
 *  Copyright 2013-2026 : The University of Western Australia                 *
 *                                                                            *
 *  This file is part of GLM (General Lake Model)                             *
 *                                                                            *
 *  GLM is free software: you can redistribute it and/or modify               *
 *  it under the terms of the GNU General Public License as published by      *
 *  the Free Software Foundation, either version 3 of the License, or         *
 *  (at your option) any later version.                                       *
 *                                                                            *
 ******************************************************************************/

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "glm.h"
#include "glm_types.h"
#include "glm_const.h"
#include "glm_globals.h"
#include "glm_util.h"
#include "glm_mixu.h"
#include "glm_flow.h"
#include "aed_csv.h"
#include "glm_wqual.h"
#include "glm_oxygenation.h"

/*----------------------------------------------------------------------------*
 * CSV bookkeeping for the direct-addition devices (Approach 2) and the        *
 * recirculation time-series (Approach 3).                                     *
 * A handle of -1 means "no CSV; use the constant values from the namelist".   *
 *----------------------------------------------------------------------------*/
typedef struct _oxy_csv_ {
    int csv;        //# CSV file handle (-1 = none)
    int time_idx;   //# time column index
    int load_idx;   //# 'load' column (mode 1, and recirc O2 mass/day)
    int flow_idx;   //# 'flow' column (mode 2 flow, and recirc flow)
    int oxy_idx;    //# 'oxy'  column (mode 2 concentration)
    int elev_idx;   //# 'elev' column (optional, dynamic height)
} oxy_csv_t;

static oxy_csv_t oxy_files[MaxInf];
static oxy_csv_t recirc_file = { -1, 0, -1, -1, -1, -1 };
static int files_initialised = FALSE;

/*----------------------------------------------------------------------------*/
static void clear_csv_handles()
{
    int i;
    if ( files_initialised ) return;
    for (i = 0; i < MaxInf; i++) {
        oxy_files[i].csv = -1;
        oxy_files[i].time_idx = 0;
        oxy_files[i].load_idx = -1;
        oxy_files[i].flow_idx = -1;
        oxy_files[i].oxy_idx  = -1;
        oxy_files[i].elev_idx = -1;
    }
    files_initialised = TRUE;
}


/******************************************************************************
 * Locate the layer whose top is at or above a given height above the bottom. *
 ******************************************************************************/
static int layer_at_height(AED_REAL height)
{
    int L;
    for (L = botmLayer; L <= surfLayer; L++)
        if ( Lake[L].Height >= height ) break;
    if ( L > surfLayer ) L = surfLayer;
    return L;
}


/******************************************************************************
 * Resolve the O2 variable index and validate the configuration.             *
 ******************************************************************************/
void init_oxygenation(void)
{
    //# Ensure CSV handles default to -1 ("no file"). Guarded no-op if
    //# oxy_open_file()/oxy_open_recirc_file() already ran.
    clear_csv_handles();

    if ( oxygenation_mode <= 0 ) return;

    //# Resolve the dissolved-oxygen WQ variable index from its name.
    if ( oxy_name == NULL ) oxy_name = "OXY_oxy";
    size_t nl = strlen(oxy_name);
    oxy_o2_idx = wq_var_index_c(oxy_name, &nl);
    if ( oxy_o2_idx < 0 ) {
        fprintf(stderr, "ERROR: oxygenation: unknown oxygen variable '%s' "
                        "(is the WQ model enabled and does it define it?)\n", oxy_name);
        exit(1);
    }

    //# Validate direct-addition devices (modes 1 & 2).
    if ( oxygenation_mode == 1 || oxygenation_mode == 2 ) {
        if ( oxy_num < 0 || oxy_num > MaxInf ) {
            fprintf(stderr, "ERROR: oxygenation: oxy_num (%d) out of range [0,%d]\n",
                    oxy_num, MaxInf);
            exit(1);
        }
    }

    //# Validate the recirculation heights (mode 3).
    if ( oxygenation_mode == 3 ) {
        if ( oxy_recirc_withdraw_height < 0.0 || oxy_recirc_return_height < 0.0 ) {
            fprintf(stderr, "ERROR: oxygenation: recirc withdraw/return heights must be >= 0\n");
            exit(1);
        }
    }
}


/******************************************************************************
 * Print a summary of the configured oxygenation system.                      *
 ******************************************************************************/
void check_oxygenation_config(void)
{
    int i;

    if ( oxygenation_mode <= 0 ) {
        printf("Oxygenation disabled\n");
        return;
    }

    printf("Oxygenation enabled (mode %d): O2 variable '%s' (index %d)\n",
           oxygenation_mode, (oxy_name != NULL) ? oxy_name : "OXY_oxy", oxy_o2_idx);

    if ( oxygenation_mode == 1 || oxygenation_mode == 2 ) {
        for (i = 0; i < oxy_num; i++) {
            if ( oxy_input_type[i] == 2 )
                printf("  device %d: input_type 2 (flow x conc) at %.2f m, flow=%.4g m3/d, conc=%.4g%s\n",
                       i, oxy_height[i], oxy_flow[i], oxy_conc[i],
                       (oxy_files[i].csv >= 0) ? " [CSV]" : "");
            else
                printf("  device %d: input_type 1 (mass rate) at %.2f m, load=%.4g /day%s\n",
                       i, oxy_height[i], oxy_load[i],
                       (oxy_files[i].csv >= 0) ? " [CSV]" : "");
        }
        if ( oxy_max > 0.0 )
            printf("  O2 concentration capped at %.4g\n", oxy_max);
    }

    if ( oxygenation_mode == 3 )
        printf("  recirculation: withdraw %.2f m -> return %.2f m, flow=%.4g m3/s, load=%.4g /day%s\n",
               oxy_recirc_withdraw_height, oxy_recirc_return_height,
               oxy_recirc_flow, oxy_recirc_add,
               (recirc_file.csv >= 0) ? " [CSV]" : "");
}


/******************************************************************************
 * Open the CSV time-series for a direct-addition device (Approach 2).        *
 ******************************************************************************/
void oxy_open_file(int idx, const char *fname, const char *timefmt)
{
    clear_csv_handles();

    if ( fname == NULL || fname[0] == '\0' ) return;  //# no CSV for this device
    if ( idx < 0 || idx >= MaxInf ) return;

    if ( (oxy_files[idx].csv = open_csv_input(fname, timefmt)) < 0 ) {
        fprintf(stderr, "ERROR: oxygenation: failed to open '%s'\n", fname);
        exit(1);
    }
    oxy_files[idx].time_idx = 0;  //# time/date is enforced to be the first column
    oxy_files[idx].load_idx = find_csv_var(oxy_files[idx].csv, "load");
    oxy_files[idx].flow_idx = find_csv_var(oxy_files[idx].csv, "flow");
    oxy_files[idx].oxy_idx  = find_csv_var(oxy_files[idx].csv, "oxy");
    oxy_files[idx].elev_idx = find_csv_var(oxy_files[idx].csv, "elev");
}


/******************************************************************************
 * Open the CSV time-series for the recirculation (Approach 3).               *
 * Columns: 'flow' (m3/s) and 'load' (O2 mass/day).                           *
 ******************************************************************************/
void oxy_open_recirc_file(const char *fname, const char *timefmt)
{
    clear_csv_handles();

    if ( fname == NULL || fname[0] == '\0' ) return;

    if ( (recirc_file.csv = open_csv_input(fname, timefmt)) < 0 ) {
        fprintf(stderr, "ERROR: oxygenation: failed to open '%s'\n", fname);
        exit(1);
    }
    recirc_file.time_idx = 0;
    recirc_file.flow_idx = find_csv_var(recirc_file.csv, "flow");
    recirc_file.load_idx = find_csv_var(recirc_file.csv, "load");
    if ( recirc_file.flow_idx < 0 )
        fprintf(stderr, "WARNING: oxygenation: 'flow' column not found in recirc CSV\n");
    if ( recirc_file.load_idx < 0 )
        fprintf(stderr, "WARNING: oxygenation: 'load' column not found in recirc CSV\n");
}


/******************************************************************************
 * Read today's CSV values for CSV-driven devices and the recirc stream.      *
 ******************************************************************************/
void read_daily_oxygenation(int julian)
{
    int i, csv;

    //# Only mode 2 (CSV-driven direct addition) reads per-device time series.
    if ( oxygenation_mode == 2 ) {
        for (i = 0; i < oxy_num; i++) {
            if ( (csv = oxy_files[i].csv) < 0 ) continue;  //# device without a CSV
            find_day(csv, oxy_files[i].time_idx, julian);

            if ( oxy_input_type[i] == 2 ) {
                if ( oxy_files[i].flow_idx >= 0 )
                    oxy_flow[i] = get_csv_val_r(csv, oxy_files[i].flow_idx);
                if ( oxy_files[i].oxy_idx >= 0 )
                    oxy_conc[i] = get_csv_val_r(csv, oxy_files[i].oxy_idx);
            } else {
                if ( oxy_files[i].load_idx >= 0 )
                    oxy_load[i] = get_csv_val_r(csv, oxy_files[i].load_idx);
            }
            if ( oxy_files[i].elev_idx >= 0 )
                oxy_height[i] = get_csv_val_r(csv, oxy_files[i].elev_idx);
        }
    }

    if ( oxygenation_mode == 3 && (csv = recirc_file.csv) >= 0 ) {
        find_day(csv, recirc_file.time_idx, julian);
        if ( recirc_file.flow_idx >= 0 )
            oxy_recirc_flow = get_csv_val_r(csv, recirc_file.flow_idx);
        if ( recirc_file.load_idx >= 0 )
            oxy_recirc_add = get_csv_val_r(csv, recirc_file.load_idx);
    }
}


/******************************************************************************
 * Approaches 1 & 2: add O2 mass to the target layer of each device.          *
 * No water is added; only the O2 concentration of the layer rises.           *
 * Returns the total O2 mass added across all devices.                        *
 ******************************************************************************/
AED_REAL do_oxygenation(AED_REAL day_fraction)
{
    int i, L;
    AED_REAL delta_mass, delta_conc, total = zero;

    if ( (oxygenation_mode != 1 && oxygenation_mode != 2) || oxy_o2_idx < 0 || Num_WQ_Vars <= 0 )
        return zero;

    for (i = 0; i < oxy_num; i++) {
        //# Mass of O2 delivered today, per the selected input mode.
        if ( oxy_input_type[i] == 2 )
            delta_mass = oxy_flow[i] * oxy_conc[i] * day_fraction;  //# flow x conc
        else
            delta_mass = oxy_load[i] * day_fraction;                //# direct mass rate

        if ( delta_mass <= zero ) continue;

        L = layer_at_height(oxy_height[i]);
        if ( Lake[L].LayerVol <= zero ) continue;

        //# Add mass without adding water: concentration rises, volume unchanged.
        delta_conc = delta_mass / Lake[L].LayerVol;
        _WQ_Vars(oxy_o2_idx, L) += delta_conc;

        //# Optional saturation cap.
        if ( oxy_max > zero && _WQ_Vars(oxy_o2_idx, L) > oxy_max )
            _WQ_Vars(oxy_o2_idx, L) = oxy_max;

        total += delta_mass;
    }
    return total;
}


/******************************************************************************
 * Approach 3 (self-contained): withdraw water at oxy_recirc_withdraw_height,  *
 * add O2, and return it at oxy_recirc_return_height. The withdrawal uses      *
 * GLM's selective-withdrawal physics (oxy_recirc_extract in glm_flow.c) and   *
 * the return re-inserts exactly the withdrawn volume, so water is conserved.  *
 * Returns the O2 mass added.                                                  *
 ******************************************************************************/
AED_REAL oxy_do_recirculation(AED_REAL day_fraction)
{
    int wqidx, L, j;
    AED_REAL want_vol, drawn, inject_density, total = zero;
    AED_REAL cap_temp = 0.0, cap_salt = 0.0;
    AED_REAL cap_wq[MaxVars];

    if ( oxygenation_mode != 3 ) return zero;

    //# Volume to recirculate this step (oxy_recirc_flow is m3/s).
    want_vol = oxy_recirc_flow * SecsPerDay * day_fraction;
    if ( want_vol <= zero ) return zero;

    //# 1. Withdraw at the withdrawal height, capturing T/S/WQ of that water.
    drawn = oxy_recirc_extract(oxy_recirc_withdraw_height, want_vol,
                               &cap_temp, &cap_salt, cap_wq);
    if ( drawn <= zero ) return zero;

    //# 2. Add O2 to the captured stream (mass/day scaled to this step).
    if ( oxy_o2_idx >= 0 && Num_WQ_Vars > 0 ) {
        cap_wq[oxy_o2_idx] += (oxy_recirc_add * day_fraction) / drawn;
        total = oxy_recirc_add * day_fraction;
    }

    //# 3. Re-inject the captured water at the return height (mass conserved).
    L = layer_at_height(oxy_recirc_return_height);
    inject_density = calculate_density(cap_temp, cap_salt);

    Lake[L].Temp = combine(Lake[L].Temp, Lake[L].LayerVol, Lake[L].Density,
                           cap_temp, drawn, inject_density);
    Lake[L].Salinity = combine(Lake[L].Salinity, Lake[L].LayerVol, Lake[L].Density,
                               cap_salt, drawn, inject_density);
    if ( Num_WQ_Vars > 0 ) {
        for (wqidx = 0; wqidx < Num_WQ_Vars; wqidx++)
            _WQ_Vars(wqidx, L) = combine_vol(_WQ_Vars(wqidx, L), Lake[L].LayerVol,
                                             cap_wq[wqidx], drawn);
    }
    Lake[L].Density = calculate_density(Lake[L].Temp, Lake[L].Salinity);
    Lake[L].LayerVol = Lake[L].LayerVol + drawn;

    //# Update cumulative volumes and re-derive layer heights.
    Lake[botmLayer].Vol1 = Lake[botmLayer].LayerVol;
    if ( surfLayer != botmLayer ) {
        for (j = (botmLayer + 1); j <= surfLayer; j++)
            Lake[j].Vol1 = Lake[j-1].Vol1 + Lake[j].LayerVol;
    }
    resize_internals(2, botmLayer);

    return total;
}
