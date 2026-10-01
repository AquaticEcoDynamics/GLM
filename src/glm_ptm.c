/******************************************************************************
 *                                                                            *
 * glm_ptm.c                                                                  *
 *                                                                            *
 * Contains Particle Tracking Model                                           *
 *                                                                            *
 * Developed by :                                                             *
 *                                                                            *
 * Copyright 2024-2026 : The University of Western Australia                  *
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
#include <math.h>
#include <time.h>       /* time */

#include "glm.h"

#include "glm_types.h"      // check if these are needed
#include "glm_const.h"
#include "glm_globals.h"

#include "glm_ptm.h"

#include "glm_util.h"
#include "glm_ncdf.h"
#include "glm_wqual.h"

            // Cell suspended in water column
#define WATER  0
            // Maybe forming a bed layer; set when hits bottom; turn off when resuspended
#define BED    1
            // Maybe forming a scum layer
#define SCUM   2
            // Maybe exiting the lake
#define EXIT   3

#define STAT   0
#define IDX2   1
#define IDX3   2
#define LAYR   3
#define FLAG   4
#define PTID   5
#define GRP    6   // which phyto group/species this particle belongs to (aed_phyto_abm.F90 GRP=7, 1-indexed)

#define MASS   0
#define DIAM   1
#define DENS   2
#define VVEL   3
#define HGHT   4

// Width of the host-owned "environment" block above (MASS..HGHT), i.e. the offset at which
// the AED particle state variables begin. Must match n_ptm_env in libaed-api/src/aed_ptm.F90,
// which allocates ptm_env as (groups, particles, n_ptm_env + n_ptm_vars) - the array this
// file aliases. Valid third indices are therefore 0 .. PTM_ENV_NVARS+Num_PTM_Vars-1.
#define PTM_ENV_NVARS 5

#define PTM_STAT_NVARS 7   // STAT,IDX2,IDX3,LAYR,FLAG,PTID,GRP - n_ptm_istat in libaed-api/src/aed_ptm.F90
#define PTM_OASIM_NBANDS 7

/* 2026-07 fix: the two PTID allocators must stay disjoint. This file hands out seed ids
 * from pg*max_particle_num + p + 1 (and a monotonic reseed counter floored at
 * num_particle_groups*max_particle_num); aed_split_particle_phyto_abm() in
 * libaed-water/src/aed_phyto_abm.F90 hands out split-born ids from PTID_SPLIT_BASE
 * upward. Keep the value below identical to the PTID_SPLIT_BASE parameter there.
 * ptm_init_glm() checks the invariant at startup rather than merely asserting it. */
#define PTID_SPLIT_BASE 1000000

#define _PTM_Stat(grp,part,var) PTM_Stat[_IDX_3d(num_particle_groups,max_particle_num,PTM_STAT_NVARS,grp,part,var)]
#define _PTM_Vars(grp,part,var) PTM_Vars[_IDX_3d(num_particle_groups,max_particle_num,(PTM_ENV_NVARS + Num_PTM_Vars),grp,part,var)]

AED_REAL get_particle_density(AED_REAL particle_density);
AED_REAL get_particle_diameter(AED_REAL particle_diameter);
AED_REAL get_settling_velocity(AED_REAL settling_velocity);
AED_REAL move_particle(AED_REAL dt_secs, AED_REAL Height, AED_REAL K_z, AED_REAL K_prime_z, AED_REAL vvel, AED_REAL rand_draw);
static AED_REAL draw_height_in_range(AED_REAL lo, AED_REAL hi);
static int ptm_next_free_slot(int pg, int max_particle_num);
static AED_REAL ptm_var_or_zero(int grp, int part, int var);
static AED_REAL ptm_state_or_zero(int grp, int part, int fidx);
static int ptm_total_particles(void);
static int ptm_group_count(int grp);
static int ptm_group_int(int *values, int grp, int fallback);
static int ptm_reseed_count(int *values, int grp, int fallback, const char *label);
static void ptm_print_reseed_config(void);
static void ptm_log_reseed_event(int pg, int active, int min_active, int target_active, int need);
static AED_REAL ptm_group_real(AED_REAL *values, int grp, AED_REAL fallback);
static int ptm_pigment_id(int grp);
static int ptm_behavior_id(int grp);
static int ptm_vertical_mode(AED_REAL vvel_m_per_s);
static AED_REAL ptm_light_preference_nm(int pigment_id);
static AED_REAL ptm_best_light_for_pigment(int pigment_id, AED_REAL light_440, AED_REAL light_490,
                                           AED_REAL light_550, AED_REAL light_565, AED_REAL light_590,
                                           AED_REAL light_620, AED_REAL light_665);
static AED_REAL ptm_particle_par(int grp, int part);
static void ptm_particle_oasim_bands(int layer, AED_REAL bands[PTM_OASIM_NBANDS], int *sample_status);
void ptm_update_state_indices(void);
static void ptm_addparticles_group(int pg, int new_particles, int max_particle_num,
                                   AED_REAL upper_height, AED_REAL lower_height);
static void ptm_addparticles_group_ex(int pg, int new_particles, int max_particle_num,
                                      AED_REAL upper_height, AED_REAL lower_height,
                                      int force_reseed_id);
static void ptm_apply_reseeding(void);
static const char *ptm_group_name(int grp);
static void ptm_group_names_attr(char *buf, size_t len);
extern void aed_sample_oasim_particle_bands(int layer, AED_REAL *bands, int nbands, int *status);
extern void aed_sample_oasim_par_kd(int layer, AED_REAL *par_center, AED_REAL *kd, int *status);
extern void aed_sample_oasim_afac(int layer, int iop, AED_REAL *afac, int *status);

/* Tier 1D: PTM state-variable index getters from the ABM module. */
extern int aed_phyto_abm_ipar_index(void);
extern int aed_phyto_abm_ibuoy_rate_index(void);
extern int aed_phyto_abm_ilightdose_index(void);
extern int aed_phyto_abm_is_act_index(void);
extern int aed_phyto_abm_ibuoy_I_index(void);
extern int aed_phyto_abm_ibuoy_stat_index(void);
extern void aed_phyto_abm_buoy_params(int grp, int *imodel, double *prm);
extern double aed_water_viscosity_c(double temp);
extern int aed_phyto_abm_ipar_scalar_depth_index(void);
extern int aed_phyto_abm_ipar_layer_center_index(void);
extern int aed_phyto_abm_ipar_depth_factor_index(void);
extern int aed_phyto_abm_ipar_scalar_status_index(void);
/* Phase A: per-band Kd sampler + index getters for the depth bands / absorbed PAR */
extern void aed_sample_oasim_band_kd(int layer, AED_REAL *kds, int nbands, int *status);
extern int aed_phyto_abm_ipar_depth_band_index(int band);
extern int aed_phyto_abm_ipar_absorbed_index(void);
/* Ranjbar NPQ: cumulative light dose CL. Indexed via a getter, not a VVEL+n offset. */
extern int aed_phyto_abm_ipar_cldose_index(void);
/* Lagrangian daily accumulators. Resolved by NAME through these accessors - never by a
 * hardcoded VVEL+N offset, which is what silently corrupted particle_age_d and six of its
 * neighbours when PTM registrations were inserted mid-list. */
extern int aed_phyto_abm_idose_scal_l_index(void);
extern int aed_phyto_abm_idose_eff_l_index(void);
extern int aed_phyto_abm_idprodC_l_index(void);
extern int aed_phyto_abm_idday_index(void);

/* Water Research 2022 colony dynamics: per-particle effective colony size. */
extern int aed_phyto_abm_ipar_fcol_index(void);

/* Tier 1D: cached offsets for particle-state variables (1-based Fortran -> 0-based C). */
static int ptm_idx_par = -1;
static int ptm_idx_buoy_rate = -1;
static int ptm_idx_lightdose = -1;
static int ptm_idx_s_act = -1;
static int ptm_idx_buoy_I = -1;      /* irradiance the integrator actually used */
static int ptm_idx_buoy_stat = -1;   /* 0 = OASIM, 1 = Lake[].Light fallback     */
/* Wallace & Hamilton (Limnol. Oceanogr. 44:273, 1999; J. Plankton Res. 22:1127, 2000) lagged
 * buoyancy regulation. Integrated HERE, at the particle sub-step, and not in aed_phyto_abm.F90,
 * because AED runs once per hydrodynamic step (3600 s) while the response time being modelled is
 * 20 min (1200 s): at hourly sampling the lag relaxes 95% per step and the mechanism vanishes.
 * The spec, with the papers reproduced, is analysis/test_buoyancy_wh.py. */
#define MAX_PTM_GROUPS_ 16          /* guarded against num_particle_groups at init */
static int      buoy_model[MAX_PTM_GROUPS_] = {0};
static AED_REAL buoy_prm[MAX_PTM_GROUPS_][11];
static int      buoy_any = 0;
static AED_REAL *buoy_par_center = NULL;   /* per-layer OASIM PAR, cached once per dt  */
static AED_REAL *buoy_kd = NULL;           /* per-layer OASIM Kd,  cached once per dt  */
static int      *buoy_stat = NULL;         /* per-layer: 0 = OASIM, 1 = Lake[].Light fallback */
static AED_REAL  buoy_I_dark = 16.41;      /* min I_dark over active groups; set at init      */
static AED_REAL *buoy_h0 = NULL;           /* height when AED last wrote ip_par_scalar_depth  */
static AED_REAL *buoy_scale = NULL;        /* AED PAR / Lake-based PAR, measured at h0        */
static int       buoy_h0_n = 0;
static AED_REAL  buoy_sfc_last = 0.0;      /* surface PAR at the previous step */
static AED_REAL  buoy_sfc_ratio = 1.0;     /* this step / previous step        */
static int       buoy_cache_n = 0;
static int       buoy_is_day = 0;          /* surface light state, for the dawn reset */
static int       buoy_dawn_reset = 0;      /* set for the step on which dawn occurs   */
static int ptm_idx_par_scalar_depth = -1;
static int ptm_idx_par_layer_center = -1;
static int ptm_idx_par_depth_factor = -1;
static int ptm_idx_par_scalar_status = -1;
static int ptm_idx_par_depth_band[PTM_OASIM_NBANDS] = {-1,-1,-1,-1,-1,-1,-1};
static int ptm_idx_par_absorbed = -1;
static int ptm_idx_cldose = -1;
static int ptm_idx_dose_scal_l = -1;
static int ptm_idx_dose_eff_l = -1;
static int ptm_idx_dprodC_l = -1;
static int ptm_idx_dday = -1;

static int ptm_idx_fcol = -1;
/* wavelengths in the fixed OASIM order used by aed_sample_oasim_particle_bands */
static const int ptm_band_nm[PTM_OASIM_NBANDS] = {440,490,550,565,590,620,665};

#define PTM_STATE_VAR(fortran_idx) (PTM_ENV_NVARS + ((fortran_idx) - 1))

/* 2026-09-10: the LEGACY PTM state block, resolved from AED instead of hardcoded.
 *
 * These 18 used to be reached by positional offsets - VVEL+3..VVEL+21 when writing the
 * particle outputs, PTM_ENV_NVARS+5..+12 when seeding BGC state. Every one of them was
 * correct, but only because nothing had been registered ahead of them in
 * aed_phyto_abm_define. When something was, the offsets silently repointed and 7 particle
 * outputs were corrupted with no error - the reads stay in bounds, so no bounds check and
 * no compiler warning can see it. Appending new variables made the symptom go away while
 * leaving the coupling in place; this removes the coupling.
 *
 * ORDER IS THE CONTRACT and it mirrors aed_phyto_abm_legacy_ptm_indices in
 * aed_phyto_abm.F90. It is NOT the registration order - that is the thing that changes.
 * A new ABM variable is APPENDED at both ends, never inserted.
 */
enum {
    PAM_L_TEM = 0, PAM_L_NO3, PAM_L_NH4, PAM_L_FRP,
    PAM_L_C,       PAM_L_N,   PAM_L_P,   PAM_L_CHL,
    PAM_L_NUM,     PAM_L_CDIV, PAM_L_TOPT, PAM_L_LNALPHACHL,
    PAM_L_AGE,     PAM_L_BIRTH, PAM_L_CUMPAR, PAM_L_CUMC,
    PAM_L_NDIV,    PAM_L_NPQ,
    PAM_N_LEGACY_C
};
static int ptm_idx_legacy[PAM_N_LEGACY_C];
extern void aed_phyto_abm_legacy_ptm_indices(int n, int *idx);

/* Fortran index of one legacy variable, or 0 when the ABM did not register it. */
#define PTM_LEG_F(which)  (ptm_idx_legacy[which])
/* Its slot in the per-particle state array. Only valid when PTM_LEG_F(which) >= 1. */
#define PTM_LEG(which)    PTM_STATE_VAR(PTM_LEG_F(which))

/* Read one legacy variable, returning 0.0 when it is not registered. Without the guard an
 * unresolved index of 0 would read PTM_STATE_VAR(0) = PTM_ENV_NVARS - 1, which is a valid
 * slot holding an unrelated environment variable: a wrong number rather than an error. */
static AED_REAL ptm_legacy_var(int grp, int part, int which)
{
    if ( PTM_LEG_F(which) < 1 ) return 0.0;
    return ptm_var_or_zero(grp, part, PTM_LEG(which));
}

/*============================================================================*/

//CONSTANTS
int num_particle_grp=1;
int *init_particle_num_by_group = NULL;
int *particle_reseed_min_by_group = NULL;
int *particle_reseed_target_by_group = NULL;
AED_REAL *init_depth_min_by_group = NULL;
AED_REAL *init_depth_max_by_group = NULL;
AED_REAL *particle_density_by_group = NULL;
AED_REAL *particle_diameter_by_group = NULL;
AED_REAL *settling_velocity_by_group = NULL;
char **particle_group_names = NULL;
char **particle_phyto_links = NULL;
char **particle_pigment_type = NULL;
int particle_group_names_n = 0;    // list lengths, set by glm_init.c
int particle_pigment_type_n = 0;
AED_REAL init_depth_min=0.0;
AED_REAL init_depth_max=2.0;
AED_REAL ptm_time_step=60.0;   // PTM substep duration in SECONDS; must evenly
                               // divide the host timestep dt (see do_ptm_update())
AED_REAL ptm_diffusivity=1e-6;
/* Buoyancy density is updated every ptm_buoy_substep'th particle sub-step. At the defaults
 * (ptm_time_step 7.5 s, substep 16) that is every 120 s, i.e. 10 samples per 20 min response
 * time. Convergence is checked by running 8 / 16 / 32. */
CINTEGER ptm_buoy_substep = 16;
/* 2026-09-10 (E5) RESTART-EXACT RNG. The particle model draws from libc rand(), whose
 * internal state cannot be read out portably, so a resumed run re-seeded and its random
 * stream diverged from the uninterrupted run's from the first step. Every draw now goes
 * through glm_rand(), which counts calls; the count is written to the restart file and a
 * resume re-seeds with particle_random_seed and replays that many draws, which restores
 * the generator's state exactly without changing the generator. Same numbers, same order,
 * for every run that never restarts. */
extern CINTEGER particle_random_seed;   /* defined below, with the other namelist defaults */
static long long ptm_rng_calls = 0;
static int glm_rand(void) { ptm_rng_calls++; return rand(); }
double ptm_rng_call_count(void) { return (double)ptm_rng_calls; }
void ptm_rng_restore(double n)
{
    long long i, m = (long long)n;
    srand((unsigned int)particle_random_seed);
    ptm_rng_calls = 0;
    for (i = 0; i < m; i++) glm_rand();
}
/* RB 2026-09-15: the buoyancy light history - the dawn hysteresis state and the previous step's surface PAR -
 * lives in statics from one step to the next, so the restart file must carry it. The start-up values (night, no
 * previous light) are right only for a start or resume in darkness: a daytime resume detected a false dawn on its
 * first step and zeroed every buoy_model-1 particle's light dose (tools/test_restart_continuity.py, 17:00 split). */
void ptm_buoy_light_state_get(double *is_day, double *sfc_last)
{ *is_day = (double)buoy_is_day; *sfc_last = (double)buoy_sfc_last; }
void ptm_buoy_light_state_set(double is_day, double sfc_last)
{ buoy_is_day = (is_day > 0.5) ? 1 : 0; buoy_sfc_last = (AED_REAL)sfc_last; }
/* Optional diagnostic-probe maintenance. Disabled by default. When enabled,
 * groups below particle_reseed_min are topped back up to particle_reseed_target
 * using the same depth ranges as initial seeding. This preserves Option A:
 * probes are replenished as instruments and still do not perturb Eulerian AED. */
CINTEGER particle_reseed_enabled = 0;
CINTEGER particle_reseed_min = 0;
CINTEGER particle_reseed_target = 0;
static int ptm_reseed_log_count[128] = {0};
static int ptm_reseed_log_suppressed[128] = {0};
static int ptm_reseed_min_warned[128] = {0};
static int ptm_reseed_target_warned[128] = {0};

// Boundary condition types
#define BC_CLAMP      1
#define BC_REFLECTIVE 2

// VARIABLES
LOGICAL sed_deactivation = FALSE;
CINTEGER num_particle_groups = 1;

/* Seed for the RNG used by particle placement and the vertical random walk, and (via the
 * matching BIND(C) declaration in glm_api_aed.F90) for the Fortran intrinsic RNG used by
 * the ABM's stochastic mortality and trait mutation - so one namelist knob drives both.
 *
 * A fixed non-zero default makes runs REPRODUCIBLE out of the box. Previously glm_main.c
 * seeded from time(NULL) before the namelist was even read, while the Fortran RNG was
 * never seeded at all: the model was half wall-clock-random and half frozen, and no two
 * runs were comparable. Set to 0 in &particles to seed from the wall clock instead, which
 * is only wanted when deliberately generating ensemble members. */
CINTEGER particle_random_seed = 20250723;
// Particle boundary conditions (1 = clamping, 2 = reflective)
int upper_boundary_cond = BC_CLAMP;   // surface
// lower_boundary_cond also governs what happens when the area-ratio settling
// probability triggers (the particle "hits a2" - see the (a1-a2)/a1 check in
// do_ptm_update()). CLAMP (default): today's behavior - the particle is
// flagged BED but its height is left where it already is (a2's height).
// REFLECTIVE: the particle instead bounces back to prev_height, undoing this
// substep's downward step, same distance it fell.
int lower_boundary_cond = BC_CLAMP;   // bottom

/*============================================================================*/

/******************************************************************************
 *                                                                            *
 *    Draw a height in [lo, hi], weighted by dMphLevelArea WITHIN that        *
 *    range, for redistributing a particle inside a mixing zone. A physically *
 *    thick mixed layer can span a real change in lake cross-sectional area   *
 *    from its top to its bottom (e.g. a layer straddling the thermocline in  *
 *    a basin that narrows with depth); drawing uniformly in HEIGHT within    *
 *    the layer implicitly assumes constant particle density per unit height  *
 *    regardless of that area change, i.e. treats the layer as a box. This    *
 *    instead biases the draw toward heights with more cross-sectional area,  *
 *    via the same inverse-CDF technique, restricted to [lo, hi] instead of   *
 *    the whole water column - and re-weighting the partial 0.1 m bins at     *
 *    both the lo and hi ends by how much of each is actually inside         *
 *    [lo, hi], not just a single partial bin at one end.                    *
 *                                                                            *
 ******************************************************************************/
static AED_REAL draw_height_in_range(AED_REAL lo, AED_REAL hi)
{
    static AED_REAL *cdf = NULL;
    static AED_REAL *bnd = NULL;
    int i_lo, i_hi, n, k;
    AED_REAL total = 0.0;
    AED_REAL u, bin_lo_cdf, frac, result;

    if (hi <= lo) return lo;

    if (cdf == NULL) {
        cdf = malloc(Nmorph * sizeof(AED_REAL));
        bnd = malloc((Nmorph + 1) * sizeof(AED_REAL));
    }

    i_lo = (int)(lo * 10.0);
    i_hi = (int)(hi * 10.0);
    if (i_lo < 0) i_lo = 0;
    if (i_hi >= Nmorph) i_hi = Nmorph - 1;
    if (i_lo > i_hi) i_lo = i_hi;
    n = i_hi - i_lo + 1;

    // Bin boundaries: bnd[0]=lo, bnd[n]=hi, interior ones on the 0.1 m grid - so
    // the first and last bins are only as wide as the part of [lo,hi] they cover.
    bnd[0] = lo;
    for (k = 1; k < n; k++) bnd[k] = (AED_REAL)(i_lo + k) / 10.0;
    bnd[n] = hi;

    for (k = 0; k < n; k++) {
        AED_REAL w = dMphLevelArea[i_lo + k];   // area per FULL 0.1 m band
        if (w < 0.0) w = 0.0;
        w *= (bnd[k+1] - bnd[k]) * 10.0;        // scale to this bin's actual (sub-)width
        if (w < 0.0) w = 0.0;
        total += w;
        cdf[k] = total;
    }

    if (total <= 0.0) {
        // Degenerate (zero-area range): fall back to uniform, same as the old behavior.
        u = (AED_REAL)glm_rand() / RAND_MAX;
        return lo + u * (hi - lo);
    }

    u = ((AED_REAL)glm_rand() / RAND_MAX) * total;
    for (k = 0; k < n; k++)
        if (u <= cdf[k]) break;
    if (k >= n) k = n - 1;

    bin_lo_cdf = (k == 0) ? 0.0 : cdf[k-1];
    frac = (cdf[k] > bin_lo_cdf) ? (u - bin_lo_cdf) / (cdf[k] - bin_lo_cdf) : 0.0;
    result = bnd[k] + frac * (bnd[k+1] - bnd[k]);

    if (result < lo) result = lo;
    if (result > hi) result = hi;
    return result;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *    Free-slot search within one particle group's pool (slots 0 ..            *
 *    max_particle_num-1 of group pg; each group keeps its own cursor).       *
 *                                                                            *
 *    A slot is "free" when STAT==0 && FLAG==EXIT(3). Every place that sets   *
 *    STAT=0 pairs it with FLAG=EXIT(3) EXCEPT do_ptm_update()'s              *
 *    sed_deactivation branch, which leaves FLAG=BED and so deliberately      *
 *    keeps those slots out of circulation.                                   *
 *                                                                            *
 *    The particle array is its own bookkeeping: there is no separate free    *
 *    list to keep in sync with it, so the two cannot disagree. (A FIFO queue *
 *    shared with aed_phyto_abm.F90 over BIND(C) was tried and removed - it   *
 *    made allocation O(1) but coupled the two languages, and forced the      *
 *    Fortran death paths to recover a global slot number from PTID because   *
 *    their loop index is layer-local, which caused a real bug.)              *
 *                                                                            *
 *    ptm_search_cursor resumes the scan where the previous one stopped,      *
 *    wrapping once. It matters because allocation happens in bursts -        *
 *    init_particle_num slots in a single startup call, and one per new       *
 *    particle on every inflow step - and restarting from 0 each time would   *
 *    make a burst of K cost O(K*max_particle_num). It is only ever a hint:   *
 *    every candidate is still tested against STAT/FLAG below, so a stale     *
 *    cursor cannot return a slot that is not genuinely free.                 *
 *                                                                            *
 ******************************************************************************/
static int ptm_search_cursor[MAX_PTM_GROUPS_] = {0};   // rolling start position per group - see below

static int ptm_next_free_slot(int pg, int max_particle_num)
{
    int k, s;
    int *cur = &ptm_search_cursor[(pg >= 0 && pg < MAX_PTM_GROUPS_) ? pg : 0];

    if (max_particle_num <= 0) return -1;
    if (*cur < 0 || *cur >= max_particle_num)
        *cur = 0;

    for (k = 0; k < max_particle_num; k++) {
        s = *cur;
        (*cur)++;
        if (*cur >= max_particle_num) *cur = 0;
        if (_PTM_Stat(pg,s,STAT) == 0 && _PTM_Stat(pg,s,FLAG) == EXIT)
            return s;
    }
    return -1;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *    This routine initialises the PTM data structure, and sets the initial   *
 *    particle properties and position                                        *
 *                                                                            *
 *                                                                            *
 *                                                                            *
 *    NOTES:                                                                  *
 *      all arrays begin with subscript of zero                               *
 *                                                                            *
 ******************************************************************************/
void ptm_init_glm()
{
//LOCALS
    int p;

    AED_REAL upper_height;
    AED_REAL lower_height;

//  int bla;
    int pg;

/*----------------------------------------------------------------------------*/
//BEGIN

    num_particle_groups = num_particle_grp;
    if (num_particle_groups < 1) num_particle_groups = 1;

    /* 2026-07 fix: enforce the PTID disjointness invariant instead of asserting it in a
     * comment. Seed ids here run up to num_particle_groups*max_particle_num; split-born
     * ids in aed_phyto_abm.F90 start above PTID_SPLIT_BASE. If the seeded range reaches
     * that base the two allocators collide and particle identities become ambiguous. */
    if ((long)num_particle_groups * (long)max_particle_num >= (long)PTID_SPLIT_BASE) {
        fprintf(stderr, "     ERROR: num_particle_grp (%d) * max_particle_num (%d) = %ld "
                        "reaches the split-particle PTID base (%d).\n",
                num_particle_groups, max_particle_num,
                (long)num_particle_groups * (long)max_particle_num, PTID_SPLIT_BASE);
        fprintf(stderr, "            Seeded and split-born particle IDs would collide; "
                        "reduce max_particle_num or num_particle_grp.\n");
        exit(1);
    }

    //NEW allocate AED ptm data structures, and GLM pointers
    if (do_particle_bgc) {
        api_set_glm_ptm(&num_particle_groups,&max_particle_num);
    } else if (Num_PTM_Vars == 0) {
        PTM_Stat = calloc(num_particle_groups * max_particle_num * PTM_STAT_NVARS, sizeof(int));
        PTM_Vars = calloc(num_particle_groups * max_particle_num * PTM_ENV_NVARS, sizeof(AED_REAL));
        if (PTM_Stat == NULL || PTM_Vars == NULL) {
            fprintf(stderr, "     ERROR: unable to allocate physical-only PTM arrays\n");
            exit(1);
        }
    } else {
        api_set_glm_ptm(&num_particle_groups,&max_particle_num);   //_WQ_SET_GLM_PTM
    }

    //printf("_PTM_Stat(0,5000-1,0)  %d \n" ,_PTM_Stat(0,5000-1,0));
    //printf("_PTM_Stat(0,5000-1,1)  %d \n"  ,_PTM_Stat(0,5000-1,1));
    //printf("_PTM_Stat(0,5000-1,2)  %d \n"  ,_PTM_Stat(0,5000-1,2));
    //printf("_PTM_Stat(0,5000-1,3)  %d \n"  ,_PTM_Stat(0,5000-1,3));
    //printf("_PTM_Stat(0,5000,0)  %d \n"  ,_PTM_Stat(0,5000,0));
    //printf("_PTM_Stat(0,5000,1)  %d \n"  ,_PTM_Stat(0,5000,1));
    //printf("_PTM_Stat(0,5000,2)  %d \n"  ,_PTM_Stat(0,5000,2));
    //printf("_PTM_Stat(0,5000,3)  %d \n"  ,_PTM_Stat(0,5000,3));
    //printf("_PTM_Stat(0,5000,4)  %d \n" ,_PTM_Stat(0,5000,4));
    //printf("_PTM_Stat(0,0,0)  %d \n"  ,_PTM_Stat(0,0,0));
    //printf("_PTM_Stat(0,1,1)  %d \n"  ,_PTM_Stat(0,1,1));
    //printf("_PTM_Stat(0,2,2)  %d \n"  ,_PTM_Stat(0,2,2));
    //printf("_PTM_Stat(0,3,3)  %d \n"  ,_PTM_Stat(0,3,3));

    //printf("_PTM_Vars(0,0,0)  %f \n"  ,_PTM_Vars(0,0,0));
    //printf("_PTM_Vars(0,0,0)  %f \n"  ,_PTM_Vars(0,0,1));
    //printf("_PTM_Vars(0,0,0)  %f \n"  ,_PTM_Vars(0,0,2));
    //printf("_PTM_Vars(0,2,2)  %f \n"  ,_PTM_Vars(0,2,2));
    //NEW initialise (integer) status array for all particles to 0
    for (pg = 0; pg < num_particle_groups; pg++) {
      for (p = 0; p < max_particle_num; p++) {
         _PTM_Stat(pg,p,STAT) = 0;     // 1 idx_stat
         _PTM_Stat(pg,p,IDX2) = 0;     // 2 idx_2
         _PTM_Stat(pg,p,IDX3) = 0;     // 3 idx_3
         _PTM_Stat(pg,p,LAYR) = 0;     // 4 idx_layer
         _PTM_Stat(pg,p,FLAG) = 3;     // 5 flag
         _PTM_Stat(pg,p,PTID) = -9999; // 6 particle id
         _PTM_Vars(pg,p,MASS) = 0.0;   //
         _PTM_Vars(pg,p,DIAM) = 0.0;
         _PTM_Vars(pg,p,DENS)  = 0.0;
         _PTM_Vars(pg,p,VVEL)  = 0.0;
         // 2026-07 fix: HGHT is deliberately NOT zeroed here. A previous change set it
         // to 0.0 on the (false) premise that the -9999 sentinel leaked into height
         // comparisons - it does not: do_ptm_update, ptm_layershift, ptm_update_layerid,
         // ptm_removeparticles and the inflow redistribution all guard on
         // _PTM_Stat(pg,p,STAT) > 0, and unused slots have STAT == 0. Meanwhile
         // ptm_write_glm writes EVERY slot, so zeroing HGHT turned unused slots from an
         // obvious -9999 sentinel into "height 0.0" (i.e. the lake bed) - a plausible
         // value that silently pollutes any plot not filtering on particle_status.
      }
    }

    for (pg = 0; pg < num_particle_groups; pg++) {
        AED_REAL group_depth_min = ptm_group_real(init_depth_min_by_group, pg, init_depth_min);
        AED_REAL group_depth_max = ptm_group_real(init_depth_max_by_group, pg, init_depth_max);

        upper_height = Lake[surfLayer].Height - group_depth_min;
        if (upper_height > Lake[surfLayer].Height) {
            fprintf(stderr, "     ERROR: Particles cannot be initialized above the surface of the lake");
            exit(1);
        }
        lower_height = Lake[surfLayer].Height - group_depth_max;
        if (lower_height < 0 || lower_height > upper_height) {
            fprintf(stderr, "     ERROR: Particle group %d has invalid init_depth_min/init_depth_max\n", pg + 1);
            exit(1);
        }

        ptm_addparticles_group(pg, ptm_group_count(pg), max_particle_num,
                               upper_height, lower_height);
    }

    ptm_update_layerid();     // assign layers to active particles
    ptm_print_reseed_config();
    ptm_update_state_indices();   // ABM state indices, buoyancy parameters (no-op without the ABM)

}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *    This routine updates the particle positions, based on random motions    *
 *    and settling/migration                                                  *
 *                                                                            *
 ******************************************************************************/
/*----------------------------------------------------------------------------*/
/* Cache OASIM per-layer PAR and Kd for this hydrodynamic step. Both are constant over the step
 * (OASIM last ran inside wq_do_glm), so caching is exact rather than an approximation, and it
 * keeps the inner loop down to a single exp().                                               */
/* The layer a particle's height is in NOW, by the rule ptm_update_layerid applies (first layer
 * whose top is above the height; at or above the surface top: the surface layer). Read-only:
 * the writer uses it so particle_layer matches the H/NS written beside it. */
static int ptm_audit_layer(int pg, int p)
{
    int i;
    if ( _PTM_Vars(pg,p,HGHT) >= Lake[surfLayer].Height ) return surfLayer;
    for (i = botmLayer; i < NumLayers; i++)
        if ( _PTM_Vars(pg,p,HGHT) < Lake[i].Height ) return i;
    return surfLayer;
}

static void ptm_buoy_cache_light(void)
{
    int i, status;
    AED_REAL pc, kd;

    if (buoy_cache_n < NumLayers) {
        buoy_par_center = realloc(buoy_par_center, NumLayers * sizeof(AED_REAL));
        buoy_kd         = realloc(buoy_kd,         NumLayers * sizeof(AED_REAL));
        buoy_stat       = realloc(buoy_stat,       NumLayers * sizeof(int));
        if (buoy_par_center == NULL || buoy_kd == NULL || buoy_stat == NULL) {
            fprintf(stderr, "PTM BUOYANCY: out of memory caching light for %d layers\n", NumLayers);
            exit(1);
        }
        buoy_cache_n = NumLayers;
    }
    /* Detect dawn from the SURFACE light. B99's I_a is the dose accumulated since dawn, so it
     * resets once per day: a particle that received no light has no ballast to burn. Holding it
     * indefinitely (the previous version) let a particle that sank into the dark keep
     * de-ballasting at the last sunny day's rate for days, driving it onto rho_min. */
    {
        /* D3: use the plumbed I_dark (buoy_I_dark, the minimum over active groups) rather than a
         * literal, and apply HYSTERESIS so cloud flickering across the threshold at dawn or dusk
         * cannot re-trigger the reset several times in one day. Day starts above I_dark and does
         * not end until the light falls below half of it. */
        int was = buoy_is_day;
        AED_REAL sfc = Lake[surfLayer].Light * 0.45;
        if (sfc > buoy_I_dark)            buoy_is_day = 1;
        else if (sfc < 0.5 * buoy_I_dark) buoy_is_day = 0;
        buoy_dawn_reset = (buoy_is_day && !was) ? 1 : 0;

        /* Temporal correction. ip_par_scalar_depth was written during the PREVIOUS wq_do_glm, using
         * that step's surface irradiance; the sun has since moved. Scale by how much the surface
         * light changed, which leaves the depth correction below to handle only the particle's own
         * movement. Without this the residual error is concentrated at dawn and dusk, where an hour
         * changes the irradiance by more than the particle's motion does. */
        if (buoy_sfc_last > 1.0e-3 && sfc > 1.0e-3) {
            buoy_sfc_ratio = sfc / buoy_sfc_last;
            if (buoy_sfc_ratio < 0.01) buoy_sfc_ratio = 0.01;
            if (buoy_sfc_ratio > 100.0) buoy_sfc_ratio = 100.0;
        } else if (sfc <= 1.0e-3) {
            buoy_sfc_ratio = 0.0;               /* night now: no light regardless of history */
        } else {
            buoy_sfc_ratio = 1.0;               /* first light after darkness; anchor unusable */
        }
        buoy_sfc_last = sfc;
    }
    for (i = 0; i < NumLayers; i++) {
        aed_sample_oasim_par_kd(i, &pc, &kd, &status);
        if (status != 0 || pc <= 0.0 || kd < 0.0 || !isfinite(pc) || !isfinite(kd)) {
            buoy_par_center[i] = Lake[i].Light * 0.45;  /* same fallback ptm_update_particle_par uses */
            buoy_kd[i] = 0.0;
            buoy_stat[i] = 1;                       /* record it -- a silent fallback is the whole
                                                     * reason the light bug went undiagnosed */
        } else {
            buoy_par_center[i] = pc;
            buoy_kd[i] = kd;
            buoy_stat[i] = 0;
        }
    }
}

/* One Wallace & Hamilton buoyancy step for a single particle. Mirrors exactly the reference
 * implementation in analysis/test_buoyancy_wh.py (function step()), which reproduces B99 Eq 7
 * analytically and their Figs 2-3.
 *
 *   B_eq = c1*I/(KI+I) - c3              in light  (B99 Eq 4, Kromkamp & Walsby 1990)
 *   B_eq = -c2*I_a - c3                  in dark   (B99 Eq 11)
 *   B    = B_eq + (B-B_eq)*exp(-dt/tau)  if B < B_eq   -- ramping up, LAGGED
 *   B    = B_eq                          otherwise     -- light fell, NO lag
 *
 * The asymmetry IS the mechanism, worth a measured 2.6x (19.36 vs 7.46 kg/m3 in the offline
 * test). Branching on B < B_eq rather than "irradiance rose since the last step" matters: under
 * B99's own constant-light experiment irradiance rises only at t=0, so an I > I_prev test would
 * apply the lag for a single step and then release it, reproducing none of Fig. 2.
 *
 * Deterministic by construction -- it must never call rand(). The random walk, the bed Bernoulli
 * draws, seeding and redistribution all share one stream, so a single extra draw would shift
 * every later trajectory and destroy the buoy_model = 0 regression.                           */
static void ptm_buoyancy_step(int pg, int p, int lay, AED_REAL dt_sec)
{
    const AED_REAL *pr = buoy_prm[pg];
    AED_REAL dt_day = dt_sec / 86400.0;
    AED_REAL h, I, scal, eff, Ia, B, beq, s_act, r, Bg, rho, fcol, d_eff, pw, mu, rho_now;

    /* D2: bounds-check the layer before using it. A stale index is possible when
     * `lay` could predate a layer merge/split since the last ptm_update_layerid, and the lake has been re-layered by
     * mixing/inflow/outflow. buoy_par_center[] is only ever grown and only refilled below
     * NumLayers, so a stale index reads a previous, deeper lake's light with no warning.
     * ptm_update_particle_par guards this (glm_ptm.c ~1311); this routine did not. */
    if (lay < botmLayer || lay >= NumLayers) return;

    /* Irradiance at the particle's CURRENT height, from the cached layer solution. */
    h = _PTM_Vars(pg,p,HGHT);
    /* D2: clamp h into its layer before extrapolating. Unclamped, exp(kd*(h-MeanHeight)) with
     * kd ~ 1/m and a multi-metre offset reaches e^10 ~ 2e4, which does not produce a NaN -- it
     * flips the branch, so a particle in genuine darkness takes the LIGHT branch and ballasts at
     * up to c1. Mirrors the clamp ptm_update_particle_par already performs. */
    {
        AED_REAL top = Lake[lay].Height;
        AED_REAL bot = (lay > botmLayer) ? Lake[lay - 1].Height : 0.0;
        if (h > top) h = top;
        if (h < bot) h = bot;
    }
    /* MEASURED BUG, fixed here: using buoy_par_center[lay] as the magnitude gave a median
     * buoy_I/particle_par of 0.066 -- 15x too little light -- because aed_sample_oasim_par_kd() is
     * called from do_ptm_update, which runs BEFORE wq_do_glm solves OASIM for this step. The
     * fallback was only 0.3% of samples, so the data was valid but from the wrong time.
     * Anchor instead on the PAR AED wrote for this particle at the end of the previous step, and
     * use the cached kd only for the small correction from h0 to the current h. */
    {
        /* Lake[].Light is CURRENT-time and already depth-resolved (fix_radiation runs at
         * glm_model.c:1058, before do_ptm_update at :1077), so using it for the shape and the
         * timing removes BOTH extrapolations the previous version needed - the particle's
         * displacement over the step and the sun's motion over the hour. buoy_scale carries AED's
         * spectral/OASIM calibration, measured for this particle at the start of the step. The only
         * extrapolation left is within a layer, exactly as ptm_update_particle_par does. */
        I = Lake[lay].Light * 0.45
            * exp(buoy_kd[lay] * (h - Lake[lay].MeanHeight))
            * buoy_scale[pg * max_particle_num + p];
    }

    /* Carry the spectral / NPQ weighting AED applied at the last hourly call, so the sub-step
     * light differs from the biology's only by depth, never by basis. */
    /* D4: scal is a PAR in W/m2, O(1)-O(1000). A 1e-12 floor admitted ratios of ~1e14 that the
     * isfinite test below happily passes. Require a physically meaningful magnitude and bound the
     * ratio -- it is a spectral/NPQ weighting, so it belongs near 1. */
    scal = _PTM_Vars(pg,p,PTM_STATE_VAR(ptm_idx_par_scalar_depth));
    eff  = _PTM_Vars(pg,p,PTM_STATE_VAR(ptm_idx_par));
    if (scal > 1.0e-3 && eff > 0.0) {
        AED_REAL w = eff / scal;
        if (w < 0.1) w = 0.1;
        if (w > 10.0) w = 10.0;
        I *= w;
    }
    if (!isfinite(I) || I < 0.0) I = 0.0;
    if (ptm_idx_buoy_I > 0)                       /* record the irradiance ACTUALLY used */
        _PTM_Vars(pg,p,PTM_STATE_VAR(ptm_idx_buoy_I)) = I;
    if (ptm_idx_buoy_stat > 0)
        _PTM_Vars(pg,p,PTM_STATE_VAR(ptm_idx_buoy_stat)) = (AED_REAL) buoy_stat[lay];

    /* Light dose: B99 Eq 11 drives overnight de-ballasting with I_a, "the integrated value of
     * irradiance accumulated since dawn" -- the cell's MEMORY of the day's light, which it burns
     * off in the dark. The EMA is therefore updated ONLY while the particle is in light, and HELD
     * through darkness. An EMA that also decays at night (the first version of this code) forgets
     * the dose by dawn, collapsing the dark term to -c3 and, measured, shrinking the realised
     * density span instead of driving the diel cycle. */
    rho_now = _PTM_Vars(pg,p,DENS);
    /* D1: the dawn reset is NOT applied here. This routine runs sub_steps/ptm_buoy_substep = 30
     * times per hydrodynamic step, so resetting here zeroed Ia before each of the 30 EMA
     * increments -- a measured ~28x suppression of the dawn dose. It is now applied exactly once,
     * in do_ptm_update, before the sub-step loop. */
    Ia = _PTM_Vars(pg,p,PTM_STATE_VAR(ptm_idx_lightdose));
    if (I > pr[5])                            /* build only in light ... */
        Ia += (I - Ia) * (1.0 - exp(-dt_day / (pr[6] / 24.0)));
                                              /* ... and hold, unchanged, through the night */

    beq = (I <= pr[5]) ? (-pr[1] * Ia - pr[2])
                       : (pr[0] * I / (pr[3] + I) - pr[2]);
    /* Visser, Passarge & Mur (1997) Eq 3: the ballast pool also decays first-order, dC/dt = -c2*C,
     * c2 = 2.64e-7 /s = 0.0228 /day. In density terms -k*(rho - rho_min). Both source papers cite
     * this term. It is what makes the dark side WORK in a deep lake:
     *   - light-independent, so a particle that leaves the photic zone can still shed ballast
     *     (with only W&H's -c2*Ia - c3 it cannot, and measured, both groups sank to the bed and
     *     stayed: cell-weighted median height 6.85 m -> 0.02 m)
     *   - proportional to the ballast carried, so it vanishes as rho -> rho_min and cannot pin a
     *     cell on the floor (a flat c3 = 3.833 did exactly that: 24% and 43% of samples clamped)
     * pr[10] is k in /day. */
    beq -= pr[10] * (rho_now - pr[7]);

    if (!isfinite(Ia) || Ia < 0.0) Ia = 0.0;
    B = _PTM_Vars(pg,p,PTM_STATE_VAR(ptm_idx_buoy_rate));
    if (B < beq)   /* exponential integrator: exact for piecewise-constant beq and
                    * unconditionally stable, so no dt/tau combination can blow up */
        B = beq + (B - beq) * exp(-dt_day / (pr[4] / 1440.0));
    else
        B = beq;

    /* Same dormancy gate as aed_phyto_abm.F90: without it a dormant particle freezes into a
     * terminal density trap (documented there). s_act is published by AED, never re-derived. */
    s_act = _PTM_Vars(pg,p,PTM_STATE_VAR(ptm_idx_s_act));
    if (s_act < 0.0) s_act = 0.0; else if (s_act > 1.0) s_act = 1.0;
    r = pr[9];
    if (r < 0.0) r = 0.0; else if (r > 1.0) r = 1.0;   /* D6: clamp as s_act already is */
    Bg = (B > 0.0) ? (B * s_act) : (B * (r + (1.0 - r) * s_act));

    /* D5: validate BEFORE committing anything. The previous version returned here after Ia and B
     * had already been written, so the "leave state untouched" claim was false and a non-finite B
     * persisted into the next call. */
    rho = rho_now + Bg * dt_day;
    if (!isfinite(rho) || !isfinite(B)) return;
    if (pr[7] <= pr[8]) {                     /* D6: inverted bounds would pin every particle */
        if (rho > pr[8]) rho = pr[8];
        if (rho < pr[7]) rho = pr[7];
    }
    _PTM_Vars(pg,p,PTM_STATE_VAR(ptm_idx_lightdose)) = Ia;
    _PTM_Vars(pg,p,PTM_STATE_VAR(ptm_idx_buoy_rate)) = B;
    _PTM_Vars(pg,p,DENS) = rho;

    /* Stokes, identical in form to aed_phyto_abm.F90 (_MOB_STOKES_). DIAM is already metres
     * (AED overwrites it with ESD*1e-6). Positive VVEL = upward. */
    fcol = _PTM_Vars(pg,p,PTM_STATE_VAR(ptm_idx_fcol));
    if (!(fcol >= 1.0)) fcol = 1.0;
    d_eff = _PTM_Vars(pg,p,DIAM) * fcol;
    pw = Lake[lay].Density;
    mu = aed_water_viscosity_c(Lake[lay].Temp);
    if (mu > 1.0e-9 && isfinite(pw))
        _PTM_Vars(pg,p,VVEL) = -9.807 * d_eff * d_eff * (rho - pw) / (18.0 * mu);
}

/*----------------------------------------------------------------------------*/


void do_ptm_update()
{
//LOCALS
    int p, tt, ij1, ij2, sub_steps, pg, layr, lbe, v;
    AED_REAL K_z, K_local, K_above, K_below, K_prime_z, dz, rand_draw, settling_efficiency_substep, ptm_dt_rem;
    float rand_float, prob, prev_height, x1, x2, y1, y2, a1, a2;

/*----------------------------------------------------------------------------*/
//BEGIN
    if (ptm_time_step <= 0.0) {
        fprintf(stderr, "     ERROR: ptm_time_step must be > 0 seconds (namelist &particles)\n");
        exit(1);
    }
    ptm_dt_rem = fmod(dt, ptm_time_step);
    if (ptm_dt_rem > 1.0e-6 && (ptm_time_step - ptm_dt_rem) > 1.0e-6) {
        fprintf(stderr, "     ERROR: host timestep dt=%.6f s is not evenly divisible "
                         "by ptm_time_step=%.6f s (namelist &particles)\n", dt, ptm_time_step);
        exit(1);
    }

    sub_steps = (int)(dt / ptm_time_step + 0.5);   // round-to-nearest; exact per the check above
    settling_efficiency_substep = 1.0 - exp(-settling_efficiency * (ptm_time_step / 86400.0));

    // LAYR was set at the END of the previous call and the column has been restructured since
    // (surface fluxes, check_layer_thickness, check_layer_stability, the mixer, the daily block).
    // Refresh first, so every read below starts from the layer each particle is actually in.
    ptm_update_layerid();

    if (buoy_any) {
        ptm_buoy_cache_light();
        /* Capture the height at which AED last evaluated this particle's PAR. ip_par_scalar_depth
         * was written at the end of the previous step by ptm_update_particle_par, i.e. at exactly
         * this height, so it is the correct anchor for the within-step depth correction. */
        {
            int need = num_particle_groups * max_particle_num;
            if (buoy_h0_n < need) {
                buoy_h0 = realloc(buoy_h0, need * sizeof(AED_REAL));
                buoy_scale = realloc(buoy_scale, need * sizeof(AED_REAL));
                if (buoy_h0 == NULL || buoy_scale == NULL) {
                    fprintf(stderr, "PTM BUOYANCY: out of memory for h0 buffer\n");
                    exit(1);
                }
                buoy_h0_n = need;
            }
            for (pg = 0; pg < num_particle_groups; pg++)
                for (p = 0; p < max_particle_num; p++) {
                    int idx0 = pg * max_particle_num + p;
                    int l0 = _PTM_Stat(pg,p,LAYR);
                    AED_REAL hh0 = _PTM_Vars(pg,p,HGHT);
                    AED_REAL aed_par, lake_par;
                    buoy_h0[idx0] = hh0;
                    buoy_scale[idx0] = 1.0;
                    if (buoy_model[pg] != 1) continue;
                    if (l0 < botmLayer || l0 >= NumLayers) continue;
                    /* Calibrate the Lake-based estimate against what AED actually gave this
                     * particle, at the same position, so its spectral/OASIM basis is preserved
                     * while the TIMING and the depth shape come from the current-step field. */
                    aed_par  = _PTM_Vars(pg,p,PTM_STATE_VAR(ptm_idx_par_scalar_depth));
                    lake_par = Lake[l0].Light * 0.45
                             * exp(buoy_kd[l0] * (hh0 - Lake[l0].MeanHeight));
                    if (aed_par > 1.0e-3 && lake_par > 1.0e-3) {
                        buoy_scale[idx0] = aed_par / lake_par;
                        if (buoy_scale[idx0] < 0.05) buoy_scale[idx0] = 0.05;
                        if (buoy_scale[idx0] > 20.0) buoy_scale[idx0] = 20.0;
                    }
                }
        }
        /* The dawn reset happens exactly ONCE per step, here, not inside ptm_buoyancy_step
         * (which runs sub_steps/ptm_buoy_substep times per step and would zero Ia before each
         * of the EMA increments - a measured ~28x suppression of the dawn dose). */
        if (buoy_dawn_reset && ptm_idx_lightdose > 0) {
            for (pg = 0; pg < num_particle_groups; pg++) {
                if (buoy_model[pg] != 1) continue;
                for (p = 0; p < max_particle_num; p++)
                    if (_PTM_Stat(pg,p,STAT) > 0)
                        _PTM_Vars(pg,p,PTM_STATE_VAR(ptm_idx_lightdose)) = 0.0;
            }
        }
    }

    for (tt = 0; tt < sub_steps; tt++) {
      for (pg = 0; pg < num_particle_groups; pg++) {
        for (p = 0; p < max_particle_num; p++) {
          if (_PTM_Stat(pg,p,STAT)>0) {

            layr = _PTM_Stat(pg,p,LAYR);
            // ptm_update_layerid() keeps LAYR current every sub-step; the clamp only guards a
            // slot whose layer index predates a layer merge/split since that call.
            if (layr < botmLayer) layr = botmLayer;
            if (layr > surfLayer) layr = surfLayer;

            /* Wallace & Hamilton buoyancy, at ptm_buoy_substep resolution. Runs BEFORE the move
             * so the density and VVEL used to move the particle this sub-step are the ones just
             * computed. nstep handles a final partial interval exactly, so the integrated time
             * always equals the elapsed time. */
            if (buoy_any && buoy_model[pg] == 1 && (tt % ptm_buoy_substep) == 0) {
                int nleft = sub_steps - tt;
                int nstep = (nleft < ptm_buoy_substep) ? nleft : ptm_buoy_substep;
                ptm_buoyancy_step(pg, p, layr, nstep * ptm_time_step);
            }

            K_z = Lake[layr].Epsilon;
            if (K_z < ptm_diffusivity) K_z = ptm_diffusivity;

            if(layr == surfLayer){
                // No layer above surfLayer to sample (Lake[surfLayer+1] would be out of
                // bounds), so there is no gradient to compute - zero is the correct
                // K_prime_z here, not a special case. The `continue` this replaced,
                // though, skipped past move_particle() below entirely: every surface-layer
                // particle received NO diffusion, NO settling, and NEVER reached the
                // BED/SCUM boundary checks further down, for the whole run. That silently
                // froze the surface population and excluded it from settling regardless
                // of species, physiology, or vertical velocity.
                K_prime_z = 0;
                K_local = K_z;
            } else {
                K_above = Lake[layr+1].Epsilon;
                if (K_above < ptm_diffusivity) K_above = ptm_diffusivity;
                lbe = (layr-1 >= botmLayer) ? layr-1 : layr;
                K_below = Lake[lbe].Epsilon;
                if (K_below < ptm_diffusivity) K_below = ptm_diffusivity;
                dz = Lake[layr+1].MeanHeight - Lake[lbe].MeanHeight;   // central spacing
                if (fabs(dz) < 1.0e-6) dz = (dz < 0.0) ? -1.0e-6 : 1.0e-6;
                K_prime_z = (K_above - K_below) / dz;                 // signed dK/dz (m/s)
                K_local = K_z + K_prime_z * (_PTM_Vars(pg,p,HGHT) - Lake[layr].MeanHeight);
                if (K_local < ptm_diffusivity) K_local = ptm_diffusivity;
            }

            // Capture current height of particle to calculate probability of settling below
            prev_height = _PTM_Vars(pg,p,HGHT);

            // Update particle position based on diffusivity and vert velocity
            _PTM_Stat(pg,p,FLAG)= WATER;

            // Every draw goes through glm_rand() so the stream can be replayed on restart.
            rand_draw = -1.0 + 2.0*((AED_REAL)glm_rand() / (AED_REAL)RAND_MAX);
            _PTM_Vars(pg,p,HGHT) = move_particle(ptm_time_step,_PTM_Vars(pg,p,HGHT), K_local, K_prime_z, _PTM_Vars(pg,p,VVEL), rand_draw);

            // An overshoot below the bed: mirrored back into the water column when the lower
            // boundary is reflective (clamping every overshoot to one coordinate builds a
            // delta function at the boundary), otherwise clamped to the bed.
            if(_PTM_Vars(pg,p,HGHT) < 0.0){
                if (lower_boundary_cond == BC_REFLECTIVE)
                    _PTM_Vars(pg,p,HGHT) = -_PTM_Vars(pg,p,HGHT);
                else
                    _PTM_Vars(pg,p,HGHT) = 0.0;
            }

            if (prev_height > _PTM_Vars(pg,p,HGHT)){

                // Get area at previous particle height
                x1 = prev_height * 10.0;
                y1 = x1 - (int)(x1 / 1.0) * 1.0;
                ij1 = (int)(x1 - y1) - 1;
                if(ij1 >= Nmorph){
                    y1 = y1 + (float)(ij1 - Nmorph);
                    ij1 = Nmorph - 1;
                } else if (ij1 < 0 ) ij1 = 0;
                a1 = MphLevelArea[ij1] + y1 * dMphLevelArea[ij1];

                // Get area at depth of current (already boundary-corrected) particle height
                x2 = _PTM_Vars(pg,p,HGHT) * 10.0;
                y2 = x2 - (int)(x2 / 1.0) * 1.0;
                ij2 = (int)(x2 - y2) - 1;
                if(ij2 >= Nmorph){
                    y2 = y2 + (float)(ij2 - Nmorph);
                    ij2 = Nmorph - 1;
                } else if (ij2 < 0 ) ij2 = 0;
                a2 = MphLevelArea[ij2] + y2 * dMphLevelArea[ij2];

                // Calculate proportional difference between two areas.
                prob = (a1 - a2) / a1;

                // Bernoulli draw to determine if particle should be assigned as BED;
                // if triggered, a second draw tests settling_efficiency_substep (this
                // call's per-substep equivalent of the per-day settling_efficiency rate).
                rand_float = (float)glm_rand() / (float)RAND_MAX;
                if(rand_float < prob){
                    rand_float = (float)glm_rand() / (float)RAND_MAX;
                    if(rand_float < settling_efficiency_substep){
                        _PTM_Stat(pg,p,FLAG) = BED;
                        if(sed_deactivation){
                            _PTM_Stat(pg,p,STAT) = 0;
                            // Zero the whole ABM state block (age, birth day, cumulative
                            // C/N/P/Chl, division count, etc.), not just STAT. A deactivated
                            // slot is exactly what ptm_next_free_slot() looks for, so leaving
                            // its old life history behind means whatever new particle reuses
                            // this slot inherits it, and can also trip the biology's
                            // "already seeded" guard (aed_phyto_abm re-seeds an empty state).
                            for (v = 0; v < Num_PTM_Vars; v++)
                                _PTM_Vars(pg,p,PTM_ENV_NVARS + v) = 0.0;
                            // Deactivated: stop visiting this particle for the remainder of
                            // this call's substeps (the STAT>0 guard above skips it from here on).
                            continue;
                        }
                    }else{
                        if (lower_boundary_cond == BC_REFLECTIVE) {
                            // Bounce back to prev_height - undo this substep's downward
                            // step, the same distance the particle fell, instead of
                            // resting at the height (a2) that triggered settling.
                            _PTM_Vars(pg,p,HGHT) = prev_height;
                        }
                    }
                }
            }

            // Apply the upper (lake surface) boundary condition.
            //   upper_boundary_cond = 1 (clamp):     height set to the surface height (SCUM).
            //   upper_boundary_cond = 2 (reflect):   height mirrored about the surface height;
            //                                        SCUM only if still at/above the surface.
            if(_PTM_Vars(pg,p,HGHT)>Lake[surfLayer].Height){
                if (upper_boundary_cond == BC_REFLECTIVE) {
                    // Mirror about the surface and remain in the water column.
                    _PTM_Vars(pg,p,HGHT) = 2.0*Lake[surfLayer].Height - _PTM_Vars(pg,p,HGHT);
                    if(_PTM_Vars(pg,p,HGHT) < 0.0) _PTM_Vars(pg,p,HGHT) = 0.0;
                    if(_PTM_Vars(pg,p,HGHT) >= Lake[surfLayer].Height){   // overshoot larger than the column
                        _PTM_Stat(pg,p,FLAG)= SCUM;
                        _PTM_Vars(pg,p,HGHT)=Lake[surfLayer].Height;
                    }
                } else { /* BC_CLAMP (default) */
                    // Determine if particle should be assigned as SCUM
                    _PTM_Stat(pg,p,FLAG)= SCUM;                  // Maybe forming a scum layer
                    _PTM_Vars(pg,p,HGHT)=Lake[surfLayer].Height;
                }
            }
          }
        }
      }
      ptm_update_layerid();
    }

    ptm_apply_reseeding();    // optional diagnostic-probe maintenance (off by default)
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *    This routine redistributes any particles within the provided depth      *
 *    range, used to capture the effect of layer mixing                       *
 *                                                                            *
 ******************************************************************************/
void ptm_redistribute(AED_REAL upper_height, AED_REAL lower_height)
{
//LOCALS
    int p;

    int pg;

/*----------------------------------------------------------------------------*/
//BEGIN
    // Check for active particles in the height range
    for (pg = 0; pg < num_particle_groups; pg++) {
      for (p = 0; p < max_particle_num; p++) {
        if (_PTM_Stat(pg,p,STAT)>0) {
            if (_PTM_Vars(pg,p,HGHT)>=lower_height && _PTM_Vars(pg,p,HGHT)<=upper_height ) {
                // Particle is in the mixing zone, so re-position, weighted by
                // cross-sectional area within [lower_height, upper_height] rather
                // than uniform: a thick mixing zone isn't a box, and a height with
                // more lake area should hold more of the well-mixed particles, not
                // the same density per unit height as a narrower part of the zone.
                _PTM_Vars(pg,p,HGHT) = draw_height_in_range(lower_height, upper_height);
            }
        }
      }
    }
    ptm_update_layerid();     // assign layers to active particles
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *    This routine adds new particles within the provided depth               *
 *    range                                                                   *                                                                        *
 *                                                                            *
 ******************************************************************************/
static void ptm_addparticles_group(int pg, int new_particles, int max_particle_num,
                                   AED_REAL upper_height, AED_REAL lower_height)
{
    ptm_addparticles_group_ex(pg, new_particles, max_particle_num,
                              upper_height, lower_height, 0);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


static void ptm_addparticles_group_ex(int pg, int new_particles, int max_particle_num,
                                      AED_REAL upper_height, AED_REAL lower_height,
                                      int force_reseed_id)
{
//LOCALS
    int p, n, v;
    static int next_reseed_ptid = 0;   // monotonic id source for reseeded slots

    AED_REAL height_range;

/*----------------------------------------------------------------------------*/
//BEGIN
    // Get vertical range in the water column that mixed
    height_range = upper_height - lower_height;

    // For each new particle, initialise their properties and height.
    // ptm_next_free_slot() carries a rolling cursor across calls (see its header
    // comment above), so allocating new_particles of them costs about one pass
    // over the pool in total rather than one pass per particle.
    for (n = 0; n < new_particles; n++) {
        p = ptm_next_free_slot(pg, max_particle_num);
        if (p < 0) {
            printf("ptm_addparticles(): WARNING group %d has no available particle slots; increase max_particle_num\n", pg + 1);
            break;
        }
        _PTM_Stat(pg,p,STAT) = 1;
        _PTM_Stat(pg,p,FLAG) = 0;
        _PTM_Stat(pg,p,GRP)  = pg + 1;   // 1-based, as aed_phyto_abm reads it
        _PTM_Vars(pg,p,MASS) = 1.0;
        _PTM_Vars(pg,p,DIAM) = get_particle_diameter(ptm_group_real(particle_diameter_by_group, pg, particle_diameter));
        _PTM_Vars(pg,p,DENS) = get_particle_density(ptm_group_real(particle_density_by_group, pg, particle_density));
        _PTM_Vars(pg,p,VVEL) = get_settling_velocity(ptm_group_real(settling_velocity_by_group, pg, settling_velocity));
        // The ABM state block starts EMPTY: aed_phyto_abm seeds an active particle whose
        // C/N/P/Chl/num are zero from its own parameter table on the next biology step, so
        // the host never duplicates those values.
        for (v = 0; v < Num_PTM_Vars; v++)
            _PTM_Vars(pg,p,PTM_ENV_NVARS + v) = 0.0;

        // Assign PTID: ids run 1..N over all groups (slot p of group pg is
        // pg*max_particle_num + p + 1), then continue monotonically for reseeded slots,
        // so no two particles of the run ever share an id.
        if(force_reseed_id){
           if (next_reseed_ptid < num_particle_groups*max_particle_num)
               next_reseed_ptid = num_particle_groups*max_particle_num;
           _PTM_Stat(pg,p,PTID) = ++next_reseed_ptid;
        } else if(_PTM_Stat(pg,p,PTID) < 0){
           _PTM_Stat(pg,p,PTID) = pg*max_particle_num + p + 1;
        } else {
           if (next_reseed_ptid < num_particle_groups*max_particle_num)
               next_reseed_ptid = num_particle_groups*max_particle_num;
           _PTM_Stat(pg,p,PTID) = ++next_reseed_ptid;
        }

        // Assign particles initial height: the whole rand() value scaled to [0,1)
        // (mean 0.5, RAND_MAX+1 levels) rather than rand() % 100.
        {
            double random_double = (double)glm_rand() / ((double)RAND_MAX + 1.0);
            random_double = random_double * height_range;           // scale unit random to requested range
            _PTM_Vars(pg,p,HGHT) = lower_height + random_double;     // set particle height
        }
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


static void ptm_apply_reseeding(void)
{
    int pg, p, active, min_active, target_active, need;
    AED_REAL group_depth_min, group_depth_max, upper_height, lower_height;

    if (!particle_reseed_enabled) return;

    for (pg = 0; pg < num_particle_groups; pg++) {
        active = 0;
        for (p = 0; p < max_particle_num; p++) {
            if (_PTM_Stat(pg,p,STAT) > 0) active++;
        }

        min_active = ptm_reseed_count(particle_reseed_min_by_group, pg, particle_reseed_min,
                                      "particle_reseed_min");
        target_active = ptm_reseed_count(particle_reseed_target_by_group, pg, particle_reseed_target,
                                         "particle_reseed_target");
        if (target_active <= 0) target_active = ptm_group_count(pg);
        if (min_active <= 0) continue;
        if (target_active <= min_active) target_active = min_active;
        if (target_active > max_particle_num) target_active = max_particle_num;
        if (active >= min_active) continue;

        group_depth_min = ptm_group_real(init_depth_min_by_group, pg, init_depth_min);
        group_depth_max = ptm_group_real(init_depth_max_by_group, pg, init_depth_max);
        upper_height = Lake[surfLayer].Height - group_depth_min;
        lower_height = Lake[surfLayer].Height - group_depth_max;
        if (upper_height > Lake[surfLayer].Height) upper_height = Lake[surfLayer].Height;
        if (lower_height < 0.0) lower_height = 0.0;
        if (lower_height > upper_height) continue;

        need = target_active - active;
        if (need > 0) {
            ptm_log_reseed_event(pg, active, min_active, target_active, need);
            ptm_addparticles_group_ex(pg, need, max_particle_num,
                                      upper_height, lower_height, 1);
        }
    }

    ptm_update_layerid();
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


void ptm_addparticles(int new_particles, int max_particle_num, AED_REAL upper_height,
                      AED_REAL lower_height)
{
    int pg;
    for (pg = 0; pg < num_particle_groups; pg++)
        ptm_addparticles_group(pg, new_particles, max_particle_num, upper_height, lower_height);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *    This routine removes new particles from a specified layer, based        *
 *    on the proportion of layer volume that is removed through outflow       *
 *                                                                            *
 *                                                                            *
 ******************************************************************************/
void ptm_removeparticles(int layer_id, AED_REAL delta_vol, AED_REAL layer_vol, int max_particle_num)
{
//LOCALS
    AED_REAL layer_prop, rand_float;
    int p, pg, v;

/*----------------------------------------------------------------------------*/
//BEGIN
    // For each particle, draw from Bernoulli distribution to see whether removed from layer
    layer_prop = delta_vol / layer_vol;
    for (pg = 0; pg < num_particle_groups; pg++) {
      for (p = 0; p < max_particle_num; p++) {
        if(_PTM_Stat(pg,p,STAT) == 1 && _PTM_Stat(pg,p,LAYR) == layer_id){
            rand_float = ((AED_REAL)glm_rand())/RAND_MAX;
            if(rand_float <= layer_prop){
                // If particle leaves through outflow, reset completely.
                //
                // STAT==0 && FLAG==3 (EXIT) is exactly the free-slot test
                // ptm_addparticles() looks for when seeding a new particle, so this
                // slot WILL be reused. The comment here always said "reset completely",
                // but the resets themselves were commented out below - so the departed
                // particle's whole ABM life history (age, birth day, cumulative C fixed,
                // division count, C/N/P/Chl state) survived into whatever new particle
                // got the slot next, and a leftover C>0 && num>0 could make the biology's
                // own seeding guard skip re-initialising the new occupant entirely.
                // Zero the ABM block by loop (not by naming individual variables, which
                // is exactly the mistake the commented-out lines below made - listing
                // MASS/DIAM/DENS/VVEL while leaving every biology variable untouched).
                // STAT=0 / FLAG=EXIT(3) is all that is needed to return the slot to
                // circulation - ptm_next_free_slot() scans for exactly that.
                _PTM_Stat(pg,p,STAT) = 0;
                _PTM_Stat(pg,p,FLAG) = 3;
                for (v = 0; v < Num_PTM_Vars; v++)
                    _PTM_Vars(pg,p,PTM_ENV_NVARS + v) = 0.0;
            }
        }
      }
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *    This routine raises or lowers particle vertical positions due to        *
 *    changes in layering associated with inflows/outflows                    *
 *                                                                            *
 *                                                                            *
 ******************************************************************************/
void ptm_layershift(AED_REAL shift_height, AED_REAL shift_amount)
{
//LOCALS
    int p,pg;

    AED_REAL upper_height;
    AED_REAL lower_height;

/*----------------------------------------------------------------------------*/
//BEGIN

    /* 2026-08-29 `shift_amount` was used BOTH as the upper bound of the affected range and as the
     * displacement applied. The sole caller is glm_flow.c:
     *     ptm_layershift(0.0, Lake[surfLayer].Height - height_start);
     * whose own comment reads "ASSUMING SHIFT IS ALL LAYERS" - i.e. the intent is that every
     * particle in the column moves with the rising free surface. As written only particles within
     * `shift_amount` metres of the bottom were lifted (a few cm for a typical inflow), so the
     * ensemble drifted steadily DOWNWARD relative to the surface over a run; and when the level
     * FELL, shift_amount < 0 made upper_height < lower_height so no particle moved at all.
     * The affected range is the whole water column; only the displacement is shift_amount.
     * Inert for these enclosures - num_inflows = 0, so this is never called. */
    lower_height = shift_height;
    upper_height = Lake[surfLayer].Height;

    // Check for active particles in the impacted height range
    for (pg = 0; pg < num_particle_groups; pg++) {
      for (p = 0; p < max_particle_num; p++) {
        if (_PTM_Stat(pg,p,STAT)>0) {
            if (_PTM_Vars(pg,p,HGHT)>lower_height && _PTM_Vars(pg,p,HGHT)<upper_height ) {
                // Particle is in the impacted zone, so re-position (lift or drop)
                _PTM_Vars(pg,p,HGHT) = _PTM_Vars(pg,p,HGHT) + shift_amount;   // adjust particle height
            }
        }
      }
    }

    ptm_update_layerid();     // assign layers to active particles
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *    B2h 2026-09-16: a water removal lowered the free surface - the          *
 *    particles it did not export follow it down.                             *
 *                                                                            *
 *    do_single_outflow accounts the export exactly once, per layer, through   *
 *    ptm_removeparticles, decrements the layer volumes and then lowers every  *
 *    layer height with resize_internals(2, botmLayer). It moved no surviving  *
 *    particle, so a particle near the old surface kept its absolute height    *
 *    while the surface dropped beneath it, and the next record showed it      *
 *    above the water (H2 BOUNDS: base 1, outlet_hi 58, outlet_lo 57).         *
 *                                                                            *
 *    The shift is DEPTH-PRESERVING, not a clamp: a particle that sat at the   *
 *    surface stays at the surface, one 1 mm below stays 1 mm below. Clamping  *
 *    every stranded particle to one coordinate would build a delta function   *
 *    at the boundary - the defect the 2026-08-29 bed comment above rejects    *
 *    on measured grounds, and the reason the bed and surface both reflect.    *
 *                                                                            *
 *    Only HGHT changes. Nothing is removed (the export already happened) and  *
 *    no STAT, FLAG, num, c or ABM slot is touched, so particle number and     *
 *    represented biomass are conserved by construction and the light and      *
 *    ballast HISTORY (cldose, lightdose, buoy_rate, fcol, s_act, age, birth,  *
 *    ndiv, the C/N/P/Chl quotas) is carried unchanged - a geometric boundary  *
 *    correction must not rewrite a cell's physiological memory.               *
 *                                                                            *
 *    No rand() call, deliberately: the random walk, the bed Bernoulli draws,  *
 *    seeding and redistribution share one stream whose call count is restart  *
 *    state (ptm_rng_calls), so a single extra draw would shift every later    *
 *    trajectory. A deterministic shift needs none.                            *
 *                                                                            *
 *    Interior particles are NOT remapped into volume coordinates. Only those  *
 *    stranded by the boundary are corrected; a volume-coordinate remap is a   *
 *    broader transport-model change and is deliberately not attempted here.   *
 *                                                                            *
 ******************************************************************************/
void ptm_follow_surface_drop(AED_REAL old_height)
{
//LOCALS
    int p, pg, moved = 0, floored = 0, above = 0;   /* B2j 2026-09-16: `above` counts particles
                                                    * already out of bounds BEFORE the water left */
    AED_REAL new_height, drop, h;   /* B2k 2026-09-16: `worst` removed - B2j deleted its only
                                     * reader, leaving a dead store. See B2j comment below. */
    static int warned = 0;   /* warn once per run, as aed_phyto_abm's split routine does */

/*----------------------------------------------------------------------------*/
//BEGIN
    if ( !ptm_sw ) return;
    new_height = Lake[surfLayer].Height;
    drop = old_height - new_height;
    /* a rise, no change, or a non-finite height: nothing for this routine to do. Written as a positive test
     * so a NaN drop falls through it rather than being compared into a branch. */
    if ( !(drop > 0.0) ) return;

    for (pg = 0; pg < num_particle_groups; pg++) {
      for (p = 0; p < max_particle_num; p++) {
        if (_PTM_Stat(pg,p,STAT) > 0 && _PTM_Vars(pg,p,HGHT) > new_height) {
            /* B2j: above the surface BEFORE any water left means the PREVIOUS step failed to
             * re-anchor it - the transport defect this routine's warning exists to surface. */
            if (_PTM_Vars(pg,p,HGHT) > old_height) above++;
            h = _PTM_Vars(pg,p,HGHT) - drop;
            if (h < 0.0) { h = 0.0; floored++; }   /* the drop exceeded the particle's own height */
            _PTM_Vars(pg,p,HGHT) = h;
            moved++;
        }
      }
    }

    if (moved) {
        /* the layer index must be fresh before the next read: every group is on Stokes settling, so a
         * particle's velocity comes from its own density against the AMBIENT density of the layer it is in */
        ptm_update_layerid();
        /* B2j 2026-09-16: report the ANOMALIES, not the routine case. The previous condition compared
         * `worst` against the surface layer's thickness, but worst == drop by construction (h = HGHT -
         * drop), so it only asked whether the withdrawal exceeded one layer - which is ordinary: the H2
         * outlet removes 8.64 m3/day, the level falls 0.1334 m and the surface layer is 0.1234 m thick,
         * so it tripped at 1.10x on a normal day with zero particles floored. A warning that fires on
         * normal operation trains the reader to ignore it, which is worse than silence. The two states
         * that ARE anomalous: a particle already above the surface before the water left (the previous
         * step failed to re-anchor it), or a drop exceeding a particle's own height. */
        if (!warned && (above > 0 || floored > 0)) {
            fprintf(stderr, "PTM SURFACE DROP: %d particle(s) re-anchored after the surface fell %.6g m; "
                            "%d were ALREADY above the surface before the removal and %d were floored at "
                            "the bed. Reported once per run.\n", moved, drop, above, floored);
            warned = 1;
        }
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *    This routine sets the layer id for each particle based on its height    *
 *                                                                            *
 ******************************************************************************/
void ptm_update_layerid()
{
//LOCALS
    int p, pg, layr;
    AED_REAL hght;

/*----------------------------------------------------------------------------*/
//BEGIN
    for (pg = 0; pg < num_particle_groups; pg++) {
      for (p = 0; p < max_particle_num; p++) {
        if (_PTM_Stat(pg,p,STAT)>0) {

            hght = _PTM_Vars(pg,p,HGHT);
            layr = _PTM_Stat(pg,p,LAYR);

            // Hunt from the particle's own previous layer instead of rescanning from
            // botmLayer every call: called every substep now, and a particle's height
            // moves by only one move_particle() step (bounded diffusion + settling
            // velocity * dt_secs) between calls, so it typically hasn't left its
            // previous layer at all. Clamp first in case NumLayers changed (layer
            // merge/split via check_layer_thickness()/check_layer_stability()/
            // do_deep_mixing()) since LAYR was last set - the walk below still
            // converges to the correct layer either way, just in more steps.
            if (layr < botmLayer)      layr = botmLayer;
            else if (layr > surfLayer) layr = surfLayer;

            // Particle sank below its bracket: walk down.
            while (layr > botmLayer && hght < Lake[layr-1].Height) layr--;
            // Particle rose above its bracket: walk up. surfLayer absorbs the clamp
            // case (do_ptm_update's upper_boundary_cond = BC_CLAMP sets HGHT exactly
            // to Lake[surfLayer].Height, which is never strictly less than itself) -
            // any particle not caught below belongs, by construction, in surfLayer.
            while (layr < surfLayer && hght >= Lake[layr].Height) layr++;

            _PTM_Stat(pg,p,LAYR) = layr;
            /* IDX3 is the cell the Fortran side bins this particle into, and it
             * is consumed 1-BASED there (aed_ptm.F90 guards on cell >= 1), while
             * C layer indices are 0-based (botmLayer 0 vs 1 in glm.h). Passing a
             * raw layr shifted every particle one layer down and made bottom-layer
             * particles (layr = 0) fail the >= 1 guard entirely, so they were never
             * binned and never had their physiology or environment updated.
             * LAYR stays 0-based - it indexes Lake[] on this side.            */
            _PTM_Stat(pg,p,IDX3) = layr + 1;
            _PTM_Stat(pg,p,IDX2) = 1;
        }
      }
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


static int ptm_total_particles(void)
{
    return num_particle_groups * max_particle_num;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


static int ptm_group_count(int grp)
{
    if (init_particle_num_by_group != NULL) return init_particle_num_by_group[grp];
    return init_particle_num;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


static int ptm_group_int(int *values, int grp, int fallback)
{
    if (values != NULL) return values[grp];
    return fallback;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


static int ptm_reseed_count(int *values, int grp, int fallback, const char *label)
{
    int raw = ptm_group_int(values, grp, fallback);
    int idx = (grp >= 0 && grp < 128) ? grp : 127;
    int *warned = (strcmp(label, "particle_reseed_min") == 0)
                  ? ptm_reseed_min_warned
                  : ptm_reseed_target_warned;

    if (raw < 0) {
        if (!warned[idx]) {
            fprintf(stderr, "     WARNING: %s for particle group %d is %d; using 0\n",
                    label, grp + 1, raw);
            warned[idx] = 1;
        }
        return 0;
    }
    if (raw > max_particle_num) {
        if (!warned[idx]) {
            fprintf(stderr, "     WARNING: %s for particle group %d is %d; capping to max_particle_num=%d\n",
                    label, grp + 1, raw, max_particle_num);
            warned[idx] = 1;
        }
        return max_particle_num;
    }
    return raw;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


static void ptm_print_reseed_config(void)
{
    int pg, min_active, target_active;

    if (!particle_reseed_enabled) return;

    fprintf(stderr, "     PTM diagnostic reseeding enabled: ");
    for (pg = 0; pg < num_particle_groups; pg++) {
        min_active = ptm_reseed_count(particle_reseed_min_by_group, pg, particle_reseed_min,
                                      "particle_reseed_min");
        target_active = ptm_reseed_count(particle_reseed_target_by_group, pg, particle_reseed_target,
                                         "particle_reseed_target");
        if (target_active <= 0) target_active = ptm_group_count(pg);
        if (target_active <= min_active) target_active = min_active;
        if (target_active > max_particle_num) target_active = max_particle_num;
        fprintf(stderr, "%sgroup %d min=%d target=%d",
                (pg == 0) ? "" : "; ", pg + 1, min_active, target_active);
    }
    fprintf(stderr, "\n");
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


static void ptm_log_reseed_event(int pg, int active, int min_active, int target_active, int need)
{
    int idx = pg;

    if (idx < 0 || idx >= 128) idx = 127;
    ptm_reseed_log_count[idx]++;

    if (ptm_reseed_log_count[idx] <= 5 || (ptm_reseed_log_count[idx] % 250) == 0) {
        if (ptm_reseed_log_suppressed[idx] > 0) {
            printf("ptm_reseed(): group %d suppressed %d repeated diagnostic-probe top-up messages\n",
                   pg + 1, ptm_reseed_log_suppressed[idx]);
            ptm_reseed_log_suppressed[idx] = 0;
        }
        printf("ptm_reseed(): group %d active %d below min %d; target %d; adding %d diagnostic probes\n",
               pg + 1, active, min_active, target_active, need);
    } else {
        ptm_reseed_log_suppressed[idx]++;
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


static AED_REAL ptm_group_real(AED_REAL *values, int grp, AED_REAL fallback)
{
    if (values != NULL) return values[grp];
    return fallback;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


static int ptm_vertical_mode(AED_REAL vvel_m_per_s)
{
    /* Convention (see move_particle(): updated_height += vvel*del_t, and
     * depth = surfHeight - height): POSITIVE vvel raises the particle, i.e.
     * upward/shallower. These labels were inverted until 2026-07, so buoyant
     * particles were reported as sinking and vice versa. */
    const AED_REAL eps = 1.0e-12;
    if (vvel_m_per_s > eps) return 1;  /* upward / buoyant  */
    if (vvel_m_per_s < -eps) return 3; /* downward / sinking */
    return 2;                          /* neutral */
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


static AED_REAL ptm_light_preference_nm(int pigment_id)
{
    switch (pigment_id) {
    case 1: return 620.0; /* phycocyanin */
    case 2: return 665.0; /* chlorophyll-a; 440 nm is the second major band */
    case 3: return 490.0; /* fucoxanthin */
    case 4: return 565.0; /* phycoerythrin */
    default: return 0.0;
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


static AED_REAL ptm_best_light_for_pigment(int pigment_id, AED_REAL light_440, AED_REAL light_490,
                                           AED_REAL light_550, AED_REAL light_565, AED_REAL light_590,
                                           AED_REAL light_620, AED_REAL light_665)
{
    switch (pigment_id) {
    case 1:
        return light_620;
    case 2:
        return (light_440 > light_665) ? light_440 : light_665;
    case 3:
        return (light_490 > light_550) ? light_490 : light_550;
    case 4:
        return (light_565 > light_590) ? light_565 : light_590;
    default:
        return 0.0;
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/* ptm_pigment_weight() lived here: a 7-entry hard-coded pigment response table applied to
 * radiant energy, with no wavelength-to-photon conversion and no connection to the 33-band
 * a*(lambda) spectra OASIM already carries in oasim.inc. It is replaced by
 * aed_sample_oasim_afac(), which reads OASIM's own spectral absorption factor. Deleted rather
 * than kept behind a flag: a second, disagreeing absorption model is a trap, not an option. */
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


static AED_REAL ptm_particle_par(int grp, int part)
{
    /* B2 2026-09-15: only ptm_write_glm calls this - the layer the height is in now, as the writer uses */
    int layer = (_PTM_Stat(grp,part,STAT) > 0) ? ptm_audit_layer(grp,part) : _PTM_Stat(grp,part,LAYR);
    if (ptm_idx_par < 0) ptm_update_state_indices();
    if (_PTM_Stat(grp, part, STAT) > 0) {
        AED_REAL par = ptm_state_or_zero(grp, part, ptm_idx_par);
        if (par != 0.0) return par;
    }
    if (layer >= botmLayer && layer < NumLayers) return Lake[layer].Light * 0.45;
    return 0.0;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


static void ptm_particle_oasim_bands(int layer, AED_REAL bands[PTM_OASIM_NBANDS], int *sample_status)
{
    int i, status = 2;
    for (i = 0; i < PTM_OASIM_NBANDS; i++) bands[i] = 0.0;
    aed_sample_oasim_particle_bands(layer, bands, PTM_OASIM_NBANDS, &status);
    if (sample_status != NULL) *sample_status = status;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


void ptm_update_state_indices(void)
{
    static int reported_absent = 0;
    int b;

    ptm_idx_par = aed_phyto_abm_ipar_index();
    if (Num_PTM_Vars < 1 || ptm_idx_par < 1) {
        /* No aed_phyto_abm state registered (physical-only particles, or the ABM is not in
         * the model list): every biology read below returns 0 and the buoyancy model is off.
         * ptm_idx_par = 0 also stops the lazy re-resolution (it tests for < 0). */
        ptm_idx_par = 0;
        ptm_idx_par_scalar_depth = ptm_idx_par_layer_center = ptm_idx_par_depth_factor = 0;
        ptm_idx_par_scalar_status = ptm_idx_par_absorbed = ptm_idx_cldose = 0;
        ptm_idx_dose_scal_l = ptm_idx_dose_eff_l = ptm_idx_dprodC_l = ptm_idx_dday = ptm_idx_fcol = 0;
        ptm_idx_buoy_rate = ptm_idx_lightdose = ptm_idx_s_act = ptm_idx_buoy_I = ptm_idx_buoy_stat = 0;
        for (b = 0; b < PTM_OASIM_NBANDS; b++) ptm_idx_par_depth_band[b] = 0;
        for (b = 0; b < PAM_N_LEGACY_C; b++) ptm_idx_legacy[b] = 0;
        buoy_any = 0;
        if (!reported_absent) {
            fprintf(stderr, "     PTM: no aed_phyto_abm particle state registered; "
                            "particle biology outputs are written as 0\n");
            reported_absent = 1;
        }
        return;
    }
    ptm_idx_par_scalar_depth = aed_phyto_abm_ipar_scalar_depth_index();
    ptm_idx_par_layer_center = aed_phyto_abm_ipar_layer_center_index();
    ptm_idx_par_depth_factor = aed_phyto_abm_ipar_depth_factor_index();
    ptm_idx_par_scalar_status = aed_phyto_abm_ipar_scalar_status_index();
    for (b = 0; b < PTM_OASIM_NBANDS; b++)
        ptm_idx_par_depth_band[b] = aed_phyto_abm_ipar_depth_band_index(b + 1);
    ptm_idx_par_absorbed = aed_phyto_abm_ipar_absorbed_index();
    ptm_idx_cldose = aed_phyto_abm_ipar_cldose_index();
    ptm_idx_dose_scal_l = aed_phyto_abm_idose_scal_l_index();
    ptm_idx_dose_eff_l = aed_phyto_abm_idose_eff_l_index();
    ptm_idx_dprodC_l = aed_phyto_abm_idprodC_l_index();
    ptm_idx_dday = aed_phyto_abm_idday_index();

    ptm_idx_fcol = aed_phyto_abm_ipar_fcol_index();

    /* The legacy block. All 18 are registered unconditionally by aed_phyto_abm_define, so
     * either the whole block resolves or the ABM is absent and none of it does. A PARTIAL
     * block is neither, and means a registration was made conditional without the C side
     * being updated - report it once and loudly, because the failure it would otherwise
     * cause is silent wrong numbers. */
    {
        int i, nset = 0;
        for (i = 0; i < PAM_N_LEGACY_C; i++) ptm_idx_legacy[i] = 0;
        aed_phyto_abm_legacy_ptm_indices(PAM_N_LEGACY_C, ptm_idx_legacy);
        for (i = 0; i < PAM_N_LEGACY_C; i++) if (ptm_idx_legacy[i] >= 1) nset++;
        if (nset != 0 && nset != PAM_N_LEGACY_C) {
            fprintf(stderr, "PTM LEGACY INDEX BLOCK INCOMPLETE: %d of %d resolved.\n"
                            "  A legacy ABM variable is no longer registered unconditionally;"
                            " glm_ptm.c must be updated to match.\n",
                    nset, PAM_N_LEGACY_C);
            exit(1);
        }
    }

    /* Wallace & Hamilton lagged buoyancy: state indices and per-group parameters. */
    ptm_idx_buoy_rate = aed_phyto_abm_ibuoy_rate_index();
    ptm_idx_lightdose = aed_phyto_abm_ilightdose_index();
    ptm_idx_s_act     = aed_phyto_abm_is_act_index();
    ptm_idx_buoy_I    = aed_phyto_abm_ibuoy_I_index();
    ptm_idx_buoy_stat = aed_phyto_abm_ibuoy_stat_index();
    buoy_any = 0;
    if (num_particle_groups > MAX_PTM_GROUPS_) {
        fprintf(stderr, "PTM BUOYANCY: num_particle_groups %d exceeds MAX_PTM_GROUPS_ %d\n",
                num_particle_groups, MAX_PTM_GROUPS_);
        exit(1);
    }
    for (int g = 0; g < num_particle_groups; g++) {
        int m = 0;
        aed_phyto_abm_buoy_params(g + 1, &m, &buoy_prm[g][0]);   /* Fortran side is 1-based */
        buoy_model[g] = m;
        if (m == 1) buoy_any = 1;
    }
    if (buoy_any) {
        if (ptm_idx_buoy_rate < 1 || ptm_idx_buoy_rate > Num_PTM_Vars ||
            ptm_idx_lightdose < 1 || ptm_idx_lightdose > Num_PTM_Vars ||
            ptm_idx_s_act     < 1 || ptm_idx_s_act     > Num_PTM_Vars) {
            fprintf(stderr, "PTM BUOYANCY INDEX ERROR: %d %d %d (Num_PTM_Vars=%d)\n",
                    ptm_idx_buoy_rate, ptm_idx_lightdose, ptm_idx_s_act, Num_PTM_Vars);
            exit(1);
        }
        if (ptm_buoy_substep < 1) ptm_buoy_substep = 1;
        /* D3: the dawn detector is global but I_dark is per-group; use the smallest active value so
         * no group's light branch can open before the day has been declared. */
        buoy_I_dark = 1.0e30;
        for (int g = 0; g < num_particle_groups; g++)
            if (buoy_model[g] == 1 && buoy_prm[g][5] < buoy_I_dark) buoy_I_dark = buoy_prm[g][5];
        if (!(buoy_I_dark > 0.0) || buoy_I_dark > 1.0e29) buoy_I_dark = 16.41;
        for (int g = 0; g < num_particle_groups; g++) if (buoy_model[g] == 1)
            fprintf(stderr, "PTM BUOYANCY: group %d Wallace-Hamilton ON  "
                    "c1=%.3f c2=%.5f c3=%.6f KI=%.2f tau=%.1fmin Idark=%.2f rho=[%.1f,%.1f] "
                    "r_deball=%.3f k_resp=%.4f\n",
                    g + 1, buoy_prm[g][0], buoy_prm[g][1], buoy_prm[g][2], buoy_prm[g][3],
                    buoy_prm[g][4], buoy_prm[g][5], buoy_prm[g][7], buoy_prm[g][8],
                    buoy_prm[g][9], buoy_prm[g][10]);
    }

    if (ptm_idx_par != 1) {
        fprintf(stderr, "PTM STATE INDEX ERROR: ip_par is %d, expected 1\n", ptm_idx_par);
        exit(1);
    }
    if (ptm_idx_par_scalar_depth < 1 || ptm_idx_par_scalar_depth > Num_PTM_Vars ||
        ptm_idx_par_layer_center < 1 || ptm_idx_par_layer_center > Num_PTM_Vars ||
        ptm_idx_par_depth_factor < 1 || ptm_idx_par_depth_factor > Num_PTM_Vars ||
        ptm_idx_par_scalar_status < 1 || ptm_idx_par_scalar_status > Num_PTM_Vars) {
        fprintf(stderr, "PTM STATE INDEX ERROR: invalid Tier 1D indices %d %d %d %d (Num_PTM_Vars=%d)\n",
                ptm_idx_par_scalar_depth, ptm_idx_par_layer_center, ptm_idx_par_depth_factor, ptm_idx_par_scalar_status, Num_PTM_Vars);
        exit(1);
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


void ptm_update_particle_par(void)
{
    int pg, p, layer;
    AED_REAL h, mean_h, par_center, kd, dz_bot, dz_top, factor, par_scalar;
    int status;

    /* 2026-09-11 (P2): aed_do_glm calls this every step whether or not the particle model is
     * on; with ptm_sw = .false. the particle arrays were never allocated and the loop below
     * dereferenced them - segmentation fault on the first step, in the paper binary too
     * (found by the runtime toggle matrix). No-op without particles. */
    if (!ptm_sw) return;                     /* particles intentionally off: nothing to do */
    if (PTM_Vars == NULL || PTM_Stat == NULL || num_particle_groups <= 0) {
        /* particles ENABLED but their arrays are missing: a required calculation cannot
         * run, and silently skipping it would hide a broken initialisation. */
        fprintf(stderr, "FATAL: ptm_update_particle_par: particles are enabled (ptm_sw) but the "
                        "particle arrays are not allocated (PTM_Vars=%p PTM_Stat=%p groups=%d)\n",
                        (void*)PTM_Vars, (void*)PTM_Stat, num_particle_groups);
        exit(1);
    }

    if (ptm_idx_par < 0) ptm_update_state_indices();
    if (ptm_idx_par < 1) return;             /* no ABM state to update */

    for (pg = 0; pg < num_particle_groups; pg++) {
        for (p = 0; p < max_particle_num; p++) {
            _PTM_Vars(pg, p, PTM_STATE_VAR(ptm_idx_par_scalar_status)) = 1.0; /* default fallback */

            if (_PTM_Stat(pg, p, STAT) <= 0) continue;

            layer = _PTM_Stat(pg, p, LAYR);
            if (layer < botmLayer || layer >= NumLayers) continue;

            h = _PTM_Vars(pg, p, HGHT);
            /* clamp particle height to layer bounds */
            dz_top = Lake[layer].Height;
            dz_bot = (layer > botmLayer) ? Lake[layer - 1].Height : 0.0;
            if (h > dz_top) h = dz_top;
            if (h < dz_bot) h = dz_bot;

            mean_h = Lake[layer].MeanHeight;

            aed_sample_oasim_par_kd(layer, &par_center, &kd, &status);

            if (status != 0 || par_center <= 0.0 || kd < 0.0 || !isfinite(par_center) || !isfinite(kd)) {
                par_center = Lake[layer].Light * 0.45;
                kd = 0.0;
                factor = 1.0;
                par_scalar = par_center;
                _PTM_Vars(pg, p, PTM_STATE_VAR(ptm_idx_par_scalar_status)) = 1.0;
            } else {
                factor = exp(kd * (h - mean_h));
                par_scalar = par_center * factor;
                _PTM_Vars(pg, p, PTM_STATE_VAR(ptm_idx_par_scalar_status)) = 0.0;
            }

            _PTM_Vars(pg, p, PTM_STATE_VAR(ptm_idx_par_scalar_depth)) = par_scalar;
            _PTM_Vars(pg, p, PTM_STATE_VAR(ptm_idx_par_layer_center)) = par_center;
            _PTM_Vars(pg, p, PTM_STATE_VAR(ptm_idx_par_depth_factor)) = factor;

            /* ---- Phase A/B: 7-band depth attenuation + absorbed-equivalent PAR ----
             * Computed HERE, in the physics update, not in the NetCDF write path, so the
             * band values stored as PTM state are provably the ones biology consumed.
             * Sampling at output time would write numbers the model never used - the
             * fabricated-vs-measured ambiguity G9 exists to prevent.
             *
             * absorbed = pigment-weighted mean band irradiance, renormalised onto the
             * scalar depth PAR. Mode 2 must change the spectral WEIGHTING of the light,
             * not its magnitude; without the renormalisation it would silently rescale
             * total light and be indistinguishable from a growth-rate change.
             */
            if (ptm_idx_par_absorbed >= 1 && ptm_idx_par_depth_band[0] >= 1) {
                AED_REAL bands[PTM_OASIM_NBANDS], kds[PTM_OASIM_NBANDS];
                AED_REAL absorbed = 0.0;
                int bstat = 2, kstat = 2, b;

                ptm_particle_oasim_bands(layer, bands, &bstat);
                aed_sample_oasim_band_kd(layer, kds, PTM_OASIM_NBANDS, &kstat);

                /* The 7 bands are still attenuated to the particle's depth and still written
                 * out - they are a real diagnostic of the spectrum at that depth. They no
                 * longer decide what the organism absorbs. */
                for (b = 0; b < PTM_OASIM_NBANDS; b++) {
                    AED_REAL eb, kb;
                    if (bstat == 0 && isfinite(bands[b]) && bands[b] >= 0.0) {
                        /* per-band Kd when available, else bulk scalar Kd, else none:
                         * degrade gracefully rather than zeroing the band */
                        kb = (kstat == 0) ? kds[b] : kd;
                        eb = bands[b] * exp(kb * (h - mean_h));
                        if (!isfinite(eb) || eb < 0.0) eb = 0.0;
                    } else {
                        eb = 0.0;
                    }
                    _PTM_Vars(pg, p, PTM_STATE_VAR(ptm_idx_par_depth_band[b])) = eb;
                }

                /* Absorbed light now comes from OASIM's own spectral absorption factor.
                 *
                 * afac = <a*(lambda) weighted by the LOCAL spectrum> / <a*(lambda) flat>,
                 * integrated over all 33 PAR bands with proper nm weights, from the same
                 * chlorophyll-specific spectra that drive shading. With afac_quanta it is
                 * weighted by PHOTON flux, Phi = I*lambda/(h c), rather than by energy.
                 *
                 * This replaces a hard-coded 7-entry weight table that was disconnected from
                 * those spectra, addressed its bands by POSITION in lambda_out, and combined
                 * them with an unweighted mean despite spacings of 15 to 60 nm.
                 *
                 * Still dimensionless and still normalised to 1.0 for a flat spectrum, so
                 * mode 2 continues to change the spectral WEIGHTING of the light and not its
                 * magnitude - no recalibration of the growth parameters is implied.
                 *
                 * Group -> IOP mapping is the identity, matching spectral_iop = 1,2,3,4 in
                 * &aed_phyto_abm. A failed sample returns afac = 1.0, leaving light unchanged.
                 */
                {
                    AED_REAL afac = 1.0;
                    int astat = 2;
                    aed_sample_oasim_afac(layer, pg + 1, &afac, &astat);
                    if (astat != 0 || !isfinite(afac) || afac <= 0.0) {
                        static int afac_warned = 0;
                        if (!afac_warned) {
                            fprintf(stderr, "     PTM: OASIM spectral absorption factor OAS_afac_iop%d is not available "
                                            "(status %d): particle absorbed PAR = scalar PAR (afac = 1) for this run. "
                                            "It needs a libaed-light that publishes the OAS_afac_iop<n> diagnostics.\n", pg + 1, astat);
                            afac_warned = 1;
                        }
                        afac = 1.0;
                    }
                    absorbed = par_scalar * afac;
                }
                _PTM_Vars(pg, p, PTM_STATE_VAR(ptm_idx_par_absorbed)) =
                    (isfinite(absorbed) && absorbed > 0.0) ? absorbed : 0.0;
            }
        }
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/* Group names and pigment types come from the &particles lists when given
 * (particle_group_names, particle_pigment_type; glm_init.c records their lengths), with a
 * generic fallback so a run that names nothing still gets distinct labels. */
static const char *ptm_group_name(int grp)
{
    static char buf[32];
    if (particle_group_names != NULL && grp >= 0 && grp < particle_group_names_n &&
        particle_group_names[grp] != NULL && particle_group_names[grp][0] != '\0')
        return particle_group_names[grp];
    snprintf(buf, sizeof(buf), "group%d", grp + 1);
    return buf;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/* Pigment id, matching the legend written on particle_pigment_id:
 *   1 phycocyanin, 2 chlorophyll_a, 3 fucoxanthin, 4 phycoerythrin; otherwise grp+1. */
static int ptm_pigment_id(int grp)
{
    if (particle_pigment_type != NULL && grp >= 0 && grp < particle_pigment_type_n &&
        particle_pigment_type[grp] != NULL) {
        const char *s = particle_pigment_type[grp];
        if (strcmp(s, "phycocyanin") == 0)   return 1;
        if (strcmp(s, "chlorophyll_a") == 0) return 2;
        if (strcmp(s, "fucoxanthin") == 0)   return 3;
        if (strcmp(s, "phycoerythrin") == 0) return 4;
    }
    return grp + 1;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/* Behaviour id: the interpreted vertical strategy of the pigment class (same legend order). */
static int ptm_behavior_id(int grp)
{
    return ptm_pigment_id(grp);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


static void ptm_group_names_attr(char *buf, size_t len)
{
    int pg;
    size_t used = 0;
    if (len == 0) return;
    buf[0] = '\0';
    for (pg = 0; pg < num_particle_groups; pg++) {
        const char *name = ptm_group_name(pg);
        int n = snprintf(buf + used, len - used, "%s%s", (pg == 0) ? "" : ",", name);
        if (n < 0 || (size_t)n >= len - used) {
            buf[len-1] = '\0';
            return;
        }
        used += (size_t)n;
    }
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *   This routine moves a particle by a random-walk diffusion step plus a     *
 *   directional sinking/floating step (vvel * dt_secs)                       *
 *                                                                            *
 ******************************************************************************/
AED_REAL move_particle(AED_REAL dt_secs, AED_REAL Height, AED_REAL K_z, AED_REAL K_prime_z, AED_REAL vvel, AED_REAL rand_draw)
{
//LOCALS

    AED_REAL updated_height;
    AED_REAL K_half;
    const AED_REAL r_var = 1.0/3.0;   // variance of a uniform[-1,1] draw

/*----------------------------------------------------------------------------*/
//BEGIN

    // rand_draw: a uniform draw in [-1,1], supplied by the caller (do_ptm_update) so that
    // every draw of the run passes through one counted stream (restart replay).

    // Visser (1997) random-walk correction for depth-varying diffusivity:
    //   z(t+dt) = z(t) + K'(z)*dt + R * sqrt(2 * K(z + 0.5*K'(z)*dt) * dt / r)
    // K'(z) (K_prime_z, the signed central dK/dHeight from do_ptm_update()) is the LOCAL
    // RATE OF CHANGE of diffusivity and K_z is the diffusivity interpolated to the
    // particle's height. With K locally linear the half-step evaluation is
    //   K(z + 0.5*K'*dt) ~ K_z + 0.5*K'^2*dt
    K_half = K_z + 0.5 * K_prime_z * K_prime_z * dt_secs;
    if (K_half < ptm_diffusivity) K_half = K_z;

    updated_height = Height + K_prime_z * dt_secs
                   + rand_draw * sqrt((2.0 * K_half * dt_secs) / r_var);   // random walk

    updated_height = updated_height + vvel * dt_secs;   // account for sinking/floating;
                                                        // vvel is per second so needs to
                                                        // be multiplied by dt_secs, the
                                                        // duration of this random walk
                                                        // substep

    return updated_height;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *        This routine returns the settling velocity for a particle in m/s    *
 *                                                                            *
 ******************************************************************************/
AED_REAL get_settling_velocity(AED_REAL settling_velocity)
{
    return settling_velocity / 86400;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *        This routine returns the density for a particle                     *
 *                                                                            *
 ******************************************************************************/
AED_REAL get_particle_density(AED_REAL particle_density)
{
    return particle_density;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 *        This routine returns the diameter for a particle                    *
 *                                                                            *
 ******************************************************************************/
AED_REAL get_particle_diameter(AED_REAL particle_diameter)
{
    return particle_diameter;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
static AED_REAL ptm_var_or_zero(int grp, int part, int var)
{
    if (var < (PTM_ENV_NVARS + Num_PTM_Vars)) return _PTM_Vars(grp, part, var);
    return 0.0;
}

/* One ABM state variable by its 1-based Fortran index; 0.0 when it is not registered
 * (index 0 would otherwise alias the last environment slot, HGHT). */
static AED_REAL ptm_state_or_zero(int grp, int part, int fidx)
{
    if (fidx < 1) return 0.0;
    return ptm_var_or_zero(grp, part, PTM_STATE_VAR(fidx));
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


static int dose_scal_l_id, dose_eff_l_id, dprodC_l_id, dday_id;
static int h_id, m_id, d_id, dn_id, vv_id, par_id, tem_id, no3_id, nh4_id, frp_id, c_id, n_id, pho_id, chl_id, num_id, cdiv_id, topt_id, lnalphachl_id, npq_id, cldose_id, fcol_id, stat_id, flag_id, ptid_id, grp_id;
/* Life-history accumulators */
static int age_id, birth_id, cumpar_id, cumpar_mol_id, cumc_id, ndiv_id;
static int group_id, pigment_id, phyto_link_id, particle_layer_id, oasim_status_id, light_440_id, light_490_id, light_550_id, light_565_id, light_590_id, light_620_id, light_665_id, spectral_score_id;
static int depth_m_id, delta_depth_m_id, vertical_mode_id, behavior_id, light_preference_nm_id, best_light_id;
static int par_scalar_depth_id, par_layer_center_id, par_depth_factor_id, par_scalar_status_id;
static int buoy_rate_id = -1, lightdose_id = -1, s_act_id = -1;
static int buoy_I_id = -1, buoy_stat_id = -1;
static int par_depth_band_id[PTM_OASIM_NBANDS], par_absorbed_id;
static int set_no_p = -1;
static size_t start[2],edges[2];
static AED_REAL *previous_particle_depth = NULL;
static int previous_particle_depth_count = 0;
/* 2026-09-10 (E5c): true once previous_particle_depth holds a real 'last written' depth -
 * after the first output record of this run, or after a restart restored it. Replaces the
 * `set_no_p == 0` test so a resumed run reports the displacement since the previous record
 * of the run it continues, instead of a spurious 0. Same values for a run that never
 * restarts: set_no_p == 0 and !prev_depth_valid coincide there. */
static int prev_depth_valid = 0;
int  ptm_outmem_len(void) { return previous_particle_depth_count; }
int  ptm_outmem_valid(void) { return prev_depth_valid; }
void ptm_outmem_get(AED_REAL *buf)
{ if (previous_particle_depth != NULL && previous_particle_depth_count > 0)
      memcpy(buf, previous_particle_depth, (size_t)previous_particle_depth_count * sizeof(AED_REAL)); }
void ptm_outmem_set(const AED_REAL *buf, int n, int valid)
{
    if (n <= 0) return;
    if (previous_particle_depth_count != n) {
        free(previous_particle_depth);
        previous_particle_depth = calloc((size_t)n, sizeof(AED_REAL));
        previous_particle_depth_count = previous_particle_depth ? n : 0;
    }
    if (previous_particle_depth == NULL) return;
    memcpy(previous_particle_depth, buf, (size_t)n * sizeof(AED_REAL));
    prev_depth_valid = valid;
}


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
//void ptm_write_glm(int *ncid, int *wlev, int *nlev, int *lvl, int *point_nlevs)
void ptm_write_glm(int ncid, int max_particle_num)
{
//LOCALS
    int p, pg, idx, total_particles, pigment;
    AED_REAL base_par, sum_weights, sampled_bands[PTM_OASIM_NBANDS];
    int sample_status;
    AED_REAL *p_height, *depth_m, *delta_depth_m, *mass, *diam, *density, *vvel, *par, *tem, *no3, *nh4, *frp, *c, *n, *pho, *chl, *num, *cdiv, *topt, *lnalphachl;
    AED_REAL *age, *birth, *cumpar, *cumpar_mol, *cumc, *ndiv, *npq, *cldose, *fcol;
    AED_REAL *dose_scal_l, *dose_eff_l, *dprodC_l, *dday;

    AED_REAL *light_440, *light_490, *light_550, *light_565, *light_590, *light_620, *light_665, *spectral_score, *light_preference_nm, *best_light;
    AED_REAL *par_scalar_depth, *par_layer_center, *par_depth_factor;
    AED_REAL *buoy_rate = NULL, *lightdose = NULL, *s_act = NULL;
    AED_REAL *buoy_I = NULL, *buoy_statv = NULL;
    AED_REAL *par_depth_band[PTM_OASIM_NBANDS], *par_absorbed;
    int *status, *flag, *ptid, *group, *pigment_arr, *phyto_link, *particle_layer, *oasim_status, *vertical_mode, *behavior;
    int *par_scalar_status, *grp;

/*----------------------------------------------------------------------------*/
//BEGIN

    set_no_p++;
    if (set_no_p > 0) prev_depth_valid = 1;   /* E5c: a previous record now exists */

    if (ptm_idx_par < 0) ptm_update_state_indices();

    total_particles = ptm_total_particles();
    start[1] = 0;             edges[1] = total_particles;
    start[0] = set_no_p;      edges[0] = 1;

    p_height  = malloc(total_particles*sizeof(AED_REAL));
    depth_m = malloc(total_particles*sizeof(AED_REAL));
    delta_depth_m = malloc(total_particles*sizeof(AED_REAL));
    mass  = malloc(total_particles*sizeof(AED_REAL));
    diam  = malloc(total_particles*sizeof(AED_REAL));
    density  = malloc(total_particles*sizeof(AED_REAL));
    vvel  = malloc(total_particles*sizeof(AED_REAL));
    par  = malloc(total_particles*sizeof(AED_REAL));
    tem  = malloc(total_particles*sizeof(AED_REAL));
    no3  = malloc(total_particles*sizeof(AED_REAL));
    nh4  = malloc(total_particles*sizeof(AED_REAL));
    frp  = malloc(total_particles*sizeof(AED_REAL));
    c  = malloc(total_particles*sizeof(AED_REAL));
    n  = malloc(total_particles*sizeof(AED_REAL));
    pho  = malloc(total_particles*sizeof(AED_REAL));
    chl  = malloc(total_particles*sizeof(AED_REAL));
    num  = malloc(total_particles*sizeof(AED_REAL));
    cdiv  = malloc(total_particles*sizeof(AED_REAL));
    topt  = malloc(total_particles*sizeof(AED_REAL));
    lnalphachl  = malloc(total_particles*sizeof(AED_REAL));
    age         = malloc(total_particles*sizeof(AED_REAL));
    birth       = malloc(total_particles*sizeof(AED_REAL));
    cumpar      = malloc(total_particles*sizeof(AED_REAL));
    cumpar_mol  = malloc(total_particles*sizeof(AED_REAL));
    cumc        = malloc(total_particles*sizeof(AED_REAL));
    ndiv        = malloc(total_particles*sizeof(AED_REAL));
    npq         = malloc(total_particles*sizeof(AED_REAL));
    cldose      = malloc(total_particles*sizeof(AED_REAL));
    dose_scal_l = malloc(total_particles*sizeof(AED_REAL));
    dose_eff_l = malloc(total_particles*sizeof(AED_REAL));
    dprodC_l = malloc(total_particles*sizeof(AED_REAL));
    dday = malloc(total_particles*sizeof(AED_REAL));

    fcol        = malloc(total_particles*sizeof(AED_REAL));
    light_440 = malloc(total_particles*sizeof(AED_REAL));
    light_490 = malloc(total_particles*sizeof(AED_REAL));
    light_550 = malloc(total_particles*sizeof(AED_REAL));
    light_565 = malloc(total_particles*sizeof(AED_REAL));
    light_590 = malloc(total_particles*sizeof(AED_REAL));
    light_620 = malloc(total_particles*sizeof(AED_REAL));
    light_665 = malloc(total_particles*sizeof(AED_REAL));
    spectral_score = malloc(total_particles*sizeof(AED_REAL));
    light_preference_nm = malloc(total_particles*sizeof(AED_REAL));
    best_light = malloc(total_particles*sizeof(AED_REAL));
    par_scalar_depth = malloc(total_particles*sizeof(AED_REAL));
    if (buoy_any) {
        buoy_rate = malloc(total_particles*sizeof(AED_REAL));
        lightdose = malloc(total_particles*sizeof(AED_REAL));
        s_act     = malloc(total_particles*sizeof(AED_REAL));
        buoy_I    = malloc(total_particles*sizeof(AED_REAL));
        buoy_statv= malloc(total_particles*sizeof(AED_REAL));
    }
    for (int b = 0; b < PTM_OASIM_NBANDS; b++)
        par_depth_band[b] = malloc(total_particles*sizeof(AED_REAL));
    par_absorbed = malloc(total_particles*sizeof(AED_REAL));
    par_layer_center = malloc(total_particles*sizeof(AED_REAL));
    par_depth_factor = malloc(total_particles*sizeof(AED_REAL));
    par_scalar_status = malloc(total_particles*sizeof(int));
    status  = malloc(total_particles*sizeof(int));
    flag  = malloc(total_particles*sizeof(int));
    ptid  = malloc(total_particles*sizeof(int));
    group = malloc(total_particles*sizeof(int));
    pigment_arr = malloc(total_particles*sizeof(int));
    phyto_link = malloc(total_particles*sizeof(int));
    particle_layer = malloc(total_particles*sizeof(int));
    oasim_status = malloc(total_particles*sizeof(int));
    vertical_mode = malloc(total_particles*sizeof(int));
    behavior = malloc(total_particles*sizeof(int));
    grp = malloc(total_particles*sizeof(int));

    if (previous_particle_depth_count != total_particles) {
        free(previous_particle_depth);
        previous_particle_depth = calloc(total_particles, sizeof(AED_REAL));
        previous_particle_depth_count = total_particles;
    }

    for (pg = 0; pg < num_particle_groups; pg++) {
      pigment = ptm_pigment_id(pg);
      for (p = 0; p < max_particle_num; p++) {
        idx = pg * max_particle_num + p;
        p_height[idx]       = _PTM_Vars(pg,p,HGHT);    //Particle[p].Height;                REAL
        /* B2 2026-09-15: LAYR was set by the last do_ptm_update and the column may have been restructured
         * since (the daily block's check_layer_thickness precedes every end-of-day record), so it can name a
         * different layer than the H/NS written beside it. Scan the height against the current layers, by
         * ptm_update_layerid's rule and for the particles it updates, and use that for particle_layer, the
         * OASIM sample and the temperature fallback (ptm_particle_par does the same). LAYR is not changed. */
        int cur_layer = (_PTM_Stat(pg,p,STAT) > 0) ? ptm_audit_layer(pg,p) : _PTM_Stat(pg,p,LAYR);
        mass[idx]           = _PTM_Vars(pg,p,MASS);    //Particle[p].Mass;                  REAL
        diam[idx]           = _PTM_Vars(pg,p,DIAM);    // Particle[p].Diam;                 REAL
        density[idx]        = _PTM_Vars(pg,p,DENS);    //Particle[p].Density;               REAL
        vvel[idx]           = _PTM_Vars(pg,p,VVEL)*86400;  //Particle[p].vvel;              REAL
                                         //  VVEL+1 = HGHT
         par[idx]            = ptm_particle_par(pg,p);
         par_scalar_depth[idx] = ptm_state_or_zero(pg, p, ptm_idx_par_scalar_depth);
         if (buoy_any) {
             buoy_rate[idx] = ptm_state_or_zero(pg, p, ptm_idx_buoy_rate);
             lightdose[idx] = ptm_state_or_zero(pg, p, ptm_idx_lightdose);
             s_act[idx]     = ptm_state_or_zero(pg, p, ptm_idx_s_act);
             buoy_I[idx]    = ptm_state_or_zero(pg, p, ptm_idx_buoy_I);
             buoy_statv[idx]= ptm_state_or_zero(pg, p, ptm_idx_buoy_stat);
         }
        for (int b = 0; b < PTM_OASIM_NBANDS; b++)
            par_depth_band[b][idx] = (ptm_idx_par_depth_band[b] >= 1)
                 ? ptm_state_or_zero(pg, p, ptm_idx_par_depth_band[b]) : 0.0;
        /* 2026-08: the -9999 "never computed" sentinel used to reach the file here.
         * The update path above already writes 0.0 when the absorbed irradiance cannot be
         * formed (no OASIM bands, oasim_status 3), but it does not run for a particle
         * created by splitting before its first light update, so that particle's state was
         * still at the sentinel when the record was written. Growth and NPQ were never
         * affected - aed_phyto_abm.F90:1747 only consumes this value when it is > 0 - but
         * the exported diagnostic was, and -9999 is not declared as a _FillValue, so any
         * analysis filtering only on |x| < 1e30 silently averaged it in. Apply the same
         * convention as the update path and export 0.0 for "not computed". */
        par_absorbed[idx] = (ptm_idx_par_absorbed >= 1)
                 ? ptm_state_or_zero(pg, p, ptm_idx_par_absorbed) : 0.0;
        if ( !isfinite(par_absorbed[idx]) || par_absorbed[idx] < 0.0 )
            par_absorbed[idx] = 0.0;
         par_layer_center[idx] = ptm_state_or_zero(pg, p, ptm_idx_par_layer_center);
         par_depth_factor[idx] = ptm_state_or_zero(pg, p, ptm_idx_par_depth_factor);
         par_scalar_status[idx] = (int)ptm_state_or_zero(pg, p, ptm_idx_par_scalar_status);
         tem[idx]            = ptm_legacy_var(pg,p,PAM_L_TEM);  //temp experienced by particle;      REAL
        if (tem[idx] == 0.0 && cur_layer >= botmLayer && cur_layer < NumLayers)   /* B2 */
            tem[idx] = Lake[cur_layer].Temp;
        no3[idx]            = ptm_legacy_var(pg,p,PAM_L_NO3);  //NO3 experienced by particle;       REAL
        nh4[idx]            = ptm_legacy_var(pg,p,PAM_L_NH4);  //NH4 experienced by particle;       REAL
        frp[idx]            = ptm_legacy_var(pg,p,PAM_L_FRP);  //FRP experienced by particle;       REAL
        c[idx]              = ptm_legacy_var(pg,p,PAM_L_C);  //internal particle C;               REAL
        n[idx]              = ptm_legacy_var(pg,p,PAM_L_N);  //internal particle N;               REAL
        pho[idx]            = ptm_legacy_var(pg,p,PAM_L_P);  //internal particle P;               REAL
        chl[idx]            = ptm_legacy_var(pg,p,PAM_L_CHL);  //internal particle chl;             REAL
        num[idx]            = ptm_legacy_var(pg,p,PAM_L_NUM);  //number of cells in particle;       REAL
        cdiv[idx]           = ptm_legacy_var(pg,p,PAM_L_CDIV);  //internal C threshold for division; REAL
        topt[idx]           = ptm_legacy_var(pg,p,PAM_L_TOPT);  //particle temperature optimum;      REAL
        lnalphachl[idx]     = ptm_legacy_var(pg,p,PAM_L_LNALPHACHL); //ln alpha chl of particle;          REAL
        age[idx]            = ptm_legacy_var(pg,p,PAM_L_AGE);  //age since seeding/division [d]
        birth[idx]          = ptm_legacy_var(pg,p,PAM_L_BIRTH);  //model day created [d]
        cumpar[idx]         = ptm_legacy_var(pg,p,PAM_L_CUMPAR);  //cumulative light dose [W/m2 d]
        /* mol photons m-2 = (W m-2 d) * 4.57 umol/J * 86400 s/d * 1e-6 */
        cumpar_mol[idx]     = cumpar[idx] * 4.57 * 86400.0 * 1.0e-6;
        cumc[idx]           = ptm_legacy_var(pg,p,PAM_L_CUMC);  //lifetime gross C fixed
        ndiv[idx]           = ptm_legacy_var(pg,p,PAM_L_NDIV);  //divisions undergone
        npq[idx]            = ptm_legacy_var(pg,p,PAM_L_NPQ);  //NPQ percent PAR reduction
        /* Ranjbar CL. Uses the computed index, not VVEL+22, so that inserting any
         * future ABM state variable cannot silently repoint this read. */
        cldose[idx]         = ptm_state_or_zero(pg, p, ptm_idx_cldose);
        dose_scal_l[idx] = (ptm_idx_dose_scal_l >= 1) ? ptm_state_or_zero(pg, p, ptm_idx_dose_scal_l) : 0.0;
        dose_eff_l[idx] = (ptm_idx_dose_eff_l >= 1) ? ptm_state_or_zero(pg, p, ptm_idx_dose_eff_l) : 0.0;
        dprodC_l[idx] = (ptm_idx_dprodC_l >= 1) ? ptm_state_or_zero(pg, p, ptm_idx_dprodC_l) : 0.0;
        dday[idx] = (ptm_idx_dday >= 1) ? ptm_state_or_zero(pg, p, ptm_idx_dday) : 0.0;

        fcol[idx]           = ptm_state_or_zero(pg, p, ptm_idx_fcol);
        status[idx]         = _PTM_Stat(pg,p,STAT);    //Particle[p].Status;                INT
        flag[idx]           = _PTM_Stat(pg,p,FLAG);    //Particle[p].Flag;                  INT
        ptid[idx]           = _PTM_Stat(pg,p,PTID);    //Particle[p].PTID;                  INT
        grp[idx]            = (_PTM_Stat(pg,p,GRP) >= 1) ? _PTM_Stat(pg,p,GRP) : pg + 1;   // the ABM's species field when it set one (split-born particles may not carry it), else the configured group
        group[idx]          = pg + 1;
        pigment_arr[idx]    = pigment;
        phyto_link[idx]     = pg + 1;
        particle_layer[idx] = cur_layer;   /* B2: the OASIM sample below reads this */
        depth_m[idx]        = Lake[surfLayer].Height - p_height[idx];
        if (depth_m[idx] < 0.0) depth_m[idx] = 0.0;
        vertical_mode[idx]  = ptm_vertical_mode(_PTM_Vars(pg,p,VVEL));
        behavior[idx]       = ptm_behavior_id(pg);

        light_preference_nm[idx] = ptm_light_preference_nm(pigment);

        ptm_particle_oasim_bands(particle_layer[idx], sampled_bands, &sample_status);
        oasim_status[idx] = sample_status;
        if (sample_status == 0 && status[idx] == 1) {
            light_440[idx] = sampled_bands[0];
            light_490[idx] = sampled_bands[1];
            light_550[idx] = sampled_bands[2];
            light_565[idx] = sampled_bands[3];
            light_590[idx] = sampled_bands[4];
            light_620[idx] = sampled_bands[5];
            light_665[idx] = sampled_bands[6];
        } else {
            /* No exact OASIM sample available. Previously this wrote
             * PAR * pigment_weight, a plausible-looking number indistinguishable
             * from measured OASIM data in the output file. Write the fill value
             * instead so it can never be mistaken for a spectral measurement.
             * particle_oasim_status records why. */
            base_par = NC_FILLER;
            light_440[idx] = base_par;
            light_490[idx] = base_par;
            light_550[idx] = base_par;
            light_565[idx] = base_par;
            light_590[idx] = base_par;
            light_620[idx] = base_par;
            light_665[idx] = base_par;
            if (status[idx] != 1 && sample_status == 0) oasim_status[idx] = 3;
        }
        if (light_440[idx] == NC_FILLER) {
            /* Derived quantities are meaningless without a real spectral sample. */
            spectral_score[idx] = NC_FILLER;
            best_light[idx] = NC_FILLER;
            if (previous_particle_depth != NULL && status[idx] == 1) {
                delta_depth_m[idx] = (!prev_depth_valid) ? 0.0 : depth_m[idx] - previous_particle_depth[idx];
                previous_particle_depth[idx] = depth_m[idx];
            } else {
                delta_depth_m[idx] = 0.0;
            }
            continue;
        }
        /* Spectral score is now the absorption factor the particle ACTUALLY experienced,
         * recovered from the two state variables biology used: absorbed / scalar-at-depth.
         * Previously it was recomputed here from the hard-coded weight table, so once the
         * biology path moved to OASIM's a*(lambda) this variable would have reported a
         * number the model never used - the fabricated-vs-measured trap the comment in
         * ptm_update_particle_par warns about. 1.0 means "this spectrum is neutral for this
         * organism"; >1 means the local spectrum favours it. */
        {
            AED_REAL pa = ptm_state_or_zero(pg, p, ptm_idx_par_absorbed);
            AED_REAL ps = par_scalar_depth[idx];
            spectral_score[idx] = (ps > 0.0 && isfinite(pa)) ? (pa / ps) : NC_FILLER;
        }
        (void)sum_weights;
        best_light[idx] = ptm_best_light_for_pigment(pigment, light_440[idx], light_490[idx],
                                                     light_550[idx], light_565[idx], light_590[idx],
                                                     light_620[idx], light_665[idx]);
        if (previous_particle_depth != NULL && status[idx] == 1) {
            delta_depth_m[idx] = (!prev_depth_valid) ? 0.0 : depth_m[idx] - previous_particle_depth[idx];
            previous_particle_depth[idx] = depth_m[idx];
        } else {
            delta_depth_m[idx] = 0.0;
        }
      }
    }

    check_nc_error(nc_put_vara(ncid, h_id, start, edges, p_height));
    check_nc_error(nc_put_vara(ncid, depth_m_id, start, edges, depth_m));
    check_nc_error(nc_put_vara(ncid, delta_depth_m_id, start, edges, delta_depth_m));
    check_nc_error(nc_put_vara(ncid, m_id, start, edges, mass));
    check_nc_error(nc_put_vara(ncid, d_id, start, edges, diam));
    check_nc_error(nc_put_vara(ncid, dn_id, start, edges, density));
    check_nc_error(nc_put_vara(ncid, vv_id, start, edges, vvel));
    check_nc_error(nc_put_vara(ncid, par_id, start, edges, par));
    check_nc_error(nc_put_vara(ncid, tem_id, start, edges, tem));
    check_nc_error(nc_put_vara(ncid, no3_id, start, edges, no3));
    check_nc_error(nc_put_vara(ncid, nh4_id, start, edges, nh4));
    check_nc_error(nc_put_vara(ncid, frp_id, start, edges, frp));
    check_nc_error(nc_put_vara(ncid, c_id, start, edges, c));
    check_nc_error(nc_put_vara(ncid, n_id, start, edges, n));
    check_nc_error(nc_put_vara(ncid, pho_id, start, edges, pho));
    check_nc_error(nc_put_vara(ncid, chl_id, start, edges, chl));
    check_nc_error(nc_put_vara(ncid, num_id, start, edges, num));
    check_nc_error(nc_put_vara(ncid, cdiv_id, start, edges, cdiv));
    check_nc_error(nc_put_vara(ncid, topt_id, start, edges, topt));
    check_nc_error(nc_put_vara(ncid, lnalphachl_id, start, edges, lnalphachl));
    check_nc_error(nc_put_vara(ncid, age_id, start, edges, age));
    check_nc_error(nc_put_vara(ncid, birth_id, start, edges, birth));
    check_nc_error(nc_put_vara(ncid, cumpar_id, start, edges, cumpar));
    check_nc_error(nc_put_vara(ncid, cumpar_mol_id, start, edges, cumpar_mol));
    check_nc_error(nc_put_vara(ncid, cumc_id, start, edges, cumc));
    check_nc_error(nc_put_vara(ncid, ndiv_id, start, edges, ndiv));
    check_nc_error(nc_put_vara(ncid, npq_id, start, edges, npq));
    check_nc_error(nc_put_vara(ncid, cldose_id, start, edges, cldose));
    if (ptm_idx_dose_scal_l >= 1) check_nc_error(nc_put_vara(ncid, dose_scal_l_id, start, edges, dose_scal_l));
    if (ptm_idx_dose_eff_l >= 1) check_nc_error(nc_put_vara(ncid, dose_eff_l_id, start, edges, dose_eff_l));
    if (ptm_idx_dprodC_l >= 1) check_nc_error(nc_put_vara(ncid, dprodC_l_id, start, edges, dprodC_l));
    if (ptm_idx_dday >= 1) check_nc_error(nc_put_vara(ncid, dday_id, start, edges, dday));

    check_nc_error(nc_put_vara(ncid, fcol_id, start, edges, fcol));
    check_nc_error(nc_put_vara(ncid, stat_id, start, edges, status));
    check_nc_error(nc_put_vara(ncid, flag_id, start, edges, flag));
    check_nc_error(nc_put_vara(ncid, ptid_id, start, edges, ptid));
    check_nc_error(nc_put_vara(ncid, grp_id, start, edges, grp));
    check_nc_error(nc_put_vara(ncid, group_id, start, edges, group));
    check_nc_error(nc_put_vara(ncid, pigment_id, start, edges, pigment_arr));
    check_nc_error(nc_put_vara(ncid, phyto_link_id, start, edges, phyto_link));
    check_nc_error(nc_put_vara(ncid, particle_layer_id, start, edges, particle_layer));
    check_nc_error(nc_put_vara(ncid, oasim_status_id, start, edges, oasim_status));
    check_nc_error(nc_put_vara(ncid, vertical_mode_id, start, edges, vertical_mode));
    check_nc_error(nc_put_vara(ncid, behavior_id, start, edges, behavior));
    check_nc_error(nc_put_vara(ncid, light_440_id, start, edges, light_440));
    check_nc_error(nc_put_vara(ncid, light_490_id, start, edges, light_490));
    check_nc_error(nc_put_vara(ncid, light_550_id, start, edges, light_550));
    check_nc_error(nc_put_vara(ncid, light_565_id, start, edges, light_565));
    check_nc_error(nc_put_vara(ncid, light_590_id, start, edges, light_590));
    check_nc_error(nc_put_vara(ncid, light_620_id, start, edges, light_620));
    check_nc_error(nc_put_vara(ncid, light_665_id, start, edges, light_665));
    check_nc_error(nc_put_vara(ncid, spectral_score_id, start, edges, spectral_score));
    check_nc_error(nc_put_vara(ncid, light_preference_nm_id, start, edges, light_preference_nm));
    check_nc_error(nc_put_vara(ncid, best_light_id, start, edges, best_light));

    check_nc_error(nc_put_vara(ncid, par_scalar_depth_id, start, edges, par_scalar_depth));
    if (buoy_any && buoy_rate_id >= 0) {
        check_nc_error(nc_put_vara(ncid, buoy_rate_id, start, edges, buoy_rate));
        check_nc_error(nc_put_vara(ncid, lightdose_id, start, edges, lightdose));
        check_nc_error(nc_put_vara(ncid, s_act_id,     start, edges, s_act));
        check_nc_error(nc_put_vara(ncid, buoy_I_id,    start, edges, buoy_I));
        check_nc_error(nc_put_vara(ncid, buoy_stat_id, start, edges, buoy_statv));
    }
    for (int b = 0; b < PTM_OASIM_NBANDS; b++)
        check_nc_error(nc_put_vara(ncid, par_depth_band_id[b], start, edges, par_depth_band[b]));
    check_nc_error(nc_put_vara(ncid, par_absorbed_id, start, edges, par_absorbed));
    check_nc_error(nc_put_vara(ncid, par_layer_center_id, start, edges, par_layer_center));
    check_nc_error(nc_put_vara(ncid, par_depth_factor_id, start, edges, par_depth_factor));
    check_nc_error(nc_put_vara(ncid, par_scalar_status_id, start, edges, par_scalar_status));

    free(p_height);
    free(depth_m);
    free(delta_depth_m);
    free(mass);
    free(diam);
    free(density);
    free(vvel);
    free(par);
    free(tem);
    free(no3);
    free(nh4);
    free(frp);
    free(c);
    free(n);
    free(pho);
    free(chl);
    free(num);
    free(cdiv);
    free(topt);
    free(lnalphachl);
    free(age);
    free(birth);
    free(cumpar);
    free(cumpar_mol);
    free(cumc);
    free(ndiv);
    free(npq);
    free(cldose);
    free(dose_scal_l);
    free(dose_eff_l);
    free(dprodC_l);
    free(dday);

    free(fcol);
    free(light_440);
    free(light_490);
    free(light_550);
    free(light_565);
    free(light_590);
    free(light_620);
    free(light_665);
    free(spectral_score);
    free(light_preference_nm);
    free(best_light);
    free(par_scalar_depth);
    if (buoy_rate) { free(buoy_rate); free(lightdose); free(s_act);
                     free(buoy_I); free(buoy_statv); }
    for (int b = 0; b < PTM_OASIM_NBANDS; b++) free(par_depth_band[b]);
    free(par_absorbed);
    free(par_layer_center);
    free(par_depth_factor);
    free(par_scalar_status);
    free(status);
    free(flag);
    free(ptid);
    free(group);
    free(pigment_arr);
    free(phyto_link);
    free(particle_layer);
    free(oasim_status);
    free(vertical_mode);
    free(behavior);
    free(grp);

    check_nc_error(nc_sync(ncid));
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 *                                                                            *
 ******************************************************************************/
// will also need to handle groups in write step; append group name onto
void ptm_init_glm_output(int ncid, int time_dim)
{
   int dims[2];
   char groups_attr[512];
   const char *pigment_legend = "1=phycocyanin/cyano,2=chlorophyll_a/green,3=fucoxanthin/diatom,4=phycoerythrin/planktothrix";

//
//------------------------------------------------------------------------------
//BEGIN
   define_mode_on(&ncid);   // Put NetCDF library in define mode.

   check_nc_error(nc_def_dim(ncid, "particles", ptm_total_particles(), &ptm_dim));
   fprintf(stderr, "     PTM NetCDF output active: writing %d particle slots and %d OASIM bands\n", ptm_total_particles(), PTM_OASIM_NBANDS);

   dims[1] = ptm_dim;
   dims[0] = time_dim;

   check_nc_error(nc_def_var(ncid, "particle_height", NC_REALTYPE, 2, dims, &h_id));
   set_nc_attributes(ncid, h_id, "meters", "Height of Particle" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_depth_m", NC_REALTYPE, 2, dims, &depth_m_id));
   set_nc_attributes(ncid, depth_m_id, "m", "Particle depth below lake surface" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_delta_depth_m", NC_REALTYPE, 2, dims, &delta_depth_m_id));
   set_nc_attributes(ncid, delta_depth_m_id, "m/output_step", "Particle depth change since previous NetCDF output; positive is deeper" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_mass", NC_REALTYPE, 2, dims, &m_id));
   set_nc_attributes(ncid, m_id, "grams", "Mass of Particle" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_diameter", NC_REALTYPE, 2, dims, &d_id));
   set_nc_attributes(ncid, d_id, "meters", "Diameter of Particle" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_density", NC_REALTYPE, 2, dims, &dn_id));
   set_nc_attributes(ncid, dn_id, "g/m3", "Density of Particle" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_vvel", NC_REALTYPE, 2, dims, &vv_id));
   /* Written as _PTM_Vars(..,VVEL)*86400 below, so the stored value is m/day, not m/s. */
   set_nc_attributes(ncid, vv_id, "m/day", "Settling Velocity of Particle (positive = upward)" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_par", NC_REALTYPE, 2, dims, &par_id));
   /* W/m2, matching the GLM 'par' global and GMK98_Ind_TempSizeLight's declared units. */
   set_nc_attributes(ncid, par_id, "W/m2", "particle layer PAR" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_tem", NC_REALTYPE, 2, dims, &tem_id));
   set_nc_attributes(ncid, tem_id, "degC", "particle layer temperature" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_no3", NC_REALTYPE, 2, dims, &no3_id));
   set_nc_attributes(ncid, no3_id, "mmolN", "particle layer NO3" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_nh4", NC_REALTYPE, 2, dims, &nh4_id));
   set_nc_attributes(ncid, nh4_id, "mmolN", "particle layer NH4" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_frp", NC_REALTYPE, 2, dims, &frp_id));
   set_nc_attributes(ncid, frp_id, "mmolP", "particle layer FRP" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_c", NC_REALTYPE, 2, dims, &c_id));
   set_nc_attributes(ncid, c_id, "pmolC/cell", "cell C concentration" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_n", NC_REALTYPE, 2, dims, &n_id));
   set_nc_attributes(ncid, n_id, "pmolN/cell", "cell N concentration" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_pho", NC_REALTYPE, 2, dims, &pho_id));
   set_nc_attributes(ncid, pho_id, "pmolP/cell", "cell P concentration" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_chl", NC_REALTYPE, 2, dims, &chl_id));
   set_nc_attributes(ncid, chl_id, "pgChl", "cell Chl concentration" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_num", NC_REALTYPE, 2, dims, &num_id));
   set_nc_attributes(ncid, num_id, "number", "number of cells/particle" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_cdiv", NC_REALTYPE, 2, dims, &cdiv_id));
   set_nc_attributes(ncid, cdiv_id, "pmol/cell", "cellular carbon content threshold for division" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_topt", NC_REALTYPE, 2, dims, &topt_id));
   set_nc_attributes(ncid, topt_id, "degC", "optimal temperature" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_lnalphachl", NC_REALTYPE, 2, dims, &lnalphachl_id));
   set_nc_attributes(ncid, lnalphachl_id, "(W m-2)-1(gChl molC)-1d-1", "slope of the P-I curve" PARAM_FILLVALUE);

   /* --- Life-history accumulators: make each particle trajectory reconstructable --- */
   check_nc_error(nc_def_var(ncid, "particle_age_d", NC_REALTYPE, 2, dims, &age_id));
   set_nc_attributes(ncid, age_id, "d", "Particle age since seeding or last division" PARAM_FILLVALUE);
   check_nc_error(nc_def_var(ncid, "particle_birth_day", NC_REALTYPE, 2, dims, &birth_id));
   set_nc_attributes(ncid, birth_id, "d", "Model day the particle was created" PARAM_FILLVALUE);
   check_nc_error(nc_def_var(ncid, "particle_cum_par", NC_REALTYPE, 2, dims, &cumpar_id));
   set_nc_attributes(ncid, cumpar_id, "W/m2 d", "Cumulative light dose (integral of PAR dt)" PARAM_FILLVALUE);
   check_nc_error(nc_def_var(ncid, "particle_cum_par_mol", NC_REALTYPE, 2, dims, &cumpar_mol_id));
   set_nc_attributes(ncid, cumpar_mol_id, "mol photons/m2", "Cumulative light dose, mol photon equivalent (WtouE=4.57)" PARAM_FILLVALUE);
   check_nc_error(nc_def_var(ncid, "particle_cum_c_fixed", NC_REALTYPE, 2, dims, &cumc_id));
   set_nc_attributes(ncid, cumc_id, "pmol C/particle", "Lifetime gross carbon fixed by the super-individual" PARAM_FILLVALUE);
   check_nc_error(nc_def_var(ncid, "particle_ndiv", NC_REALTYPE, 2, dims, &ndiv_id));
   set_nc_attributes(ncid, ndiv_id, "number", "Number of cell divisions undergone" PARAM_FILLVALUE);
   check_nc_error(nc_def_var(ncid, "particle_npq", NC_REALTYPE, 2, dims, &npq_id));
   set_nc_attributes(ncid, npq_id, "%", "Non-photochemical quenching state; percent effective PAR reduction" PARAM_FILLVALUE);
   check_nc_error(nc_def_var(ncid, "particle_cldose", NC_REALTYPE, 2, dims, &cldose_id));
   /* Defined ONLY when the diagnostic is registered. Defining them unconditionally left a
    * disabled run publishing four variables reading 0 everywhere, and 0 J/m2 reads as a
    * measurement rather than as "not computed". */
   if (ptm_idx_dose_scal_l >= 1) {
      check_nc_error(nc_def_var(ncid, "particle_dose_scalar_lastday", NC_REALTYPE, 2, dims, &dose_scal_l_id));
      set_nc_attributes(ncid, dose_scal_l_id, "J/m2", "Raw scalar PAR exposure accumulated over the last COMPLETE reporting day" PARAM_FILLVALUE);
      check_nc_error(nc_def_var(ncid, "particle_dose_effective_lastday", NC_REALTYPE, 2, dims, &dose_eff_l_id));
      set_nc_attributes(ncid, dose_eff_l_id, "J/m2", "Growth-effective PAR exposure accumulated over the last COMPLETE reporting day" PARAM_FILLVALUE);
      check_nc_error(nc_def_var(ncid, "particle_prodC_lastday", NC_REALTYPE, 2, dims, &dprodC_l_id));
      set_nc_attributes(ncid, dprodC_l_id, "mg C/cell", "Carbon fixed per cell accumulated over the last COMPLETE reporting day (signed)" PARAM_FILLVALUE);
      check_nc_error(nc_def_var(ncid, "particle_dose_dayindex", NC_REALTYPE, 2, dims, &dday_id));
      set_nc_attributes(ncid, dday_id, "yearday", "Reporting day the particle is currently accumulating into" PARAM_FILLVALUE);
   }

   set_nc_attributes(ncid, cldose_id, "umol/m2/s", "Ranjbar cumulative light dose CL; drains at NPQ_PRR and resets NPQ at zero (NPQ_mode=1 only)" PARAM_FILLVALUE);
   check_nc_error(nc_def_var(ncid, "particle_fcol", NC_REALTYPE, 2, dims, &fcol_id));
   set_nc_attributes(ncid, fcol_id, "-", "Effective colony size as a multiple of cell ESD; prognostic when simColony=1 (Water Research 2022 Eq 10-11)" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_status", NC_INT, 2, dims, &stat_id));
   check_nc_error(nc_put_att_text(ncid, stat_id, "long_name", strlen("Status of Particle"), "Status of Particle"));

   check_nc_error(nc_def_var(ncid, "particle_flag", NC_INT, 2, dims, &flag_id));
   check_nc_error(nc_put_att_text(ncid, flag_id, "long_name", strlen("Location Flag of Particle"), "Location Flag of Particle"));

   check_nc_error(nc_def_var(ncid, "particle_ptid", NC_INT, 2, dims, &ptid_id));
   check_nc_error(nc_put_att_text(ncid, ptid_id, "long_name", strlen("ID of Particle"), "ID of Particle"));

   check_nc_error(nc_def_var(ncid, "particle_group_id", NC_INT, 2, dims, &group_id));
   check_nc_error(nc_put_att_text(ncid, group_id, "long_name", strlen("Configured particle group ID"), "Configured particle group ID"));
   ptm_group_names_attr(groups_attr, sizeof(groups_attr));
   check_nc_error(nc_put_att_text(ncid, group_id, "group_names", strlen(groups_attr), groups_attr));

   check_nc_error(nc_def_var(ncid, "particle_pigment_id", NC_INT, 2, dims, &pigment_id));
   check_nc_error(nc_put_att_text(ncid, pigment_id, "long_name", strlen("Configured particle pigment ID"), "Configured particle pigment ID"));
   check_nc_error(nc_put_att_text(ncid, pigment_id, "pigment_id_legend", strlen(pigment_legend), pigment_legend));

   check_nc_error(nc_def_var(ncid, "particle_phyto_link_id", NC_INT, 2, dims, &phyto_link_id));
   check_nc_error(nc_put_att_text(ncid, phyto_link_id, "long_name", strlen("Configured AED phytoplankton link"), "Configured AED phytoplankton link"));

   check_nc_error(nc_def_var(ncid, "particle_layer", NC_INT, 2, dims, &particle_layer_id));
   check_nc_error(nc_put_att_text(ncid, particle_layer_id, "long_name", strlen("Zero-based GLM layer of particle"), "Zero-based GLM layer of particle"));

   check_nc_error(nc_def_var(ncid, "particle_oasim_status", NC_INT, 2, dims, &oasim_status_id));
   check_nc_error(nc_put_att_text(ncid, oasim_status_id, "long_name", strlen("OASIM particle spectral sample status"), "OASIM particle spectral sample status"));
   check_nc_error(nc_put_att_text(ncid, oasim_status_id, "legend", strlen("0=exact,1=invalid_layer,2=missing,3=inactive"), "0=exact,1=invalid_layer,2=missing,3=inactive"));

   check_nc_error(nc_def_var(ncid, "particle_vertical_mode", NC_INT, 2, dims, &vertical_mode_id));
   check_nc_error(nc_put_att_text(ncid, vertical_mode_id, "long_name", strlen("Interpreted vertical behavior"), "Interpreted vertical behavior"));
   check_nc_error(nc_put_att_text(ncid, vertical_mode_id, "legend", strlen("1=upward_buoyant,2=neutral,3=sinking"), "1=upward_buoyant,2=neutral,3=sinking"));

   check_nc_error(nc_def_var(ncid, "particle_behavior_id", NC_INT, 2, dims, &behavior_id));
   check_nc_error(nc_put_att_text(ncid, behavior_id, "long_name", strlen("Interpreted particle behavior ID"), "Interpreted particle behavior ID"));
   check_nc_error(nc_put_att_text(ncid, behavior_id, "legend", strlen("1=cyano_buoyant,2=green_neutral,3=diatom_sinking,4=planktothrix_deep_light_adapted"), "1=cyano_buoyant,2=green_neutral,3=diatom_sinking,4=planktothrix_deep_light_adapted"));

   check_nc_error(nc_def_var(ncid, "particle_light_440", NC_REALTYPE, 2, dims, &light_440_id));
   set_nc_attributes(ncid, light_440_id, "OASIM_band", "Particle OASIM direct+diffuse band irradiance at 440 nm" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_light_490", NC_REALTYPE, 2, dims, &light_490_id));
   set_nc_attributes(ncid, light_490_id, "OASIM_band", "Particle OASIM direct+diffuse band irradiance at 490 nm" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_light_550", NC_REALTYPE, 2, dims, &light_550_id));
   set_nc_attributes(ncid, light_550_id, "OASIM_band", "Particle OASIM direct+diffuse band irradiance at 550 nm" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_light_565", NC_REALTYPE, 2, dims, &light_565_id));
   set_nc_attributes(ncid, light_565_id, "OASIM_band", "Particle OASIM direct+diffuse band irradiance at 565 nm" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_light_590", NC_REALTYPE, 2, dims, &light_590_id));
   set_nc_attributes(ncid, light_590_id, "OASIM_band", "Particle OASIM direct+diffuse band irradiance at 590 nm" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_light_620", NC_REALTYPE, 2, dims, &light_620_id));
   set_nc_attributes(ncid, light_620_id, "OASIM_band", "Particle OASIM direct+diffuse band irradiance at 620 nm" PARAM_FILLVALUE);

   check_nc_error(nc_def_var(ncid, "particle_light_665", NC_REALTYPE, 2, dims, &light_665_id));
   set_nc_attributes(ncid, light_665_id, "OASIM_band", "Particle OASIM direct+diffuse band irradiance at 665 nm" PARAM_FILLVALUE);

    check_nc_error(nc_def_var(ncid, "particle_spectral_score", NC_REALTYPE, 2, dims, &spectral_score_id));
    set_nc_attributes(ncid, spectral_score_id, "OASIM_band", "Pigment-weighted particle spectral suitability from OASIM bands" PARAM_FILLVALUE);

    check_nc_error(nc_def_var(ncid, "particle_light_preference_nm", NC_REALTYPE, 2, dims, &light_preference_nm_id));
    set_nc_attributes(ncid, light_preference_nm_id, "nm", "Dominant pigment wavelength preference for this particle" PARAM_FILLVALUE);

    check_nc_error(nc_def_var(ncid, "particle_best_light", NC_REALTYPE, 2, dims, &best_light_id));
    set_nc_attributes(ncid, best_light_id, "OASIM_band", "Best available OASIM light for this particle pigment preference" PARAM_FILLVALUE);

    check_nc_error(nc_def_var(ncid, "particle_par_scalar_depth", NC_REALTYPE, 2, dims, &par_scalar_depth_id));
    if (buoy_any) {
        /* Wallace & Hamilton internal state. Defined ONLY when the model is active, so a
         * buoy_model = 0 run produces a structurally identical file and stays bit-comparable. */
        check_nc_error(nc_def_var(ncid, "particle_buoy_rate", NC_REALTYPE, 2, dims, &buoy_rate_id));
        set_nc_attributes(ncid, buoy_rate_id, "kg/m3/day", "W&H lagged rate of particle density change" PARAM_FILLVALUE);
        check_nc_error(nc_def_var(ncid, "particle_lightdose", NC_REALTYPE, 2, dims, &lightdose_id));
        set_nc_attributes(ncid, lightdose_id, "W/m2", "W&H light-dose EMA driving dark de-ballasting" PARAM_FILLVALUE);
        check_nc_error(nc_def_var(ncid, "particle_s_act", NC_REALTYPE, 2, dims, &s_act_id));
        set_nc_attributes(ncid, s_act_id, "-", "dormancy activity applied to the buoyancy flux" PARAM_FILLVALUE);
        check_nc_error(nc_def_var(ncid, "particle_buoy_I", NC_REALTYPE, 2, dims, &buoy_I_id));
        set_nc_attributes(ncid, buoy_I_id, "W/m2", "irradiance the buoyancy integrator actually used" PARAM_FILLVALUE);
        check_nc_error(nc_def_var(ncid, "particle_buoy_stat", NC_REALTYPE, 2, dims, &buoy_stat_id));
        set_nc_attributes(ncid, buoy_stat_id, "-", "buoyancy light source: 0 OASIM, 1 Lake Light fallback" PARAM_FILLVALUE);
    }
    set_nc_attributes(ncid, par_scalar_depth_id, "W/m2", "Depth-resolved scalar PAR at the particle location" PARAM_FILLVALUE);
    {   /* Phase A: depth-attenuated per-band irradiance + absorbed-equivalent PAR */
        char vn[64], vd[160]; int b;
        for (b = 0; b < PTM_OASIM_NBANDS; b++) {
            snprintf(vn, sizeof(vn), "particle_par_depth_%d", ptm_band_nm[b]);
            snprintf(vd, sizeof(vd), "Depth-attenuated band irradiance at %d nm "
                     "(as used by biology, not resampled at output time)", ptm_band_nm[b]);
            check_nc_error(nc_def_var(ncid, vn, NC_REALTYPE, 2, dims, &par_depth_band_id[b]));
            set_nc_attributes(ncid, par_depth_band_id[b], "W/m2/nm", vd PARAM_FILLVALUE);
        }
        check_nc_error(nc_def_var(ncid, "particle_par_absorbed", NC_REALTYPE, 2, dims, &par_absorbed_id));
        set_nc_attributes(ncid, par_absorbed_id, "W/m2",
            "Group-specific absorbed-equivalent PAR: scalar PAR at particle depth times "
            "OASIM's spectral absorption factor OAS_afac_iop<group>" PARAM_FILLVALUE);
    }

    check_nc_error(nc_def_var(ncid, "particle_par_layer_center", NC_REALTYPE, 2, dims, &par_layer_center_id));
    set_nc_attributes(ncid, par_layer_center_id, "W/m2", "Layer-centre scalar PAR from OASIM" PARAM_FILLVALUE);

    check_nc_error(nc_def_var(ncid, "particle_par_depth_factor", NC_REALTYPE, 2, dims, &par_depth_factor_id));
    set_nc_attributes(ncid, par_depth_factor_id, "1", "Exponential depth factor relative to layer centre" PARAM_FILLVALUE);

    check_nc_error(nc_def_var(ncid, "particle_par_scalar_status", NC_INT, 2, dims, &par_scalar_status_id));
    check_nc_error(nc_put_att_text(ncid, par_scalar_status_id, "long_name", strlen("Scalar PAR source flag"), "Scalar PAR source flag"));
    check_nc_error(nc_put_att_text(ncid, par_scalar_status_id, "legend", strlen("0=OASIM,1=fallback"), "0=OASIM,1=fallback"));

   check_nc_error(nc_def_var(ncid, "particle_grp", NC_INT, 2, dims, &grp_id));
   check_nc_error(nc_put_att_text(ncid, grp_id, "long_name", strlen("Phyto group/species of Particle"), "Phyto group/species of Particle"));

   define_mode_off(&ncid);
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
