/******************************************************************************
 * glm_lib.c - Library API for GLM (callable from Python, etc.)
 *
 * Provides glm_run_with_nml() so GLM can be invoked as a library without
 * the main() entry point. Used by pglm Python package.
 ******************************************************************************/
#include <string.h>
#include "glm.h"
#include "glm_globals.h"   /* CLOGICAL non_avg */

extern char glm_nml_file[];
extern void run_model(void);
/* libutil: resets the input CSV pool (_n_inf) so each library run starts fresh */
extern void close_all_csv_inputs(void);
void glm_finish(void);

/* Stub for glm_output.c (defined in glm_main.c for executable build) */
char *all_plots_name = NULL;

#ifdef _WIN32
__declspec(dllexport)
#endif
/* ---------------------------------------------------------------------------
 * BMI-style split (clarena): init / run / finish as separate entry points.
 *   glm_init_with_nml(nml) -> init_model()      (reads nml, allocates, WQ init)
 *   glm_run()              -> do_model[_non_avg] (time loop)
 *   glm_finish()           -> end_model()       (restart, close files, free WQ/zone arrays)
 * glm_run_with_nml() keeps its behaviour (= init + run + finish) so existing
 * callers are unaffected. A finish that frees everything is what allows the
 * library to be re-initialised in the same process for the next coupling step.
 * ------------------------------------------------------------------------- */
extern void init_model(int *jstart, int *nsave);
extern void do_model(int jstart, int nsave);
extern void do_model_non_avg(int jstart, int nsave);
extern void end_model(void);
static int _lib_jstart = 0, _lib_nsave = 0, _lib_initialised = 0;

#ifdef _WIN32
__declspec(dllexport)
#endif
void glm_init_with_nml(const char *nml_path)
{
    if (_lib_initialised) glm_finish();   /* never leak a previous run */
    if (nml_path && nml_path[0]) {
        strncpy(glm_nml_file, nml_path, 255);
        glm_nml_file[255] = '\0';
    }
    init_model(&_lib_jstart, &_lib_nsave);
    _lib_initialised = 1;
}

#ifdef _WIN32
__declspec(dllexport)
#endif
void glm_run(void)
{
    if (!_lib_initialised) return;
    if (non_avg) do_model_non_avg(_lib_jstart, _lib_nsave);
    else         do_model(_lib_jstart, _lib_nsave);
}

#ifdef _WIN32
__declspec(dllexport)
#endif
void glm_finish(void)
{
    if (!_lib_initialised) return;
    end_model();               /* writes restart.nc, then frees WQ + zone arrays */
    close_all_csv_inputs();    /* reset the input CSV pool (see note below) */
    _lib_initialised = 0;
}

#ifdef _WIN32
__declspec(dllexport)
#endif
void glm_run_with_nml(const char *nml_path)
{
    if (nml_path && nml_path[0]) {
        strncpy(glm_nml_file, nml_path, 255);
        glm_nml_file[255] = '\0';
    }
    run_model();

    /* Each library invocation is an independent run (clarena calls this once per
     * coupling step). run_model()'s end_model() closes the inflow/outflow/met
     * files, but close_csv_input() only reclaims a pool slot on LIFO closes, so
     * _n_inf never returns to 0 and ratchets up ~(num_inflows+num_outlets) per
     * step -> "Too many csv_files open" after a few steps. Reset the pool here. */
    close_all_csv_inputs();
}
