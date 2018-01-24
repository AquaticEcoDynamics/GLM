/******************************************************************************
 *                                                                            *
 * glm_init.c                                                                 *
 *                                                                            *
 * Developed by :                                                             *
 *     AquaticEcoDynamics (AED) Group                                         *
 *     School of Earth & Environment                                          *
 *     University of Western Australia                                        *
 *                                                                            *
 *     http://aed.see.uwa.edu.au/                                             *
 *                                                                            *
 * Copyright 2013 - 2016 -  The University of Western Australia               *
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

#include "glm.h"

#include "glm_types.h"
#include "glm_const.h"
#include "glm_globals.h"
#include "glm_csv.h"
#include "glm_input.h"
#include "glm_util.h"
#include "glm_layers.h"
#include "glm_wqual.h"
#include "glm_lnum.h"
#include "glm_bird.h"
#include "glm_ncdf.h"

#include <aed_time.h>
#include <namelist.h>

//#define dbgprt(...) fprintf(stderr, __VA_ARGS__)
#define dbgprt(...) /* __VA_ARGS__ */


extern int *WQ_VarsIdx;

static AED_REAL   base_elev;
static AED_REAL   crest_elev;
extern LOGICAL    seepage;
extern AED_REAL   seepage_rate;

char glm_nml_file[256] = "glm2.nml";
char wq_lib[256] = "aed2";

static void create_lake(int namlst);
static void initialise_lake(int namlst);
static int init_time(char *start, char *stop, int timefmt, int *nDays);

/*############################################################################*/

/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void init_glm(int *jstart, char *outp_dir, char *outp_fn, int *nsave)
{
    int jyear, jmonth, jday, julianday;

    /*---------------------------------------------
     * glm setup
     *-------------------------------------------*/
    char           *sim_name = NULL;
    int             max_layers;
    AED_REAL        min_layer_vol;
    AED_REAL        min_layer_thick;
    AED_REAL        max_layer_thick;
//  AED_REAL        Kw;
    extern AED_REAL Benthic_Imin;
//  AED_REAL        coef_mix_conv;
//  AED_REAL        coef_mix_eta;
//  AED_REAL        coef_mix_ct;
//  AED_REAL        coef_mix_cs;
//  AED_REAL        coef_mix_kh;
//  AED_REAL        coef_mix_hyp;
//  CLOGICAL        non_avg;
//  int             deep_mixing;
    /*-------------------------------------------*/

    /*---------------------------------------------
     * glm restart
     *-------------------------------------------*/
    char           *restart_file = NULL;
    /*-------------------------------------------*/

    /*---------------------------------------------
     * wq setup
     *-------------------------------------------*/
    char           *twq_lib = NULL;
    char           *wq_nml_file = "aed2.nml";
    int             lode_method;
    int             lsplit_factor;
//  LOGICAL         bioshade_feedback;
//  LOGICAL         repair_state;
//  CLOGICAL        mobility_off;
//  int             benthic_mode;
//  int             n_zones;
//  AED_REAL       *zone_heights = NULL;
    /*-------------------------------------------*/

    /*---------------------------------------------
     * time format
     *-------------------------------------------*/
    int             timefmt;
    char           *start = NULL;
    char           *stop  = NULL;
    AED_REAL        dt;        // timestep
    int             num_days;  // number of days to run the sim
//  AED_REAL        timezone_r;
    /*-------------------------------------------*/

    /*---------------------------------------------
     * output
     *-------------------------------------------*/
    char           *out_dir = NULL;
    char           *out_fn  = NULL;
//  LOGICAL         out_lkn;
//  int             nsave;
    int             csv_point_nlevs  = 0;
    LOGICAL        *csv_point_frombot = NULL;
    char           *csv_point_fname  = NULL;
    AED_REAL       *csv_point_at     = NULL;
    int             csv_point_nvars  = 0;
    char          **csv_point_vars   = NULL;
    char           *csv_lake_fname   = NULL;
    LOGICAL         csv_outlet_allinone = FALSE;
    char           *csv_outlet_fname = NULL;
    int             csv_outlet_nvars = 0;
    char          **csv_outlet_vars  = NULL;
    char           *csv_ovrflw_fname = NULL;
    /*-------------------------------------------*/

    /*---------------------------------------------
     * meteorology
     *-------------------------------------------*/
    LOGICAL         met_sw;          // Include surface meteorological forcing
    char           *lw_type = NULL;  // Type LW measurement (LW_IN/LW_CC/LW_NET)
    LOGICAL         rain_sw;         // Rainfall composition
    LOGICAL         snow_sw;         // Snowfall
    char           *meteo_fl = NULL; // Name of meteorology input file
//  int             lw_ind;          // type of longwave radiation - now in glm_input
//  LOGICAL         atm_stab;        // Account for non-neutral atmospheric stability
//  LOGICAL         subdaily;        //
//  AED_REAL        CD;
//  AED_REAL        CE;
//  AED_REAL        CH;
    extern AED_REAL wind_factor;
    extern AED_REAL sw_factor;
    extern AED_REAL lw_factor;
    extern AED_REAL at_factor;
    extern AED_REAL rh_factor;
    extern AED_REAL rain_factor;
    extern int      rad_mode;
    extern int      albedo_mode;
    extern int      cloud_mode;
    char           *timefmt_m = NULL;
    extern AED_REAL timezone_m;
    /*-------------------------------------------*/

    /*---------------------------------------------
     * inflow
     *-------------------------------------------*/
    int             num_inflows;
    LOGICAL        *subm_flag      = NULL;
    char          **names_of_strms = NULL;
    AED_REAL       *strm_hf_angle  = NULL;
    AED_REAL       *strmbd_slope   = NULL;
    AED_REAL       *strmbd_drag    = NULL;
    AED_REAL       *inflow_factor  = NULL;
    char          **inflow_fl      = NULL;
    int             inflow_varnum;
    char          **inflow_vars    = NULL;
    AED_REAL        coef_inf_entrain;
    char           *timefmt_i      = NULL;
    extern AED_REAL timezone_i;
    /*-------------------------------------------*/

    /*---------------------------------------------
     * outflow
     *-------------------------------------------*/
    int             num_outlet;
    LOGICAL        *flt_off_sw   = NULL;
    int            *outlet_type    = NULL;
    int             crit_O2        = -1;
    int             crit_O2_dep    = -1;
    int             crit_O2_days   = -1;
    AED_REAL       *outlet_crit    = NULL;
    char          **O2name         = NULL;
    int             O2idx          = 0;
    AED_REAL       *target_temp    = NULL;
    AED_REAL        min_lake_temp;
    LOGICAL         mix_withdraw;
    LOGICAL         coupl_oxy_sw;
    extern AED_REAL fac_range_upper;
    extern AED_REAL fac_range_lower;
    AED_REAL       *outl_elvs    = NULL;
    AED_REAL       *bsn_len_outl = NULL;
    AED_REAL       *bsn_wid_outl = NULL;
    char          **outflow_fl   = NULL;
    char           *withdrTemp_fl  = NULL;
    AED_REAL       *outflow_factor = NULL;
    char           *timefmt_o    = NULL;
    extern AED_REAL timezone_o;
    /*-------------------------------------------*/

    /*---------------------------------------------
     * snowice
     *-------------------------------------------*/
    extern AED_REAL         snow_albedo_factor;
    extern AED_REAL         snow_rho_max;
    extern AED_REAL         snow_rho_min;
    /*-------------------------------------------*/

    /*---------------------------------------------
     * fetch
     *-------------------------------------------*/
    extern LOGICAL     fetch_sw;
    extern int         fetch_ndirs;
    extern AED_REAL   *fetch_dirs;
    extern AED_REAL   *fetch_scale;
    extern AED_REAL    fetch_height;
    extern AED_REAL    fetch_porosity;
    /*-------------------------------------------*/

    /*---------------------------------------------
     * sed_heat
     *-------------------------------------------*/
    extern CLOGICAL         sed_heat_sw;
    extern AED_REAL         sed_heat_Ksoil;
    extern AED_REAL         sed_temp_depth;
    AED_REAL               *sed_temp_mean       = NULL;
    AED_REAL               *sed_temp_amplitude  = NULL;
    AED_REAL               *sed_temp_peak_doy   = NULL;
    //extern AED_REAL         sed_temp_amplitude;
    //extern AED_REAL         sed_temp_peak_doy;
    /*-------------------------------------------*/

    int i, j, k;
    int namlst;

    //==========================================================================
    NAMELIST glm_setup[] = {
          { "glm_setup",         TYPE_START,            NULL               },
          { "sim_name",          TYPE_STR,              &sim_name          },
          { "max_layers",        TYPE_INT,              &max_layers        },
          { "min_layer_vol",     TYPE_DOUBLE,           &min_layer_vol     },
          { "min_layer_thick",   TYPE_DOUBLE,           &min_layer_thick   },
          { "max_layer_thick",   TYPE_DOUBLE,           &max_layer_thick   },
          { "Kw",                TYPE_DOUBLE,           &Kw                },
          { "Benthic_Imin",      TYPE_DOUBLE,           &Benthic_Imin      },
          { "coef_mix_conv",     TYPE_DOUBLE,           &coef_mix_conv     },
          { "coef_wind_stir",    TYPE_DOUBLE,           &coef_wind_stir    },
          { "coef_mix_turb",     TYPE_DOUBLE,           &coef_mix_turb     },
          { "coef_mix_shear",    TYPE_DOUBLE,           &coef_mix_shear    },
          { "coef_mix_KH",       TYPE_DOUBLE,           &coef_mix_KH       },
          { "coef_mix_hyp",      TYPE_DOUBLE,           &coef_mix_hyp      },
          { "non_avg",           TYPE_BOOL,             &non_avg           },
          { "deep_mixing",       TYPE_INT,              &deep_mixing       },
          { NULL,                TYPE_END,              NULL               }
    };
    NAMELIST glm_restart[] = {
          { "glm_restart",       TYPE_START,            NULL               },
          { "restart_file",      TYPE_STR,              &restart_file      },
          { NULL,                TYPE_END,              NULL               }
    };
    NAMELIST wq_setup[] = {
          { "wq_setup",          TYPE_START,            NULL               },
          { "wq_lib",            TYPE_STR,              &twq_lib           },
          { "wq_nml_file",       TYPE_STR,              &wq_nml_file       },
          { "ode_method",        TYPE_INT,              &lode_method       },
          { "split_factor",      TYPE_INT,              &lsplit_factor     },
          { "bioshade_feedback", TYPE_BOOL,             &bioshade_feedback },
          { "repair_state",      TYPE_BOOL,             &repair_state      },
          { "mobility_off",      TYPE_BOOL,             &mobility_off      },
          { "benthic_mode",      TYPE_INT,              &benthic_mode      },
          { "n_zones",           TYPE_INT,              &n_zones           },
          { "zone_heights",      TYPE_DOUBLE|MASK_LIST, &zone_heights      },
          { NULL,                TYPE_END,              NULL               }
    };
    NAMELIST time[] = {
          { "time",              TYPE_START,            NULL               },
          { "timefmt",           TYPE_INT,              &timefmt           },
          { "start",             TYPE_STR,              &start             },
          { "stop",              TYPE_STR,              &stop              },
          { "dt",                TYPE_DOUBLE,           &dt                },
          { "num_days",          TYPE_INT,              &num_days          },
          { "timezone",          TYPE_DOUBLE,           &timezone_r        },
          { NULL,                TYPE_END,              NULL               }
    };
    NAMELIST output[] = {
          { "output",            TYPE_START,            NULL               },
          { "out_dir",           TYPE_STR,              &out_dir           },
          { "out_fn",            TYPE_STR,              &out_fn            },
          { "nsave",             TYPE_INT,               nsave             },
          { "csv_point_nlevs",   TYPE_INT,              &csv_point_nlevs   },
          { "csv_point_fname",   TYPE_STR,              &csv_point_fname   },
          { "csv_point_frombot", TYPE_BOOL|MASK_LIST,   &csv_point_frombot },
          { "csv_point_at",      TYPE_DOUBLE|MASK_LIST, &csv_point_at      },
          { "csv_point_nvars",   TYPE_INT,              &csv_point_nvars   },
          { "csv_point_vars",    TYPE_STR|MASK_LIST,    &csv_point_vars    },
          { "csv_lake_fname",    TYPE_STR,              &csv_lake_fname    },
          { "csv_outlet_allinone", TYPE_BOOL,           &csv_outlet_allinone},
          { "csv_outlet_fname",  TYPE_STR,              &csv_outlet_fname  },
          { "csv_outlet_nvars",  TYPE_INT,              &csv_outlet_nvars  },
          { "csv_outlet_vars",   TYPE_STR|MASK_LIST,    &csv_outlet_vars   },
          { "csv_ovrflw_fname",  TYPE_STR,              &csv_ovrflw_fname  },
          { NULL,                TYPE_END,              NULL               }
    };
    NAMELIST meteorology[] = {
          { "meteorology",       TYPE_START,            NULL               },
          { "met_sw",            TYPE_BOOL,             &met_sw            },
          { "lw_type",           TYPE_STR,              &lw_type           },
          { "rain_sw",           TYPE_BOOL,             &rain_sw           },
     //   { "snow_sw",           TYPE_BOOL,             &snow_sw           },
          { "meteo_fl",          TYPE_STR,              &meteo_fl          },
          { "subdaily",          TYPE_BOOL,             &subdaily          },
          { "atm_stab",          TYPE_BOOL,             &atm_stab          },
          { "rad_mode",          TYPE_INT,              &rad_mode          },
          { "albedo_mode",       TYPE_INT,              &albedo_mode       },
          { "cloud_mode",        TYPE_INT,              &cloud_mode        },
          { "wind_factor",       TYPE_DOUBLE,           &wind_factor       },
          { "sw_factor",         TYPE_DOUBLE,           &sw_factor         },
          { "lw_factor",         TYPE_DOUBLE,           &lw_factor         },
          { "at_factor",         TYPE_DOUBLE,           &at_factor         },
          { "rh_factor",         TYPE_DOUBLE,           &rh_factor         },
          { "rain_factor",       TYPE_DOUBLE,           &rain_factor       },
          { "CD",                TYPE_DOUBLE,           &CD                },
          { "CE",                TYPE_DOUBLE,           &CE                },
          { "CH",                TYPE_DOUBLE,           &CH                },
          { "catchrain",         TYPE_BOOL,             &catchrain         },
          { "rain_threshold",    TYPE_DOUBLE,           &rain_threshold    },
          { "runoff_coef",       TYPE_DOUBLE,           &runoff_coef       },
          { "time_fmt",          TYPE_STR,              &timefmt_m         },
          { "timezone",          TYPE_DOUBLE,           &timezone_m        },
          { NULL,                TYPE_END,              NULL               }
    };
    NAMELIST inflow[] = {
          { "inflow",            TYPE_START,            NULL               },
          { "num_inflows",       TYPE_INT,              &num_inflows       },
          { "subm_flag",         TYPE_BOOL|MASK_LIST,   &subm_flag         },
          { "names_of_strms",    TYPE_STR|MASK_LIST,    &names_of_strms    },
          { "strm_hf_angle",     TYPE_DOUBLE|MASK_LIST, &strm_hf_angle     },
          { "strmbd_slope",      TYPE_DOUBLE|MASK_LIST, &strmbd_slope      },
          { "strmbd_drag",       TYPE_DOUBLE|MASK_LIST, &strmbd_drag       },
          { "inflow_factor",     TYPE_DOUBLE|MASK_LIST, &inflow_factor     },
          { "inflow_fl",         TYPE_STR|MASK_LIST,    &inflow_fl         },
          { "inflow_varnum",     TYPE_INT,              &inflow_varnum     },
          { "inflow_vars",       TYPE_STR|MASK_LIST,    &inflow_vars       },
          { "coef_inf_entrain",  TYPE_DOUBLE,           &coef_inf_entrain  },
          { "time_fmt",          TYPE_STR,              &timefmt_i         },
          { "timezone",          TYPE_DOUBLE,           &timezone_i        },
          { NULL,                TYPE_END,              NULL               }
    };
    NAMELIST outflow[] = {
          { "outflow",           TYPE_START,            NULL               },
          { "num_outlet",        TYPE_INT,              &num_outlet        },
          { "outlet_type",       TYPE_INT|MASK_LIST,    &outlet_type       },
          { "crit_O2",           TYPE_INT,              &crit_O2           },
          { "crit_O2_dep",       TYPE_INT,              &crit_O2_dep       },
          { "crit_O2_days",      TYPE_INT,              &crit_O2_days      },
          { "outlet_crit",       TYPE_DOUBLE|MASK_LIST, &outlet_crit       },
          { "O2name",            TYPE_STR|MASK_LIST,    &O2name            },
          { "O2idx",             TYPE_INT,              &O2idx             },
          { "target_temp",       TYPE_DOUBLE|MASK_LIST, &target_temp       },
          { "min_lake_temp",     TYPE_DOUBLE,           &min_lake_temp     },
          { "fac_range_upper",   TYPE_DOUBLE,           &fac_range_upper   },
          { "fac_range_lower",   TYPE_DOUBLE,           &fac_range_lower   },
          { "mix_withdraw",      TYPE_BOOL,             &mix_withdraw      },
          { "coupl_oxy_sw",      TYPE_BOOL,             &coupl_oxy_sw      },
          { "flt_off_sw",        TYPE_BOOL|MASK_LIST,   &flt_off_sw        },
          { "outl_elvs",         TYPE_DOUBLE|MASK_LIST, &outl_elvs         },
          { "bsn_len_outl",      TYPE_DOUBLE|MASK_LIST, &bsn_len_outl      },
          { "bsn_wid_outl",      TYPE_DOUBLE|MASK_LIST, &bsn_wid_outl      },
          { "outflow_fl",        TYPE_STR|MASK_LIST,    &outflow_fl        },
          { "withdrTemp_fl",     TYPE_STR,              &withdrTemp_fl     },
          { "outflow_factor",    TYPE_DOUBLE|MASK_LIST, &outflow_factor    },
          { "seepage",           TYPE_BOOL,             &seepage           },
          { "seepage_rate",      TYPE_DOUBLE,           &seepage_rate      },
          { "time_fmt",          TYPE_STR,              &timefmt_o         },
          { "timezone",          TYPE_DOUBLE,           &timezone_o        },
          { NULL,                TYPE_END,              NULL               }
    };
    NAMELIST diffuser[] = {
          { "diffuser",          TYPE_START,            NULL               },
          { "NumDif",            TYPE_INT,              &NumDif            },
          { "diff",              TYPE_DOUBLE|MASK_LIST, &mol_diffusivity   },
          { NULL,                TYPE_END,              NULL               }
    };
    NAMELIST debugging[] = {
          { "debugging",         TYPE_START,            NULL               },
          { "disable_evap",      TYPE_BOOL,             &no_evap           },
          { NULL,                TYPE_END,              NULL               }
    };
    NAMELIST snowice[] = {
          { "snowice",           TYPE_START,            NULL               },
          { "snow_albedo_factor",TYPE_DOUBLE,           &snow_albedo_factor},
          { "snow_rho_max",      TYPE_DOUBLE,           &snow_rho_max      },
          { "snow_rho_min",      TYPE_DOUBLE,           &snow_rho_min      },
          { NULL,                TYPE_END,              NULL               }
    };
    NAMELIST fetch[] = {
          { "fetch",             TYPE_START,            NULL               },
          { "fetch_sw",          TYPE_BOOL,             &fetch_sw          },
          { "num_dir",           TYPE_INT,              &fetch_ndirs       },
          { "wind_dir",          TYPE_DOUBLE|MASK_LIST, &fetch_dirs        },
          { "fetch_scale",       TYPE_DOUBLE|MASK_LIST, &fetch_scale       },
          { "edge_height",       TYPE_DOUBLE,           &fetch_height      },
          { "edge_porosity",     TYPE_DOUBLE,           &fetch_porosity    },
          { NULL,                TYPE_END,              NULL               }
    };
    NAMELIST sed_heat[] = {
          { "sed_heat",          TYPE_START,            NULL               },
          { "sed_heat_Ksoil",    TYPE_DOUBLE,           &sed_heat_Ksoil    },
          { "sed_temp_depth",    TYPE_DOUBLE,           &sed_temp_depth    },
          { "sed_temp_mean",     TYPE_DOUBLE|MASK_LIST, &sed_temp_mean     },
          { "sed_temp_amplitude",TYPE_DOUBLE|MASK_LIST, &sed_temp_amplitude},
          { "sed_temp_peak_doy", TYPE_DOUBLE|MASK_LIST, &sed_temp_peak_doy },
          { NULL,                TYPE_END,              NULL               }
    };
/*----------------------------------------------------------------------------*/

    fprintf(stderr, "Reading config from %s\n",glm_nml_file);

    //-------------------------------------------------
    // Open the namelist file.
    if ( (namlst = open_namelist(glm_nml_file)) < 0 ) {
        fprintf(stderr,"Error opening namelist file %s\n", glm_nml_file);
        exit(1);
    }

    //-------------------------------------------------
    // Set some default values
    coef_inf_entrain = 0.;
    Kw = 0.2;

    if ( get_namelist(namlst, glm_setup) != 0 ) {
       fprintf(stderr,"Error reading the 'glm_setup' namelist from %s\n", glm_nml_file);
       exit(1);
    }

    MaxLayers = max_layers;
    VMin = min_layer_vol;
    DMin = min_layer_thick;
    DMax = max_layer_thick;
    NumLayers = 0;
    n_zones = 0;

    //-------------------------------------------------
    if ( get_namelist(namlst, glm_restart) ) {
        do_restart = FALSE;
    } else {
        do_restart = (restart_file != NULL);
    }

    wq_calc   = TRUE;

    if ( get_namelist(namlst, wq_setup) ) {
        fprintf(stderr, "No WQ config\n");
        twq_lib           = "aed2";
        wq_calc           = FALSE;
        ode_method        = 1;
        split_factor      = 1;
        bioshade_feedback = TRUE;
        repair_state      = FALSE;
        n_zones           = 0;
    } else {
        ode_method        = lode_method;
        split_factor      = lsplit_factor;
    }
    if ( twq_lib != NULL ) strncpy(wq_lib, twq_lib, 128);
    if (benthic_mode > 1 && n_zones <= 0) {
        fprintf(stderr, "benthic mode > 1 but no zones defined; reverting to benthic mode 1\n");
        benthic_mode = 1;
    }

    //-------------------------------------------------
    if ( get_namelist(namlst, time) ) {
        fprintf(stderr,"Error reading the 'time' namelist from %s\n", glm_nml_file);
        exit(1);
    }
    // set met, inflow and outflow data file timezones to be the default value.
    timezone_m = timezone_r;
    timezone_i = timezone_r;
    timezone_o = timezone_r;

    nDays = num_days;
    timestep = dt;

    printf("nDays %d timestep %f\n", nDays, timestep);

    //-------------------------------------------------
    create_lake(namlst);

    //-------------------------------------------------
    csv_point_nlevs = 0;
    csv_point_nvars = 0;
    csv_lake_fname = NULL;

    if ( get_namelist(namlst, output) ) {
        fprintf(stderr,"Error in output parameters specified");
        strcpy(outp_dir, ".");
        strcpy(outp_fn, "output");
        *nsave = 24;
    } else {
        strcpy(outp_dir, out_dir);
        strcpy(outp_fn, out_fn);
    }

    if ( csv_point_nlevs > MaxPointCSV ) { fprintf(stderr, "csv_point_nlevs must be < %d\n", MaxPointCSV); exit(1); }
    if ( csv_point_nvars > MaxCSVOutVars ) { fprintf(stderr, "csv_point_nvars must be < %d\n", MaxCSVOutVars); exit(1); }
    if ( csv_outlet_nvars > MaxCSVOutVars ) { fprintf(stderr, "csv_outlet_nvars must be < %d\n", MaxCSVOutVars); exit(1); }

    if ( csv_point_frombot == NULL ) {
        // CAB this is a potential source of a memory leak.
        csv_point_frombot = malloc(sizeof(LOGICAL)*csv_point_nlevs);
        for (i = 0; i < csv_point_nlevs; i++) csv_point_frombot[i] = TRUE;
    }

    if ( wq_calc ) {
        for (i = 0; i < csv_point_nvars; i++)
            set_csv_point_varname(i, csv_point_vars[i]);
    } else {
        /**********************************************************************
         * If we are not doing water quality calculations, ignore vars that   *
         * belong to WQ states                                                *
         **********************************************************************/
        int tn = csv_point_nvars, j = 0;
        for (i = 0; i < csv_point_nvars; i++) {
            if ( internal_var(csv_point_vars[i]) )
                set_csv_point_varname(j++, csv_point_vars[i]);
            else tn--;
        }
        csv_point_nvars = tn;
    }
    configure_csv(csv_point_nlevs, csv_point_at, csv_point_fname, csv_point_frombot, csv_point_nvars, csv_lake_fname);

    if ( wq_calc ) {
        for (i = 0; i < csv_outlet_nvars; i++)
            set_csv_outfl_varname(i, csv_outlet_vars[i]);
    } else {
        int tn = csv_outlet_nvars, j = 0;
        for (i = 0; i < csv_outlet_nvars; i++) {
            if ( internal_var(csv_outlet_vars[i]) )
                set_csv_outfl_varname(j++, csv_outlet_vars[i]);
            else tn--;
        }
        csv_outlet_nvars = tn;
    }
    configure_outfl_csv(csv_outlet_allinone, csv_outlet_fname, csv_outlet_nvars, csv_ovrflw_fname);

    //--------------------------------------------------------------------------
    // met
    wind_factor = 1.0;
    sw_factor = 1.0;
    lw_factor = 1.0;
    at_factor = 1.0;
    rh_factor = 1.0;
    rain_factor = 1.0;

    if ( get_namelist(namlst, meteorology) ) {
        fprintf(stderr,"Error reading 'meteorology' from namelist file %s\n", glm_nml_file);
        exit(1);
    }

    if ( lw_type == NULL )
        lw_ind = LW_CC;
    else if ( strcmp(lw_type, "LW_CC") == 0 )
        lw_ind = LW_CC;
    else if ( strcmp(lw_type, "LW_IN") == 0 )
        lw_ind = LW_IN;
    else if ( strcmp(lw_type, "LW_NET") == 0 )
        lw_ind = LW_NET;
    else {
        fprintf(stderr," Error in long wave type : '%s' unknown\n", lw_type);
        exit(1);
    }
    coef_wind_drag = CD;

    //--------------------------------------------------------------------------
    // snowice
    snow_albedo_factor = 1.0;
    snow_rho_max       = 300.0;
    snow_rho_min       = 50.0;

    if ( get_namelist(namlst, snowice) ) {
         snow_sw = FALSE;
         fprintf(stderr,"No snow and ice section, setting default parameters and assuming no snowfall\n");
    } else
         snow_sw = TRUE;

    //--------------------------------------------------------------------------
    // fetch
    if ( get_namelist(namlst, fetch) ) {
         fetch_sw = FALSE;
    } else
         fetch_sw = TRUE;

    //--------------------------------------------------------------------------
    // sediment heat (sed_heat)
    printf("*starting sed_heat = %10.5f\n",sed_heat_Ksoil);

    sed_heat_Ksoil     = 5.0;
    sed_temp_depth     = 0.1;
//    sed_temp_mean[0]   = 9.7;
    printf("*starting sed_heat = %10.5f\n",sed_temp_depth);
  //  sed_temp_amplitude = 2.7;
  //  sed_temp_peak_doy  = 151;
    if ( get_namelist(namlst, sed_heat) ) {
         sed_heat_sw = FALSE;
         fprintf(stderr,"No sed_heat section, turning off sediment heating\n");
    } else {
         sed_heat_sw = TRUE;
         fprintf(stderr,"Sed_heat section present, simulating sediment heating\n");
    }
    printf("*sed_temp_mean = %10.5f\n",sed_temp_mean[0]);
    printf("*sed_temp_mean = %10.5f\n",sed_temp_mean[1]);



    open_met_file(meteo_fl, snow_sw, rain_sw, timefmt_m);
    config_bird(namlst);

    //-------------------------------------------------
    for (i = 0; i < MaxInf; i++) {
        for (j = 0; j < MaxPar; j++ ) {
            Inflows[i].QDown[j] = 0.0;
            Inflows[i].QIns[j] = 0.0;
            Inflows[i].SDown[j] = 0.0;
            Inflows[i].SIns[j] = 0.0;
            Inflows[i].TDown[j] = 0.0;
            Inflows[i].TIns[j] = 0.0;
            Inflows[i].DDown[j] = 0.0;
            Inflows[i].DOld[j] = 0.0;
            Inflows[i].DIIns[j] = 0.0;
            Inflows[i].InPar[j] = 0;
            for (k = 0; k < MaxVars; k++) {
                Inflows[i].WQIns[j][k] = 0.0;
                Inflows[i].WQDown[j][k] = 0.0;
            }
            Inflows[i].WQInf[j] = 0.0;
        }
        Inflows[i].HFlow = 0.0;
        Inflows[i].TotIn = 0.0;
        Inflows[i].Dlwst = 0.0;
        Inflows[i].iCnt = 0;
        Inflows[i].NoIns = 0;
    }

    num_inflows = 0;
    if ( get_namelist(namlst, inflow) ) {
        NumInf = 0;
        if ( num_inflows == 0 )
            fprintf(stderr, "No inflow config, assuming no inflows\n");
        else {
            if ( num_inflows > MaxInf )
                fprintf(stderr, "Too many inflows specified in inflow config %d > %d\n", num_inflows, MaxInf);
            else
                fprintf(stderr, "Unknown error in inflow config\n");
            exit(1);
        }
    } else {
        if ( num_inflows > MaxInf ) {
            fprintf(stderr, "Too many inflows specified in inflow config %d > %d\n", num_inflows, MaxInf);
            exit(1);
        }
        NumInf = num_inflows;

        for (i = 0; i < NumInf; i++) {
            Inflows[i].SubmFlag = (subm_flag != NULL)?subm_flag[i]:FALSE;
            Inflows[i].Alpha = strm_hf_angle[i] * Pi/PiDeg;
            Inflows[i].Phi = strmbd_slope[i] * Pi/PiDeg;
            Inflows[i].DragCoeff = strmbd_drag[i];
            Inflows[i].Factor = inflow_factor[i];

            open_inflow_file(i, inflow_fl[i], inflow_varnum, (const char**)inflow_vars, timefmt_i);
        }
    }
    einff = coef_inf_entrain;

    //-------------------------------------------------
    seepage = FALSE; seepage_rate = 0.0;

    if ( get_namelist(namlst, outflow) ) {
        fprintf(stderr, "No outflow config, assuming no outflows\n");
        NumOut = 0;
    } else {
        if ( num_outlet > MaxOut) {
            fprintf(stderr, "Too many outlets specified in outflow config %d > %d\n", num_outlet, MaxOut);
            exit(1);
        }
        NumOut = num_outlet;
        if (num_outlet<=0) num_outlet=1; // BEWARE: num_outlet is only used now for malloc
                                         // remove this if you use it for anything else!!
        if ( flt_off_sw == NULL ) {
            flt_off_sw = malloc(sizeof(LOGICAL)*num_outlet);
            for (i = 0; i < NumOut; i++) flt_off_sw[i] = FALSE;
        } else if ( outlet_type == NULL ) {
            outlet_type = malloc(sizeof(int)*num_outlet);
            for (i = 0; i < NumOut; i++) outlet_type[i] = (flt_off_sw[i])?2:1;
        }
        for (i = 0; i < NumOut; i++) {
            // Outlet_type
            if ( outlet_type != NULL ) {
                if ( (outlet_type[i] > 0) && (outlet_type[i] <= 5) ) {
                    Outflows[i].Type = outlet_type[i];
                } else {
                    fprintf(stderr, "Wrong outlet type\n");
                    exit(1);
                }
            }
            Outflows[i].FloatOff = flt_off_sw[i];
            if (Outflows[i].Type == 2 || Outflows[i].FloatOff) {
                Outflows[i].FloatOff = TRUE;
                Outflows[i].Type = 2;
            }
            if ( Outflows[i].FloatOff ) {
                if ( (outl_elvs[i] > (crest_elev-base_elev)) || (outl_elvs[i] < 0.0) ) {
                    fprintf(stderr,
                    "Floating outflow (%124lf) above surface or deeper than lake depth (%12.4lf)\n",
                                    outl_elvs[i], crest_elev - base_elev);
                    exit(1);
                }
                Outflows[i].OLev = outl_elvs[i];  // if floating outlet make it is relative to surface
            } else {
                if ( (outl_elvs[i] > crest_elev) || (outl_elvs[i] < base_elev) ) {
                    fprintf(stderr,
                    "Outflow elevation (%124lf) above crest elevation (%12.4lf) or below base elevation (%12.4lf)\n",
                                    outl_elvs[i], crest_elev, base_elev);
                    exit(1);
                }
                Outflows[i].OLev = outl_elvs[i] - Base; // else make it relative to bottom
            }
//          Outflows[i].OLen = bsn_len_outl[i];
//          Outflows[i].OWid = bsn_wid_outl[i];
            Outflows[i].Factor = outflow_factor[i];
            if ( outlet_crit != NULL )
                Outflows[i].Hcrit  = outlet_crit[i];
            if ( target_temp != NULL )
                Outflows[i].TARGETtemp  = target_temp[i]; // if more than 1 withdrawals with their depth should work with a target temperature (like "isotherm")

            if (outflow_fl[i] != NULL) open_outflow_file(i, outflow_fl[i], timefmt_o);

        }
    }
    if ( outlet_crit != NULL ) { // only relevant if we have defined it.
        if ((crit_O2 < 0) || (crit_O2_dep < base_elev) || (crit_O2_days < 1)) {
            fprintf(stderr, "crit_O2 < 0 or crit_O2_dep < base elevation or crit_O2_days < 1\n");
            exit(1);
        }
    }
    O2crit = crit_O2;
    O2critdep = crit_O2_dep;
    O2critdays = crit_O2_days;
    MIXwithdraw = mix_withdraw;
    COUPLoxy = coupl_oxy_sw;
    MINlaketemp = min_lake_temp;

    if (withdrTemp_fl != NULL) open_withdrtemp_file(withdrTemp_fl, timefmt_o);

    //-------------------------------------------------
    for (i = 1; i < MaxDif; i++) mol_diffusivity[i] = 1.25E-09;
    mol_diffusivity[0] = 0.00000014;
    NumDif = 2;

    if ( get_namelist(namlst, diffuser) )
         fprintf(stderr,"No diffuser data, setting default values\n");

    //--------------------------------------------------------------------------

    if ( timefmt != 2 && timefmt != 3 ) {
        fprintf(stderr, "invalid time format \"%d\"\n", timefmt);
        exit(1);
    }
    if ( start == NULL ) {
        fprintf(stderr, "Start date is required\n"); exit(1);
    }
    if ( timefmt == 2 ) {
        if ( stop == NULL ) {
            fprintf(stderr, "Stop date is required for timefmt == 2\n"); exit(1);
        }
    }
/*
 *  The realloc was a problem for namelist reader because realloc
 *  may free the memory originally allocated if there was not enough to fill
 *  the request for extended memory. So instead we copy the original to
 *  a newly allocated memory block and free that ourselves when done.

    if ( stop == NULL ) { stop = malloc(40); *stop = 0; }
    else if ( strlen(stop) < 39 ) stop = realloc(stop, 40);
*/
    {   char *t = malloc(40);
        *t = 0;
        if ( stop != NULL ) strncpy(t, stop, 40);
        stop = t;
    }

    if ( timefmt != 2 ) *stop = 0;

    julianday = init_time(start, stop, timefmt, &nDays);
    free(stop);
    calendar_date(julianday, &jyear, &jmonth, &jday);
    //# Days since start of the year, jyear
    jday = julianday - julian_day(jyear, 1, 1) + 1;

    *jstart = julian_day(jyear, 1, jday);

    //--------------------------------------------------------------------------

    Num_WQ_Vars = 0;

    if ( wq_calc ) {
        size_t l = strlen(wq_nml_file);

        prime_wq(wq_lib);
        wq_init_glm(wq_nml_file, &l, &MaxLayers, &Num_WQ_Vars, &Num_WQ_Ben, &Kw); // Reads WQ namelist file
        fprintf(stderr, "Num_WQ_Vars = %d\n", Num_WQ_Vars);
        if ( Num_WQ_Vars > MaxVars ) {
            fprintf(stderr, "Sorry, this version of GLM only supports %d water quality variables\n", MaxVars);
            exit(1);
        }
    }
    NumDif += Num_WQ_Vars;

    initialise_lake(namlst);

    // This is where we could map inflow, met and csv_output vars to wq vars

    if ( ! WQ_VarsIdx ) {
        WQ_VarsIdx = malloc(sizeof(int)*inflow_varnum);
        for (j = 0; j < inflow_varnum; j++) WQ_VarsIdx[j] = -1;
    }
    if ( wq_calc ) {
        /* The first 3 vars are flow, temp and salt */
        for (j = 3; j < inflow_varnum; j++) {
            size_t k =  strlen(inflow_vars[j]);
            WQ_VarsIdx[j-3] = wq_var_index_c(inflow_vars[j], &k);
        }

        if ( benthic_mode > 1 ) {
            if ( (n_zones <= 0 || zone_heights == NULL) ) {
                fprintf(stderr, "benthic mode %d must define zones\n", benthic_mode);
                exit(1);
            } else
                wq_set_glm_zones(zone_heights, &n_zones, &Num_WQ_Vars, &Num_WQ_Ben);
        }

        for (j = 0; j < NumOut; j++) {
            if ( O2name != NULL ) {
                size_t tl = strlen(O2name[j]);
                O2idx = wq_var_index_c(O2name[j],&tl);
                if (O2idx < 0) {
                    fprintf(stderr, "Wrong oxygen name for outlet %3d ?\n",j+1); // How does it exit???
                    Outflows[j].O2idx = -1;
                } else  {
                    Outflows[j].O2idx = O2idx;
                }
            }
        }

        wq_set_glm_data(Lake, &MaxLayers, &MetData, &SurfData, &dt);
    }


    get_namelist(namlst, debugging);

    close_namelist(namlst);  // Close the glm.nml file

#if DEBUG
    debug_initialisation(0);
#endif
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 * This routine sets up the lake morphometry                                  *
 ******************************************************************************/
void create_lake(int namlst)
{
    /*---------------------------------------------
     * morphometry
     *-------------------------------------------*/
    char           *lake_name;
//  AED_REAL        Latitude;      // global variable
//  AED_REAL        Longitude;     // global variable
//  AED_REAL        base_elev;     // module variable
//  AED_REAL        crest_elev;    // module variable
    AED_REAL        bsn_len;
    AED_REAL        bsn_wid;
    int             bsn_vals;
    AED_REAL       *H = NULL;
    AED_REAL       *A = NULL;
    AED_REAL       *V = NULL;
    /*-------------------------------------------*/


    int kar;                     // first layer with a positive area
    int ksto;                    // first layer with a positive storage
#ifndef _VISUAL_C_
    // The visual c compiler doesn't like this so must malloc manually
    AED_REAL alpha_b[MaxLayers]; // interpolation coefficient for volume
    AED_REAL beta_b[MaxLayers];  // interpolation coefficient for area
#else
    AED_REAL *alpha_b;           // interpolation coefficient for volume
    AED_REAL *beta_b;            // interpolation coefficient for area
#endif
    int lanext;                  // temporary counter for interpolating area
    int lvnext;                  // temporary counter for interpolating volume
    AED_REAL x, y;
    int i,z,b;
    int ij,mi;

    NAMELIST morphometry[] = {
          { "morphometry",       TYPE_START,            NULL               },
          { "lake_name",         TYPE_STR,              &lake_name         },
          { "Latitude",          TYPE_DOUBLE,           &Latitude          },
          { "Longitude",         TYPE_DOUBLE,           &Longitude         },
          { "base_elev",         TYPE_DOUBLE,           &base_elev         },
          { "crest_elev",        TYPE_DOUBLE,           &crest_elev        },
          { "bsn_len",           TYPE_DOUBLE,           &bsn_len           },
          { "bsn_wid",           TYPE_DOUBLE,           &bsn_wid           },
          { "bsn_vals",          TYPE_INT,              &bsn_vals          },
          { "H",                 TYPE_DOUBLE|MASK_LIST, &H                 },
          { "A",                 TYPE_DOUBLE|MASK_LIST, &A                 },
          { "V",                 TYPE_DOUBLE|MASK_LIST, &V                 },
          { NULL,                TYPE_END,              NULL               }
    };
    AED_REAL h_z = 0.;

/*----------------------------------------------------------------------------*/

    base_elev = MISVAL;
    crest_elev = MISVAL;
    //-------------------------------------------------
    if ( get_namelist(namlst, morphometry) ) {
        fprintf(stderr,"Error reading the 'morphometry' namelist from %s\n", glm_nml_file);
        exit(1);
    }

    if (base_elev != MISVAL || crest_elev != MISVAL) {
        fprintf(stderr, "values for base_elev and crest_elev are no longer used\n");
    }
    if ( V != NULL ) {
        fprintf(stderr, "values for V are no longer used\n");
      //free(V);
        V = NULL;
    }

    base_elev = H[0];  crest_elev = H[bsn_vals-1];

    printf("Maximum lake depth is %f\n", crest_elev - base_elev);

    if ( (MaxLayers * DMax) < (crest_elev - base_elev) ) {
        fprintf(stderr, "Configuration Error. MaxLayers * max_layer_height < depth of the lake\n");
        exit(1);
    }

    if ( n_zones > 0 && zone_heights != NULL ) {
        if ( zone_heights[n_zones-1] < (crest_elev-base_elev) ) {
            fprintf(stderr, "WARNING last zone height is less than maximum depth\n");
            fprintf(stderr, "   adding an extra zone to compensate\n");
            zone_heights = realloc(zone_heights, (n_zones+2)*sizeof(AED_REAL));
            if ( zone_heights == NULL) {
                fprintf(stderr, "Memory error ...\n"); exit(1);
            }
            zone_heights[n_zones++] = (crest_elev-base_elev)+1;
        }

        zone_area = malloc(n_zones * sizeof(AED_REAL));
    }

    Lake = malloc(sizeof(LakeDataType)*MaxLayers);
    memset(Lake, 0, sizeof(LakeDataType)*MaxLayers);
    for (i = 0; i < MaxLayers; i++) Lake[i].ExtcCoefSW = Kw;

    Base = H[0];
    ksto = 0;
    kar = 0;

    if ( V == NULL ) V = malloc(sizeof(AED_REAL)*bsn_vals);
    V[0] = 0.;
    for (b = 1; b < bsn_vals; b++) {
        if ( (A[b] < A[b-1]) || (H[b] < H[b-1]) ) {
            fprintf(stderr, "Error. H and A in morphometry must be monotonically increasing\n");
            fprintf(stderr, "A[%d] = %f; A[%d] = %f; H[%d] = %f; H[%d] = %f\n",
                             b-1, A[b-1], b, A[b], b-1, H[b-1], b, H[b]);
            exit(1);
        }
        V[b] = V[b-1] + (  (A[b-1]+(A[b]-A[b-1])/2.0) * (H[b] - H[b-1]));
    }

    z = 0;
    for (b = 0; b < bsn_vals; b++) {
        H[b] -= Base;

        if (A[b] <= 0.0 ) kar++;
        if (H[b] <= 0.0 ) ksto++;

        /* Create the zone areas */
        if (benthic_mode > 1 && z < n_zones) {
            if ( zone_heights[z] <= H[b] ) {
                zone_area[z] = A[b];
                if ( b > 0 ) {
                    zone_area[z] += A[b-1];
                    zone_area[z] /= 2;
                }
                z++;
            }
        }
    }
    MaxArea = A[bsn_vals-1];

    /**************************************************************************
     * The model creates a refined lookup-table of depth-area-volume for later*
     * use. The maximum number of elements in the internal lookup table is    *
     * defined as Nmorph and calculated based on the highest supplied lake    *
     * depth index into storage arrays may be calculated as 10* the maximum   *
     * lake depth. Since the surface layer height is calculated after inflows *
     * and outflows, the height may be temporarily above the crest level,     *
     * and therefore 10 additional layers are included                        *
     **************************************************************************/
    Nmorph = ( ( H[bsn_vals-1] * MphInc ) + 1.0 / 1000.0 ) + 10;

    allocate_storage();

    CrestLevel = crest_elev - Base;
    LenAtCrest = bsn_len;
    WidAtCrest = bsn_wid;

#ifdef _VISUAL_C_
    alpha_b = malloc(sizeof(AED_REAL) * MaxLayers);
    beta_b = malloc(sizeof(AED_REAL) * MaxLayers);
#endif
    // Loop from the bottom to top of the provided depth points given in
    // &morphometry to calculate the bathymetric interpolation coefficients,
    // "a" and "b", at each level
    for (b = 1; b < (bsn_vals-1); b++) {
        if (V[b] > 0.0)
            alpha_b[b] = log10(V[b+1]/V[b]) / log10(H[b+1] / H[b]);

        dbgprt( " b = %2d V[b+1] = %24.18e V[b] = %24.18e ALOG10 = %24.18e\n", b, V[b+1], V[b], log10(V[b+1]/V[b]));
        dbgprt( " b = %2d H[b+1] = %24.18e H[b] = %24.18e ALOG10 = %24.18e\n", b, H[b+1], H[b], log10(H[b+1]/H[b]));
        dbgprt( " b = %2d  alpha_b[b] = %24.18e\n", b, alpha_b[b]);

        if (A[b] > 0.0)
            beta_b[b] = log10(A[b+1]/A[b]) / log10(H[b+1] / H[b]);

        dbgprt( " b = %2d A[b+1] = %24.18e A[b] = %24.18e ALOG10 = %24.18e\n", b, A[b+1], A[b], log10(A[b+1]/A[b]));
        dbgprt( " b = %2d H[b+1] = %24.18e H[b] = %24.18e ALOG10 = %24.18e\n", b, H[b+1], H[b], log10(H[b+1]/H[b]));
        dbgprt( " b = %2d beta_b[b] = %24.18e\n", b, beta_b[b]);
    }
    // The values of exponents  a and b for layer 0 are not used as the
    // area and volume are assumed to vary linearly from the
    // lake bottom to the top of the first layer
    alpha_b[0] = 1.0;
    beta_b[0] = 1.0;

    alpha_b[bsn_vals-1] = alpha_b[bsn_vals-2];
    beta_b[bsn_vals-1] = beta_b[bsn_vals-2];

    // Now prepare the refined depth-area-volume relationship using finer depth
    // increments (MphInc) to support layer interpolations later in the simulation
    // Note: kar is the first layer with a positive A
    // and  ksto is the first layer with a positive V
    b = -1;
    for (mi = 0; mi < Nmorph; mi++) {
    	h_z = (mi+1.0)/MphInc;

        while (b != (bsn_vals-2)) {
            b++;
            if (h_z < H[b+1]) break;
        }

        // Interpolate A and V for all depths below the
        // first layer that is a non-zero table entry (ie blank)
        if (b == 0) {
            lvnext = MAX(1, ksto);
            lanext = MAX(1, kar);
            MphLevelVol[mi] = V[lvnext] * h_z / H[lvnext];
            MphLevelArea[mi] = A[lanext] * h_z / H[lanext];
        } else {
            if (b < ksto)
                MphLevelVol[mi] = V[ksto] * h_z / H[ksto];
            else
                MphLevelVol[mi] = V[b] * pow((h_z / H[b]), alpha_b[b]);

            dbgprt( " mi=%2d j = %2d V = %24.18f H = %24.18f alpha_b = %24.18f\n", mi, b, V[b], H[b], alpha_b[b]);
            dbgprt( " h_z = %24.18f and result is %24.18f pow = %24.18f\n", h_z, MphLevelVol[mi], pow(h_z / H[b], alpha_b[b]));

            if (b < kar)
                MphLevelArea[mi] = A[kar] * h_z / H[kar];
            else
                MphLevelArea[mi] = A[b] * pow((h_z / H[b]), beta_b[b]);
        }
        b--;
    }

    // Calculate change between points for volume and area
    for (mi = 0; mi < Nmorph-1; mi++) {
        dMphLevelVol[mi] = (MphLevelVol[mi+1] - MphLevelVol[mi]);
        dMphLevelArea[mi] = (MphLevelArea[mi+1] - MphLevelArea[mi]);
    }
    dMphLevelVol[Nmorph-1] = dMphLevelVol[Nmorph-2];
    dMphLevelArea[Nmorph-1] = dMphLevelArea[Nmorph-2];

    // Calculate the maximum and minimum volume and volume at crest
    VMin = MphLevelVol[Nmorph-1] * VMin;
    VMax = MphInc * VMin;

    // Calculate storage at crest level, VolAtCrest
    x = CrestLevel * MphInc;
    y = AMOD(x, 1.0);
    ij = x - y;
    if (ij >= Nmorph) {
        y += (ij - Nmorph);
        ij = Nmorph;
    }
    if (ij > 0) {
        ij--;
        VolAtCrest = MphLevelVol[ij] + y * dMphLevelVol[ij];
    } else
        VolAtCrest = MphLevelVol[0];

    memcpy(MphLevelVoldash, MphLevelVol, sizeof(AED_REAL) * Nmorph);    // MphLevelVoldash = MphLevelVol;
    memcpy(dMphLevelVolda, dMphLevelVol, sizeof(AED_REAL) * Nmorph);    // dMphLevelVolda = dMphLevelVol;
    if ( V != NULL ) free(V);
//  if ( A != NULL ) free(A);
//  if ( H != NULL ) free(H);

#ifdef _VISUAL_C_
    free(alpha_b); free(beta_b);
#endif
}


/******************************************************************************
 * Routine to set the initial values for each layer                           *
 ******************************************************************************/
void initialise_lake(int namlst)
{
    /*---------------------------------------------
     * init_profiles
     *-------------------------------------------*/
    AED_REAL        lake_depth;
    int             num_heights; // support the old way
    AED_REAL       *the_heights; // support the old way
    int             num_depths;
    AED_REAL       *the_depths;
    AED_REAL       *the_temps;
    AED_REAL       *the_sals;
    int             num_wq_vars;
    char          **wq_names;
    AED_REAL       *wq_init_vals;
    /*-------------------------------------------*/

    NAMELIST init_profiles[] = {
          { "init_profiles",     TYPE_START,            NULL               },
          { "lake_depth",        TYPE_DOUBLE,           &lake_depth        },
          { "num_heights",       TYPE_INT,              &num_heights       },
          { "the_heights",       TYPE_DOUBLE|MASK_LIST, &the_heights       },
          { "num_depths",        TYPE_INT,              &num_depths        },
          { "the_depths",        TYPE_DOUBLE|MASK_LIST, &the_depths        },
          { "the_temps",         TYPE_DOUBLE|MASK_LIST, &the_temps         },
          { "the_sals",          TYPE_DOUBLE|MASK_LIST, &the_sals          },
          { "num_wq_vars",       TYPE_INT,              &num_wq_vars       },
          { "wq_names",          TYPE_STR|MASK_LIST,    &wq_names          },
          { "wq_init_vals",      TYPE_DOUBLE|MASK_LIST, &wq_init_vals      },
          { NULL,                TYPE_END,              NULL               }
    };

    int i, j, min_layers;
    int nx, np, nz;
    int *idx = NULL;

/*----------------------------------------------------------------------------*/
    //-------------------------------------------------
    // Do the initial profiles
    //-------------------------------------------------
    num_heights = 0;
    num_depths = 0;
    lake_depth = MISVAL;
    num_wq_vars = 0;
    if ( get_namelist(namlst, init_profiles) ) {
        fprintf(stderr,"Error reading initial_profiles from namelist file %s\n", glm_nml_file);
        exit(1);
    }
    if (! wq_calc) num_wq_vars = 0;

    // Initial values for the number of levels specified in the glm.nml file
    if ( num_heights != 0 ) {
        NumLayers = num_heights;
        for (i = 0; i < num_heights; i++) {
            Lake[i].Height = the_heights[i];
            Lake[i].Temp = the_temps[i];
            Lake[i].Salinity = the_sals[i];
        }

        if (the_heights[num_depths-1] > CrestLevel) {
            fprintf(stderr, "maximum height is greater than crest level\n");
            exit(1);
        }
        num_depths = num_heights;
    } else {
        // Initial values for the number of levels specified in the glm.nml file
        NumLayers = num_depths;
        if ( lake_depth == MISVAL ) {
            fprintf(stderr, "the depths format requires a lake_depth value\n");
            exit(1);
        }
        for (i = 0, j = num_depths-1; i < num_depths; i++, j--) {
            Lake[i].Height = lake_depth - (the_depths[j] - the_depths[0]);
            Lake[i].Temp = the_temps[j];
            Lake[i].Salinity = the_sals[j];
        }

        if (the_depths[num_depths-1] > lake_depth ) {
            fprintf(stderr, "last depth is greater the specified lake depth\n");
            exit(1);
        }
    }

    // First map the wq state var names to their indices
    if ( num_wq_vars > 0 ) {
        idx = malloc(sizeof(int)*num_wq_vars);
        for (j = 0; j < num_wq_vars; j++) {
            size_t k =  strlen(wq_names[j]);
            if ((idx[j] = wq_var_index_c(wq_names[j], &k)) < 0)
                fprintf(stderr, "Cannot find \"%s\" for initial value\n", wq_names[j]);
        }
    }

    if ( (j = get_nml_listlen(namlst, "init_profiles", "wq_init_vals")) != (num_wq_vars * num_depths) )
        fprintf(stderr, "WARNING: Initial profiles problem - expected %d wd_init_vals entries but got %d\n",
                                             (num_wq_vars * num_depths), j);

    // Likewise for each of the listed wq vars
    for (j = 0; j < num_wq_vars; j++) {
        if ( num_heights != 0 ) {
            for (i = 0; i < num_depths; i++)
                if ( idx[j] >= 0 ) _WQ_Vars(idx[j],i) = wq_init_vals[j*num_depths+i];
        } else {
            int k;
            for (i = 0, k = num_depths-1; i < num_depths; i++, k--)
                if ( idx[j] >= 0 ) _WQ_Vars(idx[j],i) = wq_init_vals[j*num_depths+k];
        }
    }

    min_layers = (lake_depth / DMin) - 1;
    if (min_layers > (MaxLayers-1)/2) min_layers = (MaxLayers-1)/2;
    if (min_layers > 50) min_layers = 50;
    if (min_layers < 3) min_layers = 3;

    // Now interpolate into at least min_layers
    while (NumLayers <= min_layers) {
        for (i = botmLayer; i < NumLayers; i++) {
            nx = 2 * (surfLayer - i);
            np = surfLayer - i;
            Lake[nx].Height = Lake[np].Height;
            Lake[nx].Temp = Lake[np].Temp;
            Lake[nx].Salinity = Lake[np].Salinity;
            for (j = 0; j < num_wq_vars; j++)
                if ( idx[j] >= 0 ) _WQ_Vars(idx[j],nx) = _WQ_Vars(idx[j],np);
        }
        for (i = botmLayer+1; i < NumLayers; i++) {
            nx = 2 * i - 1;
            np = 2 * i - 2;
            nz = 2 * i - 0;
            Lake[nx].Temp = (Lake[np].Temp + Lake[nz].Temp) / 2.0;
            Lake[nx].Height = (Lake[np].Height + Lake[nz].Height) / 2.0;
            Lake[nx].Salinity = (Lake[np].Salinity + Lake[nz].Salinity) / 2.0;
            for (j = 0; j < num_wq_vars; j++)
                if ( idx[j] >= 0 ) _WQ_Vars(idx[j],nx) = (_WQ_Vars(idx[j],np) + _WQ_Vars(idx[j],nz)) / 2.0;
        }
        NumLayers = 2*NumLayers - 1;
        if ( NumLayers * 2 >= MaxLayers ) break;
    }

    // And free the temporary index map
    if ( num_wq_vars > 0 ) free(idx);

    // calculate the density from the temp and salinity just read in
    for (i = botmLayer; i < NumLayers; i++)
        Lake[i].Density = calculate_density(Lake[i].Temp, Lake[i].Salinity);

}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 * The subroutine init_time() initialises the time module by reading          *
 * a namelist and take actions according to the specifications.               *
 ******************************************************************************/
#define INIT_T_STEP       1
#define INIT_T_BEGIN_END  2
#define INIT_T_BEGIN_STEP 3

static int init_time(char *start, char *stop, int timefmt, int *nDays)
{
    int jul1=0, secs1=0, jul2, secs2;
    int nsecs;

    switch (timefmt) {
        case INIT_T_STEP:
//          strcpy(start,"2000-01-01 00:00:00");
//          read_time_string(start, &jul1, &secs1);
            fprintf(stderr, "timefmt = 1 not supported\n");
            exit(1);
            break;
        case INIT_T_BEGIN_END:
            read_time_string(start, &jul1, &secs1);
            read_time_string(stop, &jul2, &secs2);

            nsecs = time_diff(jul2, secs2, jul1, secs1);

            *nDays = jul2-jul1;
            if (nsecs < 86400 && jul1 != jul2) nDays = nDays-1;
            nsecs = nsecs - 86400*(*nDays);
            break;
        case INIT_T_BEGIN_STEP:
            read_time_string(start, &jul1, &secs1);

            nsecs = (*nDays) * 86400;
            jul2  = jul1 + (*nDays);
            secs2 = (nsecs%86400);

            write_time_string(stop, jul2, secs2);
            break;
        default:
            fprintf(stderr, "Invalid time format specified\n");
            exit(1);
            break;
    }

    return jul1;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
