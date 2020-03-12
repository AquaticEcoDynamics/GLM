/******************************************************************************
 *                                                                            *
 * glm_init.c                                                                 *
 *                                                                            *
 * Developed by :                                                             *
 *     AquaticEcoDynamics (AED) Group                                         *
 *     School of Agriculture and Environment                                  *
 *     University of Western Australia                                        *
 *                                                                            *
 *     http://aquatic.science.uwa.edu.au/                                     *
 *                                                                            *
 * Copyright 2013 - 2020 -  The University of Western Australia               *
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
#include "glm_balance.h"

#include <aed_time.h>
#include <namelist.h>

#define DEFAULT_GLM_NML   "glm3.nml"
#define DEFAULT_GLM_NML_2 "glm2.nml"
#define DEFAULT_WQ_LIB    "aed2"
#define DEFAULT_WQ_NML    "aed2.nml"

//#define dbgprt(...) fprintf(stderr, __VA_ARGS__)
#define dbgprt(...) /* __VA_ARGS__ */

extern int *WQ_VarsIdx;

static AED_REAL   base_elev;
static AED_REAL   crest_elev;
static AED_REAL   max_elev;
extern LOGICAL    seepage;
extern AED_REAL   seepage_rate;

char glm_nml_file[256] = DEFAULT_GLM_NML;
char wq_lib[256] = DEFAULT_WQ_LIB;


static void create_lake(int namlst);
static void initialise_lake(int namlst);
static int init_time(const char *start, char *stop, int timefmt, int *startTOD, int *stopTOD, int *nDays);

/*############################################################################*/

/******************************************************************************
 *                                                                            *
 ******************************************************************************/
void init_glm(int *jstart, char *outp_dir, char *outp_fn, int *nsave)
{
    int jyear, jmonth, jday, julianday;
    extern int startTOD, stopTOD;

    /*-- %%NAMELIST glm_setup ------------------------------------------------*/
    char           *sim_name = NULL;
    int             max_layers = 500;
    AED_REAL        min_layer_vol = 0.1;
    AED_REAL        min_layer_thick = 0.1;
    AED_REAL        max_layer_thick = 0.1;
//  extern int      density_model;
//  extern CLOGICAL littoral_sw;
//  extern CLOGICAL non_avg;
    //==========================================================================
    NAMELIST glm_setup[] = {
          { "glm_setup",         TYPE_START,            NULL                  },
          { "sim_name",          TYPE_STR,              &sim_name             },
          { "max_layers",        TYPE_INT,              &max_layers           },
          { "min_layer_vol",     TYPE_DOUBLE,           &min_layer_vol        },
          { "min_layer_thick",   TYPE_DOUBLE,           &min_layer_thick      },
          { "max_layer_thick",   TYPE_DOUBLE,           &max_layer_thick      },
          { "density_model",     TYPE_INT,              &density_model        },
          { "littoral_sw",       TYPE_BOOL,             &littoral_sw          },
          { "non_avg",           TYPE_BOOL,             &non_avg              },
          { NULL,                TYPE_END,              NULL                  }
    };
    /*-- %%END NAMELIST ------------------------------------------------------*/

    /*-- %%NAMELIST mixing ---------------------------------------------------*/
//  extern int      surface_mixing;
//  extern AED_REAL coef_mix_conv;
//  extern AED_REAL coef_mix_eta;
//  extern AED_REAL coef_mix_ct;
//  extern AED_REAL coef_mix_cs;
//  extern AED_REAL coef_mix_KH;
//  extern AED_REAL coef_mix_hyp;
//  extern int      deep_mixing;
    //==========================================================================
    NAMELIST mixing[] = {
          { "mixing",            TYPE_START,            NULL                  },
          { "surface_mixing",    TYPE_INT,              &surface_mixing       },
          { "coef_mix_conv",     TYPE_DOUBLE,           &coef_mix_conv        },
          { "coef_wind_stir",    TYPE_DOUBLE,           &coef_wind_stir       },
          { "coef_mix_turb",     TYPE_DOUBLE,           &coef_mix_turb        },
          { "coef_mix_shear",    TYPE_DOUBLE,           &coef_mix_shear       },
          { "coef_mix_shreq",    TYPE_DOUBLE,           &coef_mix_shreq       },
          { "coef_mix_KH",       TYPE_DOUBLE,           &coef_mix_KH          },
          { "coef_mix_hyp",      TYPE_DOUBLE,           &coef_mix_hyp         },
          { "deep_mixing",       TYPE_INT,              &deep_mixing          },
          { "diff",              TYPE_DOUBLE|MASK_LIST, &mol_diffusivity      },
          { NULL,                TYPE_END,              NULL                  }
    };
    /*-- %%END NAMELIST ------------------------------------------------------*/

    /*-- %%NAMELIST wq_setup -------------------------------------------------*/
    char           *twq_lib = NULL;
    char           *wq_nml_file = DEFAULT_WQ_NML;
    int             lode_method = -1;
    int             lsplit_factor = 1;
//  LOGICAL         bioshade_feedback;
//  LOGICAL         repair_state;
//  CLOGICAL        mobility_off;
    //==========================================================================
    NAMELIST wq_setup[] = {
          { "wq_setup",          TYPE_START,            NULL                  },
          { "wq_lib",            TYPE_STR,              &twq_lib              },
          { "wq_nml_file",       TYPE_STR,              &wq_nml_file          },
          { "ode_method",        TYPE_INT,              &lode_method          },
          { "split_factor",      TYPE_INT,              &lsplit_factor        },
          { "bioshade_feedback", TYPE_BOOL,             &bioshade_feedback    },
          { "repair_state",      TYPE_BOOL,             &repair_state         },
          { "mobility_off",      TYPE_BOOL,             &mobility_off         },
          { NULL,                TYPE_END,              NULL                  }
    };
    /*-- %%END NAMELIST ------------------------------------------------------*/

    /*-- %%NAMELIST time -----------------------------------------------------*/
    int             timefmt;
    char           *start = NULL;
    char           *stop  = NULL;
    AED_REAL        dt = 0.0;      // timestep
    int             num_days = 0;  // number of days to run the sim
//  AED_REAL        timezone_r;
    //==========================================================================
    NAMELIST time[] = {
          { "time",              TYPE_START,            NULL                  },
          { "timefmt",           TYPE_INT,              &timefmt              },
          { "start",             TYPE_STR,              &start                },
          { "stop",              TYPE_STR,              &stop                 },
          { "dt",                TYPE_DOUBLE,           &dt                   },
          { "num_days",          TYPE_INT,              &num_days             },
          { "timezone",          TYPE_DOUBLE,           &timezone_r           },
          { NULL,                TYPE_END,              NULL                  }
    };
    /*-- %%END NAMELIST ------------------------------------------------------*/

    /*-- %%NAMELIST output ---------------------------------------------------*/
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
    //==========================================================================
    NAMELIST output[] = {
          { "output",            TYPE_START,            NULL                  },
          { "out_dir",           TYPE_STR,              &out_dir              },
          { "out_fn",            TYPE_STR,              &out_fn               },
          { "nsave",             TYPE_INT,               nsave                },
          { "csv_point_nlevs",   TYPE_INT,              &csv_point_nlevs      },
          { "csv_point_fname",   TYPE_STR,              &csv_point_fname      },
          { "csv_point_frombot", TYPE_BOOL|MASK_LIST,   &csv_point_frombot    },
          { "csv_point_at",      TYPE_DOUBLE|MASK_LIST, &csv_point_at         },
          { "csv_point_nvars",   TYPE_INT,              &csv_point_nvars      },
          { "csv_point_vars",    TYPE_STR|MASK_LIST,    &csv_point_vars       },
          { "csv_lake_fname",    TYPE_STR,              &csv_lake_fname       },
          { "csv_outlet_allinone", TYPE_BOOL,           &csv_outlet_allinone  },
          { "csv_outlet_fname",  TYPE_STR,              &csv_outlet_fname     },
          { "csv_outlet_nvars",  TYPE_INT,              &csv_outlet_nvars     },
          { "csv_outlet_vars",   TYPE_STR|MASK_LIST,    &csv_outlet_vars      },
          { "csv_ovrflw_fname",  TYPE_STR,              &csv_ovrflw_fname     },
          { NULL,                TYPE_END,              NULL                  }
    };
    /*-- %%END NAMELIST ------------------------------------------------------*/

    /*-- * %%NAMELIST meteorology --------------------------------------------*/
    LOGICAL         met_sw = FALSE;  // Include surface meteorological forcing
    char           *lw_type = NULL;  // Type LW measurement (LW_IN/LW_CC/LW_NET)
    LOGICAL         rain_sw = FALSE; // Rainfall composition
    LOGICAL         snow_sw = FALSE; // Snowfall
    char           *meteo_fl = NULL; // Name of meteorology input file
//  int             lw_ind;          // type of longwave radiation - now in glm_input
//  LOGICAL         atm_stab;        // Account for non-neutral atmospheric stability
//  LOGICAL         subdaily;        //
//  extern AED_REAL CD;
//  extern AED_REAL CE;
//  extern AED_REAL CH;
    extern AED_REAL salt_fall;
    extern AED_REAL wind_factor;
    extern int      fetch_mode;
    extern AED_REAL fetch_aws;
    extern AED_REAL fetch_xws;
    extern char *   fetch_fws;
    extern AED_REAL sw_factor;
    extern AED_REAL lw_factor;
    extern AED_REAL lw_offset;
    extern AED_REAL at_factor;
    extern AED_REAL at_offset;
    extern AED_REAL rh_factor;
    extern AED_REAL rain_factor;
    extern int      rad_mode;
    extern int      albedo_mode;
    extern int      cloud_mode;
    extern int      light_mode;
//  extern LOGICAL  link_solar_shade;
//  extern LOGICAL  link_rain_loss;
//  extern LOGICAL  link_bottom_drag;
    char           *timefmt_m = NULL;
    extern AED_REAL timezone_m;
    //==========================================================================
    NAMELIST meteorology[] = {
          { "meteorology",       TYPE_START,            NULL                  },
          { "met_sw",            TYPE_BOOL,             &met_sw               },
          { "lw_type",           TYPE_STR,              &lw_type              },
          { "rain_sw",           TYPE_BOOL,             &rain_sw              },
          { "salt_fall",         TYPE_DOUBLE,           &salt_fall            },
          { "meteo_fl",          TYPE_STR,              &meteo_fl             },
          { "subdaily",          TYPE_BOOL,             &subdaily             },
          { "atm_stab",          TYPE_BOOL,             &atm_stab             },
          { "rad_mode",          TYPE_INT,              &rad_mode             },
          { "albedo_mode",       TYPE_INT,              &albedo_mode          },
          { "cloud_mode",        TYPE_INT,              &cloud_mode           },
          { "fetch_mode",        TYPE_INT,              &fetch_mode           },
          { "wind_factor",       TYPE_DOUBLE,           &wind_factor          },
          { "sw_factor",         TYPE_DOUBLE,           &sw_factor            },
          { "lw_factor",         TYPE_DOUBLE,           &lw_factor            },
          { "lw_offset",         TYPE_DOUBLE,           &lw_offset            },
          { "at_factor",         TYPE_DOUBLE,           &at_factor            },
          { "at_offset",         TYPE_DOUBLE,           &at_offset            },
          { "rh_factor",         TYPE_DOUBLE,           &rh_factor            },
          { "rain_factor",       TYPE_DOUBLE,           &rain_factor          },
          { "CD",                TYPE_DOUBLE,           &CD                   },
          { "CE",                TYPE_DOUBLE,           &CE                   },
          { "CH",                TYPE_DOUBLE,           &CH                   },
          { "Aws",               TYPE_DOUBLE,           &fetch_aws            }, // (for mode 1 ) scalar
          { "Xws",               TYPE_DOUBLE,           &fetch_xws            }, // (for mode 2 ) scalar?
          { "Fws",               TYPE_STR,              &fetch_fws            }, // (for mode 3 ) not sure how to do this ...
          { "catchrain",         TYPE_BOOL,             &catchrain            },
          { "rain_threshold",    TYPE_DOUBLE,           &rain_threshold       },
          { "runoff_coef",       TYPE_DOUBLE,           &runoff_coef          },
          { "time_fmt",          TYPE_STR,              &timefmt_m            },
          { "timezone",          TYPE_DOUBLE,           &timezone_m           },
          { "link_solar_shade",  TYPE_BOOL,             &link_solar_shade     },
          { "link_rain_loss",    TYPE_BOOL,             &link_rain_loss       },
          { "link_bottom_drag",  TYPE_BOOL,             &link_bottom_drag     },
     //   { "snow_sw",           TYPE_BOOL,             &snow_sw              },
          { NULL,                TYPE_END,              NULL                  }
    };
    /*-- %%END NAMELIST ------------------------------------------------------*/

    /*-- %%NAMELIST light ----------------------------------------------------*/
    extern AED_REAL   *light_extc;
    extern AED_REAL   *energy_frac;
    extern AED_REAL    Benthic_Imin;
//  AED_REAL           Kw;
    char              *Kw_file = NULL;
    //==========================================================================
    NAMELIST light[] = {
          { "light",             TYPE_START,            NULL                  },
          { "albedo_mode",       TYPE_INT,              &albedo_mode          },
          { "albedo_mean",       TYPE_DOUBLE,           &albedo_mean          },
          { "albedo_amplitude",  TYPE_DOUBLE,           &albedo_amplitude     },
          { "light_mode",        TYPE_INT,              &light_mode           },
          { "n_bands",           TYPE_INT,              &n_bands              },
          { "light_extc",        TYPE_DOUBLE|MASK_LIST, &light_extc           },
          { "energy_frac",       TYPE_DOUBLE|MASK_LIST, &energy_frac          },
          { "Benthic_Imin",      TYPE_DOUBLE,           &Benthic_Imin         },
          { "Kw",                TYPE_DOUBLE,           &Kw                   },
          { "Kw_file",           TYPE_STR,              &Kw_file              },
          { NULL,                TYPE_END,              NULL                  }
    };
    /*-- %%END NAMELIST ------------------------------------------------------*/

    /*-- %%NAMELIST inflow ---------------------------------------------------*/
    int             num_inflows    = 0;
    LOGICAL        *subm_flag      = NULL;
    char          **names_of_strms = NULL;
    AED_REAL       *strm_hf_angle  = NULL;
    AED_REAL       *strmbd_slope   = NULL;
    AED_REAL       *strmbd_drag    = NULL;
    AED_REAL       *inflow_factor  = NULL;
    char          **inflow_fl      = NULL;
    int             inflow_varnum  = 0;
    char          **inflow_vars    = NULL;
    AED_REAL        coef_inf_entrain = 0.0;
    char           *timefmt_i      = NULL;
    extern AED_REAL timezone_i;
    //==========================================================================
    NAMELIST inflow[] = {
          { "inflow",            TYPE_START,            NULL                  },
          { "num_inflows",       TYPE_INT,              &num_inflows          },
          { "subm_flag",         TYPE_BOOL|MASK_LIST,   &subm_flag            },
          { "names_of_strms",    TYPE_STR|MASK_LIST,    &names_of_strms       },
          { "strm_hf_angle",     TYPE_DOUBLE|MASK_LIST, &strm_hf_angle        },
          { "strmbd_slope",      TYPE_DOUBLE|MASK_LIST, &strmbd_slope         },
          { "strmbd_drag",       TYPE_DOUBLE|MASK_LIST, &strmbd_drag          },
          { "inflow_factor",     TYPE_DOUBLE|MASK_LIST, &inflow_factor        },
          { "inflow_fl",         TYPE_STR|MASK_LIST,    &inflow_fl            },
          { "inflow_varnum",     TYPE_INT,              &inflow_varnum        },
          { "inflow_vars",       TYPE_STR|MASK_LIST,    &inflow_vars          },
          { "coef_inf_entrain",  TYPE_DOUBLE,           &coef_inf_entrain     },
          { "time_fmt",          TYPE_STR,              &timefmt_i            },
          { "timezone",          TYPE_DOUBLE,           &timezone_i           },
          { NULL,                TYPE_END,              NULL                  }
    };
    /*-- %%END NAMELIST ------------------------------------------------------*/

    /*-- %%NAMELIST outflow --------------------------------------------------*/
    int             num_outlet     = 0;
    LOGICAL        *flt_off_sw     = NULL;
    int            *outlet_type    = NULL;
    int             crit_O2        = -1;
    int             crit_O2_dep    = -1;
    int             crit_O2_days   = -1;
    AED_REAL       *outlet_crit    = NULL;
    char          **O2name         = NULL;
    int             O2idx          = 0;
    AED_REAL       *target_temp    = NULL;
    AED_REAL        min_lake_temp  = 0.0;
    LOGICAL         mix_withdraw   = FALSE;
    extern AED_REAL outflow_thick_limit;
    extern LOGICAL  single_layer_draw;
    LOGICAL         coupl_oxy_sw   = FALSE;
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
    //==========================================================================
    NAMELIST outflow[] = {
          { "outflow",           TYPE_START,            NULL                  },
          { "num_outlet",        TYPE_INT,              &num_outlet           },
          { "outlet_type",       TYPE_INT|MASK_LIST,    &outlet_type          },
          { "crit_O2",           TYPE_INT,              &crit_O2              },
          { "crit_O2_dep",       TYPE_INT,              &crit_O2_dep          },
          { "crit_O2_days",      TYPE_INT,              &crit_O2_days         },
          { "outlet_crit",       TYPE_DOUBLE|MASK_LIST, &outlet_crit          },
          { "O2name",            TYPE_STR|MASK_LIST,    &O2name               },
          { "O2idx",             TYPE_INT,              &O2idx                },
          { "target_temp",       TYPE_DOUBLE|MASK_LIST, &target_temp          },
          { "min_lake_temp",     TYPE_DOUBLE,           &min_lake_temp        },
          { "fac_range_upper",   TYPE_DOUBLE,           &fac_range_upper      },
          { "fac_range_lower",   TYPE_DOUBLE,           &fac_range_lower      },
          { "mix_withdraw",      TYPE_BOOL,             &mix_withdraw         },
          { "coupl_oxy_sw",      TYPE_BOOL,             &coupl_oxy_sw         },
          { "flt_off_sw",        TYPE_BOOL|MASK_LIST,   &flt_off_sw           },
          { "outl_elvs",         TYPE_DOUBLE|MASK_LIST, &outl_elvs            },
          { "bsn_len_outl",      TYPE_DOUBLE|MASK_LIST, &bsn_len_outl         },
          { "bsn_wid_outl",      TYPE_DOUBLE|MASK_LIST, &bsn_wid_outl         },
          { "outflow_fl",        TYPE_STR|MASK_LIST,    &outflow_fl           },
          { "withdrTemp_fl",     TYPE_STR,              &withdrTemp_fl        },
          { "outflow_factor",    TYPE_DOUBLE|MASK_LIST, &outflow_factor       },
          { "seepage",           TYPE_BOOL,             &seepage              },
          { "seepage_rate",      TYPE_DOUBLE,           &seepage_rate         },
          { "crest_width",       TYPE_DOUBLE,           &crest_width          },
          { "crest_factor",      TYPE_DOUBLE,           &crest_factor         },
          { "outflow_thick_limit", TYPE_DOUBLE,         &outflow_thick_limit  },
          { "single_layer_draw", TYPE_BOOL,             &single_layer_draw    },
          { "time_fmt",          TYPE_STR,              &timefmt_o            },
          { "timezone",          TYPE_DOUBLE,           &timezone_o           },
          { NULL,                TYPE_END,              NULL                  }
    };
    /*-- %%END NAMELIST ------------------------------------------------------*/

    /*-- @NAMELIST mass_balance ----------------------------------------------*/
    char           *balance_fname   = NULL;
    int             balance_varnum;
    char          **balance_vars = NULL;
    char           *timefmt_b = NULL;
    AED_REAL        timezone_b;
    //==========================================================================
    NAMELIST mass_balance[] = {
          { "mass_balance",      TYPE_START,            NULL                  },
          { "balance_file",      TYPE_STR,              &balance_fname        },
          { "balance_varnum",    TYPE_INT,              &balance_varnum       },
          { "balance_vars",      TYPE_STR|MASK_LIST,    &balance_vars         },
          { "time_fmt",          TYPE_STR,              &timefmt_b            },
          { "timezone",          TYPE_DOUBLE,           &timezone_b           },
          { NULL,                TYPE_END,              NULL                  }
    };
    /*-- %%END NAMELIST ------------------------------------------------------*/

    /*-- %%NAMELIST snowice --------------------------------------------------*/
    extern AED_REAL snow_albedo_factor;
    extern AED_REAL snow_rho_max;
    extern AED_REAL snow_rho_min;
    extern AED_REAL snow_water_equivalent;
    extern AED_REAL snow_rain_compact;
    extern AED_REAL K_ice_white;
    extern AED_REAL K_ice_blue;
    extern AED_REAL K_water;
    extern AED_REAL f_sw_wl1;
    extern AED_REAL f_sw_wl2;
    extern AED_REAL attn_ice_blue_wl1;
    extern AED_REAL attn_ice_blue_wl2;
    extern AED_REAL attn_ice_white_wl1;
    extern AED_REAL attn_ice_white_wl2;
    extern AED_REAL attn_snow_wl1;
    extern AED_REAL attn_snow_wl2;
    extern AED_REAL rho_ice_blue;
    extern AED_REAL rho_ice_white;
    extern AED_REAL min_ice_thickness;
    extern AED_REAL dt_iceon_avg;

    //==========================================================================
    NAMELIST snowice[] = {
          { "snowice",               TYPE_START,        NULL                  },
          { "snow_albedo_factor",    TYPE_DOUBLE,       &snow_albedo_factor   },
          { "snow_rho_max",          TYPE_DOUBLE,       &snow_rho_max         },
          { "snow_rho_min",          TYPE_DOUBLE,       &snow_rho_min         },
          { "snow_water_equivalent", TYPE_DOUBLE,       &snow_water_equivalent},
          { "snow_rain_compact",     TYPE_DOUBLE,       &snow_rain_compact    },
          { "K_ice_white",           TYPE_DOUBLE,       &K_ice_white          },
          { "K_ice_blue",            TYPE_DOUBLE,       &K_ice_blue           },
          { "K_water",               TYPE_DOUBLE,       &K_water              },
          { "f_sw_wl1",              TYPE_DOUBLE,       &f_sw_wl1             },
          { "f_sw_wl2",              TYPE_DOUBLE,       &f_sw_wl2             },
          { "attn_ice_blue_wl1",     TYPE_DOUBLE,       &attn_ice_blue_wl1    },
          { "attn_ice_blue_wl2",     TYPE_DOUBLE,       &attn_ice_blue_wl2    },
          { "attn_ice_white_wl1",    TYPE_DOUBLE,       &attn_ice_white_wl1   },
          { "attn_ice_white_wl2",    TYPE_DOUBLE,       &attn_ice_white_wl2   },
          { "attn_snow_wl1",         TYPE_DOUBLE,       &attn_snow_wl1        },
          { "attn_snow_wl2",         TYPE_DOUBLE,       &attn_snow_wl2        },
          { "rho_ice_blue",          TYPE_DOUBLE,       &rho_ice_blue         },
          { "rho_ice_white",         TYPE_DOUBLE,       &rho_ice_white        },
          { "min_ice_thickness",     TYPE_DOUBLE,       &min_ice_thickness    },
          { "dt_iceon_avg",          TYPE_DOUBLE,       &dt_iceon_avg         },
          { NULL,                    TYPE_END,          NULL                  }
    };
    /*-- %%END NAMELIST ------------------------------------------------------*/

    /*-- %%NAMELIST fetch ----------------------------------------------------*/
    extern LOGICAL     fetch_sw;
    extern int         fetch_ndirs;
    extern AED_REAL   *fetch_dirs;
    extern AED_REAL   *fetch_scale;
    extern AED_REAL    fetch_height;
    extern AED_REAL    fetch_porosity;
    //==========================================================================
    NAMELIST fetch[] = {
          { "fetch",             TYPE_START,            NULL                  },
          { "fetch_sw",          TYPE_BOOL,             &fetch_sw             },
          { "num_dir",           TYPE_INT,              &fetch_ndirs          },
          { "wind_dir",          TYPE_DOUBLE|MASK_LIST, &fetch_dirs           },
          { "fetch_scale",       TYPE_DOUBLE|MASK_LIST, &fetch_scale          },
          { "edge_height",       TYPE_DOUBLE,           &fetch_height         },
          { "edge_porosity",     TYPE_DOUBLE,           &fetch_porosity       },
          { NULL,                TYPE_END,              NULL                  }
    };
    /*-- %%END NAMELIST ------------------------------------------------------*/

    /*-- %%NAMELIST sediment -------------------------------------------------*/
//  int              benthic_mode;
//  int              n_zones;
    AED_REAL        *zone_heights = NULL;
    extern CLOGICAL  sed_heat_sw;
    extern int       sed_heat_model;
    extern AED_REAL  sed_heat_Ksoil;
    extern AED_REAL  sed_temp_depth;
    extern AED_REAL *sed_temp_mean;
    extern AED_REAL *sed_temp_amplitude;
    extern AED_REAL *sed_temp_peak_doy;
    extern AED_REAL *sed_reflectivity;
    extern AED_REAL *sed_roughness;
//  extern AED_REAL  sed_temp_amplitude;
//  extern AED_REAL  sed_temp_peak_doy;
    //==========================================================================
    NAMELIST sediment[] = {
          { "sediment",          TYPE_START,            NULL                  },
          { "benthic_mode",      TYPE_INT,              &benthic_mode         },
          { "n_zones",           TYPE_INT,              &n_zones              },
          { "zone_heights",      TYPE_DOUBLE|MASK_LIST, &zone_heights         },
          { "sed_reflectivity",  TYPE_DOUBLE|MASK_LIST, &sed_reflectivity     },
          { "sed_roughness",     TYPE_DOUBLE|MASK_LIST, &sed_roughness        },
          { "sed_temp_mean",     TYPE_DOUBLE|MASK_LIST, &sed_temp_mean        },
          { "sed_temp_amplitude",TYPE_DOUBLE|MASK_LIST, &sed_temp_amplitude   },
          { "sed_temp_peak_doy", TYPE_DOUBLE|MASK_LIST, &sed_temp_peak_doy    },
          { "sed_heat_Ksoil",    TYPE_DOUBLE,           &sed_heat_Ksoil       },
          { "sed_temp_depth",    TYPE_DOUBLE,           &sed_temp_depth       },
          { "sed_heat_model",    TYPE_INT,              &sed_heat_model       },
          { NULL,                TYPE_END,              NULL                  }
    };
    /*-- %%END NAMELIST ------------------------------------------------------*/

    /*-- %%NAMELIST debugging ------------------------------------------------*/
//  extern CLOGICAL dbg_mix;   //# debug output from mixer
//  extern CLOGICAL no_evap;   //# turn off evaporation
    //==========================================================================
    NAMELIST debugging[] = {
          { "debugging",         TYPE_START,            NULL                  },
          { "debug_mixer",       TYPE_BOOL,             &dbg_mix              },
          { "disable_evap",      TYPE_BOOL,             &no_evap              },
          { NULL,                TYPE_END,              NULL                  }
    };
    /*-- %%END NAMELIST ------------------------------------------------------*/

    //-------------------------------------------------
    int i, j, k;
    int namlst;

/*----------------------------------------------------------------------------*/

    //-------------------------------------------------
    // Open the namelist file.
    if ( (namlst = open_namelist(glm_nml_file)) < 0 ) {
        fprintf(stderr,"\n     ERROR opening the glm namelist file %s\n", glm_nml_file);
        if (strcmp(glm_nml_file, DEFAULT_GLM_NML) == 0) {
            fprintf(stderr, "     Trying %s\n", DEFAULT_GLM_NML_2);
            strcpy(glm_nml_file, DEFAULT_GLM_NML_2);
            if ( (namlst = open_namelist(glm_nml_file)) < 0 ) {
                fprintf(stderr,"\n     ERROR opening the glm namelist file %s\n", glm_nml_file);
                exit(1);
            }
        } else
            exit(1);
    }

    fprintf(stderr, "\n     Reading configuration from %s\n", glm_nml_file);

    //-------------------------------------------------
    // Set some default values
    coef_inf_entrain = 0.;
    Kw = 0.2;

    //-------------------------------------------------
    if ( get_namelist(namlst, glm_setup) != 0 ) {
       fprintf(stderr,"\n     ERROR reading the 'glm_setup' namelist from %s\n", glm_nml_file);
       exit(1);
    }

    //-------------------------------------------------
    for (i = 1; i < MaxDif; i++) mol_diffusivity[i] = 1.25E-09;
    mol_diffusivity[0] = 0.00000014;
    NumDif = 2;

    //-------------------------------------------------
    if ( get_namelist(namlst, mixing) ) {
        fprintf(stderr,"\n     ERROR reading the 'mixing' namelist from %s\n", glm_nml_file);
        exit(1);
    }

    MaxLayers = max_layers;
    VMin = min_layer_vol;
    DMin = min_layer_thick;
    DMax = max_layer_thick;
    NumLayers = 0;
    n_zones = 0;

    //-------------------------------------------------
    wq_calc   = TRUE;
    if ( get_namelist(namlst, wq_setup) ) {
        // fprintf(stderr, "No WQ config\n");
        twq_lib           = DEFAULT_WQ_LIB;
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

    //-------------------------------------------------
    if ( get_namelist(namlst, time) ) {
        fprintf(stderr,"\n     ERROR reading the 'time' namelist from %s\n", glm_nml_file);
        exit(1);
    }
    // set met, inflow and outflow data file timezones to be the default value.
    timezone_m = timezone_r;
    timezone_i = timezone_r;
    timezone_o = timezone_r;

    nDays = num_days;
    timestep = dt;

    if (quiet < 2) printf("\n     nDays= %d; timestep= %f (s)\n", nDays, timestep);

    //-------------------------------------------------
    create_lake(namlst);

    //-------------------------------------------------
    csv_point_nlevs = 0;
    csv_point_nvars = 0;
    csv_lake_fname = NULL;

    if ( get_namelist(namlst, output) ) {
        fprintf(stderr,"\n     ERROR in output parameters specified");
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

    if ( csv_point_frombot == NULL && csv_point_nlevs > 0) {
        // CAB this is a potential source of a memory leak.
        csv_point_frombot = calloc(csv_point_nlevs, sizeof(LOGICAL));
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
    // meteorology
    wind_factor = 1.0;
    sw_factor = 1.0;
    lw_factor = 1.0;
    at_factor = 1.0;
    at_offset = 0.0;
    rh_factor = 1.0;
    rain_factor = 1.0;

    if ( get_namelist(namlst, meteorology) ) {
        fprintf(stderr,"\n     ERROR reading 'meteorology' from namelist file %s\n", glm_nml_file);
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
        fprintf(stderr,"\n     ERROR in long-wave type : '%s' unknown\n", lw_type);
        exit(1);
    }
    coef_wind_drag = CD;
    coef_wind_chwn = CH;

    //-------------------------------------------------
    if ( get_namelist(namlst, light) ) {
        fprintf(stderr,"\n     ERROR reading the 'light' namelist from %s\n", glm_nml_file);
        exit(1);
    }
    if ( Kw_file != NULL ) open_kw_file(Kw_file, timefmt_m);

    for (i = 0; i < MaxLayers; i++) Lake[i].ExtcCoefSW = Kw;

    //--------------------------------------------------------------------------
    // snowice
    snow_albedo_factor = 1.0;
    snow_rho_max       = 300.0;
    snow_rho_min       = 50.0;

    if ( get_namelist(namlst, snowice) ) {
         snow_sw = FALSE;
         fprintf(stderr,"     No 'snowice' section, setting defaults & assuming no snowfall\n");
    } else
         snow_sw = TRUE;

    //--------------------------------------------------------------------------
    // fetch
    if ( get_namelist(namlst, fetch) )
         fetch_sw = FALSE;
    else
         fetch_sw = TRUE;

    //--------------------------------------------------------------------------
    // sediment heat
    sed_heat_Ksoil     = 5.0;
    sed_temp_depth     = 0.1;
    if ( get_namelist(namlst, sediment) ) {
        sed_heat_sw = FALSE;
        if (quiet < 2) fprintf(stderr,"     No 'sediment' section, turning off sediment heating\n");
    } else {
        sed_heat_sw = TRUE;
        if (quiet < 2 && sed_heat_model > 0 ) fprintf(stderr,"     'sediment' section present, simulating sediment heating\n");
      //  InitialTemp(int *m, const AED_REAL *depth, const AED_REAL *wv,
      //                           const AED_REAL *topTemp, const AED_REAL *botTemp,
      //                           const AED_REAL *nSPinUpDays, AED_REAL *tNew);
        if (sed_temp_mean != NULL && quiet < 2 && sed_heat_model == 1 ) {
            printf("       *sed_temp_mean[0] = %10.5f\n",sed_temp_mean[0]);
        }
    }
    if ( sed_reflectivity == NULL ) {
        int t_zones = 2;
        if ( n_zones > 1 ) t_zones = n_zones;
        sed_reflectivity = calloc(t_zones, sizeof(AED_REAL));
    }

    if ( n_zones > 0 && zone_heights != NULL ) {
        if ( zone_heights[n_zones-1] <= (max_elev-base_elev) ) {
            fprintf(stderr, "     WARNING last zone height is less than maximum depth\n");
            fprintf(stderr, "        ... adding an extra zone to compensate\n");
            zone_heights = realloc(zone_heights, (n_zones+2)*sizeof(AED_REAL));
            if ( zone_heights == NULL) {
                fprintf(stderr, "     Memory ERROR ...\n"); exit(1);
            }
            zone_heights[n_zones++] = (max_elev-base_elev)+1;
        }
        theZones = calloc(n_zones, sizeof(ZoneType));
        printf("     Sediment zones being set at %7.1f %7.1f %7.1f ... \n",zone_heights[0],zone_heights[1],zone_heights[2]);
        for (i = 0; i < n_zones; i++) theZones[i].zheight = zone_heights[i];
    }

    /**************************************************************************
     * If there are zones and these were not defined in the config they will  *
     * be NULL and access will cause segfault.                                *
     **************************************************************************/
    if ( sed_heat_sw ) {
        int t_zones = 2;
        if ( n_zones > 1 ) t_zones = n_zones;

        if (sed_roughness == NULL) {
            sed_roughness = calloc(t_zones, sizeof(AED_REAL));
        }
        if (sed_temp_mean == NULL) {
            sed_temp_mean = calloc(t_zones, sizeof(AED_REAL));
        }
        if (sed_temp_amplitude == NULL) {
            sed_temp_amplitude = calloc(t_zones, sizeof(AED_REAL));
        }
        if (sed_temp_peak_doy == NULL) {
            sed_temp_peak_doy = calloc(t_zones, sizeof(AED_REAL));
        }
    }

    if (benthic_mode > 1 && n_zones <= 0) {
        fprintf(stderr, "     NOTE: benthic_mode > 1 but no zones defined; reverting to benthic_mode 1\n");
        benthic_mode = 1;
    }

/*
fprintf(stderr, "n_zones %d\n", n_zones);
for (i = 0; i < n_zones; i++) {
    fprintf(stderr, "  sed_reflectivity[%d] = %e\n", i, sed_reflectivity[i]);
    fprintf(stderr, "  sed_temp_mean[%d] = %e\n", i, sed_temp_mean[i]);
    fprintf(stderr, "  sed_temp_amplitude[%d] = %e\n", i, sed_temp_amplitude[i]);
    fprintf(stderr, "  sed_temp_peak_doy[%d] = %e\n", i, sed_temp_peak_doy[i]);
    fprintf(stderr, "  zone_heights[%d] = %e\n", i, zone_heights[i]);
    fprintf(stderr, "  sed_roughness[%d] = %e\n", i, sed_roughness[i]);
}
*/

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
            fprintf(stderr, "     No 'inflow' config, assuming no inflows\n");
        else {
            if ( num_inflows > MaxInf )
                fprintf(stderr, "     Too many inflows specified in inflow config %d > %d\n", num_inflows, MaxInf);
            else
                fprintf(stderr, "     Unknown ERROR in inflow config\n");
            exit(1);
        }
    } else {
        if ( num_inflows > MaxInf ) {
            fprintf(stderr, "     ERROR: Too many inflows specified in 'inflow' config %d > %d\n", num_inflows, MaxInf);
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
        fprintf(stderr, "     No 'outflow' config, assuming no outflows\n");
        NumOut = 0;
    } else {
        if ( num_outlet > MaxOut) {
            fprintf(stderr, "     ERROR: Too many outlets specified in 'outflow' config %d > %d\n", num_outlet, MaxOut);
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
                    fprintf(stderr, "     ERROR: Wrong outlet type\n");
                    exit(1);
                }
            }
            Outflows[i].FloatOff = flt_off_sw[i];
            if (Outflows[i].Type == 2 || Outflows[i].FloatOff) {
                Outflows[i].FloatOff = TRUE;
                Outflows[i].Type = 2;
            }
            if ( Outflows[i].FloatOff ) {
                if ( (outl_elvs[i] > (max_elev-base_elev)) || (outl_elvs[i] < 0.0) ) {
                    fprintf(stderr,
                    "     ERROR: Depth of floating outflow (%124lf) is above lake surface or deeper than lake depth (%12.4lf)\n",
                                    outl_elvs[i], max_elev - base_elev);
                    exit(1);
                }
                Outflows[i].OLev = outl_elvs[i];  // if floating outlet make it is relative to surface
            } else {
                if ( (outl_elvs[i] > crest_elev) || (outl_elvs[i] < base_elev) ) {
                    fprintf(stderr,
                    "     ERROR: Outflow elevation (%124lf) above crest elevation (%12.4lf) or below base elevation (%12.4lf)\n",
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
            fprintf(stderr, "     ERROR: crit_O2 < 0 or crit_O2_dep < base elevation or crit_O2_days < 1\n");
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

    if ( timefmt != 2 && timefmt != 3 ) {
        fprintf(stderr, "     ERROR: invalid time format \"%d\"\n", timefmt);
        exit(1);
    }
    if ( start == NULL ) {
        fprintf(stderr, "     ERROR: Start date is required\n"); exit(1);
    }
    if ( timefmt == 2 ) {
        if ( stop == NULL ) {
            fprintf(stderr, "     ERROR: Stop date is required if timefmt == 2\n"); exit(1);
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

    julianday = init_time(start, stop, timefmt, &startTOD, &stopTOD, &nDays);
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
        fprintf(stderr, "     WQ plugin active: included Num_WQ_Vars = %d\n", Num_WQ_Vars);
        if ( Num_WQ_Vars > MaxVars ) {
            fprintf(stderr, "     ERROR: Sorry, this version of GLM only supports %d water quality variables\n", MaxVars);
            exit(1);
        }
    }
    NumDif += Num_WQ_Vars;

    initialise_lake(namlst);

    //--------------------------------------------------------------------------
    if ( !get_namelist(namlst, mass_balance) ) {
        //# this needs to happen after wq stuff has been initialised
        open_balance(out_dir, balance_fname, balance_varnum, (const char**)balance_vars, timefmt_b);
    }
    //--------------------------------------------------------------------------

    // This is where we could map inflow, met and csv_output vars to wq vars

    if ( ! WQ_VarsIdx ) {
        WQ_VarsIdx = calloc(inflow_varnum, sizeof(int));
    }
    if ( wq_calc ) {
        /* The first 3 vars are flow, temp and salt */
        for (j = 3; j < inflow_varnum; j++) {
            size_t k =  strlen(inflow_vars[j]);
            WQ_VarsIdx[j-3] = wq_var_index_c(inflow_vars[j], &k);
        }

        if ( benthic_mode > 1 ) {
            if ( (n_zones <= 0 || zone_heights == NULL) ) {
                fprintf(stderr, "     benthic_mode %d must define zones\n", benthic_mode);
                exit(1);
            } else {
                wq_set_glm_zones(theZones, &n_zones, &Num_WQ_Vars, &Num_WQ_Ben);
            }
        }

        for (j = 0; j < NumOut; j++) {
            if ( O2name != NULL ) {
                size_t tl = strlen(O2name[j]);
                O2idx = wq_var_index_c(O2name[j],&tl);
                if (O2idx < 0) {
                    fprintf(stderr, "     wrong oxygen name for outlet %3d ?\n",j+1); // How does it exit???
                    Outflows[j].O2idx = -1;
                } else  {
                    Outflows[j].O2idx = O2idx;
                }
            }
        }

        wq_set_glm_data(Lake, &MaxLayers, &MetData, &SurfData, &dt,
                                   rain_factor, sw_factor, biodrag);
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

    AED_REAL *alpha_b;           // interpolation coefficient for volume
    AED_REAL *beta_b;            // interpolation coefficient for area

    int lanext;                  // temporary counter for interpolating area
    int lvnext;                  // temporary counter for interpolating volume
    AED_REAL x, y;
//  int z, b, ij, mi;
    int b, ij, mi;

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
        fprintf(stderr,"     ERROR: reading the 'morphometry' namelist from %s\n", glm_nml_file);
        exit(1);
    }

    if (base_elev != MISVAL ) {
        fprintf(stderr, "     NOTE: value for base_elev is no longer used; A[1] is assumed.\n");
    }
    if ( V != NULL ) {
        fprintf(stderr, "     NOTE: values for V are no longer used\n");
        V = NULL;
    }

    base_elev = H[0];  max_elev = H[bsn_vals-1];
    if ( crest_elev == MISVAL ) {
        fprintf(stderr, "     NOTE: values for crest_elev not provided, assuming max elevation, H[bsn]\n");
        crest_elev = H[bsn_vals-1];
    }
    if (crest_elev > max_elev) crest_elev = max_elev;

    if (quiet < 2) {
        printf("     Maximum lake depth is %f\n", max_elev - base_elev);
        printf("     Depth where flow will occur over the crest is %f\n", crest_elev - base_elev);
    }

    if ( (MaxLayers * DMax) < (max_elev - base_elev) ) {
        fprintf(stderr, "Configuration ERROR. MaxLayers * max_layer_height < depth of the lake\n");
        exit(1);
    }

    Lake = calloc(MaxLayers, sizeof(LakeDataType));

    Base = H[0];
    ksto = 0;
    kar = 0;

    V = calloc(bsn_vals, sizeof(AED_REAL));
    V[0] = 0.;
    for (b = 1; b < bsn_vals; b++) {
        if ( (A[b] < A[b-1]) || (H[b] < H[b-1]) ) {
            fprintf(stderr, "     ERROR: H and A in morphometry must be monotonically increasing\n");
            fprintf(stderr, "     A[%d] = %f; A[%d] = %f; H[%d] = %f; H[%d] = %f\n",
                             b-1, A[b-1], b, A[b], b-1, H[b-1], b, H[b]);
            exit(1);
        }
        V[b] = V[b-1] + (  (A[b-1]+(A[b]-A[b-1])/2.0) * (H[b] - H[b-1]));
    }

//  z = 0;
    for (b = 0; b < bsn_vals; b++) {
        H[b] -= Base;

        if (A[b] <= 0.0 ) kar++;
        if (H[b] <= 0.0 ) ksto++;

        // this will never run because n_zones is 0 when create_lake is called
        /* Create the zone areas */
//      if (benthic_mode > 1 && z < n_zones) {
//          if ( theZones[z].zheight <= H[b] ) {
//              theZones[z].zarea = A[b];
//              if ( b > 0 ) {
//                  theZones[z].zarea += A[b-1];
//                  theZones[z].zarea /= 2;
//              }
//              z++;
//          }
//      }
    }
    MaxArea = A[bsn_vals-1];

    /**************************************************************************
     * The model creates a refined lookup-table of depth-area-volume for later*
     * use. The maximum number of elements in the internal lookup table is    *
     * defined as Nmorph and calculated based on the highest supplied lake    *
     * depth; index into storage arrays may be calculated as 10* the maximum  *
     * lake depth. Since the surface layer height is calculated after inflows *
     * and outflows, the height may be temporarily above the max lake level,  *
     * and therefore 10 additional layers are included                        *
     **************************************************************************/
    Nmorph = ( ( H[bsn_vals-1] * MphInc ) + 1.0 / 1000.0 ) + 10;

    allocate_storage();

    CrestHeight = crest_elev - Base;
    MaxHeight = max_elev - Base;
    LenAtCrest = bsn_len;
    WidAtCrest = bsn_wid;

    alpha_b = calloc(MaxLayers, sizeof(AED_REAL));
    beta_b = calloc(MaxLayers, sizeof(AED_REAL));

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

    // Calculate storage at maximum lake level, MaxVol
    x = MaxHeight * MphInc;
    y = AMOD(x, 1.0);
    ij = x - y;
    if (ij >= Nmorph) {
        y += (ij - Nmorph);
        ij = Nmorph;
    }
    if (ij > 0) {
        ij--;
        MaxVol = MphLevelVol[ij] + y * dMphLevelVol[ij];
    } else
        MaxVol = MphLevelVol[0];

    // Calculate storage at crest/overflow level, VolAtCrest
    x = CrestHeight * MphInc;
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

    if (quiet < 2) fprintf(stderr,"     VolAtCrest= %10.5f; MaxVol= %10.5f (m3)\n", VolAtCrest, MaxVol);

    memcpy(MphLevelVoldash, MphLevelVol, sizeof(AED_REAL) * Nmorph);    // MphLevelVoldash = MphLevelVol;
    memcpy(dMphLevelVolda, dMphLevelVol, sizeof(AED_REAL) * Nmorph);    // dMphLevelVolda = dMphLevelVol;
    if ( V != NULL ) free(V);

//  if ( A != NULL ) free(A);
//  if ( H != NULL ) free(H);

    free(alpha_b); free(beta_b);
}


/******************************************************************************
 * Routine to set the initial values for each layer                           *
 ******************************************************************************/
void initialise_lake(int namlst)
{
    /*-- %%NAMELIST init_profiles --------------------------------------------*/
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
    AED_REAL        snow_thickness = 0.0;
    AED_REAL        white_ice_thickness = 0.0;
    AED_REAL        blue_ice_thickness = 0.0;
    AED_REAL        avg_surf_temp = 6.0;
    AED_REAL       *restart_variables = NULL;

    //==========================================================================
    NAMELIST init_profiles[] = {
          { "init_profiles",       TYPE_START,            NULL                },
          { "lake_depth",          TYPE_DOUBLE,           &lake_depth         },
          { "num_heights",         TYPE_INT,              &num_heights        },
          { "the_heights",         TYPE_DOUBLE|MASK_LIST, &the_heights        },
          { "num_depths",          TYPE_INT,              &num_depths         },
          { "the_depths",          TYPE_DOUBLE|MASK_LIST, &the_depths         },
          { "the_temps",           TYPE_DOUBLE|MASK_LIST, &the_temps          },
          { "the_sals",            TYPE_DOUBLE|MASK_LIST, &the_sals           },
          { "num_wq_vars",         TYPE_INT,              &num_wq_vars        },
          { "wq_names",            TYPE_STR|MASK_LIST,    &wq_names           },
          { "wq_init_vals",        TYPE_DOUBLE|MASK_LIST, &wq_init_vals       },
          { "snow_thickness",      TYPE_DOUBLE,           &snow_thickness     },
          { "white_ice_thickness", TYPE_DOUBLE,           &white_ice_thickness},
          { "blue_ice_thickness",  TYPE_DOUBLE,           &blue_ice_thickness },
          { "avg_surf_temp",       TYPE_DOUBLE,           &avg_surf_temp      },
          { "restart_variables",   TYPE_DOUBLE|MASK_LIST, &restart_variables  },
          { NULL,                  TYPE_END,              NULL                }
    };
    /*-- %%END NAMELIST ------------------------------------------------------*/

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
        fprintf(stderr,"     ERROR: reading initial_profiles from namelist file %s\n", glm_nml_file);
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

        if (the_heights[num_depths-1] > CrestHeight) {
            fprintf(stderr, "     ERROR: maximum height is greater than crest level\n");
            exit(1);
        }
        num_depths = num_heights;
    } else {
        // Initial values for the number of levels specified in the glm.nml file
        NumLayers = num_depths;
        if ( lake_depth == MISVAL ) {
            fprintf(stderr, "     ERROR: the depths format requires a lake_depth value\n");
            exit(1);
        }
        for (i = 0, j = num_depths-1; i < num_depths; i++, j--) {
            Lake[i].Height = lake_depth - (the_depths[j] - the_depths[0]);
            Lake[i].Temp = the_temps[j];
            Lake[i].Salinity = the_sals[j];
        }

        if (the_depths[num_depths-1] > lake_depth ) {
            fprintf(stderr, "     ERROR: last depth is greater the specified lake depth\n");
            exit(1);
        }
    }

    // First map the wq state var names to their indices
    if ( num_wq_vars > 0 ) {
        idx = calloc(num_wq_vars, sizeof(int));
        for (j = 0; j < num_wq_vars; j++) {
            size_t k =  strlen(wq_names[j]);
            if ((idx[j] = wq_var_index_c(wq_names[j], &k)) < 0)
                fprintf(stderr, "Cannot find \"%s\" for initial value\n", wq_names[j]);
        }
    }

    if ( (j = get_nml_listlen(namlst, "init_profiles", "wq_init_vals")) != (num_wq_vars * num_depths) )
        fprintf(stderr, "     WARNING: Initial profiles problem - expected %d wd_init_vals entries but got %d\n",
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
    
    if (restart_variables == NULL) {
        restart_variables = calloc(17, sizeof(AED_REAL));

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
    } else {
        // Use the layers exactly as provided in the nml (used in restarting)
        NumLayers = num_depths;
    }

    // And free the temporary index map
    if ( num_wq_vars > 0 ) free(idx);

    // calculate the density from the temp and salinity just read in
    for (i = botmLayer; i < NumLayers; i++)
        Lake[i].Density = calculate_density(Lake[i].Temp, Lake[i].Salinity);

    if (littoral_sw) {
        Lake[onshoreLayer].Temp = Lake[surfLayer].Temp;
        Lake[offshoreLayer].Temp = Lake[surfLayer].Temp;
    }

    SurfData.delzWhiteIce = white_ice_thickness;
    SurfData.delzBlueIce = blue_ice_thickness;

    //Initializing with non-zero snow thickness causes segmentation faults
    SurfData.delzSnow = zero;

    if (SurfData.delzBlueIce > 0.0 || SurfData.delzWhiteIce > 0.0) {
        ice = TRUE;
    }

    AvgSurfTemp = avg_surf_temp;

    if (restart_variables != NULL) {
        DepMX = restart_variables[0];
        PrevThick =  restart_variables[1]; //# mixed layer thickness from previous time step
        gPrimeTwoLayer =  restart_variables[2]; //# Reduced gravity for int wave estimate
        Energy_AvailableMix = restart_variables[3];  //# Total available energy to mix (carries over from previous timesteps)
        Mass_Epi =  restart_variables[4];//# Sigma mass of Epilimnion (surface layer after Kelvin-Helmholtz) kg
        OldSlope = restart_variables[5];
        Time_end_shear =  restart_variables[6]; //# Time left before shear cut off [hours]
        Time_start_shear =  restart_variables[7];//# Time count since start of sim for shear period start [hours]
        Time_count_end_shear =  restart_variables[8]; //# Time count since start of sim for shear period end [hours]
        Time_count_sim = restart_variables[9];  //# Time count since start of simulation [hours]
        Half_Seiche_Period = restart_variables[10];//# One half the seiche period
        Thermocline_Height =  restart_variables[11];//# Height at the top of the metalimnion [m]
        FO = restart_variables[12];
        FSUM = restart_variables[13];
        u_f = restart_variables[14];
        u0 = restart_variables[15];
        u_avg = restart_variables[16];

//      free(restart_variables);
    }
    
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


/******************************************************************************
 * The subroutine init_time() initialises the time module by reading          *
 * a namelist and take actions according to the specifications.               *
 *                                                                            *
 * timefmt [integer]                                                          *
 *      method to specify start and duration of model run                     *
 *      1: duration computed from number of time steps, MaxN (bogus start     *
 *         date used) [no longer implemented!!]                               *
 *      2: duration computed from given start and stop dates (number of time  *
 *         steps MaxN computed)                                               *
 *      3: duration computed from number of time steps, (start date as        *
 *         specified, stop date computed)                                     *
 * start   [string, format = "yyyy-mm-dd hh:mm:ss"]                           *
 *      nominal start date                                                    *
 *      This variable is used only if timefmt != 1                            *
 * stop    [string, format = "yyyy-mm-dd hh:mm:ss"]                           *
 *      nominal stop date                                                     *
 *      This variable is used only if timefmt = 2                             *
 * dt      [float, minimum = 0.001, maximum = 86400, unit = s]                *
 *      Time step for integration                                             *
 * numb_days [number of days to run the simulation ]                          *
 *      This variable is used only if timefmt != 2                            *
 *                                                                            *
 ******************************************************************************/
#define INIT_T_STEP       1
#define INIT_T_BEGIN_END  2
#define INIT_T_BEGIN_STEP 3

static int init_time(const char *start, char *stop, int timefmt, int *startTOD, int *stopTOD, int *nDays)
{
    int jul1 = 0, jul2;
    int nsecs;
    extern int nDates;

    switch (timefmt) {
        case INIT_T_STEP:
            fprintf(stderr, "     ERROR: timefmt = 1 not supported\n");
            exit(1);
            break;
        case INIT_T_BEGIN_END:
            read_time_string(start, &jul1, startTOD);
            read_time_string(stop, &jul2, stopTOD);

            nsecs = time_diff(jul2, *stopTOD, jul1, *startTOD);

            *nDays = nsecs / iSecsPerDay;
            nDates = jul2 - jul1 + 1;
            break;
        case INIT_T_BEGIN_STEP:
            read_time_string(start, &jul1, startTOD);

            nsecs = (*nDays) * iSecsPerDay;
            jul2  = jul1 + (*nDays);
            nDates = *nDays + 1;

            write_time_string(stop, jul2, *startTOD);
            *stopTOD = *startTOD;
            break;
        default:
            fprintf(stderr, "     ERROR: Invalid time format specified\n");
            exit(1);
            break;
    }

    return jul1;
}
/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/
