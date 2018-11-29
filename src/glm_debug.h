/******************************************************************************
 *                                                                            *
 * glm_debug.h                                                                *
 *                                                                            *
 * Developed by :                                                             *
 *     AquaticEcoDynamics (AED) Group                                         *
 *     School of Agriculture and Environment                                  *
 *     The University of Western Australia                                    *
 *                                                                            *
 *     http://aquatic.science.uwa.edu.au/                                     *
 *                                                                            *
 * Copyright 2018 -  The University of Western Australia                      *
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
#ifndef _GLM_DEBUG_H_
#define _GLM_DEBUG_H_

#define _DBG_BEFORE_  0
#define _DBG_LOOP_    1
#define _DBG_AFTER_   2

#define _DBG_MIXER_INIT_()  _dbg_mix_init_fields()
#define _DBG_MIXER_(d1, d2) {          \
  _dbg_mixer_s(d1, d2, Epi_botmLayer, Meta_topLayer);  \
  _dbg_mixer_a(Energy_AvailableMix);   \
  _dbg_mixer_a(Energy_RequiredMix);    \
  _dbg_mixer_a(redg);                  \
  _dbg_mixer_a(Epi_dz);                \
  _dbg_mixer_a(MeanSalt);              \
  _dbg_mixer_a(MeanTemp);              \
  _dbg_mixer_a(gPrimeTwoLayer);        \
  _dbg_mixer_a(Vol_Epi);               \
  _dbg_mixer_a(Mass_Epi);              \
  _dbg_mixer_a(Time_end_shear);        \
  _dbg_mixer_a(Time_start_shear);      \
  _dbg_mixer_a(Half_Seiche_Period);    \
  _dbg_mixer_a(Thermocline_Height);    \
  _dbg_mixer_a(u_avg);                 \
  _dbg_mixer_a(WindSpeedX);            \
  _dbg_mixer_a(Dens_Epil);             \
  _dbg_mixer_a(Epi_Thick);             \
  _dbg_mixer_a(dMdz);                  \
  _dbg_mixer_a(q_cub);                 \
  _dbg_mixer_a(LengthAtThermo);        \
  _dbg_mixer_a(Hypl_Thick);            \
  _dbg_mixer_a(Dens_Hypl);             \
  _dbg_mixer_a(Energy_Conv);           \
  _dbg_mixer_a(Energy_WindStir);       \
  _dbg_mixer_a(Energy_TotStir);        \
  _dbg_mixer_a(Energy_Deepen);         \
  _dbg_mixer_a(Energy_Shear);          \
  _dbg_mixer_a(del_u);                 \
  _dbg_mixer_a(u_avgSQ);               \
  _dbg_mixer_a(u_eff);                 \
  _dbg_mixer_a(u0_old);                \
  _dbg_mixer_a(u_avg_old);             \
  _dbg_mixer_a(deltaKH);               \
  _dbg_mixer_a(del_deltaKH);           \
  _dbg_mixer_a(GPEFFC);                \
  _dbg_mixer_a(accn);                  \
  _dbg_mixer_a(zsml_tilda);            \
  _dbg_mixer_a(Slope);                 \
  _dbg_mixer_a(VMsum);                 \
  _dbg_mixer_a(Tsum);                  \
  _dbg_mixer_a(Ssum);                  \
  _dbg_mixer_a(DepMX);                 \
  _dbg_mixer_a(PrevThick);             \
  _dbg_mixer_a(OldSlope);              \
  _dbg_mixer_a(Time_count_end_shear);  \
  _dbg_mixer_a(Time_count_sim);        \
  _dbg_mixer_a(FO);                    \
  _dbg_mixer_a(FSUM);                  \
  _dbg_mixer_a(u_f);                   \
  _dbg_mixer_a(u0);                    \
  _dbg_mixer_a(coef_mix_KH);           \
  _dbg_mixer_a(coef_mix_conv);         \
  _dbg_mixer_a(coef_wind_stir);        \
  _dbg_mixer_a(coef_mix_shear);        \
  _dbg_mixer_a(coef_mix_turb);         \
  _dbg_mixer_a(U_star);                \
  _dbg_mixer_a(U_star_sqr);            \
  _dbg_mixer_a(U_star_cub);            \
  _dbg_mixer_a(Epilimnion_Mid_Ht);     \
  _dbg_mixer_a(q_sqr);                 \
  _dbg_mixer_a(IntWaveSpeed);          \
  _dbg_mixer_a(ZeroMom);               \
  _dbg_mixer_a(FirstMom);              \
  _dbg_mixer_a(delzkm1);               \
  _dbg_mixer_a(Vol_Hypl);              \
  _dbg_mixer_a(Hypl_Mass);             \
  _dbg_mixer_e(); }

void _dbg_mix_init_fields(void);
void _dbg_mix_add_field(const char *f);
void _dbg_time(int jday, int iclock);
void _dbg_mixer_s(int d1, int d2, int ebl, int mtl);
void _dbg_mixer_a(AED_REAL e1);
void _dbg_mixer_e(void);
void _mix_dbg_on(void);
void _mix_dbg_off(void);

void _glm_dbg_on(void);
void _glm_dbg_off(void);

#if DEBUG

void _glm_dbg(const char *fmt, ...);
void _DumpLake(int where, int extra);
void LakeCheck(const char *str);

#else

//#define _dbg_time(jday, iclock)
#define _glm_dbg(...)
#define _DumpLake(where, extra)
#define LakeCheck(str)

#endif

#if 0
#undef _DBG_MIXER_
#undef _DBG_MIXER_INIT_
#define _DBG_MIXER_INIT_()
#define _DBG_MIXER_(d1, d2)
#endif

#endif
