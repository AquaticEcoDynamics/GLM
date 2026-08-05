/*
**
**  Internal constants for GSW-TEOS-10 V3.05.
*/
/******************************************************************************
 *                                                                            *
 * This code was taken from the "Gibbs-SeaWater (GSW) Oceanographic Toolbox"  *
 *                                                                            *
 *   ---------------------------------------------------------------------    *
 *                                                                            *
 * Licence for the use of the Gibbs SeaWater (GSW) Oceanographic Toolbox      *
 *                                                                            *
 * Copyright (c) 2011, SCOR/IAPSO WG127 (Scientific Committee on Oceanic      *
 * Research/ International Association for the Physical Sciences of the       *
 * Oceans, Working Group 127).                                                *
 *                                                                            *
 * All rights reserved.                                                       *
 *                                                                            *
 * Redistribution and use, in source and binary forms, without modification,  *
 * is permitted provided that the following conditions are met:               *
 *                                                                            *
 *  * Redistributions of source code must retain the above copyright notice,  *
 *    this list of conditions and the following disclaimer.                   *
 *                                                                            *
 *  * Redistributions in binary form must reproduce the above copyright       *
 *    notice, this list of conditions and the following disclaimer in the     *
 *    documentation and/or other materials provided with the distribution.    *
 *                                                                            *
 *  * Neither the name of SCOR/IAPSO WG127 nor the names of its contributors  *
 *    may be used to endorse or promote products derived from this software   *
 *    without specific prior written permission.                              *
 *                                                                            *
 *                                                                            *
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"*
 * AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE  *
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE *
 * ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE  *
 * LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR        *
 * CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF       *
 * SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS   *
 * INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN    *
 * CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)    *
 * ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE *
 * POSSIBILITY OF SUCH DAMAGE.                                                *
 *                                                                            *
 * The software is available from http://www.TEOS-10.org                      *
 *                                                                            *
 ******************************************************************************/
#ifndef GSW_CONST_H
#define GSW_CONST_H

#define GSW_TEOS10_CONSTANTS \
static const AED_REAL \
    gsw_sfac = 0.0248826675584615, \
    offset = 5.971840214030754e-1

#define GSW_SPECVOL_COEFFICIENTS \
const AED_REAL v000 =  1.0769995862e-3, \
    v001 = -6.0799143809e-5, \
    v002 =  9.9856169219e-6, \
    v003 = -1.1309361437e-6, \
    v004 =  1.0531153080e-7, \
    v005 = -1.2647261286e-8, \
    v006 =  1.9613503930e-9, \
    v010 = -3.1038981976e-4, \
    v011 =  2.4262468747e-5, \
    v012 = -5.8484432984e-7, \
    v013 =  3.6310188515e-7, \
    v014 = -1.1147125423e-7, \
    v020 =  6.6928067038e-4, \
    v021 = -3.4792460974e-5, \
    v022 = -4.8122251597e-6, \
    v023 =  1.6746303780e-8, \
    v030 = -8.5047933937e-4, \
    v031 =  3.7470777305e-5, \
    v032 =  4.9263106998e-6, \
    v040 =  5.8086069943e-4, \
    v041 = -1.7322218612e-5, \
    v042 = -1.7811974727e-6, \
    v050 = -2.1092370507e-4, \
    v051 =  3.0927427253e-6, \
    v060 =  3.1932457305e-5, \
    v100 = -1.5649734675e-5, \
    v101 =  1.8505765429e-5, \
    v102 = -1.1736386731e-6, \
    v103 = -3.6527006553e-7, \
    v104 =  3.1454099902e-7, \
    v110 =  3.5009599764e-5, \
    v111 = -9.5677088156e-6, \
    v112 = -5.5699154557e-6, \
    v113 = -2.7295696237e-7, \
    v120 = -4.3592678561e-5, \
    v121 =  1.1100834765e-5, \
    v122 =  5.4620748834e-6, \
    v130 =  3.4532461828e-5, \
    v131 = -9.8447117844e-6, \
    v132 = -1.3544185627e-6, \
    v140 = -1.1959409788e-5, \
    v141 =  2.5909225260e-6, \
    v150 =  1.3864594581e-6, \
    v200 =  2.7762106484e-5, \
    v201 = -1.1716606853e-5, \
    v202 =  2.1305028740e-6, \
    v203 =  2.8695905159e-7, \
    v210 = -3.7435842344e-5, \
    v211 = -2.3678308361e-7, \
    v212 =  3.9137387080e-7, \
    v220 =  3.5907822760e-5, \
    v221 =  2.9283346295e-6, \
    v222 = -6.5731104067e-7, \
    v230 = -1.8698584187e-5, \
    v231 = -4.8826139200e-7, \
    v240 =  3.8595339244e-6, \
    v300 = -1.6521159259e-5, \
    v301 =  7.9279656173e-6, \
    v302 = -4.6132540037e-7, \
    v310 =  2.4141479483e-5, \
    v311 = -3.4558773655e-6, \
    v312 =  7.7618888092e-9, \
    v320 = -1.4353633048e-5, \
    v321 =  3.1655306078e-7, \
    v330 =  2.2863324556e-6, \
    v400 =  6.9111322702e-6, \
    v401 = -3.4102187482e-6, \
    v402 = -6.3352916514e-8, \
    v410 = -8.7595873154e-6, \
    v411 =  1.2956717783e-6, \
    v420 =  4.3703680598e-6, \
    v500 = -8.0539615540e-7, \
    v501 =  5.0736766814e-7, \
    v510 = -3.3052758900e-7, \
    v600 =  2.0543094268e-7

#endif /* GSW_INTERNAL_CONST_H */
