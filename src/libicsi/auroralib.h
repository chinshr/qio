/* Qualcomm, Inc., San Diego, USA                                            *
 * OGI School of Science and Engineering at OHSU (OGI), Oregon, USA          *
 * Internation Computer Science Institute (ICSI), California, USA            *
 *---------------------------------------------------------------------------*
 * Copyright(C) 2002 Qualcomm, Inc., OGI, ICSI                               *
 * All rights reserved.                                                      *
 *                                                                           *
 * U.S. Patent Number 5,450,522.  Patents pending as of November 2004	     *
 * include U.S. application numbers 10/137,633 (Distributed voice	     *
 * recognition system utilizing multistream network feature processing),     *
 * 10/059,737 (System and method for computing and transmitting		     *
 * parameters in a distributed voice recognition system) and 09/714,806	     *
 * (Nonlinear mapping for feature extraction in automatic speech	     *
 * recognition).                                                             *
 *                                                                           *
 * A limited license without payment of royalty is hereby granted for        *
 * non-commercial research and development use only. Redistribution in       *
 * source code form is also permitted for non-commercial research and        *
 * development use only provided that the above copyright notice, this set   *
 * of conditions and the following disclaimer are retained. Redistribution   *
 * in binary form is also permitted  for non-commercial research and         *
 * development use only provided that the above copyright notice, this set   *
 * of conditions and the following disclaimer are reproduced in the          *
 * documentation and/or other materials provided with the distribution.  All *
 * other rights, including all commercial exploitation, are retained.        *
 *                                                                           *
 * THIS SOFTWARE IS PROVIDED BY THE AUTHORS AND CONTRIBUTORS ``AS IS'' AND   *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE     *
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE*
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE  *
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL*
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS   *
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)     *
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT*
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY *
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF    *
 * SUCH DAMAGE.                                                              *
 *---------------------------------------------------------------------------*/
 
#include "libogi/fe_util.h"
#include "libtcts/dsp.h"
#include <math.h>

#ifndef _AURORALIB_H
#define _AURORALIB_H

#ifdef  min
#undef  min
#endif
#define min(a, b)       ((a) < (b) ? (a) : (b))

#ifdef  max
#undef  max
#endif
#define max(a, b)       ((a) < (b) ? (b) : (a))

#define SQR(X) ((X)*(X))

#define debugflag 0
int fprintf_debug(char* message);

/* lda filtering functions */
void get_lda_filter(char *filterfile, muT *BCoef, muT *ACoef, int *latency, int nfilters);
void lda_filter(char *filterfile, muT *in, muT *fout);

/* low-pass filtering functions */
void get_lp_filter(muT *BCoef, muT *ACoef, int *latency);
void lp_filter(muT *in, muT *fout);

/* low-pass filtering functions */
void get_vadlp_filter(muT *BCoef, muT *ACoef, int *latency);
void vadlp_filter(muT *in, muT *fout);

/* general IIR filter
   initl is the number of samples that are used for
   initialisation */
void mx_filter(muT *BCoef, muT *ACoef, int *latency, muT *in, muT *fout, int initl);

/* 1/2 downsampling filter */
void lp_iir_filter(muT *in, muT *fout);

/* general IIR/FIR  filtering routine */
void lda_iir_filter(char *filterfile, muT *in, muT *fout);

/* DC offset filter, order 1 */
void DCOffsetFilter( float *CircBuff, long BSize, int *BPointer, long nSamples, float tap);

/* Simple energy-based VAD */
void mx_vad_combine(float *silfea, float *silfea_energy_vad, int length);
float *mx_energy_vad(muT *in);

/* Reading of silence probabilities from file */
float *read_sil_fea(char* ipfilname, long feat_length);
float *read_sil_fea_down_up(char* ipfilname, long feat_length, int down_up);    

/* use vad decisions to stop updating the mean and var during silence */
int online_normalization_with_vad(muT *in, float alpha, float bias,  muT *out, int latency, char *sil_flag);

/* utterqance-based normalization + use vad decisions to stop updating the mean and var during silence */
int normalization_with_vad(muT *in, muT *out, char *sil_flag);

/* Cat different vectors into a single vector */
int mx_cat_vectors(muT *in1, muT *in2, muT *out);

/* Average 2 vectors */
int mx_average_vectors(muT *in1, muT *in2, muT *out);

/* */
int mx_sub_vectors(muT *in, muT *out, int N);

/* computation of dynamic features (HTK compatible) : this is essentially to compute
the deltas before frame-dropping*/
int mx_dynamic_htk(muT* in, muT* out, int deltawin);
int mx_ddynamic(muT* in, muT* out, int deltawin);

/* estimates noise spectrum based on first frames */
int mx_noise_level(muT *in, muT *noise);

/* estimates noise spectrum based on first frames + simple VAD */
int mx_noise_level_energyfloor(muT *in, muT *noise, float threshold, float alpha);

/* make a single noise spectral estimate for the utterance
   based on frame-level VAD flags passed in array silflag (this function
   added by gelbart in 2002, for the nr tool) */
int mx_noise_level_energyfloor_vad(muT *in, muT *noise, float threshold, float alpha, char *silflag);

/* with guard */
int mx_noise_level2(muT *in, muT *noise);

/* simple resampling routine */
void simple_resample(muT *in, float ratio, muT *out);

#endif
