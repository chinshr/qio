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
#include <stdio.h>
#include <assert.h>
#include <math.h>

#ifndef _NOISECOMP_H
#define _NOISECOMP_H

typedef struct NCParameters {
  float filter_power;
  float overest_low;
  float overest_up;
  float overest_snrlow;
  float overest_snrup;
  int h_thresh_flag;
  float h_thresh;
  float time_alpha;
  int time_delay;
  int freq_length;
  int freq_flag;
  int start_freq_band;
  float h_scale;
  int addnoise_flag;
  float addnoise_low;
  float addnoise_up;
  float addnoise_scale;
  float *addnoise_powspec;
} SPSParams;

//Noise estimation using minima tracking
//constants for noise estimation
#define th1 0.81
#define th2 0.75 //0.92//0.75
#define th3 0.998 //0.996 //0.998
#define th4 0.99 //0.96 //0.99
#define subf 0.02
#define osub 1.4 //1 //1.4
#define fr_len 7
#define D 80
#define th5 0.6
#define D_sm 5
#define short_win 16
#define num_fr_beg 0.5 //0.5
#define NOISE_FILE_CANT_OPEN      1
#define NOISE_FILE_IO_ERROR       2
#define NOISE_FILE_SUCCES         3
#define NOISE_ESTIMATION_SUCCES   4
int minima_tracking_noise_estimation(muT *mx_P_n,muT *mx_in);

//Noise compensation main function
int noise_compensation(muT *in, muT *nse, SPSParams *ncpars, muT *out, muT *out2);

//Estimate the filter
int estimate_timesmooth_filter(muT *in,muT *nse,int frmind, SPSParams *ncpars,muT *instfilt,muT *tsfilt);
int estimate_freqsmooth_filt(int frmind, SPSParams *ncpars, muT *instfilt,muT *tsfilt, muT *finalfilt);
int estimate_freqsmooth_filt_high_samp(int frmind, SPSParams *ncpars, muT *instfilt, muT *tsfilt, muT *finalfilt);

//Noise filtering
int noise_comp_filter(muT *in, muT *filt, int frmind, muT *nsespec, SPSParams *ncpars, muT *out, float noise_power, float addnoise_power, float *clean_en, float *noisy_en);

//Carlos filter
int carlos_filter(muT *in, muT *out);


#endif
