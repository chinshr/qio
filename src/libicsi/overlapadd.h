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
#include <stdio.h>
#include <assert.h>

#ifndef _OVERLAPADD_H
#define _OVERLAPADD_H

/* Perform speech resynthesis using overlapping hanning windows.
   Ideally, we should use hanning windows and the frame shift
   should be equal to the frame length divided by a power of 2.
   in:            array of vectors of nfft/2+1 power spectral values
   phase:         array of vectors of nfft/2+1 phase values
   speech:        output vector of speech samples
   window_length: the analysis/synthesis window length in samples
   window_shift:  the analysis/synthesis window shift in samples */
int speech_synthesis(muT *in, muT *phase, muT *speech, int window_length, int window_shift);

// Converts power/phase representation to real/imaginary representation
int powphase2complex(float *pow, float *phase, COMPLEX *complexvec, int length);

// Gets real part of complex vector
int complex2real(COMPLEX *complexvec, float *realvec, int length);

// Write a 16 bits samples file
int write_wav_file(char *filename, muT *mx, int swap);

// Writes raw float wavfile; no swapping supported.
int write_wav_file_float(char *filename, muT *mx);

#endif 
