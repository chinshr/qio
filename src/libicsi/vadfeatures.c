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

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "libicsi/vadfeatures.h"

int mx_autocor_based_params(muT* in, int windowsz, int framesz, int npoints, float *window, muT* out, float preemp, int sample_rate, int high_frequency)
{
  float *powspec;
  int low_freq, high_freq;
  int a;
  int count;

  int low_period, high_period;
  float norm_autocor_var, zero_xing_rate;
  int n_params;
  int n_frames;

  float *framebuff,*r;
  int i,j,h,zc;
  float m,v;
  n_params=4;
 
  if((framebuff = malloc(npoints * sizeof(float))) == NULL) {
    return -1;
  }
  if(mx_set_size((in->d1-(framesz))/framesz, n_params, out)) {
    free(framebuff);
    return -1;
  }
  n_frames = (in->d1-(framesz))/framesz;
  powspec = (float *)malloc((npoints/2+1)*sizeof(float));
  
  r = (float *)malloc((windowsz-1)*sizeof(float));

  low_period = (int)((float)sample_rate/high_frequency);
  high_period = (int)(windowsz/2);

  low_freq = 16;
  high_freq = 100;
  a = 10;

  for(i=0; i < out->d1; i++) {

    asr_memcpy(framebuff, in->data, in->d1, (i+1)*framesz, windowsz);

    if (preemp != 0) {
      preemphasis(framebuff,windowsz,preemp);
    }
    memset(framebuff+windowsz, 0, (npoints-windowsz) * sizeof(float));
    if(window) {
      for (j=0; j<windowsz; j++) {
        framebuff[j] *= window[j];
      }
    }

    // Compute normalized autocorrelation function //
    for(h=0;h<2;h++) {
      r[h]=0;
      for(j=0;j<windowsz-h;j++){
	r[h] += framebuff[j]*framebuff[j+h];  
      }
      r[h] /= (windowsz-h);
      if (h>0) {
	r[h] /= r[0];
      }
    }
    for(h=low_period;h<high_period;h++) {
      r[h]=0;
      for(j=0;j<windowsz-h;j++){
        r[h] += framebuff[j]*framebuff[j+h];  
      }
      r[h] /= (windowsz-h);
      if (h>0) {
	r[h] /= r[0];
      }
    }

    // Compute variance of normalized autocorrelation function //
    m = 0;
    v = 0;  
    for(h=low_period;h<high_period;h++) {
      v += r[h]*r[h];
      m += r[h];
    }
    m /= (high_period-low_period);
    v /= (high_period-low_period);
    norm_autocor_var = log10(v); /* - m*m); */
    
    // Compute zero xing //
    zc=0;
    for(h=1;h<windowsz;h++){
      if ( ((framebuff[h-1]>0) && (framebuff[h]<0) ) || ((framebuff[h-1]<0) && (framebuff[h]>0) ) || (framebuff[h]==0) ) {
	zc++;
      }
    }
    zero_xing_rate = (float)zc/(windowsz-1);
    
    // Compute logFFT variance //
    mu_FAST(framebuff, npoints);
    for(j=1; j < npoints/2; j++){
      powspec[j] = log(SQR(framebuff[2 * j]) + SQR(framebuff[2 * j + 1]));
    }
    powspec[npoints / 2] = log(SQR(framebuff[1]));
    powspec[0] = log(SQR(framebuff[0]));
    v=0;
    m=0;
    for(count=low_freq-1-a;count<low_freq+a;count++) {
      m += powspec[count];
    }
    m /= (2*a+1);
    for(count=low_freq;count<high_freq;count++) {
      m = m + powspec[count+a]/(2*a+1) - powspec[count-1-a]/(2*a+1);
      v += (powspec[count]-m)*(powspec[count]-m);
    }    
    v /= (high_freq-low_freq);
    
    // Prepare output vector //
    out->data[i*out->d2] = norm_autocor_var;
    out->data[i*out->d2+1] = r[1];
    out->data[i*out->d2+2] = zero_xing_rate;  
    out->data[i*out->d2+3] = v;

  }
  
  free(powspec);
  free(framebuff);
  free(r);

  return(0);

}
