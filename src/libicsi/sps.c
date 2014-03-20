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

#include "sps.h"

int mx_sps(muT *in, muT *noise, muT *out)
{
  int i,j,k;

  if(mx_set_size(in->d1, in->d2, out)) {
    return -1;
  }
  for(i=0; i < out->d1; i++) {
    MakhoulSub(&(in->data[i * out->d2]), &(noise->data[i * out->d2]), &(out->data[i * out->d2]), 23, 2, 0.001);
  }
  return 0;
}

int mx_noiselevel(muT *in, muT *noise)
{
  int i,j,k;

  if(mx_set_size(in->d1, in->d2, noise)) {
    return -1;
  }
  for(j=0; j < in->d2; j++) {
    noise->data[j]=0;
  }
  for(j=0; j < in->d2; j++) {
    for(i=0; i < 25; i++) { 
      noise->data[j] += in->data[i * in->d2 + j];
    }
    noise->data[j] /= 25;
  }
  for(i=1; i < in->d1; i++) { 
    for(j=0; j < in->d2; j++) {		       
      noise->data[i * in->d2 + j] += noise->data[j];
    }
  }
  return 0;
}

void NoiseLevelEstimate() 
{

}

void MakhoulSub(float* input, float* noise, float* output, int nbands, float alpha0, float beta)
{

  int i,j;
  float snr;
  float energy=0;
  float noiseenergy=0;
  float alpha;

  for(i=0;i<nbands;i++) {
    noiseenergy=noiseenergy+noise[i];
    energy=energy+input[i];      
  }
  
  if(energy==0) {
    snr=0;
  }
  else {
    snr=10*log10(energy/noiseenergy);
  }
   
  if(snr>=-5 && snr <=20) {
    alpha=alpha0-(snr/(20/(alpha0-1)));
  }
  if (snr>20) {
    alpha=1;
  }
  if (snr<-5) {
    alpha=(-(alpha0-1)/20)*(-5)+alpha0;
  }
   
  //alpha=10;

  for(i=0;i<nbands;i++) {
    output[i]=input[i]-(alpha*noise[i]);
    if (output[i]<(beta*noise[i])) {
      output[i]=beta*noise[i];
    }
  }

}
