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

#include "overlapadd.h"

int speech_synthesis(muT *in, muT *phase, muT *speech, int window_length, int window_shift)
{
  int i,j;
  COMPLEX *complex_fft;
  COMPLEX *complex_ifft;
  float *real_ifft;
  int nfft;

  assert(in->d2==phase->d2);
  nfft = 2*(in->d2-1);
  complex_fft = (COMPLEX*)malloc(nfft*sizeof(COMPLEX));
  complex_ifft = (COMPLEX*)malloc(nfft*sizeof(COMPLEX));
  mx_set_size(window_shift*(in->d1-1)+window_length, 1, speech); // this put values to 0 also //
  real_ifft = (float*)malloc(nfft*sizeof(float));

  for(i=0; i < in->d1; i++) {
    // builds complex fft from power and phase spectrum
    powphase2complex(in->data+i*in->d2,phase->data+i*in->d2,complex_fft,nfft);
    // ifft
    ifft_C (complex_fft, nfft, complex_ifft, nfft);
    // get speech frame
    complex2real(complex_ifft,real_ifft,nfft);
    // overlap-add
    for(j=0; j < window_length; j++) {
      // scaling to make hanning windows sum to 1   
      speech->data[i*window_shift+j] += 2.0*(float)window_shift/(float)window_length * real_ifft[j];
    }
  }
  
  free(complex_fft);
  free(complex_ifft);
  free(real_ifft);

  return 0;
}

int powphase2complex(float *pow, float *phase, COMPLEX *complexvec, int length)
{
  int i;

  for (i=0;i<=length/2;i++) {
    complexvec[i].re=sqrt(pow[i])*cos(phase[i]);
    complexvec[i].im=sqrt(pow[i])*sin(phase[i]);
  }
  for (i=1;i<length/2;i++) {
    complexvec[length-i].re=complexvec[i].re;
    complexvec[length-i].im=-complexvec[i].im;
  }
  return 0;
}

int complex2real(COMPLEX *complexvec, float *realvec, int length)
{
  int i;

  for (i=0;i<length;i++) {
    realvec[i] = complexvec[i].re;
  }  
  return 0;
}

#define swap2(a) {int t=((char*)&a)[0]; ((char*)&a)[0]=((char*)&a)[1]; ((char*)&a)[1]=t;} 
int write_wav_file(char *filename, muT *mx, int swapflag)
{
  FILE *fp;
  long i, ret = HTK_FILE_SUCCESS;
  short temp;

  if((fp = fopen(filename, "w")) == NULL) {
    return HTK_FILE_CANT_OPEN;
  }

  for(i = 0; i < mx->d1 * mx->d2; i++) {
    temp = (short)round(mx->data[i]);
    if (swapflag) {
      swap2(temp);
    }
    if (fwrite(&temp, sizeof(short), 1, fp) != 1) {
      ret = HTK_FILE_IO_ERROR;
    }
  }
 
  fclose(fp);
  return ret;
}

// Added by gelbart to allow floating point output from
// the nr tool:
//  - writes raw float wavfile
//  - no swapping is performed 
int write_wav_file_float(char *filename, muT *mx)
{
  FILE *fp;
  long  ret = HTK_FILE_SUCCESS;

  if((fp = fopen(filename, "w")) == NULL) {
    return HTK_FILE_CANT_OPEN;
  }

  if (fwrite(mx->data, sizeof(float), mx->d1 * mx->d2, fp) 
      !=  mx->d1 * mx->d2) {
    ret = HTK_FILE_IO_ERROR;
  }
 
  fclose(fp);
  return ret;
}

