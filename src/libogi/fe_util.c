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
#include <math.h>
#include <assert.h>
#include "fe_util.h"

#define SQR(X) ((X)*(X))
#define ROUND(X) (floor((X)+0.5))
#define swap4(a) {int t=((char*)&a)[0]; ((char*)&a)[0]=((char*)&a)[3]; ((char*)&a)[3]=t; t=((char*)&a)[1]; ((char*)&a)[1]=((char*)&a)[2]; ((char*)&a)[2]=t;}
#define swap2(a) {int t=((char*)&a)[0]; ((char*)&a)[0]=((char*)&a)[1]; ((char*)&a)[1]=t;}
//#define sqr(x) ((x) * (x))


#define BUNCH_SIZE 16

void asr_memcpy(float *framebuff, float *in, int n_samples, int center, int length)
{
  int left, right, true_length;
  left = center - length/2;
  right = left + length;
  
  true_length = length;
  if (left<0) {
    true_length = right;
    left = 0;
  }
  else if (right>n_samples) {
    true_length = n_samples - left;
  }
  memcpy(framebuff, in+left, true_length * sizeof(float));
  if (true_length<length) {
    memset(framebuff+true_length, 0, (length-true_length) * sizeof(float));
  }
  
}

int mx_set_size(int rows, int cols, muT *out)
{
  out->d1 = rows;
  out->d2 = cols;
  if((out->data = (float *) calloc(out->d1*out->d2, sizeof(float))) == NULL) {
    return -1;
  }
  return 0;
}

void mx_release(muT *in) {
  free(in->data);

/* set pointer to zero after free, so that there will be a
   segmentation fault if it is dereferenced */
  in->data = 0;
  
}

//Added phase for optional speech synthesis;sunil 08/01/2001
// 2001 oct 31, dupont@icsi.berkeley.edu, changed to get frames centered at multiples of frame shift
int mx_specgram_phase(muT* in, int windowsz, int framesz, int npoints, float *window, muT* out, float preemp, muT *phase)
{
  float *framebuff;
  int i, j;
  
  assert(windowsz <= npoints);
  assert(framesz <= windowsz);
  
  if((framebuff = malloc(npoints * sizeof(float))) == NULL) {
    return -1;
  }
  if(mx_set_size((in->d1-(windowsz-framesz))/framesz, npoints/2+1, out)) {
    free(framebuff);
    return -1;
  }
  if(mx_set_size((in->d1-(windowsz-framesz))/framesz, npoints/2+1, phase)) {
    free(framebuff);
    return -1;
  }
  for(i=0; i < out->d1; i++) {
    //memcpy(framebuff, in->data+i*framesz, windowsz * sizeof(float));
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
    mu_FAST(framebuff, npoints);
    for(j=1; j < npoints/2; j++){
      out->data[i*out->d2 + j] = SQR(framebuff[2 * j]) + SQR(framebuff[2 * j + 1]);
      phase->data[i*out->d2 + j] = atan2(framebuff[2 * j + 1],framebuff[2 * j]);
    }
    out->data[i*out->d2 + npoints / 2] = SQR(framebuff[1]);
    out->data[i*out->d2] = SQR(framebuff[0]);
    
    phase->data[i*out->d2 + npoints / 2] = atan2(framebuff[1],fabs(framebuff[1]));
    phase->data[i*out->d2] = atan2(framebuff[0],fabs(framebuff[0]));
  }
  free(framebuff);
  return 0;
}

void mx_log(muT* in)
{
  int i, j;
  for(i=0; i < in->d1; i++) {
    for(j=0; j < in->d2; j++) {
      if(in->data[i * in->d2 +j] < 1/*2e-22*/) {
        in->data[i * in->d2 +j] = 1/*2e-22*/;
      }
      in->data[i * in->d2 +j] = log(in->data[i * in->d2 +j]);
    }
  }
}

int get_upper_band_energy(muT *in, int fs, muT *mx_upper_band_energy)
{
  
  int i,j;
  int idx1, idx2;
  
  if(fs==11000) {
    mx_set_size(in->d1,1,mx_upper_band_energy);
    idx1 = in->d2*8000/fs;
    for(i=0;i<in->d1;i++)
      for(j=idx1;j<in->d2;j++)
	mx_upper_band_energy->data[i] += in->data[i*in->d2+j];
  } else {
    mx_set_size(in->d1,2,mx_upper_band_energy);
    idx1 = in->d2*8000/fs;
    idx2 = in->d2*11000/fs;
    for(i=0;i<in->d1;i++) {
      for(j=idx1;j<idx2;j++)
	mx_upper_band_energy->data[i*2] += in->data[i*in->d2+j];
      for(j=idx2;j<in->d2;j++)
	mx_upper_band_energy->data[i*2+1] += in->data[i*in->d2+j];
    }
  }
  
  return(0);
}

int mx_prod(muT *in1, muT *in2, muT *out)
{
  int i,j,k;
  
  assert(in1->d2==in2->d1);
  if(mx_set_size(in1->d1, in2->d2, out)) {
    return -1;
  }
  for(i=0; i < out->d1; i++) {
    for(j=0; j < out->d2; j++) {
      out->data[i * out->d2 + j] = 0;
      for(k=0; k < in1->d2; k++) {
        out->data[i * out->d2 + j] += in1->data[i * in1->d2 + k]*in2->data[k * in2->d2 + j];
      }
    }
  }
  return 0;
}

//Added on 12/18/01 Sunil
int mx_prod_high_samp(muT *in1, int fs, muT *in2, muT *out)
{
  int i,j,k;
  
  if(mx_set_size(in1->d1, in2->d2, out)) {
    return -1;
  }
  for(i=0; i < out->d1; i++) {
    for(j=0; j < out->d2; j++) {
      out->data[i * out->d2 + j] = 0;
      for(k=0; k < in2->d1; k++) {
        out->data[i * out->d2 + j] += in1->data[i * in1->d2 + k]*in2->data[k * in2->d2 + j];
      }
    }
  }
  return 0;
}

/*******************************************************************************************/

int mx_read_vec(short *data, int len, muT *out)
{
  int i;
  if(mx_set_size(len, 1, out)) {
    return -1;
  }
  for (i=0; i < len; i++) {
    out->data[i]=data[i];
  }
  return 0;
}

int mx_read_htk_file(char *htk_filename, muT *mx, x_int32 *sampPeriod, x_int16 *parmKind, int swap)
{
  FILE *fp;
  long i;
  struct htk_header header;
  
  if((fp = fopen(htk_filename, "r")) == NULL) {
    return HTK_FILE_CANT_OPEN;
  }
  if(fread(&header, sizeof(struct htk_header), 1, fp) != 1) {
    fclose(fp);
    return HTK_FILE_IO_ERROR;
  }
  if(swap) {
    swap4(header.nSamples);
    swap4(header.sampPeriod);
    swap2(header.sampSize);
    swap2(header.parmKind);
  }
  
  if(mx_set_size(header.nSamples, header.sampSize / sizeof(float), mx)) {
    fputs("Fatal: Memory allocation error", stderr);
    exit(-1);
  }
  if(fread(mx->data, header.sampSize, header.nSamples, fp) != header.nSamples) {
    mx_release(mx);
    fclose(fp);
    return HTK_FILE_IO_ERROR;
  }
  if(swap) {
    for(i = 0; i < mx->d1 * mx->d2; i++) {
      swap4(mx->data[i]);
    }
  }
  fclose(fp);

  *sampPeriod = header.sampPeriod;
  *parmKind = header.parmKind;
  
  return HTK_FILE_SUCCESS;
}

int read_raw_file(char *filename, short **data, long *length, int swap)
{
  FILE *fp;
  long i;
  
  *data = NULL;
  if((fp = fopen(filename, "rb")) == NULL) {
    return HTK_FILE_CANT_OPEN;
  }
  if((fseek(fp, 0L, SEEK_END)) != 0) {
    return HTK_FILE_IO_ERROR;
  }
  *length = ftell(fp);
  if((fseek(fp, 0L, SEEK_SET)) != 0) {
    return HTK_FILE_IO_ERROR;
  }
  *length /= sizeof(short);
  
  //  printf("Samples: %d SampSize: %d\n", header->nSamples, header->sampSize);
  *data = (short *) malloc(sizeof(short) * *length);
  if(!*data) {
    fputs("Fatal: Memory allocation error", stderr);
    exit(-1);
  }
  if(fread(*data, sizeof(short), *length, fp) != *length) {
    free(*data);
    *data = NULL;
    fclose(fp);
    return HTK_FILE_IO_ERROR;
  }
  if(swap) {
    for(i = 0; i < *length; i++) {
      swap2((*data)[i]);
    }
  }
  fclose(fp);
  return HTK_FILE_SUCCESS;
}

int mx_write_htk_file(char *htk_filename, muT *mx, x_int32 sampPeriod, x_int16 parmKind, int swap)
{
  FILE *fp;
  struct htk_header header;
  long i, ret = HTK_FILE_SUCCESS;
  
  if((fp = fopen(htk_filename, "w")) == NULL) {
    return HTK_FILE_CANT_OPEN;
  }
  header.nSamples = mx->d1;
  header.sampSize = mx->d2 * sizeof(float);
  header.sampPeriod = sampPeriod;
  header.parmKind = parmKind;
  if(swap) {
    swap4(header.nSamples);
    swap4(header.sampPeriod);
    swap2(header.sampSize);
    swap2(header.parmKind);
  }
  if(fwrite(&header, sizeof(struct htk_header), 1, fp) != 1) {
    fclose(fp);
    return HTK_FILE_IO_ERROR;
  }
  if(swap) {
    for(i = 0; i < mx->d1 * mx->d2; i++) {
      swap4(mx->data[i]);
    }
  }
  if(fwrite(mx->data, mx->d2 * sizeof(float), mx->d1, fp) != mx->d1) {
    ret = HTK_FILE_IO_ERROR;
  }
  if(swap) {
    for(i = 0; i < mx->d1 * mx->d2; i++) {
      swap4(mx->data[i]);
    }
  }
  fclose(fp);
  return ret;
}

/*******************************************************************************************/
/* dupont@icsi.berkeley.edu 6 Aug 2001: generalized hanning from previous get_hamming_window */
void get_genhanning_window(float *window, int windowsz, float alpha)
{
  int i;
  double base_angle;
  base_angle = 2.0 * M_PI / (double)(windowsz - 1);
  
  for(i=0; i< windowsz; i++) {
    window[i] = alpha - (1.0 - alpha) * cos (i * base_angle);
  }
}

void preemphasis(float *speech, int windowsz, float preemp)
{
  int i;
  for (i=0; i<windowsz-1; i++) {
    speech[i] = speech[i+1] - preemp * speech[i]; /* forward difference */
  }
  speech[i] = (1 - preemp) * speech[i];
}

int get_dcttable(int N, int K, muT *out)
{
  float base_angle, sum;
  int i, l ;
  if(mx_set_size(N, K, out)) {
    return -1;
  }
  base_angle =  M_PI / out->d1;
  for(i=0 ; i < out->d2; i++) {
    sum=0;
    for(l=0; l < out->d1; l++) {
      out->data[l * out->d2 + i] = cos(i * (l + 0.5) * base_angle);
      sum+=SQR(out->data[l * out->d2 + i]);
    }
    sum=sqrt(sum);
    for(l=0; l < out->d1; l++) {
      out->data[l * out->d2 + i] /= sum;
    }
  }
  return 0;
}

int get_dcttable_s(int N, int K, muT *dct, muT *idct)
{
  float base_angle, sum;
  int i, l ;
  
  if(mx_set_size(N, K, dct)||mx_set_size(N, K, idct)) {
    return -1;
  }
  base_angle =  M_PI / dct->d1;
  for(i=0 ; i < dct->d2; i++) {
    sum=0;
    for(l=0; l < dct->d1; l++) {
      dct->data[l * dct->d2 + i] = cos(i * (l + 0.5) * base_angle);
      sum+=SQR(dct->data[l * dct->d2 + i]);
    }
    sum=sqrt(sum);
    for(l=0; l < dct->d1; l++) {
      dct->data[l * dct->d2 + i] /= sum;
    }
  }
  for(i=0 ; i < dct->d1; i++)
    for(l=0; l < dct->d2; l++)
      idct->data[i * dct->d2 + l] = dct->data[i+ dct->d2 * l];
  return 0;
}

int get_sldtable(char *name, int nparams, muT *out)
{
  FILE *fp;
  int i, j, k, N, K ;
  muT tmp;
  
  if (!(fp=fopen(name,"r"))) {
    printf("%s can not be opened for reading\n",name);
    exit(1);
  }
  fscanf(fp,"<%d> <%d>\n",&N,&K);
  if(K < nparams) {
    fprintf(stderr,"No of output coefficients greater than the matrix size\n");
    return -1;
  }
  printf("loading sld matrix with %d %d \n",N,K);
  if(mx_set_size(N, K, &tmp)) {
    return -1;
  }
  if(mx_set_size(N, nparams, out)) {
    return -1;
  }
  
  k=0;
  for(i=0 ; i < tmp.d1; i++) {
    for(j=0; j < tmp.d2; j++) {
      if(fscanf(fp,"%f ",&(tmp.data[i*tmp.d2+j]))==EOF) {
	printf("%s contains insufficient data \n",name);
	exit(1);
      }
      k++;
    }
  }
  //added to have parameters fewer than 15
  for(i=0 ; i < out->d1; i++)
    memcpy(out->data+i*out->d2,tmp.data+i*tmp.d2,out->d2*sizeof(float));
  
  fclose(fp);
  printf("loaded %d elements\n",k);
  return(0);
}

//added Fend on 07/15/2001 by sunil
int get_mel_filterbank(int nFilters, int nPoints, float Fs, float Fstart, float Fend, muT* out) {
  float mf_start, mf_sd2, Fc;
  int *cbin, i, k;
  
  if((cbin = (int *) malloc((nPoints + 2) * sizeof(int))) == NULL) {
    return -1;
  }
  if(mx_set_size(nPoints / 2 + 1, nFilters, out)) {
    free(cbin);
    return -1;
  }
  mf_start = Melf(Fstart);
  //mf_sd2 = Melf(Fs / 2.0);
  mf_sd2 = Melf(Fend);
  
  cbin[0] = ROUND(Fstart/Fs*nPoints);
  for(i = 1; i <= nFilters; i++) {
    Fc = InvMelf(mf_start + (mf_sd2 - mf_start) / (nFilters + 1) * i);
    cbin[i] = ROUND(Fc / Fs * nPoints);
  }
  //cbin[nFilters + 1] = ROUND(nPoints / 2);
  cbin[nFilters + 1] = ROUND(Fend/Fs*nPoints );
  
  for(k = 1; k <= nFilters; k++) {
    for(i = cbin[k - 1]; i <= cbin[k]; i++) {
      out->data[out->d2 * i + (k-1)] = (i - cbin[k - 1]) / (float) (cbin[k] - cbin[k-1]);
    }
    for(i = cbin[k] + 1; i <= cbin[k + 1]; i++) {
      out->data[out->d2 * i + (k-1)] = 1 - (i - cbin[k]) / (float) (cbin[k + 1] - cbin[k]);
    }
  }
  free(cbin);
  return 0;
}

/*******************************************************************************************/

int temp_filter(muT *in, muT *filter, muT *out)
{
  int i, j, k, l, delay;
  muT tmx;
  assert(in->d2 == filter->d1);
  
  delay=(filter->d2 - 1) / 2;
  if(mx_set_size(in->d1 + delay * 2, in->d2, &tmx)) {
    return -1;
  }
  if(mx_set_size(in->d1, in->d2, out)) {
    mx_release(&tmx);
    return -1;
  }
  /* fill the beginning */
  for(i = 0; i < delay; i++) {
    for(k = 0; k < tmx.d2; k++) {
      tmx.data[(i + 1) * tmx.d2 - k - 1] = in->data[i * in->d2 + k];
    }
  }
  /* copy the data */
  memcpy(tmx.data + tmx.d2 * delay, in->data, in->d1 * in->d2 * sizeof(float));
  /* fill the end  */
  for(i = in->d1 + delay; i < tmx.d1; i++) {
    for(k = 0; k < tmx.d2; k++) {
      j=i-2*delay;
      tmx.data[i * tmx.d2 + k] = in->data[(j + 1) * in->d2 - k - 1];
    }
  }
  for(k = 0; k < tmx.d2; k++) {
    for(i = 0; i < in->d1; i++) {
      out->data[i * out->d2 + k] = 0.0;
      for(j = 0, l = i; j < filter->d2; j++, l++) {
        out->data[i * out->d2 + k] += tmx.data[l * tmx.d2 + k]
	  * filter->data[k * filter->d2 + j];
      }
    }
  }
  mx_release(&tmx);
  return 0;
}

/* dupont@icsi.berkeley.edu 16 Apr 2001: add mean and variance estimation latency argument */
int online_normalization(muT *in, float alpha, float bias, muT *out, int lat)
{
  float mean[FEAT_VEC_SIZE] = { 23.815484, 1.037240, -0.382422, -0.203596, -0.557713, -0.051042,
                                -0.208684, 0.073762, -0.100447,  0.007481, -0.062511, -0.001553,
                                -0.123591, -0.006837, -0.140246 };
  float std_sq[FEAT_VEC_SIZE] = { SQR(5.957326), SQR(1.350153), SQR(0.992368), SQR(0.685526),
                                  SQR(0.834512), SQR(0.545422), SQR(0.478728), SQR(0.476498),
                                  SQR(0.422000), SQR(0.417962), SQR(0.351819), SQR(0.361830),
                                  SQR(0.323899), SQR(0.322991), SQR(0.287901) };
  long i, j;
  int latency = lat;
  
  if(mx_set_size(in->d1, in->d2, out)) {
    return -1;
  }
  for(i = 0; i < in->d1; i++) {
    for(j = 0; j < in->d2; j++) {
      mean[j] += (in->data[i * in->d2 + j] - mean[j]) * alpha;
      std_sq[j] += (SQR(in->data[i * in->d2 + j] - mean[j]) - std_sq[j]) * alpha;
      if (i>=latency) {
	out->data[(i-latency)*out->d2+j] = (in->data[(i-latency)*in->d2+j] - mean[j]) / (sqrt(std_sq[j]) + bias);
      }
    }
  }
  if (latency) {
    for(i = 0; i < latency; i++) {
      for(j = 0; j < in->d2; j++) {
	out->data[(in->d1-i-1)*out->d2+j] = (in->data[(in->d1-i-1)*in->d2+j] - mean[j]) / (sqrt(std_sq[j]) + bias);
      }
    }
  }
  return 0;
}

int upsample2(muT *in, muT *out)
{
  long i, j;
  
  if(mx_set_size(in->d1 * 2 - 1, in->d2, out)) {
    return -1;
  }
  for(i = 0; i < in->d1-1; i++) {
    memcpy(out->data + i * 2 * out->d2, in->data + i * out->d2, in->d2 * sizeof(float));
    for(j = 0; j < in->d2; j++) {
      out->data[(i*2+1) * out->d2 + j]  = (in->data[i * in->d2 + j]
					   +  in->data[(i+1) * in->d2 + j]) / 2;
    }
    memcpy(out->data + (in->d1-1) * 2 * out->d2,
           in->data + (in->d1-1) * out->d2, in->d2 * sizeof(float));
  }
  return 0;
}

int downsample2(muT *in, int out_dim, muT *out)
{
  long i;
  
  if(mx_set_size((in->d1 + 1) / 2, out_dim, out)) {
    return -1;
  }
  for(i = 0; i < out->d1; i++) {
    memcpy(out->data + i * out->d2, in->data + i * 2 * in->d2, out->d2 * sizeof(float));
  }
  return 0;
}

long smooth_sil_flag(int range, int point, char *slifea, long length)
{
  char *tmp;
  long i, totalNSilFrms = 0, nSilFrms = 0;
  int padlen = (range - 1) / 2;
  if((tmp = (char *) malloc(length + 2 * padlen)) == NULL) {
    return -1;
  }
  memset(tmp, 1, padlen);
  memcpy(tmp+padlen, slifea, length);
  memset(tmp+padlen+length, 1, padlen);
  
  for(i = 0; i < range-1; i++) {
    nSilFrms += tmp[i];
  }
  for(i = 0; i < length; i++) {
    nSilFrms += tmp[i+range - 1];
    totalNSilFrms += slifea[i] = (char) (nSilFrms >= (range - point));
    nSilFrms -= tmp[i];
  }
  free(tmp);
  
  /* dupont@icsi.berkeley.edu 16 Apr 2001: to have at least some frames that are not silence */
  /* if (totalNSilFrms == length) {
     totalNSilFrms = length-range;
     for(i = 0; i < range; i++) {
     slifea[i] = 0;
     }
     } */
  if (totalNSilFrms == length) {
    totalNSilFrms = length-range;
    for(i = length-range; i < length; i++) {
      slifea[i] = 0;
    }
  }
  /* drop first 10 frames anyway */
  /*if (range != 0) {
    for(i = 0; i < 10; i++) {
    if (slifea[i] = 0) { totalNSilFrms++; }
    slifea[i] = 1;
    }
    }*/
  
  return length - totalNSilFrms;
  
}

//Added on 07/31/2001 by sunil
char *sil_prob_to_flag(long length, float *silfea, float vad_threshold, int vad_rank_range, int vad_rank_point)
{
  
  char *sil_flag;
  long i,nSpeechFrames = 0;
  
  //Smoothing silence probability
  if ((sil_flag = (char*)malloc(length)) == NULL) {
    return 0;
  }
  for(i = 0; i < length; i++) {
    sil_flag[i] = (char)(silfea[i] > vad_threshold);
  }
  
  if (vad_rank_range>0) {
    nSpeechFrames = smooth_sil_flag(vad_rank_range, vad_rank_point, sil_flag, length);
  }
  
  return sil_flag;
}

/*int cut_sil_parts(muT *in, float *silfea, muT *out)
  {
  char *sil_flag;
  long i, j, nSpeechFrames = 0;
  
  if((sil_flag = (char *) malloc(in->d1)) == NULL) {
  return -1;
  }
  for(i = 0; i < in->d1; i++) {
  sil_flag[i] = silfea[i] > VAD_TRESHOLD;
  }
  nSpeechFrames = smooth_sil_flag(VAD_RANK_RANGE, VAD_RANK_POINT, sil_flag, in->d1);
  if(mx_set_size(nSpeechFrames, in->d2, out)) {
  free(sil_flag);
  return -1;
  }
  for(i = 0, j = 0; i < in->d1; i++) {
  //    printf("%d %f\n", !sil_flag[i], in->data[i * in->d2]);
  if(!sil_flag[i]) {
  memcpy(out->data + j * out->d2, in->data + i * in->d2, in->d2 * sizeof(float));
  ++j;
  }
  }
  free(sil_flag);
  return 0;
  }*/

int framedropping(muT *in, char *sil_flag_ds2, int down_up,muT *out)
{
  
  char *sil_flag;
  long i, j, nSpeechFrames = 0;
  
  if((sil_flag = (char *) malloc(in->d1+1)) == NULL) {
    return -1;
  }
  
  if (down_up) {
    //Upsample the silence flag
    for(i = 0; i < (in->d1+1)/2 ; i++) {
      sil_flag[2*i]=sil_flag_ds2[i];
      sil_flag[2*i+1]= sil_flag_ds2[i];
    }
  } else {
    memcpy(sil_flag,sil_flag_ds2,in->d1);
  }
  
  for(i = 0; i < in->d1; i++) {
    nSpeechFrames += sil_flag[i];
  }
  nSpeechFrames = in->d1-nSpeechFrames;
  
  if(mx_set_size(nSpeechFrames, in->d2, out)) {
    free(sil_flag);
    return -1;
  }
  
  for(i = 0, j = 0; i < in->d1 && j < nSpeechFrames; i++) {
    if(!sil_flag[i]) {
      memcpy(out->data + j * out->d2, in->data + i * in->d2, in->d2 * sizeof(float));
      ++j;
    }
  }
  
  free(sil_flag);
  return 0;
}

int cut_sil_parts_orig_ds2_silflag(muT *in, char *sil_flag_ds2, muT *out)
{
  
  char *sil_flag;
  long i, j, nSpeechFrames = 0;
  
  if((sil_flag = (char *) malloc(in->d1)) == NULL) {
    return -1;
  }
  
  //Upsample the silence flag
  for(i = 0; i < (in->d1+1)/2 ; i++) {
    sil_flag[2*i]=sil_flag_ds2[i];
    sil_flag[2*i+1]= sil_flag_ds2[i];
  }
  
  for(i = 0; i < in->d1; i++) {
    nSpeechFrames += sil_flag[i];
  }
  nSpeechFrames = in->d1-nSpeechFrames;
  
  if(mx_set_size(nSpeechFrames, in->d2, out)) {
    free(sil_flag);
    return -1;
  }
  
  for(i = 0, j = 0; i < in->d1 && j < nSpeechFrames; i++) {
    if(!sil_flag[i]) {
      memcpy(out->data + j * out->d2, in->data + i * in->d2, in->d2 * sizeof(float));
      ++j;
    }
  }
  
  free(sil_flag);
  return 0;
}


/*******************************************************************************************/

/*float who[] = {-2.130615, -1.826538, 0.135376, -2.932373, -1.150024, 2.651001, -1.575928, -0.032104, -0.107300, -0.198364, -2.484131,
  2.331177, 0.240723, 2.337524, 2.661743, 5.097004, 1.991821, 1.747925, -0.730347, 1.911621, 0.684814, -2.793213,
  1.527100, -0.120972, 0.149170, 0.199341, 1.975952, -2.816528, -0.586670, -2.773315, -2.775513, 4.866590};
  float wih[] = {-1.696899, -0.553101, -0.040894, 0.411377, 0.450439, 0.294556, 5.659857, 0.715942, -0.212036, 1.208374, 0.012573,
  2.086060, 0.357788, 1.322334, -2.350708, 0.785767, 2.566895, 1.210938, 3.257446, 3.999390, -1.780486, 0.987305,
  -0.345337, 2.617432, 0.484497, -3.551636, 0.496338, 0.145179, -0.279663, -3.537109, 0.088989, -2.980713, 0.541260,
  -3.997070, 5.363252, 3.763550, -2.002441, 1.088379, 0.034058, -3.987915, 1.337891, 1.710441, -0.866089, 0.229126,
  1.010132, -0.060547, 0.782349, 0.088013, 5.014454, 0.087158, 0.038208, 0.067017, 0.130737, -0.009521, 0.092773,
  4.100571, -0.025269, 0.006226, 0.041138, 0.003906, -0.067261, 0.080811, 4.055887, 0.096436, 0.005493, 0.192749,
  0.037476, 0.138428, 0.017334, 3.927449, 0.478516, -0.333496, 0.124023, 0.030396, -1.102417, 0.871948, 6.408793,
  3.993652, 2.438721, 1.240356, 0.352661, -3.798584, -1.132935, 1.041499, 2.920166, 2.028320, 3.136108, 0.901367,
  2.111816, 2.610596, 3.457172, -2.600098, 0.275757, 0.713501, -0.471680, 3.154419, -1.198975, 1.227792, -3.963257,
  -1.110474, 2.421875, 0.340088, 3.921387, 2.897461, 1.522875};
  
  float bias[] = {-68.916397, -2.685470};
  float scale[] = {0.060627, 0.218782};
  
  #define INP_LAYER_SIZE  (INP_MEMORY * INP_WIDTH)*/

#define DEFBIASUNITVAL (-1.0)

/*int ninp = INP_LAYER_SIZE;
  int nhid = HID_LAYER_SIZE;
  int nout = OUT_LAYER_SIZE;
  int width = INP_WIDTH;
  int cw = INP_MEMORY;
  float bval = DEFBIASUNITVAL;
  char *wfile = NULL;
  char *nfile = NULL;*/

int onlftrs_read(FILE *fp, float *buf, int nftrs)
{
  int i;
  for(i = 0; i < nftrs; i++) {
    if(fscanf(fp, "%f", buf+i) != 1) {
      if(i == 0 && feof(fp)) {
        return 1;
      }
      return -1;
    }
  }
  return 0;
}

void norm_ftrs(const float *src, float *dst, int width,
              const float *bias, const float *scale)
{   /* normalize according to bias and scale read from file */
  int i;
  for (i = 0; i < width; ++i) {
    *dst++ = (*scale++) * (*src++ + *bias++);
  }
}

void * lmalloc(int bytes)
{
  void *tmp;
  if ((tmp = (void*)malloc(bytes)) == NULL) {
    fprintf(stderr,"Can't allocate memory\n");
    exit(-1);
  }
  return tmp;
}

void bias_matrix(float *mat,int M,int N, float bval)
{  /* setup the 'bias' values in a matrix to the defined constant value.
      Bias units are are at the end of each row i.e. every N values,
      where N is the number of columns. */
  float *p = &mat[N]-1;
  while (M--) {
    *p = bval;
    p += N;
  }
}

#define LOGISTIC(x) (1.0/(1.0+exp(-(x))))

void applysigmoid(float *mat, int N, int stride)
{
  float *endp = &mat[stride];
  while (mat != endp) {
    float *rmat = mat;
    float *rendp = &mat[N];
    while (rmat != rendp) {
      *rmat = LOGISTIC(*rmat);
      rmat++;
    }
    mat += stride;
  }
}

void applysoftmax(float *matrix, int N, int stride)
{
  float *endp = &matrix[stride];
  float *row = matrix;
  while (row != endp) {
    float *rendp = &row[N];
    float *col = row;
    float sum = 0;
    while (col != rendp) {
      *col = exp(*col);
      sum += *col++;
    }
    sum = 1.0/sum;
    col = row;
    while (col != rendp) {
      *col++ *= sum;
    }
    row += stride;
  }
}

void bunch_layer_fwd_prop(int insize, int outsize, float *w, float *in, float *out, int addstride) {
  
  int i;
  int bunch = BUNCH_SIZE;
  
  for(i = 0; i < outsize; i++) {
    out[i] = 0;
  }
  
  mul_mfmft_mf_beta0(bunch, insize+1, outsize,
		     in, w, out,
		     insize+1, insize+1, outsize+addstride); 
  
}

void layer_fwd_prop(int insize, int outsize, float *w, float *in, float *out) {
  int i, j;
  for(i = 0; i < outsize; i++) {
    out[i] = 0;
    for(j = 0; j < insize+1; j++) {
      out[i]+=w[(insize+1)*i+j]*in[j];
    }
  }
}

void nn_forw_prop(int Ni, int Nh, int No, float *who, float *wih, float *ipats, float *hpats, float *opats, int lin)
{
  layer_fwd_prop(Ni, Nh, wih, ipats, hpats);
  applysigmoid(hpats,Nh,Nh+1);
  layer_fwd_prop(Nh, No, who, hpats, opats);
  if (!lin) {
    applysoftmax(opats,No,No);
  }
}

void bunch_nn_forw_prop(int Ni, int Nh, int No, float *who, float *wih, float *ipats, float *hpats, float *opats, int lin)
{
  int i;
  int bunch = BUNCH_SIZE;
  bunch_layer_fwd_prop(Ni, Nh, wih, ipats, hpats,1);
  for (i=0;i<bunch;i++) {
    applysigmoid(hpats+(Nh+1)*i,Nh,Nh+1);
  }
  bunch_layer_fwd_prop(Nh, No, who, hpats, opats,0);
  if (!lin) {
    for (i=0;i<bunch;i++) {
      applysoftmax(opats+No*i,No,No);
    }
  }
}


void init_mlp_forward(MLPParams *params, char *mlp_filename, char *feat_norm_filename)
{
  
  float *pos;
  int i,j;
  long dummy;
  FILE *fp;
  int nih;
  int nho;
  
  fp = fopen(feat_norm_filename,"r");
  if ( fp == NULL ) {
    fprintf(stderr,"Error in opening file -> %s\n",feat_norm_filename);
  }
  fscanf(fp,"vec %d\n",&params->width);
  params->bias = malloc(params->width*sizeof(float));
  params->scale = malloc(params->width*sizeof(float));
  for (i=0;i<params->width;i++) {
    fscanf(fp,"%f\n",&params->bias[i]);
    params->bias[i] = -params->bias[i];
    //params->bias[i] = 0;
  }
  fscanf(fp,"vec %ld\n",&dummy);
  for (i=0;i<params->width;i++) {
    fscanf(fp,"%f\n",&params->scale[i]);
    //params->scale[i] = 1;
  }
  fclose(fp);
  dummy=0;
  
  fp = fopen(mlp_filename,"r");
  fscanf(fp,"%d %d %d\n",&params->ninp,&params->nhid,&params->nout);
  //fprintf (stderr,"%d,%d,%d\n",params->ninp,params->nhid,params->nout);
  params->cw = params->ninp/params->width;
  params->bval = DEFBIASUNITVAL;
  nih = (params->ninp+1)*params->nhid;
  nho = params->nout*(params->nhid+1);
  params->wih = malloc(nih*sizeof(float));
  params->who = malloc(nho*sizeof(float));
  fscanf(fp,"weigvec %ld\n",&dummy);
  assert(dummy==params->ninp*params->nhid);
  //for (j=0;j<params->nhid;j++) {
  pos = &params->wih[0];
  for (i=0;i<params->nhid;i++) {
    for (j=0;j<params->ninp;j++) {
      fscanf(fp,"%f\n",pos++);
    }
    pos++;
  }
  fscanf(fp,"weigvec %ld\n",&dummy);
  assert(dummy==params->nout*params->nhid);
  //for (j=0;j<params->nout;j++) {
  pos = &params->who[0];
  for (i=0;i<params->nout;i++) {
    for (j=0;j<params->nhid;j++) {
      fscanf(fp,"%f\n",pos++);
    }
    pos++;
  }
  i = params->ninp;
  fscanf(fp,"biasvec %ld\n",&dummy);
  pos = &params->wih[0]+params->ninp;
  for (j=0;j<params->nhid;j++) {
    fscanf(fp,"%f\n",pos);
    *pos = - *pos;
    pos += params->ninp+1;
  }
  i = params->nhid;
  pos = &params->who[0]+params->nhid;
  fscanf(fp,"biasvec %ld\n",&dummy);
  for (j=0;j<params->nout;j++) {
    fscanf(fp,"%f\n",pos);
    *pos = - *pos;
    pos += params->nhid+1;
  }
  fclose(fp);
//  fprintf(stderr,"End init MLP\n");
}

int mlp_forward_linear(muT *feat, muT *out, MLPParams params, int feat_length, char *silflag)
{
  
  float *ipats, *hpats, *opats, *silfea, *ipp;
  int nread;
  
//  fprintf(stderr,"Begin forward MLP\n");
  
  mx_set_size((feat_length+1)/2, params.nout, out);
  
  ipats = lmalloc((params.ninp+1)*sizeof(float));
  hpats = lmalloc((params.nhid+1)*sizeof(float));
  opats = lmalloc(params.nout*sizeof(float));
  silfea = lmalloc((feat_length)*sizeof(float));
  
  bias_matrix(hpats,1,params.nhid+1,params.bval);
  bias_matrix(ipats,1,params.ninp+1,params.bval);
  
  memset(opats,0,params.nout*sizeof(float));
  
  ipp = ipats + (params.ninp+1) - (params.cw-1)*params.width - 1;
  for(nread = 0; nread < params.cw - 1; nread++) {
    memcpy(ipp, &feat->data[0], sizeof(float)*params.width);
    norm_ftrs(ipp, ipp, params.width, params.bias, params.scale);
    ipp += params.width;
  }
  for(nread = 0; nread < feat_length+(params.cw-1)/2; nread++) {
    memmove(ipats, ipats + (params.ninp+1) - (params.cw-1)*params.width - 1,
	    (params.cw-1)*params.width*sizeof(float));
    memcpy(ipats+(params.cw-1)*params.width, &feat->data[nread < feat_length ? nread*params.width : (feat_length-1)*params.width],
	   params.width*sizeof(float));
    norm_ftrs(ipats+(params.cw-1)*params.width, ipats+(params.cw-1)*params.width,
	      params.width, params.bias, params.scale);
    if (nread>=(params.cw-1)/2) {
      if (!((nread-(params.cw-1)/2)%2)) {
	if (!silflag[(nread-(params.cw-1)/2)/2]) {
	  nn_forw_prop(params.ninp, params.nhid, params.nout, params.who, params.wih, ipats, hpats, opats, 1);
	  memcpy(&out->data[((nread-(params.cw-1)/2))/2*params.nout],opats,params.nout*sizeof(float));
	}
	else {
	  memcpy(&out->data[((nread-(params.cw-1)/2))/2*params.nout],opats,params.nout*sizeof(float));
	}
      }
    }
  }
  free(ipats); free(hpats); free(opats);
//  fprintf(stderr,"End forward MLP\n");
  return 1;
}


int bunch_mlp_forward_linear(muT *feat, muT *out, MLPParams params, int feat_length, char *silflag)
{
  int speech_present;
  
  int count, i, center, current;
  int bunch = BUNCH_SIZE;
  int local_bunch;
  
  float *ipats, *hpats, *opats, *silfea;
  float *normfea;
  int nread;
  
//  fprintf(stderr,"Begin forward MLP\n");
  
  //mx_set_size((feat_length+1)/2, params.nout, out);
  mx_set_size(feat_length, params.nout, out);
  
  ipats = lmalloc(bunch*(params.ninp+1)*sizeof(float));
  hpats = lmalloc(bunch*(params.nhid+1)*sizeof(float));
  opats = lmalloc(bunch*params.nout*sizeof(float));
  silfea = lmalloc((feat_length)*sizeof(float));
  normfea = lmalloc((feat_length*params.width)*sizeof(float));
  
  bias_matrix(hpats,bunch,params.nhid+1,params.bval);
  bias_matrix(ipats,bunch,params.ninp+1,params.bval);
  
  memset(opats,0,bunch*params.nout*sizeof(float));
  
  // normalize features
  for(nread = 0; nread < feat_length; nread++) {
    norm_ftrs(&feat->data[nread*params.width],normfea+nread*params.width,
	      params.width, params.bias, params.scale);
  }
  
  // build bunches
  nread = 0;
  while (nread < feat_length) {
    
    //local_bunch=min(bunch,(feat_length-nread+1)/2);
    local_bunch=min(bunch,feat_length-nread);
    
    speech_present = 0;
    if (local_bunch>0) {
      for (count=0;count<local_bunch;count++) {
	
	//center=2*count+nread;
	center=count+nread;
	
	for(i=0; i<params.cw; i++) {
	  current=center+i-(params.cw-1)/2;
	  memcpy(ipats+count*(params.width*params.cw+1)+i*params.width, normfea+max(min(current,feat_length-1),0)*params.width,
		 params.width*sizeof(float));
	}
	
	//if ( !silflag[center/2] ) {
	if ( !silflag[center] ) {
	  speech_present = 1;
	}
	speech_present = 1;
      }
      if (speech_present) {
	bunch_nn_forw_prop(params.ninp, params.nhid, params.nout, params.who, params.wih, ipats, hpats, opats, 1);
      }
      
      //memcpy(&out->data[nread/2*params.nout],opats,local_bunch*params.nout*sizeof(float));
      memcpy(&out->data[nread*params.nout],opats,local_bunch*params.nout*sizeof(float));
    }
    
    //nread+=local_bunch*2;
    nread+=local_bunch;
    
  }
  
  free(ipats); free(hpats); free(opats);
//  fprintf(stderr,"End forward MLP\n");
  return 1;
}

int bunch_mlp_forward(muT *feat, muT *out, MLPParams params, int feat_length, char *silflag, int lin)
{
  int speech_present;
  
  int count, i, center, current;
  int bunch = BUNCH_SIZE;
  int local_bunch;
  
  float *ipats, *hpats, *opats, *silfea;
  float *normfea;
  int  nread;
  
//  fprintf(stderr,"Begin forward MLP\n");
  
  mx_set_size(feat_length, params.nout, out);
  
  ipats = lmalloc(bunch*(params.ninp+1)*sizeof(float));
  hpats = lmalloc(bunch*(params.nhid+1)*sizeof(float));
  opats = lmalloc(bunch*params.nout*sizeof(float));
  silfea = lmalloc((feat_length)*sizeof(float));
  normfea = lmalloc((feat_length*params.width)*sizeof(float));
  
  bias_matrix(hpats,bunch,params.nhid+1,params.bval);
  bias_matrix(ipats,bunch,params.ninp+1,params.bval);
  
  memset(opats,0,bunch*params.nout*sizeof(float));
  
  // normalize features
  for(nread = 0; nread < feat_length; nread++) {
    norm_ftrs(&feat->data[nread*params.width],normfea+nread*params.width,
	      params.width, params.bias, params.scale);
  }
  
  // build bunches
  nread = 0;
  while (nread < feat_length) {
    local_bunch=min(bunch,(feat_length-nread));
    speech_present = 0;
    if (local_bunch>0) {
      for (count=0;count<local_bunch;count++) {
	center=count+nread;
	for(i=0; i<params.cw; i++) {
	  current=center+i-(params.cw-1)/2;
	  memcpy(ipats+count*(params.width*params.cw+1)+i*params.width, normfea+max(min(current,feat_length-1),0)*params.width,
		 params.width*sizeof(float));
	}
	if (silflag) {
	  if ( !silflag[center] ) {
	    speech_present = 1;
	  }
	}
	else {
	  speech_present = 1;
	}
	speech_present = 1;
      }
      if (speech_present) {
	bunch_nn_forw_prop(params.ninp, params.nhid, params.nout, params.who, params.wih, ipats, hpats, opats, lin);
      }
      memcpy(&out->data[nread*params.nout],opats,local_bunch*params.nout*sizeof(float));
    }
    nread+=local_bunch;
  }
  
  free(ipats); free(hpats); free(opats);
//  fprintf(stderr,"End forward MLP\n");
  return 1;
}



/*******  G. D. Bergland and M. T. DoLan (Fortran), Paul Kube (C) *********/

/*
 * Discrete Fourier analysis routine
 * from IEEE Programs for Digital Signal Processing
 * G. D. Bergland and M. T. DoLan, original authors
 * Translated from the FORTRAN with some changes by Paul Kube
 *
 * Modified to return the power spectrum by Chuck Wooters
 *
 * Modified again by Tony Robinson (ajr@eng.cam.ac.uk) Dec 92
 *
 * Slight naming mods by N. Morgan, July 1993
 *      (fft_chuck -> fft_pow)
 *      ( calling args long ll -> long winlength)
 *      (long m -> long log2length)
 */

/*
 * This routine replaces the float vector b
 * of length n with its finite discrete fourier transform.
 * DC term is returned in b[0];
 * n/2th harmonic float part in b[1].
 * jth harmonic is returned as complex number stored as
 * b[2*j] + i b[2*j + 1]
 * (i.e., remaining coefficients are as a DPCOMPLEX vector).
 *
 */

#define PI      (3.1415926535897932)
#define PI8     (0.392699081698724)     /* PI / 8.0 */
#define RT2     (1.4142135623731)       /* sqrt(2.0) */
#define IRT2    (0.707106781186548)     /* 1.0/sqrt(2.0) */
#define signum(i) (i < 0 ? -1 : i == 0 ? 0 : 1)

static void FR2TR(int, float*, float*);
static void FR4TR(int, int, float*, float*, float*, float*);
static void FORD1(int, float*);
static void FORD2(int, float*);
static int  fastlog2(int);

int mu_FAST(float *b,int n)
{
  float fn;
  int i, in, nn, n2pow, n4pow, nthpo;
  
  n2pow = fastlog2(n);
  if(n2pow <= 0) return (0);
  nthpo = n;
  fn = nthpo;
  n4pow = n2pow / 2;
  
  /* radix 2 iteration required; do it now */
  if(n2pow % 2) {
    nn = 2;
    in = n / nn;
    FR2TR(in, b, b + in);
  }
  else nn = 1;
  
  /* perform radix 4 iterations */
  for(i = 1; i <= n4pow; i++) {
    nn *= 4;
    in = n / nn;
    FR4TR(in, nn, b, b + in, b + 2 * in, b + 3 * in);
  }
  
  /* perform inplace reordering */
  FORD1(n2pow, b);
  FORD2(n2pow, b);
  
  /* take conjugates */
  for(i = 3; i < n; i += 2) b[i] = -b[i];
  
  return 1;
}

/*
 * radix 2 subroutine
 */
static void FR2TR(int in,float *b0,float *b1)
{
  int k;
  float t;
  for(k = 0; k < in; k++) {
    t = b0[k] + b1[k];
    b1[k] = b0[k] - b1[k];
    b0[k] = t;
  }
}

/*
 * radix 4 subroutine
 */
static void FR4TR(int in,int nn,float *b0,float *b1,float *b2,float *b3)
{
  float arg, piovn, th2;
  float *b4 = b0, *b5 = b1, *b6 = b2, *b7 = b3;
  float t0, t1, t2, t3, t4, t5, t6, t7;
  float r1, r5, pr, pi;
  float c1, c2, c3, s1, s2, s3;
  
  int j, k, jj, kk, jthet, jlast, ji, jl, jr, int4;
  int L[16], L1, L2, L3, L4, L5, L6, L7, L8, L9, L10, L11, L12, L13, L14, L15;
  int j0, j1, j2, j3, j4, j5, j6, j7, j8, j9, j10, j11, j12, j13, j14;
  int k0, kl;
  
  L[1] = nn / 4;
  for(k = 2; k < 16; k++) {  /* set up L's */
    switch (signum(L[k-1] - 2)) {
    case -1:
      L[k-1]=2;
    case 0:
      L[k]=2;
      break;
    case 1:
      L[k]=L[k-1]/2;
    }
  }
  
  L15=L[1]; L14=L[2]; L13=L[3]; L12=L[4]; L11=L[5]; L10=L[6]; L9=L[7];
  L8=L[8];  L7=L[9];  L6=L[10]; L5=L[11]; L4=L[12]; L3=L[13]; L2=L[14];
  L1=L[15];
  
  piovn = PI / nn;
  ji=3;
  jl=2;
  jr=2;
  
  for(j1=2;j1<=L1;j1+=2)
    for(j2=j1;j2<=L2;j2+=L1)
      for(j3=j2;j3<=L3;j3+=L2)
	for(j4=j3;j4<=L4;j4+=L3)
	  for(j5=j4;j5<=L5;j5+=L4)
	    for(j6=j5;j6<=L6;j6+=L5)
	      for(j7=j6;j7<=L7;j7+=L6)
		for(j8=j7;j8<=L8;j8+=L7)
		  for(j9=j8;j9<=L9;j9+=L8)
		    for(j10=j9;j10<=L10;j10+=L9)
		      for(j11=j10;j11<=L11;j11+=L10)
			for(j12=j11;j12<=L12;j12+=L11)
			  for(j13=j12;j13<=L13;j13+=L12)
			    for(j14=j13;j14<=L14;j14+=L13)
			      for(jthet=j14;jthet<=L15;jthet+=L14) {
				th2 = jthet - 2;
				if(th2<=0.0) {
				  for(k=0;k<in;k++) {
				    t0 = b0[k] + b2[k];
				    t1 = b1[k] + b3[k];
				    b2[k] = b0[k] - b2[k];
				    b3[k] = b1[k] - b3[k];
				    b0[k] = t0 + t1;
				    b1[k] = t0 - t1;
				  }
				  if(nn-4>0) {
				    k0 = in*4 + 1;
				    kl = k0 + in - 1;
				    for (k=k0;k<=kl;k++) {
				      kk = k-1;
				      pr = IRT2 * (b1[kk]-b3[kk]);
				      pi = IRT2 * (b1[kk]+b3[kk]);
				      b3[kk] = b2[kk] + pi;
				      b1[kk] = pi - b2[kk];
				      b2[kk] = b0[kk] - pr;
				      b0[kk] = b0[kk] + pr;
				    }
				  }
				} else {
				  arg = th2*piovn;
				  c1 = cos(arg);
				  s1 = sin(arg);
				  c2 = c1*c1 - s1*s1;
				  s2 = c1*s1 + c1*s1;
				  c3 = c1*c2 - s1*s2;
				  s3 = c2*s1 + s2*c1;
				  
				  int4 = in*4;
				  j0=jr*int4 + 1;
				  k0=ji*int4 + 1;
				  jlast = j0+in-1;
				  for(j=j0;j<=jlast;j++) {
				    k = k0 + j - j0;
				    kk = k-1; jj = j-1;
				    r1 = b1[jj]*c1 - b5[kk]*s1;
				    r5 = b1[jj]*s1 + b5[kk]*c1;
				    t2 = b2[jj]*c2 - b6[kk]*s2;
				    t6 = b2[jj]*s2 + b6[kk]*c2;
				    t3 = b3[jj]*c3 - b7[kk]*s3;
				    t7 = b3[jj]*s3 + b7[kk]*c3;
				    t0 = b0[jj] + t2;
				    t4 = b4[kk] + t6;
				    t2 = b0[jj] - t2;
				    t6 = b4[kk] - t6;
				    t1 = r1 + t3;
				    t5 = r5 + t7;
				    t3 = r1 - t3;
				    t7 = r5 - t7;
				    b0[jj] = t0 + t1;
				    b7[kk] = t4 + t5;
				    b6[kk] = t0 - t1;
				    b1[jj] = t5 - t4;
				    b2[jj] = t2 - t7;
				    b5[kk] = t6 + t3;
				    b4[kk] = t2 + t7;
				    b3[jj] = t3 - t6;
				  }
				  jr += 2;
				  ji -= 2;
				  if(ji-jl <= 0) {
				    ji = 2*jr - 1;
				    jl = jr;
				  }
				}
			      }
}

/*
 * an inplace reordering subroutine
 */
static void FORD1(int m,float *b)
{
  int j, k = 4, kl = 2, n = 0x1 << m;
  float t;
  
  for(j = 4; j <= n; j += 2) {
    if(k - j>0) {
      t = b[j-1];
      b[j - 1] = b[k - 1];
      b[k - 1] = t;
    }
    k -= 2;
    if(k - kl <= 0) {
      k = 2*j;
      kl = j;
    }
  }	
}

/*
 *  the other inplace reordering subroutine
 */
static void FORD2(int m,float *b)
{
  float t;
  
  int n = 0x1<<m, k, ij, ji, ij1, ji1;
  
  int l[16], l1, l2, l3, l4, l5, l6, l7, l8, l9, l10, l11, l12, l13, l14, l15;
  int j1, j2, j3, j4, j5, j6, j7, j8, j9, j10, j11, j12, j13, j14;
  
  
  l[1] = n;
  for(k=2;k<=m;k++) l[k]=l[k-1]/2;
  for(k=m;k<=14;k++) l[k+1]=2;
  
  l15=l[1];l14=l[2];l13=l[3];l12=l[4];l11=l[5];l10=l[6];l9=l[7];
  l8=l[8];l7=l[9];l6=l[10];l5=l[11];l4=l[12];l3=l[13];l2=l[14];l1=l[15];
  
  ij = 2;
  
  for(j1=2;j1<=l1;j1+=2)
    for(j2=j1;j2<=l2;j2+=l1)
      for(j3=j2;j3<=l3;j3+=l2)
	for(j4=j3;j4<=l4;j4+=l3)
	  for(j5=j4;j5<=l5;j5+=l4)
	    for(j6=j5;j6<=l6;j6+=l5)
	      for(j7=j6;j7<=l7;j7+=l6)
		for(j8=j7;j8<=l8;j8+=l7)
		  for(j9=j8;j9<=l9;j9+=l8)
		    for(j10=j9;j10<=l10;j10+=l9)
		      for(j11=j10;j11<=l11;j11+=l10)
			for(j12=j11;j12<=l12;j12+=l11)
			  for(j13=j12;j13<=l13;j13+=l12)
			    for(j14=j13;j14<=l14;j14+=l13)
			      for(ji=j14;ji<=l15;ji+=l14) {
				ij1 = ij-1; ji1 = ji - 1;
				if(ij-ji<0) {
				  t = b[ij1-1];
				  b[ij1-1]=b[ji1-1];
				  b[ji1-1] = t;
				  
				  t = b[ij1];
				  b[ij1]=b[ji1];
				  b[ji1] = t;
				}
				ij += 2;
			      }
}

/*
 * int fastlog(n)
 *
 */

static int fastlog2(int n)
{
  int num_bits, power = 0;
  
  if((n < 2) || (n % 2 != 0)) return(0);
  num_bits = sizeof(int) * 8;   /* How big are ints on this machine? */
  
  while(power <= num_bits) {
    n >>= 1;
    power += 1;
    if(n & 0x01) {
      if(n > 1)	return(0);
      else return(power);
    }
  }
  return(0);
}

///////////////Added by pratibha Dec 24 for TRAPS-based system/////////////
void HammingTrap (int len_trap, muT *hamming_win) {
  int i;
  float base_angle;
  
  /* init the window */
  base_angle = 2. * M_PI / (float)(len_trap - 1);
  for (i=0; i< len_trap; i++)
    hamming_win->data[i] = 0.54 - (0.46 * cos(i * base_angle));
  return;
}

void mx_create_bandclassifier_input(muT *mx_traps_in, muT *mx_traps_out,int Band,muT *hamming_win) {
  int i , j , len;
  muT mx_buffer;
  float mn, var, aa, pp;

//  len = (mx_traps_in->d1)+50+9 ;
  len = (mx_traps_in->d1)+25+4 ;
  mx_buffer.data = (float *) calloc(len, sizeof(float));
  if ( mx_buffer.data == NULL ) 
    fprintf(stderr,"Error in allocation\n");
  mx_buffer.d1 = len;
  mx_buffer.d2 = 1 ;
  
  
  //for (i=0 ; i<41; i++) 
  for (i=0 ; i<21; i++) 
    mx_buffer.data[i]=mx_traps_in->data[0*(mx_traps_in->d2)+Band];
  
  //for (i=0;  i<9; i++)
    //mx_buffer.data[i+41]=mx_traps_in->data[(8-i)*(mx_traps_in->d2)+Band];
  for (i=0;  i<4; i++)
    mx_buffer.data[i+21]=mx_traps_in->data[(3-i)*(mx_traps_in->d2)+Band];
  
  for (i=0; i<(mx_traps_in->d1); i++) 
    mx_buffer.data[25+i]=mx_traps_in->data[i*(mx_traps_in->d2)+Band];
//    mx_buffer.data[50+i]=mx_traps_in->data[i*(mx_traps_in->d2)+Band];
  
  //for (i=0; i<9; i++)
    //mx_buffer.data[((mx_traps_in->d1)+50)+i]=mx_traps_in->data[((mx_traps_in->d1) - 1 - i)*(mx_traps_in->d2) +Band];
  for (i=0; i<4; i++)
    mx_buffer.data[((mx_traps_in->d1)+25)+i]=mx_traps_in->data[((mx_traps_in->d1) - 1 - i)*(mx_traps_in->d2) +Band];
  
//    mx_write_htk_file("Flipped_CB0", &mx_buffer, 100000, 9 ,1);
  
  for (i=0; i<mx_traps_in->d1; i++) 
    for (j=0; j<mx_traps_out->d2; j++) 
      mx_traps_out->data[i*(mx_traps_out->d2)+j] = mx_buffer.data[i+j]; 
  
  //Mean variance normalization
  
  for (i=0; i<mx_traps_in->d1; i++) {
    mn=0.0;
    var=0.0;
    for (j=0; j<mx_traps_out->d2; j++) {
      aa=mx_traps_out->data[i*(mx_traps_out->d2)+j];
      var=var+aa*aa;  
      mn=mn+aa;
    } 
    mn=mn/(mx_traps_out->d2);
    var=var/(mx_traps_out->d2);
    var=var-(mn*mn); 
    
    
    if (var > 0.0)  {
      pp=sqrt(var);
      for (j=0; j<mx_traps_out->d2; j++) {
	aa=mx_traps_out->data[i*(mx_traps_out->d2)+j]-mn;
	aa=aa/pp;
	mx_traps_out->data[i*(mx_traps_out->d2)+j] = aa * hamming_win->data[j] ;
      }
    } 
    else
      for (j=0; j<mx_traps_out->d2; j++) {
      mx_traps_out->data[i*(mx_traps_out->d2)+j] = (mx_traps_out->data[i*(mx_traps_out->d2)+j]-mn) * hamming_win->data[j] ;
    } 
    
  }   
  
  mx_release(&mx_buffer);
}

///////////////////////////////////////////////////////////////

int mx_combine_traps(muT * mx_out_traps, muT *mx_out_comb_traps, int Nbands) {
  
  int i , j, band=0;
  int totaldim;
  int mlpdim;
  int totalframes;
  
  mlpdim=mx_out_traps[band].d2;
  totaldim=mlpdim*Nbands;
  totalframes=mx_out_traps[band].d1;
  
  if(mx_set_size((mx_out_traps[0].d1),totaldim,mx_out_comb_traps))
    return -1;
  
  // No need to take log if u are linearly forwarding band classifiers outputs
  
  for (i = 0 ; i < mx_out_traps[0].d1 ; i++ ) {
    for (j = 0 ; j < mx_out_traps[0].d2 ; j++ ) {
      for (band = 0 ; band < Nbands ; band++ ) 
	  mx_out_comb_traps->data[i*totaldim+band*mlpdim+j] = mx_out_traps[band].data[i*mlpdim+j];

      }
    }
  return 0;
}
 

// Transpose for taking the inverse of PCA or spectral lda
// For reconstructing bands from cepstra , inv of DCT or PCA is needed
// (for sym matrix inverse of a matrix is nothing but transpose of the matrix)
int mx_transpose(muT *dct,muT *idct) {
  int i, j ;
  if(mx_set_size(dct->d2, dct->d1, idct)) {
    return -1;
  }
  for ( i = 0 ; i < idct->d1 ; i++ ) {
    for ( j = 0 ; j < idct->d2 ; j++ ) {
      idct->data[i*idct->d2+j] = dct->data[i + dct->d2 * j];
    }
  } 
  return 0;
}

//

void mx_create_bandclassifier_input60(muT *mx_traps_in, muT *mx_traps_out,int Band,muT *hamming_win) {
  int i , j , len;
  muT mx_buffer;
  float mn, var, aa, pp;

  len = (mx_traps_in->d1)+50+9 ;
  mx_buffer.data = (float *) calloc(len, sizeof(float));
  if ( mx_buffer.data == NULL )
    fprintf(stderr,"Error in allocation\n");
  mx_buffer.d1 = len;
  mx_buffer.d2 = 1 ;
 
 
  for (i=0 ; i<41; i++)
    mx_buffer.data[i]=mx_traps_in->data[0*(mx_traps_in->d2)+Band];
 
  for (i=0;  i<9; i++)
    mx_buffer.data[i+41]=mx_traps_in->data[(8-i)*(mx_traps_in->d2)+Band];
 
  for (i=0; i<(mx_traps_in->d1); i++)
    mx_buffer.data[50+i]=mx_traps_in->data[i*(mx_traps_in->d2)+Band];
 
  for (i=0; i<9; i++)
    mx_buffer.data[((mx_traps_in->d1)+50)+i]=mx_traps_in->data[((mx_traps_in->d1) - 1 - i)*(mx_traps_in->d2) +Band];
 
//    mx_write_htk_file("Flipped_CB0", &mx_buffer, 100000, 9 ,1);
 
  for (i=0; i<mx_traps_in->d1; i++)
    for (j=0; j<mx_traps_out->d2; j++)
      mx_traps_out->data[i*(mx_traps_out->d2)+j] = mx_buffer.data[i+j];

  //Mean variance normalization

  for (i=0; i<mx_traps_in->d1; i++) {
    mn=0.0;
    var=0.0;
    for (j=0; j<mx_traps_out->d2; j++) {
      aa=mx_traps_out->data[i*(mx_traps_out->d2)+j];
      var=var+aa*aa;
      mn=mn+aa;
    }
    mn=mn/(mx_traps_out->d2);
    var=var/(mx_traps_out->d2);
    var=var-(mn*mn);


    if (var > 0.0)  {
      pp=sqrt(var);
      for (j=0; j<mx_traps_out->d2; j++) {
        aa=mx_traps_out->data[i*(mx_traps_out->d2)+j]-mn;
        aa=aa/pp;
        mx_traps_out->data[i*(mx_traps_out->d2)+j] = aa * hamming_win->data[j] ;
      }
    }
    else
      for (j=0; j<mx_traps_out->d2; j++) {
      mx_traps_out->data[i*(mx_traps_out->d2)+j] = (mx_traps_out->data[i*(mx_traps_out->d2)+j]-mn) * hamming_win->data[j] ;
    }

  }

  mx_release(&mx_buffer);
}

