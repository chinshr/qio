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
#include "libicsi/auroralib.h"

int fprintf_debug(char* message) {
  if (debugflag) {
    fprintf(stderr,message);
  }
  return 1;
}

void get_lda_filter(char *filterfile, muT *BCoef, muT *ACoef, int *latency, int nfilters) {

  FILE *filfil;  
  int count;
  float *BCoef0;
  float *ACoef0;
  float *BCoef1;
  float *ACoef1;
  int delay0;
  int delay1;
  int nb0, na0, nb1, na1, nb, na;

  filfil=fopen(filterfile,"r");
  fscanf(filfil,"%d %d",&nb0,&na0);
  BCoef0 = (float*)lmalloc((nb0+1)*sizeof(float));
  ACoef0 = (float*)lmalloc((na0+1)*sizeof(float));
  for (count=0;count<=nb0;count++) {
    fscanf(filfil,"%f",&BCoef0[count]);
  }
  for (count=0;count<=na0;count++) {
    fscanf(filfil,"%f",&ACoef0[count]);
  }
  fscanf(filfil,"%d",&delay0);

  fscanf(filfil,"%d %d",&nb1,&na1);
  BCoef1 = (float*)lmalloc((nb1+1)*sizeof(float));
  ACoef1 = (float*)lmalloc((na1+1)*sizeof(float));
  for (count=0;count<=nb1;count++) {
    fscanf(filfil,"%f",&BCoef1[count]);
  }
  for (count=0;count<=na1;count++) {
    fscanf(filfil,"%f",&ACoef1[count]);
  }
  fscanf(filfil,"%d",&delay1);
 
  nb = max(nb0,nb1);
  na = max(na0,na1);
  nb++;
  na++;
  mx_set_size(nfilters,nb,BCoef);
  mx_set_size(nfilters,na,ACoef);
  memset(BCoef->data,0,nfilters*nb*sizeof(float));
  memset(ACoef->data,0,nfilters*na*sizeof(float));
  for (count=0;count<nfilters;count++) {
    if (count<2) {
      memcpy(&(BCoef->data[count*nb]),BCoef0,(nb0+1)*sizeof(float));
      memcpy(&(ACoef->data[count*na]),ACoef0,(na0+1)*sizeof(float));
      latency[count]=delay0;
    }
    else {    
      memcpy(&(BCoef->data[count*nb]),BCoef1,(nb1+1)*sizeof(float));
      memcpy(&(ACoef->data[count*na]),ACoef1,(na1+1)*sizeof(float));
      latency[count]=delay1;
    }
  }
}

void lda_filter(char *filterfile, muT *in, muT *fout) {
  muT BCoef, ACoef;
  int *latency;
  latency = malloc(in->d2*sizeof(int));
  get_lda_filter(filterfile, &BCoef, &ACoef, latency, in->d2);
  //Changed to 28 - Sunil 01/04/2002
  mx_filter(&BCoef, &ACoef, latency, in, fout, 28);
}

void get_lp_filter(muT *BCoef, muT *ACoef, int *latency) {
  int nb, na;
  float BCoef0[6]={0.08895249559651,0.27489720439596,0.45087105994837,0.45087105994837,0.27489720439596,0.08895249559651};
  float ACoef0[6]={1.00000000000000,-0.20466646013718,0.75433478131523,-0.01516974960443,0.09118094189891,0.00376200640914};
  nb = 6;
  na = 6;
  mx_set_size(1,nb,BCoef);
  mx_set_size(1,na,ACoef);
  memcpy(BCoef->data,BCoef0,nb*sizeof(float));
  memcpy(ACoef->data,ACoef0,na*sizeof(float));
  *latency = 1;
}

void lp_filter(muT *in, muT *fout) {
  muT BCoef, ACoef;
  int latency;
  get_lp_filter(&BCoef, &ACoef, &latency);
  //Changed from 10 to 3 - SUnil 01/04/2002
  mx_filter(&BCoef, &ACoef, &latency, in, fout, 3);
}

void get_vadlp_filter(muT *BCoef, muT *ACoef, int *latency) {
  int nb, na;
  float BCoef0[6]={0.5};
  float ACoef0[6]={1.00000000000000,-0.5};
  nb = 1;
  na = 2;
  mx_set_size(1,nb,BCoef);
  mx_set_size(1,na,ACoef);
  memcpy(BCoef->data,BCoef0,nb*sizeof(float));
  memcpy(ACoef->data,ACoef0,na*sizeof(float));
  *latency = 0;
}

void vadlp_filter(muT *in, muT *fout) {
  muT BCoef, ACoef;
  int latency;
  get_vadlp_filter(&BCoef, &ACoef, &latency);
  //Changed from 10 to 0 - Sunil 01/04/2002
  mx_filter(&BCoef, &ACoef, &latency, in, fout, 4);
}

void mx_filter(muT *BCoef, muT *ACoef, int *latency, muT *in, muT *fout, int initl) {

  int r,c,i,k,l,delay;
  float *indata, *bufferdata;
  muT buffer;

  float mm;
  float *x;
  float *y;
  int trim, r0;
  int pos;

  int initialize_length;

  float *ZI;
  int nb, na, nz;
  
  initialize_length = initl;
  nb=BCoef->d2-1;
  na=ACoef->d2-1;
  nz=max(nb,na)+1;
  c=in->d2;
  ZI = (float*)lmalloc(nz*c*sizeof(float));
  memset(ZI,0,nz*c*sizeof(float));
  delay = 0;
  for(k=0;k<BCoef->d1;k++) {
    delay = max(delay, latency[k]);
  }
  r=in->d1+delay+initialize_length;
  mx_set_size(r,c,&buffer);
  mx_set_size(in->d1,in->d2,fout);

  indata = (float *) in->data;
  bufferdata = (float *) buffer.data;

  // fill the beginning //
  if (initialize_length) {
    for(k=0,l=c-1;k<c;k++,l--) {
      mm=0.0;
      for(i=0;i<10;i++) {
	mm+=indata[i*c+k];
      }
      mm/=10.0;
      for(i=0;i<initialize_length;i++) {
	//bufferdata[i*c+k]=mm;
	//Changed to repeating the first frame - Sunil 01/04/2002
	bufferdata[i*c+k]=indata[k];
      }
    }
  }

  // copy the data //
  memcpy(&(bufferdata[initialize_length*c]),indata,in->d1*in->d2*sizeof(float));
 
  // fill the end  //
  for(i=in->d1+initialize_length;i<r;i++) {
    for(k=0,l=c-1;k<c;k++,l--) {
      bufferdata[i*c+k]=indata[((in->d1)-1)*c+k];
    }
  }

  x = (float *)malloc(r*sizeof(float));
  y = (float *)malloc(r*sizeof(float));

  for(k=0;k<c;k++) {

    //put data
    for(i=0;i<r;i++) {
      x[i]=bufferdata[i*c+k];
    }

    if (BCoef->d1==1) {
      pos=0;
    }
    else {
      pos=k;
    }
    trim=latency[pos];
    r0 = in->d1+trim+initialize_length;
    filter(&(BCoef->data[pos*(nb+1)]), nb, &(ACoef->data[pos*(na+1)]), na, x, r0, y, &(ZI[k*nz]));
    
    //trim
    for(i=0;i<in->d1;i++) {
      fout->data[i*c+k]=y[i+trim+initialize_length];
    }    
  }

  free(x);
  free(y);
  free(ZI);

  mx_release(&buffer);
  return;

}

void lp_iir_filter(muT *in, muT *fout) {

  float BCoef0[6]={0.08895249559651,0.27489720439596,0.45087105994837,0.45087105994837,0.27489720439596,0.08895249559651};
  float ACoef0[6]={1.00000000000000,-0.20466646013718,0.75433478131523,-0.01516974960443,0.09118094189891,0.00376200640914};

  int initialize_length = 10;

  float mm;

  int r,c,i,k,l,delay;
  float *indata, *bufferdata;
  muT buffer;

  float *x;
  float *y;
  int trim, r0;

  float ZI0[6];

  int delay0 = 1;
  delay = 1;

  r=in->d1+delay+initialize_length;
  r0 = in->d1+delay0+initialize_length;
  c=in->d2;
  mx_set_size(r,c,&buffer);
  mx_set_size(in->d1,in->d2,fout);

  indata=(float *) in->data;
  bufferdata=(float *) buffer.data;

  // fill the beginning //
  if (initialize_length) {
    for(k=0,l=c-1;k<c;k++,l--) {
      mm=0.0;
      for(i=0;i<10;i++) {
		mm+=indata[i*c+k];
      }
      mm/=10.0;
      for(i=0;i<initialize_length;i++) {
		bufferdata[i*c+k]=mm;
      }
    }
  }

  // copy the data //
  memcpy(&(bufferdata[initialize_length*c]),indata,in->d1*in->d2*sizeof(float));
 
  // fill the end  //
  for(i=in->d1+initialize_length;i<r;i++) {
    for(k=0,l=c-1;k<c;k++,l--) {
      bufferdata[i*c+k]=indata[((in->d1)-1)*c+k];
    }
  }

  x = (float *)malloc(r*sizeof(float));
  y = (float *)malloc(r*sizeof(float));

  for(k=0;k<c;k++) {

    //put data
    for(i=0;i<r;i++) {
      x[i]=bufferdata[i*c+k];
    }

    //clean memory
    for(i=0;i<6;i++) {
      ZI0[i]=0;
    }

	trim=delay0;
	filter(BCoef0, 5, ACoef0, 5, x, r0, y, ZI0);
    
    //trim
    for(i=0;i<in->d1;i++) {
      fout->data[i*c+k]=y[i+trim+initialize_length];
    }
    
  }

  free(x);
  free(y);

  mx_release(&buffer);
  return;

}

void lda_iir_filter(char *filterfile, muT *in, muT *fout) {

  int r,c,i,k,l,delay;
  float *indata, *bufferdata;
  muT buffer;

  float *x;
  float *y;
  int trim, r0, r1;

  float mm;

  FILE *filfil;
  
  int count;

  int initialize_length;

  float *BCoef0;
  float *ACoef0;
  float *BCoef1;
  float *ACoef1;
  float *ZI0;
  float *ZI1;

  int delay0;
  int delay1;

  int nb0, na0, nb1, na1, nz;

  filfil=fopen(filterfile,"r");
  fscanf(filfil,"%d %d",&nb0,&na0);
  BCoef0 = (float*)lmalloc((nb0+1)*sizeof(float));
  ACoef0 = (float*)lmalloc((na0+1)*sizeof(float));
  for (count=0;count<=nb0;count++) {
    fscanf(filfil,"%f",&BCoef0[count]);
  }
  for (count=0;count<=na0;count++) {
    fscanf(filfil,"%f",&ACoef0[count]);
  }
  fscanf(filfil,"%d",&delay0);

  fscanf(filfil,"%d %d",&nb1,&na1);
  BCoef1 = (float*)lmalloc((nb1+1)*sizeof(float));
  ACoef1 = (float*)lmalloc((na1+1)*sizeof(float));
  for (count=0;count<=nb1;count++) {
    fscanf(filfil,"%f",&BCoef1[count]);
  }
  for (count=0;count<=na1;count++) {
    fscanf(filfil,"%f",&ACoef1[count]);
  }
  fscanf(filfil,"%d",&delay1);

  nz=max(nb0,nb1)+1;
  ZI0 = (float*)lmalloc(nz*sizeof(float));
  ZI1 = (float*)lmalloc(nz*sizeof(float));

  delay = max(delay0, delay1);


  initialize_length = 20;

  r=in->d1+delay+initialize_length;
  r0 = in->d1+delay0+initialize_length;
  r1 = in->d1+delay1+initialize_length;
  c=in->d2;
  mx_set_size(r,c,&buffer);
  mx_set_size(in->d1,in->d2,fout);

  indata=(float *) in->data;
  bufferdata=(float *) buffer.data;

  // fill the beginning //
  if (initialize_length) {
    for(k=0,l=c-1;k<c;k++,l--) {
      mm=0.0;
      for(i=0;i<10;i++) {
		mm+=indata[i*c+k];
      }
      mm/=10.0;
      for(i=0;i<initialize_length;i++) {
		bufferdata[i*c+k]=mm;
      }
    }
  }

  // copy the data //
  memcpy(&(bufferdata[initialize_length*c]),indata,in->d1*in->d2*sizeof(float));
 
  // fill the end  //
  for(i=in->d1+initialize_length;i<r;i++) {
    for(k=0,l=c-1;k<c;k++,l--) {
      bufferdata[i*c+k]=indata[((in->d1)-1)*c+k];
    }
  }

  x = (float *)malloc(r*sizeof(float));
  y = (float *)malloc(r*sizeof(float));

  for(k=0;k<c;k++) {

    //put data
    for(i=0;i<r;i++) {
      x[i]=bufferdata[i*c+k];
    }

    //clean memory
    for(i=0;i<nz;i++) {
      ZI0[i]=0; ZI1[i]=0;
    }

    if (k < 2) {
      trim=delay0;
      filter(BCoef0, nb0, ACoef0, na0, x, r0, y, ZI0);
    } else {    
      trim=delay1;
      filter(BCoef1, nb1, ACoef1, na1, x, r1, y, ZI1);
    }
    
    //trim
    for(i=0;i<in->d1;i++) {
      fout->data[i*c+k]=y[i+trim+initialize_length];
    }
    
  }

  free(x);
  free(y);
 
  free(BCoef0); free(BCoef1); free(ACoef0); free(ACoef1); free(ZI0); free(ZI1); fclose(filfil);

  mx_release(&buffer);
  return;

}

void DCOffsetFilter( float *CircBuff, long BSize, int *BPointer, long nSamples, float tap ) {
  
  long i;

  for ( i=0; i<nSamples; i++ )
    {
      /* y[n]=x[n]-x[n-1]+(1-tap)*y[n-1] */
      CircBuff[(*BPointer+1)%BSize]=CircBuff[(*BPointer+2)%BSize]-
	CircBuff[(*BPointer+1)%BSize]+(1.0-tap)*CircBuff[*BPointer];
      *BPointer=(*BPointer+1)%BSize;
    }
}

void mx_vad_combine(float *silfea, float *silfea_energy_vad, int length) {
  int i;
  for (i=0;i<length;i++) {
    silfea[i] = max (silfea[i],silfea_energy_vad[i]);
  }
}

float *mx_energy_vad(muT *in) {

  float *silfea;
  int i,j;
  float *e;
  float local_alpha;
  int init_frames;
  float mean_e;
  float snr_threshold;
  int n_speech;
  int speech;
  int latency;

  int totspeech;

  latency = 4;
  speech = 0;
  n_speech = 0;
  snr_threshold = 20;
  mean_e = 0;
  init_frames = 10;
  e = (float*)lmalloc((in->d1+latency)*sizeof(float));
  silfea = (float*)lmalloc(((in->d1+1)/2)*sizeof(float));

  /* compute log energy */
  for (i=0;i<in->d1;i++) {
    
    e[i] = 0;
    for (j=0;j<in->d2;j++) {
      e[i] += in->data[i * in->d2 + j];
    }
    e[i] = 16.0/log(2.0)*log((e[i]+64)/64);
  }
  for (i=0;i<latency;i++) {
    e[in->d1+i] = e[in->d1-1];
  }
  
  /* compute vad flag */  
  totspeech = 0;
  for (i=0;i<in->d1;i++) {

    if (i < init_frames) {
      local_alpha = 1 - 1/(i+1);
    }
    else {
      local_alpha = 0.98;
    }
    
    if ( ( (e[i+latency]-mean_e) < snr_threshold) || (i < init_frames) ) {
      mean_e += (1-local_alpha)*(e[i+latency]-mean_e);
      if (mean_e < 90) {
		mean_e = 90;
      }
    }

    if (  (e[i+latency]-mean_e) > snr_threshold ) {
      n_speech++;
      if (n_speech > 4) {
		speech = 1;
      }
    }
    else {
      n_speech = 0;
    }

    if (!(i%2)) {
      silfea[i/2] = (1-speech);
      totspeech = totspeech+speech;
    }
  }

  if (totspeech<10) {
	for (i=0;i<10;i++) {
	  if (!(i%2)) {
		silfea[i/2] = 0;
	  }
	}
  }

  return(silfea);
}

float *read_sil_fea(char* ipfilname, long feat_length) {

  int nread;
  float *silfea; float dum;
  int idum;
  FILE *ipfil;
  silfea = (float*)lmalloc((feat_length)*sizeof(float));
  ipfil = fopen(ipfilname,"r");

  nread = 0;
  while ( ( nread < feat_length ) && (fscanf(ipfil,"%d %d %f",&idum, &idum, &silfea[nread]) != EOF ) ) {
	fscanf(ipfil,"%d %d %f",&idum, &idum, &dum);
	nread++;
  }
  for (nread=nread; nread<feat_length; nread++) {
    silfea[nread] = 1.0;
  }
  fclose(ipfil);
  return (silfea);
}

//Added 12/12/01 Sunil
float *read_sil_fea_down_up(char* ipfilname, long feat_length, int down_up) {

  int nread;
  float *silfea; float dum;
  int idum;
  FILE *ipfil;
  silfea = (float*)lmalloc((feat_length)*sizeof(float));
  ipfil = fopen(ipfilname,"r");

  nread = 0;
  while ((nread < feat_length) && (fscanf(ipfil,"%d %d %f",&idum, &idum, &silfea[nread])!= EOF)) {
	if(down_up) {
	  fscanf(ipfil,"%d %d %f",&idum, &idum, &dum);
	} else {
	  nread++;
	  fscanf(ipfil,"%d %d %f",&idum, &idum, &silfea[nread]);
	}
	nread++;
  }
  for (nread=nread; nread<feat_length; nread++) {
    silfea[nread] = 1.0;
  }
  fclose(ipfil);
  return (silfea);
}

int online_normalization_with_vad(muT *in, float alpha, float bias,  muT *out, int latency, char *sil_flag)
{
  
  //16k dec 2001
  //float mean[] = {23.747812,-0.988758,-1.036022,-0.138392,0.780827,0.373286,0.172459,
  //			  0.196979,0.356768,0.051367,-0.056195,-0.141021,0.149337,0.115492,
  //			  0.141460,4.744705,4.323703};
  //float std_sq[] = {SQR(7.383385),SQR(2.002117),SQR(1.419817),SQR(1.021199),SQR(0.883927),
  //				SQR(0.650637),SQR(0.561664),SQR(0.545514),SQR(0.532843),SQR(0.488156),
  //				SQR(0.462210),SQR(0.425797),SQR(0.398226),SQR(0.347861),SQR(0.312360),
  //				SQR(1.448336),SQR(1.129791)};

  //8k amsterdam
  float mean[] = { 23.815484, 1.037240, -0.382422, -0.203596, -0.557713, -0.051042, 
		   -0.208684, 0.073762, -0.100447,  0.007481, -0.062511, -0.001553, 
		   -0.123591, -0.006837, -0.140246,
                   4.744705, 4.323703 }; 
  float std_sq[] = { SQR(5.957326), SQR(1.350153), SQR(0.992368), SQR(0.685526), 
		     SQR(0.834512), SQR(0.545422), SQR(0.478728), SQR(0.476498), 
		     SQR(0.422000), SQR(0.417962), SQR(0.351819), SQR(0.361830), 
		     SQR(0.323899), SQR(0.322991), SQR(0.287901),
                     SQR(1.448336),SQR(1.129791)}; 

  //8k + 2x 16k features - 01/10/02 Sunil
  /*float mean[] = {25.303995,-0.824313,-1.181598,-0.179810,0.315798,0.241731,
				  -0.034465, 0.245547,0.262070,0.136384,-0.135733,-0.113913,
				  0.024102,0.081258,0.018602 ,4.744705,4.323703};
  float std_sq[] = {SQR(7.449790),SQR(2.538877),SQR(1.912556),SQR(1.321070),
					SQR(1.117010),SQR(0.833661),SQR(0.782323),SQR(0.692861),
					SQR(0.654425),SQR(0.623721),SQR(0.554693),SQR(0.513114),
					SQR(0.479340),SQR(0.408876),SQR(0.385325),SQR(1.448336),SQR(1.129791) };*/

  long i, j;

  if(mx_set_size(in->d1, in->d2, out)) {
    return -1;
  }

  for(i = 0; i < in->d1; i++) {
    for(j = 0; j < in->d2; j++) {
      //Update mean and var only during speech
      if(!sil_flag[i]) {
		mean[j] += (in->data[i * in->d2 + j] - mean[j]) * alpha;
		std_sq[j] += (SQR(in->data[i * in->d2 + j] - mean[j]) - std_sq[j]) * alpha;
      }
      if (i>=latency) {
		out->data[(i-latency)*out->d2+j] = (in->data[(i-latency)*in->d2+j] - mean[j]) / ((float)sqrt((double)std_sq[j]) + bias);
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

int normalization_with_vad(muT *in, muT *out, char *sil_flag)
{
  int latency = 0;
  long i, j;
  float mean, variance;
  int count;

  if(mx_set_size(in->d1, in->d2, out)) {
    return -1;
  }

  for(j = 0; j < in->d2; j++) {

    mean = 0;
    variance = 0;
    count = 0;

    for(i = 0; i < in->d1; i++) {
      //Update mean and var only during speech
      if(!sil_flag[i]) {
		count++;
		mean += in->data[i * in->d2 + j];
		variance += in->data[i * in->d2 + j]*in->data[i * in->d2 + j];
      }
    }
    mean /= count;
    variance /= count;
    variance -= mean*mean;

    for(i = 0; i < in->d1; i++) {
      out->data[(i-latency)*out->d2+j] = (in->data[(i-latency)*in->d2+j] - mean) / sqrt(variance) ;
    }
  
  }
  return 0;
}

int mx_cat_vectors(muT *in1, muT *in2, muT *out)
{
  
  int count,count2,count3;
  int n_params, n_frames;

  if (in1->d1 != in2->d1) {
    fprintf(stderr,"matrices do not have the same number of frames");
    exit(1);
  }
  //n_frames = min(in1->d1,in2->d1);
  n_frames = in1->d1;

  n_params = in1->d2+in2->d2;
  if(mx_set_size(n_frames,n_params,out)) {
    return -1;
  }
  
  for (count=0;count<n_frames;count++) {
    count3=0;
    for (count2=0;count2<in1->d2;count2++) {
      out->data[count*n_params+count3]=in1->data[count*in1->d2+count2];
      count3++;
    }
    for (count2=0;count2<in2->d2;count2++) {
      out->data[count*n_params+count3]=in2->data[count*in2->d2+count2];
      count3++;
    }
  }
  
  return 0;

}

int mx_average_vectors(muT *in1, muT *in2, muT *out)
{
  
  int count,count2;
  int n_params, n_frames;

  if (in1->d1 != in2->d1) {
    fprintf(stderr,"matrices do not have the same number of frames");
    exit(1);
  }
  if (in1->d2 != in2->d2) {
    fprintf(stderr,"matrices do not have the same number of elements");
    exit(1);
  }
  n_frames = in1->d1;
  n_params = in1->d2;

  if(mx_set_size(n_frames,n_params,out)) {
    return -1;
  }
  
  for (count=0;count<n_frames;count++) {
    for (count2=0;count2<n_params;count2++) {
      out->data[count*n_params+count2]=
		0.5*(in1->data[count*n_params+count2]+in2->data[count*n_params+count2]);
    }
  }

  return 0;
}

int mx_sub_vectors(muT *in, muT *out, int N)
{
  
  int count,count2,count3,count4;
  int n_params, n_frames;

  n_frames = in->d1;
  n_params = 3*N;

  if(mx_set_size(n_frames,n_params,out)) {
    return -1;
  }
  
  for (count=0;count<n_frames;count++) {
    count3=0;
    count4=0;
    for (count2=0;count2<in->d2/3;count2++) {
      if (count2<N) {
		out->data[count*n_params+count3]=in->data[count*in->d2+count4];
		count3++;
      }
      count4++;
    }
    for (count2=0;count2<in->d2/3;count2++) {
      if (count2<N) {
		out->data[count*n_params+count3]=in->data[count*in->d2+count4];
		count3++;
      }
      count4++;
    }
    for (count2=0;count2<in->d2/3;count2++) {
      if (count2<N) {
		out->data[count*n_params+count3]=in->data[count*in->d2+count4];
		count3++;
      }
      count4++;
    }
  }

  return 0;
}

int mx_dynamic_htk(muT* in, muT* out, int deltawin) {
  
  int i,j,k;
  float norm;
  
  mx_set_size(in->d1,in->d2,out);
  
  norm = 0;
  for (k=1;k<=deltawin;k++) {
    norm += 2*k*k;
  }
  
  for(i=0; i < in->d1; i++) {
    
    for(j=0; j < in->d2; j++) {
      out->data[out->d2 * i + j] = 0;
      for (k=1;k<=deltawin;k++) {
		out->data[out->d2 * i + j] += k * in->data[out->d2 * min(i+k,in->d1-1) + j]
		  - k * in->data[out->d2 * max(i-k,0) + j];
      } 
      out->data[out->d2 * i + j] /= norm;
      
    }
    
  }

  return 0;
}

int mx_ddynamic(muT* in, muT* out, int deltawin) {
  
  int i,j,k;
  float norm;
  float a, b;
  int deltawinall;
  float sum;

  mx_set_size(in->d1,in->d2,out);
  
  norm = 1;
  
  sum = 0;
  for (k=1;k<=deltawin;k++) {
    sum += k*k;
  }
  deltawinall = 2*deltawin+1;
  b = - (float)deltawinall / ( 2*sum - (float)deltawinall*deltawin*deltawin );
  a = 1.0 - b*deltawin*deltawin;

  for(i=0; i < in->d1; i++) {
    
    for(j=0; j < in->d2; j++) {
      out->data[out->d2 * i + j] = 0;
      for (k=-deltawin;k<=deltawin;k++) {
		out->data[out->d2 * i + j] +=
		  ( a + b * k * k) * in->data[out->d2 * max(min(i+k,in->d1-1),0) + j];
      } 
      out->data[out->d2 * i + j] /= norm;
      
    }
    
  }

  return 0;
}

#define MIN_NOISE_FRAMES 15

int mx_noise_level(muT *in, muT *noise)
{
  int i,j;

  if(mx_set_size(in->d1, in->d2, noise)) {
    return -1;
  }

  for (i=0; i < MIN_NOISE_FRAMES; i++)
    for (j=0; j < in->d2; j++)
      noise->data[j] += in->data[i*in->d2+j];

  for (j=0; j < in->d2; j++)
    noise->data[j] /= MIN_NOISE_FRAMES;
  for (i=1; i < in->d1; i++)
    for (j=0; j < in->d2; j++)
      noise->data[i*in->d2+j] = noise->data[j];

  return 0;

}    

int mx_noise_level_energyfloor(muT *in, muT *noise, float threshold, float alpha)
{
  int i,j;
  float local_alpha;
  float noise_energy;
  float frame_energy;
  
  if(mx_set_size(in->d1, in->d2, noise)) {
    return -1;
  }

  for (i=0; i < in->d1; i++) {
    if ( i < MIN_NOISE_FRAMES ) {
      local_alpha = 1.0 - 1.0/((float)i+1.0);
    }
    else {
      frame_energy = 0;
      for (j=0; j < in->d2; j++) {
		frame_energy += in->data[i*in->d2+j];
      }
      if (frame_energy < noise_energy * threshold) {
		local_alpha = alpha;
      }
      else {
		local_alpha = 1.0;
      }
    }
    noise_energy = 0;
    for (j=0; j < in->d2; j++) {
      if (i>0) {
		noise->data[i*in->d2+j] = exp ( log(noise->data[(i-1)*in->d2+j]) + (1.0-local_alpha)*(log(1.0+in->data[i*in->d2+j])-log(noise->data[(i-1)*in->d2+j])) );
      }
      else {
		noise->data[i*in->d2+j] = (1.0-local_alpha)*(1.0+in->data[i*in->d2+j]);
      }
      noise_energy += noise->data[i*in->d2+j];
    }
  }
  
  return 0;

}    

/* Added by gelbart April 2 2002, for the nr tool.

   mx_noise_level_energyfloor_vad
   
   This version of mx_noise_level_energyfloor takes as an argument a
   vector of binary silence flags (1 for nonspeech frames, 0 for
   speech frames).  There must be as many silence flags as the number
   of frames in the utterance.  This function will calculate a single
   noise estimate over all the frames labeled nonspeech.

*/
int mx_noise_level_energyfloor_vad(muT *in, muT *noise, float threshold, float alpha, char *silflag)
{
  int i,j;
  float local_alpha;
  int silframe_count;
  
  if(mx_set_size(in->d1, in->d2, noise)) {
    return -1;
  }

  /* Count number of silence frames seen so far */
  silframe_count = 0; 

  /* For each frame */
  for (i=0; i < in->d1; i++) {

    /* If this was judged a silence frame */
    if (silflag[i]) {
      local_alpha = 1.0 - 1.0/((float)silframe_count+1.0);
    }
    else {
      local_alpha = 1.0;
    }

    for (j=0; j < in->d2; j++) {
      
      /* If this is the first silence frame */
      if (silflag[i] && (silframe_count==0)) {
	/* Initialize the noise estimate */
	noise->data[i*in->d2+j] = 1.0+in->data[i*in->d2+j];
      }
      
      /* Else, as long as we have seen at least one silence frame so far
	 to initialize the noise estimate */
      else if (silframe_count>0) {
	 noise->data[i*in->d2+j] = exp ( log(noise->data[(i-1)*in->d2+j]) + (1.0-local_alpha)*(log(1.0+in->data[i*in->d2+j])-log(noise->data[(i-1)*in->d2+j])) );
      }
      
    }
    
    if (silflag[i]) {
      silframe_count++;
    }
  }
  
  if (silframe_count==0) {
    fprintf(stderr, "WARNING: No nonspeech frames in input file ; cannot get noise spectral estimate \r\n");  // Noise spectral estimate defaults to 0.0 since mx_set_size uses calloc, which sets all bits to zero (same as setting to 0.0 for IEEE floating point)
  }

  /* Now just set all noise estimates to the last noise estimate.  (In
     mx_noise_level_energyfloor_vad, we don't care about causality of the
     noise estimation.) */
  for (i=0; i < in->d1-1; i++) {
    for (j=0; j < in->d2; j++) {
      noise->data[i*in->d2+j] = noise->data[(in->d1-1)*in->d2+j];
    }
  }
   
  return 0;
}    


int mx_noise_level2(muT *in, muT *noise)
{
  int i,j;
  float global[MIN_NOISE_FRAMES];
  float mean, variance;
  int actual_noise_frames;

  if(mx_set_size(in->d1, in->d2, noise)) {
    return -1;
  }
  
  variance = 10;
  actual_noise_frames = MIN_NOISE_FRAMES;
  while (variance>3) {
    mean = 0;
    variance = 0;
    for (i=0; i < actual_noise_frames; i++) {
      global[i] = 0;
      for (j=0; j < in->d2; j++) {
		global[i] += in->data[i*in->d2+j];
      }
      mean += 10*log10(global[i]);
      variance += 10*log10(global[i])*10*log10(global[i]);
    }
    mean /= actual_noise_frames;
    variance /= actual_noise_frames;
    variance -= mean*mean;
    variance = sqrt(variance);
    actual_noise_frames--;
  }

  actual_noise_frames++;
  printf ("frames: %d\n",actual_noise_frames);
  for (i=0; i < actual_noise_frames; i++)
    for (j=0; j < in->d2; j++)
      noise->data[j] += in->data[i*in->d2+j];

  for (j=0; j < in->d2; j++)
    noise->data[j] /= actual_noise_frames;

  for (i=1; i < in->d1; i++) 
    for (j=0; j < in->d2; j++)
      noise->data[i*in->d2+j] = noise->data[j];
  
  return 0;

}

void simple_resample(muT *in, float ratio, muT *out) {
  int i,j;
  int left,right;
  int outframes = (1.0/ratio)*in->d1;
  float pos = 0;
  mx_set_size(outframes,in->d2,out);
  for (i=0;i<outframes;i++) {
    left = floor(pos);
    right = min(left+1,in->d1-1);
    for (j=0;j<in->d2;j++) {
      out->data[i*in->d2+j] = in->data[left*in->d2+j] + (pos-left)*(in->data[right*in->d2+j]-in->data[left*in->d2+j]);
    }
    pos += ratio;
  }
}
