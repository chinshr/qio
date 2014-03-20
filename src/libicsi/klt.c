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
 
#include <stdarg.h>
#include "libicsi/klt.h"


static void  errorExit(char *format, ...)
{
  va_list args;
  va_start(args, format);
  vfprintf(stderr, format, args);
  fprintf(stderr, "\n");
  fflush(stderr);
  va_end(args);

  exit(1);
}


void readStats(FILE*f, size_t N, int ascii, double *cor, double *means, double *vecs, double *vals)
{
  size_t i,j;

  for (i=0;i<N;i++) {
    if (!ascii) {
      for (j=0;j<N;j++) 
	if (!fread(&cor[i*N+j],sizeof(cor[0]),1,f))
	  errorExit("input stat file eof error, cor(%d,%d)\n",i,j);
      if (!fread(&means[i],sizeof(means[0]),1,f))
	  errorExit("input stat file eof error, mean(%d)\n",i);
      for (j=0;j<N;j++) 
	if (!fread(&vecs[i*N+j],sizeof(vecs[0]),1,f))
	  errorExit("input stat file eof error, vecs(%d,%d)\n",i,j);
      if (!fread(&vals[i],sizeof(vals[0]),1,f))
	  errorExit("input stat file eof error, vals(%d)\n",i);
    } else {
      for (j=0;j<N;j++) 
	if (!fscanf(f,"%lf",&cor[i*N+j]))
	  errorExit("input stat file eof error, cor(%d,%d)\n",i,j);
      if (!fscanf(f,"%lf",&means[i]))
	  errorExit("input stat file eof error, mean(%d)\n",i);
      for (j=0;j<N;j++) 
	if (!fscanf(f,"%lf",&vecs[i*N+j]))
	  errorExit("input stat file eof error, vecs(%d,%d)\n",i,j);
      if (!fscanf(f,"%lf",&vals[i]))
	  errorExit("input stat file eof error, vals(%d)\n",i);
    }
  }
  // Could check for eof condition to be
  // true here.
}

void readStats2(FILE*f,
	  const size_t N,
	  const int ascii,
	  double *vecs,
	  double *vals) 
{
  size_t i,j;

  if (!ascii) {
    for (i=0;i<N;i++)
      if (!fread(&vals[i],sizeof(double),1,f))
	errorExit("input stat file eof error, vals(%d)\n",i);
    for (i=0;i<N;i++)
      for (j=0;j<N;j++) 
	if (!fread(&vecs[i*N+j],sizeof(double),1,f))
	  errorExit("input stat file eof error, vecs(%d,%d)\n",i,j);
  } else {
    for (i=0;i<N;i++)
      if (!fscanf(f,"%lf",&vals[i]))
	errorExit("input stat file eof error, vals(%d)\n",i);

    for (i=0;i<N;i++)
      for (j=0;j<N;j++) 
	if (!fscanf(f,"%lf",&vecs[i*N+j]))
	  errorExit("input stat file eof error, vecs(%d,%d)\n",i,j);

  }

  // Could check for eof condition to be
  // true here.
}

void mx_perframenorm(muT* in, muT* out)
{
  int i,j;
  float sum=0;

  mx_set_size(in->d1, in->d2, out);

  for(i=0; i < in->d1; i++) {
    sum = 0;
    for(j=0; j < in->d2; j++) {
      sum += in->data[i*in->d2+j];
    }
    sum /= in->d2;
    for(j=0; j < in->d2; j++) {
      out->data[i*in->d2+j] = in->data[i*in->d2+j] - sum;
    }
  }
}

void mx_correctspeechsilence(muT* in, muT* out, float *silfea, int silstate)
{
  int i,j;
  float sum=0;
  float psvad, pstandem;

  mx_set_size(in->d1, in->d2, out);

  for(i=0; i < in->d1; i++) {

    // Input features are already downsampled , silfea is not  so downsample silfea
    psvad = silfea[2*i];
    pstandem = exp(in->data[i*in->d2+silstate]);
    sum = 0;
    for(j=0; j < in->d2; j++) {
      sum += exp(in->data[i*in->d2+j]);
    }
    pstandem /= sum;

    psvad=min(0.9999,max(0.0001,psvad));
    pstandem=min(0.9999,max(0.0001,pstandem));

    //printf("%f %f  ",psvad,pstandem);

    for(j=0; j < in->d2; j++) {
      if (j==silstate) {
	out->data[i*in->d2+j] = in->data[i*in->d2+j] + 0.5*log(psvad/pstandem);
      }
      else {
	out->data[i*in->d2+j] = in->data[i*in->d2+j] + 0.5*log((1.0-psvad)/(1.0-pstandem));
      }
    }

  }

}

void mx_addcorrectspeechsilence(muT* in, muT* out, float *silfea)
{
  int i,j;
  float psvad;

  mx_set_size(in->d1, in->d2+1, out);

  for(i=0; i < in->d1; i++) {

    psvad = silfea[i];
    psvad=min(0.9999,max(0.0001,psvad));

    for(j=0; j < in->d2; j++) {
      out->data[i*out->d2+j] = in->data[i*in->d2+j];
    }
    out->data[i*out->d2+j] = log(psvad);

  }

}


void mx_klt_forward(muT* in, muT* out, char *in_st_filename, int unity_variance, int ascii, int contextw, int n_outputs, int ldaform)
{
  
  FILE *in_st_fp;

    size_t buf_size = in->d1;
    int n_frames = buf_size;
    const size_t n_ftrs = in->d2;
    const size_t n_cftrs= n_ftrs*contextw; // number of features when
	                                       // incorporating context
    float *ftr_buf;
    float *ftr_buf_p;
    float *cftr_buf;

    // mean vector E[X}
    double *ftr_means;
    double *ftr_means_p;
    double *ftr_means_endp;
    double *ftr_cov;
    double *ftr_eigenvecs;
    double *ftr_eigenvals;

    // Initialize the above declared arrays
    size_t i,j;

    // Other variables
    double *valp;
    double *vecp;
    double *vecpp;
    double valp_inv;

    float *oftr_buf;
    float *oftr_buf_p;
    double *ftr_dbuf_src;
    double *ftr_dbuf_src_p;
    double *ftr_dbuf_src_endp;
    double *ftr_dbuf_dst;
    double *ftr_dbuf_dst_p;   
    float *ftr_buf_base;
    float *cftr_buf_base;

    // Open stat file 
    in_st_fp = fopen(in_st_filename,"r");

    // Allocate mem
    mx_set_size(buf_size,n_outputs,out);
    cftr_buf = (float*)malloc(buf_size*n_cftrs*sizeof(float));
    ftr_means = (double*)malloc(n_cftrs*sizeof(double));
    ftr_means_endp = ftr_means+n_cftrs;
 
    memset(ftr_means,0,n_cftrs*sizeof(double));

    if (in_st_fp) {
      // read in the matrix containing the data.
      ftr_cov = (double*)malloc(n_cftrs*n_cftrs*sizeof(double));
      ftr_eigenvecs = (double*)malloc(n_cftrs*n_cftrs*sizeof(double));
      ftr_eigenvals = (double*)malloc(n_cftrs*sizeof(double));
      if (!ldaform) {
	readStats(in_st_fp,n_cftrs,ascii,
		   ftr_cov,ftr_means,
		   ftr_eigenvecs,ftr_eigenvals);
      }
      else {
	readStats2(in_st_fp,n_cftrs,ascii,
		   ftr_eigenvecs,ftr_eigenvals);
      }
    }
    else {
      fprintf (stderr,"No KLT statistics found");
    }

    // at this point we no longer need the following
    free(ftr_cov);

    if (unity_variance) {
      // multily in the eigenvalues into the eigenvectors
      valp = ftr_eigenvals;
      vecp = ftr_eigenvecs;
      for (i=0;i<n_cftrs;i++) {
	vecpp = vecp;
	if(*valp<SMALL) {
	  *valp=SMALL;
	}
	valp_inv = sqrt(1.0/(*valp));;
	for (j=0;j<n_cftrs;j++) {
	  (*vecpp) *= valp_inv;
	  vecpp+=n_cftrs;
	}
	valp++;
	vecp++;
      }
    }

    // allocate space performing the orthogonalization.
    //float *oftr_buf = new float[buf_size * n_outputs];
    ftr_dbuf_src = (double*)malloc(buf_size*n_cftrs*sizeof(double));
    ftr_dbuf_dst = (double*)malloc(buf_size*n_cftrs*sizeof(double));

    //ftr_buf should contain the data TODO
    ftr_buf = &(in->data[0]);

    //create contextualized buffer here
    ftr_buf_base = ftr_buf;
    cftr_buf_base = cftr_buf;
    for (i=0;i<buf_size-contextw+1;i++) {
      for (j=0;j<n_cftrs;j++) {
	*cftr_buf_base=*ftr_buf_base;
	ftr_buf_base++;
	cftr_buf_base++;
	
      }
      ftr_buf_base-=(contextw-1)*n_ftrs;
    }
    // done creating contextualized buffer
    
    // While subtracting off the means, convert the features to 
    // doubles and do the normalization there.
    ftr_dbuf_src_p = ftr_dbuf_src;
    ftr_dbuf_src_endp = ftr_dbuf_src + (n_frames-contextw+1)*n_cftrs;
    ftr_buf_p  = cftr_buf;
    while (ftr_dbuf_src_p != ftr_dbuf_src_endp) {
      ftr_means_p = ftr_means;
      while (ftr_means_p != ftr_means_endp)
	*ftr_dbuf_src_p++ = (double) *ftr_buf_p++ - *ftr_means_p++;
    }   
    
    // Normalize the features with matrix multiply.
    // mul_mdmd_md(const int M, const int K, const int N, 
    //             const double *const A, const double *const B, double *const C, 
    //             const int Astride, const int Bstride, const int Cstride)    
    mul_mdmd_md(n_frames-contextw+1,n_cftrs,n_cftrs,
		ftr_dbuf_src,ftr_eigenvecs,ftr_dbuf_dst,
		n_cftrs,n_cftrs,n_cftrs);    
    
    // should set out here TODO frame: n_frames-contextw+1 are computed
    oftr_buf = &(out->data[0]);

    ftr_dbuf_dst_p = ftr_dbuf_dst;
    oftr_buf_p = oftr_buf;
    for (i=0;i<buf_size-contextw+1;i++) {
      for (j=0;j<n_outputs;j++) {
	*oftr_buf_p++ = (float)ftr_dbuf_dst_p[j];
      }
      ftr_dbuf_dst_p += n_cftrs;
    }
        
    //free oftr_buf;
    free(ftr_dbuf_src);
    free(ftr_dbuf_dst);
    
    //free ftr_buf;
    free(cftr_buf);
    free(ftr_means);
    free(ftr_eigenvecs);
    free(ftr_eigenvals);

    fclose(in_st_fp);

}
