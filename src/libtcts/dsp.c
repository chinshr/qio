/*  This source file is a subset of Thierry Dutoit's dsplib.c library
    of Matlab-like signal processing routines for scalars, vectors, and
    matrices.  It is based on revision 3.02 of the library.
 
    Author:  Thierry DUTOIT
  
    Copyright (c) 1991,92,93,94 Faculte Polytechnique de Mons -  TCTS labs */

#include "libtcts/dsp.h"

/*--------------------------------------------------------------------------*/
int error(char* message) {
  fprintf(stderr,message);
  return 1;
  }

/*--------------------------------------------------------------------------*/
void filter(FLOAT* BCoef, int NumOrder, FLOAT* ACoef, int DenOrder,
            FLOAT* x, int NPoints, FLOAT* y, FLOAT* ZI)
{
  int i;
  int j;
  FLOAT* vi;
  FLOAT* ai;
  FLOAT* bi;
  if (NumOrder>DenOrder)
    for (i=0;i<NPoints;i++) {
      vi=ZI;
      ai=ACoef+1;
      bi=BCoef;
      *y=*(vi++)+*(bi++) * *x;
      for (j=1;j<NumOrder;j++) {
	*(vi-1)=*vi+*bi * *x;
	vi++;
	bi++;
      }
      *(vi-1)=*bi * *x;
      vi=ZI;
      for (j=0;j<DenOrder;j++) {
	*vi-=*ai * *y;
	vi++;
	ai++;
      }
      x++; y++;
    }
  else
    for (i=0;i<NPoints;i++) {
      vi=ZI;
      ai=ACoef+1;
      bi=BCoef;
      *y=*(vi++)+*(bi++) * *x;
      for (j=1;j<=NumOrder;j++) {
	*(vi-1)=*vi+*bi * *x;
	vi++;
	bi++;
      }
      for (j=NumOrder+1;j<DenOrder;j++) {
	*(vi-1)=*vi;
	vi++;
      }
      *(vi-1)=0;
      vi=ZI;
      for (j=0;j<DenOrder;j++) {
	*vi-=*ai * *y;
	vi++;
	ai++;
      }
      x++; y++;
    }
}

/*--------------------------------------------------------------------------*/
void conj_CV(COMPLEX* ComplexVector, int Npts, COMPLEX* Result)
{
   int i;
   for (i=0;i<Npts;i++) (Result++)->im=-(ComplexVector++)->im;
}

/*--------------------------------------------------------------------------*/
/* Internal Variables */
int CurrentNFFT = 0;   /* Current Value of NFFT */
int FFTType = 0; /* Last FFT computed : _Real or _Complex */
FLOAT* cs; /* pointers on cosine table */

/*--------------------------------------------------------------------------*/
#ifdef __TARGETPC__
void fft (COMPLEX* data, FLOAT* cossin, int N);
/*   ASM routine to compute a complex FFT with the Split-Radix method.
 */
#else
void fft (COMPLEX* data, FLOAT* cossin, int N)
  /*   ANSI C routine to compute a complex FFT with the Split-Radix method.
*/
{
  int M;
  int i0;
  int i;
  int j;
  int k;
  int pas;
  int id;
  int n2;
  int n4;
  COMPLEX r1;
  COMPLEX r2;
  COMPLEX* tmp0;
  COMPLEX* tmp1;
  COMPLEX* tmp2;
  COMPLEX* tmp3;
  FLOAT* tmp;
  FLOAT s3;
  FLOAT co;
  FLOAT so;
  FLOAT ct;
  FLOAT st;

  M=(int)(log(N)/log(2.0));
  n2=N;
  pas=1;
  for (k=1;k<M;k++)
    {
      n4=n2/4;
      /* papillon pour Wø */
      i0=0;
      tmp0=data;
      id=2*n2;
      do   {
	do   {
	  tmp1=tmp0+n4;
	  tmp2=tmp1+n4;
	  tmp3=tmp2+n4;
	  r1.re=tmp0->re-tmp2->re;
	  r1.im=tmp0->im-tmp2->im;
	  tmp0->re+=tmp2->re;
	  tmp0->im+=tmp2->im;
	  r2.re=tmp1->re-tmp3->re;
	  r2.im=tmp1->im-tmp3->im;
	  tmp1->re+=tmp3->re;
	  tmp1->im+=tmp3->im;
	  tmp2->re=r1.re+r2.im;
	  tmp2->im=r1.im-r2.re;
	  tmp3->re=r1.re-r2.im;
	  tmp3->im=r2.re+r1.im;
	  i0+=id;
	  tmp0+=id;
	}   while (i0 < N-1);
	i0=2*id-n2;
	tmp0=&(data[i0]);
	id*=4;
      }   while (i0 < N-1);
      /* papillon g‚n‚ral */
      i=-1;
      for (j=1;j<n4;j++)
	{
	  i+=pas;
	  tmp=&(cossin[4*i]);
	  co=*(tmp++);
	  so=*(tmp++);
	  ct=*(tmp++);
	  st=*(tmp++);
	  i0=j;
	  tmp0=&(data[j]);
	  id=2*n2;
	  do   {
            do   {
	      tmp1=tmp0+n4;
	      tmp2=tmp1+n4;
	      tmp3=tmp2+n4;
	      r1.re=tmp0->re-tmp2->re;
	      r1.im=tmp0->im-tmp2->im;
	      tmp0->re+=tmp2->re;
	      tmp0->im+=tmp2->im;
	      r2.re=tmp1->re-tmp3->re;
	      r2.im=tmp1->im-tmp3->im;
	      tmp1->re+=tmp3->re;
	      tmp1->im+=tmp3->im;
	      s3=r1.re-r2.im;
	      r1.re+=r2.im;
	      r2.im=r1.im-r2.re;
	      r2.re+=r1.im;
	      tmp2->re=r1.re*co+r2.im*so;
	      tmp2->im=r2.im*co-r1.re*so;
	      tmp3->re=s3*ct+r2.re*st;
	      tmp3->im=r2.re*ct-s3*st;
	      i0+=id;
	      tmp0+=id;
	    }   while (i0 < N-1);
            i0=2*id-n2+j;
            tmp0=&(data[i0]);
            id*=4;
	  }   while (i0 < N-1);
	}
      n2/=2;
      pas*=2;
    }
  i0=0;
  tmp0=data;
  id=4;
  do   {
    do   {
      tmp1=tmp0+1;
      r1=*tmp0;
      tmp0->re+=tmp1->re;
      tmp1->re=r1.re-tmp1->re;
      tmp0->im+=tmp1->im;
      tmp1->im=r1.im-tmp1->im;
      i0+=id;
      tmp0+=id;
    }   while (i0<N-1);
    i0=2*id-2;
    tmp0=&(data[i0]);
    id*=4;
  }   while (i0<N-1);
  j=0;
  tmp0=data;
  tmp1=data;
  for (i=0;i<N-1;i++)
    {
      if (i<j)
	{
	  r1=*tmp0;
	  *tmp0=*tmp1;
	  *tmp1=r1;
	}
      k=N/2;
      while (k <= j) {j-=k; tmp1-=k; k/=2;}
      j+=k;
      tmp1+=k;
      tmp0++;
    }
}
#endif
/*--------------------------------------------------------------------------*/
void initialisation_CFFT(int N, FLOAT* cossin)
  /*   Initializes the trigonometric table used by the complex FFT routine above.
   N = FFT length (=max number of samples to be used)
*/
{
#define PI 3.14159265358979323846
  int i;
  FLOAT theta;
  FLOAT arg;
  theta=2*PI/N;
  arg=0;
  for (i=0;i<((int)(N/4)-1);i++) {
    arg+=theta;
    *(cossin++)=cos(arg);
    *(cossin++)=sin(arg);
    *(cossin++)=cos(3*arg);
    *(cossin++)=sin(3*arg);
  }
#undef PI
}
/*--------------------------------------------------------------------------*/
void ifft_C (COMPLEX* ComplexFFT, int NPoints, COMPLEX* ComplexSignal, int NFFT)
{
  int i;
  if (NFFT!=CurrentNFFT)   {
    if (CurrentNFFT != 0) free(cs);
    CurrentNFFT=NFFT;
    cs=(FLOAT*) calloc(NFFT+1,sizeof(FLOAT));
    if (cs==NULL) fprintf(stderr,"Error in fft_C");
    initialisation_CFFT(NFFT,cs);
    FFTType=sizeof(COMPLEX);
  }
  if (FFTType!=sizeof(COMPLEX))   {
    initialisation_CFFT(NFFT,cs);
    FFTType=sizeof(COMPLEX);
  }
  if (NPoints>NFFT) NPoints=NFFT;
  memcpy(ComplexSignal,ComplexFFT,NPoints*2*sizeof(FLOAT));
  conj_CV(ComplexSignal,NPoints,ComplexSignal);
  memset(&(ComplexSignal[NPoints]),0,(NFFT-NPoints)*2*sizeof(FLOAT));
  fft(ComplexSignal,cs,NFFT);
  for (i=0;i<NFFT;i++) {
    ComplexSignal->im/=-NFFT;
    ComplexSignal->re/=NFFT;
    ComplexSignal++;
  }
}
