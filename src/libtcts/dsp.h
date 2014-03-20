#ifndef  __DSP_H
#define  __DSP_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#ifdef __cplusplus
extern "C" {
#endif

#define __FLOAT32__ 

#ifndef NULL
  #define NULL 0
#endif  

#ifdef __FLOAT32__
#   define FLOAT float
#else
#   define FLOAT double
#endif 

typedef struct {
  FLOAT re;
  FLOAT im;
} COMPLEX;

#define round(x) ((long int)(2.0*(x))-(long int) (x))

void filter(FLOAT* BCoef, int NumOrder, FLOAT* ACoef, int DenOrder,
            FLOAT* x, int NPoints, FLOAT* y, FLOAT* ZI);
/* Computes the output y of a general FIR-IIR digital filter when
 signal x (length=NPoints) is presented at the input.
 The filter realization is the Transposed Direct Form II.
 Its transfer function is :

          B^[0] + B^[1] z^-1   + ... + B^[NumOrder] z^-NumOrder
    T(z)=-------------------------------------------------------
            1   + A^[1] z^-1   + ... + A^[DenOrder] z^-DenOrder

 ZI is the vector of internal variables of the filter.
 It should have a minimum dimension = max(NumOrder,DenOrder)
 It is set to the final conditions on return.
 WARNING : DenOrder and NumOrder shouldn't be both set to zero
        at the same time.
 MATLAB compatibility :  total with "y=filter(b,a,x,zi)"
*/


  void ifft_C (COMPLEX* ComplexFFT, int NPoints, COMPLEX* ComplexSignal, int NFFT);  


#ifdef __cplusplus
}
#endif

#endif
