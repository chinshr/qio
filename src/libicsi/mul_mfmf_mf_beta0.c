/*
**
** PHiPAC Matrix-Matrix Code for the operation:
**    C = A*B
**
** Automatically Generated by mm_gen ($Revision: 1.1 $) using the command:
**    ../mm_gen/mm_gen -cb 1 2 10 -cb 50 50 2 -beta 0 -file mul_mfmf_mf_beta0.c -rout mul_mfmf_mf_beta0 
**
** Run '../mm_gen/mm_gen -help' for help.
**
** Generated on: Tuesday December 12 1995, 18:11:56 PST
** Created by: Jeff Bilmes <bilmes@cs.berkeley.edu>
**             http://www.icsi.berkeley.edu/~bilmes/phipac
**
**
** Usage:
**    mul_mfmf_mf_beta0(const int M, const int K, const int N, const float *const A, const float *const B, float *const C, const int Astride, const int Bstride, const int Cstride)
** where
**  A is an MxK matrix
**  B is an KxN matrix
**  C is an MxN matrix
**  Astride is the number of entries between the start of each row of A
**  Bstride is the number of entries between the start of each row of B
**  Cstride is the number of entries between the start of each row of C
**
**
** "Copyright (c) 1995 The Regents of the University of California.  All
** rights reserved."  Permission to use, copy, modify, and distribute
** this software and its documentation for any purpose, without fee, and
** without written agreement is hereby granted, provided that the above
** copyright notice and the following two paragraphs appear in all copies
** of this software.
**
** IN NO EVENT SHALL THE UNIVERSITY OF CALIFORNIA BE LIABLE TO ANY PARTY FOR
** DIRECT, INDIRECT, SPECIAL, INCIDENTAL, OR CONSEQUENTIAL DAMAGES ARISING OUT
** OF THE USE OF THIS SOFTWARE AND ITS DOCUMENTATION, EVEN IF THE UNIVERSITY OF
** CALIFORNIA HAS BEEN ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
**
** THE UNIVERSITY OF CALIFORNIA SPECIFICALLY DISCLAIMS ANY WARRANTIES,
** INCLUDING, BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
** AND FITNESS FOR A PARTICULAR PURPOSE.  THE SOFTWARE PROVIDED HEREUNDER IS
** ON AN "AS IS" BASIS, AND THE UNIVERSITY OF CALIFORNIA HAS NO OBLIGATION TO
** PROVIDE MAINTENANCE, SUPPORT, UPDATES, ENHANCEMENTS, OR MODIFICATIONS.
**
*/

#define LOAD1x1(c00,C,Cstride) \
{\
   float * _cp = C; \
   c00 = _cp[0]; \
}

#define STORE1x1(c00,C,Cstride) \
{\
   float *_cp = C; \
   _cp[0] = c00; \
}

/* Fixed M,K,N = 1,1,1 fully-unrolled matrix matrix multiply. */
#define mul_mf1x1mf1x1_mf1x1(c00,A0,B,Bstride) \
{ \
   register float _b0; \
   register float _a0; \
   \
   \
   _b0 = B[0]; \
   B += Bstride; \
   _a0 = A0[0]; \
   c00 += _a0*_b0; \
}


/* Fixed M,K,N = 1,2,1 fully-unrolled matrix matrix multiply. */
#define mul_mf1x2mf2x1_mf1x1(c00,A0,B,Bstride) \
{ \
   register float _b0; \
   register float _a0; \
   \
   \
   _b0 = B[0]; \
   B += Bstride; \
   _a0 = A0[0]; \
   c00 += _a0*_b0; \
   \
   _b0 = B[0]; \
   B += Bstride; \
   _a0 = A0[1]; \
   c00 += _a0*_b0; \
}


#define LOAD1x2(c00,c01,C,Cstride) \
{\
   float * _cp = C; \
   c00 = _cp[0]; c01 = _cp[1]; \
}

#define STORE1x2(c00,c01,C,Cstride) \
{\
   float *_cp = C; \
   _cp[0] = c00; _cp[1] = c01; \
}

/* Fixed M,K,N = 1,1,2 fully-unrolled matrix matrix multiply. */
#define mul_mf1x1mf1x2_mf1x2(c00,c01,A0,B,Bstride) \
{ \
   register float _b0,_b1; \
   register float _a0; \
   \
   \
   _b0 = B[0]; _b1 = B[1]; \
   B += Bstride; \
   _a0 = A0[0]; \
   c00 += _a0*_b0; c01 += _a0*_b1; \
}


/* Fixed M,K,N = 1,2,2 fully-unrolled matrix matrix multiply. */
#define mul_mf1x2mf2x2_mf1x2(c00,c01,A0,B,Bstride) \
{ \
   register float _b0,_b1; \
   register float _a0; \
   \
   \
   _b0 = B[0]; _b1 = B[1]; \
   B += Bstride; \
   _a0 = A0[0]; \
   c00 += _a0*_b0; c01 += _a0*_b1; \
   \
   _b0 = B[0]; _b1 = B[1]; \
   B += Bstride; \
   _a0 = A0[1]; \
   c00 += _a0*_b0; c01 += _a0*_b1; \
}


#define LOAD1x4(c00,c01,c02,c03,C,Cstride) \
{\
   float * _cp = C; \
   c00 = _cp[0]; c01 = _cp[1]; c02 = _cp[2]; c03 = _cp[3]; \
}

#define STORE1x4(c00,c01,c02,c03,C,Cstride) \
{\
   float *_cp = C; \
   _cp[0] = c00; _cp[1] = c01; _cp[2] = c02; _cp[3] = c03; \
}

/* Fixed M,K,N = 1,1,4 fully-unrolled matrix matrix multiply. */
#define mul_mf1x1mf1x4_mf1x4(c00,c01,c02,c03,A0,B,Bstride) \
{ \
   register float _b0,_b1,_b2,_b3; \
   register float _a0; \
   \
   \
   _b0 = B[0]; _b1 = B[1]; _b2 = B[2]; _b3 = B[3]; \
   B += Bstride; \
   _a0 = A0[0]; \
   c00 += _a0*_b0; c01 += _a0*_b1; c02 += _a0*_b2; c03 += _a0*_b3; \
}


/* Fixed M,K,N = 1,2,4 fully-unrolled matrix matrix multiply. */
#define mul_mf1x2mf2x4_mf1x4(c00,c01,c02,c03,A0,B,Bstride) \
{ \
   register float _b0,_b1,_b2,_b3; \
   register float _a0; \
   \
   \
   _b0 = B[0]; _b1 = B[1]; _b2 = B[2]; _b3 = B[3]; \
   B += Bstride; \
   _a0 = A0[0]; \
   c00 += _a0*_b0; c01 += _a0*_b1; c02 += _a0*_b2; c03 += _a0*_b3; \
   \
   _b0 = B[0]; _b1 = B[1]; _b2 = B[2]; _b3 = B[3]; \
   B += Bstride; \
   _a0 = A0[1]; \
   c00 += _a0*_b0; c01 += _a0*_b1; c02 += _a0*_b2; c03 += _a0*_b3; \
}


#define LOAD1x8(c00,c01,c02,c03,c04,c05,c06,c07,C,Cstride) \
{\
   float * _cp = C; \
   c00 = _cp[0]; c01 = _cp[1]; c02 = _cp[2]; c03 = _cp[3]; c04 = _cp[4]; c05 = _cp[5]; c06 = _cp[6]; c07 = _cp[7]; \
}

#define STORE1x8(c00,c01,c02,c03,c04,c05,c06,c07,C,Cstride) \
{\
   float *_cp = C; \
   _cp[0] = c00; _cp[1] = c01; _cp[2] = c02; _cp[3] = c03; _cp[4] = c04; _cp[5] = c05; _cp[6] = c06; _cp[7] = c07; \
}

/* Fixed M,K,N = 1,1,8 fully-unrolled matrix matrix multiply. */
#define mul_mf1x1mf1x8_mf1x8(c00,c01,c02,c03,c04,c05,c06,c07,A0,B,Bstride) \
{ \
   register float _b0,_b1,_b2,_b3,_b4,_b5,_b6,_b7; \
   register float _a0; \
   \
   \
   _b0 = B[0]; _b1 = B[1]; _b2 = B[2]; _b3 = B[3]; _b4 = B[4]; _b5 = B[5]; _b6 = B[6]; _b7 = B[7]; \
   B += Bstride; \
   _a0 = A0[0]; \
   c00 += _a0*_b0; c01 += _a0*_b1; c02 += _a0*_b2; c03 += _a0*_b3; c04 += _a0*_b4; c05 += _a0*_b5; c06 += _a0*_b6; c07 += _a0*_b7; \
}


/* Fixed M,K,N = 1,2,8 fully-unrolled matrix matrix multiply. */
#define mul_mf1x2mf2x8_mf1x8(c00,c01,c02,c03,c04,c05,c06,c07,A0,B,Bstride) \
{ \
   register float _b0,_b1,_b2,_b3,_b4,_b5,_b6,_b7; \
   register float _a0; \
   \
   \
   _b0 = B[0]; _b1 = B[1]; _b2 = B[2]; _b3 = B[3]; _b4 = B[4]; _b5 = B[5]; _b6 = B[6]; _b7 = B[7]; \
   B += Bstride; \
   _a0 = A0[0]; \
   c00 += _a0*_b0; c01 += _a0*_b1; c02 += _a0*_b2; c03 += _a0*_b3; c04 += _a0*_b4; c05 += _a0*_b5; c06 += _a0*_b6; c07 += _a0*_b7; \
   \
   _b0 = B[0]; _b1 = B[1]; _b2 = B[2]; _b3 = B[3]; _b4 = B[4]; _b5 = B[5]; _b6 = B[6]; _b7 = B[7]; \
   B += Bstride; \
   _a0 = A0[1]; \
   c00 += _a0*_b0; c01 += _a0*_b1; c02 += _a0*_b2; c03 += _a0*_b3; c04 += _a0*_b4; c05 += _a0*_b5; c06 += _a0*_b6; c07 += _a0*_b7; \
}


#define LOAD1x10(c00,c01,c02,c03,c04,c05,c06,c07,c08,c09,C,Cstride) \
{\
   float * _cp = C; \
   c00 = _cp[0]; c01 = _cp[1]; c02 = _cp[2]; c03 = _cp[3]; c04 = _cp[4]; c05 = _cp[5]; c06 = _cp[6]; c07 = _cp[7]; c08 = _cp[8]; c09 = _cp[9]; \
}

#define STORE1x10(c00,c01,c02,c03,c04,c05,c06,c07,c08,c09,C,Cstride) \
{\
   float *_cp = C; \
   _cp[0] = c00; _cp[1] = c01; _cp[2] = c02; _cp[3] = c03; _cp[4] = c04; _cp[5] = c05; _cp[6] = c06; _cp[7] = c07; _cp[8] = c08; _cp[9] = c09; \
}

/* Fixed M,K,N = 1,1,10 fully-unrolled matrix matrix multiply. */
#define mul_mf1x1mf1x10_mf1x10(c00,c01,c02,c03,c04,c05,c06,c07,c08,c09,A0,B,Bstride) \
{ \
   register float _b0,_b1,_b2,_b3,_b4,_b5,_b6,_b7,_b8,_b9; \
   register float _a0; \
   \
   \
   _b0 = B[0]; _b1 = B[1]; _b2 = B[2]; _b3 = B[3]; _b4 = B[4]; _b5 = B[5]; _b6 = B[6]; _b7 = B[7]; _b8 = B[8]; _b9 = B[9]; \
   B += Bstride; \
   _a0 = A0[0]; \
   c00 += _a0*_b0; c01 += _a0*_b1; c02 += _a0*_b2; c03 += _a0*_b3; c04 += _a0*_b4; c05 += _a0*_b5; c06 += _a0*_b6; c07 += _a0*_b7; c08 += _a0*_b8; c09 += _a0*_b9; \
}


/* Fixed M,K,N = 1,2,10 fully-unrolled matrix matrix multiply. */
#define mul_mf1x2mf2x10_mf1x10(c00,c01,c02,c03,c04,c05,c06,c07,c08,c09,A0,B,Bstride) \
{ \
   register float _b0,_b1,_b2,_b3,_b4,_b5,_b6,_b7,_b8,_b9; \
   register float _a0; \
   \
   \
   _b0 = B[0]; _b1 = B[1]; _b2 = B[2]; _b3 = B[3]; _b4 = B[4]; _b5 = B[5]; _b6 = B[6]; _b7 = B[7]; _b8 = B[8]; _b9 = B[9]; \
   B += Bstride; \
   _a0 = A0[0]; \
   c00 += _a0*_b0; c01 += _a0*_b1; c02 += _a0*_b2; c03 += _a0*_b3; c04 += _a0*_b4; c05 += _a0*_b5; c06 += _a0*_b6; c07 += _a0*_b7; c08 += _a0*_b8; c09 += _a0*_b9; \
   \
   _b0 = B[0]; _b1 = B[1]; _b2 = B[2]; _b3 = B[3]; _b4 = B[4]; _b5 = B[5]; _b6 = B[6]; _b7 = B[7]; _b8 = B[8]; _b9 = B[9]; \
   B += Bstride; \
   _a0 = A0[1]; \
   c00 += _a0*_b0; c01 += _a0*_b1; c02 += _a0*_b2; c03 += _a0*_b3; c04 += _a0*_b4; c05 += _a0*_b5; c06 += _a0*_b6; c07 += _a0*_b7; c08 += _a0*_b8; c09 += _a0*_b9; \
}


/* Fixed M,N = 50,20, Arbitrary K L0-blocked matrix matrix multiply. */
static void
mul_mfmf_mf_beta0_l1_arb_k(int K, const float *const A, const float *const B, float *const C, const int Astride, const int Bstride, const int Cstride)
{
   const float *a0,*b0;
   float *c0;
   const float *ap0_0;
   const float *bp0;
   float *cp0;
   const int A_sbs_stride = Astride*1;
   const int C_sbs_stride = Cstride*1;
   const int k_marg_el = K & 1;
   const int k_norm = K - k_marg_el;
   float *const c0_endp = C+50*Cstride;
   register float c00,c01,c02,c03,c04,c05,c06,c07,c08,c09;
   for (c0=C,a0=A; c0!= c0_endp; c0+=C_sbs_stride,a0+=A_sbs_stride) {
      const float* const ap0_endp = a0 + k_norm;
      float* const cp0_endp = c0 + 20;
      for (b0=B,cp0=c0; cp0!=cp0_endp; b0+=10,cp0+=10) {
         ap0_0 = a0;
         bp0=b0;
         LOAD1x10(c00,c01,c02,c03,c04,c05,c06,c07,c08,c09,cp0,Cstride);
         for (; ap0_0!=ap0_endp; ap0_0+=2) {
            mul_mf1x2mf2x10_mf1x10(c00,c01,c02,c03,c04,c05,c06,c07,c08,c09,ap0_0,bp0,Bstride);
         }
         if (k_marg_el & 0x1) {
            mul_mf1x1mf1x10_mf1x10(c00,c01,c02,c03,c04,c05,c06,c07,c08,c09,ap0_0,bp0,Bstride);
         }
         STORE1x10(c00,c01,c02,c03,c04,c05,c06,c07,c08,c09,cp0,Cstride);
      }
   }
}

/* Arbitrary M,K,N L0-blocked matrix matrix multiply. */
static void
mul_mfmf_mf_beta0_l1_arb_all(const int M, const int K, const int N, const float *const A, const float *const B, float *const C, const int Astride, const int Bstride, const int Cstride)
{
   const float *a0,*b0;
   float *c0;
   const float *ap0_0;
   const float *bp0;
   float *cp0;
   const int A_sbs_stride = Astride*1;
   const int C_sbs_stride = Cstride*1;
   const int k_marg_el = K & 1;
   const int k_norm = K - k_marg_el;
   const int m_marg_el = M & 0;
   const int m_norm = M - m_marg_el;
   const int n_marg_el = N % 10;
   const int n_norm = N - n_marg_el;
   float *const c0_endp = C+m_norm*Cstride;
   register float c00,c01,c02,c03,c04,c05,c06,c07,c08,c09;
   for (c0=C,a0=A; c0!= c0_endp; c0+=C_sbs_stride,a0+=A_sbs_stride) {
      const float* const ap0_endp = a0 + k_norm;
      float* const cp0_endp = c0 + n_norm;
      for (b0=B,cp0=c0; cp0!=cp0_endp; b0+=10,cp0+=10) {
         ap0_0 = a0;
         bp0=b0;
         LOAD1x10(c00,c01,c02,c03,c04,c05,c06,c07,c08,c09,cp0,Cstride);
         for (; ap0_0!=ap0_endp; ap0_0+=2) {
            mul_mf1x2mf2x10_mf1x10(c00,c01,c02,c03,c04,c05,c06,c07,c08,c09,ap0_0,bp0,Bstride);
         }
         if (k_marg_el & 0x1) {
            mul_mf1x1mf1x10_mf1x10(c00,c01,c02,c03,c04,c05,c06,c07,c08,c09,ap0_0,bp0,Bstride);
         }
         STORE1x10(c00,c01,c02,c03,c04,c05,c06,c07,c08,c09,cp0,Cstride);
      }
   }
   for (c0=C,a0=A; c0!= c0_endp; c0+=C_sbs_stride,a0+=A_sbs_stride) {
      const float* const ap0_endp = a0 + k_norm;
      b0 = B+n_norm;
      cp0 = c0+n_norm;
      if (n_marg_el & 0x8) {
         ap0_0 = a0;
         bp0=b0;
         LOAD1x8(c00,c01,c02,c03,c04,c05,c06,c07,cp0,Cstride);
         for (; ap0_0!=ap0_endp; ap0_0+=2) {
            mul_mf1x2mf2x8_mf1x8(c00,c01,c02,c03,c04,c05,c06,c07,ap0_0,bp0,Bstride);
         }
         if (k_marg_el & 0x1) {
            mul_mf1x1mf1x8_mf1x8(c00,c01,c02,c03,c04,c05,c06,c07,ap0_0,bp0,Bstride);
         }
         STORE1x8(c00,c01,c02,c03,c04,c05,c06,c07,cp0,Cstride);
         b0 += 8;
         cp0 += 8;
      }
      if (n_marg_el & 0x4) {
         ap0_0 = a0;
         bp0=b0;
         LOAD1x4(c00,c01,c02,c03,cp0,Cstride);
         for (; ap0_0!=ap0_endp; ap0_0+=2) {
            mul_mf1x2mf2x4_mf1x4(c00,c01,c02,c03,ap0_0,bp0,Bstride);
         }
         if (k_marg_el & 0x1) {
            mul_mf1x1mf1x4_mf1x4(c00,c01,c02,c03,ap0_0,bp0,Bstride);
         }
         STORE1x4(c00,c01,c02,c03,cp0,Cstride);
         b0 += 4;
         cp0 += 4;
      }
      if (n_marg_el & 0x2) {
         ap0_0 = a0;
         bp0=b0;
         LOAD1x2(c00,c01,cp0,Cstride);
         for (; ap0_0!=ap0_endp; ap0_0+=2) {
            mul_mf1x2mf2x2_mf1x2(c00,c01,ap0_0,bp0,Bstride);
         }
         if (k_marg_el & 0x1) {
            mul_mf1x1mf1x2_mf1x2(c00,c01,ap0_0,bp0,Bstride);
         }
         STORE1x2(c00,c01,cp0,Cstride);
         b0 += 2;
         cp0 += 2;
      }
      if (n_marg_el & 0x1) {
         ap0_0 = a0;
         bp0=b0;
         LOAD1x1(c00,cp0,Cstride);
         for (; ap0_0!=ap0_endp; ap0_0+=2) {
            mul_mf1x2mf2x1_mf1x1(c00,ap0_0,bp0,Bstride);
         }
         if (k_marg_el & 0x1) {
            mul_mf1x1mf1x1_mf1x1(c00,ap0_0,bp0,Bstride);
         }
         STORE1x1(c00,cp0,Cstride);
      }
   }
}

/* Fixed M,K,N = 50,100,20 L0-blocked matrix matrix multiply. */
static void
mul_mfmf_mf_beta0_l1(const float *const A, const float *const B, float *const C, const int Astride, const int Bstride, const int Cstride)
{
   const float *a0,*b0;
   float *c0;
   const float *ap0_0;
   const float *bp0;
   float *cp0;
   const int A_sbs_stride = Astride*1;
   const int C_sbs_stride = Cstride*1;
   float *const c0_endp = C+50*Cstride;
   register float c00,c01,c02,c03,c04,c05,c06,c07,c08,c09;
   for (c0=C,a0=A; c0!= c0_endp; c0+=C_sbs_stride,a0+=A_sbs_stride) {
      const float* const ap0_endp = a0 + 100;
      float* const cp0_endp = c0 + 20;
      for (b0=B,cp0=c0; cp0!=cp0_endp; b0+=10,cp0+=10) {
         ap0_0 = a0;
         bp0=b0;
         LOAD1x10(c00,c01,c02,c03,c04,c05,c06,c07,c08,c09,cp0,Cstride);
         for (; ap0_0!=ap0_endp; ap0_0+=2) {
            mul_mf1x2mf2x10_mf1x10(c00,c01,c02,c03,c04,c05,c06,c07,c08,c09,ap0_0,bp0,Bstride);
         }
         STORE1x10(c00,c01,c02,c03,c04,c05,c06,c07,c08,c09,cp0,Cstride);
      }
   }
}

void
mul_mfmf_mf_beta0(const int M, const int K, const int N, const float *const A, const float *const B, float *const C, const int Astride, const int Bstride, const int Cstride)
{
   /* Code for L1-blocked routine. */
   int m2,k2,n2;
   const float *a2,*b2;
   float *c2;
   const float *ap2,*bp2;
   float *cp2;
   {
      float *cprb,*cpre,*cp,*cpe;
      cpre = C + M*Cstride;
      for (cprb = C; cprb != cpre; cprb += Cstride) {
         cpe = cprb + N;
         for (cp = cprb; cp != cpe; cp++) {
            *cp = 0.0;
         }
      }
   }
   if (M < 51 && K < 101 && N < 21) {
      mul_mfmf_mf_beta0_l1_arb_all(M,K,N,A,B,C,Astride,Bstride,Cstride);
      return;
   }
   for (m2=0; m2<=M-50; m2+=50) {
      c2 = C + m2*Cstride;
      a2 = A + m2*Astride;
      for (n2=0,b2=B,cp2=c2; n2<=N-20; n2+=20,b2+=20,cp2+=20) {
         for (k2=0,bp2=b2,ap2=a2; k2<=K-100; k2+=100,bp2+=100*Bstride,ap2+=100) {
            mul_mfmf_mf_beta0_l1(ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
         if (k2 < K) {
            mul_mfmf_mf_beta0_l1_arb_k(K-k2,ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
      }
      if (n2 < N) {
         for (k2=0,bp2=b2,ap2=a2; k2<=K-100; k2+=100,bp2+=100*Bstride,ap2+=100) {
            mul_mfmf_mf_beta0_l1_arb_all(50,100,N-n2,ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
         if (k2 < K) {
            mul_mfmf_mf_beta0_l1_arb_all(50,K-k2,N-n2,ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
      }
   }
   if (m2 < M) {
      c2 = C + m2*Cstride;
      a2 = A + m2*Astride;
      for (n2=0,b2=B,cp2=c2; n2<=N-20; n2+=20,b2+=20,cp2+=20) {
         for (k2=0,bp2=b2,ap2=a2; k2<=K-100; k2+=100,bp2+=100*Bstride,ap2+=100) {
            mul_mfmf_mf_beta0_l1_arb_all(M-m2,100,20,ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
         if (k2 < K) {
            mul_mfmf_mf_beta0_l1_arb_all(M-m2,K-k2,20,ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
      }
      if (n2 < N) {
         for (k2=0,bp2=b2,ap2=a2; k2<=K-100; k2+=100,bp2+=100*Bstride,ap2+=100) {
            mul_mfmf_mf_beta0_l1_arb_all(M-m2,100,N-n2,ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
         if (k2 < K) {
            mul_mfmf_mf_beta0_l1_arb_all(M-m2,K-k2,N-n2,ap2,bp2,cp2,Astride,Bstride,Cstride);
         }
      }
   }
}