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
 
#ifndef _FE_UTIL_H
#define _FE_UTIL_H

/* The following defines x_int32 and x_int16, which are meant to be
   integers of size 32 bits and 16 bits.  If the definitions are not
   correct on your system, edit them accordingly.  We did not simply
   use the definitions from the C99 stdint.h because not all compilers
   currently have stdint.h.  We use the x_ prefix to avoid conflict
   with stdint.h on systems where it does exist. */
typedef int x_int32;
typedef short x_int16;

#define MAX(a,b)        ((a)>(b) ? (a) : (b))
#define MIN(a,b)        ((a)<(b) ? (a) : (b))
#ifdef  min
#undef  min
#endif
#define min(a, b)       ((a) < (b) ? (a) : (b))
#ifdef  max
#undef  max
#endif
#define max(a, b)       ((a) < (b) ? (b) : (a))      

typedef struct muT_ {
  int d1, d2;              /* dimensions: row and column sizes */
  float *data;             /* data used to distribute the memory */
} muT;

int mu_FAST(float *b,int n);

/* This functions copies length samples centered at center.
   Pads with 0 if this window goes outside the range of n_samples */
void asr_memcpy(float *framebuff, float *in, int n_samples, int center, int length);

void mx_write_file(muT *m, char *name);
int mx_read_vec(short *data, int len, muT *out);
int mx_set_size(int rows, int cols, muT *out);
void mx_release(muT *in);
void get_genhanning_window(float *window, int windowsz, float alpha);
void preemphasis(float *speech, int windowsz, float preemp);
int mx_specgram_phase(muT* in, int windowsz, int framesz, int npoints, float *window, muT* out, float preemp, muT *phase);
void mx_log(muT* in);
int mx_prod(muT *in1, muT *in2, muT *out);
int mx_prod_high_samp(muT *in1, int fs, muT *in2, muT *out);
int get_upper_band_energy(muT *mx_in, int fs, muT *mx_upper_band_energy);
int get_dcttable(int N, int K, muT *out);
int get_dcttable_s(int N, int K, muT *dct, muT *idct);
int get_sldtable(char *name, int nparams, muT *out);
int get_mel_filterbank(int nFilters, int nPoints, float Fs, float Fstart, float Fend, muT* out);
int temp_filter(muT *in, muT *filter, muT *fout);
int mx_pow2(muT* in, muT* inIm, muT* out);
int read_raw_file(char  *filename, short **data, long  *length, int   swap);
int mx_range(muT *in, int rfirst, int rstep, int rlast, int cfirst, int cstep, int clast, muT *out);
int mx_add(muT* in1, muT* in2, muT* out);

int mx_read_htk_file(char              *htk_filename,
                     muT               *mx,
                     x_int32           *sampPeriod,
                     x_int16           *parmKind,
                     int               swap);

int mx_write_htk_file(char              *htk_filename,
                      muT               *mx,
                      x_int32           sampPeriod,
                      x_int16           parmKind,
                      int               swap);

int mx_transpose(muT *in, muT *out);
int online_normalization(muT *in, float alpha, float bias, muT *out, int lat);
int upsample2(muT *in, muT *out);
int downsample2(muT *in, int out_dim, muT *out);
char *sil_prob_to_flag(long length, float *silfea, float vad_threshold, int vad_rank_range, int vad_rank_point);
//int cut_sil_parts(muT *in, float *silfea, muT *out);
int cut_sil_parts_orig_ds2_silflag(muT *in, char *sil_flag_ds2, muT *out);
int framedropping(muT *in, char *sil_flag_ds2, int down_up,muT *out);
void * lmalloc(int bytes);
     
#define Melf(x)       (2595.0 * log10(1.0 + (x)/700.0))
#define InvMelf(x)    ((pow(10.0, (x)/2595.0) - 1.0)*700)

#define FEAT_VEC_SIZE 15

#define HID_LAYER_SIZE 15
#define OUT_LAYER_SIZE  2
#define INP_WIDTH       2
#define INP_MEMORY      3

#define VAD_RANK_RANGE 11
#define VAD_RANK_POINT  2
//#define VAD_RANK_RANGE 21
//#define VAD_RANK_POINT  4
#define VAD_TRESHOLD  0.5

struct htk_header {
  x_int32 nSamples;
  x_int32 sampPeriod;
  x_int16 sampSize;
  x_int16 parmKind;
} header, header_out;

/*typedef float t_fea_vec[FEAT_VEC_SIZE];*/

#define HTK_FILE_SUCCESS        0
#define HTK_FILE_CANT_OPEN      1
#define HTK_FILE_IO_ERROR       2
#define HTK_FILE_INVALID_HEADER 3

typedef struct MLPParams_ {
  int ninp; // number of inputs
  int nhid; // number of hidden units
  int nout; // number of outputs
  int cw; // number of context windows
  int width; // feature vector length
  float bval; // bias value
  float *wih; // input-hidden matrix
  float *who; // hidden-output matrix
  float *bias; // bias and scale vectors for feature noralization
  float *scale; // ---
} MLPParams;

void init_mlp_forward(MLPParams *params, char *mlp_filename, char *feat_norm_filename);

int mlp_forward_linear(muT *feat, muT *out, MLPParams params, int feat_length, char *silflag);

int mlp_forward(muT *feat, muT *out, MLPParams params, int feat_length, char *silflag, int sil);

extern void mul_mfmft_mf_beta0 (const int M, const int K, const int N, 
				const float *const A, 
				const float *const B, 
				float *const C, 
				const int Astride, 
				const int Bstride, 
				const int Cstride);

extern void mul_mfmf_mf_beta0 (const int M, const int K, const int N, 
			       const float *const A, 
			       const float *const B, 
			       float *const C, 
			       const int Astride, 
			       const int Bstride, 
			       const int Cstride);

extern void mul_mftmf_mf_alpha (const int M, const int K, const int N, 
				const float *const A, 
				const float *const B, 
				float *const C, 
				const int Astride, 
				const int Bstride, 
				const int Cstride, 
				const float alpha);

void bunch_layer_fwd_prop(int insize, int outsize, float *w, float *in, float *out, int addstride);
void bunch_nn_forw_prop(int Ni, int Nh, int No, float *who, float *wih, float *ipats, float *hpats, float *opats, int lin);
int bunch_mlp_forward_linear(muT *feat, muT *out, MLPParams params, int feat_length, char *silflag);
int bunch_mlp_forward(muT *feat, muT *out, MLPParams params, int feat_length, char *silflag, int lin);
void mx_create_bandclassifier_input60(muT *mx_traps_in, muT *mx_traps_out,int Band,muT *hamming_win);
void HammingTrap (int len_trap, muT *hamming_win);
int mx_combine_traps(muT * mx_out_traps, muT *mx_out_comb_traps, int Nbands);
#endif
