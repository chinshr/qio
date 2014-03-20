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

#include <sys/stat.h>
#include <stdio.h>
#include <errno.h>
#include <stdlib.h>
#include <assert.h>
#include <string.h>
#include <sys/types.h>

// OGI includes //
#include "libogi/fe_util.h"
// ICSI includes //
#include "libicsi/vadfeatures.h"
#include "libicsi/auroralib.h"
#include "libicsi/noisecomp.h"
#include "libicsi/overlapadd.h"

#define ON_ALPHA 0.1
#define ON_BIAS  1

//Mx structures
muT mx_dct;
muT mx_idct;
muT mx_dct23by23;
muT mx_idct23by23;
muT mx_dct15by15;
muT mx_idct15by15;

muT mx_in;
muT mx_out;
muT mx_cleannoisyratio;
muT mx_noise;
muT mx_upper_band_energy;
muT mx_mfb;
muT mx_delta;
muT mx_deltadelta;
muT mx_phase;
muT mx_synth;
muT mx_silfea;
muT mx_silfea_out;
muT mx_klt;
muT mx_straight;
muT mx_mlpout;
muT mx_autocor_based;
muT mx_in_meta;
muT mx_samples;
muT mx_inmlp;
muT mx_invad;
muT mx_vaddct;

//TRAPs global variable
muT mx_trapsmlpin;
muT mx_traps_in1;
muT mx_traps_in2;
muT mx_trapsmlpout;
muT mx_in_band_traps ;
muT mx_out_band_traps[15] ;
muT mx_concate_comb_traps ;
muT mx_out_final_traps ;
muT mx_out_final_traps_klt ;
muT mx_traps_recons_cb;
muT mx_traps_baseline;
muT hamming_win ;

int silence_state;

int from_straight = 15;
int from_tandem = 15;

typedef struct CommandArgs_ {
  
  //Windowing & FFT
  int dc_offset;
  int window_length;
  int window_shift;
  int fftlen;
  
  //ON-LINE NOISE ESTIMATION
  float noisest_threshold;
  float noisest_alpha;
  
  //WIENER FILTER
  int sps;
  char *sps_type;
  int sps_type_int;
  float sps_overest_low;
  float sps_overest_up;
  float sps_h_thresh;
  int sps_h_thresh_flag;
  float sps_filter_power;
  float sps_time_alpha;
  int sps_time_delay;
  int sps_freq_length;
  int sps_freq_flag;
  float sps_h_scale;
  int sps_addnoise_flag;
  float sps_addnoise_low;
  float sps_addnoise_up;
  float sps_addnoise_scale;
  int sps_hirschnoise;
  int sps_synth;
  
  //MEL FILTERBANK
  int start_freq;
  int end_freq;
  int n_filters;
  int mel;
  
  //UPPER-BANDS FEATURES
  int upper_flag;
  
  //VOICE ACTIVITY DETECTION
  char *vad_file;
  int vad;
  float vadthresh_fd;
  float vadthresh_oln;
  char *vad_weights;
  char *vad_norm;
  int fd;
  int vrr;
  int vrp;
  int vad_shift;
  
  //FEATURES DOWNSAMPLING/UPSAMPLING
  int down_up;
  int down_filter;
  
  //TEMPORAL LDA, SPECTRAL LDA and DCT
  int serverside_lda;
  char *lda_filters_file;
  char *spec_basis_file;
  int n_params;
  int dct;
  int featsize;
  
  //MEAN and VARIANCE NORMALIZATION
  float alpha;
  int onl_latency;
  int onl_init;
  char *onl_hist;
  
  //DYNAMIC FEATURES
  int delta_order;
  
  //NEURAL NET FEATURES
  int prop2;
  int traps;
  
  //META FEATURES
  int use_meta_features;
  int use_meta_features_fortandem;
  
  //IN/OUT
  char *input_file;
  char *output_file;
  int fs;
  int swapin;
  int swapout;
  
  //Misc
  int threshold;
  char *noise_input_file;
  
} CommandArgs;

// global variables to be used in server side
float *silfea;
char *silflag;
char *silflag_oln;
void server_tandem(CommandArgs command_args);
void server_straight(CommandArgs command_args);
void server_traps(CommandArgs command_args);

void help_command (CommandArgs *command_args, char *argv) {
  fprintf (stderr, "\r\nUSAGE:");
  fprintf (stderr, "   %s [options] -i input -o output\r\n", argv);
  fprintf (stderr, "\r\nwhere:\r\n");
  fprintf (stderr, "     -h                Print this help                                              \r\n\n");
  
  
  fprintf (stderr, "WINDOWING & FFT\n");
  fprintf (stderr, "     -O 0/1            DC offset filter (%d)                             \r\n", command_args->dc_offset);
  fprintf (stderr, "     -Length <int>     Window length in ms (%d)                             \r\n", command_args->window_length);
  fprintf (stderr, "     -Shift  <int>     Window shift in ms (%d)                             \r\n", command_args->window_shift);
  fprintf (stderr, "     -FFTLen  <int>    FFT window length (%d)                             \r\n\n", command_args->fftlen);
  
  
  fprintf (stderr, "ON-LINE NOISE ESTIMATION\n");
  fprintf (stderr, "     -NEthresh float   Threshold for noise estimation (%f)  \r\n", command_args->noisest_threshold);
  fprintf (stderr, "     -NEalpha float    Alpha for noise estimation (%f)  \r\n\n", command_args->noisest_alpha);
  
  
  fprintf (stderr, "WIENER FILTER\n");
  fprintf (stderr, "     -S 0/1            Noise compensation flag (%d)                             \r\n", command_args->sps);
  fprintf (stderr, "     -Stype string     Noise compensation type: QIO_FFT, QIO_MEL... (%s)                            \r\n", command_args->sps_type);   
  fprintf (stderr, "     -Sover_l float    Noise compensation lower overestimation factor (%f)                            \r\n", command_args->sps_overest_low);
  fprintf (stderr, "     -Sover_u float    Noise compensation upper overestimation factor (%f)                            \r\n", command_args->sps_overest_up);
  fprintf (stderr, "     -Shthresh1 float  Noise compensation threshold on transfert function (%f)                            \r\n", command_args->sps_h_thresh);
  fprintf (stderr, "     -Shthresh1_f 0/1  Noise compensation flag for bin-dependent threshold (%d)                            \r\n", command_args->sps_h_thresh_flag);
  fprintf (stderr, "     -Sfiltpow float   Noise compensation filter power (%f)                            \r\n", command_args->sps_filter_power);
  
  fprintf (stderr, "     -Stime_a float    Noise compensation time smoothing alpha (%f)                            \r\n", command_args->sps_time_alpha);
  fprintf (stderr, "     -Stime_d int      Noise compensation time smoothing latency (%d)                            \r\n", command_args->sps_time_delay);
  fprintf (stderr, "     -Sfreq_l int      Noise compensation frequency smoothing window (%d)                            \r\n", command_args->sps_freq_length);
  
  fprintf (stderr, "     -Sfreq_f int      Noise compensation flag for geometric mean frequency smoothing (%d)                            \r\n", command_args->sps_freq_flag);
  
  fprintf (stderr, "     -Shthresh2 float  Noise compensation threshold (scale factor on noise) (%f)                            \r\n", command_args->sps_h_scale);
  
  fprintf (stderr, "     -Snoise_f 0/1     Noise compensation flag for constant addition (%d)                            \r\n", command_args->sps_addnoise_flag);
  fprintf (stderr, "     -Snoise_l float   Noise compensation constant addition lower threshold (%f)                            \r\n", command_args->sps_addnoise_low);
  fprintf (stderr, "     -Snoise_u float   Noise compensation constant addition upper threshold (%f)                            \r\n", command_args->sps_addnoise_up);
  fprintf (stderr, "     -Snoise_s float   Noise compensation constant addition scale factor (%f)                            \r\n", command_args->sps_addnoise_scale);
  
  fprintf (stderr, "     -Shirschnoise 0/1 Noise estimation using Hirsch method (%d)                            \r\n", command_args->sps_hirschnoise);
  
  fprintf (stderr, "     -Ssynth 0/1       Noise compensation resynthesizes signal after cleaning\n                       -> output sample file filename = outputfilename.clean.raw (%d)                            \r\n\n", command_args->sps_synth);
  
  
  fprintf (stderr, "MEL FILTERBANK\n");
  fprintf (stderr, "     -Mel 0/1          Mel Filter analysis (%d) \r\n", command_args->mel);
  fprintf (stderr, "     -M <int>          Number of mel-filters (%d) \r\n", command_args->n_filters);
  fprintf (stderr, "     -CL frequency     Low frequency cutoff for mel bands (%d)   \r\n", command_args->start_freq);
  fprintf (stderr, "     -CH frequency     High frequency cutoff for mel bands (%d)  \r\n\n", command_args->end_freq);
  
  fprintf (stderr, "UPPER-BANDS FEATURES\n");
  fprintf (stderr, "     -UPPERflag 0/1    Flag for computing upper-bands energies if sample rate > 8000 (%d) \r\n\n", command_args->upper_flag);
  
  fprintf (stderr, "VOICE ACTIVITY DETECTION\n");
  fprintf (stderr, "     -s vad            Silence probabilities filename (ascii file from feacat on a pfile) \r\n");
  fprintf (stderr, "     -VADweights file  File containing the VAD weights\n                       (%s)\r\n", command_args->vad_weights);
  fprintf (stderr, "     -VADnorm file     File containing the normalization coefficients for the VAD inputs\n                       (%s)\r\n", command_args->vad_norm);
  fprintf (stderr, "                       If '-s' is not passed and '-VADweights void -VADnorm void', all silence probabilities are 0.\r\n");
  fprintf (stderr, "     -fd 0/1           Frame drop flag (%d)       \r\n",command_args->fd);
  fprintf (stderr, "     -VADfd float      VAD threshold for frame dropping (%f)  \r\n", command_args->vadthresh_fd);
  fprintf (stderr, "     -VADoln float     VAD threshold for online normalization (%f)  \r\n", command_args->vadthresh_oln);
  fprintf (stderr, "                       If not provided, do not drop any frames.\r\n");
  fprintf (stderr, "     -V 0/1            Use VAD  (%d)\r\n", command_args->vad);  
  fprintf (stderr, "     -VRR <int>        VAD rank range    (%d)       \r\n",command_args->vrr);
  fprintf (stderr, "     -VRP <int>        VAD rank point    (%d)       \r\n",command_args->vrp);
  fprintf (stderr, "     -VADshift <int>   VAD asynchrony in frames (>0 to reduce latency) with respect to features (%d)       \r\n\n",command_args->vad_shift);
  
  fprintf (stderr, "FEATURES DOWNSAMPLING/UPSAMPLING\n");
  fprintf (stderr, "     -D 0/1            Use downsampling/upsampling (%d)              \r\n", command_args->down_up);
  fprintf (stderr, "     -Dfilter 0/1      Use downsampling filter - you may want to set this to 0 if the LDA filter already cuts a lot of the high mod. freqs. (%d) \r\n\n", command_args->down_filter);
  
  
  fprintf (stderr, "TEMPORAL LDA, SPECTRAL LDA and DCT\n");
  fprintf (stderr, "     -LDAserver 0/1    Temporal LDA done on the server-side (%d)\r\n",command_args->serverside_lda);  
  fprintf (stderr, "     -F filters        Filename of the lda filters coefficients (use the provided void filters if you don't want any filtering)\n                       (%s)\r\n", command_args->lda_filters_file);
  fprintf (stderr, "     -SF filters       Filename of the spectral basis vectors\n                       \r\n");
  fprintf (stderr, "     -DCT 0/1          Apply spectral basis (%d)             \r\n", command_args->dct);
  fprintf (stderr, "     -NF  number       Number of coefficients of the output vector (%d)             \r\n", command_args->featsize);
  fprintf (stderr, "     -N number         Number of DCT Components (%d)             \r\n\n", command_args->n_params);
  
  
  fprintf (stderr, "MEAN and VARIANCE NORMALIZATION\n");
  fprintf (stderr, "     -A alpha          Value of the alpha parameter for on-line normalisation, assuming a 50Hz rate (%f)  \r\n", command_args->alpha);
  fprintf (stderr, "     -L latency        Mean and variance estimates latency (%d)             \r\n", command_args->onl_latency);
  fprintf (stderr, "     -OnlInit 0/1      Initialisation of on-line mean (%d)                             \r\n", command_args->onl_init);
  fprintf (stderr, "     -OnlHist file     History filename for online mean and variance normalization (none)                     \r\n\n");
  
  
  fprintf (stderr, "DYNAMIC FEATURES\n");
  fprintf (stderr, "     -Delta 0/1/2      0:static, 1:static+delta, 2:static+delta+doubledelta (%d)             \r\n\n", command_args->delta_order);
  
  
  fprintf (stderr, "NEURAL NET FEATURES\n");
  fprintf (stderr, "     -prop2 0/1        Flag to compute additional TANDEM stream (%d)             \r\n", command_args->prop2);
  fprintf (stderr, "     -TRAPS 0/1        Flag to compute additional TRAPS stream (%d)             \r\n\n", command_args->traps);
  
  
  fprintf (stderr, "META FEATURES\n");
  fprintf (stderr, "     -meta 0/1         Flag to compute additional META stream (%d)             \r\n", command_args->use_meta_features);
  fprintf (stderr, "     -metatandem 0/1   Flag to compute additional META stream for tandem (%d)             \r\n\n", command_args->use_meta_features_fortandem);
  
  
  fprintf (stderr, "IN/OUT\n");
  fprintf (stderr, "     -i input          Input filename REQUIRED (raw 16 bits, 8000Hz format) \r\n");
  fprintf (stderr, "     -fs 8000/11000/16000 Input sampling frequency (%d) \r\n",command_args->fs);
  fprintf (stderr, "     -o output         Output filename REQUIRED (htk features USER type) \r\n");
  fprintf (stderr, "     -swapin 0/1       Byte swap input  (%d)       \r\n",command_args->swapin);
  fprintf (stderr, "     -swapout 0/1      Byte swap output (%d)       \r\n\n",command_args->swapout);
  
  fprintf (stderr, "THINGS WE SHOULD PROBABLY GET RID OF\n");
  fprintf (stderr, "     -n noise          Input noise filename (raw 16 bits, 8000Hz format) \r\n");
  fprintf (stderr, "     -T 0/1            Noise addition (%d)              \r\n\n", command_args->threshold);
}

enum{QIO_FFT,QIO_MEL};

//#define FFTLEN 256
//#define SAMP_FREQ 8000
char temp2[1024];

void init_command (CommandArgs *command_args) {
  char *temp;
  
  //Windowing & FFT
  command_args->window_length = 25;
  command_args->window_shift = 10;
  command_args->dc_offset = 1;
  command_args->fftlen = 256;
  
  //ON-LINE NOISE ESTIMATION
  command_args->noisest_threshold = 2.0;
  command_args->noisest_alpha = 0.99;
  
  //WIENER FILTER
  command_args->sps = 1;
  command_args->sps_type = strdup("QIO_FFT");
  command_args->sps_type_int = QIO_FFT;
  command_args->sps_overest_low = 3.125;
  command_args->sps_overest_up = 1.25;
  command_args->sps_h_thresh = 0.01;
  command_args->sps_h_thresh_flag = 0;
  command_args->sps_filter_power = 2.0;
  command_args->sps_time_alpha = 0.1;
  command_args->sps_time_delay = 2;
  command_args->sps_freq_length = 21;
  command_args->sps_freq_flag = 0;
  command_args->sps_h_scale = 0.001;
  command_args->sps_addnoise_flag = 0;
  command_args->sps_addnoise_low = 0.001;
  command_args->sps_addnoise_up = 0.1;
  command_args->sps_addnoise_scale = 0.002;
  command_args->sps_hirschnoise = 0;
  command_args->sps_synth = 0;
  
  //MEL FILTERBANK
  command_args->start_freq = 64;
  command_args->end_freq = 4000;
  command_args->n_filters = 23;
  command_args->mel = 1;
  
  //UPPER-BANDS FEATURES
  command_args->upper_flag = 0;
  
  //VOICE ACTIVITY DETECTION
  temp = getenv("AURORACALC");
  if (!temp) {
    printf("The AURORACALC variable must be set to the root directory where the software is installed!\n");
    exit(1);
  }
  strcpy(temp2,temp);
  command_args->vad_weights = strdup(strcat(temp2,"/parameters/vad/net.tim-fin-tic-it-spn-rand.54i+50h+2o.0-delay-wiener+dct+lpf.wts.head"));
  strcpy(temp2,temp);
  command_args->vad_norm = strdup(strcat(temp2,"/parameters/vad/tim-fin-tic-it-spn-rand.0-delay-wiener+dct+lpf.norms")) ;
  command_args->vad_file = 0;
  command_args->vadthresh_fd = 0.5;
  command_args->vadthresh_oln = 0.5;
  command_args->fd = 1;
  command_args->vrr = 11;
  command_args->vrp = 2;
  command_args->vad = 1;
  command_args->vad_shift = 1;
  
  //FEATURES DOWNSAMPLING/UPSAMPLING
  command_args->down_up = 1;
  command_args->down_filter = 0; // no downsampling filter because it is integrated with the LDA
  
  //TEMPORAL LDA, SPECTRAL LDA and DCT
  command_args->serverside_lda = 0;
  strcpy(temp2,temp);
  command_args->lda_filters_file = strdup(strcat(temp2,"/parameters/lda/30-causal-noise_clean+ds_filter.fil")); 
  strcpy(temp2,temp);
  command_args->spec_basis_file = 0;
  command_args->n_params = 15;
  command_args->dct = 1;
  command_args->featsize = 45;
  
  //MEAN and VARIANCE NORMALIZATION
  command_args->alpha = 0.02;
  command_args->onl_latency = 0;
  command_args->onl_hist = 0;
  command_args->onl_init = 0;
  
  //DYNAMIC FEATURES
  command_args->delta_order = 2;
  
  //NEURAL NET FEATURES
  command_args->prop2 = 0;
  command_args->traps = 0;
  
  //META FEATURES
  command_args->use_meta_features = 0;
  command_args->use_meta_features_fortandem = 0;
  
  //IN & OUT
  command_args->fs = 8000;
  command_args->swapin = 1;
  command_args->swapout = 1;
  
  //Misc
  command_args->threshold = 0;
  command_args->noise_input_file = 0;
}

int parse_command (int argc, char *argv[], CommandArgs *command_args)
{
  int count = 1;
  FILE *fp_out;
  FILE *fp_in;
  
  if (argc==0) {
    help_command(command_args, argv[0]);
    return 0;
  }
  
  while (argc) {
    if (strcmp (argv[count], "-h") == 0) {
      help_command(command_args, argv[0]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-UPPERflag") == 0) {
      --argc;
      ++count;
      command_args->upper_flag = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-fs") == 0) {
      --argc;
      ++count;
      command_args->fs = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-swapin") == 0) {
      --argc;
      ++count;
      command_args->swapin = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-swapout") == 0) {
      --argc;
      ++count;
      command_args->swapout = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-D") == 0) {
      command_args->down_up = 0;
      --argc;
      ++count;
      command_args->down_up = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-fd") == 0) {
      --argc;
      ++count;
      command_args->fd = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-FFTLen") == 0) {
      --argc;
      ++count;
      command_args->fftlen = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-VRR") == 0) {
      --argc;
      ++count;
      command_args->vrr = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-VRP") == 0) {
      --argc;
      ++count;
      command_args->vrp = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-VADshift") == 0) {
      --argc;
      ++count;
      command_args->vad_shift = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-Dfilter") == 0) {
      --argc;
      ++count;
      command_args->down_filter = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-T") == 0) {
      --argc;
      ++count;
      command_args->threshold = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-CL") == 0) {
      --argc;
      ++count;
      command_args->start_freq = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-CH") == 0) {
      --argc;
      ++count;
      command_args->end_freq = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-V") == 0) {
      --argc;
      ++count;
      command_args->vad = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-M") == 0) {
      --argc;
      ++count;
      command_args->n_filters = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-Mel") == 0) {
      --argc;
      ++count;
      command_args->mel = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-O") == 0) {
      --argc;
      ++count;
      command_args->dc_offset = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-Length") == 0) {
      --argc;
      ++count;
      command_args->window_length = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-Shift") == 0) {
      --argc;
      ++count;
      command_args->window_shift = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-S") == 0) {
      --argc;
      ++count;
      command_args->sps = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-Stype") == 0) {
      --argc;
      ++count;
      command_args->sps_type = strdup(argv[count]);
      --argc;
      ++count;
      if (strcmp(command_args->sps_type,"QIO_FFT")==0) {
	command_args->sps_type_int = QIO_FFT;
      }
      else if (strcmp(command_args->sps_type,"QIO_MEL")==0) {
	command_args->sps_type_int = QIO_MEL;
      }
      /*else if (strcmp(command_args->sps_type,"ERIC")==0) {
	command_args->sps_type_int = ERIC;
	fprintf (stderr, "Spectral subtraction type (%s) not implemented",command_args->sps_type);
	return 0;
        }*/
      else {
	fprintf (stderr, "Spectral subtraction type (%s) not known",command_args->sps_type);
	return 0;
      }
    }
    else if (strcmp (argv[count], "-Sover_l") == 0) {
      --argc;
      ++count;
      command_args->sps_overest_low = atof(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-Sover_u") == 0) {
      --argc;
      ++count;
      command_args->sps_overest_up = atof(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-Shthresh1") == 0) {
      --argc;
      ++count;
      command_args->sps_h_thresh = atof(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-Shthresh1_f") == 0) {
      --argc;
      ++count;
      command_args->sps_h_thresh_flag = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-Sfiltpow") == 0) {
      --argc;
      ++count;
      command_args->sps_filter_power = atof(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-Stime_a") == 0) {
      --argc;
      ++count;
      command_args->sps_time_alpha = atof(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-Stime_d") == 0) {
      --argc;
      ++count;
      command_args->sps_time_delay = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-Sfreq_l") == 0) {
      --argc;
      ++count;
      command_args->sps_freq_length = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-Sfreq_f") == 0) {
      --argc;
      ++count;
      command_args->sps_freq_flag = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-Shthresh2") == 0) {
      --argc;
      ++count;
      command_args->sps_h_scale = atof(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-Snoise_f") == 0) {
      --argc;
      ++count;
      command_args->sps_addnoise_flag = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-Snoise_l") == 0) {
      --argc;
      ++count;
      command_args->sps_addnoise_low = atof(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-Snoise_u") == 0) {
      --argc;
      ++count;
      command_args->sps_addnoise_up = atof(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-Snoise_s") == 0) {
      --argc;
      ++count;
      command_args->sps_addnoise_scale = atof(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-Shirschnoise") == 0) {
      --argc;
      ++count;
      command_args->sps_hirschnoise = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-Ssynth") == 0) {
      --argc;
      ++count;
      command_args->sps_synth = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-A") == 0) {
      --argc;
      ++count;
      command_args->alpha = atof(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-N") == 0) {
      --argc;
      ++count;
      command_args->n_params = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-NF") == 0) {
      --argc;
      ++count;
      command_args->featsize = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-DCT") == 0) {
      --argc;
      ++count;
      command_args->dct = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-Delta") == 0) {
      --argc;
      ++count;
      command_args->delta_order = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-prop2") == 0) {
      --argc;
      ++count;
      command_args->prop2 = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-meta") == 0) {
      --argc;
      ++count;
      command_args->use_meta_features = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-metatandem") == 0) {
      --argc;
      ++count;
      command_args->use_meta_features_fortandem = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-L") == 0) {
      --argc;
      ++count;
      command_args->onl_latency = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-OnlInit") == 0) {
      --argc;
      ++count;
      command_args->onl_init = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-OnlHist") == 0) {
      --argc;
      ++count;
      command_args->onl_hist = strdup(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-VADfd") == 0) {
      --argc;
      ++count;
      command_args->vadthresh_fd = atof(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-VADoln") == 0) {
      --argc;
      ++count;
      command_args->vadthresh_oln = atof(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-VADweights") == 0) {
      --argc;
      ++count;
      command_args->vad_weights=strdup(argv[count]);
      if ( !strcmp(command_args->vad_weights,"void") ) {
	command_args->vad_weights = 0;
      }
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-VADnorm") == 0) {
      --argc;
      ++count;
      command_args->vad_norm=strdup(argv[count]);
      if ( !strcmp(command_args->vad_norm,"void") ) {
	command_args->vad_norm = 0;
      }
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-NEthresh") == 0) {
      --argc;
      ++count;
      command_args->noisest_threshold = atof(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-NEalpha") == 0) {
      --argc;
      ++count;
      command_args->noisest_alpha = atof(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-F") == 0) {
      --argc;
      ++count;
      command_args->lda_filters_file=strdup(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-LDAserver") == 0) {
      command_args->serverside_lda = 0;
      --argc;
      ++count;
      command_args->serverside_lda = atoi(argv[count]);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-SF") == 0) {
      --argc;
      ++count;
      command_args->spec_basis_file=strdup(argv[count]);
      if ( !strcmp(command_args->spec_basis_file,"void") ) {
	command_args->spec_basis_file = 0;
      }
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-i") == 0) {
      --argc;
      ++count;
      command_args->input_file=strdup(argv[count]);
      fp_in = fopen (command_args->input_file, "rb");
      if (fp_in == NULL) {
	fprintf (stderr, "ERROR:   Can not open file '%s' !\r\n",
		 command_args->input_file);
	return 0;
      }
      fclose(fp_in);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-n") == 0) {
      --argc;
      ++count;
      command_args->noise_input_file=strdup(argv[count]);
      fp_in = fopen (command_args->noise_input_file, "rb");
      if (fp_in == NULL) {
	fprintf (stderr, "ERROR:   Can not open file '%s' !\r\n",
		 command_args->noise_input_file);
	return 0;
      }
      fclose(fp_in);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-o") == 0) {
      --argc;
      ++count;
      command_args->output_file=strdup(argv[count]);
      fp_out = fopen (command_args->output_file, "wb");
      if (fp_out == NULL) {
	fprintf (stderr, "ERROR:   Can not open file '%s' !\r\n",
		 command_args->output_file);
	return 0;
      }
      fclose (fp_out);
      --argc;
      ++count;
    }
    else if (strcmp (argv[count], "-s") == 0) {
      --argc;
      ++count;
      command_args->vad_file=strdup(argv[count]);
      fp_in = fopen (command_args->vad_file,"rb");
      if (fp_in == NULL) {
	fprintf (stderr, "ERROR:   Can not open file '%s' !\r\n",
		 command_args->vad_file);
	return 0;
      }
      fclose(fp_in);
      --argc;
      ++count;
    }
    else if (argv[count][0] == '-') {
      fprintf (stderr, "WARNING:  Un-recognized flag '%s' !\r\n", argv[count]);
      --argc;
      ++count;
    }
    else {
      fprintf (stderr, "WARNING:  Un-recognized argument '%s' !\r\n", argv[count]);
      --argc;
      ++count;
      return 0;
    }
  }
  
  if (!command_args->input_file)    {
    fprintf (stderr, "ERROR:   Input file required !\r\n");
    return 0;
  }
  if (!command_args->output_file)    {
    fprintf (stderr, "ERROR:   Output file required !\r\n");
    return 0;
  }
  
  return 1;
}

void mx_copy(muT *in, muT *out) {
  mx_set_size(in->d1,in->d2,out);
  memcpy(out->data,in->data,in->d1*in->d2*sizeof(float));
}

int main(int argc, char *argv[])
{
  
  
  //"Additive noise" for noise compensation"
  float mean_fft[129] = {36456.537660, 38785.636449, 47851.540806, 235317.838367, 1778419.222305, 12632462.860100, 63156399.997220, 155031404.669803, 202147140.748431, 146738934.243948, 95992588.599821, 113611235.404278, 141717021.904693, 168642666.690051, 182684065.900844, 183638263.941444, 203166239.702228, 200924359.944168, 147524809.963875, 99435945.782690, 84531117.952443, 99309086.398686, 117435836.430613, 95760015.012036, 86367399.408929, 74057426.212751, 46474281.463184, 29792548.539962, 27270238.604906, 38758894.052241, 38966736.521314, 27736752.662543, 22896174.637884, 20665838.317357, 15010347.658873, 16940555.638509, 20738537.339959, 16577490.825954, 13569506.251274, 14653474.336457, 15410956.146418, 13521464.207814, 13934498.711918, 13982734.029577, 12550459.114353, 11390987.611689, 10525124.107672, 11357516.620520, 10918337.708379, 11111424.634228, 12625087.516228, 13266934.796348, 11686931.493795, 10959752.764603, 9002312.073062, 7326061.820102, 7376174.880554, 7654878.732369, 6078436.929397, 5819960.940375, 6256870.348651, 5786006.907036, 5170989.139378, 5031576.437902, 6120223.238114, 7878771.343185, 5503690.699174, 3740251.537677, 3192297.075345, 2935440.268260, 2867165.329350, 2561655.156023, 2441586.579030, 2518175.409656, 2524507.523845, 2457618.147759, 2594372.862476, 2884182.497147, 2800612.627279, 2614581.748001, 2681613.164150, 2907185.303300, 2832736.720152, 2580486.350937, 2278423.604086, 2157134.012783, 2002321.369784, 1777878.069167, 1941015.767542, 2432021.179021, 2327660.185958, 2090634.262948, 1801141.700318, 1582134.321163, 1526353.144699, 1353625.673403, 1195173.242468, 1208410.500565, 1208658.038029, 1156020.434271, 1054502.740807, 996259.039510, 935618.899064, 879644.974343, 793211.423336, 809329.883492, 818901.454534, 782830.376997, 722881.630907, 644103.952388, 572793.884895, 501154.583119, 434373.759655, 373640.443991, 282752.731431, 183822.951444, 118387.194079, 70715.028821, 37820.243419, 16578.954283, 5906.188261, 2209.515849, 2702.646915, 5465.961940, 8400.778187, 12117.993339, 15689.529784, 11751.882060, 6066.191738};
  float mean_mel[23] = { 8232357.887432, 147066375.526554, 410266436.733135, 297556900.793236, 422566107.750971, 569686545.601578, 354249007.690783, 345345895.322708, 214844033.516496, 126552764.692374, 76127827.275885, 69775991.857384, 65703941.849880, 63418926.712335, 56467210.629646, 36807044.026445, 30826963.370200, 19938910.259846, 20665437.149183, 16229492.962456, 11153569.356959, 7152826.220635, 1605035.182299};
  
  long i;
  int cc;
  int start;
  int count;
  long fealen;
  float *hmw;
  float *hmw_big;
  short *inputdata, *noise_inputdata;
  char output_speechfile[2048];
  SPSParams sps_params;
  char *basedirectory;
  MLPParams vad_params;
  
  // Command line stuff
  CommandArgs command_args;
  init_command(&command_args);
  if (!parse_command (argc - 1, argv, &command_args)) {
    exit(1);
  }
  
  basedirectory = getenv("AURORACALC");
  if (!basedirectory) {
    printf("The AURORACALC variable must be set to the root directory where the software is installed!\n");
    exit(1);
  }
  
  command_args.window_length = command_args.window_length*command_args.fs/1000;
  command_args.window_shift = command_args.window_shift*command_args.fs/1000;
  if (command_args.fs == 8000) {
    command_args.fftlen = 256;
  } else if(command_args.fs == 11000) {
    command_args.fftlen = 512;
  } else if(command_args.fs == 16000) {
    command_args.fftlen = 512;
  } else {
    fprintf(stderr,"Sampling frequency not supported\n");
    exit(1);
  }
  
  hmw = (float*)malloc(command_args.window_length*sizeof(float));
  hmw_big = (float*)malloc(512*sizeof(float));
  
  // Fill structure for noise compensation
  sps_params.overest_low = command_args.sps_overest_low;  
  sps_params.overest_up = command_args.sps_overest_up;  
  sps_params.overest_snrlow = 0;
  sps_params.overest_snrup = 20;
  sps_params.h_thresh = command_args.sps_h_thresh;
  sps_params.h_thresh_flag = command_args.sps_h_thresh_flag;
  sps_params.filter_power = command_args.sps_filter_power;
  sps_params.time_alpha = command_args.sps_time_alpha;
  sps_params.time_delay = command_args.sps_time_delay;
  sps_params.freq_length = command_args.sps_freq_length;
  sps_params.freq_flag = command_args.sps_freq_flag;
  sps_params.start_freq_band = ceil(command_args.start_freq*command_args.fftlen/command_args.fs);
  sps_params.h_scale = command_args.sps_h_scale;
  sps_params.addnoise_flag = command_args.sps_addnoise_flag;
  sps_params.addnoise_low = command_args.sps_addnoise_low;
  sps_params.addnoise_up = command_args.sps_addnoise_up;
  sps_params.addnoise_scale = command_args.sps_addnoise_scale;
  if (command_args.sps_type_int==QIO_FFT) {
    sps_params.addnoise_powspec = mean_fft;
  } else if (command_args.sps_type_int==QIO_MEL) {
    if (command_args.n_filters==23) {
      sps_params.addnoise_powspec = mean_mel;
    } else {
      fprintf(stderr,"Don't have a constant power spectrum for this filter-bank size. Switched to no noise addition!!!");
      sps_params.addnoise_flag = 0;
      command_args.sps_addnoise_flag = 0;
    }
  }
  if (command_args.sps_synth) {
    if ( ((float)command_args.window_length/(float)command_args.window_shift) != 2.0 ) {
      fprintf(stderr,"!!!! Windows shoulb be half overlaping, or 3/4 overlaping, for resynthesis !!!! \n");
    }
  }
  
  // DCT coeffs
  fprintf_debug("\nGetting DCT coefs\n");    
  get_dcttable(command_args.n_filters, command_args.n_params, &mx_dct);
  get_dcttable(command_args.n_filters, 6, &mx_vaddct);
  
  // Hamming window
  get_genhanning_window(hmw, command_args.window_length, 0.54);
  get_genhanning_window(hmw_big, 512, 0.54);
  
  // Mel-filter bank weights
  if (command_args.fs == 8000) {
    get_mel_filterbank(command_args.n_filters, command_args.fftlen, command_args.fs, 
		       command_args.start_freq, command_args.end_freq, &mx_mfb);
  } else if (command_args.fs == 11000) {
    get_mel_filterbank(command_args.n_filters, (int)round(command_args.fftlen*8/11), 8000, 
		       command_args.start_freq, command_args.end_freq, &mx_mfb);
  } else {
    get_mel_filterbank(command_args.n_filters, (int)round(command_args.fftlen*8/16), 8000, 
		       command_args.start_freq, command_args.end_freq, &mx_mfb);
  }
  
  // Read the noise file
  if (command_args.noise_input_file) {
    cc=read_raw_file(command_args.noise_input_file, &noise_inputdata, &fealen,command_args.swapin );
    switch(cc) {
    case HTK_FILE_CANT_OPEN:
      fprintf(stderr, "Cannot open input wave file: %s\n", command_args.noise_input_file);
      return-1;
    case HTK_FILE_IO_ERROR:
      fprintf(stderr, "Cannot read input wave file: %s\n", command_args.noise_input_file);
      return-1;
      //case HTK_FILE_SUCCESS:
    }
  }
  
  // Read the speech file
  fprintf(stderr,"%s -> %s\n", command_args.input_file, command_args.output_file);
  cc=read_raw_file(command_args.input_file, &inputdata, &fealen, command_args.swapin);
  switch(cc) {
  case HTK_FILE_CANT_OPEN:
    fprintf(stderr, "Cannot open input wave file: %s\n", command_args.input_file);
    return-1;
  case HTK_FILE_IO_ERROR:
    fprintf(stderr, "Cannot read input wave file: %s\n", command_args.input_file);
    return-1;
    //case HTK_FILE_SUCCESS:
  }
  
  // Convert to Mx form
  mx_read_vec(inputdata, fealen, &mx_samples);
  free(inputdata);
  
  // DC offset compensation for speech
  start = 0;
  if (command_args.dc_offset) {
    float tap = 0.05;
    if (command_args.fs == 8000)
      DCOffsetFilter(mx_samples.data, fealen, &start, fealen, tap);
    else if (command_args.fs == 11000)
      DCOffsetFilter(mx_samples.data, fealen, &start, fealen, tap*8/11);
    else if (command_args.fs == 16000)
      DCOffsetFilter(mx_samples.data, fealen, &start, fealen, tap*8/16);
  }
  
  // FFT for speech
  mx_specgram_phase(&mx_samples, command_args.window_length, command_args.window_shift, 
		    command_args.fftlen, hmw, &mx_out, 0.0, &mx_phase);
  mx_release(&mx_in);
  mx_in = mx_out;
  
  //Additional energy parameters
  //Changed on 01/08/2002 - Sunil
  if(command_args.fs != 8000 ) {
    get_upper_band_energy(&mx_in,command_args.fs,&mx_upper_band_energy);
    //Append them to original specgram
    mx_set_size(mx_in.d1,129,&mx_out);
    for(i=0;i<mx_in.d1;i++)
      memcpy(mx_out.data+i*mx_out.d2,mx_in.data+i*mx_in.d2,mx_out.d2*sizeof(float));
    mx_release(&mx_in);
    mx_in = mx_out;
    mx_cat_vectors(&mx_in, &mx_upper_band_energy, &mx_out);
    mx_release(&mx_in);
    mx_release(&mx_upper_band_energy);
    mx_in = mx_out;
    //Keep phase from 0-4kHz
    mx_set_size(mx_in.d1,129,&mx_out);
    for(i=0;i<mx_in.d1;i++)
      memcpy(mx_out.data+i*mx_out.d2,mx_phase.data+i*mx_phase.d2,mx_out.d2*sizeof(float));
    mx_release(&mx_phase);
    mx_phase = mx_out;
  }

  // Noise compensation on the FFT bins
  if (command_args.sps) {
    if (command_args.sps_type_int==QIO_FFT) {
      if (command_args.sps_hirschnoise) {
	// mx_noiselevel_hirsch(&mx_in, &mx_noise);
      } else {
	if (command_args.noisest_alpha==1.0) {
	  mx_noise_level(&mx_in, &mx_noise);
	}
	else {
	  mx_noise_level_energyfloor(&mx_in, &mx_noise, command_args.noisest_threshold,
				     command_args.noisest_alpha); 
	}
      }
      noise_compensation(&mx_in, &mx_noise, &sps_params, &mx_out, &mx_cleannoisyratio);
      mx_release(&mx_in);
      mx_in = mx_out;
    }
  }

  //Added on 01/08/02 Sunil
  if(command_args.fs != 8000 ) {
    //Separate the upperband features from specgram
    mx_set_size(mx_in.d1, 129, &mx_out); 
    mx_set_size(mx_in.d1, mx_in.d2-129, &mx_upper_band_energy); 
    for(i=0;i<mx_in.d1;i++) {
      memcpy(mx_out.data+i*mx_out.d2,mx_in.data+i*mx_in.d2,mx_out.d2*sizeof(float));
      memcpy(mx_upper_band_energy.data+i*mx_upper_band_energy.d2,mx_in.data+i*mx_in.d2+mx_out.d2,mx_upper_band_energy.d2*sizeof(float));
    }
    mx_release(&mx_in);
    mx_in = mx_out;
  }
  
  // Optional synthesis of speech
  if (command_args.sps_synth) {
    fprintf_debug("Synthesising speech\n");
    speech_synthesis(&mx_in,&mx_phase,&mx_synth,command_args.window_length,
		     command_args.window_shift);
    strcpy(output_speechfile, command_args.output_file);
    strcat(output_speechfile,".clean.raw");
    write_wav_file(output_speechfile,&mx_synth,command_args.swapout);
  }
    
  // Mel-warping
  if (command_args.mel) {
    mx_prod_high_samp(&mx_in,command_args.fs, &mx_mfb, &mx_out);
    mx_release(&mx_in);
    mx_in = mx_out;
  } else {
    // Write output //
    fprintf(stderr,"Number of output frames: %d\n\n", mx_in.d1);
    switch(mx_write_htk_file(command_args.output_file, &mx_in, 100000, 9 ,command_args.swapout)) {
    case HTK_FILE_CANT_OPEN:
      fprintf(stderr, "Cannot open output feature file: %s\n", command_args.output_file);
      break;
    case HTK_FILE_IO_ERROR:
      fprintf(stderr, "Cannot write to output feature file: %s\n", command_args.output_file);
      break;
    case HTK_FILE_INVALID_HEADER:
      fprintf(stderr, "Fatal error: Trying write nonsense header to output\n");
      return -1;
      //case HTK_FILE_SUCCESS:
    }
    exit (0);
  }
  
  //Log
  mx_log(&mx_in);                                    
  if( (command_args.fs != 8000) )
    mx_log(&mx_upper_band_energy);
  
  // Reading silence probabilities & energy-based VAD
  if (command_args.vad) {
    if (command_args.vad_file) {
      fprintf(stderr,"Reading silence features\n");
      silfea = read_sil_fea_down_up(command_args.vad_file, mx_in.d1,command_args.down_up);
    } else if (command_args.vad_weights && command_args.vad_norm) {
      strcpy(temp2,basedirectory);
      strcat(temp2,"/parameters/vad/vad-lpf.fil");
      //lda_filter(temp2, &mx_in, &mx_invad);
      //Changed to vadlp_filter - Sunil 01/04/2002
      vadlp_filter(&mx_in, &mx_invad);
      if (command_args.spec_basis_file) {
	// Spectral Basis
	fprintf_debug("Applying Spectral LDA\n");
	get_sldtable(command_args.spec_basis_file, 6, &mx_vaddct);
	mx_prod(&mx_invad, &mx_vaddct, &mx_out);
	mx_release(&mx_invad);
	mx_invad = mx_out;
      } else {
	// DCT
	fprintf_debug("Applying DCT\n");
	mx_prod(&mx_invad, &mx_vaddct, &mx_out);
	mx_release(&mx_invad);
	mx_invad = mx_out;
      }
      init_mlp_forward(&vad_params, command_args.vad_weights, command_args.vad_norm);
      // no silflag, softmaxed output
      bunch_mlp_forward(&mx_invad,&mx_out,vad_params,mx_invad.d1,0,0);
      silfea = (float*)lmalloc((mx_invad.d1)*sizeof(float));
      for (count=0;count<mx_invad.d1;count++)
	silfea[count] = 1.0;
      for (count=0;count<mx_invad.d1;count++) {
	if (command_args.down_up) {
	  if (!(count%2)) {
	    silfea[count/2] = mx_out.data[max(0,count-command_args.vad_shift)*mx_out.d2+1];
	  }
	}
	else {
	  silfea[count] = mx_out.data[max(0,count-command_args.vad_shift)*mx_out.d2+1];
	}
      }
    } else {
      silfea = (float*)lmalloc((mx_in.d1)*sizeof(float));
      for(i = 0; i < mx_in.d1; i++)
	silfea[i]=0;
    }  
  }
  
  // If the LDA is on the terminal, then you don't need to
  // low-pass before downsampling
  // IIR-LDA filter
  fprintf(stderr,"LDA filtering mel-bands:%s\n",command_args.lda_filters_file);  
  lda_filter(command_args.lda_filters_file, &mx_in, &mx_out);
  mx_release(&mx_in);
  mx_in = mx_out;  

  // Keep one more copy of features after LDA
  mx_copy(&mx_in,&mx_traps_in1);
  //Additional energy parameters
  if( (command_args.fs != 8000) ) {
    lda_filter(command_args.lda_filters_file, &mx_upper_band_energy, &mx_out);
    mx_release(&mx_upper_band_energy);
    mx_upper_band_energy = mx_out; 
  }
  if (command_args.dct) {
    if (command_args.spec_basis_file) {
      // Spectral Basis
      fprintf_debug("Applying Spectral LDA\n");
      get_sldtable(command_args.spec_basis_file,command_args.n_params, &mx_dct23by23);
      mx_prod(&mx_in, &mx_dct, &mx_out);
      mx_release(&mx_in);
      mx_in = mx_out;
    } else {
      // DCT
      fprintf_debug("Applying DCT\n");
      mx_prod(&mx_in, &mx_dct, &mx_out);
      mx_release(&mx_in);
      mx_in = mx_out;
    }   
  }
  //Augment cepstral features with higher-band energies
  if((command_args.fs != 8000) ) {
    mx_cat_vectors(&mx_in, &mx_upper_band_energy, &mx_out);
    mx_release(&mx_in);
    mx_release(&mx_upper_band_energy);
    mx_in = mx_out;
  } else {
	//Put two zeros
	mx_set_size(mx_in.d1,2,&mx_upper_band_energy);
	mx_cat_vectors(&mx_in, &mx_upper_band_energy, &mx_out);
	mx_release(&mx_in);
    mx_release(&mx_upper_band_energy);
    mx_in = mx_out;
  }
  if (command_args.down_up) {
    // Downsampling pre-filter
    if (command_args.down_filter) {
      fprintf_debug("Downsampling pre-filter\n");
      lp_filter(&mx_in, &mx_out);
      mx_release(&mx_in);
      mx_in = mx_out;
    }
    // Downsample by 2
    fprintf_debug("Downsampling by 2\n");
    downsample2(&mx_in, mx_in.d2, &mx_out);
    mx_release(&mx_in);
    mx_in = mx_out;
  }  

  //Pack the silence probability with features
  if (command_args.vad) {
    mx_set_size(mx_in.d1, 1, &mx_silfea);
    memcpy(mx_silfea.data, silfea, mx_in.d1 * sizeof(float));
    mx_cat_vectors(&mx_silfea, &mx_in, &mx_out);
    mx_release(&mx_silfea);
    mx_release(&mx_in);
    mx_in = mx_out;
    free(silfea);
  }
  
  // Write output //
  fprintf(stderr,"Features: %d x %d\n\n", mx_in.d1,mx_in.d2);
  switch(mx_write_htk_file(command_args.output_file, &mx_in, 100000, 9 ,command_args.swapout)) {
  case HTK_FILE_CANT_OPEN:
    fprintf(stderr, "Cannot open output feature file: %s\n", command_args.output_file);
    break;
  case HTK_FILE_IO_ERROR:
    fprintf(stderr, "Cannot write to output feature file: %s\n", command_args.output_file);
    break;
  case HTK_FILE_INVALID_HEADER:
    fprintf(stderr, "Fatal error: Trying write nonsense header to output\n");
    return -1;
    //case HTK_FILE_SUCCESS:
  }
  
  // Release memory
  mx_release(&mx_in);
  mx_release(&mx_dct);
  mx_release(&mx_mfb);
  mx_release(&mx_delta);
  mx_release(&mx_deltadelta);
  mx_release(&mx_phase);
  mx_release(&mx_synth);
  mx_release(&mx_noise);
  
  return(0);
  
}

