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
#include "libicsi/klt.h"

#define ON_ALPHA 0.1
#define ON_BIAS  1

//Mx structures
muT mx_dct;
muT mx_idct;
muT mx_dct23x23;
muT mx_idct23x23;
muT mx_dct15x15;
muT mx_idct15x15;

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
muT mx_trapsdct;
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
int from_tandem = 9;
int from_traps = 6;

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
  int vad;
  float vadthresh_fd;
  float vadthresh_oln;
  int fd;
  int vrr;
  int vrp;
  
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

// gloabl variables to be used in server side
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
  fprintf (stderr, "     -fd 0/1           Frame drop flag (%d)       \r\n",command_args->fd);
  fprintf (stderr, "     -VADfd float      VAD threshold for frame dropping (%f)  \r\n", command_args->vadthresh_fd);
  fprintf (stderr, "     -VADoln float     VAD threshold for online normalization (%f)  \r\n", command_args->vadthresh_oln);
  fprintf (stderr, "                       If not provided, do not drop any frames.\r\n");
  fprintf (stderr, "     -V 0/1            Use VAD  (%d)\r\n", command_args->vad);  
  fprintf (stderr, "     -VRR <int>        VAD rank range    (%d)       \r\n",command_args->vrr);
  fprintf (stderr, "     -VRP <int>        VAD rank point    (%d)       \r\n",command_args->vrp);
  
  fprintf (stderr, "FEATURES DOWNSAMPLING/UPSAMPLING\n");
  fprintf (stderr, "     -D 0/1            Use downsampling/upsampling (%d)              \r\n", command_args->down_up);
  fprintf (stderr, "     -Dfilter 0/1      Use downsampling filter - you may want to set this to 0 if the LDA filter already cuts a lot of the high mod. freqs. (%d) \r\n\n", command_args->down_filter);
  
  
  fprintf (stderr, "TEMPORAL LDA, SPECTRAL LDA and DCT\n");
  fprintf (stderr, "     -LDAserver 0/1    Temporal LDA done on the server-side (%d)\r\n",command_args->serverside_lda);  
  fprintf (stderr, "     -F filters        Filename of the lda filters coefficients (use the provided void filters if you don't want any filtering)\n                       (%s)\r\n", command_args->lda_filters_file);
  fprintf (stderr, "     -SF filters       Filename of the spectral basis vectors\n");
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
  command_args->vadthresh_fd = 0.5;
  command_args->vadthresh_oln = 0.5;
  command_args->fd = 1;
  command_args->vrr = 11;
  command_args->vrp = 2;
  command_args->vad = 1;
  
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
   else if (strcmp (argv[count], "-TRAPS") == 0) {
      --argc;
      ++count;
      command_args->traps = atoi(argv[count]);
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
  
  x_int32 sampPeriod;
  x_int16 parmKind;
  char *temp;
  char trapsklt[1024];
  long i;
  int cc;
  char *basedirectory;
  int nUpBands;
  
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
  
  // DCT coeffs
  fprintf_debug("\nGetting DCT coefs\n");    
  if ( command_args.upper_flag ) {
    if ( command_args.fs == 11000 ) {
      get_dcttable_s(command_args.n_filters+1, command_args.n_params, &mx_dct, &mx_idct);
    }
    else if ( command_args.fs == 16000 ) {
      get_dcttable_s(command_args.n_filters+2, command_args.n_params, &mx_dct, &mx_idct);
    }
  }
  else {
    get_dcttable_s(command_args.n_filters, command_args.n_params, &mx_dct, &mx_idct);
  }    
  get_dcttable_s(command_args.n_filters, command_args.n_filters, &mx_dct23x23, &mx_idct23x23);
  get_dcttable_s(command_args.n_params, command_args.n_params, &mx_dct15x15, &mx_idct15x15);
  get_dcttable(command_args.n_filters, 6, &mx_vaddct);
  
  ////////////////////////////
  // SERVER SIDE PROCESSING //
  ////////////////////////////
  
  // Read the speech file
  fprintf(stderr,"%s -> %s\n", command_args.input_file, command_args.output_file);
  cc=mx_read_htk_file(command_args.input_file, &mx_in, &sampPeriod, &parmKind,command_args.swapin); 
  switch(cc) {
  case HTK_FILE_CANT_OPEN:
    fprintf(stderr, "Cannot open input wave file: %s\n", command_args.input_file);
    return-1;
  case HTK_FILE_IO_ERROR:
    fprintf(stderr, "Cannot read input wave file: %s\n", command_args.input_file);
    return-1;
    //case HTK_FILE_SUCCESS:
  }
  
  //Unpack silence probabilities and cepstral coeffs
  mx_set_size(mx_in.d1, 1, &mx_silfea);
  for(i = 0; i < mx_in.d1; i++) {
    mx_silfea.data[i] = mx_in.data[i*mx_in.d2];
  }
  nUpBands = (command_args.upper_flag && command_args.fs == 16000) ? 2 :
             (command_args.upper_flag && command_args.fs == 11000) ? 1 : 0;
  
  mx_set_size(mx_in.d1, command_args.n_params+nUpBands, &mx_out);
  for(i = 0; i < mx_in.d1; i++) {
    memcpy(mx_out.data+i*mx_out.d2, mx_in.data+i*mx_in.d2+1, mx_out.d2 * sizeof(float));
  }
  mx_release(&mx_in);
  mx_in = mx_out;
  printf("%f %f\n",mx_in.data[0],mx_in.data[1]);
  // VAD probabilities processing //
  if (command_args.vad) {
    if (command_args.down_up) {
      // Up sampling features
      fprintf_debug("Up sampling Silprob\n");
      upsample2(&mx_silfea, &mx_out);
      mx_release(&mx_silfea);
      mx_silfea = mx_out;
    }
    silfea = (float*)lmalloc((mx_silfea.d1)*sizeof(float));
    for(i = 0; i < mx_silfea.d1; i++) {
      silfea[i] = mx_silfea.data[i*mx_silfea.d2];
    }
  }
  // Convert for the OLN, no median filter smoothing //
  if(command_args.vad) {
    silflag_oln = sil_prob_to_flag(mx_silfea.d1,silfea,command_args.vadthresh_oln,0,0);
    // Convert VAD probabilities to binary decisions
    silflag = sil_prob_to_flag(mx_silfea.d1,silfea,command_args.vadthresh_fd,2*command_args.vrr-1,2*command_args.vrp+1);
  }

  // upsampling of features //
  
  // mx_in has 15 ceps till here ..

  // Straight features and tandem features processing //
  server_straight(command_args);

  if (command_args.prop2) {
    //Tandem processing
    server_tandem(command_args);
  }
  
  // TRAPS features computation : pratibha , Dec 24 2001 //
  if (command_args.traps) {
    //IDCT for reconstructing critical bands from cepstra
    fprintf_debug("Applying IDCT\n");
    mx_prod(&mx_trapsdct, &mx_idct15x15, &mx_traps_recons_cb);
    mx_release(&mx_trapsdct);
    
    // Reconstructed critical bands 
    server_traps(command_args);
  }
  
  // Concatenation of different streams
  if ( command_args.traps == 1 && command_args.prop2 == 0 ) {
    mx_sub_vectors(&mx_straight, &mx_out, from_straight);
    mx_release(&mx_straight);
    mx_straight = mx_out;
    // Cat straight features + TRAPs
    mx_cat_vectors(&mx_straight, &mx_out_final_traps, &mx_out);
//    strcpy(junkname,"/net/samson/u0/pratibha/Junk_concate_noklt");
//    mx_write_htk_file(junkname, &mx_out,100000, 9 ,command_args.swapout);

    mx_release(&mx_straight);
    mx_release(&mx_out_final_traps);
    temp = getenv("AURORACALC");
    strcpy(trapsklt,temp);
    strcat(trapsklt,"/parameters/traps/total_stats_dim51_traps6_straight45");
    mx_klt_forward(&mx_out, &mx_out_final_traps_klt, trapsklt, 0, 1, 1, 51, 0);

    fprintf(stderr,"concatenate traps and feat\n");
    mx_in = mx_out_final_traps_klt;
  }
  else if (command_args.traps == 0 && command_args.prop2 == 1) {
    // Keep only from_straight*3 straight features - to end up with from_straight*3 + from_tandem
    mx_sub_vectors(&mx_straight, &mx_out, from_straight); 
    mx_release(&mx_straight); 
    mx_straight = mx_out;
    // Cat straight features + MLP features
    mx_cat_vectors(&mx_straight, &mx_mlpout, &mx_out); 
    mx_release(&mx_in); 
    mx_in = mx_out;
  }
  else if (command_args.traps == 1 && command_args.prop2 == 1) {
    // Cat straight features + TRAPs + Feature Net
    mx_sub_vectors(&mx_straight, &mx_out, from_straight);
    mx_release(&mx_straight);
    mx_straight = mx_out;
    mx_cat_vectors(&mx_straight, &mx_out_final_traps, &mx_traps_baseline);
    mx_release(&mx_out_final_traps);
    fprintf(stderr,"After concate TRAPs + Baseline\n");
    mx_cat_vectors(&mx_traps_baseline, &mx_mlpout, &mx_out);
    fprintf(stderr,"After concate TRAPs + Baseline + FeatNet\n");
    //mx_release(&mx_in);
    mx_in = mx_out;
  } 
  else {
    mx_in = mx_straight;
  }
  
  // Frame dropping
  if (command_args.fd) {
    fprintf(stderr,"Frame dropping\n");
    framedropping(&mx_in,silflag,0,&mx_out);
    mx_release(&mx_in); 
    mx_in = mx_out;
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

////////pratibha , 24th Dec 2001 ///////
void server_traps(CommandArgs command_args) {
  int B ;
  char *temp;
  char bandclassifier_param[1024];
  char bandclassifier_weight[1024];
  char bandclassifier_norm[1024];
  char mergerclassifier_norm[1024];
  char mergerclassifier_weight[1024];
  MLPParams bandmlp_params[15];
  MLPParams mergermlp_params;
  char Band[128];
  char tempb[1024],temp_w[1024],temp_n[1024];

  temp = getenv("AURORACALC");
  strcpy(tempb,temp);
  strcat(tempb,"/parameters/traps/tnn_60_1bands_crb");
  
  mx_set_size(1,60,&hamming_win) ;
  HammingTrap (60,&hamming_win);



  for (B = 0 ; B < 15 ; B++ ) {
       sprintf(Band,"%d",B);
       mx_set_size(mx_traps_recons_cb.d1, 60, &mx_in_band_traps) ;
       mx_create_bandclassifier_input60(&mx_traps_recons_cb,&mx_in_band_traps,B,&hamming_win);
       strcpy(bandclassifier_param,tempb);
       strcat(bandclassifier_param,Band);
       strcpy(temp_w,bandclassifier_param);
       strcpy(temp_n,bandclassifier_param);

       strcpy(bandclassifier_weight,temp_w);
       strcpy(bandclassifier_norm,temp_n);
       strcat(bandclassifier_weight,".weights");
       strcat(bandclassifier_norm,".norms");
       
       init_mlp_forward(&bandmlp_params[B],bandclassifier_weight,bandclassifier_norm);
       bunch_mlp_forward_linear(&mx_in_band_traps,&mx_out_band_traps[B],bandmlp_params[B],mx_in_band_traps.d1,silflag);   
       mx_release(&mx_in_band_traps);
  }

  fprintf(stderr,"Combine traps through merging net\n");
  mx_combine_traps(mx_out_band_traps,&mx_concate_comb_traps,15);
  
  for (B = 0 ; B < 15 ; B++ ) 
  mx_release(&mx_out_band_traps[B]);
	       
  // Forward pass through merging net
  strcpy(mergerclassifier_weight,temp);
  strcpy(mergerclassifier_norm,temp);
  strcat(mergerclassifier_weight,"/parameters/traps/merger_weight_60");
  strcat(mergerclassifier_norm,"/parameters/traps/merger_norm_60");
  init_mlp_forward(&mergermlp_params,mergerclassifier_weight,mergerclassifier_norm);
  bunch_mlp_forward_linear(&mx_concate_comb_traps,&mx_out_final_traps,mergermlp_params,mx_concate_comb_traps.d1,silflag);
  fprintf(stderr,"Input to Merger %d %d\n",mx_concate_comb_traps.d1,mx_concate_comb_traps.d2);
  mx_release(&mx_concate_comb_traps);
  
  fprintf(stderr,"Merger output %d %d\n",mx_out_final_traps.d1,mx_out_final_traps.d2);


 //Upsample TRAPs features : 28th Jan .. no upsampling as straight features are already upsampled befor
 // reconstruction of critical bands

  fprintf(stderr,"After TRAPs\n");
}

void server_tandem(CommandArgs command_args) {
  
  char *temp;
  char tandemnet[1024];
  char tandemnorm[1024];
  MLPParams mlp_params;
  
    // Online normalization with VAD
  if (command_args.alpha > 0.0) {
    fprintf_debug("Online normalization\n");
    if (command_args.onl_init==1) {
      fprintf(stderr,"OnlInit not suported yet");
    }
    online_normalization_with_vad(&mx_mlpout, 0.5*command_args.alpha, ON_BIAS, &mx_out, command_args.onl_latency, silflag_oln);
    mx_release(&mx_mlpout);
    mx_mlpout = mx_out;
  }


    // Calculate deltas
  if (command_args.delta_order == 2) {
    mx_ddynamic(&mx_mlpout,&mx_deltadelta,3);
  }


  fprintf(stderr,"Calculating deltas\n");
  if (command_args.delta_order >= 1) {
    mx_dynamic_htk(&mx_mlpout,&mx_delta,3);
    mx_cat_vectors(&mx_mlpout, &mx_delta, &mx_out);
    mx_release(&mx_mlpout);
    mx_mlpout = mx_out;
  }


  if (command_args.delta_order == 2) {
    mx_cat_vectors(&mx_mlpout, &mx_deltadelta, &mx_out);
    mx_release(&mx_mlpout);
    mx_mlpout = mx_out;
  }


  downsample2(&mx_mlpout, mx_mlpout.d2, &mx_out);
  mx_release(&mx_mlpout);
  mx_mlpout = mx_out;

  // MLP stream
  // Warning: the weight file must contain the number of input,
  // hidden and output units (first 3 integers in the file)
  // the norm file is the standard format
  temp = getenv("AURORACALC");
  strcpy(temp2,temp);
  strcat(temp2,"/parameters/mlp/Release_frozen_icsi56_DS0.head");
  silence_state = 54;
  strcpy(tandemnet,temp2);
  strcpy(temp2,temp);
  strcat(temp2,"/parameters/mlp/Release_frozen.norm");
  strcpy(tandemnorm,temp2);
  init_mlp_forward(&mlp_params, tandemnet, tandemnorm);
  bunch_mlp_forward_linear(&mx_mlpout,&mx_out,mlp_params,mx_mlpout.d1,silflag);
  mx_release(&mx_mlpout);
  mx_mlpout = mx_out;
  
  // correct speech/silence discrimination
  mx_correctspeechsilence(&mx_mlpout, &mx_out, silfea, silence_state);
  mx_release(&mx_mlpout);
  mx_mlpout = mx_out;
  
  // per-frame normalization
  mx_perframenorm(&mx_mlpout, &mx_out);
  mx_release(&mx_mlpout);
  mx_mlpout = mx_out;
  
  // Upsample because MLP output computed with a skip of 1 frame
  upsample2(&mx_mlpout, &mx_out);
  mx_release(&mx_mlpout);
  mx_mlpout = mx_out;

}

void server_straight(CommandArgs command_args) {
  
  int i;
  printf("%f %f\n",mx_in.data[0],mx_in.data[1]);



  //Seperate the cepstrum and upper band energy at the server side
  if((command_args.upper_flag) && (command_args.fs != 8000)) {
    mx_set_size(mx_in.d1,command_args.n_params,&mx_out); 
    mx_set_size(mx_in.d1,mx_in.d2-command_args.n_params, &mx_upper_band_energy); 
    for(i=0;i<mx_in.d1;i++) {
      memcpy(mx_out.data+i*mx_out.d2,mx_in.data+i*mx_in.d2,mx_out.d2*sizeof(float));
      memcpy(mx_upper_band_energy.data+i*mx_upper_band_energy.d2,mx_in.data+i*mx_in.d2+mx_out.d2,mx_upper_band_energy.d2*sizeof(float));
    }
    mx_release(&mx_in);
    mx_straight = mx_out;
  }
  else
  mx_straight = mx_in ;

   fprintf(stderr,"Straight Feat %d %d\n",mx_straight.d1,mx_straight.d2);

  // REMEMBER Copy for TRAPS if want to  work with downsampled traps (TRAPS with 30 point input)
  //mx_copy(&mx_straight,&mx_trapsdct);


    if (command_args.down_up) {
    // Up sampling features
    fprintf_debug("Up sampling features\n");
    upsample2(&mx_straight, &mx_out);
    mx_release(&mx_straight);
    mx_straight = mx_out;
    printf("%15.10f %15.10f\n",mx_straight.data[0],mx_straight.data[1]);
    // Upsampling post-filter
    fprintf_debug("Downsampling post-filter\n");
    lp_filter(&mx_straight, &mx_out);
    mx_release(&mx_straight);
    mx_straight = mx_out;
  }

  mx_copy(&mx_straight,&mx_trapsdct);
  mx_copy(&mx_straight,&mx_mlpout);
  
  // Append 0 to straight
  if((command_args.upper_flag) && (command_args.fs != 8000)) {
    mx_set_size(mx_straight.d1,command_args.n_filters,&mx_out);
    for(i=0;i<mx_straight.d1;i++) {
      memcpy(mx_out.data+i*mx_out.d2,mx_straight.data+i*mx_straight.d2,mx_straight.d2*sizeof(float));
    }
    mx_release(&mx_straight);
    mx_straight = mx_out;
    
    //Inv DCT on cepstral values
    fprintf_debug("Applying Inverse DCT\n");
    mx_prod(&mx_straight, &mx_idct23x23, &mx_out);
    mx_release(&mx_straight);
    mx_straight = mx_out;
    
    if (command_args.down_up) {
    // Up sampling features
    fprintf_debug("Up sampling upper_band_energy features\n");
    upsample2(&mx_upper_band_energy, &mx_out);
    mx_release(&mx_upper_band_energy);
    mx_upper_band_energy = mx_out;
    // Upsampling post-filter
    fprintf_debug("Downsampling post-filter\n");
    lp_filter(&mx_upper_band_energy, &mx_out);
    mx_release(&mx_upper_band_energy);
    mx_upper_band_energy = mx_out;
  }

    //Pack the spectral values
    fprintf_debug("Packing spectral values\n");
    mx_cat_vectors(&mx_straight, &mx_upper_band_energy, &mx_out);
    mx_release(&mx_straight);
    mx_release(&mx_upper_band_energy);
    mx_straight = mx_out;

    //DCT again
    fprintf(stderr,"Applying DCT on 0-%d KHz spectrum\n",command_args.fs/2);
    mx_prod(&mx_straight, &mx_dct, &mx_out);
    mx_release(&mx_straight);
    mx_straight = mx_out;
  }

  // Online normalization with VAD
  if (command_args.alpha > 0.0) {
    fprintf_debug("Online normalization\n");
    if (command_args.onl_init==1) {
      fprintf(stderr,"OnlInit not suported yet");
    }
    online_normalization_with_vad(&mx_straight, 0.5*command_args.alpha, ON_BIAS, &mx_out, command_args.onl_latency, silflag_oln);
    mx_release(&mx_straight); 
    mx_straight = mx_out;
  } 
  
  // Calculate deltas
  if (command_args.delta_order == 2) {
    mx_ddynamic(&mx_straight,&mx_deltadelta,3);
  }
  
  fprintf(stderr,"Calculating deltas\n");
  if (command_args.delta_order >= 1) {
    mx_dynamic_htk(&mx_straight,&mx_delta,3);
    mx_cat_vectors(&mx_straight, &mx_delta, &mx_out); 
    mx_release(&mx_straight); 
    mx_straight = mx_out;
  }
  
  if (command_args.delta_order == 2) {
    mx_cat_vectors(&mx_straight, &mx_deltadelta, &mx_out);
    mx_release(&mx_straight);
    mx_straight = mx_out;
  }

  fprintf(stderr,"Straight %d %d\n",mx_straight.d1,mx_straight.d2);
}
