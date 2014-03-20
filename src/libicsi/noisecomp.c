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

#include <stdlib.h>
#include <string.h>
#include "noisecomp.h"

//------------------------------------------------------------------------------------//
//Function: noise_compensation                                                        //
//Input: Noisy power spectrum, Noise estimate, Noise compensation parameter structure //
//Output: Estimated clean power spectrum                                              //
//------------------------------------------------------------------------------------//
int noise_compensation(muT *in, muT *nse, SPSParams *ncpars, muT *out, muT *out2)
{
  int i,j;
  int curr_frmindex;
  muT instfilt;
  muT tsfilt;
  muT finalfilt;
  float noise_power, addnoise_power;
  float clean_en, noisy_en;
  addnoise_power = 0;

  if(mx_set_size(in->d1, in->d2, &instfilt) || 
     mx_set_size(1, in->d2, &tsfilt) || 
     mx_set_size(1, in->d2, &finalfilt) || 
     mx_set_size(in->d1, in->d2, out) ||
     mx_set_size(in->d1, 1, out2) )
    {
      return -1;
    }

  //for (i=0;i<in->d2;i++) {
  //tsfilt.data[i] = pow(ncpars->h_thresh,ncpars->filter_power);
  //}
  
  if(ncpars->addnoise_flag) {
    if (ncpars->addnoise_powspec) {
      addnoise_power = 0;
      for (j=0;j<in->d2;j++) {
		addnoise_power += ncpars->addnoise_powspec[j];
      }
    }
    else {
      fprintf(stderr,"Constant addition after noise compensation requires constant power spectrum!");
      exit(1);
    }
  }

  for (i=0; i < in->d1; i++) 
    {
      //Time smooth noise reduction filter
      estimate_timesmooth_filter(in,nse,i,ncpars,&instfilt,&tsfilt);
      if (i >= ncpars->time_delay) {
		//Estimate the final freqsmooth filter after the delay
	        //estimate_freqsmooth_filt(i,ncpars,&instfilt,&tsfilt,&finalfilt);
		estimate_freqsmooth_filt_high_samp(i,ncpars,&instfilt,&tsfilt,&finalfilt);
		curr_frmindex = i - ncpars->time_delay;
		//Apply the filter
		noise_power = 0;
		for (j=0;j<in->d2;j++) {
		  noise_power += nse->data[curr_frmindex*in->d2+j];
		}
		noise_comp_filter(in,&finalfilt,curr_frmindex,nse,ncpars,out,noise_power,addnoise_power,&clean_en,&noisy_en);
		out2->data[curr_frmindex] = clean_en/noisy_en;
      }
    }
  //For the last few frames
  for (i=in->d1-ncpars->time_delay;i < in->d1; i++) {
    noise_comp_filter(in,&finalfilt,i,nse,ncpars,out,noise_power,addnoise_power,&clean_en,&noisy_en);
    out2->data[i] = clean_en/noisy_en;
  }

  return 0;
  
}

//------------------------------------------------------------------------------------//
//Function: estimate_timesmooth_filter                                                //
//Input: Noisy power spectrum, Noise estimate, Noise compensation parameter structure //
//       frame index.                                                                 //
//Output: instantaneous filter and time smoothed filter                               //
//------------------------------------------------------------------------------------//
int estimate_timesmooth_filter(muT *in,muT *nse,int frmind, SPSParams *ncpars,muT *instfilt,muT *tsfilt)
{
  int i;
  long frmptr = frmind*in->d2;
  float post_snr;
  float signal_energy=0;
  float noise_energy=0;
  float threshold;
  float over_est_factor;

  //Compute noisy signal energy and noise energy
  //for (i=0; i< in->d2; i++) {
  //Changed to only 0-4K range even for 11K and 16K sampling frequencies
  for (i=0; i< 129; i++) {
    signal_energy += in->data[frmptr+i];
    noise_energy += nse->data[frmptr+i];
  }
  //Compute posterior snr
  post_snr = 10*log10(signal_energy/noise_energy);
  
  //Compute noise over estimation factor
  if (post_snr < ncpars->overest_snrlow)
    {
      post_snr = ncpars->overest_snrlow;
    } else if (post_snr > ncpars->overest_snrup) 
    {
      post_snr = ncpars->overest_snrup;
    }
  //over_est_factor= (-1.875/20)*post_snr+3.125;
  over_est_factor = ((ncpars->overest_up - ncpars->overest_low)/ncpars->overest_snrup)*post_snr+ncpars->overest_low;

  for (i=0; i< in->d2; i++) {
    if (ncpars->h_thresh_flag) {
      //Frequency dependent threshold
      threshold = ncpars->h_thresh*noise_energy/(in->d2*nse->data[frmptr+i]);
    } else {
      //Fixed threshold
      threshold = ncpars->h_thresh;
    }
    //Gain equation
    instfilt->data[frmptr+i] = MAX(((in->data[frmptr+i] - over_est_factor*nse->data[frmptr+i])/in->data[frmptr+i]),threshold);
    //filterpower = 1 for spectral subtraction and 2 for wiener filter
    instfilt->data[frmptr+i] = (float)pow(instfilt->data[frmptr+i],ncpars->filter_power);
    //Time smooth
    tsfilt->data[i] = (1-ncpars->time_alpha)*tsfilt->data[i] + ncpars->time_alpha*instfilt->data[frmptr+i];
  } 
  return 0;
}

//------------------------------------------------------------------------------------//
//Function: estimate_freqsmooth_filter                                                //
//Input: Frame index, Noise compensation parameter structure, instantaneous filter    //
//       and time smoothed filter                                                     //
//Output: Final filter                                                                //
//------------------------------------------------------------------------------------//
int estimate_freqsmooth_filt(int frmind, SPSParams *ncpars, muT *instfilt, muT *tsfilt, muT *finalfilt)
{
  int i, j;
  int lidx;
  int uidx;
  int winsize = (ncpars->freq_length-1)/2;
  float *ts2_filt;
  ts2_filt = (float*)malloc(tsfilt->d2*sizeof(float));

  for (i=0; i< tsfilt->d2; i++) 
    //Second level of gain smoothing
    ts2_filt[i] = (instfilt->data[(frmind - ncpars->time_delay)*instfilt->d2+i]*tsfilt->data[i]) + \
      (tsfilt->data[i]*(1-tsfilt->data[i]));
    
  //Frequency smoothing
  for (i=0; i< tsfilt->d2; i++) {
    //0 = arithmetic; 1 = geometric
    finalfilt->data[i] = ncpars->freq_flag;
	//Changed to ncpars->start_freq_band to skip the initial unused bands
    lidx = MAX((i-winsize),ncpars->start_freq_band);
	//lidx = MAX((i-winsize),0);
    uidx = MIN((i+winsize),tsfilt->d2-1);
    for (j=lidx; j<=uidx; j++) {
      if(ncpars->freq_flag) {
		finalfilt->data[i] *= ts2_filt[j];
      } else {
		finalfilt->data[i] += ts2_filt[j];
      }
    }

    if(ncpars->freq_flag) {
      finalfilt->data[i] = (float)(pow(finalfilt->data[i],1/(uidx-lidx+1)));
    } else {
       finalfilt->data[i] /= (uidx-lidx+1);
    }
  }

  free (ts2_filt);
  return 0;
}


//------------------------------------------------------------------------------------//
//Function: estimate_freqsmooth_filter_high_samp                                      //
//Input: Frame index, Noise compensation parameter structure, instantaneous filter    //
//       and time smoothed filter                                                     //
//Output: Final filter                                                                //
//------------------------------------------------------------------------------------//
int estimate_freqsmooth_filt_high_samp(int frmind, SPSParams *ncpars, muT *instfilt, muT *tsfilt, muT *finalfilt)
{
  int i, j;
  int lidx;
  int uidx;
  int winsize = (ncpars->freq_length-1)/2;
  float *ts2_filt;
  ts2_filt = (float*)malloc(tsfilt->d2*sizeof(float));

  for (i=0; i< tsfilt->d2; i++) 
    //Second level of gain smoothing
    ts2_filt[i] = (instfilt->data[(frmind - ncpars->time_delay)*instfilt->d2+i]*tsfilt->data[i]) + \
      (tsfilt->data[i]*(1-tsfilt->data[i]));
    
  //Frequency smoothing 0-4K
  for (i=0; i< 129; i++) {
    //0 = arithmetic; 1 = geometric
    finalfilt->data[i] = ncpars->freq_flag;
	//Changed to ncpars->start_freq_band to skip the initial unused bands
    //lidx = MAX((i-winsize),ncpars->start_freq_band);
	lidx = MAX((i-winsize),0);
    uidx = MIN((i+winsize),128);
    for (j=lidx; j<=uidx; j++) {
      if(ncpars->freq_flag) {
		finalfilt->data[i] *= ts2_filt[j];
      } else {
		finalfilt->data[i] += ts2_filt[j];
      }
    }

    if(ncpars->freq_flag) {
      finalfilt->data[i] = (float)(pow(finalfilt->data[i],1/(uidx-lidx+1)));
    } else {
       finalfilt->data[i] /= (uidx-lidx+1);
    }
  }
  
  for (i=129; i< finalfilt->d2; i++) 
	finalfilt->data[i] = ts2_filt[i];
  
  free (ts2_filt);
  return 0;
}

//------------------------------------------------------------------------------------//
//Function: noise_comp_filter                                                         //
//Input: noisy power spectrum, filter, frame index, noise estimate, noise compensation//
//       parameters                                                                   //
//Output: Estimated clean power spectrum                                              //
//------------------------------------------------------------------------------------//
int noise_comp_filter(muT *in, muT *filt, int frmind, muT *nse, SPSParams *ncpars, muT *out, float noise_power, float addnoise_power, float *clean_en, float *noisy_en)
{

  float scale;
  int i;
  
  *clean_en = 0;
  *noisy_en = 0;

  for (i=0; i< in->d2; i++) {
    *noisy_en += in->data[frmind * in->d2 + i];
    out->data[frmind * in->d2 + i] = MAX((filt->data[i]*in->data[frmind * in->d2 + i]),ncpars->h_scale*nse->data[frmind * in->d2 + i]);
    *clean_en += out->data[frmind * in->d2 + i];
    //out->data[frmind * in->d2 + i] = in->data[frmind * in->d2 + i];
    if(ncpars->addnoise_flag) {
      scale = MIN(ncpars->addnoise_up,MAX(ncpars->addnoise_low,ncpars->addnoise_scale*noise_power/addnoise_power));
      out->data[frmind * in->d2 + i] += scale*ncpars->addnoise_powspec[i];
    }
  }
  
  return 0;
}

//------------------------------------------------------------------------------------//
//Function: minima_tracking_noise_estimation                                          //
//------------------------------------------------------------------------------------//
int minima_tracking_noise_estimation(muT *mx_P_n,muT *mx_in)
{
  muT mx_x_k, mx_P_x, mx_Q_k;
  muT mx_mfb_norm, minimum_val;
  int i, j, k, beg_num_mean, count;
  float help_fl; 

  if(mx_set_size(mx_in->d1, mx_in->d2, &mx_x_k)){
	return -1;
  }
  if(mx_set_size(mx_in->d1, mx_in->d2, &mx_P_x)){
	return -1;
  }
  if(mx_set_size(mx_in->d1, mx_in->d2, mx_P_n)){
	return -1;
  }


  //####################################
  // initialization (averaging)
  //####################################
  memcpy(mx_x_k.data, mx_in->data, mx_in->d2*sizeof(float)); //X_k(:,1)  = X_2(:,1);

  //##choose######
  //##beg_num_mean = floor((1/(1-th2))+0.5); //mean(1/(1-th2) - num of frames for beg est.
  beg_num_mean = floor((1/(1-num_fr_beg))+0.5);
  for(j = 0; j < mx_in->d2; j++) { //P_x(:,1)  = mean(X_2(:,1:round(1/(1-th2)))')';
    mx_P_x.data[j] = mx_in->data[j];
    for(i = 1; i < beg_num_mean; i++) { 
      mx_P_x.data[j] = (mx_P_x.data[j] + mx_in->data[i * mx_in->d2 + j]);
    }
    mx_P_x.data[j] = mx_P_x.data[j]/beg_num_mean;
  }

  //####################################
  // comp. of P_x and X_k (smoothing in time)
  //####################################
  for(j = 0; j < mx_in->d2; j++) { //I and II. recursive filter - smoothing in time for each band
	for(i = 1; i < mx_in->d1; i++) {
	  mx_x_k.data[i * mx_in->d2 + j] = (th1 * mx_x_k.data[(i-1) * mx_in->d2 + j]) + ((1-th1) * mx_in->data[i * mx_in->d2 + j]);
	  mx_P_x.data[i * mx_in->d2 + j] = (th2 * mx_P_x.data[(i-1) * mx_in->d2 + j]) + ((1-th2) * mx_in->data[i * mx_in->d2 + j]);
	}
  }

  //####################################
  //comp. of P_n - minima tracking in time Martin's rewriten method - using several small windows
  //####################################
  if(mx_set_size(D_sm, mx_in->d2, &minimum_val)){
	return -1;
  }
  for(j = 0; j < mx_in->d2; j++) { //initialization
	for(i = 0; i < D_sm; i++) {
	  minimum_val.data[i * mx_in->d2 + j] = mx_P_x.data[0 + j];
	}
	count = 0;
	for(i = 0; i < mx_in->d1; i++) {
	  count = count + 1;
	  if (count >= short_win) {
		for(k = 1; k < D_sm; k++) {
	      minimum_val.data[(k-1) * mx_in->d2 + j] = minimum_val.data[k * mx_in->d2 + j];
	    }
		minimum_val.data[(D_sm-1) * mx_in->d2 + j] = mx_P_x.data[i * mx_in->d2 + j];
		count = 0;
	  }
	  if(mx_P_x.data[i * mx_in->d2 + j] < minimum_val.data[(D_sm-1) * mx_in->d2 + j]){
		minimum_val.data[(D_sm-1) * mx_in->d2 + j] = mx_P_x.data[i * mx_in->d2 + j];
	  }

	  help_fl = minimum_val.data[0 * mx_in->d2 + j];  //find the minimum
	  for(k = 0; k < D_sm; k++) {
		help_fl = MIN(help_fl, minimum_val.data[k * mx_in->d2 + j]);
	  } 
	  mx_P_n->data[i * mx_in->d2 + j] = help_fl;
	}
  }
  mx_release(&minimum_val);

  //  mx_write_file(mx_P_n,"out.mx");

  //####################################
  //P_n = P_n*osub - oversubtraction - multiplication
  //####################################
  for(j = 0; j < mx_in->d2; j++) {
	for(i = 0; i < mx_in->d1; i++) {
	  mx_P_n->data[i * mx_in->d2 + j] =  mx_P_n->data[i * mx_in->d2 + j] * osub;
	}
  }
  // mx_write_file(mx_P_n,"out.mx");
 
  //########################################################
  //#######NOISE SUBTRACTION
  //########################################################
  if(mx_set_size(mx_in->d1, mx_in->d2, &mx_Q_k)){
    return -1;
  }
  //if(mx_set_size(mx_in->d1, mx_in->d2, mx_Y_k)){
  //  return -1;
  //}
  for(j = 0; j < mx_in->d2; j++) { 
    for(i = 0; i < mx_in->d1; i++) {
      mx_Q_k.data[i * mx_in->d2 + j] = 1 - sqrt(osub * mx_P_n->data[i * mx_in->d2 + j]/mx_x_k.data[i * mx_in->d2 + j]);
      if (mx_Q_k.data[i * mx_in->d2 + j]*sqrt(mx_in->data[i * mx_in->d2 + j]) <= (subf*sqrt(mx_P_n->data[i * mx_in->d2 + j]))){
		//        mx_Y_k->data[i * mx_in->d2 + j] = subf*sqrt(mx_P_n->data[i * mx_in->d2 + j]);
      }
      else{
        //mx_Y_k->data[i * mx_in->d2 + j] = mx_Q_k.data[i * mx_in->d2 + j]*sqrt(mx_in->data[i * mx_in->d2 + j]); //only magnitude spectrum
      }
    }
  }
  //compute back the power spectrum
  //for(j = 0; j < mx_in->d2; j++) { 
  //for(i = 0; i < mx_in->d1; i++) {
      //mx_Y_k->data[i * mx_in->d2 + j]  =  SQR(mx_Y_k->data[i * mx_in->d2 + j]);
  //}
  //}
  mx_release(&mx_x_k);
  mx_release(&mx_P_x);
  mx_release(&mx_Q_k);
  mx_release(&mx_mfb_norm);
  return NOISE_ESTIMATION_SUCCES;
}

