The Qualcomm-ICSI-OGI (QIO) Aurora front end
============================================

The qio/ directory contains the code and parameter files of the QIO
front end.  This is the front end described in the ICSLP 2002 paper by
Adami et al., which is available on the ICSI Speech Group publications
page.  The full front end produces 51 features including 6 TRAPS-based
features.  The version of the front end submitted to the ETSI
committee in January 2002 was different: it did not have the 6
TRAPS-based features and so only produced 45 features.  That version
is referred to as QIO-NoTRAPS in the ICSLP 2002 paper.  It is described
in depth in the document QIO-NoTRAPS.pdf (Adobe PDF format), which is
a copy of the documentation submitted to the ETSI Aurora committee.  The
document NoiseReduction.pdf goes into additional detail about the noise 
reduction component.  To compile the software, cd to the qio/src
directory and run 'make' or 'gmake'.

The QIO front end consists of the auroracalc_front and auroracalc_back
programs.  These represent terminal-side and server-side processing in
a distributed speech recognition (DSR) environment.  A third tool,
coder, was run between auroracalc_front and auroracalc_back in order
to simulate compression and quantization for transmission.  This is
the historical reason for the split of the front end into
auroracalc_front and auroracalc_back.  For copyright reasons, the
coder tool is not included in this archive.  The recognition accuracy
results in the ICSLP 2002 paper and the QIO-NoTRAPS.pdf document may
be slightly affected by the quantization performed by the coder tool.

The auroracalc_front and auroracalc_back tools have many command line
options. Running the tools with no command line arguments will produce
an options list.  Unfortunately, some of these options are poorly
tested, obsolete or ineffectual; this is a legacy of the fact that
many options were created for experimental purposes during the
development of these tools.  Also, some options included in the lists
for both the auroracalc_front and auroracalc_back tools are in fact
only applicable to one of those tools.  Please refer to the files in the
scripts/ directory for guidance about normal options for running these
tools.  To obtain only the 45 QIO-NoTRAPS features, suppress the
TRAPS-based features by running the auroracalc_back binary as
"auroracalc_back -TRAPS 0 ...".  To obtain all 51 features, enable the
TRAPS-based features by running it as "auroracalc_back -TRAPS 1 ...".
A KLT is applied to the length-51 feature vector before it is passed
to the HMM back end; this step was accidentally omitted from the
system description in the ICSLP 2002 paper.

The ETSI Aurora benchmarking rules required the front end to process
each user utterance independently.  Thus the front end does not keep
any information from one user utterance to the next.  The ETSI Aurora
rules also required that the front end have a low algorithmic latency,
or in other words, that the algorithms used in the front end do not
require more than a small delay between the arrival of input samples
and the generation of the corresponding output frames.  This has
affected the design of some features such as the mean and variance
normalization.  See QIO-NoTRAPS.pdf for an analysis of the front end's
algorithmic latency.

Standalone VAD (Voice Activity Detection) and noise reduction
=============================================================

It is sometimes useful to use the VAD and noise reduction components
of the front end without using the entire front end.  There are
programs for this purpose under the qio/ directory:
silence_flags (silence_flags.c), nr (nr.c), and drop (drop.c).

silence_flags and nr
--------------------

The program silence_flags performs VAD and outputs an ASCII file of
'1' or '0' "silence" flags (one for every frame in the input audio
file; frame step is 10 ms) which label frames as nonspeech or speech
respectively (so '1' stands for a frame judged to be nonspeech and '0'
stands for a frame judged to be speech).  (The name "silence" is a
misnomer for nonspeech since there may be noise even when there is no
speech.)

The nr program performs the Wiener filter noise reduction on an input
waveform and outputs a noise-reduced waveform.  This depends on a
noise estimate made over frames judged to be nonspeech.  In the QIO
front end the noise estimate is initialized from the beginning of each
utterance, under the assumption that each utterance starts with a
period of nonspeech, and updated using later frames of the utterance
judged to be nonspeech based on an energy threshold. (See section
2.7.1 of QIO-NoTRAPS.pdf.)  Noise estimation is done in a causal
fashion, so that each frame is noise-reduced using a noise estimate
which does not depend on future frames.  When the "-Ssilfile
<flags_filename>" command line option, this program performs noise
estimation in a different way: it reads in a VAD flags file (as
created by silence_flags.c) for each utterance and then noise-reduces
the utterance using a single noise estimate calculated over all the
frames labeled nonspeech by the silence flag file.  This alternate
style of noise estimation may be preferable if not every
utterance starts with a sufficiently long period of nonspeech (see
section 2.7.1 of QIO-NoTRAPS.pdf for information on what "sufficiently
long" means).

The Wiener filter imposes a spectral floor, but this floor is low, so
it can be partly lost in 16-bit integer output.  Past experience shows
that the use of 16-bit output from the Wiener filter can lower ASR
performance, perhaps because of the use of the logarithm function
(which is volatile near 0) in ASR feature extraction.  To output the
resynthesized audio after noise reduction in (headerless) floating
point, instead, use the command line option "-Ssynthfloat 1".  In
tests on the Aurora 2 task, using an MFCC front end with cepstral mean
subtraction and Asela Gunawardana's "complex" back-end configuration,
using 16-bit output from nr resulted in 72.9% average word accuracy
and using floating point output resulted in 76.7% average word
accuracy.

Sometimes the output file from nr will be slightly shorter
than the input file.  The difference should always be less than
the length of the frame step (10 ms).

Memory usage for nr can be a problem with large input files.  The
tool does each stage of procesing for every frame in the input file
before moving to the next stage of processing.  Peak memory usage
occurs when the input magnitude spectrogram, input phase spectrogram,
noise estimate magnitude spectrogram, smoothed noise estimate
magnitude spectrogram, and noise-reduced magnitude spectrogram are all
held in memory at once.  The memory usage for each spectrogram is
approximately

  (number of seconds of input data) * (100 frames/second) * 
     (NFFT/2 floats/frame) * (4 bytes/float) 

where NFFT = 256 at 8000 Hz sampling rate and 512 at the 11000 and
16000 Hz sampling rates.

For example, a 10-minute file at 16000 Hz would result in an
nr peak memory usage of approximately

  (600 seconds) * (100 frames/second) * (256 floats/frame) * (4
     bytes/float) * (5 spectrograms) =~ 300 megabytes

See scripts/silenceflag_script and scripts/noisered_script for
examples of the use of these tools, including the use of command line
options to set the sampling rate and big/little endianness of the
input data.  (Confusingly, there are some command line options for the
silence_flags or nr tool which are just inherited from the
auroracalc_front and auroracalc_back source code.  Not all of these
options have been tested with silence_flags and nr; some of them will
may actually have no effect when used with silence_flag or nr.)

silence_flags and nr were used for the ICSI-SRI-UW entry in the 2004
NIST Meeting Transcription evaluation, described in the ICSLP 2004
paper by Mirghafori et al.  The scripts scripts/silenceflag_script and
scripts/noisered_script are related to the ones used for that
evaluation.  Note that the VAD used by silence_flags has not been
optimized for meeting data (the training data used for the VAD MLP was from
the Aurora 2 corpus, which is TIDIGITS with noise and convolutional
distortion added, and from SpeechDatCar) and that it is
possible to use your own VAD techniques to create flag files for nr.
Florian Metze of ISL improved his nr performance on meeting data by
using nr with his own VAD which simply computed the power for all
frames of the recording, smoothed the power curve by low-pass
filtering across frames, and then marked the 5% or so frames with the
lowest power as silence.


drop
----

drop is a simple program which can be used with the output of
silence_flags to apply frame dropping (dropping of nonspeech frames)
to HTK feature files.  It takes as input an HTK feature file and an
ASCII file of flags as produced by silence_flags, and produces a new
HTK feature file with the frames labeled as nonspeech removed. (The
ICSI tool feacat, found in ICSI's SPRACHcore software package, can be
used to convert between feature file formats such as HTK, SRI, and
ICSI pfile.)

silence_flags and drop were used to implement frame-dropping for the
ICSLP 2002 paper by Kleinschmidt and Gelbart.

AURORACALC environment variable
===============================

Before using the tools (except for drop which doesn't require this)
the AURORACALC environment variable must be set to the full path of
the installation directory aurora-front-end/qio.  This is so the tools
can load parameter files (like multi-layer perceptron weights for VAD
and TRAPS). Running the tools with no arguments will produce a list of
command line options.

If there is a space inside a directory name in AURORACALC (such as in
the name of the "Documents and Settings" directory under Windows XP),
that could cause problems with scripts/auroracalc_script and perhaps
with other code as well.

VAD training
============

The vad-training/ directory contains code, script, and label files
used for training the VAD (voice activity detection) MLP (multi-layer
perceptron) using ICSI's Quicknet MLP tools.  These are only needed if
you wish to re-train the VAD.

Technical support for users
===========================

This package has been tested under Red Hat Enterprise Linux 3.0 and
4.0 using an x86 CPU.  Under Windows, the best option is probably to
compile the software using the free Cygwin tools.  This package has
been tested under Windows XP using Cygwin.  We hope this package is portable to
other platforms, but some modification might be necessary.

The file qio/src/libogi/fe_util.h defines the x_int32 and x_int16
types which are used in the definition of HTK feature file header
fields.  If it is not true on your system that sizeof(int) == 4 and
sizeof(short) == 2, you must edit the definitions of x_int32 and
x_int16.  Note that we have not tested this code on systems where it
is not true that sizeof(int) == 4 and sizeof(short) == 2, so be sure
to test after editing (see below for information about testing).

Users' questions about this software can be directed to the ICSI
Speech Group's tools forum at
http://groups.yahoo.com/group/icsi-speech-tools.  There is no
guarantee, however, that questions will be answered.

Users should also be aware that these tools were originally created
for internal use by the developers themselves, and as a result the
user documentation consists only of some README files and the command line
option information printed when the tools are run with no arguments.
Users may find it necessary to refer to provided design documentation and
source code in order to fully understand the use of these tools.

Installation verification
=========================

The test/ directory contains sample input and output files.  
The output files can be compared against after installation
as a check that installation has proceeded correctly. 

