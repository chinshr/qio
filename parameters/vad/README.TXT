Default files
-------------

The following are the default MLP weight and normalization files used
by the VAD in the QIO front end:

net.tim-fin-tic-it-spn-rand.54i+50h+2o.0-delay-wiener+dct+lpf.wts.head
tim-fin-tic-it-spn-rand.0-delay-wiener+dct+lpf.norms

Alternate files: no Wiener filter prior to VAD
-----------------------------------------------

The following are alternate weight and normalization files that may be
useful when running VAD standalone using the silence_flags tool:

net.tim-fin-tic-spn-rand.54i+50h+2o.mel-delay+dct+lpf.wts.head
tim-fin-tic-spn-rand.mel-delay+dct+lpf.norms

These alternate files were created without running the Wiener filter 
before the VAD.  The point of this is that on non-Aurora data sets, it
may not be valid to initialize the Wiener filter's noise estimate
using the assumption that files start with a period of nonspeech.
(Also, for these alternate files, the SpeechDatCar-Italian data was
omitted from the VAD training set due to its license terms.)

See ../../../scripts/silenceflag_script for examples of how to specify
the weights and norm file used by the VAD using the -VADweights and
-VADnorm command line options and how to disable the Wiener filter
using the -S option.


Alternate files: no Wiener filter prior to VAD; 20 ms window length
-------------------------------------------------------------------

The following are alternate weight and normalization files that may be
useful when running VAD standalone using silence_flags:

net.tim-fin-tic-spn-rand.54i+50h+2o.win20-mel-delay+dct+lpf.wts.head
tim-fin-tic-spn-rand.win20-mel-delay+dct+lpf.norms

These were created without running the Wiener filter before the VAD,
and using a 20 ms DFT window length (matching the nr tool) rather than
25 ms (matching the default DFT window length of the other tools
including the silence_flags tool).  Users may prefer to use these
files when calculating silence flags which will be used with the nr tool.

See ../../../scripts/silenceflag_script for examples of how to specify
the weights and norm file used by the VAD using the -VADweights and
-VADnorm command line options, how to disable the Wiener filter
running before the VAD using the -S option, and how to change the DFT
length used by silence_flags to 20 ms.

Notes
-----

You may notice that the default filenames contain "0-delay" while the
alternate filenames instead contain "delay".  Sadly, the meaning (if
any) of the "0-delay" vs. "delay" difference has been forgotten.  If
this information is recovered it should be added to this document in
the future.
