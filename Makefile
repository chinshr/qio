# THIS SOFTWARE IS PROVIDED BY THE AUTHORS AND CONTRIBUTORS ``AS IS'' AND   
# ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE     
# IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
# ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE 
# FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL 
# DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS   
# OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)
# CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT      
# LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY 
# OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF    
# SUCH DAMAGE.                                                            

opt =
CC  = gcc -Wall -g
LIBS = -lm

all : auroracalc_front auroracalc_back nr silence_flags drop clean

##Terminal Side
auroracalc_front   : auroracalc_front.o fe_util.o dsp.o auroralib.o vadfeatures.o  noisecomp.o overlapadd.o klt.o mul_mdmd_md.o mul_mdmd_md_l0g.o mul_mdmd_md_l0nf.o mul_mfmf_mf_beta0.o mul_mfmft_mf_beta0.o
  $(CC) $(opt) -o auroracalc_front auroracalc_front.o fe_util.o dsp.o auroralib.o vadfeatures.o  noisecomp.o overlapadd.o klt.o mul_mdmd_md.o mul_mdmd_md_l0g.o mul_mdmd_md_l0nf.o mul_mfmf_mf_beta0.o mul_mfmft_mf_beta0.o -I./libogi -I./libtcts -I./libicsi  $(LIBS)

auroracalc_front.o : auroracalc_front.c
  $(CC) $(opt) -c auroracalc_front.c -I./ -I./libogi -I./libicsi -I./libtcts

##Server Side
auroracalc_back   : auroracalc_back.o fe_util.o dsp.o auroralib.o vadfeatures.o  noisecomp.o overlapadd.o klt.o mul_mdmd_md.o mul_mdmd_md_l0g.o mul_mdmd_md_l0nf.o mul_mfmf_mf_beta0.o mul_mfmft_mf_beta0.o
  $(CC) $(opt) -o auroracalc_back auroracalc_back.o fe_util.o dsp.o auroralib.o vadfeatures.o   noisecomp.o overlapadd.o klt.o mul_mdmd_md.o mul_mdmd_md_l0g.o mul_mdmd_md_l0nf.o mul_mfmf_mf_beta0.o mul_mfmft_mf_beta0.o -I./libogi -I./libtcts -I./libicsi  $(LIBS)

auroracalc_back.o : auroracalc_back.c
  $(CC) $(opt) -c auroracalc_back.c -I./ -I./libogi -I./libicsi -I./libtcts

##Standalone VAD
silence_flags: silence_flags.o fe_util.o dsp.o auroralib.o vadfeatures.o   noisecomp.o overlapadd.o klt.o mul_mdmd_md.o mul_mdmd_md_l0g.o mul_mdmd_md_l0nf.o mul_mfmf_mf_beta0.o mul_mfmft_mf_beta0.o
  $(CC) $(opt) -o silence_flags silence_flags.o fe_util.o dsp.o auroralib.o vadfeatures.o   noisecomp.o overlapadd.o klt.o mul_mdmd_md.o mul_mdmd_md_l0g.o mul_mdmd_md_l0nf.o mul_mfmf_mf_beta0.o mul_mfmft_mf_beta0.o -I./libogi -I./libtcts -I./libicsi -lm

silence_flags.o : silence_flags.c
  $(CC) $(opt) -c silence_flags.c -I./ -I./libogi -I./libicsi -I./libtcts

##Standalone frame dropping
drop: drop.o fe_util.o mul_mfmft_mf_beta0.o
  $(CC) $(opt) -o drop drop.o fe_util.o mul_mfmft_mf_beta0.o -lm

drop.o : drop.c
  $(CC) $(opt) -c drop.c -I./ -I./libogi 

##Standalone noise reduction
nr   : nr.o fe_util.o dsp.o auroralib.o vadfeatures.o   noisecomp.o overlapadd.o klt.o mul_mdmd_md.o mul_mdmd_md_l0g.o mul_mdmd_md_l0nf.o mul_mfmf_mf_beta0.o mul_mfmft_mf_beta0.o
  $(CC) $(opt) -o nr nr.o fe_util.o dsp.o auroralib.o vadfeatures.o   noisecomp.o overlapadd.o klt.o mul_mdmd_md.o mul_mdmd_md_l0g.o mul_mdmd_md_l0nf.o mul_mfmf_mf_beta0.o mul_mfmft_mf_beta0.o -I./libogi -I./libtcts -I./libicsi  -lm

nr.o :nr.c
  $(CC) $(opt) -c nr.c -I./ -I./libogi -I./libicsi -I./libtcts 



##Libraries
fe_util.o  : ./libogi/fe_util.c ./libogi/fe_util.h
  $(CC) $(opt) -c ./libogi/fe_util.c -I./libogi -I.

auroralib.o  : ./libicsi/auroralib.c ./libicsi/auroralib.h
  $(CC) $(opt) -c ./libicsi/auroralib.c -I. -I./libicsi -I./libtcts

dsp.o : ./libtcts/dsp.c ./libtcts/dsp.h
  $(CC) $(opt) -c ./libtcts/dsp.c -I./



noisecomp.o : ./libicsi/noisecomp.c ./libicsi/noisecomp.h
  $(CC) $(opt) -c ./libicsi/noisecomp.c -I. -I./libicsi

vadfeatures.o : ./libicsi/vadfeatures.c ./libicsi/vadfeatures.h
  $(CC) $(opt) -c ./libicsi/vadfeatures.c -I. -I./libicsi

overlapadd.o : ./libicsi/overlapadd.c ./libicsi/overlapadd.h
  $(CC) $(opt) -c ./libicsi/overlapadd.c -I. -I./libicsi

klt.o : ./libicsi/klt.c ./libicsi/klt.h mul_mdmd_md.o mul_mdmd_md_l0g.o mul_mdmd_md_l0nf.o
  $(CC) $(opt) -c ./libicsi/klt.c -I. -I./libicsi

mul_mfmf_mf_beta0.o : ./libicsi/mul_mfmf_mf_beta0.c
  $(CC) $(opt) -c ./libicsi/mul_mfmf_mf_beta0.c -I. -I./libicsi

mul_mfmft_mf_beta0.o : ./libicsi/mul_mfmft_mf_beta0.c
  $(CC) $(opt) -c ./libicsi/mul_mfmft_mf_beta0.c -I. -I./libicsi

mul_mdmd_md.o : ./libicsi/mul_mdmd_md.c
  $(CC) $(opt) -c ./libicsi/mul_mdmd_md.c -I. -I./libicsi

mul_mdmd_md_l0g.o : ./libicsi/mul_mdmd_md_l0g.c
  $(CC) $(opt) -c ./libicsi/mul_mdmd_md_l0g.c -I. -I./libicsi

mul_mdmd_md_l0nf.o : ./libicsi/mul_mdmd_md_l0nf.c
  $(CC) $(opt) -c ./libicsi/mul_mdmd_md_l0nf.c -I. -I./libicsi

clean   :
  rm *.o 
