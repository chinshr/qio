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

opt  =
CC   = gcc -Wall -g
LIBS = -lm

MKDIR_P = mkdir -p
SRC_DIR = src
BIN_DIR = bin

all	: folders auroracalc_front auroracalc_back nr silence_flags drop clean

# Create folders
folders: ${BIN_DIR}

${BIN_DIR}:
	${MKDIR_P} ${BIN_DIR}

# Terminal Side
auroracalc_front: auroracalc_front.o fe_util.o dsp.o auroralib.o vadfeatures.o  noisecomp.o overlapadd.o klt.o mul_mdmd_md.o mul_mdmd_md_l0g.o mul_mdmd_md_l0nf.o mul_mfmf_mf_beta0.o mul_mfmft_mf_beta0.o
	$(CC) $(opt) -o $(BIN_DIR)/auroracalc_front $(BIN_DIR)/auroracalc_front.o $(BIN_DIR)/fe_util.o $(BIN_DIR)/dsp.o $(BIN_DIR)/auroralib.o $(BIN_DIR)/vadfeatures.o $(BIN_DIR)/noisecomp.o $(BIN_DIR)/overlapadd.o $(BIN_DIR)/klt.o $(BIN_DIR)/mul_mdmd_md.o $(BIN_DIR)/mul_mdmd_md_l0g.o $(BIN_DIR)/mul_mdmd_md_l0nf.o $(BIN_DIR)/mul_mfmf_mf_beta0.o $(BIN_DIR)/mul_mfmft_mf_beta0.o -I./$(SRC_DIR)/libogi -I./$(SRC_DIR)/libtcts -I./$(SRC_DIR)/libicsi  $(LIBS)

auroracalc_front.o: $(SRC_DIR)/auroracalc_front.c
	$(CC) $(opt) -c $(SRC_DIR)/auroracalc_front.c -o $(BIN_DIR)/$@ -I./$(SRC_DIR) -I./$(SRC_DIR)/libogi -I./$(SRC_DIR)/libicsi -I./$(SRC_DIR)/libtcts

# Server Side
auroracalc_back: auroracalc_back.o fe_util.o dsp.o auroralib.o vadfeatures.o  noisecomp.o overlapadd.o klt.o mul_mdmd_md.o mul_mdmd_md_l0g.o mul_mdmd_md_l0nf.o mul_mfmf_mf_beta0.o mul_mfmft_mf_beta0.o
	$(CC) $(opt) -o $(BIN_DIR)/auroracalc_back $(BIN_DIR)/auroracalc_back.o $(BIN_DIR)/fe_util.o $(BIN_DIR)/dsp.o $(BIN_DIR)/auroralib.o $(BIN_DIR)/vadfeatures.o $(BIN_DIR)/noisecomp.o $(BIN_DIR)/overlapadd.o $(BIN_DIR)/klt.o $(BIN_DIR)/mul_mdmd_md.o $(BIN_DIR)/mul_mdmd_md_l0g.o $(BIN_DIR)/mul_mdmd_md_l0nf.o $(BIN_DIR)/mul_mfmf_mf_beta0.o $(BIN_DIR)/mul_mfmft_mf_beta0.o -I./$(SRC_DIR)/libogi -I./$(SRC_DIR)/libtcts -I./$(SRC_DIR)/libicsi  $(LIBS)

auroracalc_back.o : $(SRC_DIR)/auroracalc_back.c
	$(CC) $(opt) -c $(SRC_DIR)/auroracalc_back.c -o $(BIN_DIR)/$@ -I./$(SRC_DIR)/ -I./$(SRC_DIR)/libogi -I./$(SRC_DIR)/libicsi -I./$(SRC_DIR)/libtcts

# Standalone VAD
silence_flags: silence_flags.o fe_util.o dsp.o auroralib.o vadfeatures.o   noisecomp.o overlapadd.o klt.o mul_mdmd_md.o mul_mdmd_md_l0g.o mul_mdmd_md_l0nf.o mul_mfmf_mf_beta0.o mul_mfmft_mf_beta0.o
	$(CC) $(opt) -o $(BIN_DIR)/silence_flags $(BIN_DIR)/silence_flags.o $(BIN_DIR)/fe_util.o $(BIN_DIR)/dsp.o $(BIN_DIR)/auroralib.o $(BIN_DIR)/vadfeatures.o $(BIN_DIR)/noisecomp.o $(BIN_DIR)/overlapadd.o $(BIN_DIR)/klt.o $(BIN_DIR)/mul_mdmd_md.o $(BIN_DIR)/mul_mdmd_md_l0g.o $(BIN_DIR)/mul_mdmd_md_l0nf.o $(BIN_DIR)/mul_mfmf_mf_beta0.o $(BIN_DIR)/mul_mfmft_mf_beta0.o -I./$(SRC_DIR)/libogi -I./$(SRC_DIR)/libtcts -I./$(SRC_DIR)/libicsi -lm

silence_flags.o : $(SRC_DIR)/silence_flags.c
	$(CC) $(opt) -c $(SRC_DIR)/silence_flags.c -o $(BIN_DIR)/$@ -I./$(SRC_DIR) -I./$(SRC_DIR)/libogi -I./$(SRC_DIR)/libicsi -I./$(SRC_DIR)/libtcts

# Standalone frame dropping
drop: drop.o fe_util.o mul_mfmft_mf_beta0.o
	$(CC) $(opt) -o $(BIN_DIR)/drop $(BIN_DIR)/drop.o $(BIN_DIR)/fe_util.o $(BIN_DIR)/mul_mfmft_mf_beta0.o -lm

drop.o : $(SRC_DIR)/drop.c
	$(CC) $(opt) -c $(SRC_DIR)/drop.c -o $(BIN_DIR)/$@ -I./$(SRC_DIR) -I./$(SRC_DIR)/libogi 

# Standalone noise reduction
nr: nr.o fe_util.o dsp.o auroralib.o vadfeatures.o   noisecomp.o overlapadd.o klt.o mul_mdmd_md.o mul_mdmd_md_l0g.o mul_mdmd_md_l0nf.o mul_mfmf_mf_beta0.o mul_mfmft_mf_beta0.o
	$(CC) $(opt) -o $(BIN_DIR)/nr $(BIN_DIR)/nr.o $(BIN_DIR)/fe_util.o $(BIN_DIR)/dsp.o $(BIN_DIR)/auroralib.o $(BIN_DIR)/vadfeatures.o $(BIN_DIR)/noisecomp.o $(BIN_DIR)/overlapadd.o $(BIN_DIR)/klt.o $(BIN_DIR)/mul_mdmd_md.o $(BIN_DIR)/mul_mdmd_md_l0g.o $(BIN_DIR)/mul_mdmd_md_l0nf.o $(BIN_DIR)/mul_mfmf_mf_beta0.o $(BIN_DIR)/mul_mfmft_mf_beta0.o -I./$(SRC_DIR)/libogi -I./$(SRC_DIR)/libtcts -I./$(SRC_DIR)/libicsi  -lm

nr.o: $(SRC_DIR)/nr.c
	$(CC) $(opt) -c $(SRC_DIR)/nr.c -o $(BIN_DIR)/$@ -I./$(SRC_DIR)/ -I./$(SRC_DIR)/libogi -I./$(SRC_DIR)/libicsi -I./$(SRC_DIR)/libtcts 

# Libraries
fe_util.o  : ./$(SRC_DIR)/libogi/fe_util.c ./$(SRC_DIR)/libogi/fe_util.h
	$(CC) $(opt) -c ./$(SRC_DIR)/libogi/fe_util.c -o ./$(BIN_DIR)/$@ -I./$(SRC_DIR)/libogi -I.

auroralib.o  : ./$(SRC_DIR)/libicsi/auroralib.c ./$(SRC_DIR)/libicsi/auroralib.h
	$(CC) $(opt) -c ./$(SRC_DIR)/libicsi/auroralib.c -o ./$(BIN_DIR)/$@ -I./$(SRC_DIR) -I./$(SRC_DIR)/libicsi -I./$(SRC_DIR)/libtcts

dsp.o : ./$(SRC_DIR)/libtcts/dsp.c ./$(SRC_DIR)/libtcts/dsp.h
	$(CC) $(opt) -c ./$(SRC_DIR)/libtcts/dsp.c -o ./$(BIN_DIR)/$@ -I./$(SRC_DIR)

noisecomp.o : ./$(SRC_DIR)/libicsi/noisecomp.c ./$(SRC_DIR)/libicsi/noisecomp.h
	$(CC) $(opt) -c ./$(SRC_DIR)/libicsi/noisecomp.c -o ./$(BIN_DIR)/$@ -I./$(SRC_DIR) -I./$(SRC_DIR)/libicsi

vadfeatures.o : ./$(SRC_DIR)/libicsi/vadfeatures.c ./$(SRC_DIR)/libicsi/vadfeatures.h
	$(CC) $(opt) -c ./$(SRC_DIR)/libicsi/vadfeatures.c -o ./$(BIN_DIR)/$@ -I./$(SRC_DIR) -I./$(SRC_DIR)/libicsi

overlapadd.o : ./$(SRC_DIR)/libicsi/overlapadd.c ./$(SRC_DIR)/libicsi/overlapadd.h
	$(CC) $(opt) -c ./$(SRC_DIR)/libicsi/overlapadd.c -o ./$(BIN_DIR)/$@ -I./$(SRC_DIR) -I./$(SRC_DIR)/libicsi

klt.o : ./$(SRC_DIR)/libicsi/klt.c ./$(SRC_DIR)/libicsi/klt.h mul_mdmd_md.o mul_mdmd_md_l0g.o mul_mdmd_md_l0nf.o
	$(CC) $(opt) -c ./$(SRC_DIR)/libicsi/klt.c -o ./$(BIN_DIR)/$@ -I./$(SRC_DIR) -I./$(SRC_DIR)/libicsi

mul_mfmf_mf_beta0.o : ./$(SRC_DIR)/libicsi/mul_mfmf_mf_beta0.c
	$(CC) $(opt) -c ./$(SRC_DIR)/libicsi/mul_mfmf_mf_beta0.c -o ./$(BIN_DIR)/$@ -I./$(SRC_DIR) -I./$(SRC_DIR)libicsi

mul_mfmft_mf_beta0.o : ./$(SRC_DIR)/libicsi/mul_mfmft_mf_beta0.c
	$(CC) $(opt) -c ./$(SRC_DIR)/libicsi/mul_mfmft_mf_beta0.c -o ./$(BIN_DIR)/$@ -I./$(SRC_DIR) -I./$(SRC_DIR)/libicsi

mul_mdmd_md.o : ./$(SRC_DIR)/libicsi/mul_mdmd_md.c
	$(CC) $(opt) -c ./$(SRC_DIR)/libicsi/mul_mdmd_md.c -o ./$(BIN_DIR)/$@ -I./$(SRC_DIR) -I./$(SRC_DIR)/libicsi

mul_mdmd_md_l0g.o : ./$(SRC_DIR)/libicsi/mul_mdmd_md_l0g.c
	$(CC) $(opt) -c ./$(SRC_DIR)/libicsi/mul_mdmd_md_l0g.c -o ./$(BIN_DIR)/$@ -I./$(SRC_DIR) -I./$(SRC_DIR)/libicsi

mul_mdmd_md_l0nf.o : ./$(SRC_DIR)/libicsi/mul_mdmd_md_l0nf.c
	$(CC) $(opt) -c ./$(SRC_DIR)/libicsi/mul_mdmd_md_l0nf.c -o ./$(BIN_DIR)/$@ -I./$(SRC_DIR) -I./$(SRC_DIR)/libicsi

clean 	:
	rm $(BIN_DIR)/*.o 
