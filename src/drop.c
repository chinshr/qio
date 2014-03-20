/* THIS SOFTWARE IS PROVIDED BY THE AUTHORS AND CONTRIBUTORS ``AS IS'' AND    *
 * ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE      *
 * IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE *
 * ARE DISCLAIMED.  IN NO EVENT SHALL THE AUTHORS OR CONTRIBUTORS BE LIABLE   *  
 * FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR CONSEQUENTIAL *
 * DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF SUBSTITUTE GOODS    *
 * OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS INTERRUPTION)      * 
 * HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN CONTRACT, STRICT * 
 * LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY  *  
 * OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF     *   
 * SUCH DAMAGE.                                                               *
 */

#include <stdio.h>
#include <stdlib.h>
#include "libogi/fe_util.h"
#include "libicsi/auroralib.h"

muT mx_in;
muT mx_out;

/* 

   This program reads in a feature file and a VAD flag file.
   It drops frames according to the VAD flag file and outputs
   a frame dropped feature file.  In the VAD flag file, '1' = nonspeech
   and '0' = speech.

   Setting inSwapFlag and outSwapFlag to 1 to enable of byte-swapping
   all data on input and output may be needed if using a little-endian
   machine.

*/

int main(int argc, char *argv[])
{

  FILE *silfile;
  char *silfilename;
  char *outfilename;
  char *filename;
  char *silflag;
  x_int32 sampPeriod;
  x_int16 parmKind;
  char inSwapFlag, outSwapFlag;
  
  int c;
  int n = 0; // number of frames
  int f = 0; // frame counter

  if (argc != 6) {
    printf("Usage: drop silfile feafile outfile inSwapFlag outSwapFlag\n");
    exit(1);
  }

  inSwapFlag = argv[4][0];
  outSwapFlag = argv[5][0];

  if (inSwapFlag == '0')
    inSwapFlag = 0;
  else if (inSwapFlag == '1')
    inSwapFlag = 1;
  else {
    printf("inSwapFlag must be either 1 (byte swap) or 0 (no byte swap)\n");
    exit(1);
  }

  if (outSwapFlag == '0')
    outSwapFlag = 0;
  else if (outSwapFlag == '1')
    outSwapFlag = 1;
  else {
    printf("outSwapFlag must be either 1 (byte swap) or 0 (no byte swap)\n");
    exit(1);
  }
  

  /**********************************************************
    Read silence flags ('0' or '1' for each frame) from file.
   **********************************************************/
  
  silfilename = argv[1];  
  silfile = fopen(silfilename, "r");
  if (silfile == NULL) {
    fprintf (stderr, "ERROR:   Can not open file '%s' !\r\n", silfilename);
    exit(1);
  }

  while (((c = fgetc(silfile)) != EOF) && (c != '\n')) {
    n++;
  }

  silflag = (char *) malloc( n * sizeof(char) );
  if (silflag == NULL) {
    fprintf (stderr, "ERROR:   can't alloc %d bytes\r\n", n);
    exit(1);
  }

  rewind(silfile);
  for (f = 0; f < n; f++) {
    c = fgetc(silfile);
    if ((c != '0') && (c != '1')) {
      fprintf (stderr,
	       "ERROR:   bad char '%c' (%d) in silence flag file\r\n", c, c);
      exit(1);
    }
    silflag[f] = (char) (c - '0');
  }
  fclose(silfile);

  /**********************************************************
    Read feature file
   **********************************************************/

  filename = argv[2];
  switch(mx_read_htk_file(filename, &mx_in, &sampPeriod, &parmKind,inSwapFlag)) {

  case HTK_FILE_CANT_OPEN:
    fprintf(stderr, "Cannot open input file: %s\n", filename);
    return-1;
  case HTK_FILE_IO_ERROR:
    fprintf(stderr, "Error reading input file: %s\n", filename);
    return-1;
    //case HTK_FILE_SUCCESS:
  }

  // Make sure HTK CRC checksum flag is off, since we won't
  // write a CRC checksum to the output HTK file.
  parmKind &= ~010000;  

  /**********************************************************
    Drop silence frames
   **********************************************************/

  if (n != mx_in.d1){
    fprintf(stderr,"%d frames in silence file but %d in feature file.\n", n, mx_in.d1);
    return -1;
  }
    
  fprintf(stderr,"Features: %d x %d -> ", mx_in.d1,mx_in.d2);
  framedropping(&mx_in,silflag,0,&mx_out);
  mx_release(&mx_in); 
  mx_in = mx_out;

  /**********************************************************
    Write output features
   **********************************************************/

  outfilename = argv[3];
  fprintf(stderr,"Features: %d x %d\n", mx_in.d1,mx_in.d2);
  switch(mx_write_htk_file(outfilename, &mx_in, sampPeriod, parmKind, outSwapFlag)) {
  case HTK_FILE_CANT_OPEN:
    fprintf(stderr, "Cannot open output feature file: %s\n", outfilename);
    break;
  case HTK_FILE_IO_ERROR:
    fprintf(stderr, "Cannot write to output feature file: %s\n", outfilename);
    break;
  case HTK_FILE_INVALID_HEADER:
    fprintf(stderr, "Fatal error: Trying write nonsense header to output\n");
    return -1;
    //case HTK_FILE_SUCCESS:
  }


 return 0;

}
