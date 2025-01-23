#include "csdr/ringbuffer.hpp"
#include "csdr/cw.hpp"
#include "fftw3.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

#define USE_NEIGHBORS  1 // 1: Subtract neighbors from each FFT bucket
#define USE_AVG_BOTTOM 1 // 1: Subtract average value from each bucket
#define USE_AVG_RATIO  0 // 1: Divide each bucket by average value

#define MAX_CHANNELS (sampleRate/2/100)
#define MAX_INPUT    (MAX_CHANNELS*2)
#define INPUT_STEP   (MAX_INPUT)//MAX_INPUT/4)
#define AVG_SECONDS  (3)
#define NEIGH_WEIGHT (0.5)

unsigned int sampleRate = 48000; // Input audio sampling rate
unsigned int printChars = 8;     // Number of characters to print at once
bool use16bit = false;           // TRUE: Use S16 input values (else F32)
bool showCw   = false;           // TRUE: Show dits and dahs
bool showDbg  = false;           // TRUE: Print debug data to stderr

Csdr::Ringbuffer<float> **in;
Csdr::RingbufferReader<float> **inReader;
Csdr::Ringbuffer<unsigned char> **out;
Csdr::RingbufferReader<unsigned char> **outReader;
Csdr::CwDecoder<float> **cwDecoder;

// Print output from ith decoder
void printOutput(FILE *outFile, int i, unsigned int freq, unsigned int printChars)
{
  // Must have a minimum of printChars
  int n = outReader[i]->available();
  if(n<printChars) return;

  // Print frequency
  fprintf(outFile, "%d:", freq);

  // Print characters
  unsigned char *p = outReader[i]->getReadPointer();
  for(int j=0 ; j<n ; ++j)
  {
    if((j<n-1) && strchr("TEI ", p[j+1]) && strchr("TEI ", p[j])) ++j;
    else fprintf(outFile, "%c", p[j]);
  }

  // Done printing
  outReader[i]->advance(n);
  printf("\n");
  fflush(outFile);
}

int main(int argc, char *argv[])
{
  FILE *inFile, *outFile;
  const char *inName, *outName;
  float accPower, avgPower, maxPower;
  int j, i, k, n, remains;

  // Parse input arguments
  for(j=1, inName=outName=0, inFile=stdin, outFile=stdout ; j<argc ; ++j)
  {
    if(argv[j][0]!='-')
    {
      // First two non-option arguments are filenames
      if(!inName) inName = argv[j];
      else if(!outName) outName = argv[j];
      else
      {
        fprintf(stderr, "%s: Excessive file name '%s'!\n", argv[0], argv[j]);
        return(2);
      }
    }
    else if(strlen(argv[j])!=2)
    {
      // Single-letter options only!
      fprintf(stderr, "%s: Unrecognized option '%s'!\n", argv[0], argv[j]);
      return(2);
    }
    else switch(argv[j][1])
    {
      case 'n':
        printChars = j<argc-1? atoi(argv[++j]) : printChars;
        printChars = printChars<1? 1 : printChars>32? 32 : printChars;
        break;
      case 'r':
        sampleRate = j<argc-1? atoi(argv[++j]) : sampleRate;
        sampleRate = sampleRate<8000? 8000 : sampleRate>48000? 48000 : sampleRate;
        break;
      case 'i':
        use16bit = true;
        break;
      case 'f':
        use16bit = false;
        break;
      case 'd':
        showDbg = true;
        break;
      case 'c':
        showCw = true;
        break;
      case 'h':
        fprintf(stderr, "CSDR-Based CW Skimmer by Marat Fayzullin\n");
        fprintf(stderr, "Usage: %s [options] [<infile> [<outfile>]]\n", argv[0]);
        fprintf(stderr, "  -r <rate>  -- Use given sampling rate.\n");
        fprintf(stderr, "  -n <chars> -- Number of characters to print.\n");
        fprintf(stderr, "  -i         -- Use 16bit signed integer input.\n");
        fprintf(stderr, "  -f         -- Use 32bit floating point input.\n");
        fprintf(stderr, "  -c         -- Print dits and dahs to STDOUT.\n");
        fprintf(stderr, "  -d         -- Print debug information to STDERR.\n");
        fprintf(stderr, "  -h         -- Print this help message.\n");
        return(0);
      default:
        fprintf(stderr, "%s: Unrecognized option '%s'!\n", argv[0], argv[j]);
        return(2);
    }
  }

  // Open input and output files
  inFile = inName? fopen(inName, "rb") : stdin;
  if(!inFile)
  {
    fprintf(stderr, "%s: Failed opening input file '%s'\n", argv[0], inName);
    return(1);
  }
  outFile = outName? fopen(outName, "wb") : stdout;
  if(!outFile)
  {
    fprintf(stderr, "%s: Failed opening output file '%s'\n", argv[0], outName);
    if(inFile!=stdin) fclose(inFile);
    return(1);
  }

  // Allocate FFT plan, input, and output buffers
  fftwf_complex *fftOut = new fftwf_complex[MAX_INPUT];
  short *dataIn = new short[MAX_INPUT];
  float *fftIn  = new float[MAX_INPUT];
  fftwf_plan fft = fftwf_plan_dft_r2c_1d(MAX_INPUT, fftIn, fftOut, FFTW_ESTIMATE);

  // Allocate CSDR object storage
  in        = new Csdr::Ringbuffer<float> *[MAX_CHANNELS];
  inReader  = new Csdr::RingbufferReader<float> *[MAX_CHANNELS];
  out       = new Csdr::Ringbuffer<unsigned char> *[MAX_CHANNELS];
  outReader = new Csdr::RingbufferReader<unsigned char> *[MAX_CHANNELS];
  cwDecoder = new Csdr::CwDecoder<float> *[MAX_CHANNELS];

  // Debug output gets accumulated here
  char dbgOut[MAX_CHANNELS+16];

  // Create CSDR objects
  for(j=0 ; j<MAX_CHANNELS ; ++j)
  {
    in[j]        = new Csdr::Ringbuffer<float>(sampleRate);
    inReader[j]  = new Csdr::RingbufferReader<float>(in[j]);
    out[j]       = new Csdr::Ringbuffer<unsigned char>(printChars*4);
    outReader[j] = new Csdr::RingbufferReader<unsigned char>(out[j]);
    cwDecoder[j] = new Csdr::CwDecoder<float>(sampleRate, showCw);
    cwDecoder[j]->setReader(inReader[j]);
    cwDecoder[j]->setWriter(out[j]);
  }

  // Read and decode input
  for(remains=0, avgPower=1.0 ; ; )
  {
    if(!use16bit)
    {
      n = fread(fftIn+remains, sizeof(float), MAX_INPUT-remains, inFile);
      if(n!=MAX_INPUT-remains) break;
    }
    else
    {
      n = fread(dataIn+remains, sizeof(short), MAX_INPUT-remains, inFile);
      if(n!=MAX_INPUT-remains) break;
      // Expand shorts to floats, normalizing them to [-1;1) range
      for(j=remains ; j<MAX_INPUT ; ++j)
        fftIn[j] = (float)dataIn[j] / 32768.0;
    }

    // Compute FFT
    fftwf_execute(fft);

    // Shift input data
    remains = MAX_INPUT-INPUT_STEP;
    memcpy(fftIn, fftIn+INPUT_STEP, remains*sizeof(float));

    // Go to magnitudes, compute average
    for(j=0, accPower=0.0 ; j<MAX_INPUT/2 ; ++j)
    {
      float power = sqrt(fftOut[j][0]*fftOut[j][0] + fftOut[j][1]*fftOut[j][1]);
      fftOut[j][1] = power;
      accPower += power;
    }

    // Filter out spurs
#if USE_NEIGHBORS
    fftOut[MAX_INPUT/2-1][0] = fmax(0.0, fftOut[MAX_INPUT/2-1][1] - NEIGH_WEIGHT * fftOut[MAX_INPUT/2-2][1]);
    fftOut[0][0] = fmax(0.0, fftOut[0][1] - NEIGH_WEIGHT * fftOut[1][1]);
#endif
    accPower = fftOut[0][0] + fftOut[MAX_INPUT/2-1][0];
    maxPower = fmax(fftOut[0][0], fftOut[MAX_INPUT/2-1][0]);
    for(j=1 ; j<MAX_INPUT/2-1 ; ++j)
    {
#if USE_NEIGHBORS
      fftOut[j][0] = fmax(0.0, fftOut[j][1] - 0.5 * NEIGH_WEIGHT * (fftOut[j-1][1] + fftOut[j+1][1]));
#endif
      accPower += fftOut[j][0];
      maxPower  = fmax(maxPower, fftOut[j][0]);
    }

    // Maintain rolling average over AVG_SECONDS
    accPower /= MAX_INPUT/2;
    avgPower += (accPower - avgPower) * INPUT_STEP / sampleRate / AVG_SECONDS;

    // Decode by channel
    for(j=i=k=n=0, accPower=0.0 ; j<MAX_INPUT/2 ; ++j, ++n)
    {
      float power = fftOut[j][0];

      // If accumulated enough FFT buckets for a channel...
      if(k>=MAX_INPUT/2)
      {
        // Compute signal power for a channel
#if USE_AVG_BOTTOM
        accPower  = fmax(0.0, accPower - avgPower);
#endif
#if USE_AVG_RATIO
        accPower  = fmax(1.0, accPower / avgPower);
#endif
        dbgOut[i] = accPower<0.5? '.' : '0' + round(fmax(fmin(accPower / maxPower * 10.0, 9.0), 0.0));

        // If CW input buffer can accept samples...
        if(in[i]->writeable()>=INPUT_STEP)
        {
          // Fill input buffer with computed signal power
          float *dst = in[i]->getWritePointer();
          for(int j=0 ; j<INPUT_STEP ; ++j) dst[j] = accPower;
          in[i]->advance(INPUT_STEP);

          // Process input for the current channel
          while(cwDecoder[i]->canProcess()) cwDecoder[i]->process();

          // Print output
          printOutput(outFile, i, i * sampleRate / 2 / MAX_CHANNELS, printChars);
        }

        // Start on a new channel
        accPower = 0.0;
        k -= MAX_INPUT/2;
        i += 1;
        n  = 0;
      }

      // Maximize channel signal power
      accPower = fmax(accPower, power);
      k += MAX_CHANNELS;
    }

    // Print debug information to the stderr
    dbgOut[i] = '\0';
    if(showDbg) fprintf(stderr, "%s (%.2f, %.2f)\n", dbgOut, avgPower, maxPower);
  }

  // Final printout
  for(i=0 ; i<MAX_CHANNELS ; i++)
    printOutput(outFile, i, i * sampleRate / 2 / MAX_CHANNELS, 1);

  // Close files
  if(outFile!=stdout) fclose(outFile);
  if(inFile!=stdin)   fclose(inFile);

  // Release FFTW3 resources
  fftwf_destroy_plan(fft);
  delete [] fftOut;
  delete [] fftIn;
  delete [] dataIn;

  // Release CSDR resources
  for(j=0 ; j<MAX_CHANNELS ; ++j)
  {
    delete outReader[j];
    delete out[j];
    delete cwDecoder[j];
    delete inReader[j];
    delete in[j];
  }

  // Release CSDR object storage
  delete [] in;
  delete [] inReader;
  delete [] out;
  delete [] outReader;
  delete [] cwDecoder;

  // Done
  return(0);
}
