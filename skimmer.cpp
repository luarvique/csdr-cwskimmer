#include "csdr/ringbuffer.hpp"
#include "csdr/cw.hpp"
#include "fftw3.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

#define MAX_CHANNELS (64)
#define MAX_INPUT    (MAX_CHANNELS*4)
#define INPUT_STEP   (MAX_INPUT)//MAX_INPUT/4)

Csdr::Ringbuffer<float> *in[MAX_CHANNELS];
Csdr::RingbufferReader<float> *in_reader[MAX_CHANNELS];
Csdr::Ringbuffer<unsigned char> *out[MAX_CHANNELS];
Csdr::RingbufferReader<unsigned char> *out_reader[MAX_CHANNELS];
Csdr::CwDecoder<float> *cw[MAX_CHANNELS];

unsigned int bufSize = 1024;
unsigned int sampleRate = 48000;
unsigned int charCount = 8;
bool use16bit = false;
bool showCw = false;


int main(int argc, char *argv[])
{
  float accPower, avgPower;
  int j, i, k, n, chrCount;
  int remains;

  for(j=1 ; j<argc ; ++j)
  {
    if(argv[j][0]!='-')
    {
    }
    else if(strlen(argv[j])!=2)
    {
    }
    else switch(argv[j][1])
    {
      case 'n':
        charCount = j<argc-1? atoi(argv[++j]) : charCount;
        charCount = charCount<1? 1 : charCount>32? 32 : charCount;
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
      default:
        break;
    }
  }

  for(j=0 ; j<MAX_CHANNELS ; ++j)
  {
    in[j]         = new Csdr::Ringbuffer<float>(sampleRate);
    in_reader[j]  = new Csdr::RingbufferReader<float>(in[j]);
    out[j]        = new Csdr::Ringbuffer<unsigned char>(bufSize);
    out_reader[j] = new Csdr::RingbufferReader<unsigned char>(out[j]);
    cw[j]         = new Csdr::CwDecoder<float>(sampleRate, showCw);
    cw[j]->setReader(in_reader[j]);
    cw[j]->setWriter(out[j]);
  }

  fftwf_complex *fftOut = new fftwf_complex[MAX_INPUT];
  short *dataIn = new short[MAX_INPUT];
  float *fftIn = new float[MAX_INPUT];
  fftwf_plan fft = fftwf_plan_dft_r2c_1d(MAX_INPUT, fftIn, fftOut, FFTW_ESTIMATE);

  for(remains=0, avgPower=0.0 ; ; )
  {
    if(!use16bit)
    {
      n = fread(fftIn+remains, sizeof(float), MAX_INPUT-remains, stdin);
      if(n!=MAX_INPUT-remains) break;
    }
    else
    {
      n = fread(dataIn+remains, sizeof(short), MAX_INPUT-remains, stdin);
      if(n!=MAX_INPUT-remains) break;
      // Expand shorts to floats, normalizing them to [-1;1) range
      for(j=remains ; j<MAX_INPUT ; ++j)
        fftIn[j] = (fftIn[j] + (float)(dataIn[j]) / 32768.0) / 2.0;
    }

    // Compute FFT
    fftwf_execute(fft);

    // Shift input data
    remains = MAX_INPUT-INPUT_STEP;
    memcpy(fftIn, fftIn+INPUT_STEP, remains*sizeof(float));

    char diag[MAX_CHANNELS+16];

    // Go to magnitudes, compute average
    for(j=0, accPower=0.0 ; j<MAX_INPUT/2 ; ++j)
    {
      float power = sqrt(fftOut[j][0]*fftOut[j][0] + fftOut[j][1]*fftOut[j][1]);
      fftOut[j][1] = power;
      accPower += power;
    }

    // Filter out spurs
    fftOut[MAX_INPUT/2-1][0] = 0.5 * fftOut[MAX_INPUT/2-2][1] + fftOut[MAX_INPUT/2-1][1];
    fftOut[0][0] = fftOut[0][1] + 0.5 * fftOut[1][1];
    accPower = fftOut[0][0] + fftOut[MAX_INPUT/2-1][0];
    for(j=1 ; j<MAX_INPUT/2-1 ; ++j)
    {
      fftOut[j][0] = fftOut[j][1] + 0.5 * (fftOut[j-1][1] + fftOut[j+1][1]);
      accPower += fftOut[j][0];
    }

    // Maintain rolling average
    accPower /= MAX_INPUT/2;
    avgPower  = fmax(avgPower, accPower) * 0.9999;

    // Decode by channel
    for(j=i=k=chrCount=0, n=1, accPower=0.0 ; j<MAX_INPUT/2 ; ++j, ++n, k+=MAX_CHANNELS)
    {
      float power = fftOut[j][0];

      if(k>=MAX_INPUT/2)
      {
        accPower /= n;
        accPower  = accPower<avgPower*2.0? 0.0 : 1.0;
        diag[i]   = accPower<0.5? '.' : '0' + round(fmax(fmin((accPower-0.5) * 18.0, 9.0), 0.0));

        // If filtered input buffer can accept samples...
        if(in[i]->writeable()>=INPUT_STEP)
        {
          // Fill input buffer with computed input power
          float *dst = in[i]->getWritePointer();
          for(int j=0 ; j<INPUT_STEP ; ++j) dst[j] = accPower;
          in[i]->advance(INPUT_STEP);

          // Process input for the current channel
          while(cw[i]->canProcess()) cw[i]->process();

          // Print output
          int m = out_reader[i]->available();
          if(m>=charCount)
          {
            unsigned char *p = out_reader[i]->getReadPointer();
            printf("%d:", i * sampleRate / 2 / MAX_CHANNELS);
            for(int j=0 ; j<m ; ++j)
            {
              if((j<m-1) && (p[j+1]==' ') && ((p[j]=='E') || (p[j]=='T')))
                ++j;
              else
                printf("%c", p[j]);
            }
            out_reader[i]->advance(m);
            chrCount += m;
            printf("\n");
            fflush(stdout);
          }
        }

        accPower = 0.0;
        k -= MAX_INPUT/2;
        i += 1;
        n  = 0;
      }

      accPower += power;
    }

    // Print debug information to the stderr
    diag[i] = '\0';
//    fprintf(stderr, "%s (%.2f - %d)\n", diag, avgPower, chrCount);
  }

  // Final printout
  for(i=0 ; i<MAX_CHANNELS ; i++)
  {
    // Print output
    int n = out_reader[i]->available();
    if(n)
    {
      unsigned char *p = out_reader[i]->getReadPointer();
      printf("%d:", i * sampleRate / 2 / MAX_CHANNELS);
      for(int j=0 ; j<n ; ++j)
        if((j<n-1) && (p[j+1]==' ') && ((p[j]=='E') || (p[j]=='T')))
          ++j;
        else
          printf("%c", p[j]);
      out_reader[i]->advance(n);
      printf("\n");
      fflush(stdout);
    }
  }

  // Release FFTW3 resources
  fftwf_destroy_plan(fft);
  delete fftOut;
  delete fftIn;
  delete dataIn;

  // Release CSDR resources
  for(j=0 ; j<MAX_CHANNELS ; ++j)
  {
    delete out_reader[j];
    delete out[j];
    delete cw[j];
    delete in_reader[j];
    delete in[j];
  }

  // Done
  return(0);
}

