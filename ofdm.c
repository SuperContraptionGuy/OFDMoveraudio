// compile: gcc -lm qam.c
#define _GNU_SOURCE
#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <string.h>
#include <errno.h>
#include <time.h>

#include <complex.h>
#include <fftw3.h>

#include "utilities.h"

#define WARN_UNUSED __attribute__((warn_unused_result))


/*
 *  performs a DFT based linear convolution using overlap and save method
 *  length of inputSamples should be the impulse response length plus the filter block length
 *  I'll probably use a filter block length equal to the impulse response length
 *
*/
buffered_data_return_t channelFilter(const overlap_save_buffer_double_t *inputSamples, sample_double_t *outputSample)
{
    static circular_buffer_double_t  impulseResponse = {0};                 // time domain of impulse response
    static circular_buffer_complex_t frequencyResponse = {0};               // fft of impulse response

    static circular_buffer_complex_t frequencyDomainSamples = {0};          // fft of inputSamples
    static circular_buffer_double_t  timeDomainFilteredSamples = {0};       // output of the filter, feed to outputSample one at a time

    // the fftw plans
    static fftw_plan fftwImpulseToFrequency;
    static fftw_plan fftwSamplesToFrequency;
    static fftw_plan fftwFilteredSamplesToTime;

    static int initialized = 0; // have we generated the filter taps
    static int bufferPrimed = 0;    // have we waited for enough samples to start filtering

    if(!initialized)
    {
        // intitialize the fftw related buffers

        impulseResponse.length = inputSamples->N;
        if((impulseResponse.buffer = fftw_malloc(sizeof(double) * impulseResponse.length)) == NULL)
            fprintf(stderr, "Couldn't allocate memory for impulseResponse buffer. %s\n", strerror(errno));

        timeDomainFilteredSamples.length = inputSamples->N;
        if((timeDomainFilteredSamples.buffer = fftw_malloc(sizeof(double) * timeDomainFilteredSamples.length)) == NULL)
            fprintf(stderr, "Couldn't allocate memory for impulseResponse buffer. %s\n", strerror(errno));

        frequencyResponse.length = inputSamples->N / 2 + 1; // smaller in frequency space because of real data
        //if((frequencyResponse.buffer = calloc(frequencyResponse.length, sizeof(double complex))) == NULL)
        if((frequencyResponse.buffer = fftw_malloc(sizeof(fftw_complex) * frequencyResponse.length)) == NULL)
            fprintf(stderr, "Couldn't allocate memory for frequency response buffer: %s\n", strerror(errno));

        frequencyDomainSamples.length = inputSamples->N / 2 + 1; // smaller in frequency space because of real data
        //if((frequencyDomainSamples.buffer = calloc(frequencyDomainSamples.length, sizeof(double complex))) == NULL)
        if((frequencyDomainSamples.buffer = fftw_malloc(sizeof(fftw_complex) * frequencyDomainSamples.length)) == NULL)
            fprintf(stderr, "Couldn't allocate memory for frequency response buffer: %s\n", strerror(errno));

        fftwImpulseToFrequency =      fftw_plan_dft_r2c_1d(inputSamples->N, impulseResponse.buffer, frequencyResponse.buffer, FFTW_ESTIMATE);   // for calculating frequency response from impulse response.           used once
        fftwSamplesToFrequency =      fftw_plan_dft_r2c_1d(inputSamples->N, inputSamples->M_buffer, frequencyDomainSamples.buffer, FFTW_MEASURE);  // for calculating frequency domain of input samples.                used many times
        fftwFilteredSamplesToTime =   fftw_plan_dft_c2r_1d(inputSamples->N, frequencyDomainSamples.buffer, timeDomainFilteredSamples.buffer, FFTW_MEASURE);     // for calculating time domain of filtered samples    used many times

        // fill impulse response buffer from file
        // pull the samples from a file called impulseResponse.raw
        int impulseResponseFile = open("impulseResponse.raw", O_RDONLY);
        if(impulseResponseFile != -1)   // check that the file is readable
        {
            // file format:
            //  S32_LE PCM mono samples
            //  so 4 bytes per sample
            //  generated like
            //      ffmpeg -i impulseResponse.wav -f s32le pipe:1 > impulseResponse.raw
            for(int i = 0; i < inputSamples->M; i++)
            {
                // for type punning
                union __attribute__((packed))
                {
                   int32_t value;
                   struct
                   {
                        uint8_t bytes[4];
                   };
                } readSample;

                // read 4 bytes at a time into the type punning structure
                int readBytes = {0};
                if((readBytes = read(impulseResponseFile, &readSample.bytes, sizeof(readSample.bytes))) == 0)
                {
                    //fprintf(stderr, "Reached end of impulseResponse.raw before filter was satisfied! Attempted reading sample #%d out of %d. Zero filling.\n", i, inputSamples->M);
                    //readSample.value = 0;
                }

                impulseResponse.buffer[i] = (double)readSample.value / INT32_MAX;    // convert to a double between -1 and 1
            }
            close(impulseResponseFile);

        } else {
            // if the file isn't readable
            fprintf(stderr, "unable to open data/impulseResponse.raw file for channel equalization filter.\nUsing a default degenerate impulse response.\n");

            // make up a default impulse response, ie, 1 at the start, and zeros everywhere else
            for(int i = 0; i < inputSamples->M; i++)
            {
                impulseResponse.buffer[i] = 0;
            }
            impulseResponse.buffer[0] = 1;
        }

        // set remaining values to zero
        for(int i = 0; i < inputSamples->L; i++)
        {
            impulseResponse.buffer[i + inputSamples->M] = 0;
        }

        // initialize the frequency response buffer
        fftw_execute(fftwImpulseToFrequency);

        // run once
        initialized = 1;
    }

    // wait for enough samples to come in to do the first DFT, then we can return the samples from the fft buffer one at a time
    //if(inputSamples->n < inputSamples->length)
    if(inputSamples->n < inputSamples->L - 1)
        return AWAITING_SAMPLES;

    // return a sample from the dft buffer
    // the very first sample will just be 0. could be an issue, it's becasue the fft hasn't run for the first time yet
    //outputSample->sample = timeDomainFilteredSamples.buffer[inputSamples->n % timeDomainFilteredSamples.length + M];    // output one sample from the buffer, avoiding the first M samples always
    outputSample->sample = timeDomainFilteredSamples.buffer[(inputSamples->n % inputSamples->L) + inputSamples->M];    // output one sample from the buffer, avoiding the first M samples always, and just using the last L samples
    // normalize the output because we did ifft(fft(x)), so factor is 1/N
    outputSample->sample /= inputSamples->N;
    outputSample->sampleIndex = inputSamples->n - inputSamples->L;
    outputSample->sampleRate = inputSamples->sampleRate;

    // after filling the L buffer, do an fft

    // if the input buffer has filled up again, recalculate the frequency domain, then multiply by frequency response, and do the inverse transform
    //if(inputSamples->n % inputSamples->length == 0)
    if(inputSamples->n % inputSamples->L == inputSamples->L - 1)    // go when the last sample is put into L
    {
        // to the dft on input samples
        fftw_execute(fftwSamplesToFrequency);
        // multiply by frequency response
        for(int i = 0; i < frequencyDomainSamples.length; i++)
        {
            frequencyDomainSamples.buffer[i] = frequencyDomainSamples.buffer[i] * frequencyResponse.buffer[i];
        }
        // to the inverse transform
        fftw_execute(fftwFilteredSamplesToTime);

        // then copy the end of the L buffer into M for next fft
        for(int i = 0; i < inputSamples->M; i++)
        {
            inputSamples->M_buffer[i] = inputSamples->M_next_buffer[i];
        }

    }

    return RETURNED_SAMPLE;
}

// runs the state machine for an OFDM based communication transmitter
buffered_data_return_t OFDM(int long n, sample_double_t *outputSample, OFDM_state_t *OFDMstate)
{
    if(OFDMstate->initialized == 0)
    {
        initializeOFDMstate(OFDMstate);

        // initialize the state variables
        OFDMstate->state.frame = IDLE;
        OFDMstate->state.frameStart = n;
        
        // initialize channel simulation filter
        OFDMstate->simulateNoise = 0;
        OFDMstate->simulateChannel = 0;

        OFDMstate->dataInput = fopen("inputData", "r");
        OFDMstate->bitOffset = 0;
        OFDMstate->generatedDataOutput = fopen("senderSequence", "w");


        if(OFDMstate->simulateChannel)
        {
            // initialize the overlap and save buffer
            OFDMstate->channelSimulationBuffer.M = 5912;
            OFDMstate->channelSimulationBuffer.L = OFDMstate->channelSimulationBuffer.M;
            OFDMstate->channelSimulationBuffer.length = OFDMstate->channelSimulationBuffer.L;
            OFDMstate->channelSimulationBuffer.N = OFDMstate->channelSimulationBuffer.L + OFDMstate->channelSimulationBuffer.M;
            OFDMstate->channelSimulationBuffer.insertionIndex = 0;

            OFDMstate->channelSimulationBuffer.M_buffer = fftw_malloc(sizeof(double) * OFDMstate->channelSimulationBuffer.N);    // use fftw malloc since fftw will use this buffer
            if(OFDMstate->channelSimulationBuffer.M_buffer == NULL)
                fprintf(stderr, "cahnnelSimulationBuffer failed to allocate: %s\n", strerror(errno));

            OFDMstate->channelSimulationBuffer.L_buffer = OFDMstate->channelSimulationBuffer.M_buffer + OFDMstate->channelSimulationBuffer.M;   // offset the L
            OFDMstate->channelSimulationBuffer.buffer = OFDMstate->channelSimulationBuffer.L_buffer;
            OFDMstate->channelSimulationBuffer.M_next_buffer = OFDMstate->channelSimulationBuffer.M_buffer + OFDMstate->channelSimulationBuffer.L;  // M_next_buffer is the last M numbers of L
        }


        // initialize fftw plan, using the measure option to calculate the fastest plan, could have a few seconds startup time
        OFDMstate->fftwPlan = fftw_plan_dft_c2r_1d(
                OFDMstate->OFDMsymbol.timeDomain.length,
                OFDMstate->OFDMsymbol.frequencyDomain.buffer,
                OFDMstate->OFDMsymbol.timeDomain.buffer,
                FFTW_MEASURE);


        OFDMstate->initialized = 1;
    }

    double output = 0;

    switch(OFDMstate->state.frame)
    {
        case IDLE:

            // send nothing
            output = 0;
            
            // wait for enough data in the pipe to generate a frame, or wait for enough time to pass to send a frame anyway.

            // check for exit from IDLE frame
            if(n - OFDMstate->state.frameStart >= OFDMstate->symbolPeriod * 3 - 1) // example of state change based on timing
            {
                OFDMstate->state.frame = ACTIVE;
                OFDMstate->state.frameStart = n + 1;    // starts next index
                OFDMstate->state.fieldIndex = 0;
                OFDMstate->state.field = PREAMBLE;
                OFDMstate->state.fieldStart = n + 1;
                OFDMstate->state.symbol = GUARD_PERIOD;
                OFDMstate->state.symbolStart = n + 1;
                OFDMstate->state.symbolIndex = 0;
            }
            break;

        case ACTIVE:
            switch(OFDMstate->state.field)
            {
                case PREAMBLE:

                    switch(OFDMstate->state.symbol)
                    {
                        case GUARD_PERIOD:

                            if(n - OFDMstate->state.symbolStart == 0)
                            {
                                // set new symbol for testing
                                // generate a symetric symbol, with even frequency components only for syncronization
                                // skipping DC and Niquist frequencies
                                for(int k = 1; k < OFDMstate->channels - 1; k++)
                                {
                                    // generate a random integer for constellation choices
                                    long int randomInteger;

                                    if(OFDMstate->state.symbolIndex < 2)    // generate symetric symbols for the first two symbols
                                    {
                                        if(k % 2 == 0)    // every other for repetative symbol
                                        //if(1)
                                        //if(k > 16 && k < 100 && k % 2 == 0)
                                        //if(k > 100 && k < 200 && k % 2 == 0)
                                        //if(k > 446 && k < 446+200 && k % 2 == 0)
                                        //if(k == OFDMstate->channels / 16)
                                        //if(k == 0)
                                        //if(0)
                                        {
                                            constellation_complex_t constellation = OFDMstate->constellations[0];
                                            lrand48_r(&OFDMstate->preamblePilotsPRNG, &randomInteger);
                                            OFDMstate->OFDMsymbol.frequencyDomain.buffer[k] = constellation.points[randomInteger % constellation.length];

                                        } else {
                                            OFDMstate->OFDMsymbol.frequencyDomain.buffer[k] = 0;
                                        }
                                        //OFDMstate->OFDMsymbol.frequencyDomain[k] *= (double)OFDMstate->channels / 200 / 30;   // normalization factor for fewer channels used
                                        OFDMstate->OFDMsymbol.frequencyDomain.buffer[k] *= M_SQRT2;  // scale factor for even only channels

                                        // time domain samples are now in the OFDMstate->OFDMsymbol.timeDomain array

                                    } else if(OFDMstate->state.symbolIndex % 2 == 0)    // Then generate two duplicate symbols
                                    {
                                        if(1)   // pick all subchannels
                                        //if(k > 100 && k < 200 && k % 2 == 0)
                                        //if(k > 446 && k < 446+10)
                                        {
                                            //OFDMstate->OFDMsymbol.frequencyDomain[k] =
                                            //    rand() % 2 * 2 - 1 +
                                            //    I*(rand() % 2 * 2 - 1);
                                                //0;
                                            constellation_complex_t constellation = OFDMstate->constellations[0];
                                            lrand48_r(&OFDMstate->preamblePilotsPRNG, &randomInteger);
                                            OFDMstate->OFDMsymbol.frequencyDomain.buffer[k] = constellation.points[randomInteger % constellation.length];

                                            fprintf(OFDMstate->generatedDataOutput, "n=%i k=%i %li: %lf+%lfi\n", OFDMstate->state.symbolIndex, k, randomInteger, creal(OFDMstate->OFDMsymbol.frequencyDomain.buffer[k]), cimag(OFDMstate->OFDMsymbol.frequencyDomain.buffer[k]));
                                        } else {
                                            OFDMstate->OFDMsymbol.frequencyDomain.buffer[k] =
                                                0;
                                        }
                                        //OFDMstate->currentOFDMSymbol[k] = 0;    // testing
                                        //OFDMstate->OFDMsymbol.frequencyDomain[k] *= (double)OFDMstate->channels / 10 / 30;
                                    }
                                }
                                // now do the transform
                                if(OFDMstate->state.symbolIndex%2 == 0 || OFDMstate->state.symbolIndex < 2)    // only do fft for symetric symbols and first of the duplicate symbols
                                    fftw_execute(OFDMstate->fftwPlan);
                                // time domain samples are now in the OFDMstate->OFDMsymbol.timeDomain array
                            }

                            // output samples for the guard period
                            //output = OFDMsymbolBaseband(OFDMstate->ofdmPeriod - OFDMstate->guardPeriod + (n - OFDMstate->state.symbolStart),
                                                        //OFDMstate);
                            output = OFDMstate->OFDMsymbol.timeDomain.buffer[OFDMstate->ofdmPeriod - OFDMstate->guardPeriod + (n - OFDMstate->state.symbolStart)];

                            // check for exit form GUARD_PERIOD symbol
                            if(n - OFDMstate->state.symbolStart >= OFDMstate->guardPeriod - 1)
                            {
                                OFDMstate->state.symbol = OFDM_PERIOD;
                                OFDMstate->state.symbolStart = n + 1;
                            }
                            break;
                        case OFDM_PERIOD:

                            //output = OFDMsymbolBaseband(n - OFDMstate->state.symbolStart,
                                                        //OFDMstate);
                            output = OFDMstate->OFDMsymbol.timeDomain.buffer[n - OFDMstate->state.symbolStart];

                            // check for exit
                            if(n - OFDMstate->state.symbolStart >= OFDMstate->ofdmPeriod - 1)
                            {
                                OFDMstate->state.symbol = GUARD_PERIOD;
                                OFDMstate->state.symbolStart = n + 1;
                                OFDMstate->state.symbolIndex++;
                            }
                            break;
                    }

                    // check for exit from PREAMBLE field
                    if(n - OFDMstate->state.fieldStart >= OFDMstate->symbolPeriod * 4 - 1)   // for 3 symbol preamble
                    {
                        OFDMstate->state.field = DATA;
                        OFDMstate->state.fieldStart = n + 1;
                        OFDMstate->state.symbolIndex = 0;
                    }
                    break;
                case DATA:


                    switch(OFDMstate->state.symbol)
                    {
                        case GUARD_PERIOD:

                            if(n - OFDMstate->state.symbolStart == 0)
                            {
                                // set new OFDM symbol for testing, that's one IQ pair for each subchannel
                                // a new symbol is chosen for transmission at the beginning of it's guard period.
                                // right now it just chooses a random symbol, but presumably you'd choose this based on some data input
                                // it can be modulated with any IQ method. QPSK is one idea, but you could choose any IQ constellation to encode data. it could also be
                                // a different constellation for each sub channel, useful for taking advantage of low noise subchannels without increasing
                                // error rates on noisy channels
                                /*
                                // choosing random blocks of channels to tranmit on for fun
                                static int center;
                                static int width;
                                if(rand() % 5 == 0 || OFDMstate->state.symbolIndex == 0)
                                {
                                    center = rand() % OFDMstate->channels;
                                    width = rand() % OFDMstate->channels;
                                }
                                */
                                // skip DC and Niquist
                                for(int k = 1; k < OFDMstate->channels - 1; k++)
                                {
                                    // generate a random index
                                    long int randomIntegerPilot;
                                    long int randomIntegerData;
                                    //lrand48_r(&OFDMstate->pilotsPRNG, &randomIntegerPilot);
                                    //lrand48_r(&OFDMstate->predefinedDataPRNG, &randomIntegerData);

                                    if(k % OFDMstate->pilotSymbolsPitch == 0) // transmit pilot symbols on pilot channels
                                    {
                                        // pilot symbols
                                        constellation_complex_t constellation = OFDMstate->constellations[0];
                                        lrand48_r(&OFDMstate->pilotsPRNG, &randomIntegerPilot);

                                        OFDMstate->OFDMsymbol.frequencyDomain.buffer[k] = constellation.points[randomIntegerPilot % constellation.length];
                                        fprintf(OFDMstate->generatedDataOutput, "n=%i k=%i %li: %lf+%lfi\n", OFDMstate->state.symbolIndex, k, randomIntegerPilot, creal(OFDMstate->OFDMsymbol.frequencyDomain.buffer[k]), cimag(OFDMstate->OFDMsymbol.frequencyDomain.buffer[k]));

                                    } else  // and use the rest for data channels 
                                    if(1)   // transmit on all other subchannels
                                    //if(k > 250 && k < 6000)
                                    //if(k > 250 && k < 350)
                                    //if(k > 7000 && k < 8000) // using a select number of channels to simplify the signal for testing
                                    //int startChunk = (OFDMstate->state.symbolIndex / 1 * 100) % OFDMstate->channels / 4;
                                    //if(k % (OFDMstate->channels / 4) > startChunk && k % (OFDMstate->channels / 4) < startChunk + 10)
                                    // choose random blocks
                                    //if(k < center + width / 2 && k > center - width / 2)
                                    {
                                        
                                        // picking a sequential constellation, and a random point in that constellation discluding the first entry
                                        constellation_complex_t constellation = OFDMstate->constellations[k%(OFDMstate->constellationsLength - 1) + 1];
                                        lrand48_r(&OFDMstate->predefinedDataPRNG, &randomIntegerData);
                                        OFDMstate->OFDMsymbol.frequencyDomain.buffer[k] = constellation.points[randomIntegerData % constellation.length];
                                        fprintf(OFDMstate->generatedDataOutput, "n=%i k=%i %li: %lf+%lfi\n", OFDMstate->state.symbolIndex, k, randomIntegerData, creal(OFDMstate->OFDMsymbol.frequencyDomain.buffer[k]), cimag(OFDMstate->OFDMsymbol.frequencyDomain.buffer[k]));

                                        //rescale if I'm using only a few channels for higher power 
                                        //OFDMstate->OFDMsymbol.frequencyDomain[k] *= (double)OFDMstate->channels / 10 / 30;
                                        //OFDMstate->OFDMsymbol.frequencyDomain[k] /= 3;

                                    //} else {
                                        //OFDMstate->OFDMsymbol.frequencyDomain[k] = 0;
                                    } else {
                                        OFDMstate->OFDMsymbol.frequencyDomain.buffer[k] = 0;
                                    }
                                }
                                // now do the transform
                                fftw_execute(OFDMstate->fftwPlan);
                                // time domain samples are now in the OFDMstate->OFDMsymbol.timeDomain array
                            }
                            //output = OFDMsymbolBaseband(OFDMstate->ofdmPeriod - OFDMstate->guardPeriod + (n - OFDMstate->state.symbolStart),
                                                        //OFDMstate);
                            output = OFDMstate->OFDMsymbol.timeDomain.buffer[OFDMstate->ofdmPeriod - OFDMstate->guardPeriod + (n - OFDMstate->state.symbolStart)];

                            //output = output * 0.5;

                            // check for exit
                            if(n - OFDMstate->state.symbolStart >= OFDMstate->guardPeriod - 1)
                            {
                                OFDMstate->state.symbol = OFDM_PERIOD;
                                OFDMstate->state.symbolStart = n + 1;
                            }
                            break;
                        case OFDM_PERIOD:

                            output = OFDMstate->OFDMsymbol.timeDomain.buffer[n - OFDMstate->state.symbolStart];

                            // check for exit
                            if(n - OFDMstate->state.symbolStart >= OFDMstate->ofdmPeriod - 1)
                            {
                                OFDMstate->state.symbol = GUARD_PERIOD;
                                OFDMstate->state.symbolStart = n + 1;
                                OFDMstate->state.symbolIndex++;
                            }
                            break;
                    }

                    // check for exit from DATA field
                    if(n - OFDMstate->state.fieldStart >= OFDMstate->symbolPeriod * 1000 - 1) // set number of symbols of data for example
                    {
                        OFDMstate->state.frame = IDLE;
                        OFDMstate->state.frameStart = n + 1;
                    }
                    break;
            }
            break;
    }

    //double normalizationFactor = sqrt(OFDMstate->ofdmPeriod);   // normalization factor due to the inverse transform
    double normalizationFactor = sqrt(OFDMstate->ofdmPeriod) * 10;   // normalization factor due to the inverse transform
    output /= normalizationFactor;

    // channel simulation
    sample_double_t equalizedSample;
    // decide whether to simulate the channel response or not
    if(OFDMstate->simulateChannel)
    {
        // channel simulation filter
        OFDMstate->channelSimulationBuffer.buffer[OFDMstate->channelSimulationBuffer.insertionIndex] = output;
        OFDMstate->channelSimulationBuffer.n = n;

        buffered_data_return_t returnValue = channelFilter(&OFDMstate->channelSimulationBuffer, &equalizedSample);

        OFDMstate->channelSimulationBuffer.insertionIndex = (OFDMstate->channelSimulationBuffer.insertionIndex + 1) % OFDMstate->channelSimulationBuffer.length;

        if(returnValue != RETURNED_SAMPLE)
            return AWAITING_SAMPLES;


    } else {

        // skip channel simulation
        equalizedSample.sample = output;
    }

    // add noise to output after channel filtering
    if(OFDMstate->simulateNoise)
    {
        double noiseAmplitude = 0.003;
        double randomDouble;
        drand48_r(&OFDMstate->channelNoisePRNG, &randomDouble);
        equalizedSample.sample += randomDouble * noiseAmplitude - noiseAmplitude / 2;
    }

    *outputSample = equalizedSample;

    return RETURNED_SAMPLE;
}

// this is the point where samples are generated
static double WARN_UNUSED calculateSample(int n, OFDM_state_t *OFDMstate)
{

    // I think I need to keep polling for a sample until SAMPLE RETURNED

    //double amplitudeScaler = 0.1;
    double amplitudeScaler = 1;
    /*
    if(n == 2500)
        return 1;
    return 0;
    */


    sample_double_t returnedSample;

    static buffered_data_return_t returnValue = AWAITING_SAMPLES;
    static int offset = 0;
    if(returnValue == AWAITING_SAMPLES)
    {
        // run until first sample is returned
        for(int i = n; returnValue == AWAITING_SAMPLES; i++)
        {
            // call the function as many times as needed until it returns a sample
            returnValue = OFDM(i, &returnedSample, OFDMstate);
            offset = i;
        }
        return returnedSample.sample;
    }
    // from then on just run once per sample
    OFDM(n + offset, &returnedSample, OFDMstate);
    return returnedSample.sample;

    //return OFDM(n, &OFDMstate);
    //return raisedCosQAM(n, sampleRate) * amplitudeScaler;
    //return impulse(n % (int)(44100 * 0.13 * 1.5), 0).I;
    //return singleChannelODFM_noguard(n, sampleRate) * amplitudeScaler;
}

// generates a .wav header of 44 bytes long
// length is the number of samples in the file
static int WARN_UNUSED writeHeader(int length, int fileDescriptor, int sampleRate)
{
    riff_header_t header =
    {
        .riff = "RIFF",
        // chunk size plus rest of file size I think (so exluding first 8 bytes of header)
        .size = length * 4 + sizeof(riff_header_t) - 8,
        .format = "WAVE",
        // fmt chunk
        .chunk = "fmt ",
        .length = 16,
        .type = RIFF_TYPE_PCM,
        .channels = 1,
        .sampleRate = sampleRate,
        .dataRate = 176400,
        .blockSize = 4,
        .bitsPerSample = 32,
        .data = "data",
        .chunkSize = length * 4,
    };

#if DEBUG_LEVEL >= 1
    // dummy test samples
    uint8_t dummy[length * 4];
    memset(dummy, 0, length * 4);

    FILE* hexdumpInput = popen("hexdump -C", "w");

    if (hexdumpInput == NULL)
    {
        fprintf(stderr, "Failed to open hexdump: %s\n", strerror(errno));
        goto exit;
    }

    size_t ret = fwrite(header.bytes, sizeof(riff_header_t), 1, hexdumpInput);

    if (ret != sizeof(riff_header_t))
    {
        fprintf(stderr, "Failed to write to hexdump: %s\n", strerror(errno));
        goto exit;
    }

    ret = fwrite(header.bytes, sizeof(riff_header_t), 1, stdout);

    if (ret != sizeof(riff_header_t))
    {
        fprintf(stderr, "Failed to write to stdout: %s\n", strerror(errno));
        goto exit;
    }

    ret = fwrite(dummy, length * 4, 1, stdout);

    if (ret != sizeof(riff_header_t))
    {
        fprintf(stderr, "Failed to write to stdout: %s\n", strerror(errno));
    }

exit:
    if (hexdumpInput != NULL && fork() == 0)
    {
        pclose(hexdumpInput);
        exit(0);
        return 0;
    }
#endif

    ssize_t rets = write(fileDescriptor, header.bytes, sizeof(riff_header_t));

    if (rets < 0)
    {
        return -1;
    }

    return 0;
}

static int WARN_UNUSED generateSamplesAndOutput(char* filenameInput)
{
    int retval = 0;

    FILE* aplayStdIn = NULL;
    int useAplay = 0;   // if 1, will output to aplay

    int fileDescriptor = -1;

    // audio sample rate
    // Supported sample rates from alsa-info.sh
    //     rates [0x560]: 44100 48000 96000 192000
    int sampleRate = 44100;
    // total number of samples to generate
    long length = sampleRate * 120;
    //long length = (1<<12) * 5 + sampleRate * 0.25;
    // the number of the current sample
    long n = 0;

    OFDM_state_t OFDMstate = {0};
    OFDMstate.sampleRate = sampleRate;

    // length of the file write buffer, samples times 4 bytes per sample
    const int bufferLength = 100 * 4;
    // the file write buffer, used to buffer the write calls
    uint8_t buffer[bufferLength];
    // number of bytes ready to be written out of the buffer in case we need to flush the buffer before it's full
    int bufferReadyBytes = 0;

    // Whether to send samples over stdout or to file
    int outputstd = 0;

    // User passed '-' -> use stdout
    if (filenameInput[0] == '-')
    {
        outputstd = 1;
    }

    // set up the file descriptors for the various outputs

    // setup a file descriptor for a pipe to aplay command to play the sound through the speakers
    char aplayCommandString[30] = {0};
    int len = snprintf(aplayCommandString, 30, "aplay -f S32_LE -c1 -r %i", sampleRate);

    if (len < 0)
    {
        fprintf(stderr, "Failed to write aplay string: %s\n", strerror(errno));
        retval = 2;
        goto exit;
    }
    else if (len == 30)
    {
        fprintf(stderr, "Failed to write aplay string: truncated\n");
        retval = 2;
        goto exit;
    }

    //puts(aplayCommandString);
    if(useAplay)
    {
        aplayStdIn = popen(aplayCommandString, "w");

        if (aplayStdIn == NULL)
        {
            fprintf(stderr, "Failed to open aplay: %s\n", strerror(errno));
            retval = 3;
            goto exit;
        }
    }

#if DEBUG_LEVEL >= 1
    if (outputstd == 0)
    {
        // for the hex dump of printed bytes.
        hexdumpStdIn = popen("hexdump -C", "w");

        if (hexdumpStdIn == NULL)
        {
            fprintf(stderr, "Failed to open hexdump: %s\n", strerror(errno));
            retval = 4;
            goto exit;
        }
    }
    char *plotstr = 
        "feedgnuplot "
        "--domain --dataid --lines --points --maxcurves 100000 "
        "--title \"Debug modulator\" "
        "--xlabel \"Time (n)\" --ylabel \"value\" "
        "--legend 0 \"Symbol Index\" "
        "--legend 1 \"Sample Index\" "
        "--legend 2 \"Generated Audio Sample\" "
        "--legend 3 \"Phase Offset\" "
        "--legend 4 \"Filtered I\" "
        "--legend 5 \"Filtered Q\" "
        "--legend 6 \"I\" "
        "--legend 7 \"Q\" "
    ;
    plotStdIn = popen(plotstr, "w");
#endif

    if (outputstd == 0)
    {
        // for the file writing
        char filename[80] = {0};
        len = snprintf(filename, 80, "%s.wav", filenameInput);

        if (len < 0)
        {
            fprintf(stderr, "Failed to get filename: %s\n", strerror(errno));
            retval = 5;
            goto exit;
        }
        else if (len == 80)
        {
            fprintf(stderr, "Failed to get filename: truncated\n");
            retval = 5;
            goto exit;
        }

        //puts(filename);

        fileDescriptor = open(filename, O_WRONLY | O_CREAT, S_IRUSR | S_IWUSR | S_IRGRP | S_IWGRP | S_IROTH);

        if (fileDescriptor < 0)
        {
            fprintf(stderr, "Failed to open file: %s\n", strerror(errno));
            retval = 6;
            goto exit;
        }
    }

    // first generate the header
    if (outputstd == 0)
    {
        int ret = writeHeader(length, fileDescriptor, sampleRate);

        if (ret != 0)
        {
            fprintf(stderr, "Failed to write header: %d\n", ret);
            retval = 7;
            goto exit;
        }
    }

    // calculate all the samples
    // seed the random number generator
    //srand(time(NULL));
    while(n < length)
    {
        // calculate a chunk of samples until the buffer is full or max is reached. one sample at a time, 4 bytes at a time
        for(bufferReadyBytes = 0; (bufferReadyBytes < bufferLength) && (n < length); bufferReadyBytes += 4, n++)
        {
            // the sample value used in calculations, to be normalized
            double sampleValue;
            // sample value after put into signed integer range, then split into bytes for file writing and audio output
            sample_32_converter_t normalizedSampleValue;
            // holds each individual byte as it's written out Little Endian style
            char byte;

            // get the double sample value, should be between -1 and 1
            sampleValue = calculateSample(n, &OFDMstate);
            // calculate the final signed integer to be output as a sample
            // the magnitude of the max is always one smaller than the magnitude of the min
            normalizedSampleValue.value = sampleValue * INT32_MAX;

            // split up the normalized value into individual bytes
            for(int i = 0; i < 4; i++)
            {
                // get the byte from normalized. pointer points to the adress of the first byte in normalized
                byte = normalizedSampleValue.bytes[i];
                // add byte to the buffer
                buffer[bufferReadyBytes + i] = byte;

                // send to the pipes one byte at a time since they are buffered by the OS
                if(useAplay)
                    putc(byte, aplayStdIn);

            #if DEBUG_LEVEL >= 1
                if (outputstd == 0)
                {
                #if DEBUG_LEVEL > 1
                    putc(byte, hexdumpStdIn);
                #endif
                }
                else
                {
                    putchar(byte);
                }
            #else
                if (outputstd != 0)
                {
                    putchar(byte);
                }
            #endif
            }
        }

        // write the buffer to the file bufferReadyBytes number of bytes, usually a whole buffer full at a time, until the end.
        if (outputstd == 0)
        {
            ssize_t ret = write(fileDescriptor, buffer, bufferReadyBytes);

            if (ret < 0)
            {
                fprintf(stderr, "Failed to write output buffer: %s\n", strerror(errno));
                retval = 8;
                goto exit;
            }
        }
    }

exit:
    if(aplayStdIn)
        pclose(aplayStdIn);
    if (outputstd == 0)
    {
        close(fileDescriptor);

    #if DEBUG_LEVEL > 0
        pclose(hexdumpStdIn);
    #endif
    }

#if DEBUG_LEVEL > 0
    if(plotStdIn != NULL && fork() == 0)
    {
        // it's holding onto a reference to stdout, gotta close that off
        //freopen("/dev/null", "w", stdout);
        // nope, that wasn't it
        pclose(plotStdIn);
        exit(0);
        return 0;
    }
#endif

    return retval;
}

static void usage(const char *filename)
{
    fprintf(stderr, "Provide file name as parameter. For example:\n");
    fprintf(stderr, "  Wrap the file name in quotes:\n");
    fprintf(stderr, "  ./%s \"file name\"\n", filename);
    fprintf(stderr, "Will generate a file called \"file name.wav\".\n");
}

int main(int argc, char** args)
{
    // extract filename from arguments
    if (argc < 2)
    {
        // or use stdout
        usage(args[0]);
        return 1;
    }
    else if (argc > 2)
    {
        fprintf(stderr, "Too many arguments.\n");
        usage(args[0]);
        return 1;
    }

    return generateSamplesAndOutput(args[1]);
}
