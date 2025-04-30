#include <fftw3.h>
#include <stdlib.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <stdint.h>
#include <math.h>
#include <complex.h>
#include <string.h>
#include <errno.h>


// Debug flag
#define DEBUG_LEVEL 2
typedef union __attribute__((packed))
{
    int32_t value;
    uint8_t byte[sizeof(int32_t)];
} sample_32_converter_t;    // union to help convert from bytes to integer to double

typedef struct
{
    double sample;      // double value representation of sample
    int sampleRate;     // rate of samples per second
    int sampleIndex;    // index since begining of record
} sample_double_t;

typedef struct
{
    double complex sample;      // double value representation of sample
    int sampleRate;     // rate of samples per second
    int sampleIndex;    // index since begining of record
} sample_complex_t;

typedef enum
{
    ALLIGNED = 0,
    MIDPOINT = 1,
    NODEBUG,
} dft_debug_t;

typedef enum
{
    NO_LOCK,
    SYMBOL_LOCK,
    PHASE_LOCK,
} timing_lock_t;


typedef struct
{
    FILE* waveformPlotStdin;
    FILE* fftDebuggerStdin;
    FILE* errorPlotStdin;
    FILE* IQplotStdin;
    FILE* IQvstimeStdin;
    FILE* eyeDiagramRealStdin;
    FILE* eyeDiagramImaginaryStdin;
    FILE* filterDebugStdin;
    FILE* QAMdecoderStdin;
    FILE* gardnerAlgoStdin;
    FILE* OFDMtimingSyncStdin;
    FILE* channelFilterStdin;
    FILE* OFDMdecoderStdin;
    FILE* OFDMIQStdin;
    union
    {
        struct
        {
            unsigned int waveformEnabled            : 1;
            unsigned int fftDebugEnabled            : 1;
            unsigned int errorPlotEnabled           : 1;
            unsigned int IQplotEnabled              : 1;
            unsigned int IQvsTimeEnabled            : 1;
            unsigned int eyeDiagramRealEnabled      : 1;
            unsigned int eyeDiagramImaginaryEnabled : 1;
            unsigned int filterDebugEnabled         : 1;
            unsigned int QAMdecoderEnabled          : 1;
            unsigned int gardnerAlgoEnabled         : 1;
            unsigned int OFDMtimingSyncEnabled      : 1;
            unsigned int channelFilterEnabled       : 1;
            unsigned int OFDMdecoderEnabled         : 1;
            unsigned int OFDMIQEnabled              : 1;
        };
            unsigned long int flags;
    };
} debugPlots_t;

typedef struct
{
    double carrierFrequency;
    double carrierPhase;            // probably not going to use in the end.
    double complex IQsamplingTransform;         // used to transform the IQ samples to compensate for carrier phase and amplitude mismatches (ln(amp) + i*phase)
    double k;  // number of carrier cycles per sample
    double symbolPeriod;    // number of audio samples per symbol
    double selectedIQsamples[4];    // the closest samples to ideal sample time and mid symbol sample time for gardner algorithm
} QAM_properties_t;

typedef struct
{
    double *buffer;        // a pointer to an array for storing samples
    int length;             // the length of the buffer
    int insertionIndex;     // the index where new data is inserted into the circular buffer
    int phase;              // used for some functions choosing when to do periodic calculations with the buffer
    int sampleRate;
    int n;
} circular_buffer_double_t;
typedef struct
{
    double complex *buffer;// a pointer to an array for storing samples
    int length;             // the length of the buffer
    int insertionIndex;     // the index where new data is inserted into the circular buffer
    int phase;              // used for some functions choosing when to do periodic calculations with the buffer
    int sampleRate;
    int n;                  // index of sample at insertionIndex
} circular_buffer_complex_t;

typedef struct
{
    int sampleRate;
    int channels;

    int baseband;   // if 1 than we're working with real samples, not complex

    int guardPeriod;
    int ofdmPeriod;
    int symbolPeriod;   // sum of guard and ofdm periods

    int initialized;

    circular_buffer_double_t preambleDetectorInputBuffer; // holds the last ofdmPeriod plus a little of samples for preamble detection

    circular_buffer_double_t fftwSymbolBuffer;              // holds the incoming ofdmPeriod number of samples
    circular_buffer_complex_t fftwIQbuffer;             // the decoded symbol
    fftw_plan fftwPlan;
    int ofdmPhaseOffset;                        // stores the offset in sample number to start taking cutting out OFDM symbols

    int simulateChannel;    // 0 don't simulate, 1 simulate
    circular_buffer_double_t channelSimulationBuffer;   // holds samples to be used for simulating channel impulse response
    circular_buffer_double_t channelImpulseResponse;    // holds the impulse response of the channel.

    sample_double_t autoCorrelation;        // schmidl timing signal
    sample_double_t autoCorrelationDerivative;
    circular_buffer_double_t autoCorrelationAverageBuffer;
    circular_buffer_double_t autoCorrelationDerivativeAverageBuffer;           // avarage for the derivative thing
    
    circular_buffer_double_t timingFilterInputBuffer;   // buffer of averaged autocorrelation values for timing filter

    sample_double_t preambleEnergy; // for normalizing the timing signal and automatic gain control

    struct
    {
        int long frameStart;    // estimate of the begining of the current frame of given type
        enum
        {
            SEARCHING,  // searching for the beginning of a frame, autocorrelation search
            ACTIVE,     // found a frame, currently decoding it
        } frame;

        int long fieldStart;    // estimate of teh beginning of the current field of given type
        enum
        {
            PREAMBLE,
            DATA,
        } field;

        int long symbolStart;   // offset of the beginning of the current symbol state
        enum
        {
            GUARD_PERIOD,
            OFDM_PERIOD,
        } symbol;
    } state;

} OFDM_properties_t;

typedef enum
{
    RETURNED_SAMPLE,
    AWAITING_SAMPLES,
} buffered_data_return_t;

// calculate the discrete fourier transform of an array of real values but only at frequency 'k' (k cycles per windowSize samples)
//  debugFlag is to print the right debug info for different situations
//  offset is the starting position of the window in the buffer
//  windowPhase is the offset in time the window is from the start of the audio samples, probably modulo the symbol period
//  carrierPhase is the phase offset of the carrier relative to the windowPhase
double complex dft(double* buffer,
                   int windowSize,
                   int offset,
                   int windowPhase,
                   double carrierPhase,
                   int k,
                   double* rmsOut,
                   dft_debug_t debugFlag,
                   FILE* debug_fd,
                   int debug_n)
{
#if DEBUG_LEVEL <= 1
    (void)debugFlag;
    (void)debug_fd;
    (void)debug_n;
#endif

    // compute DFT (rn just one freq (k), but when we implement OFDM, then at many harmonic(orthogonal) frequencies)
    // here is a simple quadrature detector
    double complex IQ = 0;
    double RMS = 0;

    int bufferIndex;
    double phase;
    int n;  // the sample number since start of audio recording modulo window size

    for(int i = 0; i < windowSize; i++)
    {
        // recovering the time of each sample relative to the real time modulo symbol period.
        // this is important for having a consistant carrier phase between DFTs
        n = i + windowPhase;
        // This does mean there will be a constant phase misalignment between the original carrier and the IQ demod exponential, so IQ will be rotated.
        // That will be corrected for by the carrierPhase parameter by coasta's loop

        // starts at buffer[offset] and wraps around to the beginning of the buffer
        bufferIndex = (i + offset) % windowSize;

        // phase of the complex exponential
        // phasor offsets the cos and sin waves so that their phase is alligned with the real time of the samples, offset by the carrierPhase
        phase = (double)(n) * k / windowSize + carrierPhase;

        // use a summation over the symbol period to separate orthogonal components (quadrature components) at the single frequency
        double complex wave = cexp(I*2*M_PI*phase); // generate a sample from the complex exponential
        double complex value = buffer[bufferIndex] * wave;  // multiply the complex exponential by the input sample
        IQ += value;    // integrate the result over the window

        // compute RMS amplitude for equalization -- this kinda sucks. if there is a DC bias (ie, some low frequency interference)
        // it also comes through on RMS. Gotta fix that
        RMS += pow(buffer[i], 2);


        // this debug define simplifies the function a bit if debugging is disabled
    #if DEBUG_LEVEL > 1
        switch(debugFlag)
        {
            case MIDPOINT:
            {
                // debug graph outputs
                // debugging the phase allignment
                fprintf(debug_fd, "%i %i %f %i %f %i %f\n", debug_n + i, 5, buffer[bufferIndex] + 4, 6, creal(wave) + 4, 7, cimag(wave) + 4);
                // debugging the integral
                //fprintf(debug_fd, "%i %i %f %i %f %i %f %i %f\n", debug_n + i, 11, creal(value) + 6, 12, cimag(value) + 6, 13, creal(IQ) + 6, 14, cimag(IQ) + 6);
                break;
            }

            case ALLIGNED:
            {
                // debug graph outputs
                // debugging the phase allignment
                fprintf(debug_fd, "%i %i %f %i %f %i %f\n", debug_n + i, 0, buffer[bufferIndex], 1, creal(wave), 2, cimag(wave));
                break;
            }
            case NODEBUG:
                break;
        }
    #endif
    }
    // normalization factor (do I need to divide by k?)
    IQ *= sqrt(1. / windowSize);

    // complete the RMS fomula
    RMS = sqrt(1./windowSize * RMS);

#if DEBUG_LEVEL > 1
    switch(debugFlag)
    {
        case MIDPOINT:
        {
            // debug fft plot
            fprintf(debug_fd, "%i %i %f %i %f\n", debug_n, 8, creal(IQ) + 4, 9, cimag(IQ) + 4);
            fprintf(debug_fd, "%i %i %f %i %f\n", debug_n + windowSize, 8, creal(IQ) + 4, 9, cimag(IQ) + 4);
            break;
        }

        case ALLIGNED:
        {
            // debug fft plot
            fprintf(debug_fd, "%i %i %f %i %f\n", debug_n, 3, creal(IQ), 4, cimag(IQ));
            fprintf(debug_fd, "%i %i %f %i %f\n", debug_n + windowSize, 3, creal(IQ), 4, cimag(IQ));
            break;
        }
        case NODEBUG:
            break;
    }
#endif

    *rmsOut = RMS;
    return IQ;
}

buffered_data_return_t channelFilter(const circular_buffer_double_t *inputSamples, sample_double_t *outputSample, circular_buffer_double_t *impulseResponse, int justSimulate, debugPlots_t debugPlots)
{
    // generate the filter from the inpulseResonse of the channel, hopefully rarely, let's just say once for now?
    static double *filter = NULL;
    static circular_buffer_complex_t frequencyResponse = {0};
    static circular_buffer_complex_t reciprocalFrequencyResponse = {0};
    static circular_buffer_complex_t inverseImpulseResponse = {0};
    static circular_buffer_double_t  channelSimulationBuffer = {0};

    static int initialized = 0; // have we generated the filter taps
    static int bufferPrimed = 0;    // have we waited for enough samples to start filtering
    if(!initialized)
    {
        // initialize impulse response for testing
        // pull the samples from a file called impulseResponse.raw
        int impulseResponseFile = open("data/impulseResponse.raw", O_RDONLY);
        if(impulseResponseFile == -1)
            fprintf(stderr, "unable to open data/impulseResponse.raw file for channel equalization filter.\n");
        // file format:
        //  S32_LE PCM mono samples
        //  so 4 bytes per sample
        for(int i = 0; i < impulseResponse->length; i++)
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
            int readBytes;
            if((readBytes = read(impulseResponseFile, &readSample.bytes, sizeof(readSample.bytes))) == 0)
                fprintf(stderr, "Reached end of impulseResponse.raw before filter was satisfied!\n");

            impulseResponse->buffer[i] = (double)readSample.value / INT32_MAX;    // convert to a double between -1 and 1
            //impulseResponse->buffer[i] = cos(2*M_PI * 100 * (double)i / impulseResponse->length) + cos(2*M_PI * 2 * (double)i / impulseResponse->length);
            //impulseResponse->buffer[i] = cos(2*M_PI * 100 * (double)i / impulseResponse->length);
            //impulseResponse->buffer[i] = exp((double)-i/100);
            //impulseResponse->buffer[i] = exp(-pow((double)i/100 - 5, 2));
            //impulseResponse->buffer[i] = 0;
            //if(i % (impulseResponse->length/10) == 0)
            //if(i == 0)
                //impulseResponse->buffer[i] = 1;
            //impulseResponse->insertionIndex++;
        }
        close(impulseResponseFile);

        channelSimulationBuffer.length = impulseResponse->length;
        if((channelSimulationBuffer.buffer = calloc(channelSimulationBuffer.length, sizeof(double))) == NULL)
            fprintf(stderr, "Couldn't allocate memory for channel simulation buffer: %s\n", strerror(errno));

        if(!justSimulate)
        {
            //  fourier transform the impulseResponse

            frequencyResponse.length = impulseResponse->length;
            if((frequencyResponse.buffer = calloc(frequencyResponse.length, sizeof(double complex))) == NULL)
                fprintf(stderr, "Couldn't allocate memory for frequency response buffer: %s\n", strerror(errno));
            reciprocalFrequencyResponse.length = impulseResponse->length;
            if((reciprocalFrequencyResponse.buffer = calloc(reciprocalFrequencyResponse.length, sizeof(double complex))) == NULL)
                fprintf(stderr, "Couldn't allocate memory for reciprocal frequency response buffer: %s\n", strerror(errno));
            inverseImpulseResponse.length = impulseResponse->length;
            if((inverseImpulseResponse.buffer = calloc(inverseImpulseResponse.length, sizeof(double complex))) == NULL)
                fprintf(stderr, "Couldn't allocate memory for inverse impulse response buffer: %s\n", strerror(errno));



            int N = frequencyResponse.length;
            for(int k = 0; k < N; k++)
            {
                for(int x =  0; x < N; x++)
                {
                    // dft
                    frequencyResponse.buffer[k] += cexp(-I*2*M_PI * x/N * k) * impulseResponse->buffer[x];
                }
                //  take the reciprical
                reciprocalFrequencyResponse.buffer[k] = 1. / frequencyResponse.buffer[k];
                // filter out the high frequencies? or set max value for the gain of any frequency
                //reciprocalFrequencyResponse.buffer[k] = fmax(fmin(creal(reciprocalFrequencyResponse.buffer[k]), 30), -30) + I*fmax(fmin(cimag(reciprocalFrequencyResponse.buffer[k]), 30), -30);
                // chop down the extreme values
                //if(cabs(reciprocalFrequencyResponse.buffer[k]) > 20)
                    //reciprocalFrequencyResponse.buffer[k] = 0;

                // chop off higher frequencies
                //if(abs(k - reciprocalFrequencyResponse.length / 2) < reciprocalFrequencyResponse.length / 30)
                    //reciprocalFrequencyResponse.buffer[k] = 0;
                //reciprocalFrequencyResponse.buffer[k] = frequencyResponse.buffer[k];

                //  reverse the fourier transform
                for(int x =  0; x < N; x++)
                {
                    // inverse dft (negative phase, and normalization factor)
                    inverseImpulseResponse.buffer[x] += cexp(I*2*M_PI * x/N * k) * reciprocalFrequencyResponse.buffer[k] / N;
                }
            }
        }



        initialized = 1;
    }
    if(debugPlots.channelFilterEnabled)
    {
        fprintf(debugPlots.channelFilterStdin, "%i %i %f\n",
                inputSamples->n,
                0, inputSamples->buffer[inputSamples->insertionIndex]
            );
        if(!justSimulate && inputSamples->n < impulseResponse->length)
        {
            //fprintf(debugPlots.channelFilterStdin, "%i %i %f %i %f %i %f %i %f %i %f %i %f %i %f\n",
            fprintf(debugPlots.channelFilterStdin, "%i %i %f %i %f %i %f %i %f %i %f %i %f\n",
            //fprintf(debugPlots.channelFilterStdin, "%i %i %f %i %f\n",
                    inputSamples->n,
                    1, -3 + impulseResponse->buffer[inputSamples->n % impulseResponse->length],
                    3, 13 + cabs(frequencyResponse.buffer[(inputSamples->n + frequencyResponse.length / 2) % frequencyResponse.length]),
                    4, 7 + carg(frequencyResponse.buffer[(inputSamples->n + frequencyResponse.length / 2) % frequencyResponse.length]),
                    5, 13 + cabs(reciprocalFrequencyResponse.buffer[(inputSamples->n + reciprocalFrequencyResponse.length / 2) % reciprocalFrequencyResponse.length]),
                    6, 7 + carg(reciprocalFrequencyResponse.buffer[(inputSamples->n + reciprocalFrequencyResponse.length / 2) % reciprocalFrequencyResponse.length]),
                    //7, creal(inverseImpulseResponse.buffer[(inputSamples->n + inverseImpulseResponse.length / 2) % inverseImpulseResponse.length]),
                    //8, cimag(inverseImpulseResponse.buffer[(inputSamples->n + inverseImpulseResponse.length / 2) % inverseImpulseResponse.length])
                    //3, cabs(frequencyResponse.buffer[inputSamples->n % frequencyResponse.length]),
                    //4, carg(frequencyResponse.buffer[inputSamples->n % frequencyResponse.length]),
                    //5, cabs(reciprocalFrequencyResponse.buffer[inputSamples->n % reciprocalFrequencyResponse.length]),
                    //6, carg(reciprocalFrequencyResponse.buffer[inputSamples->n % reciprocalFrequencyResponse.length]),
                    7, -12 + creal(inverseImpulseResponse.buffer[inputSamples->n % inverseImpulseResponse.length])
                    //8, cimag(inverseImpulseResponse.buffer[inputSamples->n % inverseImpulseResponse.length])
                );
        } else {
            fprintf(debugPlots.channelFilterStdin, "%i %i %f\n",
                    inputSamples->n,
                    1, -3 + impulseResponse->buffer[inputSamples->n % impulseResponse->length]
                );
        }
    }

    // wait for enough samples to come in to do the first convolution
    // I don't have to wait, I'm pulling in past samples, and before the first sample, they can be zeros
    // at least for the test where I convolve the impulse response with the input

    //if(inputSamples->n < impulseResponse->length - 1)
    //if(inputSamples->n < 0)
        //return AWAITING_SAMPLES;

    // just do a test, let's take the impulse response and convolve it with the input samples to simulate the channel
    channelSimulationBuffer.buffer[channelSimulationBuffer.insertionIndex] = 0;
    //outputSample->sample = 0;
    for(int i = 0; i < impulseResponse->length; i++)
    {
        // multiply the filter by the input samples at respective times and sum the results
        // index, start at most recent input sample and move backwards in time
        int inputIndex = (inputSamples->insertionIndex - i) % inputSamples->length;
        if(inputIndex < 0)
            inputIndex += inputSamples->length;
        int filterIndex = i;

        //outputSample->sample += inputSamples->buffer[inputIndex] * impulseResponse->buffer[filterIndex] * -inverseImpulseResponse.buffer[filterIndex];
        //outputSample->sample += inputSamples->buffer[inputIndex] * impulseResponse->buffer[filterIndex];
        channelSimulationBuffer.buffer[channelSimulationBuffer.insertionIndex] += inputSamples->buffer[inputIndex] * impulseResponse->buffer[filterIndex];
        outputSample->sample = channelSimulationBuffer.buffer[channelSimulationBuffer.insertionIndex];
        //outputSample->sample += inputSamples->buffer[inputIndex] * inverseImpulseResponse.buffer[filterIndex];
    }
    channelSimulationBuffer.n = inputSamples->n;
    outputSample->sampleIndex = inputSamples->n;
    channelSimulationBuffer.sampleRate = inputSamples->sampleRate;
    outputSample->sampleRate = inputSamples->sampleRate;

    // increment channel sim buffer index now that we're done accessing it
    channelSimulationBuffer.insertionIndex = (channelSimulationBuffer.insertionIndex + 1) % channelSimulationBuffer.length;

    // if just simulating channel, return the output now
    if(justSimulate)
    {
        return RETURNED_SAMPLE;
    }

    // now reverse the channel by convolving the inverse filter

    // wait for enough sampls to begin convolution
    if(channelSimulationBuffer.n < channelSimulationBuffer.length - 1)
        return AWAITING_SAMPLES;
    // apply inverse filter
    outputSample->sample = 0;
    for(int i = 0; i < channelSimulationBuffer.length; i++)
    {
        // index, start at oldest sample, and move forward in time up to most recent
        //int inputIndex = (channelSimulationBuffer.insertionIndex + 1 + i) % channelSimulationBuffer.length;
        // start at newest sample and move backwards in time to oldest sample
        int inputIndex = (channelSimulationBuffer.insertionIndex  - 1- i) % channelSimulationBuffer.length;
        if(inputIndex < 0)
            inputIndex += channelSimulationBuffer.length;
        // reverse the filter direction, start at filter end, and work to beginning
        //int filterIndex = inverseImpulseResponse.length - 1 - i;
        int filterIndex = i;

        outputSample->sample += channelSimulationBuffer.buffer[inputIndex] * creal(inverseImpulseResponse.buffer[filterIndex]);
        //outputSample->sample += channelSimulationBuffer.buffer[i];
    }
    //outputSample->sample = channelSimulationBuffer.buffer[channelSimulationBuffer.insertionIndex+1];
    outputSample->sampleIndex = channelSimulationBuffer.n - channelSimulationBuffer.length;
    outputSample->sampleRate = channelSimulationBuffer.sampleRate;



    //outputSample->sampleIndex = inputSamples->n - impulseResponse->length;
    //outputSample->sampleIndex = inputSamples->n;
    //outputSample->sampleIndex = inputSamples->n - 0;

    // test, bypass
    //outputSample->sample = inputSamples->buffer[inputSamples->insertionIndex];


    if(debugPlots.channelFilterEnabled)
    {
        fprintf(debugPlots.channelFilterStdin, "%i %i %f %i %f\n",
                outputSample->sampleIndex,
                //inputSamples->n,
                2, -6 + channelSimulationBuffer.buffer[channelSimulationBuffer.insertionIndex],
                9, -18 + outputSample->sample
            );

    }
    // collect enough input samples ahead of the output timepoint to start filtering
    // do one step of convolution to get one output sample
    return RETURNED_SAMPLE;
}

buffered_data_return_t raisedCosFilter(const circular_buffer_complex_t *inputSamples, sample_complex_t *outputSample, double cutoffFrequency, debugPlots_t debugPlots)
{
    // array to store time series of filter data
    static double *filter = NULL;

    if (inputSamples == NULL)
    {
        if (filter != NULL)
        {
            free(filter);
            filter = NULL;
        }
    }

    //int symbolPeriod = 64; // audio samples per symbol
    int k = 1; // cycles per period
    //int filterSides = 10;    // number of symbols to either side of current symbol to filter with raised cos
    //int filterLengthSymbols = 2 * filterSides + 1;    // length of raised cos filter in IQ symbols, ie, how many IQ samples we need to generate the current symbol
    //int filterLength = filterLengthSymbols * symbolPeriod;  // length in audio samples

    //double filterCutoffFrequency = cutoffFrequency * 1.25; // cutoff frequency of the low pass filter
    double filterCutoffFrequency = cutoffFrequency; // cutoff frequency of the low pass filter
    //double filterCutoffFrequency = 500. * 1.25; // cutoff frequency of the low pass filter

    // array to store timeseries of IQ samples
    //static double complex *IQdata;

    static int initialized = 0;         // initialize the filter kernel first, and only once
    static int bufferPrimed = 0;        // get enough samples to start processing them.
    if(!initialized)
    {
        // initialize raised cos filter data
        //filter = malloc(sizeof(double) * filterLength); // this never gets released. so, might wanna fix that TODO
        filter = malloc(sizeof(double) * inputSamples->length);  // filter kernel = size of input buffer

        for(int i = 0; i < inputSamples->length; i++)
        {
            int filterIndex = i - inputSamples->length / 2;    // should go -filterLength/2 -> 0 -> filterLength/2
            // raised cos filter math
            //double b = 0.42;    // filter parameter beta, has to do with frequency falloff and time domain fall off
            double b = 0.3;    // filter parameter beta, has to do with frequency falloff and time domain fall off

            double filterValue =
            (
                sin(2*M_PI * filterCutoffFrequency * filterIndex / inputSamples->sampleRate) /
                (2*M_PI * filterCutoffFrequency * filterIndex / inputSamples->sampleRate) *
                (cos(2*M_PI * filterCutoffFrequency * filterIndex / inputSamples->sampleRate * b)) /
                (1 - pow(4 * b * filterCutoffFrequency * filterIndex / inputSamples->sampleRate, 2))
            );

            if(!isfinite(filterValue))   // in case it's undefined, ie divide by zero case
                filterValue = sin(M_PI / 2 / b) / (M_PI / 2 / b);

            if(filterIndex == 0)
                filterValue = 1;    // the math gives a divide by zero at index 0

            filter[i] = filterValue;
        }

        // ensure initialize only runs once
        initialized = 1;
    }
    // wait for enough samples to come in
    if(!bufferPrimed)
    {
        // basically, we need to wait for future samples before we can filter the first sample. skip the remainder of this function until enough samples are collected.
        if(inputSamples->insertionIndex < inputSamples->length / 2)    // filled the buffer half way, now ready to start processing
            return AWAITING_SAMPLES;    // wait for samples

        bufferPrimed = 1;   // otherwise stop waiting and process
    }

    //int relativeSampleIndex = inputSamples.insertionIndex - inputSamples.length / 2;    // 0 is the 'current' time in the circular input buffer. negatives are into the past, positives are into the future
    //  insertionIndex should the the index of the last inserted sample, inserted by the calling function
    outputSample->sample = 0;
    outputSample->sampleRate = inputSamples->sampleRate;
    outputSample->sampleIndex = inputSamples->n - inputSamples->length / 2;   // the index is shifted by half the filter width
    for(int i = 0; i < inputSamples->length; i++)
    {
        // calculate the convoution for the sample at relative index 0

        // generate relative indexes
        int relativeSampleIndex = i - inputSamples->length / 2;     // 0 is the 'current' time in the circular input buffer. negatives are into the past, positives are into the future
        int relativeFilterIndex = -relativeSampleIndex;     // 0 is the center of the filter

        // generate absolute index positions
        int sampleIndex = (relativeSampleIndex + inputSamples->insertionIndex - inputSamples->length / 2) % inputSamples->length;

        if(sampleIndex < 0)
            sampleIndex += inputSamples->length; // wrap negative indexes back around to positive values

        int filterIndex = (relativeFilterIndex + inputSamples->length / 2) % inputSamples->length;

        if(filterIndex < 0)
            filterIndex += inputSamples->length; // wrap negative indexes back around to positive values

        // calculate the multiplication and sumation
        outputSample->sample += inputSamples->buffer[sampleIndex] * filter[filterIndex];
        if(debugPlots.filterDebugEnabled && relativeSampleIndex == 0)
        {
            fprintf(
                debugPlots.filterDebugStdin,
                "%i %i %f %i %f %i %f\n",
                outputSample->sampleIndex,
                1,
                creal(inputSamples->buffer[sampleIndex]),
                2,
                cimag(inputSamples->buffer[sampleIndex]),
                3,
                filter[outputSample->sampleIndex%inputSamples->length]
            );
        }
    }
    //printf("output: n=%i, %f+%fi\n", inputSamples.n, creal(outputSample->sample), cimag(outputSample->sample));
    //outputSample->sample /= 24.;        // bit arbitrary, gotta figure out normalization factor for raised cos filter

    // print out debug info for the filter
    //fprintf(debugPlots.filterDebugStdin, "%i %i %f %i %f %i %f\n", inputSamples.n, 1, creal(inputSamples.buffer[inputSamples.insertionIndex]), 2, cimag(inputSamples.buffer[inputSamples.insertionIndex]), 3, filter[outputSample->sampleIndex%inputSamples.length]);
    if(debugPlots.filterDebugEnabled)
        fprintf(debugPlots.filterDebugStdin, "%i %i %f %i %f\n", outputSample->sampleIndex, 6, creal(outputSample->sample), 7, cimag(outputSample->sample));

    return RETURNED_SAMPLE;
}

buffered_data_return_t averagingFilter(const circular_buffer_double_t *inputSamples, sample_double_t *outputSample, double *runningAverage)
{
    /*
    double average = 0;
    // convlutional method
    for(int i  = 0; i < inputSamples->length; i++)
    {
        average += inputSamples->buffer[i];
    }
    average /= inputSamples->length;
    outputSample->sampleIndex = inputSamples->n - inputSamples->length / 2; // center the filter so it averages before and after in a sliding window
    outputSample->sampleRate = inputSamples->sampleRate;
    outputSample->sample = average;
    return RETURNED_SAMPLE;
    */

    // iterative method
    //static double iterativeAverage = 0;
    int lateIndex = inputSamples->insertionIndex;   // newest index
    int earlyIndex = (inputSamples->insertionIndex + 1) % inputSamples->length;  // oldest index
    *runningAverage += (inputSamples->buffer[lateIndex] - inputSamples->buffer[earlyIndex]) / inputSamples->length;

    outputSample->sampleIndex = inputSamples->n - inputSamples->length / 2; // center the filter so it averages before and after in a sliding window
    outputSample->sampleRate = inputSamples->sampleRate;
    outputSample->sample = *runningAverage;

    return RETURNED_SAMPLE;
}

// detect the most likely/optimal for the channel conditions ofdm frame start time
buffered_data_return_t timingPeakDetectionFilter(const circular_buffer_double_t *inputSamples, sample_double_t *outputSample, debugPlots_t debugPlots)
{
    // width of input window should be some 4 times larger than guard interval
    sample_double_t minimumFirstHalf = {0};
    minimumFirstHalf.sample = inputSamples->buffer[(inputSamples->insertionIndex + 1) % inputSamples->length];    // equals value of first sample
    minimumFirstHalf.sampleIndex = inputSamples->n - inputSamples->length + 1;
    minimumFirstHalf.sampleRate = inputSamples->sampleRate;

    sample_double_t minimumSecondHalf = {0};
    minimumSecondHalf.sample = inputSamples->buffer[inputSamples->insertionIndex];    // equals value of last sample
    minimumSecondHalf.sampleIndex = inputSamples->n;
    minimumSecondHalf.sampleRate = inputSamples->sampleRate;

    sample_double_t maximum = {0};
    maximum.sample = inputSamples->buffer[(inputSamples->insertionIndex + 1 + inputSamples->length / 2) % inputSamples->length];    // equals value right in the middle
    maximum.sampleIndex = inputSamples->n - inputSamples->length + inputSamples->length / 2;
    maximum.sampleRate = inputSamples->sampleRate;
    /*
    // width of input window should be some 4 times larger than guard interval
    sample_double_t minimumFirstHalf = {0};
    //minimumFirstHalf.sample = inputSamples->buffer[inputSamples->insertionIndex + 1];    // equals value of first sample
    minimumFirstHalf.sample = 10;    // equals value of first sample
    //minimumFirstHalf.sampleIndex = inputSamples->n - inputSamples->length + 1;
    minimumFirstHalf.sampleIndex = -1;
    minimumFirstHalf.sampleRate = inputSamples->sampleRate;

    sample_double_t minimumSecondHalf = {0};
    //minimumSecondHalf.sample = inputSamples->buffer[inputSamples->insertionIndex + 1];    // equals value of first sample
    minimumSecondHalf.sample = 10;    // equals value of first sample
    //minimumSecondHalf.sampleIndex = inputSamples->n - inputSamples->length + 1;
    minimumSecondHalf.sampleIndex = -1;
    minimumSecondHalf.sampleRate = inputSamples->sampleRate;

    sample_double_t maximum = {0};
    //maximum.sample = inputSamples->buffer[inputSamples->insertionIndex + 1];    // equals value of first sample
    maximum.sample = -1;    // equals value of first sample
    //maximum.sampleIndex = inputSamples->n - inputSamples->length + 1;
    maximum.sampleIndex = -1;
    maximum.sampleRate = inputSamples->sampleRate;
    */
    // if first and second half minimums are below some threshold
    // if maximum overall is above some threshold
    // if the maximum occurs between the minima
    double lowerThreshold = 0.6;
    double upperThreshold = 0.7;
    
    // only run the convolution if it's likely to be close, because it's really slow
    if(minimumFirstHalf.sample < lowerThreshold && minimumSecondHalf.sample < lowerThreshold && maximum.sample > upperThreshold)
    {
        for(int i = 0; i < inputSamples->length; i++)
        {
            double value = inputSamples->buffer[(inputSamples->insertionIndex + 1 + i) % inputSamples->length];
            double index = inputSamples->n - inputSamples->length + 1 + i;

            if(value >= maximum.sample)
            {
                // maximum overall
                maximum.sample = value;
                maximum.sampleIndex = index;
            } 
            if(i < inputSamples->length / 2 && value < minimumFirstHalf.sample)
            {
                // determine minimum of first half of samples
                minimumFirstHalf.sample = value;
                minimumFirstHalf.sampleIndex = index;
            }
            if(i >= inputSamples->length / 2 && value <= minimumSecondHalf.sample)
            {
                // minimum of last half of samples
                minimumSecondHalf.sample = value;
                minimumSecondHalf.sampleIndex = index;
            }
        }
        
        outputSample->sample = 0;
        outputSample->sampleIndex = 0;
        outputSample->sampleRate = inputSamples->sampleRate;
        if(minimumSecondHalf.sample < lowerThreshold && minimumFirstHalf.sample < lowerThreshold && maximum.sample > upperThreshold && minimumFirstHalf.sampleIndex < maximum.sampleIndex&& minimumSecondHalf.sampleIndex > maximum.sampleIndex)
        {
            // the offset is the index value of the maximum
            outputSample->sample = 1;
            outputSample->sampleIndex = maximum.sampleIndex;    // index is the offset of maximum
            return RETURNED_SAMPLE; // only return a sample when preamble is detected
        }
    }
    return AWAITING_SAMPLES;
}

buffered_data_return_t gardnerAlgorithm(const circular_buffer_complex_t *inputSamples, sample_complex_t *outputSamples, double symbolPeriod, debugPlots_t debugPlots)
{
    // now do some timing alignemnt.
    //      I'm imagining a clock that counts up by a phase quantity scaled by the number of IQ samples elapsed. accumulates phase
    //      The rate of phase change will be adjusted by a PLL
    //      The zero crossing for the clock will be projected so that the two samples it falls between can be chosen.
    //      I'll keep an array of past IQ samples to choose from.
    //      The next filteredIQsample.sampleIndex where calculations will occur is calculated by the projection
    //          concept after PLL updates during said calculations

    // Choosing samples for timing lock and symbol detection
    static double symbolSamplerAccumulatedPhase = 0;
    static double symbolSamplerPhaseRate = 0;
    static int symbolSamplerNextIndex = 0;  // next index to trigger calculations, should be just after the ideal sample time.
    static int idealIQindex = 0;

    if(debugPlots.gardnerAlgoEnabled)
    {
        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %f\n", inputSamples->n, 5, creal(inputSamples->buffer[inputSamples->insertionIndex]));
        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %f\n", inputSamples->n, 6, cimag(inputSamples->buffer[inputSamples->insertionIndex]));
    }

    if(inputSamples->n < symbolSamplerNextIndex) // check if it's too early
        return AWAITING_SAMPLES;    // it's too early, wait till the right number of IQ samples has passed.


    if(symbolSamplerAccumulatedPhase == 0)
    {
        // initialize phase rate
        symbolSamplerPhaseRate = symbolPeriod; // initialize the phase rate to the idealized value
        symbolSamplerAccumulatedPhase = symbolSamplerPhaseRate + 0.25;    // trigger calculations one period from now
        symbolSamplerNextIndex = (int)ceil(symbolSamplerAccumulatedPhase);    // trigger calculations one period from now
        return AWAITING_SAMPLES;
    }

    // Gardner Algorithm
    //  this can happen just once per IQ symbol, so not every audio sample
    int postIdealIndexOffset = 0;
    int preIdealIndexOffset = -1;
    int postMidIndexOffset = postIdealIndexOffset - symbolPeriod / 2;
    int preMidIndexOffset = postMidIndexOffset - 1;

    int postIdealIndex = inputSamples->insertionIndex;
    int preIdealIndex = postIdealIndex + preIdealIndexOffset;
    int postMidIndex =  postIdealIndex + postMidIndexOffset;
    int preMidIndex =   postIdealIndex+ preMidIndexOffset;


    postIdealIndex =    postIdealIndex < 0 ? inputSamples->length + postIdealIndex : postIdealIndex;     // wrap to positive
    preIdealIndex =     preIdealIndex < 0 ? inputSamples->length + preIdealIndex : preIdealIndex;        // wrap
    postMidIndex =      postMidIndex < 0 ? inputSamples->length + postMidIndex : postMidIndex;           // wrap
    preMidIndex =       preMidIndex < 0 ? inputSamples->length + preMidIndex : preMidIndex;              // wrap

    //  Interpolation between samples
    double complex IQmidpoint = (inputSamples->buffer[preMidIndex] - inputSamples->buffer[postMidIndex]) * (symbolSamplerNextIndex - symbolSamplerAccumulatedPhase) + inputSamples->buffer[postMidIndex];
    static double complex IQideal = 0;
    double complex IQlast = IQideal;
    IQideal = (inputSamples->buffer[preIdealIndex] - inputSamples->buffer[postIdealIndex]) * (symbolSamplerNextIndex - symbolSamplerAccumulatedPhase) + inputSamples->buffer[postIdealIndex];

    // calculate error signal
    //      I either need to changed the phaseErrorEstimate to rateOfPhaseErrorEstimate, or else allow the PID loop to set the phase offset directly rather than adjusting it's derivitive
    static double symbolSamplerPhaseErrorEstimate = 0;
    static double phaseErrorEstimateArray[10];
    static int phaseErrorEstimateArraySize = 10;
    static int phaseErrorEstimateArrayIndex = 0;

    // update with new values
    double errorLast = symbolSamplerPhaseErrorEstimate;
    symbolSamplerPhaseErrorEstimate = creal((IQideal - IQlast) * conj(IQmidpoint));  // this should get us a rough sync

    // take a rolling average of the error estimate
    phaseErrorEstimateArray[phaseErrorEstimateArrayIndex] = symbolSamplerPhaseErrorEstimate;
    phaseErrorEstimateArrayIndex = (phaseErrorEstimateArrayIndex + 1) % phaseErrorEstimateArraySize;

    static double phaseErrorEstimateAverage = 0;
    double phaseErrorEstimateAverageLast = phaseErrorEstimateAverage;
    phaseErrorEstimateAverage = 0;
    for(int i = 0; i < phaseErrorEstimateArraySize; i++)
    {
        phaseErrorEstimateAverage += phaseErrorEstimateArray[i];
    }
    phaseErrorEstimateAverage /= phaseErrorEstimateArraySize;


    //double error_I = symbolSamplerPhaseErrorEstimate;  // integral is just the value before derivative
    //double error = symbolSamplerPhaseErrorEstimate - errorLast; // error is derivative of phase offset estimate
    double error_I = phaseErrorEstimateAverage;  // integral is just the value before derivative
    double error = phaseErrorEstimateAverage - phaseErrorEstimateAverageLast; // error is derivative of phase offset estimate

    // this is shitty code, should be a function
    // PID loop
    //symbolSamplerPhaseRate -= error_I * 0.6 + error * 1.5;  // tight for pure alternating I
    //symbolSamplerPhaseRate -= error_I * 0.2 + error * 0.5;  // tight for pure alternating I
    symbolSamplerPhaseRate -= error_I * 0.002 + error * 0.06;  // tight for pure alternating I
    symbolSamplerPhaseRate = fmin(fmax(symbolSamplerPhaseRate, (double)symbolPeriod / 1.5), (double)symbolPeriod * 1.5);

    if(debugPlots.gardnerAlgoEnabled)
    {
        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %i\n", inputSamples->n + postIdealIndexOffset - 1, 8, 0);
        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %f\n", inputSamples->n + postIdealIndexOffset, 8, creal(inputSamples->buffer[postIdealIndex]));
        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %i\n", inputSamples->n + postIdealIndexOffset + 1, 8, 0);

        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %i\n", inputSamples->n + postIdealIndexOffset - 1, 9, 0);
        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %f\n", inputSamples->n + postIdealIndexOffset, 9, cimag(inputSamples->buffer[postIdealIndex]));
        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %i\n", inputSamples->n + postIdealIndexOffset + 1, 9, 0);

        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %i\n", inputSamples->n + preIdealIndexOffset - 1, 10, 0);
        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %f\n", inputSamples->n + preIdealIndexOffset, 10, creal(inputSamples->buffer[preIdealIndex]));
        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %i\n", inputSamples->n + preIdealIndexOffset + 1, 10, 0);

        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %i\n", inputSamples->n + preIdealIndexOffset - 1, 11, 0);
        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %f\n", inputSamples->n + preIdealIndexOffset, 11, cimag(inputSamples->buffer[preIdealIndex]));
        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %i\n", inputSamples->n + preIdealIndexOffset + 1, 11, 0);

        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %i\n", inputSamples->n + postMidIndexOffset - 1, 12, 0);
        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %f\n", inputSamples->n + postMidIndexOffset, 12, creal(inputSamples->buffer[postMidIndex]));
        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %i\n", inputSamples->n + postMidIndexOffset + 1, 12, 0);

        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %i\n", inputSamples->n + postMidIndexOffset - 1, 13, 0);
        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %f\n", inputSamples->n + postMidIndexOffset, 13, cimag(inputSamples->buffer[postMidIndex]));
        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %i\n", inputSamples->n + postMidIndexOffset + 1, 13, 0);

        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %i\n", inputSamples->n + preMidIndexOffset - 1, 14, 0);
        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %f\n", inputSamples->n + preMidIndexOffset, 14, creal(inputSamples->buffer[preMidIndex]));
        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %i\n", inputSamples->n + preMidIndexOffset + 1, 14, 0);

        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %i\n", inputSamples->n + preMidIndexOffset - 1, 15, 0);
        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %f\n", inputSamples->n + preMidIndexOffset, 15, cimag(inputSamples->buffer[preMidIndex]));
        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %i\n", inputSamples->n + preMidIndexOffset + 1, 15, 0);

        //fprintf(debugPlots.gardnerAlgoStdin, "%f %i %f\n", symb
        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %f\n", inputSamples->n, 7, symbolSamplerPhaseErrorEstimate);
        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %f\n", inputSamples->n, 21, phaseErrorEstimateAverage);

        //fprintf(debugPlots.gardnerAlgoStdin, "%f %i %i\n", symbolSamplerAccumulatedPhase - 0.5, 16, 0);
        fprintf(debugPlots.gardnerAlgoStdin, "%f %i %f\n", symbolSamplerAccumulatedPhase, 16, creal(IQideal));
        //fprintf(debugPlots.gardnerAlgoStdin, "%f %i %i\n", symbolSamplerAccumulatedPhase + 0.5, 16, 0);

        //fprintf(debugPlots.gardnerAlgoStdin, "%f %i %i\n", symbolSamplerAccumulatedPhase - 0.5, 17, 0);
        fprintf(debugPlots.gardnerAlgoStdin, "%f %i %f\n", symbolSamplerAccumulatedPhase, 17, cimag(IQideal));
        //fprintf(debugPlots.gardnerAlgoStdin, "%f %i %i\n", symbolSamplerAccumulatedPhase + 0.5, 17, 0);

        fprintf(debugPlots.gardnerAlgoStdin, "%f %i %i\n", symbolSamplerAccumulatedPhase - symbolPeriod / 2 - 0.5, 18, 0);
        fprintf(debugPlots.gardnerAlgoStdin, "%f %i %f\n", symbolSamplerAccumulatedPhase - symbolPeriod / 2, 18, creal(IQmidpoint));
        fprintf(debugPlots.gardnerAlgoStdin, "%f %i %i\n", symbolSamplerAccumulatedPhase - symbolPeriod / 2 + 0.5, 18, 0);

        fprintf(debugPlots.gardnerAlgoStdin, "%f %i %i\n", symbolSamplerAccumulatedPhase - symbolPeriod / 2 - 0.5, 19, 0);
        fprintf(debugPlots.gardnerAlgoStdin, "%f %i %f\n", symbolSamplerAccumulatedPhase - symbolPeriod / 2, 19, cimag(IQmidpoint));
        fprintf(debugPlots.gardnerAlgoStdin, "%f %i %i\n", symbolSamplerAccumulatedPhase - symbolPeriod / 2 + 0.5, 19, 0);

        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %i\n", symbolSamplerNextIndex, 20, 0);
        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %i\n", symbolSamplerNextIndex, 20, 2);
        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %i\n", symbolSamplerNextIndex, 20, -2);
        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %i\n", symbolSamplerNextIndex, 20, 0);
    }

    if(debugPlots.IQplotEnabled)
        fprintf(debugPlots.IQplotStdin, "%f %i %f\n", creal(IQideal), symbolSamplerNextIndex % (int)symbolPeriod, cimag(IQideal));

    /*
     * I want this to work in some way to show the IQ plot with the phase offset determined by the timing algorithm
     * so that it's from the perspective of the timing plot, and we can see the effectiveness of time sync
    if(debugPlots.eyeDiagramRealEnabled)
        fprintf(debugPlots.eyeDiagramRealStdin, "%f %i %f\n",
                fmod(symbolSamplerAccumulatedPhase, symbolPeriod * 4),
                -1,
                creal(IQideal)
                );
    if(debugPlots.eyeDiagramRealEnabled)
        fprintf(debugPlots.eyeDiagramRealStdin, "%f %i %f\n",
                fmod(symbolSamplerAccumulatedPhase - symbolPeriod / 2, symbolPeriod * 4),
                -1,
                creal(IQmidpoint)
                );
    */

    // update the phase rate for symbol sampler
    symbolSamplerAccumulatedPhase += symbolSamplerPhaseRate;    // apply the phase offset
    symbolSamplerNextIndex = (int)ceil(symbolSamplerAccumulatedPhase); // determine the next index to make calculations

    outputSamples->sample = IQideal;
    outputSamples->sampleIndex = idealIQindex;
    //outputSamples->sampleRate =   // not sure actually
    idealIQindex++; // increment symbol number
    return RETURNED_SAMPLE;
}

buffered_data_return_t demodulateQAM(const sample_double_t *sample, QAM_properties_t QAMstate, debugPlots_t debugPlots)
{
    // currently does not handle samples after about half the filter size, instead just using them as future samples, but never getting to them. Should add a command to take care of the remainder I guess and cycle the filter function with zeros until it's done TODO
    static int initialized = 0;

    static circular_buffer_double_t  channelFilterBuffer = {0};
    static circular_buffer_double_t  impulseResponse = {0};
    static circular_buffer_complex_t filterInputBuffer = {0};
    static circular_buffer_complex_t timingSyncBuffer = {0};

    // deallocation of dynamic memory if the input sample is NULL
    if (sample == NULL)
    {
        if(channelFilterBuffer.buffer != NULL)
        {
            free(channelFilterBuffer.buffer);
            channelFilterBuffer.buffer = NULL;
        }
        if (filterInputBuffer.buffer != NULL)
        {
            free(filterInputBuffer.buffer);
            filterInputBuffer.buffer = NULL;
        }

        if (timingSyncBuffer.buffer != NULL)
        {
            free(timingSyncBuffer.buffer);
            timingSyncBuffer.buffer = NULL;
        }
    }

    if(!initialized)
    {
        // initialization of circular buffers

        // initialize channel equilztion filter
        channelFilterBuffer.length = 5912;
        channelFilterBuffer.insertionIndex = 0;
        channelFilterBuffer.buffer = calloc(channelFilterBuffer.length, sizeof(double));
        if(channelFilterBuffer.buffer == NULL)
            fprintf(stderr, "ChannelFIlterBuffer failed to allocate: %s\n", strerror(errno));

        // initialize impulse response
        impulseResponse.length = 5912;
        impulseResponse.insertionIndex = 0;
        impulseResponse.buffer = calloc(impulseResponse.length, sizeof(double));
        if(impulseResponse.buffer == NULL)
            fprintf(stderr, "ImpulseResponse failed to allocate: %s\n", strerror(errno));

        // initialize filter buffer
        filterInputBuffer.length = sample->sampleRate / QAMstate.carrierFrequency * 4 * (10 * 2 + 1);
        filterInputBuffer.insertionIndex = 0;
        //filterInputBuffer.buffer = malloc(sizeof(double complex) * filterInputBuffer.length);   // allocate some space for the buffer
        filterInputBuffer.buffer = calloc(filterInputBuffer.length, sizeof(double complex));   // allocate some space for the buffer, and zero it
        //printf("got %lu bytes from malloc\n", sizeof(double complex) * filterInputBuffer.length);
        if(filterInputBuffer.buffer == NULL)
            fprintf(stderr, "FilterInputBuffer failed to allocate: %s\n", strerror(errno));

        // initialize timing sync buffer
        timingSyncBuffer.length = QAMstate.symbolPeriod * 2;    // give it some extra size.
        timingSyncBuffer.insertionIndex = 0;
        timingSyncBuffer.buffer = calloc(timingSyncBuffer.length, sizeof(double complex));
        if(timingSyncBuffer.buffer == NULL)
            fprintf(stderr, "TimingSyncBuffer failed to allocate: %s\n", strerror(errno));


        // only run once
        initialized = 1;

    }


    // raw inputs
    if(debugPlots.QAMdecoderEnabled)
        fprintf(debugPlots.QAMdecoderStdin, "%i %i %f\n", sample->sampleIndex, 0, sample->sample);

    // channel equalization filter
    channelFilterBuffer.n = sample->sampleIndex;
    channelFilterBuffer.sampleRate = sample->sampleRate;
    channelFilterBuffer.phase = 0;
    channelFilterBuffer.buffer[channelFilterBuffer.insertionIndex] = sample->sample;
    sample_double_t equalizedSample;

    buffered_data_return_t returnValue = channelFilter(&channelFilterBuffer, &equalizedSample, &impulseResponse, 0, debugPlots);

    channelFilterBuffer.insertionIndex = (channelFilterBuffer.insertionIndex + 1) % channelFilterBuffer.length;

    if(returnValue != RETURNED_SAMPLE)
        return AWAITING_SAMPLES;

    if(debugPlots.QAMdecoderEnabled)
        fprintf(debugPlots.QAMdecoderStdin, "%i %i %f\n", equalizedSample.sampleIndex, 1, equalizedSample.sample);


    // IQ sampling
    //  multiply by sin and cos
    //      this step may have issues as the carrier frequency comes up to a quarter of the sample rate, since
    //      the multiplication step generates frequencies in the IQsample centered a  twice the original carrier
    //      frequency. The bandwidth of the signal on that elevated carrier may alias. I should double sample rate before this step,
    //      then reduce the sample rate to fraction of the original after filtering out that high frequency stuff.
    double complex wave = cexp(I*(2*M_PI * QAMstate.carrierFrequency * equalizedSample.sampleIndex / equalizedSample.sampleRate + QAMstate.carrierPhase));
    double complex IQsample = equalizedSample.sample * wave;

    if(debugPlots.filterDebugEnabled)
    {
        fprintf(debugPlots.filterDebugStdin, "%i %i %f\n", equalizedSample.sampleIndex, 0, equalizedSample.sample);
        fprintf(debugPlots.filterDebugStdin, "%i %i %f %i %f\n", equalizedSample.sampleIndex, 4, creal(wave), 5, cimag(wave));
    }

    if(debugPlots.gardnerAlgoEnabled)
        fprintf(debugPlots.gardnerAlgoStdin, "%i %i %f\n", equalizedSample.sampleIndex, 0, equalizedSample.sample);

    //  low pass filter
    //      convolution with a raised cos filter would probably do fine
    filterInputBuffer.n = equalizedSample.sampleIndex;
    filterInputBuffer.sampleRate = equalizedSample.sampleRate;
    filterInputBuffer.phase = 0;
    filterInputBuffer.buffer[filterInputBuffer.insertionIndex] = IQsample;
    sample_complex_t filteredIQsample;

    // filter the IQ samples
    returnValue = raisedCosFilter(&filterInputBuffer, &filteredIQsample, QAMstate.carrierFrequency * 1, debugPlots);

    filterInputBuffer.insertionIndex = (filterInputBuffer.insertionIndex + 1) % filterInputBuffer.length;

    if(returnValue != RETURNED_SAMPLE) // don't continue processing unless a sample is returned
        return AWAITING_SAMPLES;

    if(debugPlots.QAMdecoderEnabled)
        fprintf(debugPlots.QAMdecoderStdin, "%i %i %f %i %f\n", filteredIQsample.sampleIndex, 5, creal(filteredIQsample.sample), 6, cimag(filteredIQsample.sample));
    if(debugPlots.eyeDiagramRealEnabled)
        fprintf(debugPlots.eyeDiagramRealStdin, "%i %i %f\n",
                filteredIQsample.sampleIndex % ((int)QAMstate.symbolPeriod * 4), 
                filteredIQsample.sampleIndex / ((int)QAMstate.symbolPeriod * 4),
                creal(filteredIQsample.sample)
                );
    if(debugPlots.eyeDiagramImaginaryEnabled)
        fprintf(debugPlots.eyeDiagramImaginaryStdin, "%i %i %f\n",
                filteredIQsample.sampleIndex % ((int)QAMstate.symbolPeriod * 4), 
                filteredIQsample.sampleIndex / ((int)QAMstate.symbolPeriod * 4),
                cimag(filteredIQsample.sample)
                );
    if(debugPlots.IQvsTimeEnabled)
        fprintf(debugPlots.IQvstimeStdin, "%f %f\n", creal(filteredIQsample.sample), cimag(filteredIQsample.sample));

    // collect samples into a buffer
    timingSyncBuffer.buffer[timingSyncBuffer.insertionIndex] = filteredIQsample.sample;
    timingSyncBuffer.n = filteredIQsample.sampleIndex;
    sample_complex_t roughSyncedIQsample;

    returnValue = gardnerAlgorithm(&timingSyncBuffer, &roughSyncedIQsample, QAMstate.symbolPeriod, debugPlots);

    timingSyncBuffer.insertionIndex = (timingSyncBuffer.insertionIndex + 1) % timingSyncBuffer.length;

    if(returnValue != RETURNED_SAMPLE)
        return AWAITING_SAMPLES;


    // Costas loop for precise Phase locking
    //  happens every IQ symbol

    //outputSample = ;  // output the determined IQ data
    return RETURNED_SAMPLE;
 }

buffered_data_return_t demodualteOFDM( const sample_double_t *sample, OFDM_properties_t *OFDMstate, debugPlots_t debugPlots)
{
    if(OFDMstate->initialized == 0)
    {
        // initialize the state
        OFDMstate->state.frame = SEARCHING;
        OFDMstate->state.frameStart = sample->sampleIndex;

        OFDMstate->sampleRate = 44100;
        OFDMstate->guardPeriod = (1<<12);  // 2^12=1024 closest power of 2 to the impulse response length, slightly shorter
        //OFDMstate->guardPeriod = 128;
        OFDMstate->ofdmPeriod = OFDMstate->guardPeriod * 4;   // dunno what the best OFDM period is compared to the guard period. I assume longer is better for channel efficiency, but maybe it's worse for noise? don't know
        //OFDMstate->ofdmPeriod = OFDMstate->guardPeriod * 8;   // dunno what the best OFDM period is compared to the guard period. I assume longer is better for channel efficiency, but maybe it's worse for noise? don't know
        OFDMstate->symbolPeriod = OFDMstate->guardPeriod + OFDMstate->ofdmPeriod;
        OFDMstate->channels = OFDMstate->ofdmPeriod / 2 + 1;  // half due to using real symbols (ie, not modulating to higher frequency carrier wave but staying in baseband)

        // initialize buffer to hold input samples for further processing
        //OFDMstate->preambleDetectorInputBuffer.length = OFDMstate->ofdmPeriod + 1;
        OFDMstate->preambleDetectorInputBuffer.length = OFDMstate->symbolPeriod * 1.5;  // a little extra, this will also be the buffer used for FFTs after a transmission is detected
        OFDMstate->preambleDetectorInputBuffer.buffer = calloc(OFDMstate->preambleDetectorInputBuffer.length, sizeof(double));
        OFDMstate->preambleDetectorInputBuffer.sampleRate = OFDMstate->sampleRate;
        if(OFDMstate->preambleDetectorInputBuffer.buffer == NULL)
            fprintf(stderr, "couldn't allocate ofdm symbol buffer.\n");

        // initialize channel simulation buffer
        OFDMstate->simulateChannel = 0;
        if(OFDMstate->simulateChannel)
        {
            OFDMstate->channelSimulationBuffer.length = 5912;
            OFDMstate->channelSimulationBuffer.insertionIndex = 0;
            OFDMstate->channelSimulationBuffer.buffer = calloc(OFDMstate->channelSimulationBuffer.length, sizeof(double));
            if(OFDMstate->channelSimulationBuffer.buffer == NULL)
                fprintf(stderr, "cahnnelSimulationBuffer failed to allocate: %s\n", strerror(errno));

            // initialize impulse response buffer
            OFDMstate->channelImpulseResponse.length = 5912;
            OFDMstate->channelImpulseResponse.insertionIndex = 0;
            OFDMstate->channelImpulseResponse.buffer = calloc(OFDMstate->channelImpulseResponse.length, sizeof(double));
            if(OFDMstate->channelImpulseResponse.buffer == NULL)
                fprintf(stderr, "ImpulseResponse failed to allocate: %s\n", strerror(errno));
        }

        // initialize auto correlation averaging buffer
        OFDMstate->autoCorrelationAverageBuffer.length = 500;
        OFDMstate->autoCorrelationAverageBuffer.insertionIndex = 0;
        OFDMstate->autoCorrelationAverageBuffer.buffer = calloc(OFDMstate->autoCorrelationAverageBuffer.length, sizeof(double));
        if(OFDMstate->autoCorrelationAverageBuffer.buffer == NULL)
            fprintf(stderr, "failed to allocate convolutionAverage buffer. %s\n", strerror(errno));

        // initialize timing filter input buffer
        OFDMstate->timingFilterInputBuffer.length = OFDMstate->guardPeriod * 5;
        OFDMstate->timingFilterInputBuffer.insertionIndex = 0;
        OFDMstate->timingFilterInputBuffer.buffer = calloc(OFDMstate->timingFilterInputBuffer.length, sizeof(double));
        if(OFDMstate->timingFilterInputBuffer.buffer == NULL)
            fprintf(stderr, "Failed to allocate timing input buffer: %s\n", strerror(errno));

        // initialize ofdm symbol buffer
        OFDMstate->fftwSymbolBuffer.length = OFDMstate->ofdmPeriod;
        OFDMstate->fftwSymbolBuffer.insertionIndex = 0;
        OFDMstate->fftwSymbolBuffer.buffer = fftw_malloc(sizeof(double) * OFDMstate->fftwSymbolBuffer.length);
        if(OFDMstate->fftwSymbolBuffer.buffer == NULL)
            fprintf(stderr, "failed to allocate fftw symbol buffer. %s\n", strerror(errno));

        // initialize ofdm symbol buffer
        // I'm assuming that fftw_complex is the same width as double complex, which it seems to be
        OFDMstate->fftwIQbuffer.length = OFDMstate->channels;
        OFDMstate->fftwIQbuffer.insertionIndex = 0;
        OFDMstate->fftwIQbuffer.buffer = fftw_malloc(sizeof(fftw_complex) * OFDMstate->fftwIQbuffer.length);
        if(OFDMstate->fftwIQbuffer.buffer == NULL)
            fprintf(stderr, "failed to allocate fftw IQ buffer. %s\n", strerror(errno));

        printf("double complex: %li, fftw_complex: %li\n", sizeof(double complex), sizeof(fftw_complex));

        OFDMstate->fftwPlan = fftw_plan_dft_r2c_1d(OFDMstate->ofdmPeriod, OFDMstate->fftwSymbolBuffer.buffer, (fftw_complex*)OFDMstate->fftwIQbuffer.buffer, FFTW_MEASURE);

        // run once
        OFDMstate->initialized = 1;
    }

    /*
    sample_double_t equalizedSample;
    // decide whether to simulate the channel response or not
    if(OFDMstate->simulateChannel)
    {
        // channel simulation filter
        OFDMstate->channelSimulationBuffer.n = sample->sampleIndex;
        OFDMstate->channelSimulationBuffer.sampleRate = sample->sampleRate;
        OFDMstate->channelSimulationBuffer.phase = 0;
        OFDMstate->channelSimulationBuffer.buffer[OFDMstate->channelSimulationBuffer.insertionIndex] = sample->sample;

        buffered_data_return_t returnValue = channelFilter(&OFDMstate->channelSimulationBuffer, &equalizedSample, &OFDMstate->channelImpulseResponse, 1, debugPlots);

        OFDMstate->channelSimulationBuffer.insertionIndex = (OFDMstate->channelSimulationBuffer.insertionIndex + 1) % OFDMstate->channelSimulationBuffer.length;

        if(returnValue != RETURNED_SAMPLE)
            return AWAITING_SAMPLES;

    } else {
        
        // skip channel simulation
        equalizedSample = *sample;
    }
    */
    sample_double_t equalizedSample;
    equalizedSample = *sample;

    OFDMstate->preambleDetectorInputBuffer.buffer[OFDMstate->preambleDetectorInputBuffer.insertionIndex] = equalizedSample.sample;
    // add sample to the ofdm symbol window buffer
    OFDMstate->preambleDetectorInputBuffer.n = equalizedSample.sampleIndex;


    // run the state machine
    switch(OFDMstate->state.frame)
    {
        case SEARCHING:


            // auto correlation filter preample detector
            // run a correlation between the first half of the symbol in the buffer, and the second half
            /*
            double convolutionSecondHalfEnergy = 0;
            double convolutionAutocorrelation = 0;
            for(int i = 0; i < OFDMstate->ofdmPeriod / 2; i++)
            {
                //break;
                int firstHalfIndex = (OFDMstate->preambleDetectorInputBuffer.insertionIndex + 1 + i) % OFDMstate->preambleDetectorInputBuffer.length;
                int secondHalfIndex = (firstHalfIndex + OFDMstate->ofdmPeriod / 2) % OFDMstate->preambleDetectorInputBuffer.length;

                //OFDMstate->autoCorrelation.sample += OFDMstate->preambleDetectorInputBuffer.buffer[firstHalfIndex] * OFDMstate->preambleDetectorInputBuffer.buffer[secondHalfIndex];
                convolutionAutocorrelation += OFDMstate->preambleDetectorInputBuffer.buffer[firstHalfIndex] * OFDMstate->preambleDetectorInputBuffer.buffer[secondHalfIndex];
                convolutionSecondHalfEnergy += pow(OFDMstate->preambleDetectorInputBuffer.buffer[secondHalfIndex], 2);
            }
            */
            // simplified iterative version of the auto correlation
            static double iterativeAutocorrelation = 0;;
            static double secondHalfEnergy = 0;
            OFDMstate->autoCorrelation.sample = 0;
            OFDMstate->autoCorrelation.sampleIndex = OFDMstate->preambleDetectorInputBuffer.n - OFDMstate->ofdmPeriod;
            //int secondHalfSecondIndex = (OFDMstate->preambleDetectorInputBuffer.insertionIndex) % OFDMstate->preambleDetectorInputBuffer.length;
            //int secondHalfFirstIndex = (secondHalfSecondIndex - (OFDMstate->preambleDetectorInputBuffer.length / 2 - 1)) % OFDMstate->preambleDetectorInputBuffer.length;

            //int firstHalfFirstIndex = (OFDMstate->preambleDetectorInputBuffer.insertionIndex + 1) % OFDMstate->preambleDetectorInputBuffer.length;
            //int firstHalfSecondIndex = (firstHalfFirstIndex + OFDMstate->preambleDetectorInputBuffer.length / 2 - 1) % OFDMstate->preambleDetectorInputBuffer.length;
            //int secondHalfFirstIndex = (firstHalfFirstIndex + OFDMstate->preambleDetectorInputBuffer.length / 2 - 1) % OFDMstate->preambleDetectorInputBuffer.length;
            //int secondHalfSecondIndex = (secondHalfFirstIndex + OFDMstate->preambleDetectorInputBuffer.length / 2 - 1) % OFDMstate->preambleDetectorInputBuffer.length;

            int firstHalfFirstIndex = (OFDMstate->preambleDetectorInputBuffer.insertionIndex + OFDMstate->preambleDetectorInputBuffer.length - OFDMstate->ofdmPeriod) % OFDMstate->preambleDetectorInputBuffer.length;
            int firstHalfSecondIndex = (firstHalfFirstIndex + OFDMstate->ofdmPeriod / 2) % OFDMstate->preambleDetectorInputBuffer.length;
            int secondHalfFirstIndex = (firstHalfSecondIndex) % OFDMstate->preambleDetectorInputBuffer.length;
            int secondHalfSecondIndex = (secondHalfFirstIndex + OFDMstate->ofdmPeriod / 2) % OFDMstate->preambleDetectorInputBuffer.length;

            if(secondHalfFirstIndex < 0)    // fix possible negatives due to subtraction
                secondHalfFirstIndex += OFDMstate->preambleDetectorInputBuffer.length;
            double point_d =   OFDMstate->preambleDetectorInputBuffer.buffer[firstHalfFirstIndex];      // d    oldest sample
            double point_dl =  OFDMstate->preambleDetectorInputBuffer.buffer[firstHalfSecondIndex];   // d+L  middle sample
            double point_dl2 =  OFDMstate->preambleDetectorInputBuffer.buffer[secondHalfFirstIndex];   // d+L  middle sample
            double point_d2l = OFDMstate->preambleDetectorInputBuffer.buffer[secondHalfSecondIndex];         // d+2L most recent sample 
            iterativeAutocorrelation = iterativeAutocorrelation + 
                point_dl * point_d2l -
                point_d * point_dl;
            secondHalfEnergy = secondHalfEnergy +
                pow(point_d2l, 2) -
                pow(point_dl, 2);

            // I think this causes issues with quiet tones auto correlating to high values and triggering a frame detection. I don't know the solution
            OFDMstate->autoCorrelation.sample = fabs(iterativeAutocorrelation) / secondHalfEnergy; // normalization to remove average gain dependance

            if(secondHalfEnergy == 0)
                OFDMstate->autoCorrelation.sample = 0;  // eliminate divide by zero error

            OFDMstate->autoCorrelation.sample = (iterativeAutocorrelation); // enable if the autocorrelation normalization is to be ignored

            // average filtered auto correlation signal
            OFDMstate->autoCorrelationAverageBuffer.buffer[OFDMstate->autoCorrelationAverageBuffer.insertionIndex] = OFDMstate->autoCorrelation.sample;
            OFDMstate->autoCorrelationAverageBuffer.n = OFDMstate->autoCorrelation.sampleIndex;
            static double runningAverage = 0;       // for the averaging function
            sample_double_t averageValue = {0};
            averagingFilter(&OFDMstate->autoCorrelationAverageBuffer, &averageValue, &runningAverage);
            OFDMstate->autoCorrelationAverageBuffer.insertionIndex = (OFDMstate->autoCorrelationAverageBuffer.insertionIndex + 1) % OFDMstate->autoCorrelationAverageBuffer.length;


            OFDMstate->timingFilterInputBuffer.buffer[OFDMstate->timingFilterInputBuffer.insertionIndex] = averageValue.sample;
            OFDMstate->timingFilterInputBuffer.n = averageValue.sampleIndex;
            sample_double_t timingSignal = {0};
            //timingSignal.sample = 1;
            //buffered_data_return_t returnValue = AWAITING_SAMPLES;
            buffered_data_return_t returnValue = timingPeakDetectionFilter(&OFDMstate->timingFilterInputBuffer, &timingSignal, debugPlots);

            OFDMstate->timingFilterInputBuffer.insertionIndex = (OFDMstate->timingFilterInputBuffer.insertionIndex + 1) % OFDMstate->timingFilterInputBuffer.length;

            // graph the correlation function
            if(debugPlots.OFDMtimingSyncEnabled)
            {
                fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", OFDMstate->autoCorrelation.sampleIndex, 0, OFDMstate->autoCorrelation.sample);
                fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", equalizedSample.sampleIndex, 1, equalizedSample.sample);
                fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", sample->sampleIndex, 2, sample->sample);


                //fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", OFDMstate->autoCorrelationDerivative.sampleIndex, 3, OFDMstate->autoCorrelationDerivative.sample + 1);
                fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", averageValue.sampleIndex, 4, averageValue.sample);
                //fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", averageDerivative.sampleIndex, 5, averageDerivative.sample + 1);

                //fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", OFDMstate->autoCorrelation.sampleIndex + OFDMstate->ofdmPeriod / 2, 7, point_dl * point_d2l);
                //fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", OFDMstate->autoCorrelation.sampleIndex, 8, point_d * point_dl);
                //fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", OFDMstate->autoCorrelation.sampleIndex, 9, convolutionAutocorrelation);
                //fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", OFDMstate->autoCorrelation.sampleIndex, 10, convolutionSecondHalfEnergy);

                //fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", OFDMstate->autoCorrelation.sampleIndex, 11, secondHalfEnergy);
            }

            // preambe triggered
            if(returnValue == RETURNED_SAMPLE)
            {

                // plot timing mark
                if(debugPlots.OFDMtimingSyncEnabled)
                {
                    //fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %i\n", timingSignal.sampleIndex, 6, 0);
                    fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %i\n", timingSignal.sampleIndex, 6, 1);
                    //fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %i\n", timingSignal.sampleIndex, 6, 0);
                }

                // save the timing offset
                OFDMstate->ofdmPhaseOffset = timingSignal.sampleIndex;

                // draw the preamble buffer to the decoder chart
                if(debugPlots.OFDMdecoderEnabled)
                {
                    for(int i = 0; i < OFDMstate->preambleDetectorInputBuffer.length; i++)
                    {
                        int preambleBufferIndex = (OFDMstate->preambleDetectorInputBuffer.insertionIndex + 1 + i) % OFDMstate->preambleDetectorInputBuffer.length;
                        int timeIndex = i - OFDMstate->preambleDetectorInputBuffer.length + 1 + (OFDMstate->preambleDetectorInputBuffer.n - OFDMstate->ofdmPhaseOffset + OFDMstate->guardPeriod);
                        fprintf(debugPlots.OFDMdecoderStdin, "%i %i %f\n",timeIndex, 0, OFDMstate->preambleDetectorInputBuffer.buffer[preambleBufferIndex]);
                    }
                }
                OFDMstate->fftwSymbolBuffer.insertionIndex = 0; // set indexes to the last sample inserted
                OFDMstate->fftwSymbolBuffer.n = 0;
                OFDMstate->fftwSymbolBuffer.sampleRate = OFDMstate->preambleDetectorInputBuffer.sampleRate;

                // change state to active
                OFDMstate->state.frame = ACTIVE;
            }


            break;
        case ACTIVE:

            // check if there are enough samples in the preamble buffer to do the next symbol
            int correctedBufferTime = OFDMstate->preambleDetectorInputBuffer.n - OFDMstate->ofdmPhaseOffset + OFDMstate->guardPeriod;
            static int symbolsProcessed = 0;    // how many symbols have been processed so far
            int symbolIndex = correctedBufferTime / OFDMstate->symbolPeriod;   // calculate the symbol that's currently being filled in the buffer

            // draw the raw sample
            if(debugPlots.OFDMdecoderEnabled)
            {
                fprintf(debugPlots.OFDMdecoderStdin, "%i %i %f\n", correctedBufferTime, 1, OFDMstate->preambleDetectorInputBuffer.buffer[OFDMstate->preambleDetectorInputBuffer.insertionIndex]);
            }

            if(symbolIndex > symbolsProcessed)
            {
                // it's time to tear out an ofdmPeriod of samples into the fftw buffer and do the math
                for(int i = 0; i < OFDMstate->ofdmPeriod; i++)
                {
                    // determine where it is in the buffer
                    int preambleIndex = (OFDMstate->preambleDetectorInputBuffer.insertionIndex + (symbolsProcessed * OFDMstate->symbolPeriod + OFDMstate->guardPeriod + i - correctedBufferTime)) % OFDMstate->preambleDetectorInputBuffer.length;
                    if(preambleIndex < 0)
                        preambleIndex += OFDMstate->preambleDetectorInputBuffer.length;

                    // copy it in order to the fftw buffer
                    OFDMstate->fftwSymbolBuffer.buffer[i] = OFDMstate->preambleDetectorInputBuffer.buffer[preambleIndex];

                    if(debugPlots.OFDMdecoderEnabled)
                    {
                        // draw snipped out portion
                        //int plotIndex = 3 + symbolsProcessed;
                        int plotIndex = 3;
                        fprintf(debugPlots.OFDMdecoderStdin, "%i %i %f\n", symbolsProcessed * OFDMstate->symbolPeriod + i + OFDMstate->guardPeriod, plotIndex, OFDMstate->fftwSymbolBuffer.buffer[i] - 1);
                    }
                }

                // do the fft
                fftw_execute(OFDMstate->fftwPlan);
                
                // debug chart for the subchannel IQ plots
                if(debugPlots.OFDMIQEnabled)
                {
                    // draw a point for every subchannel
                    double normalizationFactor = sqrt(OFDMstate->ofdmPeriod) * 2;
                    for(int k = 0; k < OFDMstate->channels; k++)
                    {
                        // organize plots in a grid from left to right, bottom to top starting at the origin, roughly square
                        int square = (int)ceil(sqrt(OFDMstate->channels));
                        double x = ((k % square));
                        double y = ((k / square));

                        // plot index for coloring preamble samples
                        int plotIndex = 0;
                        if(symbolsProcessed < 2)
                            plotIndex = symbolsProcessed + 1;
                        
                        // plot each point
                        fprintf(debugPlots.OFDMIQStdin, "%f %i %f\n", x + creal(OFDMstate->fftwIQbuffer.buffer[k]) / normalizationFactor, plotIndex, y + cimag(OFDMstate->fftwIQbuffer.buffer[k]) / normalizationFactor);
                    }
                }

                // mark that this symbol was extracted and processed
                symbolsProcessed++;
            }

            break;
    }

    // increment insertion index of input buffer
    OFDMstate->preambleDetectorInputBuffer.insertionIndex = (OFDMstate->preambleDetectorInputBuffer.insertionIndex + 1) % OFDMstate->preambleDetectorInputBuffer.length;

    return AWAITING_SAMPLES;
}

int main(void)
{
    // length of each symbol in samples, and so the fft window.
    // This is the period of time all the orthogonal symbols will be integrated over
#define SYMBOL_PERIOD 32
    // the OFDM channel number, how many cycles per symbol
    int k = 1;
    // buffer is the length of the symbol period, so that symbols are orthogonal
    double sampleBuffer[SYMBOL_PERIOD] = {0.0};

    //int guardPeriod = 4./1000 * 44100;      // 4ms guard period for room echos
    //int guardPeriod = 0;
    //int totalPeriod = SYMBOL_PERIOD + guardPeriod;
    sample_32_converter_t sampleConvert;
    sample_double_t sample;

    // the offset of sample buffer or fftwindow in time, to get symbol time sync
    int windowPhase = 0;
    //double windowPhaseReal = (double)rand() / RAND_MAX * SYMBOL_PERIOD; // the floating point window phase, to be quantized into windowPhase
    double windowPhaseReal = 0;

    // to help the if statements that choose when to take an OFDM DFT sample not get stuck in a loop when the window offset changes
    int tookSampleAt = 0;

    // prepare file descriptors for output files and streams
    int retval = 0;

    debugPlots_t debugPlots;    // holds all the file descriptors for plots
    char *iq_plot_buffer = NULL;

    // process arguments
    //  none right now

    // for live plotting, pipe to feedgnuplot
    const char *plot =
        "feedgnuplot "
        "--domain --lines --points "
        "--title \"Raw signal Time domain\" "
        "--xlabel \"Time (microphone sample #\" --ylabel \"value\" "
        "--legend 0 \"Signal\" "
        "--legend 1 \"equalization factor\" "
    ;
    //const char *plot =
        //"hexdump -C "
        //"tee testoutput.txt"
    //;

    // using it to plot the time domain signal
    debugPlots.waveformPlotStdin = popen(plot, "w");

    if (debugPlots.waveformPlotStdin == NULL)
    {
        fprintf(stderr, "Failed to create waveform plot: %s\n", strerror(errno));
        retval = 1;
        goto exit;
    }

    const char *debugPlot =
        "feedgnuplot "
        "--domain --dataid --lines --points "
        "--title \"FFT debuger graph\" "
        "--xlabel \"Time (microphone sample #\" --ylabel \"value\" "
        "--legend 0 \"input samples\" "
        "--legend 1 \"fft real\" "
        "--legend 2 \"fft imaginary\" "
        "--legend 3 \"I decision\" "
        "--legend 4 \"Q decision\" "
        "--legend 5 \"input samples mid\" "
        "--legend 6 \"fft real mid\" "
        "--legend 7 \"fft imaginary mid\" "
        "--legend 8 \"I decision mid\" "
        "--legend 9 \"Q decision mid\" "
        "--legend 11 \"samp*real\" "
        "--legend 12 \"samp*imag\" "
        "--legend 13 \"integral real\" "
        "--legend 14 \"integral imag\" "
        "--legend 10 \"phase error signal\" "
        "--legend 15 \"window phase\" "
    ;

    // using it to plot the time domain signal
    debugPlots.fftDebuggerStdin = popen(debugPlot, "w");

    if (debugPlots.fftDebuggerStdin == NULL)
    {
        fprintf(stderr, "Failed to create fft debug plot: %s\n", strerror(errno));
        retval = 2;
        goto exit;
    }

    const char *errorPlot =
        "feedgnuplot "
        "--domain --lines --points "
        "--title \"Time Domain error signal\" "
        "--legend 0 \"estimated phase offset\" "
        "--legend 1 \"rolling average error signal\" "
        "--legend 2 \"PI filtered signal\" "
        "--legend 3 \"fft window phase shift\" "
        "--legend 4 \"real target phase offset\" "
    ;

    // using it to plot the time domain signal
    debugPlots.errorPlotStdin = popen(errorPlot, "w");

    if (debugPlots.errorPlotStdin == NULL)
    {
        fprintf(stderr, "Failed to create error plot: %s\n", strerror(errno));
        retval = 3;
        goto exit;
    }

    //FILE* IQplotStdin = popen("feedgnuplot --domain --points --title \"IQ plot\" --unset key", "w");

    // Allocate buffer for call to feedgnuplot for IQ plot
    iq_plot_buffer = malloc(120 + (50 * SYMBOL_PERIOD));

    // Check if allocation failed
    if (iq_plot_buffer == NULL)
    {
        fprintf(stderr, "Failed to allocate buffer for IQ plot: %s\n", strerror(errno));
        retval = 4;
        goto exit;
    }

    int stringLength = 0;

    stringLength += snprintf(iq_plot_buffer, 120,
        "feedgnuplot "
        "--dataid --domain --points --maxcurves %i "
        "--title \"IQ plot\" "
        "--xlabel \"I\" --ylabel \"Q\" ",
        SYMBOL_PERIOD * 2 + 1
    );

    if (stringLength < 0)
    {
        fprintf(stderr, "Printing iq plot buffer failed: %s\n", strerror(errno));
        retval = 5;
        goto exit;
    }
    else if (stringLength == 120)
    {
        fprintf(stderr, "Printing iq plot buffer failed: truncated");
        retval = 5;
        goto exit;
    }

    for(int i = 0; i < SYMBOL_PERIOD; i++)
    {
        stringLength += snprintf(iq_plot_buffer + stringLength, 50,
            "--legend %i \"Phase offset %i samples\" ",
            i,
            i
        );

        if (stringLength < 0)
        {
            fprintf(stderr, "Printing iq plot buffer failed: %s\n", strerror(errno));
            retval = 5;
            goto exit;
        }
        else if (stringLength == 50)
        {
            fprintf(stderr, "Printing iq plot buffer failed: truncated");
            retval = 5;
            goto exit;
        }
    }

    debugPlots.IQplotStdin = popen(iq_plot_buffer, "w");

    // Free buffer, even if popen failed
    free(iq_plot_buffer);
    iq_plot_buffer = NULL;

    if (debugPlots.IQplotStdin == NULL)
    {
        fprintf(stderr, "Failed to create IQ plot: %s\n", strerror(errno));
        retval = 6;
        goto exit;
    }

    // for ploting IQ values over time to hopefully obtain an error function
    char eyeDiagramPlotReal[300] = {0};
    sprintf(
        eyeDiagramPlotReal,
        "feedgnuplot "
        "--domain --dataid --lines --points --maxcurves %i "
        "--title \"Eye Diagram, Real part\" "
        "--xlabel \"Time (Audio sample #)\" --ylabel \"I\" ",
        3000);
    debugPlots.eyeDiagramRealStdin = popen(eyeDiagramPlotReal, "w");

    char eyeDiagramPlotImaginary[300] = {0};
    sprintf(
        eyeDiagramPlotImaginary,
        "feedgnuplot "
        "--domain --dataid --lines --points --maxcurves %i "
        "--title \"Eye Diagram, Imaginary part\" "
        "--xlabel \"Time (Audio sample #)\" --ylabel \"Q\" ",
        3000);
    debugPlots.eyeDiagramImaginaryStdin = popen(eyeDiagramPlotImaginary, "w");

    // using it to plot the time domain signal
    char IQvstimeplot[300] = {0};
    sprintf(
        IQvstimeplot,
        "feedgnuplot "
        "--domain --lines --points "
        "--title \"continuous IQ plot\" "
        "--xlabel \"I\" --ylabel \"Q\" "
        );
    debugPlots.IQvstimeStdin = popen(IQvstimeplot, "w");

    if (debugPlots.IQvstimeStdin == NULL)
    {
        fprintf(stderr, "Failed to create IQ vs Time plot: %s\n", strerror(errno));
        retval = 7;
        goto exit;
    }

    char *filterDebugPlot =
        "feedgnuplot "
        "--domain --dataid --lines --points "
        "--title \"Filter debugging plot\" "
        "--xlabel \"Time (Audio sample #)\" --ylabel \"value\" "
        "--legend 0 \"Original Audio samples\" "
        "--legend 1 \"Original I\" "
        "--legend 2 \"Original Q\" "
        "--legend 3 \"Filter array\" "
        "--legend 4 \"IQ sampler internal exponential real part\" "
        "--legend 5 \"IQ sampler internal exponential imaginary part\" "
        "--legend 6 \"Filtered I\" "
        "--legend 7 \"Filtered Q\" "

    ;
    debugPlots.filterDebugStdin = popen(filterDebugPlot, "w");
    if(debugPlots.filterDebugStdin == NULL)
    {
        fprintf(stderr, "Failed to create filter debug plot: %s\n", strerror(errno));
        retval = 8;
        goto exit;
    }

    char *QAMdemodulatePlot =
        "feedgnuplot "
        "--domain --dataid --lines --points "
        "--title \"QAM demodulation debug plot\" "
        "--xlabel \"Time (sample index)\" --ylabel \"value\" "
        "--legend 0 \"Original audio samples\" "
        "--legend 1 \"Equalized audio samples\" "
        "--legend 5 \"I\" "
        "--legend 6 \"Q\" "
        "--legend 7 \"Phase error signal\" "
    ;
    debugPlots.QAMdecoderStdin = popen(QAMdemodulatePlot, "w");
    if(debugPlots.QAMdecoderStdin == NULL)
    {
        fprintf(stderr, "Failed to create QAM decoder debug plot: %s\n", strerror(errno));
        retval = 9;
        goto exit;
    }

    char *gardenerAlgoPlot =
        "feedgnuplot "
        "--domain --dataid --lines --points "
        "--title \"Gardener algorithm debug plot\" "
        "--xlabel \"Time (sample index)\" --ylabel \"value\" "
        "--legend 0 \"Original audio samples\" "
        "--legend 5 \"I\" "
        "--legend 6 \"Q\" "
        "--legend 7 \"Phase error signal\" "
        "--legend 21 \"Phase error signal\" "
        "--legend 8 \"Post Ideal I\" "
        "--legend 9 \"Post Ideal Q\" "
        "--legend 10 \"Pre ideal I\" "
        "--legend 11 \"Pre ideal Q\" "
        "--legend 12 \"Post Mid I\" "
        "--legend 13 \"Post Mid Q\" "
        "--legend 14 \"Pre Mid I\" "
        "--legend 15 \"Pre Mid Q\" "
        "--legend 16 \"Interpolated ideal I\" "
        "--legend 17 \"Interpolated ideal Q\" "
        "--legend 18 \"Interpolated mid I\" "
        "--legend 19 \"Interpolated mid Q\" "
        "--legend 20 \"sample time Ceil\" "
    ;
    debugPlots.gardnerAlgoStdin = popen(gardenerAlgoPlot, "w");
    if(debugPlots.gardnerAlgoStdin == NULL)
    {
        fprintf(stderr, "Failed to create gardner algo debug plot: %s\n", strerror(errno));
        retval = 10;
        goto exit;
    }

    char *channelFilterPlot =
        "feedgnuplot "
        "--domain --dataid --lines --points "
        "--title \"channel filter debug plot\" "
        "--xlabel \"sample #\" --ylabel \"Value\" "
        "--legend 0 \"input samples\" "
        "--legend 1 \"impulse Response\" "
        "--legend 2 \"channel simulation\" "
        "--legend 9 \"channel correction\" "
        "--legend 3 \"frequency response magnitude\" "
        "--legend 4 \"frequency response phase\" "
        "--legend 5 \"reciprocal frequency response magnitude\" "
        "--legend 6 \"reciprocal frequency response phase\" "
        "--legend 7 \"inverse impulse response magnitude\" "
        "--legend 8 \"inverse impulse response phase\" "
    ;
    debugPlots.channelFilterStdin = popen(channelFilterPlot, "w");
    if(debugPlots.channelFilterStdin == NULL)
    {
        fprintf(stderr, "Failed to create channel filter plot: %s\n", strerror(errno));
        retval = 11;
        goto exit;
    }

    char *OFDMtimingSyncPlot =
        "feedgnuplot "
        "--domain --dataid --lines --points "
        "--title \"OFDM timing sync\" "
        "--xlabel \"sample number\" --ylabel \"value\" "
        "--legend 0 \"preamble auto correlation\" "
        "--legend 1 \"Channel simulated samples\" "
        "--legend 2 \"original samples\" "
        "--legend 3 \"preamble autor correlation derivative\" "
        "--legend 4 \"preamble auto correlation average\" "
        "--legend 5 \"preamble auto correlation derivative average\" "
        "--legend 6 \"OFDM Transmission found, at given offset\" "
        "--legend 7 \"d+l * d+2l\" "
        "--legend 8 \"d * d+l\" "
        "--legend 9 \"convolution auto correlation\" "
        "--legend 10 \"convolution second half energy\" "
        "--legend 11 \"second half energy\" "

    ;
    debugPlots.OFDMtimingSyncStdin = popen(OFDMtimingSyncPlot, "w");
    if(debugPlots.OFDMtimingSyncStdin == NULL)
    {
        fprintf(stderr, "Failed to create OFDM timing sync plot: %s\n", strerror(errno));
        retval = 11;
        goto exit;
    }

    char *OFDMdecoderPlot =
        "feedgnuplot "
        "--domain --dataid --lines "
        "--title \"OFDM symbol decoder\" "
        "--xlabel \"sample number\" --ylabel \"value\" "
        "--legend 0 \"Preamble detection raw samples\" "
        "--legend 1 \"further raw samples\" "
    ;
    debugPlots.OFDMdecoderStdin = popen(OFDMdecoderPlot, "w");
    if(debugPlots.OFDMtimingSyncStdin == NULL)
    {
        fprintf(stderr, "Failed to create OFDM decoder plot: %s\n", strerror(errno));
        retval = 12;
        goto exit;
    }

    char *OFDMIQPlot =
        "feedgnuplot "
        "--domain --dataid --points "
        "--title \"IQ plots for every subchannel\" "
        "--xlabel \"I\" --ylabel \"Q\" "
        "--legend 0 \"Data samples\" "
        "--legend 1 \"Preamble 1 samples\" "
        "--legend 2 \"Preamble 2 samples\" "

    ;
    debugPlots.OFDMIQStdin = popen(OFDMIQPlot, "w");
    if(debugPlots.OFDMIQStdin == NULL)
    {
        fprintf(stderr, "Failed to create OFDM IQ plot: %s\n", strerror(errno));
        retval = 13;
        goto exit;
    }
    debugPlots.flags = 0;   // reset all the flags

    // set some debug flags
    //debugPlots.waveformEnabled = 1;
    //debugPlots.QAMdecoderEnabled = 1;
    //debugPlots.filterDebugEnabled = 1;
    //debugPlots.gardnerAlgoEnabled = 1;
    //debugPlots.eyeDiagramRealEnabled = 1;
    //debugPlots.eyeDiagramImaginaryEnabled = 1;
    //debugPlots.channelFilterEnabled = 1;
    debugPlots.OFDMtimingSyncEnabled = 1;
    debugPlots.OFDMdecoderEnabled = 1;
    debugPlots.OFDMIQEnabled = 1;


    // while there is data to recieve, not end of file -> right now just a fixed number of 2000
    //for(int audioSampleIndex = 0; audioSampleIndex < SYMBOL_PERIOD * 600; audioSampleIndex++)
    for(int audioSampleIndex = 0; audioSampleIndex < 44100 * 60; audioSampleIndex++)
    //for(int audioSampleIndex = 0; audioSampleIndex < SYMBOL_PERIOD * 2000; audioSampleIndex++)
    {
        // recieve data on stdin, signed 32bit integer
        for(size_t i = 0; i < sizeof(sampleConvert.value); i++)
        {
            // get the bytes, little endian, of signed 32bit integer, using the struct to do type punning
            sampleConvert.byte[i] = getchar();
        }


        // convert to a double ranged from -1 to 1
        sample.sample = (double)sampleConvert.value / INT32_MAX;
        sample.sampleRate = 44100;
        sample.sampleIndex = audioSampleIndex;

        if(debugPlots.waveformEnabled)
            fprintf(debugPlots.waveformPlotStdin, "%i %f\n", sample.sampleIndex, sample.sample);

        /*
        QAM_properties_t QAMstate;
        QAMstate.carrierFrequency = (double)sample.sampleRate / SYMBOL_PERIOD;
        QAMstate.carrierPhase = 0;
        QAMstate.k = k;
        //QAMstate.symbolPeriod = (int)(sample.sampleRate / QAMstate.carrierFrequency);
        QAMstate.symbolPeriod = SYMBOL_PERIOD;
        demodulateQAM(&sample, QAMstate, debugPlots);
        */

        static OFDM_properties_t OFDMstate = {0};
        demodualteOFDM(&sample, &OFDMstate, debugPlots);


        //OFDM_properties_t OFDMdemodulateState = {SYMBOL_PERIOD, &sampleBuffer, k};
        //demodulateOFDM(n, sampleDouble, &OFDMdemodulateState, &debugPlots);

    }


exit:
    if (iq_plot_buffer != NULL)
    {
        free(iq_plot_buffer);
        iq_plot_buffer = NULL;
    }

    // close the streams, allowing the plots to render and open a window, then wait for them to terminate in separate threads.
    if((debugPlots.IQvstimeStdin != NULL) && (fork() == 0))
    {
        pclose(debugPlots.IQvstimeStdin);
        return 0;
    }

    if((debugPlots.IQplotStdin != NULL) && (fork() == 0))
    {
        pclose(debugPlots.IQplotStdin);
        return 0;
    }

    if((debugPlots.errorPlotStdin != NULL) && (fork() == 0))
    {
        pclose(debugPlots.errorPlotStdin);
        return 0;
    }

    if((debugPlots.waveformPlotStdin != NULL) && (fork() == 0))
    {
        pclose(debugPlots.waveformPlotStdin);
        return 0;
    }

    if((debugPlots.fftDebuggerStdin != NULL) && (fork() == 0))
    {
        pclose(debugPlots.fftDebuggerStdin);
        return 0;
    }
    if((debugPlots.eyeDiagramRealStdin != NULL) && (fork() == 0))
    {
        pclose(debugPlots.eyeDiagramRealStdin);
        return 0;
    }
    if((debugPlots.eyeDiagramImaginaryStdin != NULL) && (fork() == 0))
    {
        pclose(debugPlots.eyeDiagramImaginaryStdin);
        return 0;
    }
    if((debugPlots.filterDebugStdin != NULL) && (fork() == 0))
    {
        pclose(debugPlots.filterDebugStdin);
        return 0;
    }
    if((debugPlots.gardnerAlgoStdin != NULL) && (fork() == 0))
    {
        pclose(debugPlots.gardnerAlgoStdin);
        return 0;
    }
    if((debugPlots.OFDMtimingSyncStdin != NULL) && (fork() == 0))
    {
        pclose(debugPlots.OFDMtimingSyncStdin);
        return 0;
    }
    if((debugPlots.channelFilterStdin != NULL) && (fork() == 0))
    {
        pclose(debugPlots.channelFilterStdin);
        return 0;
    }
    if((debugPlots.OFDMdecoderStdin != NULL) && (fork() == 0))
    {
        pclose(debugPlots.OFDMdecoderStdin);
        return 0;
    }

    // after all the data is recieved, generate a plot of IQ
    return retval;
}
