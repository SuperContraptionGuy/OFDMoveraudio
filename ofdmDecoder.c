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
#include <time.h>

#include "utilities.h"

buffered_data_return_t interpolateSample(const circular_buffer_double_t *inputSamples, sample_double_t *outputSample, OFDM_state_t *OFDMstate, debugPlots_t debugPlots)
{
    // convert from output sample time to input sample time
    //double output_clock = OFDMstate->samplerAccumulatedPhase; // the fractional number on samples accumulated up to the current input sample
    //int output_sample = OFDMstate->sample;
    //int input_clock = inputSamples->n;

    if(debugPlots.OFDMinterpolatorEnabled)
    {
        fprintf(debugPlots.OFDMinterpolatorStdin, "%i %i %f\n", inputSamples->n, 0, inputSamples->buffer[inputSamples->insertionIndex]);
    }
    // determine if we've recieved enough samples, or need to wait for more
    if(inputSamples->n <= OFDMstate->samplerAccumulatedPhase + ceil((double)inputSamples->length / 2) - 1)    // if last input index is at least half the interpolation filter width in the future, then we can process samples until it's not (The architecture doesn't really allow us to generate multiple output samples for one input sample, so we'll just generate once. This assumes that the input sample rate is always higher than the output sample rate)
        return AWAITING_SAMPLES;

    // if so, use the kernel to interpolate between input samples to generate an output sample
    outputSample->sample = 0;
    for(int i = 0; i < inputSamples->length; i++)
    {
        // multiply by impulse response at given offset, and sum
        int inputSampleIndex = (inputSamples->insertionIndex + 1 + i) % inputSamples->length;  // the index in the buffer of the current input sample being used
        int inputSampleTime = inputSamples->n - inputSamples->length + 1 + i;   // the time of the current input sample
        double filterCoordinate = (double)inputSampleTime - OFDMstate->samplerAccumulatedPhase;    // the position on the interpolation filter to multiply by the input sample

        double inputSample = inputSamples->buffer[inputSampleIndex];
        double filterValue = filterCoordinate < -1 || filterCoordinate > 1 ?    // piece wise function for the 'tent' filter kernel
                                0 :
                                filterCoordinate < 0 ?
                                    filterCoordinate + 1 :
                                    -filterCoordinate + 1;

        outputSample->sample += inputSample * filterValue;  // convolve the filter kernel with the input samples to get one output sample.
                                                            // can't use an actual convolution algo since the output samples' phases will likely change
                                                            // for every output sample as the sample rate ratio is adjusted
    }
    // graph the interpolator debug plot
    if(debugPlots.OFDMinterpolatorEnabled)
    {
        fprintf(debugPlots.OFDMinterpolatorStdin, "%f %i %f\n", OFDMstate->samplerAccumulatedPhase, 1, outputSample->sample);
    }

    // output a sample
    outputSample->sampleRate = OFDMstate->sampleRate;   // the assumed sample rate, though it will be slightly different depending on the  sampling ratio
    outputSample->sampleIndex = OFDMstate->sample;

    // increase the output sample clock's phase
    //OFDMstate->resamplingRatio = ((double)OFDMstate->sampleRate / (1 + 25./1000000))/ inputSamples->sampleRate;
    OFDMstate->resamplingRatio = ((double)OFDMstate->sampleRate * (1 - OFDMstate->samplingFrequencyOffsetEstimate))/ inputSamples->sampleRate;
    OFDMstate->samplerAccumulatedPhase += 1 / OFDMstate->resamplingRatio;   // add one clock duration to the sampler time
    OFDMstate->sample++;

    // graph the offset used on the SFO estimator plot
    if(debugPlots.OFDMsfoEstimatorEnabled && outputSample->sampleIndex % (OFDMstate->symbolPeriod / 100) == 0)
    {
        fprintf(debugPlots.OFDMsfoEstimatorStdin, "%f %i %f\n", ((double)outputSample->sampleIndex - OFDMstate->ofdmPhaseOffset - OFDMstate->ofdmPeriod) / OFDMstate->symbolPeriod, 2, OFDMstate->samplingFrequencyOffsetEstimate * pow(10, 6));
    }

    return RETURNED_SAMPLE;
    // on major issue with this topology is that sometimes I'll need to output two samples after recieving only one.
    // but only if the input sample rate is lower than the output sample rate
    // I'm able to output no samples after recieving some, but not more than one sample. And the function isn't called again until a new input
    // sample is ready. I think the main reason for this issue is that my entire processing function tree is initiated for each recieved sample from the sound card, so effectively it's only called at the sampling frequency, but I may need to call it slightly more often, or a lot more often depending on the resampling ratio
    //
    // One other solution is to simply have an output sample buffer in addition to an input sample buffer,
    // with a value attached that specifies how many new values there are, or some other marker
    // so the consumer knows what samples need to be processed. or else I just have to run
    // the consumer multiple times for each new sample

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
buffered_data_return_t timingPeakDetectionFilter(const circular_buffer_double_t *inputSamples, sample_double_t *outputSample, sample_double_t *windowAverage, OFDM_state_t *OFDMstate, debugPlots_t debugPlots)
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

    sample_double_t thresholdFirstHalf = {0};
    thresholdFirstHalf.sample = inputSamples->buffer[(inputSamples->insertionIndex + 1) % inputSamples->length];    // equals value of first sample
    thresholdFirstHalf.sampleIndex = inputSamples->n - inputSamples->length + 1;
    thresholdFirstHalf.sampleRate = inputSamples->sampleRate;

    sample_double_t thresholdSecondHalf = {0};
    thresholdSecondHalf.sample = inputSamples->buffer[inputSamples->insertionIndex];    // equals value of last sample
    thresholdSecondHalf.sampleIndex = inputSamples->n;
    thresholdSecondHalf.sampleRate = inputSamples->sampleRate;

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
    //double lowerThreshold = 0.2;
    //double upperThreshold = 0.25;
    //double lowerThreshold = windowAverage == 0 ? 0.1 : windowAverage->sample / 3;
    //double upperThreshold = windowAverage == 0 ? 0.1 : windowAverage->sample * 1.5;
    //double lowerThreshold = OFDMstate->receivedRMS.sample == 0 ? 0.5 : OFDMstate->receivedRMS.sample;
    //double upperThreshold = OFDMstate->receivedRMS.sample == 0 ? 1 : OFDMstate->receivedRMS.sample * 2;
    //double upperThreshold = lowerThreshold * 1.25;
    double lowerThreshold = 0.001;
    double upperThreshold = 0.4;
    
    if(debugPlots.OFDMtimingSyncEnabled && inputSamples->n % 500 == 0)
    {
        fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", inputSamples->n, 9, lowerThreshold);
        fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", inputSamples->n, 10, upperThreshold);
    }
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

        double threshold = lowerThreshold;
        //double threshold = maximum.sample * 0.1;
        //threshold = threshold < lowerThreshold ? lowerThreshold : threshold;
        for(int i = 0; i < inputSamples->length; i++)
        {
            double value = inputSamples->buffer[(inputSamples->insertionIndex + 1 + i) % inputSamples->length];
            double index = inputSamples->n - inputSamples->length + 1 + i;
            // get closest to the threshold in the first half
            //double threshold = (lowerThreshold + upperThreshold) / 2;
            if(index < maximum.sampleIndex && fabs(value - threshold) < fabs(thresholdFirstHalf.sample - threshold))
            {
                // determine minimum of first half of samples
                thresholdFirstHalf.sample = value;
                thresholdFirstHalf.sampleIndex = index;
            }
            // get the closest to the threshold in the second half
            if(index >= maximum.sampleIndex && fabs(value - threshold) < fabs(thresholdSecondHalf.sample - threshold))
            {
                // minimum of last half of samples
                thresholdSecondHalf.sample = value;
                thresholdSecondHalf.sampleIndex = index;
            }
        }
        
        outputSample->sample = 0;
        outputSample->sampleIndex = 0;
        outputSample->sampleRate = inputSamples->sampleRate;
        if(minimumSecondHalf.sample < lowerThreshold && minimumFirstHalf.sample < lowerThreshold && maximum.sample > upperThreshold && minimumFirstHalf.sampleIndex < maximum.sampleIndex&& minimumSecondHalf.sampleIndex > maximum.sampleIndex)
        {
            // the offset is the index value of the maximum
            outputSample->sample = 1;
            //outputSample->sampleIndex = maximum.sampleIndex;    // index is the offset of maximum
            int centerOfPlateu = (thresholdFirstHalf.sampleIndex + thresholdSecondHalf.sampleIndex) / 2;    // head a bit ahead of the center of the plateu, near the end of it ideally
            outputSample->sampleIndex = centerOfPlateu + OFDMstate->guardPeriod / 2 * 3 / 4;    // head a bit ahead of the center of the plateu, near the end of it ideally
            //outputSample->sampleIndex = centerOfPlateu;    // head a bit ahead of the center of the plateu, near the end of it ideally
            // graph the correlation function
            if(debugPlots.OFDMtimingSyncEnabled)
            {
                fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", thresholdFirstHalf.sampleIndex, 7, threshold);
                fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", centerOfPlateu - OFDMstate->guardPeriod / 2, 7, maximum.sample);
                fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", centerOfPlateu, 7, maximum.sample);
                fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", centerOfPlateu + OFDMstate->guardPeriod / 2, 7, maximum.sample);
                fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", thresholdSecondHalf.sampleIndex, 7, threshold);
                fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", outputSample->sampleIndex, 6, threshold);
                fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", outputSample->sampleIndex, 6, maximum.sample);
                fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", outputSample->sampleIndex, 6, 1.);

            }
            return RETURNED_SAMPLE; // only return a sample when preamble is detected
        }
    }
    return AWAITING_SAMPLES;
}

buffered_data_return_t demodualteOFDM( const sample_double_t *sample, OFDM_state_t *OFDMstate, debugPlots_t debugPlots)
{
    if(OFDMstate->initialized == 0)
    {
        initializeOFDMstate(OFDMstate);

        OFDMstate->dataOutput = fopen("receiverSequence", "w");

        // initialize the state
        OFDMstate->state.frame = SEARCHING;
        OFDMstate->state.frameStart = sample->sampleIndex;

        // initialize buffer to hold input samples for further processing
        initializeCircularBuffer_double(&OFDMstate->preambleDetectorInputBuffer, OFDMstate->symbolPeriod * 1.5, OFDMstate->sampleRate);

        OFDMstate->disableSFOestimation = 0;
        srand(time(NULL));
        if(OFDMstate->disableSFOestimation)
        {
            //OFDMstate->samplingFrequencyOffsetEstimate = 48/pow(10,6);
            //OFDMstate->samplingFrequencyOffsetEstimate = (rand()%120 - 60)/pow(10, 6); // initial assumed sampling error. can be set to an inital value for testing the estimator
            //OFDMstate->samplingFrequencyOffsetResidual = -OFDMstate->samplingFrequencyOffsetEstimate;
            OFDMstate->samplingFrequencyOffsetEstimate = 0 / pow(10,6);
            OFDMstate->samplingFrequencyOffsetResidual = 0;
        } else {
            //OFDMstate->samplingFrequencyOffsetEstimate = (rand()%120 - 60)/pow(10, 6); // initial assumed sampling error. can be set to an inital value for testing the estimator
            OFDMstate->samplingFrequencyOffsetEstimate = 0/pow(10,6);
            //OFDMstate->samplingFrequencyOffsetEstimate = 0;
            OFDMstate->samplingFrequencyOffsetResidual = 0;
        }
        OFDMstate->pilotSymbolsLength = OFDMstate->channels / OFDMstate->pilotSymbolsPitch;
        OFDMstate->pilotSymbols = malloc(sizeof(complex double) * OFDMstate->pilotSymbolsLength);

        OFDMstate->initialChannelEstimate = malloc(sizeof(complex double) * OFDMstate->channels);    // estimate generated during preamble
        for(int k = 0; k < OFDMstate->channels; k++)
            OFDMstate->initialChannelEstimate[k] = 1;

        OFDMstate->channelEstimate = malloc(sizeof(complex double) * OFDMstate->channels);
        for(int k = 0; k < OFDMstate->channels; k++)
            OFDMstate->channelEstimate[k] = 1;

        OFDMstate->pilotChannelEstimate = malloc(sizeof(complex double) * OFDMstate->pilotSymbolsLength);
        for(int i = 0; i < OFDMstate->pilotSymbolsLength; i++)
            OFDMstate->pilotChannelEstimate[i] = 1;

        // initialize channel simulation buffer
        OFDMstate->simulateChannel = 0;
        if(OFDMstate->simulateChannel)
        {
            initializeOverlapAndSaveBuffer(&OFDMstate->channelSimulationBuffer, 5912);
            // initialize impulse response buffer
            initializeCircularBuffer_double(&OFDMstate->channelImpulseResponse, 5912, 0);
        }

        initializeCircularBuffer_double(&OFDMstate->sampleInterpolatorBuffer, 2, 0);
        initializeCircularBuffer_double(&OFDMstate->autoCorrelationAverageBuffer, 500, 0);
        initializeCircularBuffer_double(&OFDMstate->timingFilterInputBuffer, OFDMstate->symbolPeriod * 1.1, 0);
        // initialize a buffer to store the second preamble symbol, ie, the first of two identical symbols for sample frequency offset estimation
        initializeCircularBuffer_complex(&OFDMstate->sfoFirstSymbol, OFDMstate->channels, 0);

        //initializeCircularBuffer_double(&OFDMstate->OFDMsymbol.timeDomain, OFDMstate->ofdmPeriod, 0);
        // I'm assuming that fftw_complex is the same width as double complex, which it seems to be
        //initializeCircularBuffer_fftw_complex(&OFDMstate->OFDMsymbol.frequencyDomain, OFDMstate->channels, 0);
        // make two buffers to hold the results of DFTs for the first and second half of the symmetric preamble/pilot symbol for sample frequency offset detection
        initializeCircularBuffer_fftw_complex(&OFDMstate->IQrateDetectorFirstHalf.frequencyDomain, OFDMstate->ofdmPeriod / 4 + 1, 0);
        initializeCircularBuffer_fftw_complex(&OFDMstate->IQrateDetectorSecondHalf.frequencyDomain, OFDMstate->ofdmPeriod / 4 + 1, 0);
        // the IQ rate detector time domains point to the first and second half of the time sample buffer
        OFDMstate->IQrateDetectorFirstHalf.timeDomain.buffer = &OFDMstate->OFDMsymbol.timeDomain.buffer[0];
        OFDMstate->IQrateDetectorFirstHalf.timeDomain.length = OFDMstate->OFDMsymbol.timeDomain.length / 2;
        OFDMstate->IQrateDetectorSecondHalf.timeDomain.buffer = &OFDMstate->OFDMsymbol.timeDomain.buffer[OFDMstate->OFDMsymbol.timeDomain.length / 2];
        OFDMstate->IQrateDetectorSecondHalf.timeDomain.length = OFDMstate->OFDMsymbol.timeDomain.length / 2;

        // generate fftw plans
        OFDMstate->fftwPlan = fftw_plan_dft_r2c_1d(
                OFDMstate->OFDMsymbol.timeDomain.length,
                OFDMstate->OFDMsymbol.timeDomain.buffer, 
                (fftw_complex*)OFDMstate->OFDMsymbol.frequencyDomain.buffer, 
                FFTW_MEASURE);
        OFDMstate->fftwRateDetectorFirstHalfPlan = fftw_plan_dft_r2c_1d(
                OFDMstate->IQrateDetectorFirstHalf.timeDomain.length,
                OFDMstate->IQrateDetectorFirstHalf.timeDomain.buffer,
                (fftw_complex*)OFDMstate->IQrateDetectorFirstHalf.frequencyDomain.buffer, 
                FFTW_MEASURE);
        OFDMstate->fftwRateDetectorSecondHalfPlan = fftw_plan_dft_r2c_1d(
                OFDMstate->IQrateDetectorSecondHalf.timeDomain.length,
                OFDMstate->IQrateDetectorSecondHalf.timeDomain.buffer,
                (fftw_complex*)OFDMstate->IQrateDetectorSecondHalf.frequencyDomain.buffer, 
                FFTW_MEASURE);
        // FYI, I never deallocate the buffers or destroy the plans. Could be bad I guess
        // fftw manual says deallocate the arrays with fftw_free() and destroy the plans with fftw_destroy_plan()

        // run once
        OFDMstate->initialized = 1;
    }

    sample_double_t equalizedSample;
    equalizedSample = *sample;

    // run the sample interpolator
    //
    //
    sample_double_t retimedSample;
    //retimedSample = equalizedSample;

    OFDMstate->sampleInterpolatorBuffer.n = equalizedSample.sampleIndex;
    OFDMstate->sampleInterpolatorBuffer.sampleRate = equalizedSample.sampleRate;
    OFDMstate->sampleInterpolatorBuffer.phase = 0;
    OFDMstate->sampleInterpolatorBuffer.buffer[OFDMstate->sampleInterpolatorBuffer.insertionIndex] = equalizedSample.sample;

    buffered_data_return_t returnValue = interpolateSample(&OFDMstate->sampleInterpolatorBuffer, &retimedSample, OFDMstate, debugPlots);
    //buffered_data_return_t returnValue = RETURNED_SAMPLE;
    //retimedSample = equalizedSample;

    OFDMstate->sampleInterpolatorBuffer.insertionIndex = (OFDMstate->sampleInterpolatorBuffer.insertionIndex + 1) % OFDMstate->sampleInterpolatorBuffer.length;

    if(returnValue != RETURNED_SAMPLE)
        return AWAITING_SAMPLES;

    //retimedSample = equalizedSample;

    // then move on to the preamble detector
    OFDMstate->preambleDetectorInputBuffer.buffer[OFDMstate->preambleDetectorInputBuffer.insertionIndex] = retimedSample.sample;
    // add sample to the ofdm symbol window buffer
    OFDMstate->preambleDetectorInputBuffer.n = retimedSample.sampleIndex;

    // run general channel statistics for autocorrelation and RMS determination

    // run the state machine
    switch(OFDMstate->state.frame)
    {
        case SEARCHING:


            // auto correlation filter preample detector
            // simplified iterative version of the auto correlation
            static double iterativeAutocorrelation = 0;
            static double secondHalfEnergy = 0;
            OFDMstate->autoCorrelation.sample = 0;
            OFDMstate->autoCorrelation.sampleIndex = OFDMstate->preambleDetectorInputBuffer.n - OFDMstate->ofdmPeriod;
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

            OFDMstate->receivedRMS.sample = sqrt(secondHalfEnergy / OFDMstate->preambleDetectorInputBuffer.length);
            OFDMstate->receivedRMS.sampleRate = 0;
            OFDMstate->receivedRMS.sampleIndex = 0;

            // I think this causes issues with quiet tones auto correlating to high values and triggering a frame detection. I don't know the solution
            //OFDMstate->autoCorrelation.sample = fabs(iterativeAutocorrelation) / secondHalfEnergy; // normalization to remove average gain dependance
            //if(secondHalfEnergy == 0)
                //OFDMstate->autoCorrelation.sample = 0;  // eliminate divide by zero error
            // uncomment to disable the shitty auto gain
            //OFDMstate->autoCorrelation.sample = (iterativeAutocorrelation / 10); // enable if the autocorrelation normalization is to be ignored
            //OFDMstate->autoCorrelation.sample = sqrt(fabs(iterativeAutocorrelation) / ((double)OFDMstate->preambleDetectorInputBuffer.length / 2)); // enable if the autocorrelation normalization is to be ignored
            double iterAuto_2 = pow(iterativeAutocorrelation, 2);
            double secHaEn_2 = pow(secondHalfEnergy, 2);
            OFDMstate->autoCorrelation.sample = secHaEn_2 == 0 ? 0 : iterAuto_2/ secHaEn_2; // enable if the autocorrelation normalization is to be ignored

            // average filtered auto correlation signal
            OFDMstate->autoCorrelationAverageBuffer.buffer[OFDMstate->autoCorrelationAverageBuffer.insertionIndex] = OFDMstate->autoCorrelation.sample;
            OFDMstate->autoCorrelationAverageBuffer.n = OFDMstate->autoCorrelation.sampleIndex;
            static double runningAverage = 0;       // for the averaging function
            sample_double_t averageValue = {0};
            averagingFilter(&OFDMstate->autoCorrelationAverageBuffer, &averageValue, &runningAverage);
            OFDMstate->autoCorrelationAverageBuffer.insertionIndex = (OFDMstate->autoCorrelationAverageBuffer.insertionIndex + 1) % OFDMstate->autoCorrelationAverageBuffer.length;

            OFDMstate->timingFilterInputBuffer.buffer[OFDMstate->timingFilterInputBuffer.insertionIndex] = averageValue.sample;
            OFDMstate->timingFilterInputBuffer.n = averageValue.sampleIndex;

            static double autocorrelationAverage = 0;
            sample_double_t windowAverage = {0};
            averagingFilter(&OFDMstate->timingFilterInputBuffer, &windowAverage, &autocorrelationAverage);

            sample_double_t timingSignal = {0};
            //timingSignal.sample = 1;
            //buffered_data_return_t returnValue = AWAITING_SAMPLES;
            buffered_data_return_t returnValue = timingPeakDetectionFilter(&OFDMstate->timingFilterInputBuffer, &timingSignal, &windowAverage, OFDMstate, debugPlots);

            OFDMstate->timingFilterInputBuffer.insertionIndex = (OFDMstate->timingFilterInputBuffer.insertionIndex + 1) % OFDMstate->timingFilterInputBuffer.length;

            // graph the correlation function
            if(debugPlots.OFDMtimingSyncEnabled && retimedSample.sampleIndex % 500 == 0)
            {
                fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", retimedSample.sampleIndex, 1, retimedSample.sample);
                //fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", averageValue.sampleIndex, 8, windowAverage.sample);
                fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", averageValue.sampleIndex, 4, averageValue.sample);
                fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", averageValue.sampleIndex, 11, OFDMstate->receivedRMS.sample);
                //fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", averageValue.sampleIndex, 12, iterAuto_2);
                //fprintf(debugPlots.OFDMtimingSyncStdin, "%i %i %f\n", averageValue.sampleIndex, 13, secHaEn_2);
            }

            // preambe triggered
            if(returnValue == RETURNED_SAMPLE)
            {

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
                OFDMstate->OFDMsymbol.timeDomain.insertionIndex = 0; // set indexes to the last sample inserted
                OFDMstate->OFDMsymbol.timeDomain.n = 0;
                OFDMstate->OFDMsymbol.timeDomain.sampleRate = OFDMstate->preambleDetectorInputBuffer.sampleRate;

                // change state to active
                OFDMstate->state.frame = ACTIVE;
                OFDMstate->state.field = PREAMBLE;
                OFDMstate->state.symbolIndex = 0;
                OFDMstate->state.processedSymbols = 0;      // how many symbols have been processed so far total
            }


            break;
        case ACTIVE:


            // check if there are enough samples in the preamble buffer to do the next symbol
            int correctedBufferTime = OFDMstate->preambleDetectorInputBuffer.n - OFDMstate->ofdmPhaseOffset + OFDMstate->guardPeriod;
            int symbolIndex = correctedBufferTime / OFDMstate->symbolPeriod;   // calculate the symbol that's currently being filled in the buffer

            // draw the raw samples for debug purposes
            if(debugPlots.OFDMdecoderEnabled)
            {
                fprintf(debugPlots.OFDMdecoderStdin, "%i %i %f\n", correctedBufferTime, 1, OFDMstate->preambleDetectorInputBuffer.buffer[OFDMstate->preambleDetectorInputBuffer.insertionIndex]);
            }

            // check if there are enough samples in the input buffers to process the next sample
            // ie, is the calculated symbolIndex greater than the last processed symbol index
            if(symbolIndex > OFDMstate->state.processedSymbols)
            {
                // it's time to tear out an ofdmPeriod of samples into the fftw buffer and do the math
                for(int i = 0; i < OFDMstate->ofdmPeriod; i++)
                {
                    // determine where it is in the buffer
                    int preambleIndex = (OFDMstate->preambleDetectorInputBuffer.insertionIndex + (OFDMstate->state.processedSymbols * OFDMstate->symbolPeriod + OFDMstate->guardPeriod + i - correctedBufferTime)) % OFDMstate->preambleDetectorInputBuffer.length;
                    if(preambleIndex < 0)
                        preambleIndex += OFDMstate->preambleDetectorInputBuffer.length;

                    // copy it in order to the fftw buffer
                    // correct for estimated RMS
                    OFDMstate->OFDMsymbol.timeDomain.buffer[i] = OFDMstate->preambleDetectorInputBuffer.buffer[preambleIndex] / OFDMstate->receivedRMS.sample;

                    if(debugPlots.OFDMdecoderEnabled)
                    {
                        // draw snipped out portion
                        //int plotIndex = 3 + OFDMstate->state.processedSymbols;
                        int plotIndex = 3;
                        fprintf(debugPlots.OFDMdecoderStdin, "%li %i %f\n", OFDMstate->state.processedSymbols * OFDMstate->symbolPeriod + i + OFDMstate->guardPeriod, plotIndex, OFDMstate->OFDMsymbol.timeDomain.buffer[i] - 7);
                    }
                }

                // do the fft
                fftw_execute(OFDMstate->fftwPlan);
                // correct for channel estimate
                for(int k = 0; k < OFDMstate->channels; k++)
                {
                    // normalize the dft and correct for channel effects
                    //OFDMstate->OFDMsymbol.frequencyDomain.buffer[k] /= OFDMstate->channelEstimate[k] * sqrt(OFDMstate->ofdmPeriod);
                    OFDMstate->OFDMsymbol.frequencyDomain.buffer[k] /= sqrt(OFDMstate->ofdmPeriod);
                }

                // debug chart for the subchannel IQ plots
                if(debugPlots.OFDMrawIQEnabled)
                {
                    // draw a point for every subchannel
                    for(int k = 0; k < OFDMstate->channels; k++)
                    {
                        // plot index for coloring preamble samples
                        int plotIndex = 0;
                        if(OFDMstate->state.field == PREAMBLE && OFDMstate->state.symbolIndex < 2)   // first preamble symbol
                            plotIndex = 1;
                        else if(OFDMstate->state.field == PREAMBLE && OFDMstate->state.symbolIndex < 4)   // first preamble symbol
                            plotIndex = 2;
                        else if(OFDMstate->state.field == DATA && OFDMstate->state.symbolIndex < 10)   // first preamble symbol
                            plotIndex = 3;

                        plotPointPerSubchannel(
                                debugPlots.OFDMrawIQStdin,
                                OFDMstate->OFDMsymbol.frequencyDomain.buffer[k],
                                k, 
                                8, 
                                OFDMstate->channels,
                                plotIndex);
                    }
                }


                // then do different things with the data depending on the field
                switch(OFDMstate->state.field)
                {
                    case PREAMBLE:

                        // use the symbols of the preamble to estimate sampling frequency offset and do channel estimation

                        // random integer for reconstructing signal
                        long int randomInteger;

                        // do SFO estimation on the preamble samples by running the FFT on the first half and the second half independently
                        // this is for symetric symbols (just the first two)
                        if(OFDMstate->state.symbolIndex < 2)
                        {
                            fftw_execute(OFDMstate->fftwRateDetectorFirstHalfPlan);
                            fftw_execute(OFDMstate->fftwRateDetectorSecondHalfPlan);

                            double samplingFrequencyOffsetEstimate = 0;
                            double normalizationFactor = 0;
                            // calculate sampling frequency offset (skip DC and highest frequency in the DFTs)
                            //for(int i = 1; i < OFDMstate->channels / 2 - 1; i++)
                            for(int k = 1; k < OFDMstate->channels - 1; k++)
                            {
                                if(k % 2 == 0)
                                {
                                    int i = k / 2;
                                    constellation_complex_t constellation = OFDMstate->constellations[0];
                                    lrand48_r(&OFDMstate->preamblePilotsPRNG, &randomInteger);
                                    complex double expectedIQ = constellation.points[randomInteger % constellation.length];

                                    // estimate sampling frequency offset
                                    // does not take into account the length of a guard interval inbetween
                                    // result is in samplingFrequencyOffsetEstimateratio of sampling frequency offset
                                    samplingFrequencyOffsetEstimate +=
                                    (
                                        carg(
                                            OFDMstate->IQrateDetectorFirstHalf.frequencyDomain.buffer[i] /
                                            OFDMstate->IQrateDetectorSecondHalf.frequencyDomain.buffer[i]
                                        ) / (2 * M_PI * i) * i
                                    );

                                    normalizationFactor += i;
                                }
                            }
                            samplingFrequencyOffsetEstimate /= normalizationFactor;

                            if(debugPlots.OFDMsfoEstimatorEnabled)
                            {
                                fprintf(
                                    debugPlots.OFDMsfoEstimatorStdin,
                                    "%i %i %f\n",
                                    OFDMstate->state.symbolIndex,
                                    0,
                                    OFDMstate->samplingFrequencyOffsetEstimate * pow(10, 6)
                                );

                                fprintf(debugPlots.OFDMsfoEstimatorStdin, "%i %i %f\n", OFDMstate->state.symbolIndex, 1, samplingFrequencyOffsetEstimate * pow(10,6));
                            }
                            // adjust the offset estimate
                            if(!OFDMstate->disableSFOestimation)
                                OFDMstate->samplingFrequencyOffsetEstimate += samplingFrequencyOffsetEstimate * 0.9;// * 0.7;// *  0.70;
                            // there is an issue with this correction, since the offset may not take effect before the next OFDM due to buffering
                        }
                        // then do SFO estimation again on the on the next two symbol
                        // two or three steps should be enough to get close, then we'll use pilot symbols to track the phases and correct for residual SFO
                        // using phase corrections only rather than resampling corrections
                        // we'll also use the third preamble symbol for channel estimation, so estimation will happen before the fine SFO correction, should be ok?

                        // the next two conditions are for full length repeated symbols
                        // I'm not doing anything with the first preamble symbol right now
                        // for the second, store the decoded IQ values into a buffer
                        else if(OFDMstate->state.symbolIndex % 2 == 0)
                        {
                            for(int k = 1; k < OFDMstate->channels - 1; k++)    // skip DC and Niquist
                            {
                                constellation_complex_t constellation = OFDMstate->constellations[0];
                                lrand48_r(&OFDMstate->preamblePilotsPRNG, &randomInteger);
                                complex double expectedIQ = constellation.points[randomInteger % constellation.length];
                                complex double receivedIQ = OFDMstate->OFDMsymbol.frequencyDomain.buffer[k];

                                // copy values
                                OFDMstate->sfoFirstSymbol.buffer[k] = receivedIQ;

                                // measure the channel response on each corrected symbol
                                OFDMstate->initialChannelEstimate[k] = receivedIQ / expectedIQ;
                                OFDMstate->channelEstimate[k] = OFDMstate->initialChannelEstimate[k];

                                /*
                                fprintf(
                                    OFDMstate->dataOutput,
                                    "n=%i k=%i %li: %lf+%lfi : %lf+%lfi : %lf+%lfi abs-> %lf %li\n",
                                    OFDMstate->state.symbolIndex,
                                    k,
                                    randomInteger,
                                    creal(expectedIQ),
                                    cimag(expectedIQ),
                                    creal(receivedIQ),
                                    cimag(receivedIQ),
                                    creal(OFDMstate->initialChannelEstimate[k]),
                                    cimag(OFDMstate->initialChannelEstimate[k]),
                                    cabs(OFDMstate->initialChannelEstimate[k]),
                                    OFDMstate->state.processedSymbols
                                );
                                */
                            }
                        }
                        // and on the third, do some calculations to estimate the sampling frequency offset
                        else if(OFDMstate->state.symbolIndex % 2 == 1)
                        {
                            double samplingFrequencyOffsetEstimate = 0;
                            double normalizationFactor = 0;
                            // calculate sampling frequency offset (skip DC and highest frequency
                            for(int k = 1; k < OFDMstate->channels - 1; k++)
                            {
                                // estimate sampling frequency offset
                                // takes into account the additional time for phasing due to the guard period, which was not included in Sliskovic2001, due to my usage of two independant complete OFDM symbols with a guard period in between
                                // result is in samplingFrequencyOffsetEstimateratio of sampling frequency offset
                                samplingFrequencyOffsetEstimate += carg(OFDMstate->sfoFirstSymbol.buffer[k] / OFDMstate->OFDMsymbol.frequencyDomain.buffer[k]) / (2 * M_PI * k * ((double)OFDMstate->guardPeriod / OFDMstate->ofdmPeriod + 1)) * k;
                                normalizationFactor += k;
                            }
                            samplingFrequencyOffsetEstimate /= normalizationFactor;


                            if(debugPlots.OFDMsfoEstimatorEnabled)
                            {
                                fprintf(debugPlots.OFDMsfoEstimatorStdin, "%i %i %f\n", OFDMstate->state.symbolIndex, 0, OFDMstate->samplingFrequencyOffsetEstimate * pow(10, 6));
                                fprintf(debugPlots.OFDMsfoEstimatorStdin, "%i %i %f\n", OFDMstate->state.symbolIndex, 1, samplingFrequencyOffsetEstimate * pow(10,6));
                            }


                            // adjust the offset estimate
                            samplingFrequencyOffsetEstimate *= 1;// weight
                            if(!OFDMstate->disableSFOestimation)
                                OFDMstate->samplingFrequencyOffsetEstimate += samplingFrequencyOffsetEstimate;

                        }


                        // check if it's time to exit the preamble
                        if(OFDMstate->state.symbolIndex == 3)
                        {
                            OFDMstate->state.field = DATA;
                            OFDMstate->state.symbolIndex = -1;   // reset the symbol index for the data field
                        }
                        break;

                    case DATA:

                        double samplingFrequencyOffsetEstimate = 0;
                        double normalizationFactor = 0;

                        // for each channel, estimate SFO and the channel, for pilots
                        for(int k = 1; k < OFDMstate->channels - 1; k++)
                        {

                            // calculating original data and pilot symbols
                            long int randomIntegerPilot;
                            //long int randomIntegerData;
                            //lrand48_r(&OFDMstate->predefinedDataPRNG, &randomIntegerPilot);
                            //lrand48_r(&OFDMstate->predefinedDataPRNG, &randomIntegerData);

                            // grab pilot symbols and calculate an estimate for sampling frequency offset
                            if(k % OFDMstate->pilotSymbolsPitch == 0)   // pilot channels
                            {
                                if(OFDMstate->state.symbolIndex == 0)
                                {
                                    // first pilot symbol
                                    //for(int i = 1; i < OFDMstate->pilotSymbolsLength + 1; i++)
                                    //for(int k = 1; k < OFDMstate->channels - 1; k++)
                                    {
                                        //if(k%OFDMstate->pilotSymbolsPitch == 0)
                                        {
                                            int i = k / OFDMstate->pilotSymbolsPitch;   // index into pilot symbol arrays
                                            constellation_complex_t constellation = OFDMstate->constellations[0];
                                            lrand48_r(&OFDMstate->pilotsPRNG, &randomIntegerPilot);
                                            complex double expectedIQ = constellation.points[randomIntegerPilot % constellation.length];    // multiply by the BPSK expected pilot to remove randomness
                                            complex double receivedIQ = OFDMstate->OFDMsymbol.frequencyDomain.buffer[k];
                                            complex double currentPilotSymbol = receivedIQ * expectedIQ;    // multiply by the BPSK expected pilot to remove randomness
                                            // store pilot symbol for next symbol calculations
                                            OFDMstate->pilotSymbols[i] = currentPilotSymbol;
                                            //fprintf(OFDMstate->dataOutput, "n=%i k=%i %li: %lf+%lfi : %lf+%lfi : %lf+%lfi %li\n", OFDMstate->state.symbolIndex, k, randomIntegerPilot, creal(expectedIQ), cimag(expectedIQ), creal(receivedIQ), cimag(receivedIQ), creal(currentPilotSymbol), cimag(currentPilotSymbol), OFDMstate->state.processedSymbols);
                                        }
                                    }
                                } else if(OFDMstate->state.symbolIndex > 0)
                                {
                                    // for each channel skip DC and highest frequency
                                    //for(int k = 1; k < OFDMstate->channels - 1; k++)
                                    {
                                    //for(int i = 1; i < OFDMstate->pilotSymbolsLength + 1; i++)
                                        //if(k % OFDMstate->pilotSymbolsPitch == 0)
                                        //for(int i = OFDMstate->pilotSymbolsLength / 5; i < OFDMstate->pilotSymbolsLength * 3 / 5; i++)
                                        {
                                            // calculate sampling frequency offset 
                                            // calculate channel index from pilot index
                                            //int k = i * OFDMstate->pilotSymbolsPitch;
                                            int i = k / OFDMstate->pilotSymbolsPitch;   // index into pilot symbol arrays


                                            // estimate sampling frequency offset
                                            // takes into account the additional time for phasing due to the guard period, which was not included in Sliskovic2001, due to my usage of two independant complete OFDM symbols with a guard period in between
                                            // result is in samplingFrequencyOffsetEstimateratio of sampling frequency offset
                                            constellation_complex_t constellation = OFDMstate->constellations[0];
                                            lrand48_r(&OFDMstate->pilotsPRNG, &randomIntegerPilot);
                                            complex double expectedIQ = constellation.points[randomIntegerPilot % constellation.length];    // multiply by the BPSK expected pilot to remove randomness
                                            complex double receivedIQ = OFDMstate->OFDMsymbol.frequencyDomain.buffer[k];
                                            complex double currentPilotSymbol= receivedIQ * expectedIQ;    // multiply by the BPSK expected pilot to remove randomness
                                            /*
                                            fprintf(
                                                    OFDMstate->dataOutput,
                                                    "n=%i k=%i %li: %lf+%lfi : %lf+%lfi : %lf+%lfi %li\n",
                                                    OFDMstate->state.symbolIndex,
                                                    k,
                                                    randomIntegerPilot,
                                                    creal(expectedIQ),
                                                    cimag(expectedIQ),
                                                    creal(receivedIQ),
                                                    cimag(receivedIQ),
                                                    creal(currentPilotSymbol),
                                                    cimag(currentPilotSymbol),
                                                    OFDMstate->state.processedSymbols
                                                   );
                                               */

                                            // estimate channel for pilots
                                            OFDMstate->pilotChannelEstimate[i] = receivedIQ / expectedIQ / OFDMstate->initialChannelEstimate[k];


                                            // estimate SFO
                                            double samplingFrequencyOffsetEstimateChange = 
                                                carg(OFDMstate->pilotSymbols[i] / (currentPilotSymbol)) 
                                                / (2 * M_PI * k * ((double)OFDMstate->guardPeriod / OFDMstate->ofdmPeriod + 1))
                                                * k;

                                            // move pilot value to the old buffer
                                            OFDMstate->pilotSymbols[i] = currentPilotSymbol;

                                            // testing the idea of excluding noisy pilot channels
                                            // according to the amplitude of it's channel estimation
                                            // Should update this in the future to use the calculated SNR of
                                            // a particular pilot channel
                                            if(cabs(OFDMstate->channelEstimate[k]) < 0.5)
                                                continue;

                                            samplingFrequencyOffsetEstimate += samplingFrequencyOffsetEstimateChange;
                                            normalizationFactor += k;
                                        }
                                    }

                                }
                            } else { // data channels.I think this needs to be moved to it's own for loop, so the equalizer can have up to date info from pilots for this symbol
                                     // doing nothing here
                            }
                        }
                        if(OFDMstate->state.symbolIndex > 0)    // calculations for SFO estimation using pilot channels
                        {
                            // finish calculations for samplingFrequencyOffset
                            samplingFrequencyOffsetEstimate = normalizationFactor == 0 ? 0 : samplingFrequencyOffsetEstimate / normalizationFactor;

                            if(debugPlots.OFDMsfoEstimatorEnabled)
                                fprintf(debugPlots.OFDMsfoEstimatorStdin, 
                                        "%i %i %f\n"
                                        "%i %i %f\n"
                                        "%i %i %f\n",
                                        OFDMstate->state.symbolIndex + 4, 3, OFDMstate->samplingFrequencyOffsetResidual * pow(10, 6),
                                        OFDMstate->state.symbolIndex + 4, 4, samplingFrequencyOffsetEstimate * pow(10,6),
                                        OFDMstate->state.symbolIndex + 4, 5, (OFDMstate->samplingFrequencyOffsetResidual + OFDMstate->samplingFrequencyOffsetEstimate) * pow(10, 6));

                            // decaying weight, with a minimum floor
                            double initialValue = 1;
                            double decayConstant = 10;   // time in samples until the weight falls to 1/e of initial value
                            double minimum = 0.01;      // minimum value
                            //double minimum = initialValue;      // disables the exponential decay bit
                            double weight = initialValue * exp(-OFDMstate->state.symbolIndex / decayConstant);
                            weight = weight < minimum ? minimum : weight;
                            // adjust the offset estimate
                            if(!OFDMstate->disableSFOestimation)
                            {
                                // store the offset into the residual for phase corrections since initial channel estimation
                                OFDMstate->samplingFrequencyOffsetResidual += samplingFrequencyOffsetEstimate * weight;// *  0.70;
                                                                                                                       // correct the resampler rate for future symbols
                                OFDMstate->samplingFrequencyOffsetEstimate += samplingFrequencyOffsetEstimate * weight;// *  0.70;
                            }
                        }

                        double errorRate = 0;
                        // now process each data channel
                        for(int k = 1; k < OFDMstate->channels - 1; k++)
                        {
                            // correct for channel equalization
                            // no pilot info available for first symbol
                            if(OFDMstate->state.symbolIndex != 0)
                            {
                                // update the channelEstimate using information from the pilots
                                int pilotBeforeIndex = floor((double)k / OFDMstate->pilotSymbolsPitch);
                                pilotBeforeIndex = pilotBeforeIndex < 1 ? 1 : pilotBeforeIndex;
                                double pilotBeforeWeight = 1 - (double)(k % OFDMstate->pilotSymbolsPitch - 1) / OFDMstate->pilotSymbolsPitch;
                                int pilotAfterIndex = ceil((double)k / OFDMstate->pilotSymbolsPitch);
                                pilotAfterIndex = pilotAfterIndex > OFDMstate->pilotSymbolsLength - 2 ? OFDMstate->pilotSymbolsLength - 2 : pilotAfterIndex;
                                double pilotAfterWeight = 1 - pilotBeforeWeight;

                                // linear interpolation between nearby pilot symbol channel estimates
                                complex double newEstimate =
                                    pilotBeforeWeight * OFDMstate->pilotChannelEstimate[pilotBeforeIndex] +
                                    pilotAfterWeight * OFDMstate->pilotChannelEstimate[pilotAfterIndex];
                                newEstimate *= OFDMstate->initialChannelEstimate[k];

                                // limit corrections to high SNR pilots
                                //if(cabs(OFDMstate->pilotChannelEstimate[pilotBeforeIndex]) > 0.5 && 
                                        //cabs(OFDMstate->pilotChannelEstimate[pilotAfterIndex]) > 0.5)
                                if(1)   // use all pilots
                                {
                                    if(OFDMstate->state.symbolIndex == 0)
                                    {
                                        OFDMstate->channelEstimate[k] = newEstimate;
                                    } else {
                                        OFDMstate->channelEstimate[k] = newEstimate;
                                    }
                                }
                            }

                            // correct for channel estimate, only for the data channels
                            OFDMstate->OFDMsymbol.frequencyDomain.buffer[k] /= OFDMstate->channelEstimate[k];

                            // plot corrected IQ data
                            if(debugPlots.OFDMdataIQEnabled)
                            {
                                plotPointPerSubchannel(
                                        debugPlots.OFDMdataIQStdin,
                                        OFDMstate->OFDMsymbol.frequencyDomain.buffer[k],
                                        k, 
                                        sqrt(2) * 4, 
                                        OFDMstate->channels,
                                        0);
                            }
                            // make sure it's not a pilot channel, otherwise the random number gen gets out of sync
                            if(k % OFDMstate->pilotSymbolsPitch != 0)
                            {
                                long int randomIntegerData;

                                // run a discriminator to classify the recieved symbols
                                // TODO

                                constellation_complex_t constellation = OFDMstate->constellations[k%(OFDMstate->constellationsLength - 1) + 1];
                                complex double receivedIQ = OFDMstate->OFDMsymbol.frequencyDomain.buffer[k];
                                // I'll try an easy inefficient algorithm, just test for nearest constellation point
                                int decodedConstellationIndex = quantizeIQsample(&constellation, receivedIQ);


                                if(0)
                                {
                                    // checking the accuracy of the channel estimation given known data
                                    // picking a sequential constellation, and a random point in that constellation discluding the first entry
                                    lrand48_r(&OFDMstate->predefinedDataPRNG, &randomIntegerData);
                                    int expectedConstellationIndex = randomIntegerData % constellation.length;
                                    complex double expectedIQ = constellation.points[expectedConstellationIndex];
                                    //fprintf(OFDMstate->dataOutput, "n=%i k=%i %li: %lf+%lfi : %lf+%lfi %li\n", OFDMstate->state.symbolIndex, k, randomIntegerData, creal(expectedIQ), cimag(expectedIQ), creal(receivedIQ), cimag(receivedIQ), OFDMstate->state.processedSymbols);
                                    if(decodedConstellationIndex != expectedConstellationIndex)
                                        errorRate += 1. / OFDMstate->channels;
                                        /*
                                        fprintf(OFDMstate->dataOutput,
                                                "n=%i k=%i integer=%li expected vs received Index=\t%i,%i\n",
                                                OFDMstate->state.symbolIndex,
                                                k,
                                                randomIntegerData,
                                                expectedConstellationIndex,
                                                decodedConstellationIndex);
                                        */
                                    if(expectedIQ != 0) // make sure we can actually make an estimate, expected symbol isn't at the origin
                                    {
                                        complex double estimatedEqualizerError = receivedIQ / expectedIQ;

                                        // debug chart for the subchannel IQ plots
                                        if(debugPlots.OFDMequalizerErrorIQEnabled)
                                        {
                                            plotPointPerSubchannel(
                                                    debugPlots.OFDMequalizerErrorIQStdin,
                                                    estimatedEqualizerError,
                                                    k, 
                                                    4, 
                                                    OFDMstate->channels,
                                                    k);
                                        }
                                    }
                                } else {
                                    reverseHuffmanTree(OFDMstate, &constellation, decodedConstellationIndex);
                                }
                            }

                            // plot channel estimate for all subchannels
                            if(debugPlots.OFDMequalizerIQEnabled)
                            {
                                plotPointPerSubchannel(
                                        debugPlots.OFDMequalizerIQStdin,
                                        OFDMstate->channelEstimate[k],
                                        k, 
                                        4, 
                                        OFDMstate->channels,
                                        k + OFDMstate->channels);
                            }
                        }
                        /*
                        fprintf(OFDMstate->dataOutput,
                                "n=%i errorRate=%f\n",
                                OFDMstate->state.symbolIndex,
                                errorRate);
                                */



                        break;
                }


                // mark that this symbol was extracted and processed, increment the symbol index
                OFDMstate->state.processedSymbols++;
                OFDMstate->state.symbolIndex++;
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

    // choose debug plots
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
    //debugPlots.OFDMdecoderEnabled = 1;
    //debugPlots.OFDMrawIQEnabled = 1;
    debugPlots.OFDMsfoEstimatorEnabled = 1;
    //debugPlots.OFDMinterpolatorEnabled = 1;
    debugPlots.OFDMequalizerErrorIQEnabled = 1;
    debugPlots.OFDMequalizerIQEnabled = 1;
    debugPlots.OFDMdataIQEnabled = 1;

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
    if(debugPlots.waveformEnabled)
    {
        debugPlots.waveformPlotStdin = popen(plot, "w");

        if (debugPlots.waveformPlotStdin == NULL)
        {
            fprintf(stderr, "Failed to create waveform plot: %s\n", strerror(errno));
            retval = 1;
            goto exit;
        }
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

    if(debugPlots.fftDebugEnabled)
    {
        debugPlots.fftDebuggerStdin = popen(debugPlot, "w");

        if (debugPlots.fftDebuggerStdin == NULL)
        {
            fprintf(stderr, "Failed to create fft debug plot: %s\n", strerror(errno));
            retval = 2;
            goto exit;
        }
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

    if(debugPlots.errorPlotEnabled)
    {
        debugPlots.errorPlotStdin = popen(errorPlot, "w");

        if (debugPlots.errorPlotStdin == NULL)
        {
            fprintf(stderr, "Failed to create error plot: %s\n", strerror(errno));
            retval = 3;
            goto exit;
        }
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

    if(debugPlots.IQplotEnabled)
    {
        debugPlots.IQplotStdin = popen(iq_plot_buffer, "w");


        if (debugPlots.IQplotStdin == NULL)
        {
            fprintf(stderr, "Failed to create IQ plot: %s\n", strerror(errno));
            retval = 6;
            goto exit;
        }
    }
    // Free buffer, even if popen failed
    free(iq_plot_buffer);
    iq_plot_buffer = NULL;


    // for ploting IQ values over time to hopefully obtain an error function
    char eyeDiagramPlotReal[300] = {0};
    sprintf(
        eyeDiagramPlotReal,
        "feedgnuplot "
        "--domain --dataid --lines --points --maxcurves %i "
        "--title \"Eye Diagram, Real part\" "
        "--xlabel \"Time (Audio sample #)\" --ylabel \"I\" ",
        3000);
    if(debugPlots.eyeDiagramRealEnabled)
        debugPlots.eyeDiagramRealStdin = popen(eyeDiagramPlotReal, "w");

    char eyeDiagramPlotImaginary[300] = {0};
    sprintf(
        eyeDiagramPlotImaginary,
        "feedgnuplot "
        "--domain --dataid --lines --points --maxcurves %i "
        "--title \"Eye Diagram, Imaginary part\" "
        "--xlabel \"Time (Audio sample #)\" --ylabel \"Q\" ",
        3000);
    if(debugPlots.eyeDiagramImaginaryEnabled)
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
    if(debugPlots.IQvsTimeEnabled)
    {
        debugPlots.IQvsTimeStdin = popen(IQvstimeplot, "w");

        if (debugPlots.IQvsTimeStdin == NULL)
        {
            fprintf(stderr, "Failed to create IQ vs Time plot: %s\n", strerror(errno));
            retval = 7;
            goto exit;
        }
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
    if(debugPlots.filterDebugEnabled)
    {
        debugPlots.filterDebugStdin = popen(filterDebugPlot, "w");
        if(debugPlots.filterDebugStdin == NULL)
        {
            fprintf(stderr, "Failed to create filter debug plot: %s\n", strerror(errno));
            retval = 8;
            goto exit;
        }
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
    if(debugPlots.QAMdecoderEnabled)
    {
        debugPlots.QAMdecoderStdin = popen(QAMdemodulatePlot, "w");
        if(debugPlots.QAMdecoderStdin == NULL)
        {
            fprintf(stderr, "Failed to create QAM decoder debug plot: %s\n", strerror(errno));
            retval = 9;
            goto exit;
        }
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
    if(debugPlots.gardnerAlgoEnabled)
    {
        debugPlots.gardnerAlgoStdin = popen(gardenerAlgoPlot, "w");
        if(debugPlots.gardnerAlgoStdin == NULL)
        {
            fprintf(stderr, "Failed to create gardner algo debug plot: %s\n", strerror(errno));
            retval = 10;
            goto exit;
        }
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
    if(debugPlots.channelFilterEnabled)
    {
        debugPlots.channelFilterStdin = popen(channelFilterPlot, "w");
        if(debugPlots.channelFilterStdin == NULL)
        {
            fprintf(stderr, "Failed to create channel filter plot: %s\n", strerror(errno));
            retval = 11;
            goto exit;
        }
    }

    char *OFDMtimingSyncPlot =
        "feedgnuplot "
        "--domain --dataid --lines --points "
        "--title \"OFDM timing sync\" "
        "--xlabel \"sample number\" --ylabel \"value\" "
        "--legend 0 \"preamble auto correlation\" "
        "--legend 1 \"Channel simulated samples\" "
        "--legend 2 \"original samples\" "
        "--legend 6 \"OFDM Transmission found, at given offset\" "
        "--legend 7 \"Plateu Estimation Points\" "
        "--legend 8 \"symbol wide average\" "
        "--legend 9 \"lower threshold\" "
        "--legend 10 \"upper threshold\" "
        "--legend 11 \"Recieved Signal RMS\" "

    ;
    if(debugPlots.OFDMtimingSyncEnabled)
    {
        debugPlots.OFDMtimingSyncStdin = popen(OFDMtimingSyncPlot, "w");
        if(debugPlots.OFDMtimingSyncStdin == NULL)
        {
            fprintf(stderr, "Failed to create OFDM timing sync plot: %s\n", strerror(errno));
            retval = 11;
            goto exit;
        }
    }

    char *OFDMdecoderPlot =
        "feedgnuplot "
        "--domain --dataid --lines "
        "--title \"OFDM symbol decoder\" "
        "--xlabel \"sample number\" --ylabel \"value\" "
        "--legend 0 \"Preamble detection raw samples\" "
        "--legend 1 \"further raw samples\" "
    ;
    if(debugPlots.OFDMdecoderEnabled)
    {
        debugPlots.OFDMdecoderStdin = popen(OFDMdecoderPlot, "w");
        if(debugPlots.OFDMtimingSyncStdin == NULL)
        {
            fprintf(stderr, "Failed to create OFDM decoder plot: %s\n", strerror(errno));
            retval = 12;
            goto exit;
        }
    }

    char *OFDMIQPlot =
        "feedgnuplot "
        "--domain --dataid --points "
        "--title \"IQ plots for every subchannel\" "
        "--xlabel \"I\" --ylabel \"Q\" "
        "--legend 0 \"Data samples\" "
        "--legend 1 \"Preamble 1 samples: autocorrelation\" "
        "--legend 2 \"Preamble 2,3 samples: sample rate estimator\" "
        "--legend 3 \"pre lockin Data samples\" "

    ;
    if(debugPlots.OFDMrawIQEnabled)
    {
        debugPlots.OFDMrawIQStdin = popen(OFDMIQPlot, "w");
        //debugPlots.OFDMrawIQStdin = fopen("IQrawPlotOutput", "w");
        if(debugPlots.OFDMrawIQStdin == NULL)
        {
            fprintf(stderr, "Failed to create OFDM IQ plot: %s\n", strerror(errno));
            retval = 13;
            goto exit;
        }
    }

    char *OFDMinterpolatorPlot =
        "feedgnuplot "
        "--domain --dataid --lines --points "
        "--title \"OFDM sample interpolator\" "
        "--xlabel \"sample time (inputSampleTimed)\" --ylabel \"value\" "
        "--legend 0 \"Input Samples\" "
        "--legend 1 \"Output Samples\" "
    ;
    if(debugPlots.OFDMinterpolatorEnabled)
    {
        debugPlots.OFDMinterpolatorStdin= popen(OFDMinterpolatorPlot, "w");
        if(debugPlots.OFDMinterpolatorStdin == NULL)
        {
            fprintf(stderr, "Failed to create OFDM timing sync plot: %s\n", strerror(errno));
            retval = 11;
            goto exit;
        }
    }

    char *OFDMsfoEstimatorPlot =
        "feedgnuplot "
        "--domain --dataid --lines --points "
        "--title \"OFDM sampling frequency offset estimator\" "
        "--xlabel \"OFDM preamble symbol number\" --ylabel \"estimated frequency error (Hz ppm)\" "
        "--legend 0 \"Estimated Offset (used for resampler)\" "
        "--legend 1 \"Detected Estimated Offset\" "
        "--legend 2 \"SFO Estimation actually used for the given samples\" "
        "--legend 3 \"residual Offset (used for phase correction only)\" "
        "--legend 4 \"Detected residual Offset\" "
        "--legend 5 \"Total SFO correction (resampler and phase residual phase correction combined)\" "
    ;
    if(debugPlots.OFDMsfoEstimatorEnabled)
    {
        debugPlots.OFDMsfoEstimatorStdin = popen(OFDMsfoEstimatorPlot, "w");
        if(debugPlots.OFDMsfoEstimatorStdin == NULL)
        {
            fprintf(stderr, "Failed to create OFDM SFO estimator plot: %s\n", strerror(errno));
            retval = 11;
            goto exit;
        }
    }

    char *OFDMEqualizerErrorIQPlot =
        "feedgnuplot "
        "--domain --dataid --points --lines "
        "--maxcurves 9000 "
        "--title \"IQ plots for every subchannel's equalizer error\" "
        "--xlabel \"I\" --ylabel \"Q\" "
        "--legend 0 \"samples\" "

    ;
    if(debugPlots.OFDMequalizerErrorIQEnabled)
    {
        debugPlots.OFDMequalizerErrorIQStdin = popen(OFDMEqualizerErrorIQPlot, "w");
        //debugPlots.OFDMrawIQStdin = fopen("IQrawPlotOutput", "w");
        if(debugPlots.OFDMequalizerErrorIQStdin == NULL)
        {
            fprintf(stderr, "Failed to create OFDM equalizer IQ plot: %s\n", strerror(errno));
            retval = 14;
            goto exit;
        }
    }

    char *OFDMEqualizerIQPlot =
        "feedgnuplot "
        "--domain --dataid --points --lines "
        "--maxcurves 9000 "
        "--title \"IQ plots for every subchannel's equalizer parameter\" "
        "--xlabel \"I\" --ylabel \"Q\" "
        "--legend 0 \"samples\" "

    ;
    if(debugPlots.OFDMequalizerIQEnabled)
    {
        debugPlots.OFDMequalizerIQStdin = popen(OFDMEqualizerIQPlot, "w");
        //debugPlots.OFDMrawIQStdin = fopen("IQrawPlotOutput", "w");
        if(debugPlots.OFDMequalizerIQStdin == NULL)
        {
            fprintf(stderr, "Failed to create OFDM equalizer IQ plot: %s\n", strerror(errno));
            retval = 14;
            goto exit;
        }
    }

    char *OFDMdataIQPlot =
        "feedgnuplot "
        "--domain --dataid --points "
        "--title \"equalizer corrected IQ plots for every subchannel\" "
        "--xlabel \"I\" --ylabel \"Q\" "
        "--legend 0 \"Data samples\" "

    ;
    if(debugPlots.OFDMdataIQEnabled)
    {
        debugPlots.OFDMdataIQStdin = popen(OFDMdataIQPlot, "w");
        //debugPlots.OFDMrawIQStdin = fopen("IQrawPlotOutput", "w");
        if(debugPlots.OFDMdataIQStdin == NULL)
        {
            fprintf(stderr, "Failed to create OFDM data IQ plot: %s\n", strerror(errno));
            retval = 15;
            goto exit;
        }
    }

    // I might want to slightly over sample the signal for the benefit of interpolation.
    // basically this is taking advantage of ffmpeg's interpolation filters, simplifying my program's architecture
    // I could internalize the call to ffmpeg, or write my own filter that can handle increasing sample rates rather than just decreasing sample rates
    //     rates [0x560]: 44100 48000 96000 192000
    //int sampleRate = 44100;
    int sampleRate = 192000;
    //int upconvertionSampleRate = 192000;

    OFDM_state_t OFDMstate = {0};
    initializeOFDMstate(&OFDMstate);
    OFDMstate.sampleRate = 44100;   // set the main sample rate for the decoder, but the incomming samples could have a higher rate

    // while there is data to recieve, not end of file -> right now just a fixed number of 2000
    //for(int audioSampleIndex = 0; audioSampleIndex < SYMBOL_PERIOD * 600; audioSampleIndex++)
    for(int audioSampleIndex = 0; audioSampleIndex < sampleRate * OFDMstate.duration; audioSampleIndex++)
    //for(int audioSampleIndex = 0; audioSampleIndex < SYMBOL_PERIOD * 2000; audioSampleIndex++)
    {
        // recieve data on stdin, signed 32bit integer
        for(size_t i = 0; i < sizeof(sampleConvert.value); i++)
        {
            // get the bytes, little endian, of signed 32bit integer, using the struct to do type punning
            sampleConvert.bytes[i] = getchar();
        }


        // convert to a double ranged from -1 to 1
        sample.sample = (double)sampleConvert.value / INT32_MAX;
        sample.sampleRate = sampleRate;
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

        demodualteOFDM(&sample, &OFDMstate, debugPlots);


        //OFDM_state_t OFDMdemodulateState = {SYMBOL_PERIOD, &sampleBuffer, k};
        //demodulateOFDM(n, sampleDouble, &OFDMdemodulateState, &debugPlots);

    }


exit:
    if (iq_plot_buffer != NULL)
    {
        free(iq_plot_buffer);
        iq_plot_buffer = NULL;
    }

    // close the streams, allowing the plots to render and open a window, then wait for them to terminate in separate threads.
    if((debugPlots.IQvsTimeStdin != NULL) && (fork() == 0))
    {
        pclose(debugPlots.IQvsTimeStdin);
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
    if((debugPlots.OFDMinterpolatorStdin != NULL) && (fork() == 0))
    {
        pclose(debugPlots.OFDMinterpolatorStdin);
        return 0;
    }

    // after all the data is recieved, generate a plot of IQ
    return retval;
}
