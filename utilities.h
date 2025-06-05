#include <stdlib.h>
#include <stdint.h>
#include <fftw3.h>


// enums
//

typedef enum
{
    RIFF_TYPE_PCM = 1,
} riff_type_t;

typedef enum
{
    RETURNED_SAMPLE,
    AWAITING_SAMPLES,
} buffered_data_return_t;

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

// used for constellation constellations
typedef enum
{
    OPHIUCHUS,
    URSA_MINOR,
    ARA,
    LUPUS,
    LEO,
    URSA_MAJOR,
    PUPPIS,
    ANTLIA,
    SAGITTARIUS,
    LIBRA,
    CIRCINUS,
    CYGNUS,
    SCORPIO,
    CAMELOPARDALIS,
    GEMINI,
    CANIS_MINOR,
    CRUX,
    BOOTES,
    CARINA,
    TAURUS,
    CANCER,
    CASSIOPEIA,
    LYRA,
    AQUARIUS,
    CENTAURUS,
    PISCES,
    VIRGO,
    APUS,
    ORION,
    CAELUM,
    ARIES,
    AQUILA,
    CANIS_MAJOR,
    AURIGA,
    CAPRICORN,
    MONOCEROS,
    CANES_VENATICI,
    ANDROMEDA,
    CORONA_BOREALIS,
    CONSTELLATION_MAX,
} constellation;


// typedef structs
//

typedef struct
{
    FILE* waveformPlotStdin;
    FILE* fftDebuggerStdin;
    FILE* errorPlotStdin;
    FILE* IQplotStdin;
    FILE* IQvsTimeStdin;
    FILE* eyeDiagramRealStdin;
    FILE* eyeDiagramImaginaryStdin;
    FILE* filterDebugStdin;
    FILE* QAMdecoderStdin;
    FILE* gardnerAlgoStdin;
    FILE* OFDMtimingSyncStdin;
    FILE* channelFilterStdin;
    FILE* OFDMdecoderStdin;
    FILE* OFDMrawIQStdin;
    FILE* OFDMinterpolatorStdin;
    FILE* OFDMsfoEstimatorStdin;
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
            unsigned int OFDMrawIQEnabled           : 1;
            unsigned int OFDMinterpolatorEnabled    : 1;
            unsigned int OFDMsfoEstimatorEnabled    : 1;

        };
            unsigned long int flags;
    };
} debugPlots_t;


typedef struct __attribute__((packed)) riff_header
{
    union
    {
        struct
        {
            char riff[4];               // 0
            uint32_t size;              // 4
            char format[4];             // 8
            char chunk[4];              // 12
            uint32_t length;            // 16
            uint16_t type;              // 20
            uint16_t channels;          // 22
            uint32_t sampleRate;        // 24
            uint32_t dataRate;          // 28
            uint16_t blockSize;         // 32
            uint16_t bitsPerSample;     // 34
            char data[4];               // 36
            uint32_t chunkSize;         // 40
        };
        uint8_t bytes[44];
    };
} riff_header_t;

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
    _Complex double sample;      // double value representation of sample
    int sampleRate;     // rate of samples per second
    int sampleIndex;    // index since begining of record
} sample_complex_t;

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
    _Complex double *buffer;// a pointer to an array for storing samples
    int length;             // the length of the buffer
    int insertionIndex;     // the index where new data is inserted into the circular buffer
    int phase;              // used for some functions choosing when to do periodic calculations with the buffer
    int sampleRate;
    int n;                  // index of sample at insertionIndex
} circular_buffer_complex_t;

typedef struct
{
    circular_buffer_complex_t frequencyDomain;              // holds the incoming ofdmPeriod number of samples
    circular_buffer_double_t timeDomain;             // the decoded symbol
} fftw_OFDM_buffer_t;

typedef struct
{
    double *M_buffer;   // start of the block
    double *L_buffer;   // offset M units into the block
    double *buffer;     // same as L_buffer
    double *M_next_buffer;  // offset (M + L) - M = L into the block (M units from the end)
    int insertionIndex;     // insertion index for circular buffer starting at L_buffer with length L
    int length;             // length of circular buffer L length of L
    int L;                  // same as length
    int M;                  // length of M_buffer
    int N;                  // sum L + M
    int n;                  // time index of the last sample inserted in L
    int sampleRate;
} overlap_save_buffer_double_t;

// this is a bit lame. should use complex datatype
typedef struct
{
    double InPhase;
    double Quadrature;
} iqsample_t;

typedef struct __attribute__((packed))
{
    union
    {
        int32_t value;
        uint8_t bytes[sizeof(int32_t)];
    };
} int32_to_bytes_t;

// a structure for splitting out one bit at a time from a byte
// for implementing a hoffman tree code scheme
typedef struct
{
    uint8_t *bytes;
    int insertionIndex; // index ehere the newest data is in the buffer
    int readIndex;      // index of the current byte being read
    int bitIndex;       // index of the bit within the current byte being used
    uint8_t nextBit;
} circular_buffer_byte_t;

typedef struct huffman_tree_t
{
    char isLeaf;    // is it a leaf with an IQ value, or does it have children

    struct huffman_tree_t *child[2];    // the two children

    int constellationIndex;   // the index of the constellation point 

} huffman_tree_t;

typedef struct
{
    _Complex double *points;   // a pointer to an array of points in the constellation
    int length;               // the number of points in the constellation

    huffman_tree_t* huffmanTree;    // a binary tree used for encoding a bit stream into constellation points

} constellation_complex_t;

typedef struct
{
    double carrierFrequency;
    double carrierPhase;            // probably not going to use in the end.
    _Complex double IQsamplingTransform;         // used to transform the IQ samples to compensate for carrier phase and amplitude mismatches (ln(amp) + i*phase)
    double k;  // number of carrier cycles per sample
    double symbolPeriod;    // number of audio samples per symbol
    double selectedIQsamples[4];    // the closest samples to ideal sample time and mid symbol sample time for gardner algorithm
} QAM_properties_t;

typedef struct
{
    // Common between sender and reciever
    //
    int sampleRate;
    int channels;

    // carrier frequenncy?
    // if 1 than we're working with real samples, not complex
    int baseband;   // if set to 1, IQ sampels are taken as real so OFDM symbols can be transmitted at baseband. Halves the number of channels. to even channels only

    int guardPeriod;
    int ofdmPeriod;     // just the OFDM symbol itslef
    int symbolPeriod;   // includes guard period and ODFM period, sum of guard and ofdm periods

    int initialized;

    int pilotSymbolsPitch;  // spacing from one pilot subchannel to the next
    _Complex double* pilotSymbols;   // an array for storing one time generated pilot symbols
    int pilotSymbolsLength; // length of the pilotSymbols array

    // fftw constructs
    fftw_plan fftwPlan;

    // random number generators
    struct drand48_data preamblePilotsPRNG;
    struct drand48_data pilotsPRNG;
    struct drand48_data predefinedDataPRNG;
    struct drand48_data channelNoisePRNG;

    // Sender relevant
    //

    FILE* dataInput;    // pipe to read data from that will be sent of the channel
    char inputByte;     // the current working byte
    int bitOffset;      // hold the offset in bits from the beginning of the current byte pulled from the dataInput file
    FILE* generatedDataOutput;  // print the randomly generated data that was encoded in a file

    // array length of channels representing the current entire OFDM symbol
    //double complex *currentOFDMSymbol;
    
    // channel simulation stuff
    int simulateNoise;
    int simulateChannel;    // 0 don't simulate, 1 simulate
    overlap_save_buffer_double_t channelSimulationBuffer;   // holds samples to be used for simulating channel impulse response
    circular_buffer_double_t channelImpulseResponse;    // holds the impulse response of the channel.

    fftw_OFDM_buffer_t OFDMsymbol;
    //circular_buffer_double_t fftwSymbolBuffer;              // holds the incoming ofdmPeriod number of samples
    //circular_buffer_complex_t fftwIQbuffer;             // the decoded symbol

    fftw_OFDM_buffer_t IQrateDetectorFirstHalf;
    fftw_OFDM_buffer_t IQrateDetectorSecondHalf;
    //circular_buffer_complex_t fftwIQrateDetectorFirstHalf;          // for holding a half ofdmPeriod DFT used in sample rate error detection
    //circular_buffer_complex_t fftwIQrateDetectorSecondHalf;          // for holding a half ofdmPeriod DFT used in sample rate error detection


    // constellations array
    constellation_complex_t *constellations;
    int constellationsLength;

    struct
    {
        int long frameStart;    // estimate of the begining of the current frame of given type
        int fieldIndex;     // index of the field within the frame
        enum
        {
            IDLE,
            SEARCHING,
            ACTIVE,
        } frame;    // a series of fields

        int long fieldStart;    // estimate of teh beginning of the current field of given type
        int symbolIndex;     // index of symbol in the current field

        enum
        {
            // Active field types

            // estimating Sampling Frequency Offset (SFO) and Channel (amplitude phase per subchannel)
            PREAMBLE,   // for frame detection and timing synchronization 
            DATA,       // data field, also contains pilot symbols interspersed

        } field;    // a series of OFDM symbols

        int long symbolStart;
        int long processedSymbols;
        enum
        {
            // all otehr symbol states
            GUARD_PERIOD,   // cyclic prefix for obsorbing channel effects (ISI)
            OFDM_PERIOD,  // transmission of the actual OFDM symbol information
            

        } symbol;   // an atomic OFDM symbol
    } state;

    // Reciever relevant
    //

    circular_buffer_double_t preambleDetectorInputBuffer; // holds the last ofdmPeriod plus a little of samples for preamble detection

    fftw_plan fftwRateDetectorFirstHalfPlan;
    fftw_plan fftwRateDetectorSecondHalfPlan;
    int ofdmPhaseOffset;                        // stores the offset in sample number to start taking cutting out OFDM symbols

    double samplingFrequencyOffsetEstimate; // the estimated frequency offset in hertz
    double samplingFrequencyOffsetResidual; // remaining frequency offset after resampling corrections, used for phase corrections
    double resamplingRatio; // ratio to multiply by the incoming samples' rates by to time the interpolations.
    int long sample;    // the current sample number, used by the sample interpolator
    double samplerAccumulatedPhase; // the time of the sample at index sample
    circular_buffer_double_t sampleInterpolatorBuffer;  // for resampling

    int disableSFOestimation;   // flag to disable any sfo corrections
    circular_buffer_complex_t sfoFirstSymbol;   // sampling frequency offset detector buffer to hold the previous IQ samples

    _Complex double* channelEstimate;  // equalizer parameters

    sample_double_t autoCorrelation;        // schmidl timing signal
    sample_double_t autoCorrelationDerivative;
    circular_buffer_double_t autoCorrelationAverageBuffer;
    circular_buffer_double_t autoCorrelationDerivativeAverageBuffer;           // avarage for the derivative thing
    
    circular_buffer_double_t timingFilterInputBuffer;   // buffer of averaged autocorrelation values for timing filter

    sample_double_t preambleEnergy; // for normalizing the timing signal and automatic gain control

} OFDM_state_t;


// functions
//

void initializeOFDMstate(OFDM_state_t*);
//void initializeCircularBuffer_complex(circular_buffer_complex_t*, int);
void initializeCircularBuffer_fftw_complex(circular_buffer_complex_t*, int, int);
void initializeCircularBuffer_complex(circular_buffer_complex_t*, int, int);
void initializeCircularBuffer_double(circular_buffer_double_t*, int, int);
void initializeOverlapAndSaveBuffer(overlap_save_buffer_double_t*, int);
void generateHuffmanTree(constellation_complex_t*);
_Complex double traverseHuffmanTree(OFDM_state_t*, constellation_complex_t*);

// some functions to generate IQ streams with different properties
iqsample_t impulse(int, int);
iqsample_t alternateI(int);
iqsample_t alternateQ(int);
iqsample_t randomQAM(int);
iqsample_t randomQAM_withPreamble(int, int);
iqsample_t sequentialIQ(int, int);

// some IQ generators for the OFDM implimentation
// generate a linear constellation on the I axis of 'levels' number of points
void fftw_ASK(int, constellation_complex_t*);
void fftw_PSK(int, constellation_complex_t*);
// generate a square QAM constelation grid with side length 'square'
//  for example, 16QAM would be square=4
void fftw_squareQAM(int, constellation_complex_t*);
// a hexagonal array of points with the length of two adjacent sides as edgeLength
void fftw_hexQAM(int, constellation_complex_t*);

// uses the stars of orion for the constellation
// coordinates pulled from starfetch json files in ~/Games
void fftw_starQAM(constellation, constellation_complex_t*);
