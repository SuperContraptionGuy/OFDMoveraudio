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

typedef struct sample_double
{
    double sample;      // double value representation of sample
    int sampleRate;     // rate of samples per second
    int sampleIndex;    // index since begining of record
} sample_double_t;

typedef struct sample_complex
{
    _Complex double sample;      // double value representation of sample
    int sampleRate;     // rate of samples per second
    int sampleIndex;    // index since begining of record
} sample_complex_t;

typedef struct circular_buffer_double
{
    double *buffer;        // a pointer to an array for storing samples
    int length;             // the length of the buffer
    int insertionIndex;     // the index where new data is inserted into the circular buffer
    int phase;              // used for some functions choosing when to do periodic calculations with the buffer
    int sampleRate;
    int n;
} circular_buffer_double_t;

typedef struct circular_buffer
{
    _Complex double *buffer;// a pointer to an array for storing samples
    int length;             // the length of the buffer
    int insertionIndex;     // the index where new data is inserted into the circular buffer
    int phase;              // used for some functions choosing when to do periodic calculations with the buffer
    int sampleRate;
    int n;                  // index of sample at insertionIndex
} circular_buffer_complex_t;

typedef struct overlap_save_buffer_double
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
typedef struct iqsample
{
    double InPhase;
    double Quadrature;
} iqsample_t;

typedef struct __attribute__((packed)) int32_to_bytes
{
    union
    {
        int32_t value;
        uint8_t bytes[sizeof(int32_t)];
    };
} int32_to_bytes_t;

// a structure for splitting out one bit at a time from a byte
// for implementing a hoffman tree code scheme
typedef struct circular_buffer_byte
{
    uint8_t *bytes;
    int insertionIndex; // index ehere the newest data is in the buffer
    int readIndex;      // index of the current byte being read
    int bitIndex;       // index of the bit within the current byte being used
    uint8_t nextBit;
} circular_buffer_byte_t;

typedef struct huffman_tree
{
    char isLeaf;    // is it a leaf with an IQ value, or does it have children

    struct huffman_tree *child[2];    // the two children

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
    int sampleRate;
    int channels;

    // carrier frequenncy?
    int baseband;   // if set to 1, IQ sampels are taken as real so OFDM symbols can be transmitted at baseband. Halves the number of channels. to even channels only

    int guardPeriod;
    int symbolPeriod;   // includes guard period and ODFM period
    int ofdmPeriod;     // just the OFDM symbol itslef

    int initialized;

    FILE* dataInput;    // pipe to read data from that will be sent of the channel
    char inputByte;     // the current working byte
    int bitOffset;      // hold the offset in bits from the beginning of the current byte pulled from the dataInput file
    FILE* generatedDataOutput;  // print the randomly generated data that was encoded in a file

    int pilotSymbolsPitch;  // spacing from one pilot subchannel to the next
    _Complex double* pilotSymbols;   // an array for storing one time generated pilot symbols

    // array length of channels representing the current entire OFDM symbol
    //double complex *currentOFDMSymbol;
    
    // channel simulation stuff
    int simulateNoise;
    int simulateChannel;    // 0 don't simulate, 1 simulate
    overlap_save_buffer_double_t channelSimulationBuffer;   // holds samples to be used for simulating channel impulse response
    circular_buffer_double_t channelImpulseResponse;    // holds the impulse response of the channel.

    // fftw constructs
    fftw_plan fftwPlan;
    struct
    {
        fftw_complex    *frequencyDomain;   // array of size two doubel arrays  size floor(ofdmPeriod / 2) + 1.     FYI this input array will be clobbered by fftw_execute calls
        double          *timeDomain;        // array of doubles size ofdmPeriod.
    } OFDMsymbol;

    // constellations array
    constellation_complex_t *constellations;
    int constellationsLength;

    struct
    {
        int long frameStart;
        int fieldIndex;     // index of the field within the frame
        enum
        {
            IDLE,
            ACTIVE,
        } frame;    // a series of fields

        int long fieldStart;
        int symbolIndex;    // index of the symbol within the field
        enum
        {
            // Active field types
            PREAMBLE,   // for frame detection and timing synchronization 
            DATA,       // data field, also contains pilot symbols interspersed

        } field;    // a series of OFDM symbols

        int long symbolStart;
        enum
        {
            // all otehr symbol states
            GUARD_PERIOD,   // cyclic prefix for obsorbing channel effects (ISI)
            OFDM_PERIOD,  // transmission of the actual OFDM symbol information
            

        } symbol;   // an atomic OFDM symbol


    } state;
} OFDM_state_t;


// functions
//

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
