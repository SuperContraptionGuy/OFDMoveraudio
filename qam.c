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

#define WARN_UNUSED __attribute__((warn_unused_result))

#define DEBUG_LEVEL 0
#if DEBUG_LEVEL >= 1
    FILE* hexdumpStdIn = NULL;
    FILE* plotStdIn = NULL;
#endif

typedef enum
{
    RIFF_TYPE_PCM = 1,
} riff_type_t;

typedef struct __attribute__((packed))
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

typedef enum
{
    RETURNED_SAMPLE,
    AWAITING_SAMPLES,
} buffered_data_return_t;

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

typedef struct
{
    fftw_complex *points;   // a pointer to an array of points in the constellation
    int length;             // the number of points in the constellation

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



// some functions to generate IQ streams with different properties
iqsample_t impulse(int symbolIndex, int impulseTime)
{
    iqsample_t sample = {0, 0};
    if(symbolIndex == impulseTime)
        sample.InPhase = 1;
    return sample;
}

iqsample_t alternateI(int symbolIndex)
{
    iqsample_t sample = {symbolIndex % 2 * 2 - 1, 0};
    //iqsample_t sample = {symbolIndex % 2, 0};
    return sample;
}

iqsample_t alternateQ(int symbolIndex)
{
    iqsample_t sample = {0, symbolIndex % 2 * 2 - 1};
    //iqsample_t sample = {symbolIndex % 2, 0};
    return sample;
}

iqsample_t randomQAM(int square)
{
    // square is the number of states of I and of Q, total states is square squared
    // I and Q values returned between -1 and 1
    iqsample_t sample = {0, 0};
    sample.InPhase = (double)(rand() % square) / (square - 1) * 2 - 1;
    sample.Quadrature = (double)(rand() % square) / (square - 1) * 2 - 1;
    return sample;
}

iqsample_t randomQAM_withPreamble(int symbolIndex, int square)
{
    iqsample_t sample = {0, 0};
    srand((unsigned int)symbolIndex); // to get a stable random sequence everytime a specific symbol is requested
    int preamble = 150; // length of preamble in symbols
    if (symbolIndex < preamble)
    {
        // add a preamble that's easy to get rough time sync to
        sample = alternateI(symbolIndex);
    } else { // choose a new random QAM IQ value at start of every total period
        // then start sending random data
        sample = randomQAM(square);
    }
    return sample;
}

iqsample_t sequentialIQ(int symbolIndex, int square)
{
    iqsample_t symbol;
    // sequentially hit all the IQ values in order in the constelation defined by power
    symbol.InPhase = (double)(symbolIndex % square) / (square - 1) * 2 - 1;
    symbol.Quadrature = (double)(symbolIndex / square % square) / (square - 1) * 2 - 1;
    return symbol;
}

// some IQ generators for the OFDM implimentation

// generate a linear constellation on the I axis of 'levels' number of points
void fftw_ASK(int levels, constellation_complex_t *constellation)
{
    if(levels < 2)
        levels = 1;

    constellation->points = malloc(sizeof(*constellation->points) * levels);
    constellation->length = levels;

    if(levels == 1)
    {
        constellation->points[0] = 0;
        return;
    }
    double rms = 0;
    for(int i = 0; i < levels; i++)
    {
        constellation->points[i] = (double)i / (levels - 1) * 2.0 - 1;
        rms += cabs(constellation->points[i]);
    }
    // calculate the RMS
    rms = rms == 0 ? 1 : rms / levels * M_SQRT1_2;
    for(int i = 0; i < levels; i++)
    {
        // correct the RMS power to 1
        constellation->points[i] /= rms;
    }
    return;
}

void fftw_PSK(int levels, constellation_complex_t *constellation)
{
    if(levels < 2)
        levels = 1;

    constellation->points = malloc(sizeof(*constellation->points) * levels);
    constellation->length = levels;

    if(levels == 1)
    {
        constellation->points[0] = 0;
        return;
    }

    double angle = 2*M_PI / levels;

    double rms = 0;
    for(int i = 0; i < levels; i++)
    {
        constellation->points[i] = cos(angle * i) + sin(angle * i)*I;
        rms += cabs(constellation->points[i]);
    }
    rms = rms == 0 ? 1 : rms / levels * M_SQRT1_2;
    for(int i = 0; i < levels; i++)
        constellation->points[i] /= rms;
    return;
}

// generate a square QAM constelation grid with side length 'square'
//  for example, 16QAM would be square=4
void fftw_squareQAM(int square, constellation_complex_t *constellation)
{

    if(square < 2)
        square = 1;


    int levels = pow(square, 2);
    constellation->points = malloc(sizeof(fftw_complex) * levels);
    constellation->length = levels;

    if(levels == 1)
    {
        constellation->points[0] = 0;
        return;
    }


    double rms = 0;
    for(int i = 0; i < levels; i++)
    {
        constellation->points[i] = 
            ((double)(i % square) / (square - 1) * 2. - 1. +
            I*((double)(i / square) / (square - 1) * 2. - 1));     // random (square)QAM, square^2 constellation points
        rms += cabs(constellation->points[i]);
    }
    rms = rms == 0 ? 1 : rms / levels * M_SQRT1_2;
    for(int i = 0; i < levels; i++)
        constellation->points[i] /= rms;
    return;
}

// a hexagonal array of points with the length of two adjacent sides as edgeLength
void fftw_hexQAM(int gridSize, constellation_complex_t *constellation)
{
    if(gridSize < 1)
        gridSize = 1;
    // the points will be layed out with three basis vectors over three two dimensional parallelagram shaped fields

    int edgeLength = ceil((double)gridSize / 2);    // convert length of two sides to length of one side.
    // determin if the pattern has a constellation point at the origin or not. also changes some offsets and number of points to draw
    int hasCenter = 0;
    int n_field;    // number of points in each field depends on if the grid is centered or offset
    if(gridSize % 2 == 1)
    //if(1)
    {
        hasCenter = 1;
        // total number of points is the number of points in one parallelogram field (not counting one connecting edge) times 3 plus the center point
        n_field = (edgeLength - 1) * edgeLength;
    } else {
        n_field = edgeLength * edgeLength;
    }

    int n_total = n_field * 3;
    if(hasCenter)
        n_total += 1;

    // allocate array for constellation points
    constellation->points = malloc(sizeof(*constellation->points) * n_total);
    constellation->length = n_total;

    double rms = 0;
    for(int index = 0; index < n_total; index++)
    {

        if(index == n_total - 1 && hasCenter) // the last value
        {
            constellation->points[index] = 0;   // return the center point if there is one and the index hit it
        } else {

            // the index of the field to use
            int field = index / n_field;
            // the index of the point within a field
            int subIndex = index % n_field;
            //int x = subIndex / edgeLength;
            int x = subIndex / edgeLength;
            if(hasCenter)
                x += 1; // skip one row by starting at 1 instead of 0 if there is a center
            int y = subIndex % edgeLength;

            // the angles of the two basis vectors to use
            double angle1 = 2 * M_PI / 3 * field;
            double angle2 = angle1 + 2 * M_PI / 3;
            // the two basis vectors to use for constructing a field
            fftw_complex basis1 = cos(angle1) + sin(angle1)*I;
            fftw_complex basis2 = cos(angle2) + sin(angle2)*I;

            // the angle of the offset vector
            double offsetAngle = angle1 + 2*M_PI/12;
            double offsetLength = sin(2*M_PI/12) / sin(2*M_PI/3);
            fftw_complex offsetVector = offsetLength*cos(offsetAngle) + offsetLength*sin(offsetAngle)*I;

            // calculate point in the field's grid
            fftw_complex point = (basis1 * x + basis2 * y);
            //double normalizationFactor = (double)gridSize / 2;
            if(hasCenter)
                constellation->points[index] = point / (edgeLength - 1);
            else
                // normalization factor is a bit different for a non-centered hex grid
                constellation->points[index] = (point + offsetVector) / sqrt(pow(edgeLength-1, 2) + pow(offsetLength, 2) - 2*(edgeLength-1)*(offsetLength)*cos(M_PI - 2*M_PI/12));
        }

        rms += cabs(constellation->points[index]);
    }
    rms = rms == 0 ? 1 : rms / n_total * M_SQRT1_2;
    for(int i = 0; i < n_total; i++)
        constellation->points[i] /= rms;
    return;
}

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

// uses the stars of orion for the constellation
// coordinates pulled from starfetch json files in ~/Games
void fftw_starQAM(constellation name, constellation_complex_t *constellation)
{
    const fftw_complex ophiuchus[10] =
    {
        7  + 2  * I,
        12 + 3  * I,
        16 + 4  * I,
        6  + 5  * I,
        18 + 5  * I,
        5  + 6  * I,
        17 + 6  * I,
        15 + 8  * I,
        4  + 9  * I,
        12 + 9  * I,
    };

    const fftw_complex ursa_minor[7] =
    {
        13 + 2  * I,
        11 + 3  * I,
        10 + 5  * I,
        11 + 7  * I,
        8  + 8  * I,
        14 + 9  * I,
        11 + 10 * I,
    };

    const fftw_complex ara[6] =
    {
        4  + 2  * I,
        16 + 3  * I,
        16 + 6  * I,
        9  + 6  * I,
        19 + 9  * I,
        6  + 10 * I,
    };

    const fftw_complex lupus[9] =
    {
        6  + 1  * I,
        11 + 2  * I,
        2  + 3  * I,
        8  + 4  * I,
        12 + 4  * I,
        16 + 5  * I,
        11 + 6  * I,
        17 + 8  * I,
        12 + 10 * I,
    };

    const fftw_complex leo[9] =
    {
        17 + 2  * I,
        19 + 3  * I,
        14 + 4  * I,
        14 + 5  * I,
        7  + 6  * I,
        17 + 6  * I,
        8  + 8  * I,
        18 + 8  * I,
        3  + 9  * I,
    };

    const fftw_complex ursa_major[7] =
    {
        5  + 2  * I,
        8  + 3  * I,
        9  + 4  * I,
        11 + 6  * I,
        9  + 8  * I,
        17 + 9  * I,
        13 + 10 * I,
    };

    const fftw_complex puppis[9] =
    {
        7  + 1  * I,
        5  + 2  * I,
        10 + 2  * I,
        11 + 3  * I,
        14 + 6  * I,
        6  + 7  * I,
        12 + 8  * I,
        18 + 8  * I,
        16 + 10 * I,
    };

    const fftw_complex antlia[4] =
    {
        15 + 3  * I,
        8  + 5  * I,
        19 + 6  * I,
        4  + 8  * I,
    };

    const fftw_complex sagittarius[8] =
    {
        12 + 3  * I,
        6  + 4  * I,
        8  + 4  * I,
        14 + 5  * I,
        4  + 6  * I,
        17 + 6  * I,
        6  + 7  * I,
        13 + 9  * I,
    };

    const fftw_complex libra[7] =
    {
        12 + 2  * I,
        18 + 3  * I,
        7  + 5  * I,
        4  + 6  * I,
        17 + 6  * I,
        12 + 8  * I,
        13 + 9  * I,
    };

    const fftw_complex circinus[5] =
    {
        8  + 2  * I,
        12 + 2  * I,
        11 + 5  * I,
        11 + 7  * I,
        12 + 9  * I,
    };

    const fftw_complex cygnus[8] =
    {
        18 + 2  * I,
        9  + 3  * I,
        15 + 4  * I,
        11 + 5  * I,
        7  + 7  * I,
        15 + 7  * I,
        4  + 8  * I,
        18 + 9  * I,
    };

    const fftw_complex scorpio[16] =
    {
        17 + 1  * I,
        18 + 2  * I,
        15 + 3  * I,
        14 + 4  * I,
        18 + 4  * I,
        13 + 5  * I,
        18 + 5  * I,
        5  + 7  * I,
        6  + 7  * I,
        11 + 7  * I,
        4  + 8  * I,
        10 + 8  * I,
        3  + 9  * I,
        4  + 10 * I,
        7  + 10 * I,
        9  + 10 * I,
    };

    const fftw_complex camelopardalis[10] =
    {
        1  + 1  * I,
        5  + 2  * I,
        9  + 3  * I,
        12 + 5  * I,
        6  + 6  * I,
        17 + 7  * I,
        7  + 8  * I,
        18 + 8  * I,
        9  + 10 * I,
        20 + 10 * I,
    };

    const fftw_complex gemini[16] =
    {
        7  + 1  * I,
        7  + 2  * I,
        14 + 2  * I,
        4  + 3  * I,
        11 + 3  * I,
        15 + 3  * I,
        19 + 3  * I,
        6  + 5  * I,
        4  + 6  * I,
        8  + 6  * I,
        15 + 6  * I,
        4  + 9  * I,
        8  + 9  * I,
        14 + 9  * I,
        17 + 9  * I,
        18 + 10 * I,
    };

    const fftw_complex canis_minor[2] =
    {
        18 + 4  * I,
        3  + 5  * I,
    };

    const fftw_complex crux[4] =
    {
        10 + 2  * I,
        16 + 4  * I,
        5  + 5  * I,
        13 + 9  * I,
    };

    const fftw_complex bootes[10] =
    {
        7  + 2  * I,
        12 + 3  * I,
        4  + 4  * I,
        12 + 5  * I,
        8  + 6  * I,
        13 + 8  * I,
        8  + 9  * I,
        17 + 9  * I,
        7  + 10 * I,
        18 + 10 * I,
    };

    const fftw_complex carina[9] =
    {
        12 + 1  * I,
        3  + 2  * I,
        8  + 3  * I,
        5  + 5  * I,
        19 + 5  * I,
        12 + 7  * I,
        19 + 8  * I,
        7  + 9  * I,
        11 + 10 * I,
    };

    const fftw_complex taurus[11] =
    {
        11 + 2  * I,
        5  + 3  * I,
        12 + 4  * I,
        8  + 5  * I,
        10 + 6  * I,
        11 + 6  * I,
        14 + 8  * I,
        9  + 9  * I,
        19 + 9  * I,
        10 + 10 * I,
        19 + 10 * I,
    };

    const fftw_complex cancer[5] =
    {
        7  + 2  * I,
        10 + 4  * I,
        11 + 5  * I,
        17 + 7  * I,
        8  + 8  * I,
    };

    const fftw_complex cassiopeia[5] =
    {
        4  + 3  * I,
        7  + 5  * I,
        12 + 5  * I,
        20 + 6  * I,
        16 + 8  * I,
    };

    const fftw_complex lyra[6] =
    {
        12 + 2  * I,
        16 + 3  * I,
        12 + 4  * I,
        7  + 5  * I,
        11 + 8  * I,
        6  + 9  * I,
    };

    const fftw_complex aquarius[13] =
    {
        10 + 2  * I,
        11 + 2  * I,
        12 + 3  * I,
        15 + 3  * I,
        4  + 4  * I,
        6  + 5  * I,
        17 + 5  * I,
        12 + 6  * I,
        7  + 7  * I,
        5  + 8  * I,
        12 + 8  * I,
        4  + 9  * I,
        19 + 9  * I,
    };

    const fftw_complex centaurus[11] =
    {
        5  + 1  * I,
        12 + 1  * I,
        1  + 3  * I,
        3  + 3  * I,
        8  + 3  * I,
        10 + 3  * I,
        8  + 5  * I,
        17 + 6  * I,
        10 + 7  * I,
        8  + 9  * I,
        3  + 10 * I,
    };

    const fftw_complex pisces[18] =
    {
        8  + 1  * I,
        9  + 1  * I,
        8  + 2  * I,
        8  + 3  * I,
        6  + 5  * I,
        4  + 7  * I,
        19 + 8  * I,
        6  + 9  * I,
        7  + 9  * I,
        12 + 9  * I,
        13 + 9  * I,
        16 + 9  * I,
        18 + 9  * I,
        21 + 9  * I,
        2  + 10 * I,
        3  + 10 * I,
        18 + 10 * I,
        20 + 10 * I,
    };

    const fftw_complex virgo[9] =
    {
        10 + 2  * I,
        12 + 4  * I,
        2  + 5  * I,
        21 + 5  * I,
        6  + 6  * I,
        14 + 6  * I,
        18 + 6  * I,
        10 + 7  * I,
        6  + 9  * I,
    };

    const fftw_complex apus[4] =
    {
        4  + 4  * I,
        10 + 5  * I,
        7  + 6  * I,
        19 + 7  * I,
    };

    const fftw_complex orion[8] =
    {
        11 + 2  * I,
        7  + 3  * I,
        14 + 4  * I,
        9  + 6  * I,
        11 + 6  * I,
        13 + 6  * I,
        15 + 8  * I,
        8  + 9  * I,
    };

    const fftw_complex caelum[4] =
    {
        15 + 2  * I,
        11 + 4  * I,
        4  + 5  * I,
        2  + 8  * I,
    };

    const fftw_complex aries[4] =
    {
        5  + 3  * I,
        15 + 5  * I,
        17 + 7  * I,
        16 + 8  * I,
    };

    const fftw_complex aquila[7] =
    {
        12 + 1  * I,
        10 + 2  * I,
        7  + 3  * I,
        20 + 3  * I,
        2  + 5  * I,
        11 + 6  * I,
        10 + 9  * I,
    };

    const fftw_complex canis_major[9] =
    {
        13 + 1  * I,
        10 + 2  * I,
        11 + 3  * I,
        14 + 4  * I,
        20 + 5  * I,
        8  + 7  * I,
        6  + 9  * I,
        1  + 10 * I,
        10 + 10 * I,
    };

    const fftw_complex auriga[6] =
    {
        5  + 2  * I,
        14 + 2  * I,
        16 + 4  * I,
        3  + 6  * I,
        15 + 8  * I,
        10 + 10 * I,
    };

    const fftw_complex capricorn[11] =
    {
        3  + 3  * I,
        19 + 3  * I,
        5  + 4  * I,
        8  + 4  * I,
        12 + 4  * I,
        18 + 4  * I,
        17 + 5  * I,
        5  + 6  * I,
        7  + 7  * I,
        14 + 7  * I,
        12 + 8  * I,
    };

    const fftw_complex monoceros[9] =
    {
        19 + 3  * I,
        20 + 4  * I,
        17 + 5  * I,
        21 + 5  * I,
        2  + 6  * I,
        12 + 6  * I,
        6  + 9  * I,
        18 + 9  * I,
        21 + 9  * I,
    };

    const fftw_complex canes_venatici[8] =
    {
        20 + 2  * I,
        11 + 3  * I,
        17 + 5  * I,
        7  + 6  * I,
        5  + 8  * I,
        14 + 8  * I,
        2  + 9  * I,
        10 + 9  * I,
    };

    const fftw_complex andromeda[6] =
    {
        19 + 1  * I,
        6  + 3  * I,
        9  + 4  * I,
        17 + 4  * I,
        12 + 6  * I,
        4  + 9  * I,
    };

    const fftw_complex corona_borealis[7] =
    {
        15 + 3  * I,
        3  + 5  * I,
        20 + 5  * I,
        5  + 7  * I,
        18 + 7  * I,
        10 + 8  * I,
        13 + 8  * I,
    };

    fftw_complex normalizationFactor = 20 + 10 * I;
    const fftw_complex *stars = NULL;
    size_t length = 0;

    switch (name % CONSTELLATION_MAX)
    {
        case OPHIUCHUS:
        {
            stars = ophiuchus;
            length = sizeof(ophiuchus) / sizeof(ophiuchus[0]);
            break;
        }

        case URSA_MINOR:
        {
            stars = ursa_minor;
            length = sizeof(ursa_minor) / sizeof(ursa_minor[0]);
            break;
        }

        case ARA:
        {
            stars = ara;
            length = sizeof(ara) / sizeof(ara[0]);
            break;
        }

        case LUPUS:
        {
            stars = lupus;
            length = sizeof(lupus) / sizeof(lupus[0]);
            break;
        }

        case LEO:
        {
            stars = leo;
            length = sizeof(leo) / sizeof(leo[0]);
            break;
        }

        case URSA_MAJOR:
        {
            stars = ursa_major;
            length = sizeof(ursa_major) / sizeof(ursa_major[0]);
            break;
        }

        case PUPPIS:
        {
            stars = puppis;
            length = sizeof(puppis) / sizeof(puppis[0]);
            break;
        }

        case ANTLIA:
        {
            stars = antlia;
            length = sizeof(antlia) / sizeof(antlia[0]);
            break;
        }

        case SAGITTARIUS:
        {
            stars = sagittarius;
            length = sizeof(sagittarius) / sizeof(sagittarius[0]);
            break;
        }

        case LIBRA:
        {
            stars = libra;
            length = sizeof(libra) / sizeof(libra[0]);
            break;
        }

        case CIRCINUS:
        {
            stars = circinus;
            length = sizeof(circinus) / sizeof(circinus[0]);
            break;
        }

        case CYGNUS:
        {
            stars = cygnus;
            length = sizeof(cygnus) / sizeof(cygnus[0]);
            break;
        }

        case SCORPIO:
        {
            stars = scorpio;
            length = sizeof(scorpio) / sizeof(scorpio[0]);
            break;
        }

        case CAMELOPARDALIS:
        {
            stars = camelopardalis;
            length = sizeof(camelopardalis) / sizeof(camelopardalis[0]);
            break;
        }

        case GEMINI:
        {
            stars = gemini;
            length = sizeof(gemini) / sizeof(gemini[0]);
            break;
        }

        case CANIS_MINOR:
        {
            stars = canis_minor;
            length = sizeof(canis_minor) / sizeof(canis_minor[0]);
            break;
        }

        case CRUX:
        {
            stars = crux;
            length = sizeof(crux) / sizeof(crux[0]);
            break;
        }

        case BOOTES:
        {
            stars = bootes;
            length = sizeof(bootes) / sizeof(bootes[0]);
            break;
        }

        case CARINA:
        {
            stars = carina;
            length = sizeof(carina) / sizeof(carina[0]);
            break;
        }

        case TAURUS:
        {
            stars = taurus;
            length = sizeof(taurus) / sizeof(taurus[0]);
            break;
        }

        case CANCER:
        {
            stars = cancer;
            length = sizeof(cancer) / sizeof(cancer[0]);
            break;
        }

        case CASSIOPEIA:
        {
            stars = cassiopeia;
            length = sizeof(cassiopeia) / sizeof(cassiopeia[0]);
            break;
        }

        case LYRA:
        {
            stars = lyra;
            length = sizeof(lyra) / sizeof(lyra[0]);
            break;
        }

        case AQUARIUS:
        {
            stars = aquarius;
            length = sizeof(aquarius) / sizeof(aquarius[0]);
            break;
        }

        case CENTAURUS:
        {
            stars = centaurus;
            length = sizeof(centaurus) / sizeof(centaurus[0]);
            break;
        }

        case PISCES:
        {
            stars = pisces;
            length = sizeof(pisces) / sizeof(pisces[0]);
            break;
        }

        case VIRGO:
        {
            stars = virgo;
            length = sizeof(virgo) / sizeof(virgo[0]);
            break;
        }

        case APUS:
        {
            stars = apus;
            length = sizeof(apus) / sizeof(apus[0]);
            break;
        }

        case ORION:
        {
            stars = orion;
            length = sizeof(orion) / sizeof(orion[0]);
            break;
        }

        case CAELUM:
        {
            stars = caelum;
            length = sizeof(caelum) / sizeof(caelum[0]);
            break;
        }

        case ARIES:
        {
            stars = aries;
            length = sizeof(aries) / sizeof(aries[0]);
            break;
        }

        case AQUILA:
        {
            stars = aquila;
            length = sizeof(aquila) / sizeof(aquila[0]);
            break;
        }

        case CANIS_MAJOR:
        {
            stars = canis_major;
            length = sizeof(canis_major) / sizeof(canis_major[0]);
            break;
        }

        case AURIGA:
        {
            stars = auriga;
            length = sizeof(auriga) / sizeof(auriga[0]);
            break;
        }

        case CAPRICORN:
        {
            stars = capricorn;
            length = sizeof(capricorn) / sizeof(capricorn[0]);
            break;
        }

        case MONOCEROS:
        {
            stars = monoceros;
            length = sizeof(monoceros) / sizeof(monoceros[0]);
            break;
        }

        case CANES_VENATICI:
        {
            stars = canes_venatici;
            length = sizeof(canes_venatici) / sizeof(canes_venatici[0]);
            break;
        }

        case ANDROMEDA:
        {
            stars = andromeda;
            length = sizeof(andromeda) / sizeof(andromeda[0]);
            break;
        }

        case CORONA_BOREALIS:
        {
            stars = corona_borealis;
            length = sizeof(corona_borealis) / sizeof(corona_borealis[0]);
            break;
        }

        default:
        {
            fprintf(stderr, "unknown constellation type: %d/n", name);
            return;
        }
    }

    // calculate the average position of the constellation to balance the
    // weight of the phase components so hopefully you don't get clipping of the waveform
    // also calculate the normalization factor so the RMS power is 1

    constellation->points = malloc(sizeof(fftw_complex) * length);
    constellation->length = length;

    if(length < 2)
    {
        constellation->points[0] = 0;
        return;
    }

    fftw_complex average = 0;
    for(int i = 0; i < length; i++)
        average += stars[i];
    average /= length;

    // choose a random point in the selected constellation
    int index = rand() % length;
    double RMS = 0;
    for(int index = 0; index < length; index++)
    {
        constellation->points[index] = (creal(stars[index] - average) / creal(normalizationFactor) + cimag(stars[index] - average) / cimag(normalizationFactor) * -I);
        RMS += cabs(constellation->points[index]);
    }
    RMS = RMS == 0 ? 1 : RMS / length * M_SQRT1_2;

    for(int i = 0; i < length; i++)
        constellation->points[i] /= RMS;
    return;
    //return 0;
    //return stars[index] / cimag(normalizationFactor) * 2 - 1 - I;
}


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

    /*
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
    */

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

double raisedCosQAM(int n, int sampleRate)
{
    int carrierPeriod = 32;
    int k = 1; // cycles per period
    //double carrierFrequency = 5000;
    double carrierFrequency = (double)sampleRate / carrierPeriod;
    //int carrierFrequency = sampleRate / carrierPeriod;

    //int symbolPeriod = 64; // audio samples per symbol
    //double symbolPeriod = sampleRate / carrierFrequency * k; // audio samples per symbol
    // array to store time series of filter data
    static double *filter = NULL;
    // array to store timeseries of IQ samples
    static iqsample_t *IQdata = NULL;

    if (n < 0)
    {
        if (filter != NULL)
        {
            free(filter);
            filter = NULL;
        }

        if (IQdata != NULL)
        {
            free(IQdata);
            IQdata = NULL;
        }

        return NAN;
    }

    // This is a sorta bug. the fact that everything is based on these integers is an issue. I think it causes discretization of carrier
    // frequency and symbol periods, which means the output frequency is not what you put in. for example, 5000 Hz input turns out to be
    // about 5500 Hz actual
    int symbolPeriod = (int)(sampleRate / carrierFrequency * k); // audio samples per symbol
    int filterSides = 10;    // number of symbols to either side of current symbol to filter with raised cos
    int filterLengthSymbols = 2 * filterSides + 1;    // length of raised cos filter in IQ symbols, ie, how many IQ samples we need to generate the current symbol
    int filterLength = filterLengthSymbols * symbolPeriod;  // length in audio samples


    // phase offset that cycles through sequentially all phase offsets
    //int phaseOffset = n *  2 / 1200 % symbolPeriod;
    int phaseOffset = 0;
    
    /*
    static int phaseOffset = -1;
    if (phaseOffset == -1)
        phaseOffset = rand() % symbolPeriod;
    */
    /*
    static int phaseOffset = -1;
    if (phaseOffset == -1 || (n % (carrierPeriod * 30)) == 0)
        phaseOffset = rand() % symbolPeriod;
    */

    int originalN = n;  // for plotting
    n -= phaseOffset;
    n = n < 0 ? 0 : n - phaseOffset, 0;

    // concept:
    //  generate the raised cos filter data once
    //  generate IQ samples ahead of time, just in time
    //  generate audio samples based on those two pieces of info

    static int initialized = 0;
    if(!initialized)
    {
        // initialize raised cos filter data
        filter = malloc((unsigned int)(sizeof(double) * filterLength)); // this never gets released. so, might wanna fix that TODO
        for(int i = 0; i < filterLength; i++)
        {
            int filter_symbolIndex = i - filterLength / 2;    // should go -filterLength/2 -> 0 -> filterLength/2
            // raised cos filter math
            double b = 0.42;    // filter parameter beta
            double filterValue = sin(M_PI * filter_symbolIndex / symbolPeriod) / (M_PI * filter_symbolIndex / symbolPeriod) * (cos(M_PI * b * filter_symbolIndex / symbolPeriod)) / (1 - pow(2 * b * filter_symbolIndex / symbolPeriod, 2));
            if(!isfinite(filterValue))   // in case it's undefined, ie divide by zero case
                filterValue = (double)symbolPeriod / 2 / b;
            if(filter_symbolIndex == 0)
                filterValue = 1;    // the math gives a divide by zero at index 0

            filter[i] = filterValue;
        }

        // generate enough IQ samples for first audio sample
        IQdata = malloc((unsigned int)(sizeof(iqsample_t) * filterLengthSymbols));  // TODO this is never freed. These should prob be passed in as a paramter and freed somewhere in the larger scope
        for(int i = 0; i < filterLengthSymbols; i++)
        {
            int IQIndex = i - filterLengthSymbols / 2;    // shoud go from -filterlengthsymbols / 2 -> 0 -> filterlengthsymbols / 2
            iqsample_t sample = {0, 0};
            // fill the circular buffer in the right order so that index 0 turns out to be the first IQ sample
            if(IQIndex < 0)
            {
                IQIndex = filterLengthSymbols + IQIndex;    // make positive
                IQdata[IQIndex % filterLengthSymbols] = sample; // the negative time samples are 0
            } else {
                //IQdata[IQIndex % filterLengthSymbols] = impulse(IQIndex, 20);
                //IQdata[IQIndex % filterLengthSymbols] = alternateI(IQIndex);
                //IQdata[IQIndex % filterLengthSymbols] = alternateQ(IQIndex);
                //IQdata[IQIndex % filterLengthSymbols] = randomQAM(2);
                //IQdata[IQIndex % filterLengthSymbols] = sequentialIQ(IQIndex, 4);
                IQdata[IQIndex % filterLengthSymbols] = randomQAM_withPreamble(IQIndex, 2);
            }
        }

        // ensure initialize only runs once
        initialized = 1;
    }

    // current symbol number
    int symbolIndex = n / symbolPeriod;
    int sampleIndex = n % symbolPeriod; // index of each sample in a symbol, where as n is increasing for the whole signal
    //printf("phaseOffset: %i\tn: %i\tsymbolIndex: %i\tsampleIndex: %i\n", phaseOffset, n, symbolIndex, sampleIndex);
    //int sampleIndex = fmod(n, symbolPeriod); // index of each sample in a symbol, where as n is increasing for the whole signal
    // circular buffer index
    int IQsampleIndex = symbolIndex % filterLengthSymbols;    // IQdata[IQsampleIndex] is current IQ sample, indexes above are future IQdata and below past IQ data both wrapping around until the are filterLengthSymbols / 2 away from current sample index

    // generate the next future IQ sample when entering a new symbol period
    //IQdata[(IQsampleIndex + filterLengthSymbols / 2) % filterLengthSymbols] = alternateI(symbolIndex + filterLengthSymbols / 2);
    if(sampleIndex == 0)
    {
        int IQindex = (IQsampleIndex + filterLengthSymbols / 2) % filterLengthSymbols;
        //IQdata[IQindex] = impulse(symbolIndex + filterLengthSymbols / 2, 20);
        //IQdata[IQindex] = alternateI(symbolIndex + filterLengthSymbols / 2);
        //IQdata[IQindex] = alternateQ(symbolIndex + filterLengthSymbols / 2);
        //IQdata[IQindex] = randomQAM(2);
        //IQdata[IQindex] = sequentialIQ(symbolIndex + filterLengthSymbols / 2, 4);
        IQdata[IQindex] = randomQAM_withPreamble(symbolIndex + filterLengthSymbols / 2, 2);
        IQindex = IQsampleIndex;
    #if DEBUG_LEVEL > 0
        fprintf(plotStdIn, "%i %i %f %i %f\n", originalN, 6, IQdata[IQsampleIndex].InPhase, 7, IQdata[IQsampleIndex].Quadrature);
        fprintf(plotStdIn, "%i %i %i\n", originalN, 8, IQindex);
    #endif
    }

    // add up raised cos contributions from all samples in the IQdata array
    iqsample_t filteredIQsample = {0, 0};
    for(int i = 0; i < filterLengthSymbols; i++)
    //int i = symbolIndex % filterLengthSymbols;
    //i = filterLengthSymbols - i - 1;
    {
        int filterIndex = (filterLengthSymbols - i - 1) * symbolPeriod + sampleIndex;    // pick a filter index
        int IQIndex = (IQsampleIndex + i - filterLengthSymbols / 2) % filterLengthSymbols;
        if(IQIndex < 0)
            IQIndex = filterLengthSymbols + IQIndex;    // make sure index is positive and wraps backwards
        filteredIQsample.InPhase += filter[filterIndex] * IQdata[IQIndex].InPhase;
        filteredIQsample.Quadrature += filter[filterIndex] * IQdata[IQIndex].Quadrature;
        //filteredIQsample.I += filter[filterIndex];
        //filteredIQsample.I = IQdata[IQIndex].I;
        //filteredIQsample.I += IQdata[IQIndex].I / filterLengthSymbols;
    #if DEBUG_LEVEL > 0
        //fprintf(plotStdIn, "%i %i %f\n", originalN + i, 8 + symbolIndex, IQdata[IQIndex].I);
    #endif
    }
    // normalization for raised cos filter, prob not correct
    //filteredIQsample.I /= 2;
    //filteredIQsample.Q /= 2;
    filteredIQsample.InPhase *= M_SQRT1_2;    // normalization factor for QAM constellation with max I and Q of 1, converts to max magnitude of 1
    filteredIQsample.Quadrature *= M_SQRT1_2;

    //return filter[symbolIndex%filterLengthSymbols*symbolPeriod+sampleIndex];
    //return IQdata[(IQsampleIndex + n%filterLengthSymbols - filterLengthSymbols / 2) %filterLengthSymbols].I;
    //return IQdata[IQsampleIndex].I;
    //return filteredIQsample.I;
    //return (double)i / filterLengthSymbols;
    //return filter[n%filterLength];

    //if(n % symbolPeriod < filterLengthSymbols)
        //return IQdata[(IQsampleIndex + sampleIndex) % filterLengthSymbols].I;
    //return -0.5;

    //filteredIQsample.I = 1;
    //filteredIQsample.Q = 0;

    double audioSample =
    (
        (filteredIQsample.InPhase) * cos(2.0 * M_PI * sampleIndex * k / symbolPeriod) +
        (filteredIQsample.Quadrature) * sin(2.0 * M_PI * sampleIndex * k / symbolPeriod)
    ) / 2.0;    // normalization factor for sum of sin+cos of amp 1 not exceeding 1
#if DEBUG_LEVEL > 0
    fprintf(plotStdIn, "%i %i %i %i %i %i %f %i %i\n", originalN, 0, symbolIndex, 1, sampleIndex, 2, audioSample, 3, phaseOffset);
    fprintf(plotStdIn, "%i %i %f %i %f\n", originalN, 4, filteredIQsample.InPhase, 5, filteredIQsample.Quadrature);
#endif

   return audioSample;

}

// runs the state machine for an OFDM based communication transmitter
buffered_data_return_t OFDM(int long n, sample_double_t *outputSample, OFDM_state_t *OFDMstate)
{
    if(OFDMstate->initialized == 0)
    {
        // initialize the state variables
        OFDMstate->state.frame = IDLE;
        OFDMstate->state.frameStart = n;

        //OFDMstate->sampleRate = 44100;
        //OFDMstate->guardPeriod = 0.13 * OFDMstate->sampleRate;  // measured impulse response of room lasts like a bit over a tenth of a second would be nice to have this adjusted on the fly
        OFDMstate->guardPeriod = (1<<12);  // 2^12=1024 closest power of 2 to the impulse response length, slightly shorter
        //OFDMstate->guardPeriod = 128;   // faster testing
        //OFDMstate->guardPeriod = 0.001 * OFDMstate->sampleRate;  // measured impulse response of room lasts like a bit over a tenth of a second would be nice to have this adjusted on the fly
        OFDMstate->ofdmPeriod = OFDMstate->guardPeriod * 4;   // dunno what the best OFDM period is compared to the guard period. I assume longer is better for channel efficiency, but maybe it's worse for noise? don't know
        //OFDMstate->ofdmPeriod = OFDMstate->guardPeriod * 8;   // dunno what the best OFDM period is compared to the guard period. I assume longer is better for channel efficiency, but maybe it's worse for noise? don't know
        OFDMstate->symbolPeriod = OFDMstate->guardPeriod + OFDMstate->ofdmPeriod;
        OFDMstate->channels = OFDMstate->ofdmPeriod / 2 + 1;  // half due to using real symbols (ie, not modulating to higher frequency carrier wave but staying in baseband) ie niquist
        
        // initialize channel simulation filter
        OFDMstate->simulateNoise = 1;
        OFDMstate->simulateChannel = 1;


        if(OFDMstate->simulateChannel)
        {
            // initialize the overlap and save buffer
            OFDMstate->channelSimulationBuffer.M = 5912;
            OFDMstate->channelSimulationBuffer.L = OFDMstate->channelSimulationBuffer.M;
            OFDMstate->channelSimulationBuffer.length = OFDMstate->channelSimulationBuffer.L;
            OFDMstate->channelSimulationBuffer.N = OFDMstate->channelSimulationBuffer.L + OFDMstate->channelSimulationBuffer.M;
            OFDMstate->channelSimulationBuffer.insertionIndex = 0;
            //OFDMstate->channelSimulationBuffer.buffer = calloc(OFDMstate->channelSimulationBuffer.length, sizeof(double));
            OFDMstate->channelSimulationBuffer.M_buffer = fftw_malloc(sizeof(double) * OFDMstate->channelSimulationBuffer.N);    // use fftw malloc since fftw will use this buffer
            if(OFDMstate->channelSimulationBuffer.M_buffer == NULL)
                fprintf(stderr, "cahnnelSimulationBuffer failed to allocate: %s\n", strerror(errno));
            OFDMstate->channelSimulationBuffer.L_buffer = OFDMstate->channelSimulationBuffer.M_buffer + OFDMstate->channelSimulationBuffer.M;   // offset the L
            OFDMstate->channelSimulationBuffer.buffer = OFDMstate->channelSimulationBuffer.L_buffer;
            OFDMstate->channelSimulationBuffer.M_next_buffer = OFDMstate->channelSimulationBuffer.M_buffer + OFDMstate->channelSimulationBuffer.L;  // M_next_buffer is the last M numbers of L
        }

        // initialize fftw arrays
        OFDMstate->OFDMsymbol.frequencyDomain = (fftw_complex*)malloc(sizeof(fftw_complex) * OFDMstate->channels);
        if(OFDMstate->OFDMsymbol.frequencyDomain == NULL)
            fprintf(stderr, "couldn't allocate fftw input array\n");

        OFDMstate->OFDMsymbol.timeDomain = (double*)malloc(sizeof(double) * OFDMstate->ofdmPeriod);
        if(OFDMstate->OFDMsymbol.timeDomain == NULL)
            fprintf(stderr, "couldn't allocate fftw output array\n");

        // initialize fftw plan, using the measure option to calculate the fastest plan, could have a few seconds startup time
        OFDMstate->fftwPlan = fftw_plan_dft_c2r_1d(OFDMstate->ofdmPeriod, OFDMstate->OFDMsymbol.frequencyDomain, OFDMstate->OFDMsymbol.timeDomain, FFTW_MEASURE);

        // initialize constellations
        OFDMstate->constellationsLength = 7;
        OFDMstate->constellations = malloc(sizeof(constellation_complex_t) * OFDMstate->constellationsLength);
        for(int i = 0; i < OFDMstate->constellationsLength; i++)
        {
            fftw_squareQAM(pow(2, i), &OFDMstate->constellations[i]);
        }

        OFDMstate->initialized = 1;
    }

    double output = 0;

    switch(OFDMstate->state.frame)
    {
        case IDLE:

            // send nothing
            output = 0;
            
            // wait for enough data in the pipe to generate a frame, or wait for enough time to pass to send a frame anyway.

            // transmit a burst of noise at high power to help synchronizer detect future ISI?
            //if(n - OFDMstate->state.frameStart >= OFDMstate->symbolPeriod - OFDMstate->guardPeriod - 1)
            if(0)
            {
                double noiseAmplitude = 0.1 * OFDMstate->ofdmPeriod;  // there is a normalization factor
                output += (double)rand() / ((double)RAND_MAX / (noiseAmplitude)) - (noiseAmplitude / 2);
            }

            // check for exit from IDLE frame
            if(n - OFDMstate->state.frameStart >= OFDMstate->symbolPeriod * 1 - 1) // example of state change based on timing
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

                    output = 0.70; // testing only

                    switch(OFDMstate->state.symbol)
                    {
                        case GUARD_PERIOD:

                            if(n - OFDMstate->state.symbolStart == 0)
                            {
                                // set new symbol for testing
                                // generate a symetric symbol, with even frequency components only for syncronization
                                for(int k = 0; k < OFDMstate->channels; k++)
                                {
                                    if(OFDMstate->state.symbolIndex == 0)
                                    {
                                        if(k % 2 == 0)
                                        //if(k > 16 && k < 100 && k % 2 == 0)
                                        //if(k > 100 && k < 200 && k % 2 == 0)
                                        //if(k > 446 && k < 446+200 && k % 2 == 0)
                                        //if(k == OFDMstate->channels / 16)
                                        //if(k == 0)
                                        //if(0)
                                        {
                                            OFDMstate->OFDMsymbol.frequencyDomain[k] =
                                                rand() % 2 * 2 - 1 +
                                                I*(rand() % 2 * 2 - 1);
                                                //I;    // testing
                                        } else {
                                            //OFDMstate->currentOFDMSymbol[k] = 0;
                                            OFDMstate->OFDMsymbol.frequencyDomain[k] = 0;
                                        }

                                        //OFDMstate->OFDMsymbol.frequencyDomain[k] *= (double)OFDMstate->channels / 200 / 30;

                                    } else if(OFDMstate->state.symbolIndex == 1)
                                    {
                                        if(1)
                                        //if(k > 100 && k < 200 && k % 2 == 0)
                                        //if(k > 446 && k < 446+10)
                                        {
                                            OFDMstate->OFDMsymbol.frequencyDomain[k] =
                                                rand() % 2 * 2 - 1 +
                                                I*(rand() % 2 * 2 - 1);
                                                //0;
                                        } else {
                                            OFDMstate->OFDMsymbol.frequencyDomain[k] =
                                                0;
                                        }
                                        //OFDMstate->currentOFDMSymbol[k] = 0;    // testing
                                        //OFDMstate->OFDMsymbol.frequencyDomain[k] *= (double)OFDMstate->channels / 10 / 30;
                                    }
                                }
                                // now do the transform
                                fftw_execute(OFDMstate->fftwPlan);
                                // time domain samples are now in the OFDMstate->OFDMsymbol.timeDomain array
                            }

                            // output samples for the guard period
                            //output = OFDMsymbolBaseband(OFDMstate->ofdmPeriod - OFDMstate->guardPeriod + (n - OFDMstate->state.symbolStart),
                                                        //OFDMstate);
                            output = OFDMstate->OFDMsymbol.timeDomain[OFDMstate->ofdmPeriod - OFDMstate->guardPeriod + (n - OFDMstate->state.symbolStart)];

                            if(OFDMstate->state.symbolIndex == 0)
                            {
                                output *= M_SQRT2;  // scale factor for even only channels
                                //output = 0; // trying no guard period for the first symbol of preamble, to help find the precise sync point
                            }
                            //output = output * 0.5;

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
                            output = OFDMstate->OFDMsymbol.timeDomain[n - OFDMstate->state.symbolStart];

                            if(OFDMstate->state.symbolIndex == 0)
                                output *= M_SQRT2;
                            //output = output;

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
                    if(n - OFDMstate->state.fieldStart >= OFDMstate->symbolPeriod * 2 - 1)   // for 2 symbol preamble
                    {
                        OFDMstate->state.field = DATA;
                        OFDMstate->state.fieldStart = n + 1;
                        OFDMstate->state.symbolIndex = 0;
                    }
                    break;
                case DATA:

                    output = 1; // testing only

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
                                static int center;
                                static int width;
                                if(rand() % 5 == 0 || OFDMstate->state.symbolIndex == 0)
                                {
                                    center = rand() % OFDMstate->channels;
                                    width = rand() % OFDMstate->channels;
                                }
                                for(int k = 0; k < OFDMstate->channels; k++)
                                {
                                    if(1)   // transmit on all subchannels
                                    //if(k > 250 && k < 6000)
                                    //if(k > 250 && k < 350)
                                    //if(k > 7000 && k < 8000) // using a select number of channels to simplify the signal for testing
                                    //int startChunk = (OFDMstate->state.symbolIndex / 1 * 100) % OFDMstate->channels / 4;
                                    //if(k % (OFDMstate->channels / 4) > startChunk && k % (OFDMstate->channels / 4) < startChunk + 10)
                                    // choose random blocks
                                    //if(k < center + width / 2 && k > center - width / 2)
                                    {
                                        
                                        // picking a sequential constellation, and a random point in that constellation
                                        constellation_complex_t constellation = OFDMstate->constellations[k%OFDMstate->constellationsLength];
                                        OFDMstate->OFDMsymbol.frequencyDomain[k] = constellation.points[rand() % constellation.length];
                                        //rescale if I'm using only a few channels for higher power 
                                        //OFDMstate->OFDMsymbol.frequencyDomain[k] *= (double)OFDMstate->channels / 10 / 30;
                                        //OFDMstate->OFDMsymbol.frequencyDomain[k] /= 3;

                                    } else {
                                        OFDMstate->OFDMsymbol.frequencyDomain[k] = 0;
                                    }
                                }
                                // now do the transform
                                fftw_execute(OFDMstate->fftwPlan);
                                // time domain samples are now in the OFDMstate->OFDMsymbol.timeDomain array
                            }
                            //output = OFDMsymbolBaseband(OFDMstate->ofdmPeriod - OFDMstate->guardPeriod + (n - OFDMstate->state.symbolStart),
                                                        //OFDMstate);
                            output = OFDMstate->OFDMsymbol.timeDomain[OFDMstate->ofdmPeriod - OFDMstate->guardPeriod + (n - OFDMstate->state.symbolStart)];

                            //output = output * 0.5;

                            // check for exit
                            if(n - OFDMstate->state.symbolStart >= OFDMstate->guardPeriod - 1)
                            {
                                OFDMstate->state.symbol = OFDM_PERIOD;
                                OFDMstate->state.symbolStart = n + 1;
                            }
                            break;
                        case OFDM_PERIOD:

                            //output = OFDMsymbolBaseband(n - OFDMstate->state.symbolStart,
                                                        //OFDMstate);
                            output = OFDMstate->OFDMsymbol.timeDomain[n - OFDMstate->state.symbolStart];
                            //output = output;

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

    /*
    // impulse train for testing the channel simulator
    output = 0;
    if(n % 10000 == 0)
        output = 1;
    */

    // channel simulation
    sample_double_t equalizedSample;
    // decide whether to simulate the channel response or not
    if(OFDMstate->simulateChannel)
    {
        // channel simulation filter
        //OFDMstate->channelSimulationBuffer.insertionIndex = n % OFDMstate->channelSimulationBuffer.length;
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
        equalizedSample.sample += (double)rand() / ((double)RAND_MAX / (noiseAmplitude)) - (noiseAmplitude / 2);
    }

    *outputSample = equalizedSample;

    return RETURNED_SAMPLE;
}

// really this is generating a single OFDM channel without guard periods,
// missing some critical infrustructure
static double singleChannelODFM_noguard(int n, int sampleRate)
{
    int symbolPeriod = 256;
    int guardPeriod = 4./1000 * sampleRate;     // I've found that the echos in a room last for about 3ms, durring that period, the symbol is phase offset and otherwise changed due to the last symbol and the transition between symbols
    //int guardPeriod = 0;    // have to disable for QAM, ie, set to 0, instead use a raised cosine filter for ISI combat
    int totalPeriod = symbolPeriod + guardPeriod;
    int k = 16;      // this is effectively the OFDM channel number, how many cycles per sample period

    // generating offsets in time to test the frame time syncronizer in qamDecoder
    // phase offset that cycles through sequentially all phase offsets
    //int phaseOffset = n *  2 / 2000 % symbolPeriod;
    // no phase offset
    //int phaseOffset = 0;
    // random phase offset
    static int phaseOffset = -1;
    if (phaseOffset == -1)
        phaseOffset = rand() % symbolPeriod;

    // apply phase offset
    n += phaseOffset;

    iqsample_t sample;  // The IQ sample currently being worked on, temporal resolution to audio sample
    double audioSample; // the final audio sample
    int symbolIndex = n / totalPeriod;  // symbol number indexed from 0, temporal resolution to IQ sample

    // random IQ in constelation defined by power
    int power = 2;  // log base2 of number of symbols. number of symbols should also be a perfect square
    int symbols = (int)pow(2, power);
    int square = (int)sqrt(symbols);

    // get IQ samples from the IQ sampling function
    sample = randomQAM_withPreamble(symbolIndex, square);

    // amplitude adjustment
    //double totalAmplitude = 0.01;
    double totalAmplitude = 1;
    //double totalAmplitude = 0.5;

    // random noise injection to IQ value
    double randomness = 0.0;
    double randI = ((double)rand() / RAND_MAX * 2 - 1) * randomness;
    double randQ = ((double)rand() / RAND_MAX * 2 - 1) * randomness;

    // implementing guard period
    // turns out this fucks up the gardner algorithm. ONLY for OFDM, not for QAM. Will use this later probs. disable by setting guard period to 0
    // time index of each symbol to resolution of audio sample
    int symbolStep = n % totalPeriod - guardPeriod;   // should be guardPeriod -> symbolPeriod -> 0 -> symbolPeriod as n increases from 0 -> totalPeriod
    audioSample =
    (
        (sample.InPhase + randI) * cos(2.0 * M_PI * symbolStep * k / symbolPeriod) +
        (sample.Quadrature + randQ) * sin(2.0 * M_PI * symbolStep * k / symbolPeriod)
    ) / 2.0 * sqrt(2.0) * totalAmplitude;

   return audioSample;
}

// this is the point where samples are generated
static double WARN_UNUSED calculateSample(int n, int sampleRate)
{

    // I think I need to keep polling for a sample until SAMPLE RETURNED

    //double amplitudeScaler = 0.1;
    double amplitudeScaler = 1;
    /*
    if(n == 2500)
        return 1;
    return 0;
    */


    static OFDM_state_t OFDMstate = {0};
    OFDMstate.sampleRate = sampleRate;
    sample_double_t returnedSample;

    static buffered_data_return_t returnValue = AWAITING_SAMPLES;
    static int offset = 0;
    if(returnValue == AWAITING_SAMPLES)
    {
        // run until first sample is returned
        for(int i = n; returnValue == AWAITING_SAMPLES; i++)
        {
            // call the function as many times as needed until it returns a sample
            returnValue = OFDM(i, &returnedSample, &OFDMstate);
            offset = i;
        }
        return returnedSample.sample;
    }
    // from then on just run once per sample
    OFDM(n + offset, &returnedSample, &OFDMstate);
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


    int fileDescriptor = -1;

    // audio sample rate
    // Supported sample rates from alsa-info.sh
    //     rates [0x560]: 44100 48000 96000 192000
    int sampleRate = 44100;
    // total number of samples to generate
    long length = sampleRate * 120;
    // the number of the current sample
    long n = 0;

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
    aplayStdIn = popen(aplayCommandString, "w");

    if (aplayStdIn == NULL)
    {
        fprintf(stderr, "Failed to open aplay: %s\n", strerror(errno));
        retval = 3;
        goto exit;
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
    srand(time(NULL));
    while(n < length)
    {
        // calculate a chunk of samples until the buffer is full or max is reached. one sample at a time, 4 bytes at a time
        for(bufferReadyBytes = 0; (bufferReadyBytes < bufferLength) && (n < length); bufferReadyBytes += 4, n++)
        {
            // the sample value used in calculations, to be normalized
            double sampleValue;
            // sample value after put into signed integer range, then split into bytes for file writing and audio output
            int32_to_bytes_t normalizedSampleValue;
            // holds each individual byte as it's written out Little Endian style
            char byte;

            // get the double sample value, should be between -1 and 1
            sampleValue = calculateSample(n, sampleRate);
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
