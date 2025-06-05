#include <stdlib.h>
#include <stdint.h>
#include <stdio.h>
#include <unistd.h>
#include <fcntl.h>
#include <math.h>
#include <string.h>
#include <time.h>
#include <errno.h>

#include <complex.h>
#include <fftw3.h>

#include "utilities.h"


void initializeOFDMstate(OFDM_state_t *OFDMstate)
{
    // initialize the state

    OFDMstate->guardPeriod = (1<<12);  // 2^12=1024 closest power of 2 to the impulse response length, slightly shorter
                                       //OFDMstate->guardPeriod = 128;
    OFDMstate->ofdmPeriod = OFDMstate->guardPeriod * 4;   // dunno what the best OFDM period is compared to the guard period. I assume longer is better for channel efficiency, but maybe it's worse for noise? don't know
                                                          //OFDMstate->ofdmPeriod = OFDMstate->guardPeriod * 8;   // dunno what the best OFDM period is compared to the guard period. I assume longer is better for channel efficiency, but maybe it's worse for noise? don't know
    OFDMstate->symbolPeriod = OFDMstate->guardPeriod + OFDMstate->ofdmPeriod;
    OFDMstate->channels = OFDMstate->ofdmPeriod / 2 + 1;  // half due to using real symbols (ie, not modulating to higher frequency carrier wave but staying in baseband)
    OFDMstate->pilotSymbolsPitch = 25;
    
    // initialize fftw arrays for both recieve and transmit
    initializeCircularBuffer_fftw_complex(&OFDMstate->OFDMsymbol.frequencyDomain, OFDMstate->channels, 0);
    initializeCircularBuffer_double(&OFDMstate->OFDMsymbol.timeDomain, OFDMstate->ofdmPeriod, 0);

    // initialize constellations and their huffman trees
    OFDMstate->constellationsLength = 16;
    OFDMstate->constellations = malloc(sizeof(constellation_complex_t) * OFDMstate->constellationsLength);
    if(OFDMstate->constellations == NULL)
        fprintf(stderr, "constellations failed to allocate: %s\n", strerror(errno));
    for(int i = 0; i < OFDMstate->constellationsLength; i++)
    {
        //fftw_squareQAM(i, &OFDMstate->constellations[i]);
        if(i==0)    // first constellation is used for the preamble
        {
            fftw_ASK(2, &OFDMstate->constellations[i]);
        } else {

            // for data channels
            fftw_squareQAM(i, &OFDMstate->constellations[i]);
            /*
               if(i / 16 == 0)
               fftw_hexQAM(i / 16, &OFDMstate->constellations[i % 16]);
               else if(i / 16 == 1)
               fftw_squareQAM(i / 16, &OFDMstate->constellations[i % 16]);
               else if(i / 16 == 1)
               fftw_PSK(i / 16, &OFDMstate->constellations[i % 16]);
               else if(i / 16 == 1)
               fftw_ASK(i / 16, &OFDMstate->constellations[i % 16]);
               */
        }

        // generate the huffman tree for the number of points in the constellation
        generateHuffmanTree(&OFDMstate->constellations[i]);

    }

    // initialize random number generators
    srand(1);
    srand48_r(rand(), &OFDMstate->preamblePilotsPRNG);
    srand48_r(rand(), &OFDMstate->pilotsPRNG);
    srand48_r(rand(), &OFDMstate->predefinedDataPRNG);
    srand48_r(rand(), &OFDMstate->channelNoisePRNG);

}

void initializeCircularBuffer_fftw_complex(circular_buffer_complex_t *buf, int length, int sampleRate)
{
    buf->length = length;
    buf->insertionIndex = 0;
    buf->phase = 0;
    buf->n = 0;
    buf->sampleRate = sampleRate;
    buf->buffer = fftw_malloc(sizeof(fftw_complex) * buf->length);
    if(buf->buffer == NULL)
        fprintf(stderr, "failed to allocate a complex circlular buffer. %s\n", strerror(errno));
}
void initializeCircularBuffer_complex(circular_buffer_complex_t *buf, int length, int sampleRate)
{
    initializeCircularBuffer_fftw_complex(buf, length, sampleRate);
}

void initializeCircularBuffer_double(circular_buffer_double_t *buf, int length, int sampleRate)
{
    buf->length = length;
    buf->insertionIndex = 0;
    buf->buffer = fftw_malloc(sizeof(double) * buf->length);
    if(buf->buffer == NULL)
        fprintf(stderr, "failed to allocate a double circlular buffer. %s\n", strerror(errno));
}

void initializeOverlapAndSaveBuffer(overlap_save_buffer_double_t *buf, int length)
{
    buf->length = length;
    buf->insertionIndex = 0;
    buf->buffer = calloc(buf->length, sizeof(double));
    if(buf->buffer == NULL)
        fprintf(stderr, "cahnnelSimulationBuffer failed to allocate: %s\n", strerror(errno));
}



void generateHuffmanTree(constellation_complex_t *constellation)
{
    // generate a huffman tree based on the number of points in the given constellation object
    // and add a pointer to the tree in the constellation object
    // assuming all points must be as equaprobable as possible (won't be possible except powers of 2 number of points)


    // linked list
    typedef struct listElem
    {
        struct listElem *next;
        huffman_tree_t *huffmanTree;
    } listElem;

    // VLA
    listElem roots[constellation->length];
    huffman_tree_t *huffmanList = malloc(sizeof(huffman_tree_t) * constellation->length);
    for(int i = 0; i < constellation->length; i++)
    {
        // if it's not the last in the list, reference the next in the linked list
        if(i != constellation->length - 1)
            roots[i].next = &roots[i+1];
        else
            roots[i].next = 0;

        roots[i].huffmanTree = &huffmanList[i];
        roots[i].huffmanTree->isLeaf = 1;
        roots[i].huffmanTree->constellationIndex = i;
        // clear the children pointers
        roots[i].huffmanTree->child[0] = 0;
        roots[i].huffmanTree->child[1] = 0;
    }

    int rootsNumber = constellation->length;
    listElem *root = &roots[0];
    listElem *tail = &roots[constellation->length - 1];
    while(rootsNumber > 1)
    {
        // each loop, combine the huffman nodes from the top two items of the list under a new huffman node
        // then add the new huffman node to the end of the linked list.

        listElem *oldRoot = root;   // first item
        listElem *newRoot = root->next->next;// third item down the list

        // the two huffman nodes of interest
        huffman_tree_t *huffman1 = oldRoot->huffmanTree;
        huffman_tree_t *huffman2 = oldRoot->next->huffmanTree;

        // we'll reuse the oldRoot linked list element for the new huffman node
        oldRoot->huffmanTree = malloc(sizeof(huffman_tree_t));  // create a new huffman node, adding the two found nodes as children
        oldRoot->huffmanTree->isLeaf = 0;
        oldRoot->huffmanTree->child[0] = huffman1;
        oldRoot->huffmanTree->child[1] = huffman2;

        // clear the constellation index field
        oldRoot->huffmanTree->constellationIndex = -1;

        // add the oldRoot to the end of the linked list
        tail->next = oldRoot;
        tail = oldRoot;
        tail->next = 0;
        // we'll use the new root to hold the newly created huffman tree node

        // see if this is the last run
        if(rootsNumber > 2)
            root = newRoot; // set the new root
        else
            root = oldRoot; // last root is the list item containing the newly created huffman node

        // clear out the data for the unused linked list item
        // we loose one list item each loop
        rootsNumber--;  // subtract one from the length of the linked list
    }

    // at the end of the loop, we should have a huffman tree constructed and pointed to by the root liste item
    constellation->huffmanTree = root->huffmanTree;
}

complex double traverseHuffmanTree(OFDM_state_t *OFDMstate, constellation_complex_t *constellation)
{
    // advance the bit pointer and consume bytes from dataInput file descriptor until a huffman leaf is reached
    huffman_tree_t *node = constellation->huffmanTree;
    char bitmask = 1;
    OFDMstate->inputByte = 85;
    while(!node->isLeaf)
    {
        bitmask = 1 << OFDMstate->bitOffset;

        node = node->child[(bitmask & OFDMstate->inputByte) > 0];
        if(OFDMstate->bitOffset >= 7)
        {
            // get a new byte
            OFDMstate->bitOffset = 0;
        } else {
            OFDMstate->bitOffset++;
        }
        //OFDMstate->bitOffset = OFDMstate->bitOffset == 7 ? OFDMstate->bitOffset = 0 : OFDMstate->bitOffset+1;
    }

    return constellation->points[node->constellationIndex];
}

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
