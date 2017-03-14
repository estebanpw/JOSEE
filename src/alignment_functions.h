#include <cstdio>
#include <inttypes.h>
#include <cstdlib>
#include <float.h>
#include <cstdint>
#include "structs.h"
/*
    Returns the score for a nucl
*/

int valOfNucl(char c);

/*
    Calculates NW table with two rows and stores a cellpath of scores, identities, gaps and starting and ending positions
*/

struct cell NWscore2rows(char * X, uint64_t Xstart, uint64_t Xend, char * Y, uint64_t Ystart, uint64_t Yend, int iGap, int eGap, struct cell * mc, struct cell * f0, struct cell * f1);