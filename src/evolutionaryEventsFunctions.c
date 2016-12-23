#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <inttypes.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "structs.h"
#include "commonFunctions.h"


void map_frags_to_genomes(unsigned char ** map_table, struct FragFile * frags, uint64_t n_frags){

	uint64_t i, j, from, to, seq;
	//For all frags
	for(i=0;i<n_frags;i++){

		seq = frags[i].seqX;
		from = frags[i].xStart;
		to = frags[i].xEnd;

		//Map coordinates in frag for seqX, which is always forward
		for(j=from;j<to;j++){
			if(map_table[seq][j] == NOFRAG || map_table[seq][j] == COVERFRAG){
				if(j == from) map_table[seq][j] = OPENFRAG;
				else if(j == to) map_table[seq][j] = CLOSEFRAG;
				else map_table[seq][j] = COVERFRAG;
			}
		}

		//Map coordinates in frag for seqY, which might be reversed
		//Remember RAMGECKO coordinates are global respective to forward and with Ystart > Yend when reversed

		seq = frags[i].seqY;
		from = frags[i].yEnd;
		to = frags[i].yStart;

		for(j=from;j<to;j++){
			if(map_table[seq][j] == NOFRAG || map_table[seq][j] == COVERFRAG){
				if(j == from) map_table[seq][j] = OPENFRAG;
				else if(j == to) map_table[seq][j] = CLOSEFRAG;
				else map_table[seq][j] = COVERFRAG;
			}
		}	
	}

	//At this point all coordinates have been mapped for the current fragments
}