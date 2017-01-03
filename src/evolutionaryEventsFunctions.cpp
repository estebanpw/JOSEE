#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <inttypes.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "structs.h"
#include "commonFunctions.h"
#include "comparisonFunctions.h"

void map_frags_to_genomes(unsigned char ** map_table, struct FragFile * frags, uint64_t n_frags){

	uint64_t i, j, from, to, seq;
	//For all frags
	for(i=0;i<n_frags;i++){

		seq = frags[i].seqX;
		from = frags[i].xStart;
		to = frags[i].xEnd;

		//Map coordinates in frag for seqX, which is always forward
		for(j=from;j<to;j++){
			if(j == from) map_table[seq][j] = OPENFRAG;
			if(j == to) map_table[seq][j] = CLOSEFRAG;
			
			if(map_table[seq][j] == NOFRAG){
				map_table[seq][j] = COVERFRAG;
			}
		}

		//Map coordinates in frag for seqY, which might be reversed
		//Remember RAMGECKO coordinates are global respective to forward and with Ystart > Yend when reversed
		//So reversed should only be switched

		if(frags[i].strand == 'f'){
			seq = frags[i].seqY;
			from = frags[i].yStart;
			to = frags[i].yEnd;
		}else{
			seq = frags[i].seqY;
			from = frags[i].yEnd;
			to = frags[i].yStart;	
		}
		

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



inline void copyFragWithNewCoordinates(struct FragFile * destination, struct FragFile * source, uint64_t xStart, uint64_t yStart, uint64_t xEnd, uint64_t yEnd, uint64_t len){
	destination->diag = source->diag;
    destination->xStart = xStart;
    destination->yStart = yStart;
    destination->xEnd = xEnd;
    destination->yEnd = yEnd;
    destination->length = len;
    destination->seqX = source->seqX;
    destination->seqY = source->seqY;
    destination->strand = source->strand;
}

struct FragFile * trim_fragments_and_map(unsigned char ** map_table, struct FragFile * frags, uint64_t * n_frags, uint64_t min_len, Sequence * sequences){

	struct FragFile new_frag;
	struct FragFile * list_new_frags;
	uint64_t list_reallocs = 1;
	uint64_t new_frags_number = 0;
	uint64_t size_fragment = sizeofFragment(); //To not compute it every time

	//Allocate memory
	list_new_frags = (struct FragFile *) malloc(INIT_TRIM_FRAGS*sizeofFragment());
	if(list_new_frags == NULL) terror("Could not allocate memory for list of new fragments in trimming");

	//Start trimming process
	uint64_t i, jX, jY, fromX, fromY, toX, toY, seqX, seqY, frag_len;
	char strand;
	uint64_t cur_new_len;

	for(i=0; i<*n_frags; i++){
		
		//Copy frag values
		fromX = frags[i].xStart; 
		toX = frags[i].xEnd; 
		

		fromY = frags[i].yStart;
		toY = frags[i].yEnd;	
		strand = frags[i].strand;

		

		seqX = frags[i].seqX; seqY = frags[i].seqY;

		frag_len = frags[i].length;
		cur_new_len = 1;
		jX = fromX+1;
		if(strand == 'f') jY = fromY+1; else jY = fromY-1;

		while(jX < toX){
			//Check how long until there is a break (by starting or ending of frag)
			while(cur_new_len < frag_len && jX < sequences[seqX].len && jY < sequences[seqY].len && jY >= 0){

				if(map_table[seqX][jX] != COVERFRAG) break;
				if(map_table[seqY][jY] != COVERFRAG) break;

				cur_new_len++; jX++;
				if(strand == 'f'){ jY++; }else{ if(jY > 0) jY--; else break;} //To scape the buffer overflow of uints64
			}

			//At this point, jX and jY hold the ending coordinates of the new fragment
			//And fromX and fromY hold the starting coordinates
			if(cur_new_len >= min_len){ //Filtering

				//The fragment must be snipped out and saved
				copyFragWithNewCoordinates(&new_frag, &frags[i], fromX, fromY, jX, jY, cur_new_len);
				memcpy(&list_new_frags[new_frags_number], &new_frag, size_fragment);
				new_frags_number++;
				//Check if we need to realloc the list of new frags
				if(new_frags_number == list_reallocs*INIT_TRIM_FRAGS){
					list_reallocs++;
					list_new_frags = (struct FragFile *) realloc(list_new_frags, list_reallocs*INIT_TRIM_FRAGS*size_fragment);
					if(list_new_frags == NULL) terror("Could not realloc fragments on the trimming process");
				}

				//And set the mapping grid to the new values
				map_table[seqX][jX] = 3;
				map_table[seqY][jY] = 3;
				if(jX+1 < sequences[seqX].len && jX+1 < toX) map_table[seqX][jX+1] = 1;
				if(strand == 'f'){
					if(jY+1 < sequences[seqY].len && jY+1 < toY) map_table[seqY][jY+1] = 1;	
				}else{
					if(jY > 0 && jY-1 < fromY) map_table[seqY][jY-1] = 1;
				}
				
				//Set the fromX and fromY to 1 (start frag) again in case this is not the first time we split

			}
			//If you are here, either the fragment was too short, or was written correctly or we are at the end of the frag
			//Just keep going
			//Copy frag values
			
			fromX = jX+1; //One to move from an ending 3 to an opening 1
			if(strand == 'f') fromY = jY+1; else fromY = jY-1; //Same
			
			//And one more to skip the opening 1
			jX = fromX+1;
			if(strand == 'f') jY = fromY+1; else fromY = jY-1; 

			cur_new_len = 1;

			//End of outside while
		}


	}

	*n_frags = new_frags_number;
	return list_new_frags;

}


void blocks_from_grid(unsigned char ** map_table, struct FragFile * frags, uint64_t * n_frags, Sequence * sequences){
	
}
