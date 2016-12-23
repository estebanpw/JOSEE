#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <inttypes.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "structs.h"

void terror(const char *s) {
    printf("ERR**** %s ****\n", s);
    exit(-1);
}

char buffered_fgetc(char *buffer, uint64_t *pos, uint64_t *read, FILE *f) {
    if (*pos >= READBUF) {
        *pos = 0;
        memset(buffer, 0, READBUF);
        *read = fread(buffer, 1, READBUF, f);
    }
    *pos = *pos + 1;
    return buffer[*pos-1];
}

void load_fragments_local(FILE * fragsfile, uint64_t * n_frags, Sequence * sequences, struct FragFile * loaded_frags){
    
    struct FragFile temp_frag;
    uint64_t unused_len, total_frags;

    //Compute number of fragments in file
    fseeko(fragsfile, 0L, SEEK_END);
    total_frags = ftello(fragsfile) - 2*(sizeof(uint64_t)); //Remove the headers
    fseeko(fragsfile, 0L, SEEK_SET);

    total_frags = total_frags/sizeof(struct FragFile); //Divide by size of frag to get the number of fragments
    
    //Allocate memory for all frags
    loaded_frags = (struct FragFile *) malloc(total_frags * sizeof(struct FragFile));
    if(loaded_frags == NULL) terror("Could not allocate heap for all fragments");


    //Skip headers
    readSequenceLength(&unused_len, fragsfile);
    readSequenceLength(&unused_len, fragsfile);

    //To keep track of current frag
    *n_frags = 0;

    while(!feof(fragsfile)){
        readFragment(&temp_frag, fragsfile);
        
        //Transform coordinates to local
        temp_frag.xStart = temp_frag.xStart - sequences[temp_frag.seqX].acum;
        temp_frag.xEnd = temp_frag.xEnd - sequences[temp_frag.seqX].acum;

        temp_frag.yEnd = temp_frag.yEnd - sequences[temp_frag.seqY].acum;
        temp_frag.yEnd = temp_frag.yEnd - sequences[temp_frag.seqY].acum;

        //Copy temp fragment into array
        memcpy(&loaded_frags[*n_frags], &temp_frag, sizeof(struct FragFile));

        *n_frags = *n_frags + 1;
    }
}