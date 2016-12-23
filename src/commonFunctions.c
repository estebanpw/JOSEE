#include <stdio.h>
#include <stdlib.h>
#include <sys/time.h>
#include <inttypes.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "structs.h"
#include "comparisonFunctions.h"

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

void load_sequences_descriptors(Sequence ** sequences, FILE * lengths_file){

    uint64_t n_files;

    //Calculate number of sequences according to size of lengths file
    fseeko(lengths_file, 0L, SEEK_END);
    n_files = ftello(lengths_file)/sizeof(uint64_t);
    fseeko(lengths_file, 0L, SEEK_SET);

    
    //Allocate heap for sequences struct to hold lengths and ids
    Sequence * st = (Sequence *) malloc(n_files*sizeofSequence());

    if(st == NULL) terror("Could not allocate memory for sequence descriptors");

    //Load sequence data into sequences descriptors
    uint64_t i=0, acum = 0;
    while(i<n_files){
        st[i].id = i;
        st[i].acum = acum;
        if(1 != fread(&st[i].len, sizeof(uint64_t), 1, lengths_file)) terror("Wrong number of sequences or sequence file corrupted");
        acum += st[i].len;
        fprintf(stdout, "[INFO] Sequence %"PRIu64" has length %"PRIu64"\n", i, st[i].len);
        i++;
    }

    //Copy pointer direction to allocated heap
    *sequences = st;

}

void load_fragments_local(FILE * fragsfile, uint64_t * n_frags, Sequence * sequences, struct FragFile ** loaded_frags){
    
    struct FragFile temp_frag;
    uint64_t unused_len, total_frags;

    //Compute number of fragments in file
    fseeko(fragsfile, 0L, SEEK_END);
    total_frags = ftello(fragsfile) - 2*(sizeof(uint64_t)); //Remove the headers
    fseeko(fragsfile, 0L, SEEK_SET);

    //Divide by size of frag to get the number of fragments
    //Plus one because it might have padding, thus rounding up to bottom and missing 1 struct
    total_frags = 1 + total_frags/sizeofFragment(); 
    
    //Allocate memory for all frags
    struct FragFile * temp_frags_array = (struct FragFile *) malloc(total_frags * sizeofFragment());
    if(temp_frags_array == NULL) terror("Could not allocate heap for all fragments");

    fprintf(stdout, "[INFO] There are %"PRIu64" fragments to be loaded, requiring %"PRIu64" Megabyte(s) of RAM\n", total_frags, (total_frags*sizeof(struct FragFile))/(1024*1024));

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
        memcpy(&temp_frags_array[*n_frags], &temp_frag, sizeofFragment());
        
        *n_frags = *n_frags + 1;

        if(*n_frags > total_frags){ terror("Something went wrong. More fragments than expected");}
    }

    //Copy pointer of array heap
    *loaded_frags = temp_frags_array;
}