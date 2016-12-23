#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <pthread.h>
#include <float.h>
#include "structs.h"
#include "commonFunctions.h"
#include "comparisonFunctions.h"



int main(int ac, char **av) {
    if (ac < 3) {
        terror("USE: JOSEE <multicomp.frags> <out.csv>");
    }
    //The number of sequences in the comparison, number of comparisons, accumulated length of sequences
    // to transform global coordinates to local
    uint64_t n_files, n_comps, acum;
    uint64_t i;

    //Open frags file and lengths file
    FILE * frags_file, * lengths_file;
    frags_file = fopen64(av[1], "rb");
    if(frags_file == NULL) terror("Could not open input frags file");

    //Concat .lengths to path of multifrags
    char path_lengths[READLINE];
    strcpy(path_lengths, av[1]);
    strcat(path_lengths, ".lengths");
    lengths_file = fopen64(path_lengths, "rb");
    if(lengths_file == NULL) terror("Could not open input lengths file");


    //Calculate number of sequences according to size of lengths file
    fseeko(lengths_file, 0L, SEEK_END);
    n_files = ftello(lengths_file)/sizeof(uint64_t);
    fseeko(lengths_file, 0L, SEEK_SET);
    n_comps = (n_files * (n_files-1)) / 2;

    
    //Allocate heap for sequences struct to hold lengths and ids
    Sequence * comparison_sequences = (Sequence *) malloc(n_files*sizeof(Sequence));
    if(comparison_sequences == NULL) terror("Could not allocate memory for sequence descriptors");

    //Load sequence data into sequences descriptors
    i=0; acum = 0;
    while(i<n_files){
        comparison_sequences[i].id = i;
        comparison_sequences[i].acum = acum;
        if(1 != fread(&comparison_sequences[i].len, sizeof(uint64_t), 1, lengths_file)) terror("Wrong number of sequences or sequence file corrupted");
        acum += comparison_sequences[i].len;
        fprintf(stdout, "[INFO] Sequence %"PRIu64" has length %"PRIu64"\n", i, comparison_sequences[i].len);
        i++;
    }


    


    fclose(lengths_file);



    fclose(frags_file);
    free(comparison_sequences);

    
    return 0;
}