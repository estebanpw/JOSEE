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
    
    //Iterator
    uint64_t i;
    //Number of frags loaded
    uint64_t total_frags;
    //Array of fragments to hold them all
    struct FragFile * loaded_frags = NULL;
    //Clocks to measure time
    clock_t begin, end;


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

    //The genome descriptor 
    Sequence * sequences = NULL;

    //Load lengths and substract accumulated length
    begin = clock();
    load_sequences_descriptors(&sequences, lengths_file);
    end = clock();
    fprintf(stdout, "[INFO] Loading sequence descriptors. T = %e\n", (double)(end-begin)/CLOCKS_PER_SEC);
    fclose(lengths_file);


    //Load fragments into array
    begin = clock();
    load_fragments_local(frags_file, &total_frags, sequences, &loaded_frags);
    end = clock();
    fprintf(stdout, "[INFO] Loading fragments into memory. T = %e\n", (double)(end-begin)/CLOCKS_PER_SEC);
    fclose(frags_file);


        



    
    free(sequences);

    
    return 0;
}