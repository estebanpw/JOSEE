#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <pthread.h>
#include <float.h>
#include "structs.h"
#include "commonFunctions.h"
#include "comparisonFunctions.h"
#include "evolutionaryEventsFunctions.h"

int main(int ac, char **av) {
    if (ac < 5) {
        terror("USE: JOSEE <multicomp.frags> <out.csv> <min_len_trimming> <num_trimm_iterations>");
    }
    
    //Iterator
    uint64_t i;
    //Number of frags loaded, number of sequences loaded
    uint64_t total_frags, n_files;
    //Array of fragments to hold them all, pointer to free when changing pointer
    struct FragFile * loaded_frags = NULL, * aux_pointer = NULL;
    //Clocks to measure time
    clock_t begin, end;
    //Minimum length to fragment to accept a trimming
    uint64_t min_len = (uint64_t) atoi(av[3]);
    //The genome descriptor 
    Sequence * sequences = NULL;
    //The number of iterations to trimm
    uint64_t N_ITERA = (uint64_t) atoi(av[4]);

    //Open frags file and lengths file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FILE * frags_file, * lengths_file;
    frags_file = fopen64(av[1], "rb");
    if(frags_file == NULL) terror("Could not open input frags file");

    //Concat .lengths to path of multifrags %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    char path_lengths[READLINE];
    strcpy(path_lengths, av[1]);
    strcat(path_lengths, ".lengths");
    lengths_file = fopen64(path_lengths, "rb");
    if(lengths_file == NULL) terror("Could not open input lengths file");


    //Load lengths and substract accumulated length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    n_files = load_sequences_descriptors(&sequences, lengths_file);
    end = clock();
    fprintf(stdout, "[INFO] Loading sequence descriptors. T = %e\n", (double)(end-begin)/CLOCKS_PER_SEC);
    fclose(lengths_file);


    //Load fragments into array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    load_fragments_local(frags_file, &total_frags, sequences, &loaded_frags);
    end = clock();
    fprintf(stdout, "[INFO] Loading fragments into memory. T = %e\n", (double)(end-begin)/CLOCKS_PER_SEC);
    fclose(frags_file); 


    //Initial mapping of fragments to table of genomes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    unsigned char ** map_table = (unsigned char **) calloc(n_files, sizeof(unsigned char *));
    //Allocate map table
    for(i=0; i<n_files; i++){ map_table[i] = calloc(sequences[i].len, sizeof(unsigned char)); if(map_table[i] == NULL) terror("Could not allocate map table"); }
    map_frags_to_genomes(map_table, loaded_frags, total_frags);
    end = clock();
    fprintf(stdout, "[INFO] Initial mapping completed. T = %e\n", (double)(end-begin)/CLOCKS_PER_SEC);


    //Trimming %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    for(i=0;i<N_ITERA;i++){
        aux_pointer = trim_fragments_and_map(map_table, loaded_frags, &total_frags, min_len, sequences);
        free(loaded_frags);
        loaded_frags = aux_pointer;
    }
    end = clock();
    fprintf(stdout, "[INFO] Trimming of fragments completed after %"PRIu64" iteration(s).\n       Number of final fragments: %"PRIu64". T = %e\n", N_ITERA, total_frags, (double)(end-begin)/CLOCKS_PER_SEC);
    
    
    for(i=0; i<n_files; i++){
        free(map_table[i]);
    }
    free(map_table);
    free(sequences);
    free(loaded_frags);

    
    return 0;
}