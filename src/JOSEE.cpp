#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <ctype.h>
#include <pthread.h>
#include <float.h>
#include "structs.h"
#include "commonFunctions.h"
#include "comparisonFunctions.h"
#include "evolutionaryEventsFunctions.h"

int DEBUG_ACTIVE = 0;
int HARD_DEBUG_ACTIVE = 0;

void init_args(int argc, char ** av, FILE ** multifrags, FILE ** out_file,
    uint64_t * min_len_trimming, uint64_t * min_trim_itera, char * path_frags, uint64_t * ht_size);

int main(int ac, char **av) {
    
    
    //Iterator
    uint64_t i;
    //Number of frags loaded, number of sequences loaded
    uint64_t total_frags, n_files;
    //Array of fragments to hold them all, pointer to free when changing pointer
    struct FragFile * loaded_frags = NULL, * aux_pointer = NULL;
    //Clocks to measure time
    clock_t begin, end;
    //Minimum length to fragment to accept a trimming
    uint64_t min_len = 100; //Default
    //The genome descriptor 
    Sequence * sequences = NULL;
    //The number of iterations to trimm
    uint64_t N_ITERA = 4; //Default
    //Path to the multifrags file
    char multifrags_path[512];
    multifrags_path[0] = '\0';
    //Initial hash table size (divisor of the longest sequence)
    uint64_t ht_size = 100; //Default


    //Open frags file, lengths file and output files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FILE * frags_file, * lengths_file, * out_file;
    init_args(ac, av, &frags_file, &out_file, &min_len, &N_ITERA, multifrags_path, &ht_size);

    //Concat .lengths to path of multifrags %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    char path_lengths[READLINE];
    path_lengths[0]='\0';
    strcpy(path_lengths, multifrags_path);
    strcat(path_lengths, ".lengths");
    lengths_file = fopen64(path_lengths, "rb");
    if(lengths_file == NULL) terror("Could not open input lengths file");


    //Load lengths and substract accumulated length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    n_files = load_sequences_descriptors(&sequences, lengths_file);
    end = clock();
    fprintf(stdout, "[INFO] Loaded sequence descriptors. T = %e\n", (double)(end-begin)/CLOCKS_PER_SEC);
    fclose(lengths_file);


    //Load fragments into array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    load_fragments_local(frags_file, &total_frags, sequences, &loaded_frags);
    end = clock();
    fprintf(stdout, "[INFO] Loaded fragments into memory. T = %e\n", (double)(end-begin)/CLOCKS_PER_SEC);
    fclose(frags_file); 


    //Initial mapping of fragments to table of genomes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    unsigned char ** map_table = (unsigned char **) std::calloc(n_files, sizeof(unsigned char *));
    //Allocate map table
    for(i=0; i<n_files; i++){ 
        map_table[i] = (unsigned char *) std::calloc(sequences[i].len, sizeof(unsigned char)); 
        if(map_table[i] == NULL) terror("Could not allocate map table"); 
    }
    map_frags_to_genomes(map_table, loaded_frags, total_frags);
    //Compute initial coverage
    get_coverage_from_genome_grid(map_table, sequences, n_files, min_len);
    print_sequences_data(sequences, n_files);
    end = clock();
    fprintf(stdout, "[INFO] Initial mapping completed. T = %e\n", (double)(end-begin)/CLOCKS_PER_SEC);


    //Trimming %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    for(i=0;i<N_ITERA;i++){
        aux_pointer = trim_fragments_and_map(map_table, loaded_frags, &total_frags, min_len, sequences);
        free(loaded_frags); //A new list is being allocated in the function
        loaded_frags = aux_pointer;
    }
    //Compute final coverage
    get_coverage_from_genome_grid(map_table, sequences, n_files, min_len);
    print_sequences_data(sequences, n_files);
    end = clock();
    fprintf(stdout, "[INFO] Trimming of fragments completed after %"PRIu64" iteration(s).\n       Number of final fragments: %"PRIu64". T = %e\n", N_ITERA, total_frags, (double)(end-begin)/CLOCKS_PER_SEC);

    
    //Frags to blocks conversion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    memory_pool * mp = new memory_pool(MAX_MEM_POOLS);
    uint64_t max_len_sequence = get_maximum_length(sequences, n_files);
    uint64_t coord_aux;
    hash_table * ht = new hash_table(mp, max_len_sequence/ht_size, sequences, max_len_sequence);
    for(i=0;i<total_frags;i++){
        //Switch coordinates of reversed fragments. This can only be done at the end of trimming and not meanwhile!
        if(loaded_frags[i].strand == 'r'){ coord_aux = loaded_frags[i].yStart; loaded_frags[i].yStart = loaded_frags[i].yEnd; loaded_frags[i].yEnd = coord_aux;}
        ht->insert_block(&loaded_frags[i]);
    }
    compute_order_of_blocks(ht, n_files);
    end = clock();
    fprintf(stdout, "[INFO] Insertion of fragments into hash table completed. Load factor = %e. T = %e\n", ht->get_load_factor(), (double)(end-begin)/CLOCKS_PER_SEC);
    
    
    //Generate synteny blocks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    Synteny_list * synteny_block_list = compute_synteny_list(ht, n_files, mp);
    //traverse_synteny_list(synteny_block_list);
    end = clock();
    fprintf(stdout, "[INFO] Generated synteny blocks. T = %e\n", (double)(end-begin)/CLOCKS_PER_SEC);
    
    //Start detecting evolutionary events %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    detect_evolutionary_event(synteny_block_list, n_files);
    end = clock();
    fprintf(stdout, "[INFO] Finished detecting evolutionary events. T = %e\n", (double)(end-begin)/CLOCKS_PER_SEC);
    
    // Debug %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    /*
    print_maptable_portion(map_table, 208657, 209236, 60, 1);
    print_maptable_portion(map_table, 208576, 209155, 60, 3);
    print_maptable_portion(map_table, 208576, 209155, 60, 2);
    printf("the not seen\n");
    print_maptable_portion(map_table, 208419, 208998, 60, 0);
    print_maptable_portion(map_table, 208657, 209236, 60, 1);
    getchar();
    */
    /*
    printf("Despeerate.: \n");find_fragments_from_maptable(map_table, 208419, 208998, 0, loaded_frags, total_frags);
    getchar();
    */
    /*
    Bucket * b = ht->get_value_at(208419);
    if(b == NULL) terror("Could not find requested");

    Synteny_list * sbl = find_synteny_block_from_block(sbl, &b->b);
    if(sbl != NULL) printSyntenyListNode(sbl);
    getchar();
    */

    if(DEBUG_ACTIVE){
        char write_debug[512];
        sprintf(write_debug, "%s_%"PRIu64"_%"PRIu64".trim.csv", multifrags_path, N_ITERA, min_len);
        write_maptable_to_disk(map_table, n_files, sequences, write_debug);
        
    }
    if(HARD_DEBUG_ACTIVE){
        ht->print_hash_table(2);
    }

    
    
    
    for(i=0; i<n_files; i++){
        std::free(map_table[i]);
    }
    std::free(map_table);
    std::free(sequences);
    std::free(loaded_frags);

    delete ht;
    delete mp;

    
    return 0;
}


void init_args(int argc, char ** av, FILE ** multifrags, FILE ** out_file,
    uint64_t * min_len_trimming, uint64_t * min_trim_itera, char * path_frags, uint64_t * ht_size){
    
    int pNum = 0;
    while(pNum < argc){
        if(strcmp(av[pNum], "--debug") == 0) DEBUG_ACTIVE = 1;
        if(strcmp(av[pNum], "--hdebug") == 0) HARD_DEBUG_ACTIVE = 1;
        if(strcmp(av[pNum], "--help") == 0){
            fprintf(stdout, "USAGE:\n");
            fprintf(stdout, "           JOSEE -multifrags [query] -out [results]\n");
            fprintf(stdout, "OPTIONAL:\n");
            fprintf(stdout, "           -min_len_trimming   [Integer:   0<=X] (default 100)\n");
            fprintf(stdout, "           -min_trim_itera     [Integer:   0<=X] (default 4)\n");
            fprintf(stdout, "           -hash_table_divisor [Integer:   1<=X] (default 100)\n");
            fprintf(stdout, "           --debug     Turns debug on\n");
            fprintf(stdout, "           --help      Shows the help for program usage\n");
            exit(1);
        }
        if(strcmp(av[pNum], "-multifrags") == 0){
            *multifrags = fopen64(av[pNum+1], "rb");
            strncpy(path_frags, av[pNum+1], strlen(av[pNum+1]));
            path_frags[strlen(av[pNum+1])] = '\0';
            if(multifrags==NULL) terror("Could not open multifrags file");
        }
        if(strcmp(av[pNum], "-out") == 0){
            *out_file = fopen64(av[pNum+1], "wt");
            if(out_file==NULL) terror("Could not open output file");
        }
        if(strcmp(av[pNum], "-min_len_trimming") == 0){
            *min_len_trimming = (uint64_t) atoi(av[pNum+1]);
            if(*min_len_trimming < 0) terror("Minimum trimming length must be zero or more");
        }
        if(strcmp(av[pNum], "-min_trim_itera") == 0){
            *min_trim_itera = (uint64_t) atoi(av[pNum+1]);
            if(*min_trim_itera < 0) terror("Minimum number of trimming iterations must be zero or more");
        }
        if(strcmp(av[pNum], "-hash_table_divisor") == 0){
            *ht_size = (uint64_t) atoi(av[pNum+1]);
            if(*ht_size < 1) terror("The hash table divisor must be one at least");
        }
        pNum++;
    }
    
    if(*multifrags==NULL || *out_file==NULL) terror("A frags file and an output file must be specified");
}

