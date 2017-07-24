#define __STDC_FORMAT_MACROS

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

uint64_t total_bytes_in_use = 0;
int DEBUG_ACTIVE = 0;
int HARD_DEBUG_ACTIVE = 0;

void print_all();
void init_args(int argc, char ** av, FILE ** multifrags, FILE ** out_file, FILE ** out_file_log,
    uint64_t * min_len_trimming, uint64_t * min_trim_itera, char * path_frags, uint64_t * ht_size,
    FILE ** out_blocks, FILE ** out_breakpoints, char * path_files, char * path_annotations, 
    uint32_t * kmer_size, FILE ** trim_frags_file, bool * trim_frags_file_write, char * file_out_write);
int main(int ac, char **av) {
    
    
    //Iterator
    uint64_t i;
    //Number of frags loaded, number of sequences loaded
    uint64_t total_frags, n_files;
    //Last ID of synteny lists
    uint64_t last_s_id;
    //Array of fragments to hold them all, pointer to free when changing pointer
    struct FragFile * loaded_frags = NULL, * aux_pointer = NULL;
    //Clocks to measure time
    clock_t begin, end;
    //Minimum length to fragment to accept a trimming
    uint64_t min_len = 50; //Default
    //The sequences manager to store ids, lengths, etc
    sequence_manager * seq_manager = new sequence_manager();
    //The number of iterations to trimm
    uint64_t N_ITERA = 500; //Default
    //Path to the multifrags file
    char multifrags_path[512];
    multifrags_path[0] = '\0';
    //Path to the fastas
    char fastas_path[512];
    fastas_path[0] = '\0';
    //Path to annotations
    char path_annotations[512];
    path_annotations[0] = '\0';
    //Initial hash table size (divisor of the longest sequence)
    uint64_t ht_size = 100; //Default
    //The kmer size for blocks alignment 
    uint32_t kmer_size = 16; //Default
    //Default behaviour is the trimmed frags do not already exist
    bool trim_frags_file_write = true; 
    // Fix
    char file_out_char[512];
    file_out_char[0] = '\0';


    //Open frags file, lengths file and output files %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    FILE * frags_file, * lengths_file, * out_file, * out_file_log, * out_blocks = NULL, * out_breakpoints = NULL, * trim_frags_file = NULL;
    init_args(ac, av, &frags_file, &out_file, &out_file_log, &min_len, &N_ITERA, multifrags_path, &ht_size, &out_blocks, &out_breakpoints, fastas_path, path_annotations, &kmer_size, &trim_frags_file, &trim_frags_file_write, file_out_char);

    //Concat .lengths to path of multifrags %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    char path_lengths[READLINE];
    path_lengths[0]='\0';
    strcpy(path_lengths, multifrags_path);
    strcat(path_lengths, ".lengths");
    lengths_file = fopen64(path_lengths, "rb");
    if(lengths_file == NULL) terror("Could not open input lengths file");


    //Set path to annotations (if given)
    if(path_annotations[0] != '\0') seq_manager->set_path_annotations(path_annotations);

    //Load lengths and substract accumulated length %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    n_files =  seq_manager->load_sequences_descriptors(lengths_file);
    end = clock();
    fprintf(stdout, "[INFO] Loaded sequence descriptors. T = %e\n", (double)(end-begin)/CLOCKS_PER_SEC);
    fclose(lengths_file);


    //Load fragments into array %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    load_fragments_local(frags_file, &total_frags, &loaded_frags);
    end = clock();
    print_memory_usage();
    fprintf(stdout, "[INFO] Loaded fragments into memory. T = %e\n", (double)(end-begin)/CLOCKS_PER_SEC);
    fclose(frags_file); 


    //Initial mapping of fragments to table of genomes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    unsigned char ** map_table = (unsigned char **) std::calloc(n_files, sizeof(unsigned char *));
    total_bytes_in_use += n_files*sizeof(unsigned char *);
    //Allocate map table
    for(i=0; i<n_files; i++){ 
        map_table[i] = (unsigned char *) std::calloc(seq_manager->get_sequence_by_label(i)->len, sizeof(unsigned char)); 
        total_bytes_in_use += seq_manager->get_sequence_by_label(i)->len * sizeof(unsigned char);
        if(map_table[i] == NULL) terror("Could not allocate map table"); 
    }
    map_frags_to_genomes(map_table, loaded_frags, total_frags, seq_manager);
    //Compute initial coverage
    get_coverage_from_genome_grid(map_table, seq_manager, n_files, min_len);
    seq_manager->print_sequences_data();
    end = clock();
    print_memory_usage();
    fprintf(stdout, "[INFO] Initial mapping completed. T = %e\n", (double)(end-begin)/CLOCKS_PER_SEC);

    //Trimming %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    uint64_t ratio_itera = N_ITERA/5;
    //Check if frags had already been computed and exist 
    if(trim_frags_file != NULL && trim_frags_file_write == false){
        // Load fragments trimmed
        total_bytes_in_use -= total_frags * sizeofFragment();
        if(0 == fread(&total_frags, sizeof(uint64_t), 1, trim_frags_file)) terror("Incorrect number of fragments to load");
        loaded_frags = (struct FragFile *) std::realloc(loaded_frags, total_frags*sizeofFragment());
        total_bytes_in_use += total_frags * sizeofFragment();
        if(loaded_frags == NULL) terror("Could not load existing trimmed frags");
        if(0 == fread(loaded_frags, sizeofFragment(), total_frags, trim_frags_file)) terror("Specified fragments are of size zero");
    }else{
        // Compute fragments trim
        total_bytes_in_use -= total_frags * sizeofFragment();
        for(i=0;i<N_ITERA;i++){
            if(i % ratio_itera == 0) fprintf(stdout, "[INFO] Iteration %"PRIu64"\n", i);
            aux_pointer = trim_fragments_and_map(map_table, loaded_frags, &total_frags, min_len);
            //fprintf(stdout, "FRAGS: %"PRIu64"\n", total_frags);
            free(loaded_frags); //A new list is being allocated in the function
            loaded_frags = aux_pointer;
        }
        total_bytes_in_use += total_frags * sizeofFragment();
    }
    
    
    //Compute final coverage
    get_coverage_from_genome_grid(map_table, seq_manager, n_files, min_len);
    //Write frags if it was set
    if(trim_frags_file != NULL && trim_frags_file_write == true){
        //Write number of frags and frags
        fwrite(&total_frags, sizeof(uint64_t), 1, trim_frags_file);
        fwrite(loaded_frags, sizeofFragment(), total_frags, trim_frags_file);
    }
    //seq_manager->print_sequences_data();
    end = clock();
    print_memory_usage();
    fprintf(stdout, "[INFO] Trimming of fragments completed after %"PRIu64" iteration(s).\n       Number of final fragments: %"PRIu64". T = %e\n", N_ITERA, total_frags, (double)(end-begin)/CLOCKS_PER_SEC);

    
    //Frags to blocks conversion %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    memory_pool * mp = new memory_pool(POOL_SIZE);
    uint64_t max_len_sequence = seq_manager->get_maximum_length();
    uint64_t coord_aux;
    // ram debug over here
    hash_table * ht = new hash_table(mp, max_len_sequence/ht_size, seq_manager, max_len_sequence);
    for(i=0;i<total_frags;i++){
        //Switch coordinates of reversed fragments. This can only be done at the end of trimming and not meanwhile!
        if(loaded_frags[i].strand == 'r'){ coord_aux = loaded_frags[i].yStart; loaded_frags[i].yStart = loaded_frags[i].yEnd; loaded_frags[i].yEnd = coord_aux;}
        ht->insert_block(&loaded_frags[i]);
    }
    //ht->print_hash_table(2);

    compute_order_of_blocks(ht, n_files);
    end = clock();
    //fprintf(stdout, "[INFO] Insertion of fragments into hash table completed. Load factor = %e. T = %e\n", ht->get_load_factor(), (double)(end-begin)/CLOCKS_PER_SEC);
    print_memory_usage();
    fprintf(stdout, "[INFO] Insertion of fragments into hash table completed. T = %e\n", (double)(end-begin)/CLOCKS_PER_SEC);
     
    
    
    //Generate synteny blocks %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    Synteny_list * synteny_block_list = compute_synteny_list(ht, n_files, mp, &last_s_id);
    traverse_synteny_list_and_write(synteny_block_list, n_files, "init");
    //traverse_synteny_list(synteny_block_list);
    end = clock();
    print_memory_usage();
    fprintf(stdout, "[INFO] Generated synteny blocks. T = %e\n", (double)(end-begin)/CLOCKS_PER_SEC);
    

    //Load sequences to run pairwise alignment %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    seq_manager->read_dna_sequences(fastas_path);
    //And load annotations if provided
    if(seq_manager->get_path_annotations() != NULL) seq_manager->read_annotations();
    end = clock();
    print_memory_usage();
    fprintf(stdout, "[INFO] Loaded DNA sequences. T = %e\n", (double)(end-begin)/CLOCKS_PER_SEC);

    //Write blocks and breakpoints to file %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if(out_blocks != NULL && out_breakpoints != NULL){
        begin = clock();
        ht->write_blocks_and_breakpoints_to_file(out_blocks, out_breakpoints);
        end = clock();
        fprintf(stdout, "[INFO] Wrote blocks and breakpoints to output files. T = %e\n", (double)(end-begin)/CLOCKS_PER_SEC);
    }



    //Start detecting evolutionary events %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    begin = clock();
    detect_evolutionary_event(synteny_block_list, seq_manager, kmer_size, ht, &last_s_id, out_file_log, file_out_char, mp);
    //fprintf(stdout, "\t\t After applying EE(s)\n");
    //traverse_synteny_list(synteny_block_list);
    end = clock();
    print_memory_usage();
    fprintf(stdout, "[INFO] Finished detecting evolutionary events. T = %e\n", (double)(end-begin)/CLOCKS_PER_SEC);
    


    //Write output sequences and finish %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    for(i=0;i<n_files;i++){
        fprintf(out_file, ">seq %"PRIu64" with len %"PRIu64"\n", i, seq_manager->get_sequence_by_label(i)->len);
        seq_manager->print_sequence_region(out_file, i, 0, seq_manager->get_sequence_by_label(i)->len);
    }
    

    // DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    // DEBUG %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    //ht->print_hash_table(2);

    
    
    //print_maptable_portion(map_table, 0, 194, 50, 0);
    //print_maptable_portion(map_table, 0, 231, 50, 1);
    
    /*
    print_maptable_portion(map_table, 833568  , 833983, 60, 1);
    print_maptable_portion(map_table, 834344  , 834759, 60, 2);
    print_maptable_portion(map_table, 834361  , 834776, 60, 3);
    */
    /*
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
        write_maptable_to_disk(map_table, seq_manager, write_debug);
        
    }
    if(HARD_DEBUG_ACTIVE){
        ht->print_hash_table(2);
    }

    
    
    
    for(i=0; i<n_files; i++){
        std::free(map_table[i]);
    }
    std::free(map_table);
    std::free(loaded_frags);

    fclose(out_file);
    fclose(out_file_log);
    if(out_blocks != NULL){
        fclose(out_blocks);
        fclose(out_breakpoints);
    }
    if(trim_frags_file != NULL){
        fclose(trim_frags_file);
    }
    

    delete ht;
    delete mp;
    delete seq_manager;

    
    return 0;
}

void print_all(){
    fprintf(stdout, "USAGE:\n");
    fprintf(stdout, "           JOSEE -multifrags [query] -pathfiles [fastas] -out [results]\n");
    fprintf(stdout, "OPTIONAL:\n");
    fprintf(stdout, "           -min_len_trimming   [Integer:   0<=X] (default 50)\n");
    fprintf(stdout, "           -min_trim_itera     [Integer:   0<=X] (default 500)\n");
    fprintf(stdout, "           -hash_table_divisor [Integer:   1<=X] (default 100)\n");
    fprintf(stdout, "           -kmer               [Integer:   1<=X] (default 16)\n");
    fprintf(stdout, "           -write_blocks_bps   [Path without format extension]\n");
    fprintf(stdout, "           -reuse_trim_frags   [Path to destination output file]\n");
    fprintf(stdout, "                               Notice that this will save the file if\n");
    fprintf(stdout, "                               it does not exist and load it if it does\n");
    fprintf(stdout, "           -annotations        [Path without format extension]\n");
    fprintf(stdout, "           --debug     Turns debug on\n");
    fprintf(stdout, "           --help      Shows the help for program usage\n");
}

void init_args(int argc, char ** av, FILE ** multifrags, FILE ** out_file, FILE ** out_file_log,
    uint64_t * min_len_trimming, uint64_t * min_trim_itera, char * path_frags, uint64_t * ht_size,
    FILE ** out_blocks, FILE ** out_breakpoints, char * path_files, char * path_annotations, 
    uint32_t * kmer_size, FILE ** trim_frags_file, bool * trim_frags_file_write, char * file_out_char){
    
    int pNum = 0;
    while(pNum < argc){
        if(strcmp(av[pNum], "--debug") == 0) DEBUG_ACTIVE = 1;
        if(strcmp(av[pNum], "--hdebug") == 0) HARD_DEBUG_ACTIVE = 1;
        if(strcmp(av[pNum], "--help") == 0){
            print_all();
            exit(1);
        }
        if(strcmp(av[pNum], "-multifrags") == 0){
            *multifrags = fopen64(av[pNum+1], "rb");
            strncpy(path_frags, av[pNum+1], strlen(av[pNum+1]));
            path_frags[strlen(av[pNum+1])] = '\0';
            if(multifrags==NULL) terror("Could not open multifrags file");
        }
        if(strcmp(av[pNum], "-reuse_trim_frags") == 0){
            if(exists_file(av[pNum+1]) == 0){
                *trim_frags_file = fopen64(av[pNum+1], "wb");
                *trim_frags_file_write = true;
                if(trim_frags_file == NULL) terror("Could not open trimmed frags file");
            }else{
                *trim_frags_file = fopen64(av[pNum+1], "rb");
                *trim_frags_file_write = false;
                if(trim_frags_file == NULL) terror("Could not open existing trimmed frags file");
            }
        }
        if(strcmp(av[pNum], "-pathfiles") == 0){
            strncpy(path_files, av[pNum+1], strlen(av[pNum+1]));
            path_files[strlen(av[pNum+1])] = '\0';
        }
        if(strcmp(av[pNum], "-kmer") == 0){
            *kmer_size = (uint32_t) atoi(av[pNum+1]);
            if(*kmer_size < 1) terror("Minimum kmer size must be larger than one");
        }
        if(strcmp(av[pNum], "-annotations") == 0){
            strncpy(path_annotations, av[pNum+1], strlen(av[pNum+1]));
            path_annotations[strlen(av[pNum+1])] = '\0';
        }
        if(strcmp(av[pNum], "-write_blocks_bps") == 0){
            char templine[READLINE]; templine[0] = '\0';
            strcpy(templine, av[pNum+1]);
            strcat(templine, ".blocks");
            *out_blocks = fopen64(templine, "wt");

            templine[0] = '\0';
            strcpy(templine, av[pNum+1]);
            strcat(templine, ".breakpoints");
            *out_breakpoints = fopen64(templine, "wt");
            if(*out_blocks==NULL || *out_breakpoints==NULL) terror("Could not open blocks/breakpoints file");
        }
        if(strcmp(av[pNum], "-out") == 0){
            *out_file = fopen64(av[pNum+1], "wt");
            char tmp[MAX_LINE];
            sprintf(&tmp[0], "%s.log", av[pNum+1]);
            strcpy(file_out_char, tmp); 
            *out_file_log = fopen64(tmp, "wt");
            if(out_file==NULL) terror("Could not open output file");
            if(out_file_log==NULL) terror("Could not open output log file");
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
    
    if(*multifrags==NULL || *out_file==NULL || path_files[0] == '\0'){
        print_all();
        terror("A frags file, a path to the fasta files and an output file must be specified");
    }
}

