#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <ctype.h>
#include <pthread.h>
#include <float.h>
#include "structs.h"
#include "commonFunctions.h"
#include "comparisonFunctions.h"

#define PRINT_RATE 1000


/*
    Case 0 is default (output all)
    Case 1 is overlapping
    Case 2 is centered +- N bases

*/

uint64_t getNumberOfSequences(FILE * f){
    uint64_t t_seqs = 0;
    char c = 'Z';
    while(!feof(f)){
        if(c=='\n') t_seqs++;
        c=fgetc(f);
    }
    //Rewind
    fseeko64(f, 0L, SEEK_SET);
    return t_seqs;
}

//Function to read all files from a text list
char ** read_all_vs_all_files(char * list_of_files, uint64_t * n_files){

    FILE * lf = fopen64(list_of_files, "rt");
    if(lf == NULL) terror("Could not open list of genomic files");

    *n_files = getNumberOfSequences(lf);
    uint64_t i;

    char ** all_sequences = (char **) std::malloc (*n_files*sizeof(char *));
    for(i=0;i<*n_files;i++){
        all_sequences[i] = (char *) std::malloc(READLINE*sizeof(char));
        if(all_sequences[i] == NULL) terror("Could not allocate paths to files");
    }


    i = 0;
    while(i < *n_files && !feof(lf)){
        if(fgets(all_sequences[i], READLINE, lf) > 0){
            if(all_sequences[i][0] != '\0' && all_sequences[i][0] != '\n'){
                all_sequences[i][strlen(all_sequences[i])-1] = '\0';
        
                i++;
            }
        }
        
    }
    *n_files = i;
    fclose(lf);
    return all_sequences;
}


int main(int ac, char **av) {
    if (ac < 7) {
        terror("USE: cutter <path_list> <out_dna> <out_class> <file_blocks_breakpoints> <min_len_filter> <w_mode [0/1]>\nUse working mode = 0 for classic output; 1 to output blocks from middle point and breakpoints from endings of blocks.");
    }

    //Iterators
    uint64_t i, n_files, curr_pos;

    //Load files to compare
    char ** paths_to_files = read_all_vs_all_files(av[1], &n_files);
    if(n_files < 2) terror("At least two files need to be loaded");
    for(uint64_t k=0;k<n_files;k++){ fprintf(stdout, "[INFO] File %"PRIu64": %s\n", k, paths_to_files[k]);}

    //Files to write results
    FILE * dna_out = fopen64(av[2], "wt"); if(dna_out == NULL) terror("Could not open output dna file");
    FILE * dna_class = fopen64(av[3], "wt"); if(dna_class == NULL) terror("Could not open class output file");
    char blocks_bps[READLINE]; blocks_bps[0]='\0';
    strcpy(blocks_bps, av[4]);
    strcat(blocks_bps, ".blocks");
    fprintf(stdout, "[INFO] Opening %s\n", blocks_bps);
    FILE * blocks = fopen64(blocks_bps, "rt"); if(blocks == NULL) terror("Could not open input blocks file");
    blocks_bps[0]='\0';
    strcpy(blocks_bps, av[4]);
    strcat(blocks_bps, ".breakpoints");
    fprintf(stdout, "[INFO] Opening %s\n", blocks_bps);
    FILE * breakpoints = fopen64(blocks_bps, "rt"); if(breakpoints == NULL) terror("Could not open breakpoints file");

    //Min length to write blocks and breakpoints
    uint64_t min_len_filter = (uint64_t) atoi(av[5]);

    //Working mode
    int w_mode = atoi(av[6]);

    //The file that will open all sequence files
    FILE * current;

    //Vector to tell for sequence reallocs
    uint64_t * n_reallocs = (uint64_t *) std::calloc(n_files, sizeof(uint64_t));
    if(n_reallocs == NULL) terror("Could not allocate realloc count vector");

    //Char to hold all sequences
    char ** all_sequences = (char **) std::calloc(n_files, sizeof(char *));
    if(all_sequences == NULL) terror("Could not allocate sequences pointer");

    //Sizes for each sequence
    uint64_t * seq_sizes = (uint64_t *) std::calloc(n_files, sizeof(uint64_t));
    if(seq_sizes == NULL) terror("Could not allocate sequence sizes");

    //Read using buffered fgetc
    uint64_t idx = 0, r = 0;
    char * temp_seq_buffer = NULL;
    if ((temp_seq_buffer = (char *) std::calloc(READBUF, sizeof(char))) == NULL) {
        terror("Could not allocate memory for read buffer");
    }
    //To force reading from the buffer
    idx = READBUF + 1;
    char c;
    
    //Read sequences and load into array
    for(i=0;i<n_files;i++){
        current = fopen64(paths_to_files[i], "rt");
        all_sequences[i] = (char *) std::calloc(SEQ_REALLOC, sizeof(char));
        if(all_sequences[i] == NULL) terror("Could not allocate genome sequence");
        if(current == NULL) terror("Could not open fasta file");

        curr_pos = 0;
        idx = READBUF + 1;
        r = 0;

        c = buffered_fgetc(temp_seq_buffer, &idx, &r, current);
        while((!feof(current) || (feof(current) && idx < r))){

            if(c == '>'){
                while(c != '\n') c = buffered_fgetc(temp_seq_buffer, &idx, &r, current); //Skip id

                while(c != '>' && (!feof(current) || (feof(current) && idx < r))){ //Until next id
                    c = buffered_fgetc(temp_seq_buffer, &idx, &r, current);
                    c = toupper(c);
                    if(c >= 'A' && c <= 'Z'){
                        all_sequences[i][curr_pos++] = c;
                        if(curr_pos >= SEQ_REALLOC*n_reallocs[i]){
                            n_reallocs[i]++;
                            all_sequences[i] = (char *) std::realloc(all_sequences[i], n_reallocs[i]*SEQ_REALLOC);
                            if(all_sequences[i] == NULL) terror("Could not realloc sequence");
                        }
                    }
                }
                curr_pos++; //one for the *
            }else{
                c = buffered_fgetc(temp_seq_buffer, &idx, &r, current);
            }
            
        }
        seq_sizes[i] = curr_pos;
        fclose(current);
    }


    std::free(temp_seq_buffer);

    for(i=0;i<n_files;i++){
        std::free(paths_to_files[i]);
    }
    std::free(paths_to_files);


    //At this point all sequences are loaded
    uint64_t b_number, b_start, b_end, b_order, b_sequence, b_len;
    char header[READLINE];
    char * seq_region = (char *) std::malloc(SEQ_REALLOC*sizeof(char));
    if(seq_region == NULL) terror("Could not allocate sequence region to output");
    uint64_t seq_region_reallocs = 1;

    //Number of breakpoints counted
    uint64_t n_breakpoints = 0, n_blocks = 0;



    //Blocks later to have the ratio, since there are more breakpoints on average
    //Skip header
    if(NULL == fgets(header, READLINE, blocks)) terror("No header was read or empty blocks file");
    while(!feof(blocks)){
        //Read a line i.e. a block
        if(6 == fscanf(blocks, "%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"", &b_number, &b_sequence, &b_order, &b_start, &b_end, &b_len)){
            //Check that the region fits in vector
            if(b_len >= seq_region_reallocs*SEQ_REALLOC){
                seq_region_reallocs++;
                seq_region = (char *) std::realloc(seq_region, seq_region_reallocs*SEQ_REALLOC);
                if(seq_region == NULL) terror("Could not realloc region sequence");
            }

            switch(w_mode){
                
                case 2: //Center cutting
                {
                    if(b_len >= min_len_filter){//Only if it is long enough
                        //Copy the region to buffer
                        if(b_number % PRINT_RATE == 0) fprintf(stdout, "[INFO] On block %"PRIu64"\n", b_number);
                        uint64_t midpoint = b_start + b_len/2;
                        memcpy(&seq_region[0], all_sequences[b_sequence]+midpoint-min_len_filter, 2*min_len_filter);
                        seq_region[2*min_len_filter+1] = '\0';
                        fprintf(dna_out, "%s\n", seq_region);
                        fprintf(dna_class, "1\n"); //Its a block
                        n_blocks++;
                    }  
                }
                break;
                case 1:
                {
                    //if the block is longer than 2 times the minimum filter
                    //We want the same number of breakpoints and blocks
                    //OVERLAPPING
                    while(b_end - b_start >= 2*min_len_filter){
                        //Copy 2 times len filter
                        memcpy(&seq_region[0], all_sequences[b_sequence]+b_start+(2*min_len_filter), 2*min_len_filter);
                        seq_region[2*min_len_filter] = '\0';
                        fprintf(dna_out, "%s\n", seq_region);
                        fprintf(dna_class, "1\n"); //Its a block
                        n_blocks++;
                        //And change coordinates
                        b_start = b_start + 2*min_len_filter;
                    }
                }
                break;
                default:
                {
                    if(b_len >= min_len_filter){//Only if it is long enough
                        //Copy the region to buffer
                        if(b_number % PRINT_RATE == 0) fprintf(stdout, "[INFO] On block %"PRIu64"\n", b_number);
                        memcpy(&seq_region[0], all_sequences[b_sequence]+b_start, b_len);
                        seq_region[b_end-b_start+1] = '\0';
                        fprintf(dna_out, "%s\n", seq_region);
                        fprintf(dna_class, "1\n"); //Its a block
                        n_blocks++;
                    }  
                }
            } 
        }
    }







    //Do breakpoints first
    //Skip header
    if(NULL == fgets(header, READLINE, breakpoints)) terror("No header was read or empty breakpoints file");
    while(!feof(breakpoints)){
        //Read a line i.e. a block
        
        if(5 == fscanf(breakpoints, "%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"", &b_number, &b_sequence, &b_start, &b_end, &b_len)){
            //Check that the region fits in vector
            if(b_len >= seq_region_reallocs*SEQ_REALLOC){
                seq_region_reallocs++;
                seq_region = (char *) std::realloc(seq_region, seq_region_reallocs*SEQ_REALLOC);
                if(seq_region == NULL) terror("Could not realloc region sequence");
            }

            switch(w_mode){

                //Center cutting
                case 2:
                {
                    if(n_breakpoints < n_blocks && b_len >= min_len_filter){ //Only if it is long enough
                        //Copy the region to buffer
                        if(b_number % PRINT_RATE == 0) fprintf(stdout, "[INFO] On breakpoint %"PRIu64"\n", b_number);
                        uint64_t midpoint = b_start + b_len/2;
                        memcpy(&seq_region[0], all_sequences[b_sequence]+midpoint-min_len_filter, 2*min_len_filter);
                        seq_region[2*min_len_filter+1] = '\0';
                        fprintf(dna_out, "%s\n", seq_region);
                        fprintf(dna_class, "2\n"); //Its a breakpoint
                        n_breakpoints++;
                    }
                }
                case 1:
                {
                    //If the start of the breakpoint is MIN_LEN bases away from seq begin
                    //and the start + MIN_LEN i
                    //OVERLAPPING
                    while(n_breakpoints < n_blocks && (b_start + 2*min_len_filter) < seq_sizes[b_sequence]){
                        memcpy(&seq_region[0], all_sequences[b_sequence]+b_start, 2*min_len_filter);
                        seq_region[2*min_len_filter] = '\0';
                        fprintf(dna_out, "%s\n", seq_region);
                        fprintf(dna_class, "2\n"); //Its a breakpoint
                        n_breakpoints++;
                        b_start = b_start + 2*min_len_filter;
                    }
                }
                break;
                default:
                {
                    if(b_len >= min_len_filter){ //Only if it is long enough
                        //Copy the region to buffer
                        if(b_number % PRINT_RATE == 0) fprintf(stdout, "[INFO] On breakpoint %"PRIu64"\n", b_number);
                        memcpy(&seq_region[0], all_sequences[b_sequence]+b_start, b_len);
                        seq_region[b_end-b_start+1] = '\0';
                        fprintf(dna_out, "%s\n", seq_region);
                        fprintf(dna_class, "2\n"); //Its a breakpoint
                        n_breakpoints++;
                    }
                }
            }

            
        }
        //printf("%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\n", b_number, b_sequence, b_start, b_end, b_len);
        //getchar();
        
        
    }

    
    


    for(i=0;i<n_files;i++){
        std::free(all_sequences[i]);
    }
    std::free(all_sequences);
    std::free(n_reallocs);
    std::free(seq_region);
    std::free(seq_sizes);

    fclose(dna_class);
    fclose(dna_out);
    fclose(breakpoints);
    fclose(blocks);

    return 0;
}