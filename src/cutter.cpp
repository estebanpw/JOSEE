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

//Function to read all files from a text list
char ** read_all_vs_all_files(char * list_of_files, uint64_t * n_files, uint64_t * t_alloc){

    FILE * lf = fopen64(list_of_files, "rt");
    if(lf == NULL) terror("Could not open list of genomic files");

    *t_alloc = 1;
    *n_files = 0;
    uint64_t i;

    char ** all_sequences = (char **) std::malloc (INIT_SEQS*sizeof(char *));
    for(i=0;i<INIT_SEQS;i++){
        all_sequences[i] = (char *) std::malloc(READLINE*sizeof(char));
        if(all_sequences[i] == NULL) terror("Could not allocate paths to files");
    }



    while(!feof(lf)){
        if(fgets(all_sequences[*n_files], READLINE, lf) > 0){
            if(all_sequences[*n_files][0] != '\0' && all_sequences[*n_files][0] != '\n'){
                all_sequences[*n_files][strlen(all_sequences[*n_files])-1] = '\0';
                (*n_files)++;
            }
        }
        if(*n_files == INIT_SEQS*(*t_alloc)){
            (*t_alloc)++;
            all_sequences = (char **) std::realloc(all_sequences, (*t_alloc)*INIT_SEQS);
            if(all_sequences == NULL) terror("Could not re-allocate paths to files");
        }
    }

    fclose(lf);
    return all_sequences;
}


int main(int ac, char **av) {
    if (ac < 6) {
        terror("USE: cutter <path_list> <out_dna> <out_class> <file_blocks_breakpoints> <min_len_filter>");
    }

    //Iterators
    uint64_t i, n_files, t_alloc, curr_pos;

    //Load files to compare
    char ** paths_to_files = read_all_vs_all_files(av[1], &n_files, &t_alloc);
    if(n_files < 2) terror("At least two files need to be loaded");
    for(uint64_t k=0;k<n_files;k++){ fprintf(stdout, "File %"PRIu64": %s\n", k, paths_to_files[k]);}

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

    //The file that will open all sequence files
    FILE * current;

    //Vector to tell for sequence reallocs
    uint64_t * n_reallocs = (uint64_t *) std::calloc(n_files, sizeof(uint64_t));
    if(n_reallocs == NULL) terror("Could not allocate realloc count vector");

    //Char to hold all sequences
    char ** all_sequences = (char **) std::calloc(n_files, sizeof(char *));
    if(all_sequences == NULL) terror("Could not allocate sequences pointer");

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

        fclose(current);
    }


    std::free(temp_seq_buffer);

    for(i=0;i<t_alloc;i++){
        std::free(paths_to_files[i]);
    }
    std::free(paths_to_files);


    //At this point all sequences are loaded
    uint64_t b_number, b_start, b_end, b_order, b_sequence, b_len;
    char header[READLINE];
    char * seq_region = (char *) std::malloc(SEQ_REALLOC*sizeof(char));
    if(seq_region == NULL) terror("Could not allocate sequence region to output");
    uint64_t seq_region_reallocs = 1;

    //Skip header
    if(NULL == fgets(header, READLINE, blocks)) terror("No header was read or empty blocks file");
    while(!feof(blocks)){
        //Read a line i.e. a block
        if(6 == fscanf(blocks, "%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"", &b_number, &b_sequence, &b_order, &b_start, &b_end, &b_len)){
            //Check that the region fits in vector
            if((b_end - b_start +1) >= seq_region_reallocs*SEQ_REALLOC){
                seq_region_reallocs++;
                seq_region = (char *) std::realloc(seq_region, seq_region_reallocs*SEQ_REALLOC);
                if(seq_region == NULL) terror("Could not realloc region sequence");
            }
            if(b_end - b_start >= min_len_filter){//Only if it is long enough
                //Copy the region to buffer
                if(b_number % PRINT_RATE == 0) fprintf(stdout, "[INFO] On block %"PRIu64"\n", b_number);
                memcpy(&seq_region[0], all_sequences[b_sequence]+b_start, b_end-b_start);
                seq_region[b_end-b_start+1] = '\0';
                fprintf(dna_out, "%s\n", seq_region);
                fprintf(dna_class, "1\n"); //Its a block
            }   
        }
    }
    //Repeat for breakpoints
    //Skip header
    if(NULL == fgets(header, READLINE, breakpoints)) terror("No header was read or empty breakpoints file");
    while(!feof(breakpoints)){
        //Read a line i.e. a block
        
        if(5 == fscanf(breakpoints, "%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"", &b_number, &b_sequence, &b_start, &b_end, &b_len)){
            //Check that the region fits in vector
            if((b_end - b_start +1) >= seq_region_reallocs*SEQ_REALLOC){
                seq_region_reallocs++;
                seq_region = (char *) std::realloc(seq_region, seq_region_reallocs*SEQ_REALLOC);
                if(seq_region == NULL) terror("Could not realloc region sequence");
            }
            if(b_end - b_start >= min_len_filter){ //Only if it is long enough
                //Copy the region to buffer
                if(b_number % PRINT_RATE == 0) fprintf(stdout, "[INFO] On breakpoint %"PRIu64"\n", b_number);
                memcpy(&seq_region[0], all_sequences[b_sequence]+b_start, b_end-b_start);
                seq_region[b_end-b_start+1] = '\0';
                fprintf(dna_out, "%s\n", seq_region);
                fprintf(dna_class, "2\n"); //Its a breakpoint
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

    fclose(dna_class);
    fclose(dna_out);
    fclose(breakpoints);
    fclose(blocks);

    return 0;
}