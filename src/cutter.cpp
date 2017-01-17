#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <pthread.h>
#include <float.h>
#include "structs.h"
#include "commonFunctions.h"
#include "comparisonFunctions.h"

#define SEQ_REALLOC 10000000

char ** read_all_vs_all_files(char * list_of_files, uint64_t * n_files, uint64_t * t_alloc){

    FILE * lf = fopen64(list_of_files, "rt");
    if(lf == NULL) terror("Could not open list of genomic files");

    *t_alloc = 1;
    *n_files = 0;
    uint64_t i;

    char ** all_sequences = (char **) malloc (INIT_SEQS*sizeof(char *));
    for(i=0;i<INIT_SEQS;i++){
        all_sequences[i] = (char *) malloc(READLINE*sizeof(char));
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
            all_sequences = (char **) realloc(all_sequences, (*t_alloc)*INIT_SEQS);
            if(all_sequences == NULL) terror("Could not re-allocate paths to files");
        }
    }

    fclose(lf);
    return all_sequences;
}


int main(int ac, char **av) {
    if (ac < 5) {
        terror("USE: cutter <path_list> <out_dna> <out_class> <file_blocks_breakpoints>");
    }

    //Load files to compare
    uint64_t i, n_files, t_alloc, j, curr_comp = 0, k;

    char ** paths_to_files = read_all_vs_all_files(av[1], &n_files, &t_alloc);
    if(n_files < 2) terror("At least two files need to be loaded");

    //Files to write results
    FILE * dna_out = fopen64(av[2], "wt"); if(dna_out == NULL) terror("Could not open output dna file");
    FILE * dna_class = fopen64(av[3], "wt"); if(dna_class == NULL) terror("Could not open class output file");
    char blocks_bps[READLINE]; blocks_bps[0]='\0';
    strcpy(blocks_bps, av[4]);
    strcat(blocks_bps, ".blocks");
    FILE * blocks = fopen64(blocks_bps, "wt"); if(blocks == NULL) terror("Could not open input blocks file");
    blocks_bps[0]='\0';
    strcpy(blocks_bps, av[4]);
    strcat(blocks_bps, ".breakpoints");
    FILE * breakpoints = fopen64(blocks_bps, "wt"); if(breakpoints == NULL) terror("Could not open breakpoints file");

    FILE * current;

    //Char to hold all sequences
    char ** all_sequences = (char **) std::calloc(n_files, sizeof(char *));

    
    //Read sequences and load into array
    for(i=0;i<n_files;i++){
        current = fopen64(path_to_files[i], "rt");
        all_sequences[i] = (char *) std::calloc(SEQ_REALLOC, sizeof(char));
        if(all_sequences[i] == NULL) terror("Could not allocate genome sequence");
        if(current == NULL) terror("Could not open fasta file");

        curr_pos = 0;
        char c = fgetc(current);
        while(!feof(current)){

            if(c == '>'){
                while(c != '\n') c = fgetc(current); //Skip id

                while(c != '>' && !feof(current)){ //Until next id
                    c = fgetc(current);
                    if(c != toupper(c)){
                        all_sequences[i][curr_pos++] = c;
                    }
                }
                curr_pos++; //one for the *
            }else{
                c = fgetc(current);
            }
            
        }

        fclose(current);
    }

    for(i=0;i<t_alloc;i++){
        free(paths_to_files[i]);
    }
    free(paths_to_files);


    //At this point all sequences are loaded
    uint64_t b_number, b_start, b_end, b_order, b_sequence, b_len;

    while(!feof(blocks)){
        fscanf(blocks, "%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64, &b_number, &b_sequence, &b_order, &b_start, &b_end, &b_len);
    }


    for(i=0;i<n_files;i++){
        free(all_sequences[i]);
    }
    free(all_sequences);

    fclose(dna_class);
    fclose(dna_out);
    fclose(breakpoints);
    fclose(blocks);

    return 0;
}