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



int main(int ac, char **av) {
    if (ac < 3) {
        terror("USE: event_generator <original> <n_sequences>");
    }

    //Iterators
    uint64_t i, n_files, curr_pos;
    

    //Files to read/write results
    FILE * original = fopen64(av[1], "rt"); if(original == NULL) terror("Could not open input sequence file");
    FILE * out_mod;
    n_files = (uint64_t) atoi(av[2]);
    
    //Vector to tell for sequence reallocs
    uint64_t * n_reallocs = (uint64_t *) std::calloc(n_files, sizeof(uint64_t));
    if(n_reallocs == NULL) terror("Could not allocate realloc count vector");

    //Char to hold all sequences
    char ** all_sequences = (char **) std::calloc(n_files, sizeof(char *));
    if(all_sequences == NULL) terror("Could not allocate sequences pointer");
    all_sequences[0] = (char *) std::malloc(SEQ_REALLOC*sizeof(char));
    if(all_sequences[0] == NULL) terror("Could not allocate initial sequence array");
    n_reallocs[0] = 1;

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
    curr_pos = 0;

    //Read original sequence
    c = buffered_fgetc(temp_seq_buffer, &idx, &r, original);
    while((!feof(original) || (feof(original) && idx < r))){
        if(c == '>'){
            while(c != '\n') c = buffered_fgetc(temp_seq_buffer, &idx, &r, original); //Skip id

            while(c != '>' && (!feof(original) || (feof(original) && idx < r))){ //Until next id
                c = buffered_fgetc(temp_seq_buffer, &idx, &r, original);
                c = toupper(c);
                if(c >= 'A' && c <= 'Z'){
                    all_sequences[0][curr_pos++] = c;
                    if(curr_pos >= SEQ_REALLOC*n_reallocs[0]){
                        n_reallocs[0]++;
                        all_sequences[0] = (char *) std::realloc(all_sequences[0], n_reallocs[0]*SEQ_REALLOC);
                        if(all_sequences[0] == NULL) terror("Could not realloc sequence");
                    }
                }
            }
            curr_pos++; //one for the *
        }else{
            c = buffered_fgetc(temp_seq_buffer, &idx, &r, original);
        }
    }
    seq_sizes[0] = curr_pos;
    
    for(i=0;i<seq_sizes[0];i++){
        printf("%c", all_sequences[0][i]);
    }
    printf("\n");

    //The original sequence is loaded
    







    //Free everything

    for(i=0;i<n_files;i++){
        if(all_sequences[i] != NULL) std::free(all_sequences[i]);
    }
    std::free(seq_sizes);
    std::free(temp_seq_buffer);
    std::free(n_reallocs);
    fclose(original);
    //fclose(out_mod);

    return 0;
}