#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <ctype.h>
#include <pthread.h>
#include <float.h>
#include "structs.h"
#include "evolution.h"
#include "commonFunctions.h"
#include "comparisonFunctions.h"

#define PRINT_RATE 70

void set_base_name(char * s, char * d);

int main(int ac, char **av) {
    if (ac < 4) {
        terror("USE: event_generator <original> <n_sequences> <n_itera>");
    }

    //Iterators
    uint64_t i, j, k, n_files, curr_pos, n_itera;
    

    //Files to read/write results
    FILE * original = fopen64(av[1], "rt"); if(original == NULL) throw "Could not open input sequence file";
    FILE * out_mod;
    n_files = (uint64_t) atoi(av[2]);
    n_itera = (uint64_t) atoi(av[3]);
    
    //Vector to tell for sequence reallocs
    uint64_t * n_reallocs = (uint64_t *) std::calloc(n_files, sizeof(uint64_t));
    if(n_reallocs == NULL) throw "Could not allocate realloc count vector";

    //Char to hold all sequences
    a_sequence * all_sequences = (a_sequence *) std::calloc(n_files, sizeofASequence());
    if(all_sequences == NULL) throw "Could not allocate sequences pointer";
    for(i=0;i<n_files;i++){
        all_sequences[i].s = (char *) std::malloc(SEQ_REALLOC*sizeof(char));
        if(all_sequences[i].s == NULL) throw "Could not allocate initial sequence array";
        n_reallocs[i] = 1;
    }
    
    //Notice: position 0 holds the master sequence

    //Sizes for each sequence
    uint64_t * seq_sizes = (uint64_t *) std::calloc(n_files, sizeof(uint64_t));
    if(seq_sizes == NULL) throw "Could not allocate sequence sizes";


    //Read using buffered fgetc
    uint64_t idx = 0, r = 0;
    char * temp_seq_buffer = NULL;
    if ((temp_seq_buffer = (char *) std::calloc(READBUF, sizeof(char))) == NULL) {
        throw "Could not allocate memory for read buffer";
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
                    all_sequences[0].s[curr_pos++] = c;
                    if(curr_pos >= SEQ_REALLOC*n_reallocs[0]){
                        n_reallocs[0]++;
                        all_sequences[0].s = (char *) std::realloc(all_sequences[0].s, n_reallocs[0]*SEQ_REALLOC);
                        if(all_sequences[0].s == NULL) terror("Could not realloc sequence");
                    }
                }
            }
            curr_pos++; //one for the *
        }else{
            c = buffered_fgetc(temp_seq_buffer, &idx, &r, original);
        }
    }
    //Initialize all sequences equally
    for(i=0;i<n_files;i++) seq_sizes[i] = curr_pos;
    for(i=1;i<n_files;i++) memcpy(&all_sequences[i].s[0], &all_sequences[0].s[0], curr_pos*sizeof(char));
    
    
    //The original sequence is loaded
    
    //Create evolution processes
    mutation ** mutation_proc = (mutation **) std::malloc((n_files-1) * sizeof(mutation *));
    if(mutation_proc == NULL) throw "Could not allocate mutation processes";

    //Attach processes
    uint64_t seed;
    for(i=0;i<n_files-1;i++){
        for(j=0;j<=strlen(av[1])/(n_files-1);j++){
            //Dont even initialize the seed
            seed += (uint64_t) av[1][j];
        }
        mutation_proc[i] = new mutation(((long double)1/seq_sizes[i+1]), &all_sequences[i+1], &seq_sizes[i+1], seed);
    }

    //Apply evolution iteratively
    fprintf(stdout, "[INFO] Running on %"PRIu64" iterations and %"PRIu64" sequences.\n", n_itera, n_files);
    for(i=0;i<n_itera;i++){
        for(j=0;j<n_files-1;j++){
            mutation_proc[0]->step();
            mutation_proc[1]->step();
            mutation_proc[2]->step();
        }
    }

    //Find new seq maximum len
    uint64_t m_len = 0;
    for(i=0;i<n_files;i++) if(seq_sizes[i] > m_len) m_len = seq_sizes[i];

    //Compare seqs
    
    i = 0;
    while(i<m_len){
        if(i % PRINT_RATE == 0){            
            for(k=0;k<n_files;k++){

                fprintf(stdout, "[%"PRIu64"]\t", k);
                j = i;
                while(j<i+PRINT_RATE){
                    if(j < seq_sizes[k]) fprintf(stdout, "%c", all_sequences[k].s[j]);
                    j++;
                }  
                fprintf(stdout, "\n");
            }
        }
        i++;
    }
    

    //Write everything to disk
    char _path[READLINE];
    for(i=1;i<n_files;i++){
        set_base_name(av[1], _path);
        sprintf(_path, "%s%"PRIu64".fasta", _path, i);
        out_mod = fopen64(_path, "wt");
        if(out_mod == NULL) throw "Could not open output files";

        fprintf(out_mod, ">%s\n", _path);
        all_sequences[i].s[seq_sizes[i]] = '\0';
        k = 0;
        for(j=0;j<seq_sizes[i];j++){
            if(j % PRINT_RATE == 0){
                fprintf(out_mod, "%*.*s", 0, PRINT_RATE, &all_sequences[i].s[k]);
                fprintf(out_mod, "\n");
                k+=PRINT_RATE;
            }
        }
        

        fclose(out_mod);
        
    }

    //Free everything

    for(i=0;i<n_files;i++) std::free(all_sequences[i].s);
    for(i=0;i<n_files-1;i++) delete mutation_proc[i];
    std::free(seq_sizes);
    std::free(temp_seq_buffer);
    std::free(n_reallocs);
    std::free(mutation_proc);
    std::free(all_sequences);
    fclose(original);

    return 0;
}

void set_base_name(char * s, char * d){
    uint64_t i;
    memcpy(&d[0], &s[0], strlen(s));
    for(i=strlen(s);i>0;i--){
        if(s[i] == '.') break;
    }
    d[i] = '\0';
}