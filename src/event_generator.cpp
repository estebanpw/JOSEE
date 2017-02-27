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
#define DUP_SIZE 100
#define INS_SIZE 100
#define DEL_SIZE 100
#define INV_SIZE 100

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

    //Probabilities multiplier
    srand(time(NULL));
    uint64_t seed = seed + rand();
    std::default_random_engine generator(seed);
    std::uniform_real_distribution<long double> d_r_unif(0.90, 1.25); //Empirically chosen
    uint64_t event_size = seq_sizes[0]/10;

    //Create evolution processes
    dna_mutation ** mutation_proc = (dna_mutation **) std::malloc((n_files-1) * sizeof(dna_mutation *));
    if(mutation_proc == NULL) throw "Could not allocate mutation processes";

    dna_duplication ** duplication_proc = (dna_duplication **) std::malloc((n_files-1) * sizeof(dna_duplication *));
    if(duplication_proc == NULL) throw "Could not allocate duplication processes";

    dna_insertion ** insertion_proc = (dna_insertion **) std::malloc((n_files-1) * sizeof(dna_insertion *));
    if(insertion_proc == NULL) throw "Could not allocate insertion processes";

    dna_deletion ** deletion_proc = (dna_deletion **) std::malloc((n_files-1) * sizeof(dna_deletion *));
    if(deletion_proc == NULL) throw "Could not allocate deletion processes";

    dna_inversion ** inversion_proc = (dna_inversion **) std::malloc((n_files-1) * sizeof(dna_inversion *));
    if(inversion_proc == NULL) throw "Could not allocate inversion processes";

    dna_transposition ** transposition_proc = (dna_transposition **) std::malloc((n_files-1) * sizeof(dna_transposition *));
    if(transposition_proc == NULL) throw "Could not allocate transposition processes";

    //Attach processes
    
    long double p_mut = ((long double) 10);
    long double p_dup = ((long double) n_itera*d_r_unif(generator));
    long double p_ins = ((long double) n_itera*d_r_unif(generator));
    long double p_del = ((long double) n_itera*d_r_unif(generator));
    long double p_inv = ((long double) n_itera*d_r_unif(generator));
    long double p_tra = ((long double) n_itera*d_r_unif(generator));

    long double s_size_init = (long double)1/seq_sizes[0];

    fprintf(stdout, "[INFO] Using probabilities %Le, %Le, %Le, %Le, %Le, %Le\n", p_mut, p_dup, p_ins, p_del, p_inv, p_tra);
    for(i=0;i<n_files-1;i++){
        for(j=0;j<=strlen(av[1])/(n_files-1);j++){
            //Dont even initialize the seed
            seed += (uint64_t) av[1][j];
        }
        printf("[INFO] Seed: %"PRIu64"\n", seed);
        mutation_proc[i] = new dna_mutation(s_size_init/p_mut, &all_sequences[i+1], &seq_sizes[i+1], seed);
        duplication_proc[i] = new dna_duplication(s_size_init/p_dup, &all_sequences[i+1], &seq_sizes[i+1], event_size, seed);
        insertion_proc[i] = new dna_insertion(s_size_init/p_ins, &all_sequences[i+1], &seq_sizes[i+1], event_size, seed);
        deletion_proc[i] = new dna_deletion(s_size_init/p_del, &all_sequences[i+1], &seq_sizes[i+1], event_size, seed);
        inversion_proc[i] = new dna_inversion(s_size_init/p_inv, &all_sequences[i+1], &seq_sizes[i+1], event_size, seed);
        transposition_proc[i] = new dna_transposition(s_size_init/p_tra, &all_sequences[i+1], &seq_sizes[i+1], event_size, seed);
    }

    

    //Apply evolution iteratively
    
    uint64_t step_size = n_itera/10;
    uint64_t acum = 0;
    for(i=0;i<n_itera;i++){
        //fprintf(stdout, "[INFO] Using probabilities %Le, %Le, %Le, %Le, %Le, %Le\n", s_size_init/p_mut, s_size_init/p_dup, s_size_init/p_ins, s_size_init/p_del, s_size_init/p_inv, s_size_init/p_tra);
        if(i % step_size == 0){ acum+= 10; fprintf(stdout, "..%"PRIu64"%%", acum); fflush(stdout); }
        for(j=0;j<n_files-1;j++){
        
            
            mutation_proc[j]->step();
            duplication_proc[j]->step();
            insertion_proc[j]->step();
            deletion_proc[j]->step();
            inversion_proc[j]->step();
            transposition_proc[j]->step();
            
            // Update probabilities (because of length change)
            s_size_init = (long double)1/seq_sizes[j+1];
            mutation_proc[j]->set_p(s_size_init/p_mut);
            duplication_proc[j]->set_p(s_size_init/p_dup);
            insertion_proc[j]->set_p(s_size_init/p_ins);
            deletion_proc[j]->set_p(s_size_init/p_del);
            inversion_proc[j]->set_p(s_size_init/p_inv);
            transposition_proc[j]->set_p(s_size_init/p_tra);
            
        }
    }
    fprintf(stdout, "\n");

    //Find new seq maximum len
    uint64_t m_len = 0;
    for(i=0;i<n_files;i++) if(seq_sizes[i] > m_len) m_len = seq_sizes[i];

    //Compare seqs
    /*
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
    */

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
    std::free(duplication_proc);
    std::free(insertion_proc);
    std::free(deletion_proc);
    std::free(inversion_proc);
    std::free(transposition_proc);
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