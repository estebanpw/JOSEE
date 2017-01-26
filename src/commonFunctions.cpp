#include <stdio.h>
#include <cstdlib>
#include <sys/time.h>
#include <inttypes.h>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "structs.h"
#include "comparisonFunctions.h"


void terror(const char *s) {
    printf("ERR**** %s ****\n", s);
    exit(-1);
}

char buffered_fgetc(char *buffer, uint64_t *pos, uint64_t *read, FILE *f) {
    if (*pos >= READBUF) {
        *pos = 0;
        memset(buffer, 0, READBUF);
        *read = fread(buffer, 1, READBUF, f);
    }
    *pos = *pos + 1;
    return buffer[*pos-1];
}

void load_fragments_local(FILE * fragsfile, uint64_t * n_frags, struct FragFile ** loaded_frags){
    
    struct FragFile temp_frag;
    uint64_t unused_len, total_frags;

    //Compute number of fragments in file
    fseeko(fragsfile, 0L, SEEK_END);
    total_frags = ftello(fragsfile) - 2*(sizeof(uint64_t)); //Remove the headers
    fseeko(fragsfile, 0L, SEEK_SET);

    //Divide by size of frag to get the number of fragments
    //Plus one because it might have padding, thus rounding up to bottom and missing 1 struct
    total_frags = 1 + total_frags/sizeofFragment(); 
    
    //Allocate memory for all frags
    struct FragFile * temp_frags_array = (struct FragFile *) std::malloc(total_frags * sizeofFragment());
    if(temp_frags_array == NULL) terror("Could not allocate heap for all fragments");

    fprintf(stdout, "[INFO] There are %"PRIu64" fragments to be loaded, requiring %"PRIu64" Megabyte(s) of RAM\n", total_frags, (total_frags*sizeof(struct FragFile))/(1024*1024));

    //Skip headers
    readSequenceLength(&unused_len, fragsfile);
    readSequenceLength(&unused_len, fragsfile);

    //To keep track of current frag
    *n_frags = 0;


    while(!feof(fragsfile)){
        readFragment(&temp_frag, fragsfile);
        
        //Transform coordinates to local
        //Actually, it is not required anymore.
        
        /*if(prevx != temp_frag.seqX || prevy != temp_frag.seqY){
            printf("Frag: (%"PRIu64", %"PRIu64") coord: (%"PRIu64", %"PRIu64") diag:: %"PRId64"\n", temp_frag.seqX, temp_frag.seqY, temp_frag.xStart, temp_frag.yStart, temp_frag.diag);
            getchar();
        }
        prevx = temp_frag.seqX; prevy = temp_frag.seqY;*/

        
        //printf("SeqX SeqY: (%"PRIu64", %"PRIu64")\n", temp_frag.seqX, temp_frag.seqY);
        //printf("Frags. (%"PRIu64", %"PRIu64")\n", sequences[temp_frag.seqX].acum, sequences[temp_frag.seqY].acum);        

        temp_frag.xStart = temp_frag.xStart;
        temp_frag.xEnd = temp_frag.xEnd;
        

        temp_frag.yStart = temp_frag.yStart;
        temp_frag.yEnd = temp_frag.yEnd;


        //temp_frag.xStart = temp_frag.xStart - sequences[temp_frag.seqX].acum;
        //temp_frag.xEnd = temp_frag.xEnd - sequences[temp_frag.seqX].acum; 

        //temp_frag.yStart = temp_frag.yStart - sequences[temp_frag.seqY].acum;
        //temp_frag.yEnd = temp_frag.yEnd - sequences[temp_frag.seqY].acum;


        //Copy temp fragment into array
        memcpy(&temp_frags_array[*n_frags], &temp_frag, sizeofFragment());
        
        *n_frags = *n_frags + 1;

        if(*n_frags > total_frags){ terror("Something went wrong. More fragments than expected");}
    }

    //Copy pointer of array heap
    *loaded_frags = temp_frags_array;
}

void get_coverage_from_genome_grid(unsigned char ** maptable, sequence_manager * seq_manager, uint64_t n_seqs, uint64_t min_len_without_breaks){
    uint64_t i, j;

    uint64_t sum_bases;
    uint64_t current;

    for(i=0;i<n_seqs;i++){
        current = 0;
        sum_bases = 0;
        for(j=0;j<seq_manager->get_sequence_by_label(i)->len;j++){
            if(maptable[i][j] != COVERFRAG){
                if(current >= min_len_without_breaks){
                    sum_bases += current;
                }
                current = 0;
            }else{
                current++;
            }
        }
        seq_manager->get_sequence_by_label(i)->coverage = (100*sum_bases)/seq_manager->get_sequence_by_label(i)->len;
    }
}


void write_maptable_to_disk(unsigned char ** maptable, sequence_manager * seq_manager, const char * out_file_path){
    uint64_t i, j;
    FILE * out_file;
    
    char output_file[512];

    for(i=0; i<seq_manager->get_number_of_sequences(); i++){

        output_file[0]='\0';
        sprintf(output_file, "%s_%d", out_file_path, (int) i);
        out_file = fopen64(output_file, "wt");
        if(out_file == NULL) terror("Could not open output debug file");

        for(j=0; j<seq_manager->get_sequence_by_label(i)->len; j++){
            fprintf(out_file, "%u", maptable[i][j]);
            if(j != 0 && j % 50 == 0) fprintf(out_file, "\n");
        }
        fprintf(out_file, "\n");

        fclose(out_file);
    }

}

void print_maptable_portion(unsigned char ** maptable, uint64_t from, uint64_t to, uint64_t rate, uint64_t seq){
    uint64_t i, acu = 0;
    static uint64_t CALL = 0;
    fprintf(stdout, "===========(%"PRIu64" -> %"PRIu64")@%"PRIu64"\n", from, to, seq);
    for(i=from;i<to;i++){
        fprintf(stdout, "%u", maptable[seq][i]);
        if(acu != 0 && acu % rate == 0) fprintf(stdout, "\n");
        acu++;
    }
    fprintf(stdout, "\n===========CALL: %"PRIu64"\n", CALL);
    CALL++;
}

void traverse_synteny_list(Synteny_list * sbl){
    Synteny_list * ptr_sbl = sbl;
    Synteny_block * ptr_sb;
    while(ptr_sbl != NULL){
        fprintf(stdout, "SBL:\n");
        ptr_sb = ptr_sbl->sb;
        while(ptr_sb != NULL){
            fprintf(stdout, "\t");printBlock(ptr_sb->b);
            ptr_sb = ptr_sb->next;
        }
        ptr_sbl = ptr_sbl->next;
        getchar();
    }
}

Synteny_list * find_synteny_block_from_block(Synteny_list * sbl, Block * b){
    Synteny_list * ptr_sbl = sbl;
    Synteny_block * ptr_sb;
    while(ptr_sbl != NULL){
        ptr_sb = ptr_sbl->sb;
        while(ptr_sb != NULL){
            if(isBlockEqualTo(b, ptr_sb->b) == 1){
                return ptr_sbl;
            }
            ptr_sb = ptr_sb->next;
        }
        ptr_sbl = ptr_sbl->next;
    }
    return NULL;
}

void find_fragments_from_maptable(unsigned char ** maptable, uint64_t start, uint64_t end, uint64_t seq, struct FragFile * frags, uint64_t n_frags){
    uint64_t i = 0;
    struct FragFile target;
    target.xStart = start;
    target.xEnd = end;
    target.seqX = seq;

    while(i<n_frags){
        if(isFragmentEqualTo(&target, &frags[i]) == 1){
            printFragment(&frags[i]);
            getchar();
        }
        i++;
    }
}


