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

int compare_ranges(Annotation * a, Annotation * b){
    if(a->start <= b->end && b->start <= a->end) return 0; //Overlap
    if(a->start < b->start) return -1;
    return 1;
}

int compare_two_annotations(Annotation * a, Annotation * b){
    if(a->start == b->start) return 0;
    if(a->start < b->start) return -1;
    return 1;
}

Annotation * binary_search_annotations(uint64_t start, uint64_t end, Annotation * anot, uint64_t n_annots){

   uint64_t first = 0;
   uint64_t last = n_annots - 1;
   uint64_t middle = (first+last)/2;
   int compare;
   Annotation aux;
   aux.start = start;
   aux.end = end;
   
   while (first <= last) {

       compare = compare_ranges(&aux, &anot[middle]);
       if (compare == 0) return &anot[middle];
       if (compare<0) last = middle - 1;
       else first = middle + 1;
       middle = (first + last)/2;
   }
   return NULL; // Not found
}

void quick_sort_annotations(Annotation * array, uint64_t x, uint64_t y) {

    Annotation pivot, aux;
    uint64_t x1, y1;

    memcpy(&pivot, &array[(x+y)/2], sizeofAnnotation());
    x1 = x;
    y1 = y;

    do{

        while (compare_two_annotations(&pivot, &array[x1]) > 0) x1++;
        while (compare_two_annotations(&pivot, &array[y1]) < 0) y1--;
        if (x1 < y1) { 
            memcpy(&aux, &array[x1], sizeofAnnotation());
            memcpy(&array[x1], &array[y1], sizeofAnnotation());
            memcpy(&array[y1], &aux, sizeofAnnotation());
            x1++;
            y1--;
        }
        else if (x1 == y1) x1++;
    } while (x1 <= y1);

    if (x < y1) quick_sort_annotations(array, x, y1);
    if (x1 < y) quick_sort_annotations(array, x1, y);
}

uint64_t quick_pow4(uint32_t n){
    static uint64_t pow4[33]={1L, 4L, 16L, 64L, 256L, 1024L, 4096L, 16384L, 65536L, 
    262144L, 1048576L, 4194304L, 16777216L, 67108864L, 268435456L, 1073741824L, 4294967296L, 
    17179869184L, 68719476736L, 274877906944L, 1099511627776L, 4398046511104L, 17592186044416L, 
    70368744177664L, 281474976710656L, 1125899906842624L, 4503599627370496L, 18014398509481984L, 
    72057594037927936L, 288230376151711744L, 1152921504606846976L, 4611686018427387904L};
    return pow4[n];
}

uint64_t quick_pow4byLetter(uint32_t n, const char c){
    static uint64_t pow4_G[33]={2*1L, 2*4L, 2*16L, 2*64L, 2*256L, 2*1024L, 2*4096L, 2*16384L, 2*65536L, 
    (uint64_t)2*262144L, (uint64_t)2*1048576L,(uint64_t)2*4194304L, (uint64_t)2*16777216L, (uint64_t)2*67108864L, (uint64_t)2*268435456L, (uint64_t)2*1073741824L, (uint64_t)2*4294967296L, 
    (uint64_t)2*17179869184L, (uint64_t)2*68719476736L, (uint64_t)2*274877906944L, (uint64_t)2*1099511627776L, (uint64_t)2*4398046511104L, (uint64_t)2*17592186044416L, 
    (uint64_t)2*70368744177664L, (uint64_t)2*281474976710656L, (uint64_t)2*1125899906842624L, (uint64_t)2*4503599627370496L, (uint64_t)2*18014398509481984L, 
    (uint64_t)2*72057594037927936L, (uint64_t) 2*288230376151711744L, (uint64_t) 2*1152921504606846976L, (uint64_t) 2*4611686018427387904L};
    
    static uint64_t pow4_T[33]={3*1L, 3*4L, 3*16L, 3*64L, 3*256L, 3*1024L, 3*4096L, 3*16384L, 3*65536L, 
    (uint64_t)3*262144L, (uint64_t) 3*1048576L, (uint64_t)3*4194304L, (uint64_t)3*16777216L, (uint64_t)3*67108864L, (uint64_t)3*268435456L, (uint64_t)3*1073741824L, (uint64_t)3*4294967296L, 
    (uint64_t)3*17179869184L, (uint64_t)3*68719476736L, (uint64_t)3*274877906944L, (uint64_t)3*1099511627776L, (uint64_t)3*4398046511104L, (uint64_t)3*17592186044416L, 
    (uint64_t)3*70368744177664L, (uint64_t)3*281474976710656L, (uint64_t)3*1125899906842624L, (uint64_t)3*4503599627370496L, (uint64_t)3*18014398509481984L, 
    (uint64_t)3*72057594037927936L, (uint64_t) 3*288230376151711744L, (uint64_t) 3*1152921504606846976L, (uint64_t) 3*4611686018427387904L};
    
    if(c == 'A') return 0;
    if(c == 'C') return quick_pow4(n);
    if(c == 'G') return pow4_G[n];
    if(c == 'T') return pow4_T[n];
    return 0;
}

uint64_t hashOfWord(const char * word, uint32_t k){
    
    uint64_t value = 0, jIdx;
    for(jIdx=0;jIdx<k;jIdx++){
        value += quick_pow4byLetter(k-(jIdx+1), word[jIdx]);
    }
    return value;
    
}

inline int64_t compare_letters(char a, char b){
    if(a != 'N') return (a == b) ? POINT : -POINT;
    return -POINT;
}

int overlapped_words(uint64_t xstart, uint64_t xend, uint64_t ystart, uint64_t yend){
    if(xstart <= yend && ystart <= xend) return 0; //Overlap
    if(xstart < ystart) return -1;
    return 1;
}

void alignment_from_hit(sequence_manager * seq_man, Word * a, Word * b, Quickfrag * qf, uint64_t kmer_size){

    //@Important: a->pos and b->pos should be ending of the hit + 1

    /*
    printf("From hits\n");
    seq_man->print_sequence_region(a->genome->id, a->pos - kmer_size, a->pos);
    seq_man->print_sequence_region(b->genome->id, b->pos - kmer_size, b->pos);
    */

    int64_t curr_pos_a = (int64_t) a->pos;
    int64_t curr_pos_b = (int64_t) b->pos;
    int64_t final_end_a = (int64_t) a->pos - 1, final_start_a = final_end_a - kmer_size + 1, final_start_b = curr_pos_b - kmer_size;
    int64_t score_right = (int64_t) kmer_size * POINT;
    int64_t score_left = score_right;
    int64_t high_left = score_left, high_right = score_right;
    int64_t start_block_a = (int64_t) a->b->start;
    int64_t end_block_a = (int64_t) a->b->end;
    int64_t start_block_b = (int64_t) b->b->start;
    int64_t end_block_b = (int64_t) b->b->end;
    qf->t_len = kmer_size;
    uint64_t idents = kmer_size;
    uint64_t final_idents = 0;

    int keep_going = 1;

    //Forward search
    while(keep_going == 1){
        
        
        if(score_right > 0 && curr_pos_a < end_block_a && curr_pos_b < end_block_b){
            if(curr_pos_a  > end_block_a ||  curr_pos_b > end_block_b) break;
            if(compare_letters(a->b->genome->seq[curr_pos_a], b->b->genome->seq[curr_pos_b]) == POINT){ score_right+=POINT; idents++; }else{ score_right-=POINT;}
            if(high_right <= score_right){
                final_end_a = curr_pos_a;
                high_right = score_right;
                final_idents = idents;
            }
            curr_pos_a++;
            curr_pos_b++;
        }else{
            keep_going = 0;
        }
    }

    keep_going = 1;
    curr_pos_a = a->pos - kmer_size - 1;
    curr_pos_b = b->pos - kmer_size - 1;

    score_left = high_right;

    //Backward search
    while(keep_going == 1){
        
        if(score_left > 0 && curr_pos_a >= start_block_a && curr_pos_b >= start_block_b){
            if(curr_pos_a < start_block_a || curr_pos_b < start_block_b ) break;
            if(compare_letters(a->b->genome->seq[curr_pos_a], b->b->genome->seq[curr_pos_b]) == POINT){ score_left+=POINT; idents++; }else{ score_left-=POINT;}
            if(high_left <= score_left){
                final_start_a = curr_pos_a;
                final_start_b = curr_pos_b;
                high_left = score_left;
                final_idents = idents;
            }
            curr_pos_a--;
            curr_pos_b--;
        }else{
            keep_going = 0;
        }
    }

    

    qf->t_len = final_end_a - final_start_a + 1;
    qf->sim = (final_idents / (long double) qf->t_len)*100;
    qf->x_start = final_start_a;
    qf->y_start = final_start_b;
    qf->diag = (int64_t)qf->x_start - (int64_t)qf->y_start;

    /*
    printf("Aligned\n");
    seq_man->print_sequence_region(a->genome->id, final_start_a, final_end_a);
    seq_man->print_sequence_region(b->genome->id, final_start_b, final_start_b+qf->t_len);
    getchar();
    */

}



void read_words_from_synteny_block_and_align(sequence_manager * seq_man, Synteny_list * sbl, uint32_t kmer_size, dictionary_hash * dhw, Quickfrag ** qfmat, unsigned char ** qfmat_state){

    //@Important: A single hit between two sequences should be enough to align two blocks

    //Erase what we had for previous alignments
    dhw->clear();
    uint64_t i,j;
    for(i=0;i<seq_man->get_number_of_sequences();i++){
        for(j=0;j<seq_man->get_number_of_sequences();j++) qfmat_state[i][j] = 0;
    }
    Synteny_list * sbl_ptr = sbl;
    Synteny_block * sb_ptr;
    //To keep track of where we are
    uint64_t advanced_steps;

    //Kmer reading
    char curr_kmer[kmer_size];
    uint64_t kmer_index = 0;
    char c;

    //To make things clearer
    Quickfrag * qf;
    Quickfrag aligned_qf;
    int64_t curr_diag;
    Word align_word;

    //And to speed things up
    uint64_t precomputed_sizeofQuickfrag = sizeofQuickfrag();
    Wordbucket * hit;
    Wordbucket ** hit_list;

    sb_ptr = sbl_ptr->sb;

    //For all blocks in the synteny
    while(sb_ptr != NULL){

        //Get next word of current block
        advanced_steps = sb_ptr->b->start;
        while(advanced_steps < sb_ptr->b->end){

            //Get nucleotide
            c = sb_ptr->b->genome->seq[advanced_steps];
            advanced_steps++;

            if(c == 'A' || c == 'C' || c == 'T' || c == 'G'){
                curr_kmer[kmer_index] = c;
                kmer_index++;
            }else{
                kmer_index = 0;
            }

            //Check if we have a kmer big enough
            if(kmer_index == kmer_size){

                //printf("Putting %s at %"PRIu64"\n", curr_kmer, advanced_steps);

                //Insert word in dictionary 
                //printf("%s\n", curr_kmer);
                hit_list = dhw->put_and_hit(curr_kmer, 'f', advanced_steps, sb_ptr->b);

                for(i=0;i<dhw->get_candidates_number();i++){
                    //For every hit candidate
                    hit = hit_list[i];

                    if(hit != NULL){
                        //printf("Got hit and the state is %u. The synteny level is %"PRIu64"\n", (int)qfmat_state[sb_ptr->b->genome->id][hit->w.b->genome->id], sbl->synteny_level);
                        //printf("saux_id: %"PRIu64" and hit genome id: %"PRIu64"\n", sb_ptr->b->genome->id, hit->w.b->genome->id);
                        //getchar();
                        //In the sake of clarity
                        qf = &qfmat[sb_ptr->b->genome->id][hit->w.b->genome->id];

                        //We have a hit we should try an alignment
                        //only if the hit is not overlapping
                        //or it is overlapping but on different diagonal
                        if(qfmat_state[sb_ptr->b->genome->id][hit->w.b->genome->id] == 0){
                            //printf("BINGO!!!!!!!!!!!!!\n");
                            
                            //There is no frag yet, so try first one
                            align_word.pos = advanced_steps;
                            align_word.strand = 'f';
                            align_word.b = sb_ptr->b;
                            alignment_from_hit(seq_man, &align_word, &hit->w, &aligned_qf, kmer_size);
                            printQuickfrag(&aligned_qf);
                            qfmat_state[sb_ptr->b->genome->id][hit->w.b->genome->id] = 1;
                            qfmat_state[hit->w.b->genome->id][sb_ptr->b->genome->id] = 1;
                            memcpy(qf, &aligned_qf, precomputed_sizeofQuickfrag);
                            memcpy(&qfmat[hit->w.b->genome->id][sb_ptr->b->genome->id], &aligned_qf, precomputed_sizeofQuickfrag);

                        }else{
                            //printf("on the else track\n");
                            //getchar();
                            //There is already a frag, check for overlapping and diagonal
                            if(overlapped_words(qf->x_start, qf->x_start+qf->t_len, advanced_steps-kmer_size-1, advanced_steps-1) != 0){
                                //If it is not overlapped 
                                curr_diag = (int64_t) advanced_steps - (int64_t) hit->w.pos;
                                if(curr_diag != qf->diag){
                                    //We can try new alignment
                                    align_word.pos = advanced_steps;
                                    align_word.strand = 'f';
                                    align_word.b->genome = sb_ptr->b->genome;
                                    alignment_from_hit(seq_man, &align_word, &hit->w, &aligned_qf, kmer_size);
                                    printQuickfrag(&aligned_qf);
                                    //printf("HAPPENS AS WELL\n");
                                    //getchar();
                                    //getchar();

                                    //Only copy if new alignment is better 

                                    if(aligned_qf.sim > qf->sim){
                                        memcpy(qf, &aligned_qf, precomputed_sizeofQuickfrag);
                                        memcpy(&qfmat[hit->w.b->genome->id][sb_ptr->b->genome->id], &aligned_qf, precomputed_sizeofQuickfrag);
                                    }
                                }
                            }
                        }
                        
                    }
                }

                

                //Insert reversed word in dictionary

                //Displace 
                memmove(&curr_kmer[0], &curr_kmer[1], kmer_size-1);
                kmer_index--;
            }

        }

        //Advance block
        sb_ptr = sb_ptr->next;
        kmer_index = 0;
    }
    printSyntenyBlock(sbl->sb);
    printQuickFragMatrix(qfmat, qfmat_state, seq_man->get_number_of_sequences());
    getchar();
    
}

