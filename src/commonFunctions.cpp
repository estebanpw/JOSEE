#include <stdio.h>
#include <cstdlib>
#include <sys/time.h>
#include <inttypes.h>
#include <ctype.h>
#include <cfloat>
#include <string.h>
#include <math.h>
#include "structs.h"
#include "comparisonFunctions.h"
#include "alignment_functions.h"


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
        
        //printFragment(&temp_frag); getchar();
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
    //bool forward = true;
    while(ptr_sbl != NULL){
        fprintf(stdout, "SBL:\n");
        ptr_sb = ptr_sbl->sb;
        while(ptr_sb != NULL){
            fprintf(stdout, "\t");printBlock(ptr_sb->b);
            ptr_sb = ptr_sb->next;
        }
        /*
        if(ptr_sbl->next == NULL){
            forward = !forward;
            printf("And changing direction!!!!!!!\n");
        }
        
        if(forward){
            ptr_sbl = ptr_sbl->next;
        }else{
            ptr_sbl = ptr_sbl->prev;
        }
        */
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
    qf->x = a->b->genome;
    qf->y = b->b->genome;

    /*
    printf("Aligned\n");
    seq_man->print_sequence_region(a->genome->id, final_start_a, final_end_a);
    seq_man->print_sequence_region(b->genome->id, final_start_b, final_start_b+qf->t_len);
    getchar();
    */

}


inline char complement(char c){
    
    switch(c){
    case ('A'): return 'T';
    break;
    case ('C'): return 'G';
    break;
    case ('G'): return 'C';
    break;
    case ('T'): return 'A';
    break;
    }
    terror("Unrecognized nucleotide in hit");
    return 0; //I know it will never get here, but the compiler des not.

}


void alignment_from_hit_reverse(sequence_manager * seq_man, Word * a, Word * b, Quickfrag * qf, uint64_t kmer_size){

    //@Important: a->pos and b->pos should be ending of the hit + 1
    //@Important: a will be reversed
    //@Important: Everything is the same except how the access is done

    /*
    printf("From hits\n");
    seq_man->print_sequence_region(a->genome->id, a->pos - kmer_size, a->pos);
    seq_man->print_sequence_region(b->genome->id, b->pos - kmer_size, b->pos);
    */

    int64_t curr_pos_a = (int64_t) a->pos - kmer_size - 1;
    int64_t curr_pos_b = (int64_t) b->pos;
    int64_t final_end_a = (int64_t) curr_pos_a, final_start_a = a->pos - 1, final_start_b = curr_pos_b - kmer_size;
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


    //printf("From: %"PRId64" to %"PRId64"\n", start_block_a, end_block_a);
    //printf("currspos. %"PRId64"\n", curr_pos_a);
    //getchar();

    //Forward search
    while(keep_going == 1){
        
        
        if(score_right > 0 && curr_pos_a >= start_block_a && curr_pos_b < end_block_b){
            if(curr_pos_a  < start_block_a ||  curr_pos_b > end_block_b) break;
            if(compare_letters(complement(a->b->genome->seq[curr_pos_a]), b->b->genome->seq[curr_pos_b]) == POINT){ score_right+=POINT; idents++; }else{ score_right-=POINT;}
            if(high_right <= score_right){
                final_end_a = curr_pos_a;
                high_right = score_right;
                final_idents = idents;
            }
            curr_pos_a--;
            curr_pos_b++;
        }else{
            keep_going = 0;
        }
    }

    keep_going = 1;
    curr_pos_a = a->pos;
    curr_pos_b = b->pos - kmer_size - 1;

    score_left = high_right;

    //printf("From: %"PRId64" to %"PRId64"\n", start_block_a, end_block_a);
    //printf("currspos. %"PRId64"\n", curr_pos_a);
    //printf("Aktuelle: %"PRId64" and %"PRId64"\n", final_start_a, final_end_a);
    //getchar();

    //Backward search
    while(keep_going == 1){
        
        if(score_left > 0 && curr_pos_a < end_block_a && curr_pos_b >= start_block_b){
            if(curr_pos_a >= end_block_a || curr_pos_b < start_block_b ) break;
            if(compare_letters(complement(a->b->genome->seq[curr_pos_a]), b->b->genome->seq[curr_pos_b]) == POINT){ score_left+=POINT; idents++; }else{ score_left-=POINT;}
            if(high_left <= score_left){
                final_start_a = curr_pos_a;
                final_start_b = curr_pos_b;
                high_left = score_left;
                final_idents = idents;
            }
            curr_pos_a++;
            curr_pos_b--;
        }else{
            keep_going = 0;
        }
    }

    //printf("Ending\n");
    //printf("Aktuelle: %"PRId64" and %"PRId64"\n", final_start_a, final_end_a);
    //getchar();

    qf->t_len = final_start_a - final_end_a + 1;
    qf->sim = (final_idents / (long double) qf->t_len)*100;
    qf->x_start = final_end_a;
    qf->y_start = final_start_b;
    qf->diag = (int64_t)qf->x_start - (int64_t)qf->y_start;
    qf->x = a->b->genome;
    qf->y = b->b->genome;
    /*
    printf("Aligned\n");
    seq_man->print_sequence_region(a->genome->id, final_start_a, final_end_a);
    seq_man->print_sequence_region(b->genome->id, final_start_b, final_start_b+qf->t_len);
    getchar();
    */

}

void inplace_dna_switch(char *a, char *b, uint64_t l1, uint64_t l2, uint64_t t){
    uint64_t i;
    char aux;
    for(i=0;i<t;i++){
        aux = a[i+l1];
        a[i+l1] = b[i+l2];
        b[i+l2] = aux;
    }
}

void inplace_reverse_and_complement(char *d, uint64_t l){
    uint64_t i;
    char c;
    for(i=0;i<l/2;i++){
        c = d[i];
        d[i] = d[l-i-1];
        d[l-i-1] = c;
    }
    for(i=0;i<l;i++){
        d[i] = complement(d[i]);
    }
}


inline void strrev(char *p, char *d, uint32_t k){
    uint32_t i;
    char c;
    for(i=0;i<k;i++){
    	c = p[k-i-1];
    	switch(c){
    	case ('A'): c=('T');
    	break;
    	case ('C'): c=('G');
    	break;
    	case ('G'): c=('C');
    	break;
    	case ('T'): c=('A');
    	break;
    	}
        d[i] = c; 
    }
  
}

void align_region(sequence_manager * seq_man, Synteny_block * unlinked_sb, uint32_t kmer_size, dictionary_hash * dhw, Quickfrag * result){
    //Restore dictionary
    dhw->clear();

    //Kmer reading
    char curr_kmer[kmer_size], rev_kmer[kmer_size];
    uint64_t kmer_index = 0;
    char c;
    uint64_t precomputed_sizeofQuickfrag = sizeofQuickfrag();
    bool first = true;

    //To get hits
    Wordbucket * hit;
    Wordbucket ** hit_list;

    //Current pos in sequence
    uint64_t advanced_steps;
    uint64_t i;

    //Results of alignment
    Quickfrag aligned_qf;
    int64_t curr_diag;
    Word align_word;

    //Current unlinked block
    Synteny_block * unlinked = unlinked_sb;

    //For all sequences in the unlinked synteny block
    while(unlinked != NULL){
        while(advanced_steps < unlinked->b->end){
            //Get nucleotide
            c = unlinked->b->genome->seq[advanced_steps];
            advanced_steps++;

            if(c == 'A' || c == 'C' || c == 'T' || c == 'G'){
                curr_kmer[kmer_index] = c;
                kmer_index++;
            }else{
                kmer_index = 0;
            }

            //Check if we have a kmer big enough
            if(kmer_index == kmer_size){

                //Insert word in dictionary 
                hit_list = dhw->put_and_hit(curr_kmer, 'f', advanced_steps, unlinked->b);

                for(i=0;i<dhw->get_candidates_number();i++){
                    //For every hit candidate
                    hit = hit_list[i];

                    if(hit != NULL){
                        
                        //We have a hit we should try an alignment
                        //only if the hit is not overlapping
                        //or it is overlapping but on different diagonal
                        if(first){
                            
                            //There is no frag yet, so try first one
                            align_word.pos = advanced_steps;
                            align_word.strand = 'f';
                            align_word.b = unlinked->b;
                            if(hit->w.strand == 'f'){
                                alignment_from_hit(seq_man, &align_word, &hit->w, result, kmer_size);
                            }else{
                                alignment_from_hit_reverse(seq_man, &hit->w, &align_word, result, kmer_size);
                            }
                            first = false;
                        }else{

                            //There is already a frag, check for overlapping and diagonal
                            if(overlapped_words(result->x_start, result->x_start+result->t_len, advanced_steps-kmer_size-1, advanced_steps-1) != 0){
                                //If it is not overlapped 
                                curr_diag = (int64_t) advanced_steps - (int64_t) hit->w.pos;
                                if(curr_diag != result->diag){
                                    //We can try new alignment
                                    align_word.pos = advanced_steps;
                                    align_word.strand = 'f';
                                    align_word.b->genome = unlinked->b->genome;
                                    if(hit->w.strand == 'f'){
                                        alignment_from_hit(seq_man, &align_word, &hit->w, &aligned_qf, kmer_size);
                                    }else{
                                        alignment_from_hit_reverse(seq_man, &hit->w, &align_word, &aligned_qf, kmer_size);
                                    }
                                    
                                    //Only copy if new alignment is better 

                                    if(aligned_qf.sim > result->sim){
                                        memcpy(result, &aligned_qf, precomputed_sizeofQuickfrag);
                                    }
                                }
                            }
                        }
                        
                    }
                }

                //Insert reversed word in dictionary
                strrev(curr_kmer, rev_kmer, kmer_size);

                hit_list = dhw->put_and_hit(rev_kmer, 'r', advanced_steps, unlinked->b);

                for(i=0;i<dhw->get_candidates_number();i++){
                    //For every hit candidate
                    hit = hit_list[i];

                    if(hit->w.strand == 'r') continue; //Skip since a 'r' - 'r' word is the same as 'f' - 'f'

                    if(hit != NULL){
                        
                        //We have a hit we should try an alignment
                        //only if the hit is not overlapping
                        //or it is overlapping but on different diagonal
                        if(first){
                            
                            //There is no frag yet, so try first one
                            align_word.pos = advanced_steps;
                            align_word.strand = 'r';
                            align_word.b = unlinked->b;
                            alignment_from_hit_reverse(seq_man, &align_word, &hit->w, result, kmer_size);
                            //printQuickfrag(&aligned_qf);
                            memcpy(result, &aligned_qf, precomputed_sizeofQuickfrag);
                            first = false;
                        }else{
                            //printf("on the else track\n");
                            //getchar();
                            //There is already a frag, check for overlapping and diagonal
                            if(overlapped_words(result->x_start, result->x_start+result->t_len, advanced_steps-kmer_size-1, advanced_steps-1) != 0){
                                //If it is not overlapped 
                                curr_diag = (int64_t) advanced_steps - (int64_t) hit->w.pos;
                                if(curr_diag != result->diag){
                                    //We can try new alignment
                                    align_word.pos = advanced_steps;
                                    align_word.strand = 'r';
                                    align_word.b->genome = unlinked->b->genome;
                                    alignment_from_hit_reverse(seq_man, &align_word, &hit->w, &aligned_qf, kmer_size);
                                    
                                    if(aligned_qf.sim > result->sim){
                                        memcpy(result, &aligned_qf, precomputed_sizeofQuickfrag);
                                    }
                                }
                            }
                        }
                        
                    }
                }
                //Displace 
                memmove(&curr_kmer[0], &curr_kmer[1], kmer_size-1);
                memmove(&rev_kmer[0], &rev_kmer[1], kmer_size-1);
                kmer_index--;
            }

        }
        //Advance block
        unlinked = unlinked->next;
        kmer_index = 0;
    }
    

    
    

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
    char curr_kmer[kmer_size], rev_kmer[kmer_size];
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
                            if(hit->w.strand == 'f'){
                                alignment_from_hit(seq_man, &align_word, &hit->w, &aligned_qf, kmer_size);
                            }else{
                                alignment_from_hit_reverse(seq_man, &hit->w, &align_word, &aligned_qf, kmer_size);
                            }
                            
                            //printQuickfrag(&aligned_qf);
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
                                    if(hit->w.strand == 'f'){
                                        alignment_from_hit(seq_man, &align_word, &hit->w, &aligned_qf, kmer_size);
                                    }else{
                                        alignment_from_hit_reverse(seq_man, &hit->w, &align_word, &aligned_qf, kmer_size);
                                    }
                                    //printQuickfrag(&aligned_qf);
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
                strrev(curr_kmer, rev_kmer, kmer_size);

                hit_list = dhw->put_and_hit(rev_kmer, 'r', advanced_steps, sb_ptr->b);

                for(i=0;i<dhw->get_candidates_number();i++){
                    //For every hit candidate
                    hit = hit_list[i];

                    if(hit->w.strand == 'r') continue; //Skip since a 'r' - 'r' word is the same as 'f' - 'f'

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
                            align_word.strand = 'r';
                            align_word.b = sb_ptr->b;
                            alignment_from_hit_reverse(seq_man, &align_word, &hit->w, &aligned_qf, kmer_size);
                            //printQuickfrag(&aligned_qf);
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
                                    align_word.strand = 'r';
                                    align_word.b->genome = sb_ptr->b->genome;
                                    alignment_from_hit_reverse(seq_man, &align_word, &hit->w, &aligned_qf, kmer_size);
                                    //printQuickfrag(&aligned_qf);
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




                //Displace 
                memmove(&curr_kmer[0], &curr_kmer[1], kmer_size-1);
                memmove(&rev_kmer[0], &rev_kmer[1], kmer_size-1);
                kmer_index--;
            }

        }

        //Advance block
        sb_ptr = sb_ptr->next;
        kmer_index = 0;
    }
    fprintf(stdout, "Synteny level: %"PRIu64"\n", sbl->synteny_level);
    printSyntenyBlock(sbl->sb);
    for(i=0;i<seq_man->get_number_of_sequences();i++){
        for(j=0;j<seq_man->get_number_of_sequences();j++){
            if(qfmat_state[i][j] == 1) qfmat[i][j].sim = 100.0 - qfmat[i][j].sim;
        }
    }
    printQuickFragMatrix(qfmat, qfmat_state, seq_man->get_number_of_sequences());
    //getchar();
    
}

//void read_words_from_synteny_block_and_align(sequence_manager * seq_man, Synteny_list * sbl, uint32_t kmer_size, dictionary_hash * dhw, Quickfrag ** qfmat, unsigned char ** qfmat_state)

void fill_quickfrag_matrix_NW(sequence_manager * seq_man, char * seq_for_reverse, Synteny_list * sbl, Quickfrag ** qfmat, unsigned char ** qfmat_state, int iGap, int eGap, struct cell * mc, struct cell * f0, struct cell * f1){
    Synteny_block * sb_ptr = sbl->sb;
    Synteny_block * sb_ptr_intern;
    struct cell alignment;
    struct cell alignment_reverse;
    Quickfrag qf;

    uint64_t i,j;
    for(i=0;i<seq_man->get_number_of_sequences();i++){
        for(j=0;j<seq_man->get_number_of_sequences();j++) qfmat_state[i][j] = 0;
    }

    while(sb_ptr != NULL){
        sb_ptr_intern = sb_ptr->next;
        while(sb_ptr_intern != NULL){
            printf("aligning :"); printBlock(sb_ptr->b); printBlock(sb_ptr_intern->b);
            alignment = NWscore2rows(&sb_ptr->b->genome->seq[sb_ptr->b->start], 0, sb_ptr->b->end - sb_ptr->b->start, &sb_ptr_intern->b->genome->seq[sb_ptr_intern->b->start], 0, sb_ptr_intern->b->end - sb_ptr_intern->b->start, iGap, eGap, mc, f0, f1);
            
            memcpy(&seq_for_reverse[sb_ptr_intern->b->start], &sb_ptr_intern->b->genome->seq[sb_ptr_intern->b->start], sb_ptr_intern->b->end - sb_ptr_intern->b->start);
            inplace_reverse_and_complement(&seq_for_reverse[sb_ptr_intern->b->start], sb_ptr_intern->b->end - sb_ptr_intern->b->start);
            alignment_reverse = NWscore2rows(&sb_ptr->b->genome->seq[sb_ptr->b->start], 0, sb_ptr->b->end - sb_ptr->b->start, &seq_for_reverse[sb_ptr_intern->b->start], 0, sb_ptr_intern->b->end - sb_ptr_intern->b->start, iGap, eGap, mc, f0, f1);
            
           
            if(alignment_reverse.ident > alignment.ident) alignment = alignment_reverse;
            qf.x_start = alignment.xs;
            qf.y_start = alignment.ys;
            qf.t_len = alignment.xe - alignment.xs;
            qf.diag = 0; // Not needed
            qf.sim = 100 - (long double)(alignment.ident * 100) / (long double) qf.t_len;
            qf.x = sb_ptr->b->genome;
            qf.y = sb_ptr_intern->b->genome;

            printf("Best is:\n");
            printCell(&alignment);

            memcpy(&qfmat[sb_ptr->b->genome->id][sb_ptr_intern->b->genome->id], &qf, sizeofQuickfrag());
            memcpy(&qfmat[sb_ptr_intern->b->genome->id][sb_ptr->b->genome->id], &qf, sizeofQuickfrag());
            qfmat_state[sb_ptr->b->genome->id][sb_ptr_intern->b->genome->id] = 1;
            qfmat_state[sb_ptr_intern->b->genome->id][sb_ptr->b->genome->id] = 1;
            sb_ptr_intern = sb_ptr_intern->next;
        }
        sb_ptr = sb_ptr->next;
    }
    printQuickFragMatrix(qfmat, qfmat_state, seq_man->get_number_of_sequences());
}

Slist * UPGMA_joining_clustering(Quickfrag ** M, double ** submat, unsigned char ** qfmat_state, uint64_t N, memory_pool * mp){


    //printQuickFragMatrix(M, qfmat_state, 2);

    uint64_t i, j, i_p, j_p;
    double f_min = -1;
    uint64_t i_min = 0, j_min = 1;
    uint64_t nodes_in_dendro = 1;
    
    Slist ** dendrogram = (Slist **) mp->request_bytes(N*sizeof(Slist*));

    uint64_t * dendro_lengths = (uint64_t *) mp->request_bytes(N*sizeof(uint64_t));

    unsigned char * skip_i = (unsigned char *) mp->request_bytes(N*sizeof(unsigned char));
    for(i=0;i<N;i++){
        dendrogram[i] = NULL;
        skip_i[i] = 0;
        dendro_lengths[i] = 1;
    }
    

    //Generate submatrix, since some rows or columns will be null
    i_p = 0; j_p = 0;
    for(i=0;i<N;i++){
        j_p = 0;
        for(j=0;j<N;j++){
            if(qfmat_state[i][j] == 1){
                // Copy value
                if(i_p == j_p) j_p++;
                submat[i_p][j_p] = M[i][j].sim;
                if(dendrogram[i_p] == NULL){
                    dendrogram[i_p] = (Slist *) mp->request_bytes(sizeofSlist());
                    if(M[i][j].x->id == i) dendrogram[i_p]->s = M[i][j].x;
                    if(M[i][j].y->id == i) dendrogram[i_p]->s = M[i][j].y;
                    dendrogram[i_p]->next = NULL;
                    
                }
                j_p++;
            }
        }
        
        if(j_p > 0) i_p++;
    }

    /*
    for(i=0;i<i_p;i++){
        printf("ID: %"PRIu64"\n", dendrogram[i]->s->id);
    }
    */

    j_p = i_p;
    while(nodes_in_dendro < i_p){

        printUnstatedDoubleMatrix(submat, i_p, skip_i);

        //Find smallest
        for(i=0;i<i_p;i++){
            //Jump those that have been removed
            if(skip_i[i] == 1) continue;

            for(j=0;j<i;j++){
                //Jump removed & Find minimum
                if(skip_i[j] == 0 && (f_min == -1 || f_min >= submat[i][j])){
                    f_min = submat[i][j];
                    i_min = i;
                    j_min = j;
                }
                
                
            }
        }

        //Join minimums
        //printf("Min on %"PRIu64", %"PRIu64"\n", i_min, j_min);

        
        if( dendrogram[j_min]->next == NULL || dendrogram[i_min]->next == NULL){ 
            //Add one makes a tuple
            Slist * aux = dendrogram[i_min];
            while(aux->next != NULL){
                aux = aux->next;
            }

            aux->next = dendrogram[j_min];
            //last_d->next = dendrogram[j_min];
            //last_d = dendrogram[j_min];
            dendro_lengths[i_min]++;
        }else{
            Slist * null_seq = (Slist *) mp->request_bytes(sizeofSlist());
            null_seq->s = NULL; 
            null_seq->next = dendrogram[j_min];

            Slist * aux = dendrogram[i_min];
            while(aux->next != NULL){
                aux = aux->next;
            }
            aux->next = null_seq;
            dendro_lengths[i_min] += dendro_lengths[j_min];                
            
        }
        

        

        //Remove existing rows and cols
        skip_i[j_min] = 1;

        //Update new distances
        for(i=0;i<i_p;i++){
            if(skip_i[i] == 0){
                submat[i_min][i] = ((submat[i_min][i] + submat[j_min][i]) * dendro_lengths[i_min]) /2;
                submat[i][i_min] = submat[i_min][i];
            }
        }
        /*
        for(i=0;i<i_p;i++){
            if(skip_i[i] == 0){
                submat[j_min][i] = ((submat[i_min][i] + submat[j_min][i]) * dendro_lengths[i_min]) /2;
            }
        }
        */
        
        nodes_in_dendro++;
        

        f_min = -1; //Restart minimum search
        
    }

    //Print clusters
    printf("Phylogenetic Clustering: \n");printDendrogramList(dendrogram[i_min]);
    printf("\n");


    return dendrogram[i_min]; 
    
}

void find_event_location(Slist * dendrogram, Event e, void * data, bool * genomes_affected){
    Slist * d = dendrogram;
    Slist * prev, * next;
    
    while(d != NULL){
        //if(d->s == NULL) printf(" ] "); else printf(" %"PRIu64" ", d->s->id);
        prev = d;
        d = d->next;
        if(d != NULL) next = d->next;

        if(prev != NULL && d != NULL && next != NULL){
            if(prev->s != NULL && d->s != NULL && next->s != NULL){
                switch(e){
                    case inversion: 
                    {
                        //check for strands
                        
                        strand_matrix * sm_B = (strand_matrix *) data;
                        printf("%u %u %u %u\n", sm_B->get_strands(next->s->id, d->s->id), sm_B->get_strands(d->s->id, prev->s->id), sm_B->get_strands(next->s->id, prev->s->id), sm_B->get_strands(prev->s->id, d->s->id));
                        sm_B->print_strand_matrix();
                        if(sm_B->get_strands(next->s->id, d->s->id) == sm_B->get_strands(d->s->id, prev->s->id)){
                            //The reversion happened in 'prev'
                            genomes_affected[prev->s->id] = true;
                        }                        
                        if(sm_B->get_strands(next->s->id, prev->s->id) == sm_B->get_strands(prev->s->id, d->s->id)){
                            //The reversion happened in 'd'
                            genomes_affected[d->s->id] = true;
                        }

                    }
                    break;
                }
            }else{
                //Collapse
                // TODO                
                
            }
        }
    }
}

int compare_distances_indel(const void * a, const void * b){
    
    if( *(uint64_t *) a > *(uint64_t *) b) return 1;
    if(*(uint64_t *) a == *(uint64_t *) b) return 0;
    return -1;
}

long double median_from_vector(uint64_t * v, uint64_t l){
    if(l % 2 != 0){
        //its not pair
        return (long double) v[l/2];
    }else{
        //Pair
        return ((long double) v[(l-1)/2] + (long double) v[(l-1)/2 + 1])/2;
    }
}

