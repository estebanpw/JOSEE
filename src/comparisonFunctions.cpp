#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <ctype.h>
#include <arpa/inet.h>
#include <pthread.h>
#include <math.h>
#include "structs.h"
#include "commonFunctions.h"



void endianessConversion(char *source, char *target, int numberOfBytes) {
    int i, j;
    for (i = numberOfBytes - 1; i >= 0; i--) {
        j = numberOfBytes - 1 - i;
        target[j] = source[i];
    }
}


/**
 * Function to read a fragment from the specified file
 */
void readFragment(struct FragFile *frag, FILE *f) {
    char tmpArray[sizeof(long double)];

    if (htons(1) == 1) {
        //big endian
        if (fread(&frag->diag, sizeof(int64_t), 1, f) != 1) {
            if (feof(f))return;
            terror("Error reading the HSP diagonal");
        }
        if (fread(&frag->xStart, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP xStart");
        }
        if (fread(&frag->yStart, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP yStart");
        }
        if (fread(&frag->xEnd, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP xEnd");
        }
        if (fread(&frag->yEnd, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP yEnd");
        }
        if (fread(&frag->length, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP length");
        }
        if (fread(&frag->ident, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP identities");
        }
        if (fread(&frag->score, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP score");
        }
        if (fread(&frag->similarity, sizeof(float), 1, f) != 1) {
            terror("Error reading the HSP similarity");
        }
        if (fread(&frag->seqX, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP seqX");
        }
        if (fread(&frag->seqY, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP seqY");
        }
        if (fread(&frag->block, sizeof(int64_t), 1, f) != 1) {
            terror("Error reading the HSP block");
        }
        frag->strand = fgetc(f);
        if (fread(&frag->evalue, sizeof(long double), 1, f) != 1) {
            terror("Error reading the HSP evalue");
        }
    } else {
        //little endian
        if (fread(tmpArray, sizeof(int64_t), 1, f) != 1) {
            if (feof(f))return;
            terror("Error reading the HSP diagonal");
        }
        endianessConversion(tmpArray, (char *) (&frag->diag), sizeof(int64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP xStart");
        }
        endianessConversion(tmpArray, (char *) (&frag->xStart), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP yStart");
        }
        endianessConversion(tmpArray, (char *) (&frag->yStart), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP xEnd");
        }
        endianessConversion(tmpArray, (char *) (&frag->xEnd), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP yEnd");
        }
        endianessConversion(tmpArray, (char *) (&frag->yEnd), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP length");
        }
        endianessConversion(tmpArray, (char *) (&frag->length), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP identity");
        }
        endianessConversion(tmpArray, (char *) (&frag->ident), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP score");
        }
        endianessConversion(tmpArray, (char *) (&frag->score), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(float), 1, f) != 1) {
            terror("Error reading the HSP float");
        }
        endianessConversion(tmpArray, (char *) (&frag->similarity), sizeof(float));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP seqX");
        }
        endianessConversion(tmpArray, (char *) (&frag->seqX), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading the HSP seqY");
        }
        endianessConversion(tmpArray, (char *) (&frag->seqY), sizeof(uint64_t));
        if (fread(tmpArray, sizeof(int64_t), 1, f) != 1) {
            terror("Error reading the HSP block");
        }
        endianessConversion(tmpArray, (char *) (&frag->block), sizeof(int64_t));
        frag->strand = fgetc(f);
        if (fread(tmpArray, sizeof(long double), 1, f) != 1) {
            terror("Error reading the HSP evalue");
        }
        endianessConversion(tmpArray, (char *) (&frag->evalue), sizeof(long double));
    }
}

/**
 * Function to write a fragment to the specified file
 */
void writeFragment(struct FragFile *frag, FILE *f) {
    char tmpArray[sizeof(long double)];
    if (htons(1) == 1) {
        //Big endian
        fwrite(&frag->diag, sizeof(int64_t), 1, f);
        fwrite(&frag->xStart, sizeof(uint64_t), 1, f);
        fwrite(&frag->yStart, sizeof(uint64_t), 1, f);
        fwrite(&frag->xEnd, sizeof(uint64_t), 1, f);
        fwrite(&frag->yEnd, sizeof(uint64_t), 1, f);
        fwrite(&frag->length, sizeof(uint64_t), 1, f);
        fwrite(&frag->ident, sizeof(uint64_t), 1, f);
        fwrite(&frag->score, sizeof(uint64_t), 1, f);
        fwrite(&frag->similarity, sizeof(float), 1, f);
        fwrite(&frag->seqX, sizeof(uint64_t), 1, f);
        fwrite(&frag->seqY, sizeof(uint64_t), 1, f);
        fwrite(&frag->block, sizeof(int64_t), 1, f);
        fputc(frag->strand, f);
        fwrite(&frag->evalue, sizeof(long double), 1, f);
    } else {
        //Little endian
        endianessConversion((char *) (&frag->diag), tmpArray, sizeof(int64_t));
        fwrite(tmpArray, sizeof(int64_t), 1, f);
        endianessConversion((char *) (&frag->xStart), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->yStart), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->xEnd), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->yEnd), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->length), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->ident), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->score), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->similarity), tmpArray, sizeof(float));
        fwrite(tmpArray, sizeof(float), 1, f);
        endianessConversion((char *) (&frag->seqX), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->seqY), tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
        endianessConversion((char *) (&frag->block), tmpArray, sizeof(int64_t));
        fwrite(tmpArray, sizeof(int64_t), 1, f);
        fputc(frag->strand, f);
        endianessConversion((char *) (&frag->evalue), tmpArray, sizeof(long double));
        fwrite(tmpArray, sizeof(long double), 1, f);
    }
}

/**
 * Function to read the sequence length
 */
void readSequenceLength(uint64_t *length, FILE *f) {
    char tmpArray[8];
    if (htons(1) == 1) {
        //big endian
        if (fread(length, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading sequence length");
        }
    } else {
        //little endian
        if (fread(tmpArray, sizeof(uint64_t), 1, f) != 1) {
            terror("Error reading sequence length");
        }
        endianessConversion(tmpArray, (char *) length, sizeof(uint64_t));
    }
}

/**
 * Function to write the sequence length
 */
void writeSequenceLength(uint64_t *length, FILE *f) {
    char tmpArray[8];
    if (htons(1) == 1) {
        //big endian
        fwrite(length, sizeof(uint64_t), 1, f);
    } else {
        //little endian
        endianessConversion((char *) length, tmpArray, sizeof(uint64_t));
        fwrite(tmpArray, sizeof(uint64_t), 1, f);
    }
}

uint64_t sizeofFragment() {
    return 9 * sizeof(uint64_t) + 2 * sizeof(int64_t) + 1 * sizeof(float) + 1 * sizeof(char) + 1 * sizeof(long double);
}

uint64_t sizeofSequence() {
    return 3 * sizeof(uint64_t) + sizeof(uint32_t) + sizeof(char *) + sizeof(uint64_t);
}

uint64_t sizeofBlock(){
    return  4*sizeof(uint64_t) + sizeof(Sequence *) + sizeof(Frags_list *) + 1*sizeof(unsigned char) + sizeof(Synteny_list *) + 2*sizeof(Block *);
}

uint64_t sizeofFrags_list(){
    return sizeof(struct FragFile *) + sizeof(struct frags_list *);
}

uint64_t sizeofBucket(){
    return sizeofBlock() + sizeof(struct bucket *) + sizeof(Frags_list *);
}

uint64_t sizeofSyntenyBlock(){
    return sizeof(Block *) + sizeof(Synteny_block *);
}

uint64_t sizeofSyntenyList(){
    return sizeof(Synteny_block *) + 2*sizeof(Synteny_list *) + 2*sizeof(uint64_t);
}

uint64_t sizeofAnnotation(){
    return 2*sizeof(uint64_t) + sizeof(char) + sizeof(char *);
}

uint64_t sizeofWord(){
    return 2*sizeof(uint64_t) + sizeof(Block *) + sizeof(char);
}

uint64_t sizeofWordbucket(){
    return sizeofWord() + sizeof(Wordbucket *);
}

uint64_t sizeofQuickfrag(){
    return 3*sizeof(uint64_t) + sizeof(long double) + sizeof(int64_t) + 2*sizeof(Sequence *);
}

uint64_t sizeofSlist(){
    return sizeof(Sequence *) + sizeof(Slist *);
}

uint64_t sizeofE_inversion(){
    return sizeof(Block *);
}

uint64_t sizeofE_duplication(){
    return 2*sizeof(Block *);
}

uint64_t sizeofRearrangement(){
    return 4*sizeof(uint64_t) + sizeof(unsigned char) + 2*sizeof(int64_t);
}

int isFragmentEqualTo(struct FragFile * a, struct FragFile * b){
    if(a->seqX != b->seqX) return 0;
    if(a->xStart != b->xStart) return 0;
    if(a->xEnd != b->xEnd) return 0;
    return 1;
}

int isBlockEqualTo(Block * a, Block * b){
    if(a->start == b->start && a->end == b->end && a->genome->id == b->genome->id) return 1;
    return 0;
}

int isBlockEqualToWithOrder(Block * a, Block * b){
    if(a->start == b->start && a->end == b->end && a->genome->id == b->genome->id && a->order == b->order) return 1;
    return 0;
}

int idNotInList(Frags_list * fl, struct FragFile * f){
    Frags_list * ptr = fl;
    while(ptr != NULL){
        
        if(ptr->f->seqX == f->seqX && ptr->f->seqY == f->seqY) return 0;
        
        
        ptr = ptr->next;
    }
    return 1;
}

void printFragment(struct FragFile * f){
    
    fprintf(stdout, "FRAG::(%"PRIu64", %"PRIu64") to (%"PRIu64", %"PRIu64"): [%"PRIu64"]-[%"PRIu64"] %c LEN[%"PRIu64"]\n", f->xStart, f->yStart, f->xEnd, f->yEnd, f->seqX, f->seqY, f->strand, f->length);
}

void printBlock(Block * b){
    fprintf(stdout, "BLOCK::(%"PRIu64", %"PRIu64"): order[%"PRIu64"] len[%"PRIu64"] @genome[%"PRIu64"] #%"PRIu64"\n", b->start, b->end, b->order, b->end-b->start+1, b->genome->id, b->id);
}

void printFragsFromBlock(Block * b){
    Frags_list * fl = b->f_list;
    while(fl != NULL){
        printf("\t");printFragment(fl->f);
        fl = fl->next;
    }
}

void printBlockWriteMode(Block * b){
    fprintf(stdout, "%"PRIu64";%"PRIu64";%"PRIu64";%"PRIu64";%"PRIu64"\n", b->start, b->end, b->order, b->end-b->start+1, b->genome->id);
}

void printSyntenyBlock(Synteny_block * b){
    Synteny_block * ptr = b;
    while(ptr != NULL){
        //printBlockWriteMode(ptr->b);
        printBlock(ptr->b);
        //printFragsFromBlock(ptr->b);
        ptr = ptr->next;
    }
}

void printSyntenyListNode(Synteny_list * sbl){
    Synteny_block * sb_ptr = sbl->sb;
    while(sb_ptr != NULL){
        printFragsFromBlock(sb_ptr->b);
        sb_ptr = sb_ptr->next;
    }
    getchar();
}

void printAnnotation(Annotation * a){
    fprintf(stdout, "GENE: (%"PRIu64", %"PRIu64") %c : %s\n", a->start, a->end, a->strand, a->product);
}

void printQuickfrag(Quickfrag * qf){
    printf("[%"PRIu64", %"PRIu64"] l:%"PRIu64", s:%Le, d=%"PRId64"\n",qf->x_start, qf->y_start, qf->t_len, qf->sim, qf->diag);
}

void printQuickFragMatrix(Quickfrag ** qfmat, unsigned char ** qfmat_state, uint64_t n_seqs){
    uint64_t i,j;
    printf(" \t");
    for(i=0;i<n_seqs;i++){
        printf("%"PRIu64"\t", i);
    }
    printf("\n");
    for(i=0;i<n_seqs;i++){
        printf("%"PRIu64"\t", i);
        for(j=0;j<n_seqs;j++){
            if(qfmat_state[i][j] == 1){
                printf("[%.1f] ", (float)qfmat[i][j].sim);
            }else{
                printf("[ ** ] ");
            }
            
        }
        printf("\n");
    }
}

void printUnstatedDoubleMatrix(double ** qfmat, uint64_t n_seqs, unsigned char * skip_i){
    uint64_t i,j;
    printf(" \t");
    for(i=0;i<n_seqs;i++){
        printf("%"PRIu64"\t", i);
    }
    printf("\n");
    for(i=0;i<n_seqs;i++){
        printf("%"PRIu64"\t", i);
        for(j=0;j<n_seqs;j++){

            if(skip_i[i] == 1 || skip_i[j] == 1){
                printf("[ ** ]");
            }else{
                if(i==j) printf("[ ** ]"); else printf("[%.1f] ", (float)qfmat[i][j]);
            }

            
            
            
        }
        printf("\n");
    }
}

void printDendrogramList(Slist * dendrogram){
    Slist * d = dendrogram;
    uint64_t safe = 0;
    while(d != NULL){
        if(d->s == NULL) printf(" ] "); else printf(" %"PRIu64" ", d->s->id);
        d = d->next;
        safe++;
        if(safe > 20){ printf("Had to break\n"); break; }
    }
}

void printInvolvedGenomes(uint64_t * genomes_counters, uint64_t n_sequences){
    for(uint64_t i=0;i<n_sequences;i++){
        printf("%"PRIu64"\t", i);
    }
    printf("\n");
    for(uint64_t i=0;i<n_sequences;i++){
        printf("%"PRIu64"\t", genomes_counters[i]);
    }
}

void printDebugBlockOrderByGenome(Synteny_list * sl, uint64_t genome_id){
    if(sl != NULL){
        Synteny_block * sb_ptr = sl->sb;
        while(sb_ptr != NULL){

            if(sb_ptr->b->genome->id == genome_id){
                Block * b_ptr = sb_ptr->b;
                bool change = false;
                while(b_ptr != NULL){
                    printf("%"PRIu64", ", b_ptr->order);
                    if(b_ptr->next == NULL){
                        change = true;
                    }
                    if(!change){
                        b_ptr = b_ptr->next;
                    }else{
                        b_ptr = b_ptr->prev;
                    }
                    
                    
                    
                }
                
                
                return;
            }
            sb_ptr = sb_ptr->next;
        }
    }
    
}


void printRearrangement(rearrangement * r){
    printf("REARRANGEMENT:: MC->%"PRId64" MO->%"PRId64" APPLIES IN->[%"PRIu64", %"PRIu64"] ENDS:%"PRIu64" AFFECTS->%"PRIu64" LIFE:%d\n", r->mod_coordinates, r->mod_order, r->b1_id, r->b2_id, r->ends_at, r->affects_who, (int)r->type);
}