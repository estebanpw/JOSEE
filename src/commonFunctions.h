#ifndef COMMON_FUNCTIONS_H
#define COMMON_FUNCTIONS_H
#include "structs.h"
/**
 * Print the error message 's' and exit(-1)
 */
void terror(const char *s);

/**
 * Function to read char by char buffered from a FILE
 */
char buffered_fgetc(char *buffer, uint64_t *pos, uint64_t *read, FILE *f);

/*
	Load all fragments and transform global coordinates to local
	n_frags should keep the number of fragments read
*/
void load_fragments_local(FILE * fragsfile, uint64_t * n_frags, struct FragFile ** loaded_frags);


/*
	Computes the coverage per sequence
*/
void get_coverage_from_genome_grid(unsigned char ** maptable, sequence_manager * seq_manager, uint64_t n_seqs, uint64_t min_len_without_breaks);


/*
	Writes all rows of a maptable to a file path that is created
	This function is intended for debug
*/
void write_maptable_to_disk(unsigned char ** maptable, sequence_manager * seq_manager, const char * out_file_path);

/*
	Prints a section of the maptable
	This function is intended for debug
*/
void print_maptable_portion(unsigned char ** maptable, uint64_t from, uint64_t to, uint64_t rate, uint64_t seq);

/*
	Prints the blocks that compose a synteny block
*/
void traverse_synteny_list(Synteny_list * sbl);

/*
	Travels the synteny list until finding the synteny block to which the
	given block corresponds to
*/
Synteny_list * find_synteny_block_from_block(Synteny_list * sbl, Block * b);

/*
	Looks for fragments in the maptable that start at the given positions in the given sequence
	WARNING: Performs linear search.
*/
void find_fragments_from_maptable(unsigned char ** maptable, uint64_t start, uint64_t end, uint64_t seq, struct FragFile * frags, uint64_t n_frags);

/*
	Compares the range of two annotations
*/
int compare_ranges(Annotation * a, Annotation * b);

/*
	Compares two annotations by starting coordinates
*/
int compare_two_annotations(Annotation * a, Annotation * b);

/*
	Binary searches an annotation array looking for ranges that overlap
*/
Annotation * binary_search_annotations(uint64_t start, uint64_t end, Annotation * anot, uint64_t n_annots);

/*
	Sorts an array of of annotations
*/
void quick_sort_annotations(Annotation * array, uint64_t x, uint64_t y);

/**
 *	Compute the power of n to for with a lookup table
 */
uint64_t quick_pow4(uint32_t n);

/*
 *	Compute the power of 0,1,2,3 to the power of n quickly by lookup table
*/
uint64_t quick_pow4byLetter(uint32_t n, const char c);

/**
 *	Compute the hash of a word of letters of length k (in ascii)
 */
uint64_t hashOfWord(const char * word, uint32_t k);

/**
 *	Compares two nucleotides and returns either POINT or -POINT
 */
inline int64_t compare_letters(char a, char b);

/**
 *	Checks if two regions are overlapped
 */
int overlapped_words(uint64_t xstart, uint64_t xend, uint64_t ystart, uint64_t yend);

/**
 *	Ungapped alignment in hit from words a and b
 */
void alignment_from_hit(sequence_manager * seq_man, Word * a, uint64_t start_block_a, uint64_t end_block_a, Word * b, uint64_t start_block_b, uint64_t end_block_b, Quickfrag * qf, uint64_t kmer_size);

/**
 *	Ungapped alignment in hit from words a (REVERSED) and b
 */
void alignment_from_hit_reverse(sequence_manager * seq_man, Word * a, Word * b, Quickfrag * qf, uint64_t kmer_size);

/*
	Aligns two unlinked blocks (i.e. take coordinates and put them as blocks in a synteny block)
*/
void align_region(sequence_manager * seq_man, Synteny_block * unlinked_sb, uint32_t kmer_size, dictionary_hash * dhw, Quickfrag * result);
/*
	Produces the dictionary and aligns the blocks in a synteny block
	The qfmat stores the best fragments
*/
void read_words_from_synteny_block_and_align(sequence_manager * seq_man, Synteny_list * sbl, uint32_t kmer_size, dictionary_hash * dhw, Quickfrag ** qfmat, unsigned char ** qfmat_state);

/*
	Same as read_words_and_align but for NW 
*/
void fill_quickfrag_matrix_NW(sequence_manager * seq_man, char * seq_for_reverse, Synteny_list * sbl, Quickfrag ** qfmat, unsigned char ** qfmat_state, int iGap, int eGap, struct cell * mc, struct cell * f0, struct cell * f1);
/*
	Complements a nucleotide 
*/
inline char complement(char c);

/*
	Swaps two DNA segments 
*/
void inplace_dna_switch(char *a, char *b, uint64_t l1, uint64_t l2, uint64_t t);
/*
	Reverses and complements a sequence in place
*/
void inplace_reverse_and_complement(char *d, uint64_t l);
/*
	Computes the reverse complementary of 'p' and returns it in 'd' for lenght 'k'
*/
inline void strrev(char *p, char *d, uint32_t k);

/*
	Computes a phylogenetic clustering based on a distance matrix using the UPGMA method (assumes equidistances)
	returns the id of the genome where the reversion happened
	or -1 if it cant be decided
*/
Slist * UPGMA_joining_clustering(Quickfrag ** M, double ** submat, unsigned char ** qfmat_state, uint64_t N, memory_pool * mp);

/*
	Using the clustered phylogenetic tree, finds where did the event happen (if possible)
*/
void find_event_location(Slist * dendrogram, Event e, void * data, bool * genomes_affected);
/*
	Compare function for distances to sort
*/
int compare_distances_indel(const void * a, const void * b);

/*
	Computes the median from a sorted vector
*/
long double median_from_vector(uint64_t * v, uint64_t l);

#endif /* COMMON_FUNCTIONS_H */
