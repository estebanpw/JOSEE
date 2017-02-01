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

#endif /* COMMON_FUNCTIONS_H */
