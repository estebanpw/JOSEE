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
	Load sequences lengths and compute their accumulated values
	Returns the number of sequences loaded
*/
uint64_t load_sequences_descriptors(Sequence ** sequences, FILE * lengths_file);

/*
	Load all fragments and transform global coordinates to local
	n_frags should keep the number of fragments read
*/
void load_fragments_local(FILE * fragsfile, uint64_t * n_frags, Sequence * sequences, struct FragFile ** loaded_frags);

/*
	Computes the longest sequence
*/
uint64_t get_maximum_length(Sequence * sequences, uint64_t n_seqs);

/*
	Computes the coverage per sequence
*/
void get_coverage_from_genome_grid(unsigned char ** maptable, Sequence * sequences, uint64_t n_seqs, uint64_t min_len_without_breaks);

/*
	Shows length and coverage per sequence
*/
void print_sequences_data(Sequence * sequences, uint64_t n_seqs);

/*
	Writes all rows of a maptable to a file path that is created
	This function is intended for debug
*/
void write_maptable_to_disk(unsigned char ** maptable, uint64_t n_seqs, Sequence * sequences, const char * out_file_path);

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
#endif /* COMMON_FUNCTIONS_H */
