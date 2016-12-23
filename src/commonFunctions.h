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

#endif /* COMMON_FUNCTIONS_H */
