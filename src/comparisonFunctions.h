#ifndef COMPARISON_FUNCTIONS_H
#define COMPARISON_FUNCTIONS_H

#define min(x, y)    (((x) < (y)) ? (x) : (y))


/**
 * Function to read a fragment from the specified file
 */
void readFragment(struct FragFile *frag, FILE *f);

/**
 * Function to write a fragment to the specified file
 */
void writeFragment(struct FragFile *frag, FILE *f);

/**
 * Function to read the sequence length
 */
void readSequenceLength(uint64_t *length, FILE *f);

/**
 * Function to write the sequence length
 */
void writeSequenceLength(uint64_t *length, FILE *f);

#endif /* COMPARISON_FUNCTIONS_H */
