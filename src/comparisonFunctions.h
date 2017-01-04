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

/*
	To compute sizeof fragment when not using padding
*/
uint64_t sizeofFragment();

/*
	To compute sizeof sequence when not using padding
*/
uint64_t sizeofSequence();

/*
	To compute sizeof Block when not using padding
*/
uint64_t sizeofBlock();

/*
	To compute sizeof Frags_list when not using padding
*/
uint64_t sizeofFrags_list();

/*
	To compute sizeof Bucket when not using padding
*/
uint64_t sizeofBucket();

/*
	Check if two blocks are equal
*/
int isBlockEqualTo(Block * a, Block * b);

/*
	Checks if a SEQ_ID is contained in a list of frags
*/
int idNotInList(Frags_list * fl, uint64_t seq_id);

#endif /* COMPARISON_FUNCTIONS_H */
