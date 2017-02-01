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
	To compute sizeof Synteny block when not using padding
*/
uint64_t sizeofSyntenyBlock();

/*
	To compute sizeof Synteny block list when not using padding
*/
uint64_t sizeofSyntenyList();

/*
	To compute sizeof Annotation when not using padding
*/
uint64_t sizeofAnnotation();

/*
	To compute sizeof Word when not using padding
*/
uint64_t sizeofWord();
/*
	To compute sizeof Wordbucket when not using padding
*/
uint64_t sizeofWordbucket();
/*
	Check if two fragments are equal based on sequence ids, coordinates and strand
*/
int isFragmentEqualTo(struct FragFile * a, struct FragFile * b);

/*
	Check if two blocks are equal
*/
int isBlockEqualTo(Block * a, Block * b);

/*
	Checks if a SEQ_ID is contained in a list of frags
*/
int idNotInList(Frags_list * fl, struct FragFile * f);

/*
	Prints a FragFile to stdout
*/
void printFragment(struct FragFile * f);

/*
	Prints a Block to stdout
*/
void printBlock(Block * b);

/*
	Prints a Block to stdout with its fragments
*/
void printFragsFromBlock(Block * b);

/*
	Prints a Block to stdout without tags
*/
void printBlockWriteMode(Block * b);

/*
	Prints a Synteny Block to stdout
*/
void printSyntenyBlock(Synteny_block * b);

/*
	Prints all frags that correspond to a synteny block
	for all synteny blocks in a synteny node
*/
void printSyntenyListNode(Synteny_list * sbl);

/*
	Prints an annotation
*/
void printAnnotation(Annotation * a);

#endif /* COMPARISON_FUNCTIONS_H */
