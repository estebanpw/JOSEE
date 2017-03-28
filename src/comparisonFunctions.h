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
	To compute sizeof quickfrag when not using padding
*/
uint64_t sizeofQuickfrag();

/*
	To compute sizeof Slist without padding
*/
uint64_t sizeofSlist();

/*
	To compute sizeof E_inversion without padding
*/
uint64_t sizeofE_inversion();

/*
	To compute sizeof E_duplication without padding
*/
uint64_t sizeofE_duplication();

/*
	To compute sizeof Rearrangement without padding
*/
uint64_t sizeofRearrangement();

/*
	To compute sizeof Rearrangement without padding
*/
uint64_t sizeofASequence();

/*
	To compute sizeof a Cell without padding
*/
uint64_t sizeofCell();

/*
	Size of holder without padding
*/
uint64_t sizeofHolder();

/*
	Likewise for triplets.
*/
uint64_t sizeofTriplet();
/*
	Check if two fragments are equal based on sequence ids, coordinates and strand
*/
int isFragmentEqualTo(struct FragFile * a, struct FragFile * b);

/*
	Check if two blocks are equal WITHOUT COMPARING ORDER!
*/
int isBlockEqualTo(Block * a, Block * b);

/*
	Check if two blocks are equal WITH ORDER COMPARISON
*/
int isBlockEqualToWithOrder(Block * a, Block * b);

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
	Prints a block to be used in mgvisualizator format
*/
void printBlockJoseMode(Block * b, FILE * f);
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

/*
	Print a quickfrag
*/
void printQuickfrag(Quickfrag * qf);
/*
	Prints the symmetric quickfrag matrix after a synteny block is aligned
*/
void printQuickFragMatrix(Quickfrag ** qfmat, unsigned char ** qfmat_state, uint64_t n_seqs);

/*
	Print a cell struct
*/
void printCell(struct cell * c);
/*
	Print a submatrix from a quickfrags matrix
*/
void printUnstatedDoubleMatrix(double ** qfmat, uint64_t n_seqs, unsigned char * skip_i);

/*
	Prints a dendrogram as a list
*/
void printDendrogramList(Slist * dendrogram);

/*
	Prints the number of blocks involved in a synteny list
*/
void printInvolvedGenomes(uint64_t * genomes_counters, uint64_t n_sequences);

/*
	Debugging function to print all orders in synteny lists for a genome
*/
void printDebugBlockOrderByGenome(Synteny_list * sl, uint64_t genome_id);

/*
	Prints a rearrangement structure
*/
void printRearrangement(rearrangement * r);
#endif /* COMPARISON_FUNCTIONS_H */
