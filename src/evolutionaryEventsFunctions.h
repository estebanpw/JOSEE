#include <cstdarg>
#include "structs.h"
#include "commonFunctions.h"


/*
	This function maps fragments to a table of size n genomes by their lengths
	0 -> No fragment in position
	1 -> Fragment start
	2 -> Covered by fragment
	3 -> Fragment end

	Important: The input to this function should have global coordinates already removed
*/
void map_frags_to_genomes(unsigned char ** map_table, struct FragFile * frags, uint64_t n_frags, sequence_manager * seq_manager);

/*
	Copy frags properties to new snipped fragment
*/
inline void copyFragWithNewCoordinates(struct FragFile * destination, struct FragFile * source, uint64_t xStart, uint64_t yStart, uint64_t xEnd, uint64_t yEnd, uint64_t len);


/*
	This function cuts frags into more if there are other frags that overlap partially
	min_len acts as a filter to remove short fragments
*/
struct FragFile * trim_fragments_and_map(unsigned char ** map_table, struct FragFile * frags, uint64_t * n_frags, uint64_t min_len);

/*
	Iterates through the keys in the hash table and computes and stores the
	order for each block found. 
*/
void compute_order_of_blocks(hash_table * ht, uint64_t n_seqs);

/*
	Produces the list of syteny blocks
*/
Synteny_list * compute_synteny_list(hash_table * ht, uint64_t n_seqs, memory_pool * mp);



/*
	Checks that blocks in a synteny list are consecutive in respect to their genomes
	Assumes: Synteny level + same number of genomes involved
*/
bool consecutive_block_order(uint64_t * pairs_diff, uint64_t args_count, ...);

/*
	Substracts the accumulated order offset from the blocks in given synteny lists
*/
void recompute_orders_from_offset(uint64_t * orders, uint64_t args_count, ...);
/*
	Returns TRUE if there is the same number of blocks per each genome, otherwise FALSE
*/
bool genomes_involved_in_synteny(uint64_t * genomes_counters, uint64_t n_sequences, uint64_t args_count, ...);

/*
	Returns the synteny level for any given number of synteny lists
	returns 0 if the lists dont share the synteny level
	otherwise returns the level of synteny
*/
uint64_t synteny_level_across_lists(uint64_t args_count, ...);
/*
	Concatenates three synteny blocks into one
*/
void concat_synteny_blocks(Synteny_list * A, Synteny_list * B, Synteny_list * C);


/*
	Computes the strand matrix for a synteny block to compute reversions
	The strand_matrix is assumed to have length n*n where n is the number of sequences
*/
void generate_strand_matrix(Synteny_block * sb, char ** strand_matrix);

/*
	Detects candidates for evolutionary events
*/
void detect_evolutionary_event(Synteny_list * sbl, sequence_manager * seq_man, uint32_t kmer_size);