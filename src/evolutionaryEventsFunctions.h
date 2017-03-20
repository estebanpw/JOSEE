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
Synteny_list * compute_synteny_list(hash_table * ht, uint64_t n_seqs, memory_pool * mp, uint64_t * last_s_id);

/*
	Computes the distance between blocks in two synteny lists
*/
void distance_between_blocks(uint64_t * distances, Synteny_list * A, Synteny_list * B);
/*
	Checks that blocks in a synteny list are consecutive in respect to their genomes
	Assumes: Synteny level + same number of genomes involved
*/
bool consecutive_block_order(uint64_t * pairs_diff, uint64_t args_count, ...);

/*
	Tries to separate the blocks of two synteny blocks into two groups by order
	i.e. the order differences between blocks in A and blocks in B can only take two values
	(diff 1 and diff 2) then a transposition can exist
	If result is true, then cons_order_T1 and cons_order_T2 hold the separated groups respectively
*/
bool consecutive_block_order_except_one(uint64_t * pairs_diff, uint64_t n_sequences, Block ** cons_order_T1, Block ** cons_order_T2, uint64_t args_count, ...);

/*
	Checks whether the separated groups from T1 and T2 have the same genomes
	A pointer to a block of the synteny list to be retrieved is returned
*/
Block * compare_order_clusters(Block ** cons_order_A_B_T1, Block ** cons_order_A_B_T2, Block ** cons_order_B_C_T1, Block ** cons_order_B_C_T2, uint64_t n_sequences);
/*
	Substracts the accumulated order offset from the blocks in given synteny lists
*/
void recompute_orders_from_offset(uint64_t * orders, uint64_t args_count, ...);
/*
	Returns TRUE if there is the same number of blocks per each genome, otherwise FALSE
*/
bool genomes_involved_in_synteny(uint64_t * genomes_counters, uint64_t n_sequences, uint64_t args_count, ...);

/*
	Applies event operations to a synteny list
*/
void apply_queue_operation(rearrangement * _r, Synteny_list * sl);
/*
	Returns the synteny level for any given number of synteny lists
	returns 0 if the lists dont share the synteny level
	otherwise returns the level of synteny
*/
uint64_t synteny_level_across_lists(uint64_t args_count, ...);

/*
	Remove DNA seq from an insertion 
*/
void remove_insertion_DNA(Block * a, Block * b, Block * c, uint64_t diff);

/*
	Remove DNA seq from a deletion
*/
void remove_deletion_DNA(Block * a, Block * b, Block * c, uint64_t diff);

/*
	Returns a list indicating which block breakpoints should be trimmed
	0 -> not present 
	1 -> no trimm 
	2,3 -> insertion / deletion
*/
void handle_indels(Synteny_list * A, Synteny_list * B, Synteny_list * C, uint64_t * indel_distance, uint64_t n_sequences, uint64_t * genomes_block_count, uint64_t * indel_kept, uint64_t * indel_type, events_queue * operations_queue, uint64_t * t_insertions, uint64_t * t_deletions);

/*
	Concatenates three synteny blocks into one
*/
void concat_synteny_blocks(Synteny_list ** A, Synteny_list ** B, Synteny_list ** C);


/*
	Computes the strand matrix for a synteny block to compute reversions
	The strand_matrix is assumed to have length n*n where n is the number of sequences
*/
void generate_strand_matrix(Synteny_block * sb, char ** strand_matrix);

/*
	Reverses a duplication and adds the operation to the event queue
*/
void reverse_duplication(Synteny_list * A, Synteny_list * B, Synteny_list * C, Block * dup, hash_table * ht, events_queue * operations_queue, uint64_t last_s_id);

/*
	Reverses a reversion by changing the sequence in place (the DNA) and the strand of fragments that involve the genomes
*/
void reverse_reversion(Synteny_list * B, sequence_manager * sm, bool * genome_ids_affected);

/*
	Reverses a transposition by displacing the affected blocks to the midpoint of the original syntenys
*/
bool reverse_tranposition(Synteny_list * A, Synteny_list * B, Synteny_list * C, Synteny_list * K1, Synteny_list * K2, Block ** blocks_to_move, Block ** blocks_to_stay, uint64_t n_sequences, events_queue * operations_queue);

/*
	Applies an operation queue for all blocks at once
*/
void apply_operation(Block * b, int64_t coordinates, int64_t order, uint64_t range1, uint64_t range2);
/*
	Detects candidates for evolutionary events
*/
void detect_evolutionary_event(Synteny_list * sbl, sequence_manager * seq_man, uint32_t kmer_size, hash_table * blocks_ht, uint64_t * last_s_id);

