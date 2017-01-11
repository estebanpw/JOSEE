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
void map_frags_to_genomes(unsigned char ** map_table, struct FragFile * frags, uint64_t n_frags);

/*
	Copy frags properties to new snipped fragment
*/
inline void copyFragWithNewCoordinates(struct FragFile * destination, struct FragFile * source, uint64_t xStart, uint64_t yStart, uint64_t xEnd, uint64_t yEnd, uint64_t len);


/*
	This function cuts frags into more if there are other frags that overlap partially
	min_len acts as a filter to remove short fragments
*/
struct FragFile * trim_fragments_and_map(unsigned char ** map_table, struct FragFile * frags, uint64_t * n_frags, uint64_t min_len, Sequence * sequences);

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
	Checks if a truple of consecutive blocks is a reversion or not
*/
void has_reversion_in_truple(Bucket * a, Bucket * b, Bucket * c);

/*
	Detects if a reversion took place
*/
void detect_evolutionary_event(hash_table * ht, uint64_t max_len_seq);