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