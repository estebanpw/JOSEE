#include <stdio.h>
#include <cstdlib>
#include <sys/time.h>
#include <inttypes.h>
#include <cstdarg>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "structs.h"
#include "commonFunctions.h"
#include "comparisonFunctions.h"



void map_frags_to_genomes(unsigned char ** map_table, struct FragFile * frags, uint64_t n_frags, sequence_manager * seq_manager){

	uint64_t i, j, from, to, seq;
	//For all frags


	for(i=0;i<n_frags;i++){

		//printFragment(&frags[i]);
		seq = frags[i].seqX;
		from = frags[i].xStart;
		to = frags[i].xEnd;

		seq_manager->get_sequence_by_label(seq)->n_frags++;

		//Map coordinates in frag for seqX, which is always forward
		map_table[seq][from] = OPENFRAG;
		for(j=from+1;j<to;j++){
			if(map_table[seq][j] == NOFRAG){
				map_table[seq][j] = COVERFRAG;
			}
		}
		map_table[seq][to] = CLOSEFRAG;
		
		//printFragment(&frags[i]);
		//Map coordinates in frag for seqY, which might be reversed
		//Remember RAMGECKO coordinates are global respective to forward and with Ystart > Yend when reversed
		//So reversed should only be switched

		seq = frags[i].seqY;
		if(frags[i].strand == 'f'){
			from = frags[i].yStart;
			to = frags[i].yEnd;
		}else{			
			from = frags[i].yEnd;
			to = frags[i].yStart;	
		}

		seq_manager->get_sequence_by_label(seq)->n_frags++;

		map_table[seq][from] = OPENFRAG;
		for(j=from+1;j<to;j++){
			if(map_table[seq][j] == NOFRAG) map_table[seq][j] = COVERFRAG;
		}
		map_table[seq][to] = CLOSEFRAG;	

		

	}
	//At this point all coordinates have been mapped for the current fragments
	
}



inline void copyFragWithNewCoordinates(struct FragFile * destination, struct FragFile * source, uint64_t xStart, uint64_t yStart, uint64_t xEnd, uint64_t yEnd, uint64_t len){
    destination->xStart = xStart;
    destination->xEnd = xEnd;

	destination->yStart = yStart;
    destination->yEnd = yEnd;	

    destination->length = len;
    destination->seqX = source->seqX;
    destination->seqY = source->seqY;
    destination->strand = source->strand;
}

struct FragFile * trim_fragments_and_map(unsigned char ** map_table, struct FragFile * frags, uint64_t * n_frags, uint64_t min_len){

	//For debug only
	static int64_t itera = -1;
	itera++;
	uint64_t t_len = 0;

	struct FragFile new_frag;
	struct FragFile * list_new_frags;
	uint64_t list_reallocs = 1;
	uint64_t new_frags_number = 0;
	uint64_t size_fragment = sizeofFragment(); //To not compute it every time

	
	//Allocate memory
	list_new_frags = (struct FragFile *) std::malloc(INIT_TRIM_FRAGS*sizeofFragment());
	if(list_new_frags == NULL) terror("Could not allocate memory for list of new fragments in trimming");

	//Start trimming process
	uint64_t i, jX, jY, fromX, fromY, toX, toY, seqX, seqY;
	char strand;
	uint64_t cur_new_len;

	for(i=0; i<*n_frags; i++){
		
		
		//printf("Working on frag: %"PRIu64"\n", i);

		//Copy frag values
		fromX = frags[i].xStart; 
		toX = frags[i].xEnd + 1;  //Because coordinates are including [x,y]
		
		
		strand = frags[i].strand;
		fromY = frags[i].yStart;
		toY = frags[i].yEnd;
		/*
		if(strand == 'f'){
			toY += 1; 
		}else{
			if(toY > 0) toY -= 1;
		}
		*/

		jX = fromX;
		jY = fromY;	

		seqX = frags[i].seqX; seqY = frags[i].seqY;
		cur_new_len = 1;


		
		

		//(seqX == 0 && seqY == 1) &&
		while((jX <= toX && ( (strand == 'f' && jY <= toY) || (strand == 'r' && jY >= toY)))){
			//Check how long until there is a break (by starting or ending of frag)
			//Increase one to skip starting OPENFRAG
			
			jX++;
			cur_new_len++;
			if(strand == 'f'){ jY++; }else{ if(jY > 0) jY--; else break;} //To scape the buffer overflow of uints64
			
			//while((jX <= toX && ( (strand == 'f' && jY <= toY) || (strand == 'r' && jY >= toY))) && map_table[seqX][jX] == COVERFRAG && map_table[seqY][jY] == COVERFRAG){
			while((jX <= toX && ( (strand == 'f' && jY <= toY) || (strand == 'r' && jY >= toY))) && map_table[seqX][jX] == COVERFRAG && map_table[seqY][jY] == COVERFRAG){
				jX++;
				if(strand == 'f'){ jY++; }else{ if(jY > 0) jY--; else break;} //To scape the buffer overflow of uints64
				cur_new_len++;
			}
			

			//At this point, jX and jY hold the ending coordinates of the new fragment
			//And fromX and fromY hold the starting coordinates

			/*
			printf("jX: %"PRIu64"; from: %"PRIu64", to: %"PRIu64"\n", jX, fromX, toX);
			printf("jY: %"PRIu64"; from: %"PRIu64", to: %"PRIu64"\n", jY, fromY, toY);
			getchar();
			*/

			
			
			
			/*
			if(itera == 8){
				printFragment(&frags[i]);			
				if(strand == 'r'){
					printf("(%u - %u) %c@ [%"PRIu64", %"PRIu64"] up to [%"PRIu64", %"PRIu64"] aprox: (%"PRIu64") Frag: %"PRIu64"\n", map_table[seqX][jX], map_table[seqY][jY], strand, jX, jY, toX, toY, jX-fromX, frags[i].diag);
					printf("BEFORE\n");
					print_maptable_portion(map_table, fromX, toX+2, 50, seqX);
					print_maptable_portion(map_table, toY, fromY+2, 50, seqY);
				}else{
					printf("(%u - %u) %c@ [%"PRIu64", %"PRIu64"] up to [%"PRIu64", %"PRIu64"] aprox: (%"PRIu64") Frag: %"PRIu64"\n", map_table[seqX][jX], map_table[seqY][jY], strand, jX, jY, toX, toY, jX-fromX, frags[i].diag);
					printf("BEFORE\n");
					print_maptable_portion(map_table, fromX, toX+2, 50, seqX);
					print_maptable_portion(map_table, fromY, toY+2, 50, seqY);
				}
			}
			*/
			

			map_table[seqX][jX] = CLOSEFRAG;
			if(strand == 'f') map_table[seqY][jY] = CLOSEFRAG; else map_table[seqY][jY] = CLOSEFRAG;
			map_table[seqX][fromX] = OPENFRAG;
			if(strand == 'f') map_table[seqY][fromY] = OPENFRAG; else map_table[seqY][fromY] = OPENFRAG;

			

			/*
			if(itera == 8){
				if(strand == 'r'){
					printf("AFter\n");
					print_maptable_portion(map_table, fromX, toX+2, 50, seqX);
					print_maptable_portion(map_table, toY, fromY+2, 50, seqY);
					getchar();
				}
				
				else{
					printf("AFter\n");
					print_maptable_portion(map_table, fromX, toX+2, 50, seqX);
					print_maptable_portion(map_table, fromY, toY+2, 50, seqY);
					getchar();
				}
			}
			*/
			
			if(cur_new_len >= min_len){ //Filtering

				//DEBUG
				t_len += cur_new_len;
				
				
				//The fragment must be snipped out and saved
				copyFragWithNewCoordinates(&new_frag, &frags[i], fromX, fromY, jX, jY, cur_new_len);
				memcpy(&list_new_frags[new_frags_number], &new_frag, size_fragment);
				/*
				if(frags[i].length < cur_new_len){
					printf("See this???\n"); 
					printFragment(&list_new_frags[new_frags_number]);
					printFragment(&frags[i]);
					printf("It was supposed to be: \n");
					printf("jX: %"PRIu64"; from: %"PRIu64", to: %"PRIu64"\n", jX, fromX, toX);
					printf("jY: %"PRIu64"; from: %"PRIu64", to: %"PRIu64"\n", jY, fromY, toY);
					getchar();

				}
				*/

				new_frags_number++;
				//Check if we need to realloc the list of new frags
				if(new_frags_number == list_reallocs*INIT_TRIM_FRAGS){
					list_reallocs++;
					list_new_frags = (struct FragFile *) std::realloc(list_new_frags, list_reallocs*INIT_TRIM_FRAGS*size_fragment);
					if(list_new_frags == NULL) terror("Could not realloc fragments on the trimming process");
				}

			}
			

			//If you are here, either the fragment was too short, or was written correctly or we are at the end of the frag
			//Just keep going
			//Copy frag values
			fromX = jX;
			fromY = jY;

			cur_new_len = 1;

			//End of outside while
		}


	}

	//printf("f avg: %"PRIu64" t_len: %"PRIu64"\n", t_len/new_frags_number, t_len);

	*n_frags = new_frags_number;

	return list_new_frags;

}


void compute_order_of_blocks(hash_table * ht, uint64_t n_seqs){
	uint64_t i;
	Bucket * ptr;
	uint64_t * seq_orders = (uint64_t *) std::calloc(n_seqs, sizeof(uint64_t));
	if(seq_orders == NULL) terror("Could not allocate vector of orders");

	//Compute orders
	for(i=0;i<ht->get_size();i++){
		ptr = ht->get_key_at(i);
		while(ptr != NULL){

			ptr->b.order = seq_orders[ptr->b.genome->id];
			seq_orders[ptr->b.genome->id]++;

			ptr = ptr->next;
		}
	}
	//Not needed anymore
	std::free(seq_orders);
	
}

Synteny_list * compute_synteny_list(hash_table * ht, uint64_t n_seqs, memory_pool * mp){
	uint64_t i, synteny_level;
	//Pointers
	Bucket * ptr;
	Block * aux_block = NULL;
	Frags_list * flptr;
	//Avoid overcomputation
	uint64_t pre_comp_sb = sizeofSyntenyBlock();
	uint64_t pre_comp_sbl = sizeofSyntenyList();

	//Bit mask to tell if the synteny block already contains a genome
	unsigned char * had_genome_bitmask = (unsigned char *) std::malloc(n_seqs*sizeof(unsigned char));
	if(had_genome_bitmask == NULL) terror("Could not allocate bit mask");

	Synteny_list * sbl = (Synteny_list *) mp->request_bytes(pre_comp_sbl);
	Synteny_list * curr_sbl = sbl;
	Synteny_list * last_sbl = NULL;
	Synteny_block * curr_sb = NULL;


	for(i=0;i<ht->get_size();i++){
		ptr = ht->get_key_at(i);

		while(ptr != NULL){
			//Each block in the has is here
			//For each block, add the blocks linked by the fragments
			memset(had_genome_bitmask, 0, n_seqs); //Reset genome counters
			synteny_level = 0; //Restart synteny level
			flptr = ptr->b.f_list;
			while(flptr != NULL){
				//printFragment(flptr->f);
				//Each fragment in the current block
				//Only if we did not already have the genome in the frags
				aux_block = ht->get_block_from_frag(flptr->f, 0);

				/*
				if(aux_block != NULL){
					printf("FETCH:\n"); printBlock(aux_block); getchar();
				}
				*/

				//Insert frag_x
				if(aux_block != NULL && aux_block->present_in_synteny == 0){
					aux_block->present_in_synteny = 1;
					Synteny_block * aux_sb = (Synteny_block *) mp->request_bytes(pre_comp_sb);
					aux_sb->b = aux_block; //insert at the head
					aux_sb->next = curr_sb;
					curr_sb = aux_sb;
					aux_block = NULL;
					if(had_genome_bitmask[flptr->f->seqX] == 0) synteny_level++;
					had_genome_bitmask[flptr->f->seqX] = 1;
					
					//printf("\t"); printBlock(aux_sb->b);
				}

				aux_block = ht->get_block_from_frag(flptr->f, 1);

				/*
				if(aux_block != NULL){
					printf("FETCH:\n"); printBlock(aux_block); getchar();
				}
				*/

				//Insert frag_y
				if(aux_block != NULL && aux_block->present_in_synteny == 0){
					aux_block->present_in_synteny = 1;
					Synteny_block * aux_sb = (Synteny_block *) mp->request_bytes(pre_comp_sb);
					aux_sb->b = aux_block; //insert at the head
					aux_sb->next = curr_sb;
					curr_sb = aux_sb;
					aux_block = NULL;
					if(had_genome_bitmask[flptr->f->seqY] == 0) synteny_level++;
					had_genome_bitmask[flptr->f->seqY] = 1;
					
					//printf("\t"); printBlock(aux_sb->b);
				}
				
				flptr = flptr->next;
			}
			//printf("broke stnyteny ---------------------------\n");

			//No more frags to add 

			// End synteny block
			if(synteny_level > 1){
				curr_sbl->sb = curr_sb;
				curr_sbl->synteny_level = synteny_level;
				//curr_sbl->prev = last_sbl;
				
				curr_sbl->next = (Synteny_list *) mp->request_bytes(pre_comp_sbl);
				curr_sbl = curr_sbl->next;
				curr_sbl->next = NULL;
				curr_sbl->prev = last_sbl;
				last_sbl = curr_sbl;
				
				//printf("Generated\n");
			}
			if(synteny_level <= 1){
				//Since there is a minimum synteny, restore its level so that it can be used
				//curr_sb->b->present_in_synteny = 0;

				//Restore levels of all of those used
				Synteny_block * rest_ptr = curr_sbl->sb;
				while(rest_ptr != NULL){
					rest_ptr->b->present_in_synteny = 0;
					rest_ptr = rest_ptr->next;
				}
				mp->reset_n_bytes(pre_comp_sb);
				
				//printf("Failed at generatin\n");
			}
			
			curr_sb = NULL;
			//Go to next block
			ptr = ptr->next;
			
			
			//printf("broke stnyteny ---------------------------\n");
			//getchar();
		}
		//printf("broke stnyteny ---------------------------\n");
	}

	curr_sbl = NULL;
	sbl->next->prev = sbl;

	std::free(had_genome_bitmask);
	return sbl;
}

//@Assumes: Synteny level + same number of genomes involved
bool consecutive_block_order(uint64_t * pairs_diff, uint64_t args_count, ...){
	va_list sbl_args;
	va_start(sbl_args, args_count);
	Synteny_list * sl_ptr;
	Synteny_block * sb_ptr;
	
	uint64_t i;
	for(i=0;i<args_count;i++){
		sl_ptr = va_arg(sbl_args, Synteny_list *);
		if(sl_ptr != NULL){
			sb_ptr = sl_ptr->sb;
			while(sb_ptr != NULL){

				if(i == 0){
					pairs_diff[sb_ptr->b->genome->id] = sb_ptr->b->order;
				}else{
					if(sb_ptr->b->order - pairs_diff[sb_ptr->b->genome->id] != 1) return false;
					pairs_diff[sb_ptr->b->genome->id] = sb_ptr->b->order;
				}
				sb_ptr = sb_ptr->next;
			}
		}
	}

	return true;
}

void recompute_orders_from_offset(uint64_t * orders, uint64_t args_count, ...){
	va_list sbl_args;
	va_start(sbl_args, args_count);
	Synteny_list * sl_ptr;
	Synteny_block * sb_ptr;
	
	uint64_t i;
	for(i=0;i<args_count;i++){
		sl_ptr = va_arg(sbl_args, Synteny_list *);
		if(sl_ptr != NULL){
			sb_ptr = sl_ptr->sb;
			while(sb_ptr != NULL){
				
				if(sb_ptr->b->order < orders[sb_ptr->b->genome->id]){
					printf("Happening at\n");
					printBlock(sb_ptr->b);
					printf("I want to substract %"PRIu64"\n", orders[sb_ptr->b->genome->id]);
					getchar();
				}
				sb_ptr->b->order = sb_ptr->b->order - orders[sb_ptr->b->genome->id];
				sb_ptr = sb_ptr->next;
			}
		}
	}
	va_end(sbl_args);
}


uint64_t synteny_level_across_lists(uint64_t args_count, ...){
	va_list sbl_args;
	va_start(sbl_args, args_count);
	Synteny_list * sl_ptr;
	sl_ptr = va_arg(sbl_args, Synteny_list *);
	uint64_t s_level;
	if(sl_ptr != NULL) s_level = sl_ptr->synteny_level; else return 0;
	

	uint64_t i;
	for(i=1;i<args_count;i++){
		sl_ptr = va_arg(sbl_args, Synteny_list *);
		if(sl_ptr == NULL || sl_ptr->synteny_level != s_level) return 0;
	}
	va_end(sbl_args);
	return s_level;
}

bool genomes_involved_in_synteny(uint64_t * genomes_counters, uint64_t n_sequences, uint64_t args_count, ...){
	va_list sbl_args;
	va_start(sbl_args, args_count);
	Synteny_list * sl_ptr;
	Synteny_block * sb_ptr;
		
	uint64_t i;
	for(i=0;i<args_count;i++){
		sl_ptr = va_arg(sbl_args, Synteny_list *);
		sb_ptr = sl_ptr->sb;
		
		while(sb_ptr != NULL){
			genomes_counters[sb_ptr->b->genome->id]++;
			sb_ptr = sb_ptr->next;
		}

		
	}
	va_end(sbl_args);

	//printInvolvedGenomes(genomes_counters, n_sequences);

	uint64_t first = 0;
	for(i=0;i<n_sequences;i++){

		if(first == 0 && genomes_counters[i] != 0){
			first = genomes_counters[i];
		}else{
			if(genomes_counters[i] != 0 && first != genomes_counters[i]) return false;
		}
		
	}

	return true;
}

void concat_synteny_blocks(Synteny_list * A, Synteny_list * B, Synteny_list * C){
	//printf("I would like to concat\n");
	//printf("And it would look like this:\n");
	

	uint64_t i;
	Synteny_block * start_sb_ptr = A->sb;
	Synteny_block * mid_ptr = B->sb;
	Synteny_block * end_sb_ptr = C->sb;
	Frags_list * fl_A, * fl_B;

	for(i=0;i<A->synteny_level;i++){
		start_sb_ptr->b->end = end_sb_ptr->b->end;
		//Append frags lists
		//Find last pointer in A
		fl_A = start_sb_ptr->b->f_list;
		while(fl_A->next != NULL) fl_A = fl_A->next;

		//Find last pointer in B
		fl_B = mid_ptr->b->f_list;
		while(fl_B->next != NULL) fl_B = fl_B->next;

		//Append C to B, and B to A
		fl_B->next = end_sb_ptr->b->f_list;
		fl_A->next = mid_ptr->b->f_list;

		start_sb_ptr = start_sb_ptr->next;
		mid_ptr = mid_ptr->next;
		end_sb_ptr = end_sb_ptr->next;
	}
	
	//printSyntenyBlock(A->sb);
	//getchar();

	//Remove intermediate synteny block list
	A->next = C->next;
	//B cant be accessed now. Its not dangling because of the mempool.

	//All synteny lists should be updated from 

}

void detect_evolutionary_event(Synteny_list * sbl, sequence_manager * seq_man, uint32_t kmer_size){
	
	//Data structures needed

	uint64_t i;
	uint64_t n_sequences = seq_man->get_number_of_sequences();

	//Until nothing else cant be done, keep iterating
	bool stop_criteria = false; 
	bool update_pointers_after_concat = false;

	//For telling whether the same number of blocks per genome is involved in synteny
	uint64_t * genomes_block_count = (uint64_t *) std::malloc(n_sequences*sizeof(uint64_t));
	if(genomes_block_count == NULL) terror("Could not allocate genome blocks counter");

	//To recompute order for the next blocks after an event
	uint64_t * order_offsets = (uint64_t *) std::malloc(n_sequences*sizeof(uint64_t));
	uint64_t * order_offsets_after_concat = (uint64_t *) std::malloc(n_sequences*sizeof(uint64_t));
	if(order_offsets == NULL || order_offsets_after_concat == NULL) terror("Could not allocate order offsets for after-events");

	//To check that blocks are consecutive in their genome
	uint64_t * pairs_diff = (uint64_t *) std::malloc(n_sequences*sizeof(uint64_t));
	if(pairs_diff == NULL) terror("Could not allocate consecutive order of blocks array");

	//For strand matrices
	strand_matrix * sm_A, * sm_B, * sm_C, * sm_D, * sm_E;
	unsigned char ** _tmp1, ** _tmp2;
	sm_A = new strand_matrix(n_sequences);
	sm_B = new strand_matrix(n_sequences);
	sm_C = new strand_matrix(n_sequences);
	sm_D = new strand_matrix(n_sequences);
	sm_E = new strand_matrix(n_sequences);

	//For hits and frags computation
	dictionary_hash * words_dictionary = new dictionary_hash(seq_man->get_maximum_length()/TABLE_RATE, seq_man->get_maximum_length(), kmer_size);
	Quickfrag ** qfmat = (Quickfrag **) std::malloc(n_sequences*n_sequences*sizeof(Quickfrag *));
	double ** qf_submat = (double **) std::malloc(n_sequences*n_sequences*sizeof(double *));
	unsigned char ** qfmat_state = (unsigned char **) std::malloc(n_sequences*n_sequences*sizeof(unsigned char *));
	if(qfmat == NULL || qfmat_state == NULL || qf_submat == NULL) terror("Could not allocate pairwise alignment matrix (1)");
	
	for(i=0;i<n_sequences;i++){
		order_offsets_after_concat[i] = 2; //Always to use when recomputing after concat
		qfmat[i] = (Quickfrag *) std::malloc(n_sequences*sizeofQuickfrag());
		qf_submat[i] = (double *) std::malloc(n_sequences*sizeof(double));
		qfmat_state[i] = (unsigned char *) std::malloc(n_sequences*sizeof(unsigned char));
		if(qfmat[i] == NULL || qfmat_state[i] == NULL || qf_submat == NULL) terror("Could not allocate pairwsie alignment matrix (2)");
	}

	//For clustering
	memory_pool * mp = new memory_pool(1, POOL_SIZE);
	/*
	memory_pool * mp = new memory_pool(1,
	 seq_man->get_number_of_sequences()*sizeofSlist()*2 +
	  2*seq_man->get_number_of_sequences()*sizeof(unsigned char) +
	  seq_man->get_number_of_sequences()*sizeof(Slist *));
	*/

	


	//Lists of synteny blocks to address evolutionary events
	Synteny_list * A, * B = NULL, * C = NULL, * D = NULL, * E = NULL;

	//To have some statistics
	uint64_t current_step = 0, current_concats = 0, t_concats = 0;
	uint64_t t_inversions = 0, t_duplications = 0;


	while(!stop_criteria){
		
		//Display current iteration
		printf("\nOn step: %"PRIu64". Total concats: %"PRIu64", this round: %"PRIu64"\n", current_step++, t_concats, current_concats);
		printf("Total inversions: %"PRIu64". Total duplications: %"PRIu64"\n", t_inversions, t_duplications);
		current_concats = 0;
		getchar();

		//In case nothing gets done, stop iterating
		stop_criteria = true;

		//Get head of synteny list
		A = sbl; 

		//Set offset orders to zero
		memset(order_offsets, 0, n_sequences*sizeof(int64_t));

		//Copy pointers of first consecutive blocks
		if(A != NULL) B = A->next; else return;
		if(B != NULL) C = B->next; else return; //Three at least
		if(C != NULL) D = C->next;
		if(D != NULL) E = D->next;


		//Generate their strand matrices
		sm_A->add_fragment_strands(A);
		sm_B->add_fragment_strands(B);
		sm_C->add_fragment_strands(C);
		sm_D->add_fragment_strands(D);
		sm_E->add_fragment_strands(E);

		while(A != NULL && B != NULL && C != NULL){ // AT least three



			//Recompute order of the last one added to the list
			printf("Pre all\n");
			recompute_orders_from_offset(order_offsets, 1, E);
			
			//printDebugBlockOrderByGenome(E, 0);
			
			
			printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
			printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
			if(A != NULL){ printSyntenyBlock(A->sb); printf("=was A=======000000\n");}
			if(B != NULL){ printSyntenyBlock(B->sb); printf("=was B=======000000\n");}
			if(C != NULL){ printSyntenyBlock(C->sb); printf("=was C=======000000\n");}
			if(D != NULL){ printSyntenyBlock(D->sb); printf("=was D=======000000\n");}
			

			// Transpositions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			// I think only 4 synteny groups are needed ??
			if(D != NULL){
				if(synteny_level_across_lists(3, A, B, D)){
					//same synteny, check the number of genomes involved
					memset(genomes_block_count, 0, n_sequences*sizeof(uint64_t));
					if(genomes_involved_in_synteny(genomes_block_count, n_sequences, 3, A, B, D)){
						//Check that A and B have their order consecutive
						if(consecutive_block_order(pairs_diff, 2, A, B)){
							//And if D is not consecutive with B, then there must exist C with insertions
							if(!consecutive_block_order(pairs_diff, 2, B, D)){
								//But only if C does not have as much synteny level as D
								if(C->synteny_level < D->synteny_level){
									// TODO I think I have to check orders more carefully
									printf("Detected insertions\n");
									getchar();
								}
							}
						}
					}
				}
			}
			


			// Duplications %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			memset(genomes_block_count, 0, n_sequences*sizeof(uint64_t));
			//If there is not the same number of genomes involved
			if(!genomes_involved_in_synteny(genomes_block_count, n_sequences, 1, B)){
				//There are duplications in B
				printf("Stopping it\n");
				//Find those that have more synteny level
				for(i=0;i<n_sequences;i++){
					if(genomes_block_count[i] > 1){
						// Genome i has duplications
						printf("Duplications in %"PRIu64"\n", i);
						t_duplications++;
					}
				}
				//getchar();
				
			}
			
			// Inversions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if(synteny_level_across_lists(3, A, B, C) > 0){
				if(sm_A->get_strands_type() != MIXED && 
				sm_A->get_strands_type() == sm_C->get_strands_type() &&
				sm_B->get_strands_type() == MIXED){

					memset(genomes_block_count, 0, n_sequences*sizeof(uint64_t));
					if(genomes_involved_in_synteny(genomes_block_count, n_sequences, 3, A, B, C)){
						if(consecutive_block_order(pairs_diff, 3, A, B, C)){
							printf("Attention: this looks like a reversion\n");

							read_words_from_synteny_block_and_align(seq_man, B, kmer_size, words_dictionary, qfmat, qfmat_state);
							mp->reset_to(0,0);
							UPGMA_joining_clustering(qfmat, qf_submat, qfmat_state, seq_man->get_number_of_sequences(), mp);
							//getchar();
							//What do here?
							t_inversions++;
						}
					}
				}else{
					printf("Strands not qualifying for inversion\n");
				}
			}else{
				printf("Different synteny for inversion\n");
			}


			// Concatenation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			// @IMPORTANT: Currently blocks are not being removed or recomputed
			// in hash table after concatenation

			//Check they share the synteny level
			if(synteny_level_across_lists(3, A, B, C) > 0){
				//Concat synteny blocks if they have the same strand
				if(sm_A->get_strands_type() != MIXED && 
				sm_A->get_strands_type() == sm_B->get_strands_type() &&
				sm_B->get_strands_type() == sm_C->get_strands_type()){

					//Erase genome counter
					memset(genomes_block_count, 0, n_sequences*sizeof(uint64_t));
					//for(i=0;i<n_sequences;i++) genomes_block_count[i] = 0;

					//Check that there is the same number of blocks per genome involved
					if(genomes_involved_in_synteny(genomes_block_count, n_sequences, 3, A, B, C)){

						if(consecutive_block_order(pairs_diff, 3, A, B, C)){
							concat_synteny_blocks(A, B, C);
							t_concats++;
							current_concats++;
							//getchar();
							//Add offset to orders
							for(i=0;i<n_sequences;i++){
								//If genome was involved we have to add offset to the concat
								if(genomes_block_count[i] != 0) order_offsets[i] += 2; //Two because two blocks are shrinked into one
								printf("G: %"PRIu64" -> O:%"PRIu64"\n", i, order_offsets[i]);
							}
							//Update current pointers
							update_pointers_after_concat = true;
							stop_criteria = false;
						}//else{
						//	printf("Non consecutive order in blocks...\n");
						//}
					}//else{
					//	printf("Genomes involved different number ...\n"); //getchar();
					//}
				}//else{
				//	printf("Frags differ in strand...\n"); //getchar();
				//}
					
			}//else{
			//	printf("Different synteny levels...\n"); //getchar();
			//}
			

			/*
			// Piece of code for pairwise alignment between blocks
			read_words_from_synteny_block_and_align(seq_man, A, kmer_size, words_dictionary, qfmat, qfmat_state);
			mp->reset_to(0,0);
			UPGMA_joining_clustering(qfmat, qf_submat, qfmat_state, seq_man->get_number_of_sequences(), mp);
			*/
			
			if(update_pointers_after_concat){
				
				sm_A->reset();
				sm_A->add_fragment_strands(A);
				_tmp1 = sm_B->sm; //Do not lose B sm pointer
				_tmp2 = sm_C->sm; //Same for C
				sm_B->sm = sm_D->sm;
				sm_C->sm = sm_E->sm;
				sm_D->sm = _tmp1; // Acquire
				sm_E->sm = _tmp2; // Dont lose it


				//A is still A
				B = D;
				C = E;
				if(C != NULL) D = E->next; else D = NULL;
				if(D != NULL) E = D->next; else E = NULL;

				if(D != NULL){ sm_D->reset(); sm_D->add_fragment_strands(D);}
				if(E != NULL){ sm_E->reset(); sm_E->add_fragment_strands(E);}

				printf("Recomputing normals\n");
				recompute_orders_from_offset(order_offsets_after_concat, 3, B, C, D); //E is recomputed at beginning

				update_pointers_after_concat = false;

			}else{
				//Only generate the new strand matrix and pass the others
				_tmp1 = sm_A->sm; // Do not lose pointer to strand matrix
				sm_A->sm = sm_B->sm;
				sm_B->sm = sm_C->sm;
				sm_C->sm = sm_D->sm;
				sm_D->sm = sm_E->sm;
				sm_E->sm = _tmp1; // Recover strand matrix
				
				//advance pointers
				A = B;
				B = C;	
				C = D;
				D = E;

				//next iteration
				if(E != NULL) E = E->next;

				//And generate new strand matrix
				sm_E->reset(); //Need hard reset to not use fragments from last synteny
				sm_E->add_fragment_strands(E);
			}
		}
	}

	printf("\nAfter %"PRIu64" step(s). Total concats: %"PRIu64", this round: %"PRIu64"\n", current_step++, t_concats, current_concats);
	printf("Total inversions: %"PRIu64". Total duplications: %"PRIu64"\n", t_inversions, t_duplications);

	for(i=0;i<seq_man->get_number_of_sequences();i++){
		std::free(qfmat[i]);
		std::free(qfmat_state[i]);
	}
	std::free(qf_submat);
	std::free(qfmat);
	std::free(qfmat_state);
	std::free(genomes_block_count);
	std::free(order_offsets);
	std::free(order_offsets_after_concat);
	std::free(pairs_diff);
	
	delete sm_A;
	delete sm_B;
	delete sm_C;
	delete sm_D;
	delete sm_E;

	delete mp;
	delete words_dictionary;
}

