#include <stdio.h>
#include <cstdlib>
#include <sys/time.h>
#include <inttypes.h>
#include <cstdarg>
#include <iostream>
#include <stack>
#include <ctype.h>
#include <string.h>
#include <math.h>
#include "structs.h"
#include "commonFunctions.h"
#include "comparisonFunctions.h"
#include "alignment_functions.h"



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
			ptr->b.id = ptr->b.order; // Order will be modified when applying events, but the id will not
			seq_orders[ptr->b.genome->id]++;

			if(ht->last_blocks[ptr->b.genome->id] != NULL){
				ht->last_blocks[ptr->b.genome->id]->next = &ptr->b;
			}
			ptr->b.prev = ht->last_blocks[ptr->b.genome->id];
			ht->last_blocks[ptr->b.genome->id] = &ptr->b;

			ptr = ptr->next;
		}
	}
	ht->release_temp_last_blocks();
	//Not needed anymore
	std::free(seq_orders);
	
}


Synteny_list * compute_synteny_list(hash_table * ht, uint64_t n_seqs, memory_pool * mp, uint64_t * last_s_id){
	uint64_t i, synteny_level, curr_id = 0;
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

	//List to add pointers to the blocks whose fragments will have to be added
	std::stack<Block *> * blocks_to_add = new std::stack<Block *>();

	Synteny_list * sbl = (Synteny_list *) mp->request_bytes(pre_comp_sbl);
	Synteny_list * curr_sbl = sbl;
	Synteny_list * last_sbl = NULL;
	Synteny_block * curr_sb = NULL, * find_pos_in_sb = NULL, * last_to_insert = NULL;
	Block * curr_block;


	for(i=0;i<ht->get_size();i++){
		ptr = ht->get_key_at(i);

		
		curr_sb = NULL;
		while(ptr != NULL){
			//Each block is here
			
			synteny_level = 0; //Restart synteny level
			memset(had_genome_bitmask, 0, n_seqs); //Reset genome counters

			//Add current block
			blocks_to_add->push(&ptr->b);
			//For each block, add the blocks linked by the fragments
			while(blocks_to_add->size() > 0){
				//Get next block added to the stack and remove it
				curr_block = blocks_to_add->top();
				blocks_to_add->pop();			
			

				//Check its fragments
				flptr = curr_block->f_list;
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
					// curr_sb is always the head

					//Insert frag_x
					if(aux_block != NULL && aux_block->present_in_synteny == NULL){
						//Add block so that we check its fragments
						blocks_to_add->push(aux_block);
						aux_block->present_in_synteny = curr_sbl;
						Synteny_block * aux_sb = (Synteny_block *) mp->request_bytes(pre_comp_sb);
						aux_sb->b = aux_block; // Insert in order by genomes (later use!!)
						
						if(curr_sb == NULL){
							//Head is null
							curr_sb = aux_sb;
						}else{
							//Finds prev and next
							find_pos_in_sb = curr_sb;
							last_to_insert = NULL;
							
							while(find_pos_in_sb != NULL && aux_sb->b->genome->id > find_pos_in_sb->b->genome->id){
								last_to_insert = find_pos_in_sb;
								find_pos_in_sb = find_pos_in_sb->next;
							}

							aux_sb->next = find_pos_in_sb;
							if(last_to_insert == NULL){
								//its new head
								curr_sb = aux_sb;
							}else{
								last_to_insert->next = aux_sb;
							}
							
						}
						
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
					if(aux_block != NULL && aux_block->present_in_synteny == NULL){
						//Add block so that we check its fragments
						blocks_to_add->push(aux_block);
						aux_block->present_in_synteny = curr_sbl;
						Synteny_block * aux_sb = (Synteny_block *) mp->request_bytes(pre_comp_sb);
						aux_sb->b = aux_block; // Insert in order by genomes (later use!!)

						if(curr_sb == NULL){
							//Head is null
							curr_sb = aux_sb;
						}else{

							//Finds prev and next
							find_pos_in_sb = curr_sb;
							last_to_insert = NULL;
							
							while(find_pos_in_sb != NULL && aux_sb->b->genome->id > find_pos_in_sb->b->genome->id){
								last_to_insert = find_pos_in_sb;
								find_pos_in_sb = find_pos_in_sb->next;
							}

							aux_sb->next = find_pos_in_sb;
							if(last_to_insert == NULL){
								//its new head
								curr_sb = aux_sb;
							}else{
								last_to_insert->next = aux_sb;
							}
							
						}


						aux_block = NULL;
						if(had_genome_bitmask[flptr->f->seqY] == 0) synteny_level++;
						had_genome_bitmask[flptr->f->seqY] = 1;
						
						//printf("\t"); printBlock(aux_sb->b);
					}
					
					flptr = flptr->next;
				}
			}
			
			//printf("broke stnyteny ---------------------------\n");

			// End synteny block
			if(synteny_level > 1){
				curr_sbl->sb = curr_sb;
				curr_sbl->synteny_level = synteny_level;
				curr_sbl->id = curr_id;
				*last_s_id = curr_id;
				curr_id++;
				//curr_sbl->prev = last_sbl;
				
				curr_sbl->next = (Synteny_list *) mp->request_bytes(pre_comp_sbl);
				curr_sbl = curr_sbl->next;
				curr_sbl->next = NULL;
				curr_sbl->prev = last_sbl;
				last_sbl = curr_sbl;
				
				//printf("Generated\n");
			}else{
				//Since there is a minimum synteny, restore its level so that it can be used
				//curr_sb->b->present_in_synteny = 0;

				//Restore levels of all of those used
				Synteny_block * rest_ptr = curr_sbl->sb;
				while(rest_ptr != NULL){
					rest_ptr->b->present_in_synteny = NULL;
					rest_ptr = rest_ptr->next;
				}
				mp->reset_n_bytes(pre_comp_sb);
				
				//printf("Failed at generatin\n");
			}
			
			curr_sb = NULL;


			//No more frags to add 
			//Go to next block
			ptr = ptr->next;
			
			
			//printf("broke stnyteny ---------------------------\n");
			//getchar();
		}
		

		//printf("broke stnyteny ---------------------------\n");
	}

	curr_sbl = curr_sbl->prev;
	if(curr_sbl != NULL) curr_sbl->next = NULL;
	//sbl->next->prev = sbl;

	std::free(had_genome_bitmask);
	delete blocks_to_add;

	return sbl;
}



// @Assumes same synteny level between lists
// and same number of genomes involved
void distance_between_blocks(uint64_t * distances, Synteny_list * A, Synteny_list * B){
	Synteny_block * sb_A = A->sb, * sb_B = B->sb;

	while(sb_A != NULL){

		distances[sb_A->b->genome->id] = sb_B->b->start - sb_A->b->end; 

		sb_A = sb_A->next;
		sb_B = sb_B->next;
	}
}

// @Assumes: Synteny level + same number of genomes involved
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

// At this stage of development it should only be used with a pair of syntenys
bool consecutive_block_order_except_one(int64_t * pairs_diff, uint64_t n_sequences, Block ** cons_order_T1, Block ** cons_order_T2, uint64_t args_count, ...){
	va_list sbl_args;
	va_start(sbl_args, args_count);
	Synteny_list * sl_ptr;
	Synteny_block * sb_ptr;
	uint64_t i;
	int64_t diff_type_1 = 0, diff_type_2 = 0, current_diff;
	bool is_in_T1_or_T2 = false;

	for(i=0;i<args_count;i++){
		sl_ptr = va_arg(sbl_args, Synteny_list *);
		if(sl_ptr != NULL){
			sb_ptr = sl_ptr->sb;
			while(sb_ptr != NULL){

				if(i == 0){
					//This inserts the orders from A
					pairs_diff[sb_ptr->b->genome->id] = sb_ptr->b->order;
				}else{
					//Now check for different orders
					is_in_T1_or_T2 = false;
					current_diff = sb_ptr->b->order - pairs_diff[sb_ptr->b->genome->id];


					if(diff_type_1 == 0){
						//Set the difference type
						diff_type_1 = current_diff;
						//Add the block to the difference cluster type
						cons_order_T1[sb_ptr->b->genome->id] = sb_ptr->b;
						is_in_T1_or_T2 = true;
					}else{
						if(current_diff == diff_type_1){
							cons_order_T1[sb_ptr->b->genome->id] = sb_ptr->b;
							is_in_T1_or_T2 = true;
						}
					}
					if(diff_type_2 == 0 && current_diff != diff_type_1){
						diff_type_2 = current_diff;
						cons_order_T2[sb_ptr->b->genome->id] = sb_ptr->b;
						is_in_T1_or_T2 = true;
					}else{
						if(current_diff == diff_type_2){
							cons_order_T2[sb_ptr->b->genome->id] = sb_ptr->b;
							is_in_T1_or_T2 = true;
						}
					}

					if(is_in_T1_or_T2 == false) return false; //There are more difference types
					
					
				}
				sb_ptr = sb_ptr->next;
			}
		}
	}
	#ifdef VERBOSE
	printf("DIFFTYPES: %"PRId64", %"PRId64"\n", diff_type_1, diff_type_2);
	#endif
	/*
	for(i=0;i<n_sequences;i++){
		if(cons_order_T1[i] != NULL) {printf("in T1@%"PRIu64": ",i); printBlockWriteMode(cons_order_T1[i]);}
		if(cons_order_T2[i] != NULL) {printf("in T2@%"PRIu64": ",i); printBlockWriteMode(cons_order_T2[i]);}
	}
	*/
	if(diff_type_1 != 1 && diff_type_2 != 1) return false; //Avoid false positives
	if(diff_type_1 == 1){
		//Switch them so that const_order_T1 has the one with higher diff
		Block ** aux_list = cons_order_T1;
		cons_order_T1 = cons_order_T2;
		cons_order_T2 = aux_list;
	}
	/*
	for(i=0;i<n_sequences;i++){
		if(cons_order_T1[i] != NULL) {printf("in T1@%"PRIu64": ",i); printBlockWriteMode(cons_order_T1[i]);}
		if(cons_order_T2[i] != NULL) {printf("in T2@%"PRIu64": ",i); printBlockWriteMode(cons_order_T2[i]);}
	}
	*/
	return true;
}

Block * compare_order_clusters(Block ** cons_order_A_B_T1, Block ** cons_order_A_B_T2, Block ** cons_order_B_C_T1, Block ** cons_order_B_C_T2, uint64_t n_sequences){
	
	Block * the_pointer_to_retrieve_synteny = NULL;
	uint64_t either_all_null_or_none_T1, either_all_null_or_none_T2, i;
	for(i=0;i<n_sequences;i++){
		//If one is not null, the other can't be either
		either_all_null_or_none_T1 = 0; either_all_null_or_none_T2 = 0;
		if(cons_order_A_B_T1[i] != NULL) either_all_null_or_none_T1++;
		if(cons_order_B_C_T1[i] != NULL) either_all_null_or_none_T1++;
		if(cons_order_A_B_T2[i] != NULL) either_all_null_or_none_T2++;
		if(cons_order_B_C_T2[i] != NULL) either_all_null_or_none_T2++;
		
		if(either_all_null_or_none_T1 == 1) return NULL; 
		if(either_all_null_or_none_T2 == 1) return NULL; 

		if(either_all_null_or_none_T1 == 2){
			if(cons_order_A_B_T1[i]->genome->id != cons_order_B_C_T1[i]->genome->id){
				return NULL; 
			}else{
				if(the_pointer_to_retrieve_synteny == NULL){
					//With the first block that will retrieve the synteny we are looking for is enough
					//So no need for the other ones
					the_pointer_to_retrieve_synteny = cons_order_A_B_T1[i];
					#ifdef VERBOSE
					printf("enter at@@@@@@@@@: "); printBlock(the_pointer_to_retrieve_synteny);
					#endif
				}
			}
		}
		
		if(either_all_null_or_none_T2 == 2){
			if(cons_order_A_B_T2[i]->genome->id != cons_order_B_C_T2[i]->genome->id){
				return NULL;
			}else{
				
			}
		}
	}
	return the_pointer_to_retrieve_synteny;
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
					throw "Overflown in order recomputation";
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

void remove_insertion_DNA(Block * a, Block * b, Block * c, uint64_t diff){
	uint64_t half_diff = diff/2;
	//printf("########->>>\n");printBlock(a);printBlock(b);printBlock(c);
	uint64_t b_point_AB = (a->end + b->start)/2;
	uint64_t b_point_BC = (b->end + c->start)/2;

	#ifdef VERBOSE
	printf("At seq: %"PRIu64"\n", a->genome->id);
	#endif
	if(diff % 2 == 0){
		#ifdef VERBOSE
		printf("(1)Moving %"PRIu64" to %"PRIu64"\n", b_point_BC, b_point_BC-half_diff);
		#endif
		memmove(&a->genome->seq[b_point_BC-half_diff], &a->genome->seq[b_point_BC], a->genome->len - b_point_BC);
		a->genome->len -= half_diff;
		#ifdef VERBOSE
		printf("Moving %"PRIu64" to %"PRIu64" a total ch of %"PRIu64" while having len: %"PRIu64"\n", b_point_AB+half_diff, b_point_AB, a->genome->len - b_point_AB + half_diff, a->genome->len);
		#endif
		memmove(&a->genome->seq[b_point_AB], &a->genome->seq[b_point_AB+half_diff], a->genome->len - (b_point_AB + half_diff));
		a->genome->len -= half_diff;
	}else{
		#ifdef VERBOSE
		printf("(2)Moving %"PRIu64" to %"PRIu64"\n", b_point_BC, b_point_BC-half_diff);
		#endif
		//Add extra one for the integer division
		memmove(&a->genome->seq[b_point_BC-half_diff-1], &a->genome->seq[b_point_BC], a->genome->len - b_point_BC - 1);
		a->genome->len -= (half_diff + 1);
		#ifdef VERBOSE
		printf("Moving %"PRIu64" to %"PRIu64"\n", b_point_AB+half_diff, b_point_AB);
		#endif
		memmove(&a->genome->seq[b_point_AB], &a->genome->seq[b_point_AB+half_diff], a->genome->len - (b_point_AB + half_diff));
		a->genome->len -= half_diff;
	}
	//getchar();
}

void remove_deletion_DNA(Block * a, Block * b, Block * c, uint64_t diff){
	uint64_t half_diff = diff/2;
	uint64_t i;
	#ifdef VERBOSE
	printf("At seq: %"PRIu64"\n", a->genome->id);
	#endif
	if(diff % 2 == 0){
		#ifdef VERBOSE
		printf("(1()Adding %"PRIu64" to %"PRIu64"\n", b->end, b->end+half_diff);
		#endif
		memmove(&a->genome->seq[b->end+half_diff], &a->genome->seq[b->end], a->genome->len - b->end);
		a->genome->len += half_diff;
		//Fill
		for(i=0;i<half_diff;i++) a->genome->seq[i+b->end] = 'N';
		#ifdef VERBOSE
		printf("Adding %"PRIu64" to %"PRIu64"\n", b->start, b->start+half_diff);
		#endif
		memmove(&a->genome->seq[b->start+half_diff], &a->genome->seq[b->start], a->genome->len - b->start);
		a->genome->len += half_diff;
		//Fill
		for(i=0;i<half_diff;i++) a->genome->seq[i+b->start] = 'N';
		
	}else{
		#ifdef VERBOSE
		printf("(2)Adding %"PRIu64" to %"PRIu64"\n", b->end, b->end+half_diff+1);
		#endif
		memmove(&a->genome->seq[b->end+half_diff+1], &a->genome->seq[b->end], a->genome->len - (b->end + 1));
		a->genome->len += (half_diff + 1);
		//Fill
		for(i=0;i<half_diff+1;i++) a->genome->seq[i+b->end] = 'N';
		#ifdef VERBOSE
		printf("Adding %"PRIu64" to %"PRIu64"\n", b->start, b->start+half_diff);
		#endif
		memmove(&a->genome->seq[b->start+half_diff], &a->genome->seq[b->start], a->genome->len - b->start);
		a->genome->len += half_diff;
		//Fill
		for(i=0;i<half_diff;i++) a->genome->seq[i+b->start] = 'N';
	}
	//getchar();
	
}

void apply_operation(Block * b, int64_t coordinates, int64_t order, uint64_t range1, uint64_t range2){

	Block * ptr = b->next;
	while(ptr != NULL){
		if(range1 <= ptr->end && ptr->start <= range2){
			ptr->order = (uint64_t)((int64_t) ptr->order + order);
			ptr->start = (uint64_t)((int64_t) ptr->start + coordinates);
			ptr->end = (uint64_t)((int64_t) ptr->end + coordinates);
		}
		ptr = ptr->next;
	}
	if(b->prev != NULL) ptr = b->prev; else return;

	while(ptr != NULL){
		if(range1 <= ptr->end && ptr->start <= range2){
			ptr->order = (uint64_t)((int64_t) ptr->order + order);
			ptr->start = (uint64_t)((int64_t) ptr->start + coordinates);
			ptr->end = (uint64_t)((int64_t) ptr->end + coordinates);
		}
		ptr = ptr->prev;
	}
}

void handle_indels(Synteny_list * A, Synteny_list * B, Synteny_list * C, uint64_t * indel_distance, 
	uint64_t n_sequences, uint64_t * genomes_block_count, uint64_t * indel_kept, uint64_t * indel_type,
	events_queue * operations_queue, uint64_t * t_insertions, uint64_t * t_deletions){

	uint64_t indel_used = 0;	
	memset(indel_distance, 0, n_sequences*sizeof(uint64_t));
	distance_between_blocks(indel_distance, A, C);
	
	/*
	printf("Got some concat here\n");
	if(A != NULL){ printSyntenyBlock(A->sb); printf("=was A=======000000\n");}
	if(B != NULL){ printSyntenyBlock(B->sb); printf("=was B=======000000\n");}
	if(C != NULL){ printSyntenyBlock(C->sb); printf("=was C=======000000\n");}
	getchar();
	*/

	//Add only those that were involved
	indel_used = 0;
	Synteny_block * sb_ptr = A->sb;
	while(sb_ptr != NULL){
		
		if(genomes_block_count[sb_ptr->b->genome->id] != 0){ 
			indel_kept[indel_used] = indel_distance[sb_ptr->b->genome->id];
			indel_used++;
		}
		sb_ptr = sb_ptr->next;
	}
	//Sort indels
	qsort(indel_kept, indel_used, sizeof(uint64_t), compare_distances_indel);
	//Get median
	uint64_t median = (uint64_t) median_from_vector(indel_kept, indel_used);
	//Generate table of types
	sb_ptr = C->sb;
	uint64_t diff;
	memset(indel_type, 0, n_sequences*sizeof(uint64_t));
	while(sb_ptr != NULL){
		//If genome was involved we have to add offset to the concat
		if(genomes_block_count[sb_ptr->b->genome->id] != 0){ 
			if(indel_distance[sb_ptr->b->genome->id] > median){
				indel_type[sb_ptr->b->genome->id] = INSERTION;
				*t_insertions = *t_insertions + 1;
				//subtract the difference
				diff = indel_distance[sb_ptr->b->genome->id] - median;
				#ifdef VERBOSE
				printf("My Diff: %"PRIu64"\n", diff);
				#endif
				
				//Modify sequence 
				//The diff can be removed either from breakpoint A-B or B-C
				//So remove "diff" from the end of A-> and the starting of <-C
				remove_insertion_DNA(sb_ptr->b->prev->prev, sb_ptr->b->prev, sb_ptr->b, diff);

				//Modify last block
				sb_ptr->b->start -= diff;
				sb_ptr->b->end -= diff;
								

				//Add to queue
				//rearrangement _r = { -((int64_t)(diff)), 0, sb_ptr->b->end, 0xFFFFFFFFFFFFFFFF, C->id, 1, sb_ptr->b->genome->id};
				//operations_queue->insert_event(_r);
				apply_operation(sb_ptr->b, -((int64_t)(diff)), 0, sb_ptr->b->end, 0xFFFFFFFFFFFFFFFF);

			}else if(indel_distance[sb_ptr->b->genome->id] < median){
				indel_type[sb_ptr->b->genome->id] = DELETION;
				*t_deletions = *t_deletions + 1;
				diff = median - indel_distance[sb_ptr->b->genome->id];
				#ifdef VERBOSE
				printf("My Diff: %"PRIu64"\n", diff);
				#endif
				//Modify sequence 
				remove_deletion_DNA(sb_ptr->b->prev->prev, sb_ptr->b->prev, sb_ptr->b, diff);

				//Modify last block
				sb_ptr->b->start += diff;
				sb_ptr->b->end += diff;
				
				//Add to queue
				//rearrangement _r = { ((int64_t)(diff)), 0, sb_ptr->b->end, 0xFFFFFFFFFFFFFFFF, C->id, 1, sb_ptr->b->genome->id};
				//operations_queue->insert_event(_r);
				apply_operation(sb_ptr->b, ((int64_t)(diff)), 0, sb_ptr->b->end, 0xFFFFFFFFFFFFFFFF);

			}else{
				indel_type[sb_ptr->b->genome->id] = NOTHING;
			}
		}
		sb_ptr = sb_ptr->next;
	}
	#ifdef VERBOSE
	printf("Median: %"PRIu64"\n", median);
	
	for(uint64_t i=0;i<n_sequences;i++){
		printf("TYP-O[%"PRIu64"]: %"PRIu64"\n", i, indel_type[i]);
	}
	#endif
}



void concat_synteny_blocks(Synteny_list ** A, Synteny_list ** B, Synteny_list ** C){
	//printf("I would like to concat\n");
	//printf("And it would look like this:\n");
	

	uint64_t i;
	Synteny_block * start_sb_ptr = (*A)->sb;
	Synteny_block * mid_ptr = (*B)->sb;
	Synteny_block * end_sb_ptr = (*C)->sb;
	Frags_list * fl_A, * fl_B;

	for(i=0;i<(*A)->synteny_level;i++){
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

	// All involved blocks in B and C must have their synteny changed
	Synteny_block * sb_ptr = (*B)->sb;
	while(sb_ptr != NULL){
		sb_ptr->b->present_in_synteny = (*A);
		sb_ptr = sb_ptr->next;
	}
	sb_ptr = (*C)->sb;
	while(sb_ptr != NULL){
		sb_ptr->b->present_in_synteny = (*A);
		sb_ptr = sb_ptr->next;
	}

	//printSyntenyBlock(A->sb);
	//getchar();

	//Remove intermediate synteny block list
	/*
	(*A)->next = (*C)->next;
	if((*C)->next != NULL) (*A)->next->prev = (*A);
	(*B) = (*A)->next;
	if((*B) != NULL) (*C) = (*B)->next; else (*B) = NULL;
	*/
	
	//(*A)->next = (*C)->next;
	//(*C)->next->prev = (*A);
	
	/*
	Synteny_block * sb;
    uint64_t synteny_level;
    struct synteny_list * next;
    struct synteny_list * prev;
    uint64_t id;
	*/
	//B and C cant be accessed now. Its not dangling because of the mempool.

	

	/*
	printf("Situation after concat\n");
	if(A != NULL){ printSyntenyBlock((*A)->sb); printf("=was A=======000000\n");}
	if((*B) != NULL){ printSyntenyBlock((*B)->sb); printf("=was B=======000000\n");}
	if((*C) != NULL){ printSyntenyBlock((*C)->sb); printf("=was C=======000000\n");}
	getchar();
	*/

	//All synteny lists should be updated from 

}




void reverse_reversion(Synteny_list * B, sequence_manager * sm, bool * genome_ids_affected){
	uint64_t i;
	Sequence * current;
	Synteny_block * sb_ptr = B->sb;
	Frags_list * fl_ptr; 
	//Get sequence handlers for those affected
	for(i=0;i<sm->get_number_of_sequences();i++){
		if(genome_ids_affected[i] == true){
			current = sm->get_sequence_by_label(i);
			//Find the blocks in particular
			//Since they are sorted
			while(sb_ptr != NULL && sb_ptr->b->genome->id != i){
				sb_ptr = sb_ptr->next;
			}
			if(sb_ptr != NULL){
				//Change reversion here in sequence
				uint64_t start = sb_ptr->b->start;
				uint64_t end = sb_ptr->b->end;
 
				inplace_reverse_and_complement(current->seq+start, end-start);


			}else{
				throw "Could not find sequence and/or block to reversion";
			}
		}
	}
	sb_ptr = B->sb;
	while(sb_ptr != NULL){
		//Change list of fragments for all since the reversion was effected
		fl_ptr = sb_ptr->b->f_list;
		while(fl_ptr != NULL){
			if(fl_ptr->f->strand == 'r') fl_ptr->f->strand = 'f';
			fl_ptr = fl_ptr->next;
		}
		sb_ptr = sb_ptr->next;
	}
}

void reverse_duplication(Synteny_list * A, Synteny_list * B, Synteny_list * C, Block * dup, hash_table * ht, events_queue * operations_queue, uint64_t last_s_id){
	if(dup == NULL) return;
	//Modify DNA sequence by removing block and shifting char bytes
	char * dna_ptr = dup->genome->seq;
	#ifdef VERBOSE
	printf("How could be wrong: %"PRIu64", %"PRIu64", %"PRIu64", %"PRIu64"\n", dup->start, dup->end, dup->genome->len, dup->genome->len - dup->end);
	#endif
	memmove(&dna_ptr[dup->start], &dna_ptr[dup->end], dup->genome->len - dup->end);
	//Change max len 
	dup->genome->len -= (dup->end - dup->start);
	Block * who_was_next = NULL;

	//Remove references from synteny list
	Synteny_block * sb_ptr = B->sb;
	Synteny_block * last = NULL;
	while(sb_ptr != NULL){

		if(isBlockEqualToWithOrder(sb_ptr->b, dup)){
			sb_ptr = sb_ptr->next;
			if(last != NULL){
				last->next = sb_ptr;
			}else{
				//Its the head
				B->sb = B->sb->next;
			}
			
		}
		last = sb_ptr;
		if(sb_ptr != NULL) sb_ptr = sb_ptr->next;
	}

	
	//Remove block from ht
	who_was_next = dup->next;
	ht->remove_block(dup);


	//Add operation to queue
	//Coordinates and order
	#ifdef VERBOSE
	printf("ENDINGSSSSSSS %"PRIu64"\n", A->id);
	//getchar();
	#endif
	rearrangement _r = { -((int64_t)(dup->end - dup->start)), -1, dup->end, 0xFFFFFFFFFFFFFFFF, B->id, 1, dup->genome->id};
	//operations_queue->insert_event(_r);
	apply_operation(dup, -((int64_t)(dup->end - dup->start)), -1, dup->end, 0xFFFFFFFFFFFFFFFF);

	//in case there were blocks after (now without order)
	
	
	while(who_was_next->present_in_synteny == dup->present_in_synteny){
		//Remove what was added
		#ifdef VERBOSE
		printf("Gone to\n");
		#endif
		who_was_next->order = (uint64_t)((int64_t)who_was_next->order + _r.mod_order);
		who_was_next->start = (uint64_t)((int64_t)who_was_next->start + _r.mod_coordinates);
		who_was_next->end = (uint64_t)((int64_t)who_was_next->end + _r.mod_coordinates);
		who_was_next = who_was_next->next;
	}
	
	//And Same for C
	sb_ptr = C->sb;
	while(sb_ptr != NULL && sb_ptr->b->genome->id <= dup->genome->id){
		
		if(sb_ptr->b->genome->id == dup->genome->id){
			sb_ptr->b->order = (uint64_t)((int64_t)sb_ptr->b->order + _r.mod_order);
			sb_ptr->b->start = (uint64_t)((int64_t)sb_ptr->b->start + _r.mod_coordinates);
			sb_ptr->b->end = (uint64_t)((int64_t)sb_ptr->b->end + _r.mod_coordinates);
		}

		sb_ptr = sb_ptr->next;
	}
}

bool reverse_tranposition(Synteny_list * A, Synteny_list * B, Synteny_list * C, Synteny_list * K1, Synteny_list * K2, Block ** blocks_to_move, Block ** blocks_to_stay, uint64_t n_sequences, events_queue * operations_queue){
	uint64_t i;

	// Remember: T1 holds those blocks with furthest distance (i.e. contained in K1 K2)

	// Move blocks contained in A_B_T either from B to between K1 and K2
	// or from between K1 and K2 to B
	Synteny_block * sb_ptr_A, * sb_ptr_C, * sb_ptr_K1, * sb_ptr_K2;


	for(i=0;i<n_sequences;i++){
		if(blocks_to_move[i] != NULL){
			//This block has to be moved
			#ifdef VERBOSE
			printf("I will try to move %"PRIu64"\n", blocks_to_move[i]->genome->id);
			#endif

			//Check if it is possible to exchange the block (no other block in its way)

			//Check if the previous block corresponds to A or to K1
			if(blocks_to_move[i]->prev->present_in_synteny == K1){
				#ifdef VERBOSE
				printf("Chose the k1 wpath\n"); 
				#endif
				sb_ptr_A = A->sb;
				sb_ptr_C = C->sb;

				while(sb_ptr_A != NULL && sb_ptr_A->b->genome->id != blocks_to_move[i]->genome->id){
					sb_ptr_A = sb_ptr_A->next;
				}
				//Have the previous block in A
				//Find the post block in C 
				while(sb_ptr_C != NULL && sb_ptr_C->b->genome->id != blocks_to_move[i]->genome->id){
					sb_ptr_C = sb_ptr_C->next;
				}
				//Check that both are not null 
				if(sb_ptr_A != NULL && sb_ptr_C != NULL){
					//Check that there is place 
					if(sb_ptr_A->b->next == sb_ptr_C->b){
						//Sufficient room?
						uint64_t l_breach = sb_ptr_C->b->start - sb_ptr_A->b->end;
						uint64_t l_move = blocks_to_move[i]->end - blocks_to_move[i]->start;
						uint64_t d_breach_midpoint = (sb_ptr_A->b->end + (l_breach/2));
						if(l_breach >= l_move){
							//Calculate middle point where to copy
							//And exchange dna
							inplace_dna_switch(blocks_to_move[i]->genome->seq, blocks_to_move[i]->genome->seq, d_breach_midpoint-(l_move/2), blocks_to_move[i]->start, l_move);

							//Calculate if the block is moved left or right
							//bool going_backwards = (blocks_to_move[i]->start >= sb_ptr_A->b->end);

							//Now we have to change blocks coordinates
							blocks_to_move[i]->start = d_breach_midpoint-(l_move/2);
							blocks_to_move[i]->end = blocks_to_move[i]->start + l_move;
							//Update prev and next
							Block * b_comes_after = blocks_to_move[i]->next;
							Block * b_comes_before = blocks_to_move[i]->prev;
							blocks_to_move[i]->prev->next = blocks_to_move[i]->next;
							b_comes_after->prev = b_comes_before;
							sb_ptr_A->b->next = blocks_to_move[i];
							sb_ptr_C->b->prev = blocks_to_move[i];
							blocks_to_move[i]->prev = sb_ptr_A->b;
							blocks_to_move[i]->next = sb_ptr_C->b;


							//Change order 
							blocks_to_move[i]->order = sb_ptr_A->b->order + 1;
							b_comes_after->order = b_comes_after->order - 1;
							//And the one in C as well 
							blocks_to_move[i]->next->order = sb_ptr_A->b->order + 1;

							//Insert queue operation to change orders and coordinates
							//if(going_backwards){
							//rearrangement _r = { 0, (int64_t)1, blocks_to_move[i]->next->start, b_comes_after->start, A->id, 1, blocks_to_move[i]->genome->id};
							//operations_queue->insert_event(_r);
							apply_operation(blocks_to_move[i], 0, 1, blocks_to_move[i]->next->start, b_comes_after->start);
							
							//}
							/*else{
								rearrangement _r = { 0, (int64_t)-1, b_comes_after->id, sb_ptr_C->b->id, A->id, 1, blocks_to_move[i]->genome->id};
								operations_queue->insert_event(_r);
							}
							*/

							//Change frags list
							//TODO

						}else{
							#ifdef VERBOSE
							printf("Problem at transposition: breach not long enough\n");
							#endif
							return false;
						}
					}else{
						printf("Problem at transposition: there is something in between\n"); 
						//if(sb_ptr_A->b != NULL && sb_ptr_A->b->next != NULL) printBlock(sb_ptr_A->b->next); else printf("First is null\n");
						//if(sb_ptr_C->b != NULL) printBlock(sb_ptr_C->b); else printf("Second is null\n");
						//getchar();
						#ifdef VERBOSE
						printDebugBlockOrderByGenome(A, 0);
						printDebugBlockOrderByGenome(A, 1);
						printDebugBlockOrderByGenome(A, 2);
						printDebugBlockOrderByGenome(A, 3);
						#endif
						//getchar();
						return false;
					}
				}
			}else{
				#ifdef VERBOSE
				printf("Chose the A wpath\n"); //getchar();
				#endif
				//It belongs to A
				sb_ptr_K1 = K1->sb;
				sb_ptr_K2 = K2->sb;

				while(sb_ptr_K1 != NULL && sb_ptr_K1->b->genome->id != blocks_to_move[i]->genome->id){
					sb_ptr_K1 = sb_ptr_K1->next;
				}
				//Have the previous block in A
				//Find the post block in C 
				while(sb_ptr_K2 != NULL && sb_ptr_K2->b->genome->id != blocks_to_move[i]->genome->id){
					sb_ptr_K2 = sb_ptr_K2->next;
				}
				//Check that both are not null 
				if(sb_ptr_K1 != NULL && sb_ptr_K2 != NULL){
					//Check that there is place 
					if(sb_ptr_K1->b->next == sb_ptr_K2->b){
						//Sufficient room?
						uint64_t l_breach = sb_ptr_K2->b->start - sb_ptr_K1->b->end;
						uint64_t l_move = blocks_to_move[i]->end - blocks_to_move[i]->start;
						uint64_t d_breach_midpoint = (sb_ptr_K1->b->end + (l_breach/2));
						if(l_breach >= l_move){
							//Calculate middle point where to copy
							//And exchange dna
							inplace_dna_switch(blocks_to_move[i]->genome->seq, blocks_to_move[i]->genome->seq, d_breach_midpoint-(l_move/2), blocks_to_move[i]->start, l_move);

							//Calculate if the block is moved left or right
							//bool going_backwards = (blocks_to_move[i]->start >= sb_ptr_A->b->end);

							//Now we have to change blocks coordinates
							blocks_to_move[i]->start = d_breach_midpoint-(l_move/2);
							blocks_to_move[i]->end = blocks_to_move[i]->start + l_move;
							//Update prev and next
							Block * comes_after = blocks_to_move[i]->next;
							comes_after->prev = blocks_to_move[i]->prev;
							//blocks_to_move[i]->next->order -= 1;
							blocks_to_move[i]->prev->next = blocks_to_move[i]->next;
							sb_ptr_K1->b->next = blocks_to_move[i];
							sb_ptr_K2->b->prev = blocks_to_move[i];
							blocks_to_move[i]->prev = sb_ptr_K1->b;
							blocks_to_move[i]->next = sb_ptr_K2->b;


							//Change order 
							blocks_to_move[i]->order = sb_ptr_K2->b->order - 1;
							//And the one in K1 as well 
							//blocks_to_move[i]->next->order = sb_ptr_K1->b->order - 1;

							//Insert queue operation to change orders and coordinates
							//if(going_backwards){
							//rearrangement _r = { 0, (int64_t)-1, comes_after->start, sb_ptr_K1->b->start, K2->id, 1, blocks_to_move[i]->genome->id};
							//operations_queue->insert_event(_r);

							apply_operation(blocks_to_move[i], 0, -1, comes_after->start, sb_ptr_K1->b->start);

							//}
							/*else{
								rearrangement _r = { 0, (int64_t)-1, b_comes_after->id, sb_ptr_C->b->id, A->id, 1, blocks_to_move[i]->genome->id};
								operations_queue->insert_event(_r);
							}
							*/

							//Change frags list
							//TODO

						}else{
							#ifdef VERBOSE
							printf("Problem at transposition: breach not long enough\n");
							#endif
							return false;
						}
					}else{
						
						printf("Problem at transposition: there is something in between\n"); 
						//if(sb_ptr_A->b != NULL && sb_ptr_A->b->next != NULL) printBlock(sb_ptr_A->b->next); else printf("First is null\n");
						//if(sb_ptr_C->b != NULL) printBlock(sb_ptr_C->b); else printf("Second is null\n");
						//getchar();
						#ifdef VERBOSE
						printDebugBlockOrderByGenome(A, 0);
						printDebugBlockOrderByGenome(A, 1);
						printDebugBlockOrderByGenome(A, 2);
						printDebugBlockOrderByGenome(A, 3);
						#endif
						//getchar();
						return false;
					}
				}



			}

						
		}
	}
	return true;
}




void detect_evolutionary_event(Synteny_list * sbl, sequence_manager * seq_man, uint32_t kmer_size, hash_table * blocks_ht, uint64_t * last_s_id){
	
	//Data structures needed
	uint64_t i;
	uint64_t n_sequences = seq_man->get_number_of_sequences();

	//First id in the synteny_list
	//uint64_t first_s_id = sbl->id;

	//Until nothing else cant be done, keep iterating
	bool stop_criteria = false; 
	bool had_modifying_event = false;

	//For telling whether the same number of blocks per genome is involved in synteny
	uint64_t * genomes_block_count = (uint64_t *) std::malloc(n_sequences*sizeof(uint64_t));
	if(genomes_block_count == NULL) terror("Could not allocate genome blocks counter");

	//To recompute order for the next blocks after an event
	uint64_t * order_offsets = (uint64_t *) std::malloc(n_sequences*sizeof(uint64_t));
	uint64_t * order_offsets_after_concat = (uint64_t *) std::malloc(n_sequences*sizeof(uint64_t));
	if(order_offsets == NULL || order_offsets_after_concat == NULL) terror("Could not allocate order offsets for after-events");

	//To check for indels (distances between blocks in concat)
	uint64_t * indel_distance = (uint64_t *) std::malloc(n_sequences*sizeof(uint64_t));
	uint64_t * indel_kept = (uint64_t *) std::malloc(n_sequences*sizeof(uint64_t));
	uint64_t * indel_type = (uint64_t *) std::malloc(n_sequences*sizeof(uint64_t));
	if(indel_distance == NULL || indel_type == NULL || indel_kept == NULL) terror("Could not allocate indels distance vector");

	//To check which genomes have reversion
	bool * genomes_affected = (bool *) std::malloc(n_sequences*sizeof(bool));
	if(genomes_affected == NULL) terror("Could not allocate vector to keep track of inversions");

	//To check that blocks are consecutive in their genome
	uint64_t * pairs_diff = (uint64_t *) std::malloc(n_sequences*sizeof(uint64_t));
	int64_t * pairs_diff_integer = (int64_t *) std::malloc(n_sequences*sizeof(int64_t));
	Block ** cons_order_A_B_T1 = (Block **) std::calloc(n_sequences, sizeof(Block *));
	Block ** cons_order_A_B_T2 = (Block **) std::calloc(n_sequences, sizeof(Block *));
	Block ** cons_order_B_C_T1 = (Block **) std::calloc(n_sequences, sizeof(Block *));
	Block ** cons_order_B_C_T2 = (Block **) std::calloc(n_sequences, sizeof(Block *));
	if(pairs_diff == NULL || pairs_diff_integer == NULL) terror("Could not allocate consecutive order of blocks array");
	if(cons_order_A_B_T1 == NULL || cons_order_B_C_T1 == NULL) terror("Could not allocate consecutive order of block pointers array (1)");
	if(cons_order_A_B_T2 == NULL || cons_order_B_C_T2 == NULL) terror("Could not allocate consecutive order of block pointers array (2)");

	//For strand matrices
	strand_matrix * sm_A, * sm_B, * sm_C;// * sm_D, * sm_E;
	//unsigned char ** _tmp1;//, ** _tmp2;
	sm_A = new strand_matrix(n_sequences);
	sm_B = new strand_matrix(n_sequences);
	sm_C = new strand_matrix(n_sequences);
	//sm_D = new strand_matrix(n_sequences);
	//sm_E = new strand_matrix(n_sequences);

	//For hits and frags computation
	dictionary_hash * words_dictionary = new dictionary_hash((uint64_t) (seq_man->get_maximum_length()/TABLE_RATE), seq_man->get_maximum_length(), kmer_size);
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
	memory_pool * mp = new memory_pool(2*POOL_SIZE);
	/*
	memory_pool * mp = new memory_pool(1,
	 seq_man->get_number_of_sequences()*sizeofSlist()*2 +
	  2*seq_man->get_number_of_sequences()*sizeof(unsigned char) +
	  seq_man->get_number_of_sequences()*sizeof(Slist *));
	*/
	//For NW computation
	char * seq_for_reverse = (char *) std::malloc(seq_man->get_maximum_length()*sizeof(char));
	if(seq_for_reverse == NULL) terror("Could not allocate reverse sequence");
	struct cell * mc, * f0, * f1;
	mc = (struct cell *) std::malloc(seq_man->get_maximum_length()*sizeofCell());
	f0 = (struct cell *) std::malloc(seq_man->get_maximum_length()*sizeofCell());
	f1 = (struct cell *) std::malloc(seq_man->get_maximum_length()*sizeofCell());
	if(mc == NULL || f0 == NULL || f1 == NULL) terror("Could not allocate rows for NW");

	strand_matrix * trans_matrix = new strand_matrix(n_sequences);

	
	// For handling rearragement operations
	rearrangement * current_rea;
	events_queue * operations_queue = new events_queue(n_sequences);

	//Lists of synteny blocks to address evolutionary events

	//Get first block (whoever it is)
	uint64_t current_head_block = 0;
	Bucket * ptr; Block * block_ptr;


	// Now we have an initial block and we can build a synteny list
	Synteny_list * A = NULL, * B = NULL, * C = NULL; //, * D = NULL, * E = NULL;
	uint64_t current_genome_start = 0; // Label of the genome where we start
	Synteny_block * sb_ptr;

	//To have some statistics
	uint64_t current_step = 0, current_concats = 0, t_concats = 0;
	uint64_t t_inversions = 0, t_duplications = 0, t_transpositions = 0;
	uint64_t t_insertions = 0, t_deletions = 0, unsolvable_reversions = 0;


	while(!stop_criteria || current_genome_start < n_sequences){
		
		//Display current iteration
		printf("\nAfter %"PRIu64" step(s):\n\tTotal concats: \t\t%"PRIu64", this round: %"PRIu64"\n", current_step++, t_concats, current_concats);
		printf("\tTotal inversions: \t%"PRIu64" from which %"PRIu64" are unsolvable.\n\tTotal duplications: \t%"PRIu64"\n", t_inversions, unsolvable_reversions, t_duplications);
		printf("\tTotal transpositions: \t%"PRIu64"\n", t_transpositions);
		printf("\tTotal insertions: \t%"PRIu64"\n", t_insertions);
		printf("\tTotal deletions: \t%"PRIu64"\n", t_deletions);
		current_concats = 0;
		getchar();

		

		// Get head of synteny list
		// Advance pointers
		//A = sbl;
		ptr = blocks_ht->get_key_at(0);
		while(ptr == NULL || ptr->b.genome->id != current_genome_start || ptr->b.present_in_synteny == NULL){
			if(ptr != NULL){
				ptr = ptr->next;
			}else{
				ptr = blocks_ht->get_key_at(current_head_block++);
			}
			if(current_head_block >= blocks_ht->get_size()) break;

			#ifdef VERBOSE
			printf("in loop ::: %"PRIu64" \n", current_genome_start); if(ptr != NULL) { printBlock(&ptr->b);} getchar();
			#endif 
		}
		current_head_block = 0;



		
		
		#ifdef VERBOSE
		fprintf(stdout, "^^^^^^--------------;;;;;;;;;;;;;;;;;;\n");
		fprintf(stdout, "THIS ROUND I AM USING %"PRIu64"\n", current_genome_start);
		getchar();
		#endif

		if(stop_criteria) current_genome_start++;
		//In case nothing gets done, stop iterating
		stop_criteria = true;
		

		B = ptr->b.present_in_synteny;
		if(ptr->b.prev != NULL && ptr->b.prev->present_in_synteny != NULL && ptr->b.prev->present_in_synteny != B) A = ptr->b.prev->present_in_synteny;
		if(ptr->b.next != NULL && ptr->b.next->present_in_synteny != NULL && ptr->b.next->present_in_synteny != B) C = ptr->b.next->present_in_synteny;
		block_ptr = &ptr->b;

		/*
		while(block_ptr != NULL){

			block_ptr = block_ptr->next;
			printf("---------------\n");
			if(B != NULL) B = block_ptr->present_in_synteny;
			if(block_ptr->prev != NULL && block_ptr->prev->present_in_synteny != NULL) A = block_ptr->prev->present_in_synteny; else A = NULL;
			if(block_ptr->next != NULL && block_ptr->next->present_in_synteny != NULL) C = block_ptr->next->present_in_synteny; else C = NULL;
			if(A != NULL) printSyntenyBlock(A->sb);
			if(B != NULL) printSyntenyBlock(B->sb);
			if(C != NULL) printSyntenyBlock(C->sb);
			printf("--------------\n");
		}
		*/
		

		//Set offset orders to zero
		//memset(order_offsets, 0, n_sequences*sizeof(int64_t));

		//Copy pointers of first consecutive blocks
		//if(A != NULL) B = A->next; else return; // A must exist (minimum)
		//if(B != NULL) C = B->next;
		//if(C != NULL) D = C->next;
		//if(D != NULL) E = D->next;

		//Generate their strand matrices
		sm_A->reset();
		sm_B->reset();
		sm_C->reset();
		if(A != NULL) sm_A->add_fragment_strands(A);
		if(B != NULL) sm_B->add_fragment_strands(B);
		if(C != NULL) sm_C->add_fragment_strands(C);
		//sm_D->add_fragment_strands(D);
		//sm_E->add_fragment_strands(E);
		/*
		#ifdef VERBOSE
		printf("BEFORE ALL !!!!!!!!!!!!!!!!!!!!!!!!!!\n");
		printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
		printf("BEFORE$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
		if(A != NULL){ printSyntenyBlock(A->sb); printf("=was A====with %"PRIu64"===000000\n", A->id);}
		if(B != NULL){ printSyntenyBlock(B->sb); printf("=was B====with %"PRIu64"===000000\n", B->id);}
		if(C != NULL){ printSyntenyBlock(C->sb); printf("=was C====with %"PRIu64"===000000\n", C->id);}

		operations_queue->print_queue();
		#endif
		*/
		//Apply queue actions for A and B (they wont be applied at iteration restart)
		//And only if last action was not a concat (since if it was a concat there was no new synteny block!)

		/*
		if(!had_modifying_event){
			if(A != NULL){
				#ifdef VERBOSE
				printf("RESTART AT A &&&&&&&&&&&&&&&&&&&&&&\n");	
				#endif
				sb_ptr = A->sb;
				while(sb_ptr != NULL){
					current_rea = operations_queue->get_aggregated_event(sb_ptr->b, A->id);
					#ifdef VERBOSE
					printf("In block: "); printBlock(sb_ptr->b);
					#endif
					if(current_rea != NULL){
						#ifdef VERBOSE
						printf("Applying to A (%"PRIu64")...:", A->id); printRearrangement(current_rea);
						#endif
						sb_ptr->b->order = (uint64_t)((int64_t) sb_ptr->b->order + current_rea->mod_order);
						sb_ptr->b->start = (uint64_t)((int64_t) sb_ptr->b->start + current_rea->mod_coordinates);
						sb_ptr->b->end = (uint64_t)((int64_t) sb_ptr->b->end + current_rea->mod_coordinates);
					}
					sb_ptr = sb_ptr->next;
				}
			}
			if(B != NULL){
				#ifdef VERBOSE
				printf("RESTART AT B &&&&&&&&&&&&&&&&&&&&&&\n");	
				#endif
				sb_ptr = B->sb;
				while(sb_ptr != NULL){
					current_rea = operations_queue->get_aggregated_event(sb_ptr->b, B->id);
					
					#ifdef VERBOSE
					printf("In block: "); printBlock(sb_ptr->b);
					#endif
					if(current_rea != NULL){
						#ifdef VERBOSE
						printf("Applying to B (%"PRIu64")...:", B->id); printRearrangement(current_rea);
						#endif
						sb_ptr->b->order = (uint64_t)((int64_t) sb_ptr->b->order + current_rea->mod_order);
						sb_ptr->b->start = (uint64_t)((int64_t) sb_ptr->b->start + current_rea->mod_coordinates);
						sb_ptr->b->end = (uint64_t)((int64_t) sb_ptr->b->end + current_rea->mod_coordinates);
					}
					sb_ptr = sb_ptr->next;
				}
			}
		}
		*/

		//had_modifying_event = false; //It would have been applied here

		while(B != NULL){ // AT least one to detect duplications

			//printf("Traversing BEFOREEEEEEEEE QUEUE\n");
			//traverse_synteny_list(sbl);
			#ifdef VERBOSE
			printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
			printf("BEFORE$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
			if(A != NULL){ printSyntenyBlock(A->sb); printf("=was A====with %"PRIu64"===000000\n", A->id);}
			if(B != NULL){ printSyntenyBlock(B->sb); printf("=was B====with %"PRIu64"===000000\n", B->id);}
			if(C != NULL){ printSyntenyBlock(C->sb); printf("=was C====with %"PRIu64"===000000 %p \n", C->id, C);}
			
			//operations_queue->print_queue();
			#endif

			//Traverse rearrangement queue and apply modifications
			/*
			if(B != NULL){
				#ifdef VERBOSE
				printf("In the GOOD AT B &&&&&&&&&&&&&&&&&&&%"PRIu64",\n", B->id);	
				#endif
				sb_ptr = B->sb;
				while(sb_ptr != NULL){
					current_rea = operations_queue->get_aggregated_event(sb_ptr->b, B->id);
					#ifdef VERBOSE
					printf("In block: "); printBlock(sb_ptr->b);
					#endif
					if(current_rea != NULL){
						#ifdef VERBOSE
						printf("Applying to B (%"PRIu64")...:", B->id); printRearrangement(current_rea);
						#endif
						sb_ptr->b->order = (uint64_t)((int64_t) sb_ptr->b->order + current_rea->mod_order);
						sb_ptr->b->start = (uint64_t)((int64_t) sb_ptr->b->start + current_rea->mod_coordinates);
						sb_ptr->b->end = (uint64_t)((int64_t) sb_ptr->b->end + current_rea->mod_coordinates);
					}
					sb_ptr = sb_ptr->next;
				}
			}
			*/
			/*
			if(C != NULL){
				#ifdef VERBOSE
				printf("In the GOOD  AT C &&&&&&&&&&&&&&&&&&&&&& %"PRIu64",\n", C->id);	
				#endif
				sb_ptr = C->sb;
				while(sb_ptr != NULL){
					current_rea = operations_queue->get_aggregated_event(sb_ptr->b, C->id);
					#ifdef VERBOSE
					printf("In block: "); printBlock(sb_ptr->b);
					#endif
					if(current_rea != NULL){
						#ifdef VERBOSE
						printf("Applying to C (%"PRIu64")...:", C->id); printRearrangement(current_rea);
						#endif
						sb_ptr->b->order = (uint64_t)((int64_t) sb_ptr->b->order + current_rea->mod_order);
						sb_ptr->b->start = (uint64_t)((int64_t) sb_ptr->b->start + current_rea->mod_coordinates);
						sb_ptr->b->end = (uint64_t)((int64_t) sb_ptr->b->end + current_rea->mod_coordinates);
					}
					sb_ptr = sb_ptr->next;
				}
			}
			*/
			/*
			if(had_modifying_event){
				if(B != NULL){
					#ifdef VERBOSE
					printf("In the GOOD AT B &&&&&&&&&&&&&&&&&&&&&&\n");	
					#endif
					sb_ptr = B->sb;
					while(sb_ptr != NULL){
						current_rea = operations_queue->get_aggregated_event(sb_ptr->b, B->id);
						#ifdef VERBOSE
						printf("In block: "); printBlock(sb_ptr->b);
						#endif
						if(current_rea != NULL){
							#ifdef VERBOSE
							printf("Applying to B (%"PRIu64")...:", B->id); printRearrangement(current_rea);
							#endif
							sb_ptr->b->order = (uint64_t)((int64_t) sb_ptr->b->order + current_rea->mod_order);
							sb_ptr->b->start = (uint64_t)((int64_t) sb_ptr->b->start + current_rea->mod_coordinates);
							sb_ptr->b->end = (uint64_t)((int64_t) sb_ptr->b->end + current_rea->mod_coordinates);
						}
						sb_ptr = sb_ptr->next;
					}
				}
			}
			if(C != NULL){
				#ifdef VERBOSE
				printf("In the GOOD  AT C &&&&&&&&&&&&&&&&&&&&&&\n");	
				#endif
				sb_ptr = C->sb;
				while(sb_ptr != NULL){
					current_rea = operations_queue->get_aggregated_event(sb_ptr->b, C->id);
					#ifdef VERBOSE
					printf("In block: "); printBlock(sb_ptr->b);
					#endif
					if(current_rea != NULL){
						#ifdef VERBOSE
						printf("Applying to C (%"PRIu64")...:", C->id); printRearrangement(current_rea);
						#endif
						sb_ptr->b->order = (uint64_t)((int64_t) sb_ptr->b->order + current_rea->mod_order);
						sb_ptr->b->start = (uint64_t)((int64_t) sb_ptr->b->start + current_rea->mod_coordinates);
						sb_ptr->b->end = (uint64_t)((int64_t) sb_ptr->b->end + current_rea->mod_coordinates);
					}
					sb_ptr = sb_ptr->next;
				}
				
				
			}
			*/
			//getchar();
			//rearrangement r = {1,2,3,4};
			//operations_queue->insert_event(r);
			//getchar();
			


			//Recompute order of the last one added to the list (because of concatenation)
			//recompute_orders_from_offset(order_offsets, 1, C);
			//if(!had_modifying_event) recompute_orders_from_offset(order_offsets, 1, B);
			
			//printDebugBlockOrderByGenome(E, 0);
			
			
			//printf("Traversing after QUEUE\n");
			//traverse_synteny_list(sbl);
			#ifdef VERBOSE
			printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
			printf("AFTER$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
			if(A != NULL){ printSyntenyBlock(A->sb); printf("=was A====with %"PRIu64"===000000\n", A->id);}
			if(B != NULL){ printSyntenyBlock(B->sb); printf("=was B====with %"PRIu64"===000000\n", B->id);}
			if(C != NULL){ printSyntenyBlock(C->sb); printf("=was C====with %"PRIu64"===000000\n", C->id);}
			#endif
			
			//if(D != NULL){ printSyntenyBlock(D->sb); printf("=was D=======000000\n");}
			//if(E != NULL){ printSyntenyBlock(E->sb); printf("=was E=======000000\n");}
			//getchar();
			

			// Transpositions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			
			if(A != NULL && B != NULL && C != NULL && synteny_level_across_lists(3, A, B, C)){
				//same synteny, check the number of genomes involved
				memset(genomes_block_count, 0, n_sequences*sizeof(uint64_t));
				if(genomes_involved_in_synteny(genomes_block_count, n_sequences, 3, A, B, C)){
					
					//Check that A and B have their order consecutive except for the transposed
					//Same for B and C (there can only be one transposed atm)
					memset(cons_order_A_B_T1, 0x0, n_sequences*sizeof(Block *));
					memset(cons_order_A_B_T2, 0x0, n_sequences*sizeof(Block *));
					memset(cons_order_B_C_T1, 0x0, n_sequences*sizeof(Block *));
					memset(cons_order_B_C_T2, 0x0, n_sequences*sizeof(Block *));
					
					if(consecutive_block_order_except_one(pairs_diff_integer, n_sequences, cons_order_A_B_T1, cons_order_A_B_T2, 2, A, B)
					&& consecutive_block_order_except_one(pairs_diff_integer, n_sequences, cons_order_B_C_T1, cons_order_B_C_T2, 2, B, C)){
						//That is, blocks can be separated using only two order discriminators
						//Now check that they are grouped the same way in both "clusters"
						Block * retrieve_synteny = NULL;
						retrieve_synteny = compare_order_clusters(cons_order_A_B_T1, cons_order_A_B_T2, cons_order_B_C_T1, cons_order_B_C_T2, n_sequences);

						//If the synteny block retrieved is not null, 

						if(retrieve_synteny != NULL){

							//printf("USING: "); printBlock(retrieve_synteny);

							//Retrieve the syteny from the block
							Synteny_list * sl_prev = NULL, * sl_after = NULL;
							Block * aux;
							aux = retrieve_synteny->prev;// blocks_ht->get_previous_block(retrieve_synteny);
							if(aux != NULL) sl_prev = aux->present_in_synteny;
							aux = retrieve_synteny->next;//blocks_ht->get_next_block(retrieve_synteny);
							if(aux != NULL) sl_after = aux->present_in_synteny;

							/*
							printf("##############################\n");
							if(sl_prev != NULL) printSyntenyBlock(sl_prev->sb);
							printf("##############################\n");
							if(sl_after != NULL) printSyntenyBlock(sl_after->sb);
							printf("##############################\n");
							*/
							if(sl_prev != A && sl_after != C && synteny_level_across_lists(5, A, B, C, sl_prev, sl_after)){
								#ifdef VERBOSE
								printSyntenyBlock(sl_prev->sb);
								printf("!!!!!!!!!!!!\n");
								printSyntenyBlock(sl_after->sb);
								printf("Detected transposition at B: %"PRIu64"\n", B->id);
								getchar();
								#endif
								// Debug ~~ check if the chain is broken
								

								//To reverse the transposition we have to align the B synteny block
								//To find out which block moved first
						
								memset(genomes_affected, false, n_sequences*sizeof(bool));

								//read_words_from_synteny_block_and_align(seq_man, B, kmer_size, words_dictionary, qfmat, qfmat_state);
								fill_quickfrag_matrix_NW(seq_man, seq_for_reverse, B, qfmat, qfmat_state, -5, -2, mc, f0, f1);
								printQuickFragMatrix(qfmat, qfmat_state, seq_man->get_number_of_sequences());

								mp->reset_to(0,0);
								//Note: The "genomes_affected" should hold which one are the blocks that moved (i.e. genome ids)
								Slist * dendro = UPGMA_joining_clustering(qfmat, qf_submat, qfmat_state, seq_man->get_number_of_sequences(), mp);
								//Now we know which blocks moved
								//cons_A_B_T1 has the "further" blocks
								//whereas cons_A_B_T2 has the closest
								trans_matrix->reset();
								for(uint64_t w=0;w<n_sequences;w++){
									//Put group number in diagonal
									if(cons_order_A_B_T1[w] != NULL) trans_matrix->set_strands(w, w, 1);
									if(cons_order_A_B_T2[w] != NULL) trans_matrix->set_strands(w, w, 2);
								}
								//Now fill each relation
								
								
								//Heuristic: Less changes
								//Count how many blocks in the groups T1 and T2
								uint64_t in_T1 = 0, in_T2 = 0;
								for(uint64_t w=0;w<n_sequences;w++){
									if(cons_order_A_B_T1[w] != NULL) in_T1++;
									if(cons_order_A_B_T2[w] != NULL) in_T2++;
								}
								

								if(in_T1 < in_T2){
									stop_criteria = !reverse_tranposition(A, B, C, sl_prev, sl_after, cons_order_A_B_T1, cons_order_A_B_T2, n_sequences, operations_queue);
								}else{
									stop_criteria = !reverse_tranposition(A, B, C, sl_prev, sl_after, cons_order_A_B_T2, cons_order_A_B_T1, n_sequences, operations_queue);
									
								}
								
								if(stop_criteria == false){
									t_transpositions++;
								}
								
								
								//stop_criteria = false;
							}else{
								#ifdef VERBOSE
								printf("Different synteny in ALL for transposition\n");
								#endif
							}
						}else{
							#ifdef VERBOSE
							printf("Retrieved synteny block is null for transposition\n");
							#endif
						}
					}else{
						#ifdef VERBOSE
						printf("Wrong consecutive order cant be discriminated for transposition\n");
						#endif
					}
				}else{
					#ifdef VERBOSE
					printf("Sorry, genomes involved differ in transposition\n");
					#endif
				}
			}else{
				#ifdef VERBOSE
				printf("A,B,C have different synteny in transposition\n");
				#endif
			}
			
			//getchar();


			// Duplications %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			// Left and right synteny must have same synteny level
			if(A != NULL && B != NULL && C != NULL && synteny_level_across_lists(2, A, C) > 0){

				memset(genomes_block_count, 0, n_sequences*sizeof(uint64_t));
				//If there is not the same number of genomes involved
				if(!genomes_involved_in_synteny(genomes_block_count, n_sequences, 1, B)){
					//There are duplications in B
					//Find those that have more synteny level


					for(i=0;i<n_sequences;i++){
						if(genomes_block_count[i] > 1){
							// Genome i has duplications
							#ifdef VERBOSE
							printf("Duplications in %"PRIu64"\n", i);
							#endif
							// Now we would know which blocks are duplications
							// list_of_dups = who_is_dup(sbl ...)

							// And reverse it

							Synteny_block * sb_ptr_dup = B->sb;
							Block * current_dup = NULL;
							
							while(sb_ptr_dup != NULL && current_dup == NULL){

								if(sb_ptr_dup->b->genome->id == i) current_dup = sb_ptr_dup->b;
								sb_ptr_dup = sb_ptr_dup->next;
							}
							
							if(current_dup != NULL) reverse_duplication(A, B, C, current_dup, blocks_ht, operations_queue, *last_s_id);
							// UNTIL HERE

							t_duplications++;
							//getchar();
							stop_criteria = false;
						}
					}
					//getchar();
				}else{
					#ifdef VERBOSE
					printf("Genomes involved not qualifying for duplication\n");
					#endif
				}
			}else{
				#ifdef VERBOSE
				printf("Wrong synteny for duplication\n");
				#endif
			}
			
			// Inversions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if(A != NULL && B != NULL && C != NULL && synteny_level_across_lists(3, A, B, C) > 0){
				if(sm_A->get_strands_type() != MIXED && 
				sm_A->get_strands_type() == sm_C->get_strands_type() &&
				sm_B->get_strands_type() == MIXED){

					memset(genomes_block_count, 0, n_sequences*sizeof(uint64_t));
					if(genomes_involved_in_synteny(genomes_block_count, n_sequences, 3, A, B, C)){
						if(consecutive_block_order(pairs_diff, 3, A, B, C)){
							#ifdef VERBOSE
							printf("Attention: this looks like a reversion\n");
							#endif

							//Clear out array of genomes affected
							memset(genomes_affected, false, n_sequences*sizeof(bool));

							//read_words_from_synteny_block_and_align(seq_man, B, kmer_size, words_dictionary, qfmat, qfmat_state);
							//printQuickFragMatrix(qfmat, qfmat_state, seq_man->get_number_of_sequences());
							fill_quickfrag_matrix_NW(seq_man, seq_for_reverse, B, qfmat, qfmat_state, -5, -2, mc, f0, f1);
							printQuickFragMatrix(qfmat, qfmat_state, seq_man->get_number_of_sequences());
							mp->reset_to(0,0);
							Slist * dendro = UPGMA_joining_clustering(qfmat, qf_submat, qfmat_state, seq_man->get_number_of_sequences(), mp);
							find_event_location(dendro, inversion, (void *) sm_B, genomes_affected);
							//getchar();
							
							// IMPORTANT: UPGMA should modify "genomes_affected" to tell which genomes (blocks) have the reversion in B
							/*
							int make_change;
							for(uint64_t w=0;w<n_sequences;w++){
								make_change = sm_B->do_forwards_require_less_changes(w);
								if(sm_B->get_frags_forward() >= sm_B->get_frags_reverse()){
									if(make_change == -1) genomes_affected[w] = 1;
								}else{
									if(make_change == 1) genomes_affected[w] = 1;
								}
							}
							*/
							unsigned char can_be_undone = 0;
							for(uint64_t w=0;w<n_sequences;w++){
								#ifdef VERBOSE
								printf("is %"PRIu64" reversed? %d\n", w, genomes_affected[w]);
								#endif
								if(genomes_affected[w] == true) can_be_undone = 1;

							}
							if(can_be_undone == 1){
								reverse_reversion(B, seq_man, genomes_affected);
								stop_criteria = false;
							}else{
								//If it is not solvable we can stop processing
								unsolvable_reversions++;
							}
							
							//Recalculate strand matrix (in case there is a concatenation)
							sm_B->reset();
							sm_B->add_fragment_strands(B);

							//printf("AFTER\n:");printSyntenyBlock(B->sb);

							t_inversions++;
							
						}else{
							#ifdef VERBOSE
							printf("Not consecutive order for inversion\n");
							#endif
						}
					}
				}else{
					#ifdef VERBOSE
					printf("Strands not qualifying for inversion\n");
					#endif
				}
			}else{
				#ifdef VERBOSE
				printf("Different synteny for inversion\n");
				#endif
			}

			// Concatenation %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			// @IMPORTANT: Currently blocks are not being removed or recomputed
			// in hash table after concatenation

			//Check they share the synteny level
			if(A != NULL && B != NULL && C != NULL && synteny_level_across_lists(3, A, B, C) > 0){
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
							#ifdef VERBOSE
							printf("Got some concat here\n");
							#endif
							//getchar();
							
							handle_indels(A, B, C, indel_distance, n_sequences, genomes_block_count, indel_kept, indel_type, operations_queue, &t_insertions, &t_deletions);
							

							concat_synteny_blocks(&A, &B, &C);
							
							#ifdef VERBOSE
							printf("Just in case after\n");
							if(A != NULL){ printSyntenyBlock(A->sb); printf("=was A=======000000with %"PRIu64"\n", A->id);}
							if(B != NULL){ printSyntenyBlock(B->sb); printf("=was B=======000000with %"PRIu64"\n", B->id);}
							if(C != NULL){ printSyntenyBlock(C->sb); printf("=was C=======000000with %"PRIu64"\n", C->id);}
							if(C != NULL && C->next != NULL){ printSyntenyBlock(C->next->sb); printf("=was C (NEXT!!!!!)=======000000with %"PRIu64" %p\n", C->next->id, C->next);}
							getchar();
							#endif
							
							t_concats++;
							current_concats++;
							
							//Add offset to orders
							sb_ptr = A->sb;
							while(sb_ptr != NULL){
								//If genome was involved we have to add offset to the concat
								if(genomes_block_count[sb_ptr->b->genome->id] != 0){ 
									
									//Add here also the coordinate indel
									//rearrangement _r = {0, -2, sb_ptr->b->start, 0xFFFFFFFFFFFFFFFF, C->id, 1, sb_ptr->b->genome->id};
									//operations_queue->insert_event(_r);
									apply_operation(sb_ptr->b, 0, -2, sb_ptr->b->start, 0xFFFFFFFFFFFFFFFF);

								}
								sb_ptr = sb_ptr->next;
							}
							//getchar();
							//Make the machine dont stop
							had_modifying_event = true;
							stop_criteria = false;

							// Repoint
							B = A;
							A = NULL; C = NULL;

						}else{
							#ifdef VERBOSE
							printf("Non consecutive order in blocks for concat...\n");
							#endif
						}
					}else{
						#ifdef VERBOSE
						printf("Genomes involved different number concat...\n"); //getchar();
						#endif
					}
				}else{
					#ifdef VERBOSE
					printf("Frags differ in strand for concat...\n"); //getchar();
					#endif
				}	
			}else{
				#ifdef VERBOSE
				printf("Different synteny levels for concat...\n"); //getchar();
				#endif
			}

			//Advance pointers
			
			


			if(block_ptr != NULL){
				block_ptr = block_ptr->next;
				while(block_ptr != NULL && block_ptr->present_in_synteny == NULL) block_ptr = block_ptr->next; 

				if(block_ptr != NULL){
					B = block_ptr->present_in_synteny;
					if(block_ptr->prev != NULL && block_ptr->prev->present_in_synteny != NULL && block_ptr->prev->present_in_synteny != B) A = block_ptr->prev->present_in_synteny; else A = NULL;
					if(block_ptr->next != NULL && block_ptr->next->present_in_synteny != NULL && block_ptr->next->present_in_synteny != B) C = block_ptr->next->present_in_synteny; else C = NULL;
					sm_A->reset();
					sm_B->reset();
					sm_C->reset();
					if(A != NULL){ sm_A->add_fragment_strands(A); }
					if(B != NULL){ sm_B->add_fragment_strands(B); }
					if(C != NULL){ sm_C->add_fragment_strands(C); }
				}else{
					//We reached to end of genome, get next
					A = NULL; B = NULL; C = NULL;
				}
			}else{
				stop_criteria = true;
			} 
			#ifdef VERBOSE
			printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
			printf("AFTER --------- CALLING THE ADVANCE$\n");
			if(A != NULL){ printSyntenyBlock(A->sb); printf("=was A=======00with %"PRIu64"\n", A->id);}
			if(B != NULL){ printSyntenyBlock(B->sb); printf("=was B=======00with %"PRIu64"\n", B->id);}
			if(C != NULL){ printSyntenyBlock(C->sb); printf("=was C=======00with %"PRIu64"\n", C->id);}
			if(C != NULL && C->next != NULL){ printSyntenyBlock(C->next->sb); printf("=was C (NEXT!!!!!)=======000000with %"PRIu64"\n", C->next->id);}
			getchar();
			#endif	

			/*
			if(had_modifying_event){
				A = sbl; sm_A->reset(); sm_A->add_fragment_strands(A);
				sm_B->reset(); sm_C->reset();
				if(A != NULL){ B = A->next; sm_B->add_fragment_strands(B); }
				if(B != NULL){ C = B->next; sm_C->add_fragment_strands(C); }
			}else{
				if(C != NULL && C->next != NULL){
					C = C->next; sm_C->reset(); sm_C->add_fragment_strands(C);
					goto _out_l;
				}else if(B != NULL && B->next != NULL){
					B = B->next; sm_B->reset(); sm_B->add_fragment_strands(B);
					if(B != NULL){ C = B->next; sm_C->add_fragment_strands(C);}
					goto _out_l;
				}else{

					//#ifdef VERBOSE
					printf("Here, advancing A with id: %"PRIu64"\n", A->id);
					//printSyntenyBlock(A->sb);
					//#endif
					if(A != NULL){
						A = A->next; sm_A->reset(); sm_A->add_fragment_strands(A);
						if(A != NULL){ B = A->next; sm_B->add_fragment_strands(B); }
						if(B != NULL){ C = B->next; sm_C->add_fragment_strands(C); }
						goto _out_l;
					}
				}
			}
			*/

			/*
			if(!had_modifying_event) A = A->next;

			sm_A->reset();
			sm_A->add_fragment_strands(A);
			
			if(A != NULL){
				if(!had_modifying_event){
					B = A->next;
					sm_B->reset();
					sm_B->add_fragment_strands(B);
				}
				
			}else{
				B = NULL;
			}
			if(B != NULL){
				C = B->next;
				sm_C->reset();
				sm_C->add_fragment_strands(C);
			}else{
				C = NULL;
			}
			*/
			/*
			_out_l:
			;
			*/
			//printf("skipping\n"); // To enable label
			
		}
	}

	printf("\nAfter %"PRIu64" step(s):\n\tTotal concats: \t\t%"PRIu64", this round: %"PRIu64"\n", current_step++, t_concats, current_concats);
	printf("\tTotal inversions: \t%"PRIu64" from which %"PRIu64" are unsolvable.\n\tTotal duplications: \t%"PRIu64"\n", t_inversions, unsolvable_reversions, t_duplications);
	printf("\tTotal transpositions: \t%"PRIu64"\n", t_transpositions);
	printf("\tTotal insertions: \t%"PRIu64"\n", t_insertions);
	printf("\tTotal deletions: \t%"PRIu64"\n", t_deletions);
	getchar();


	for(i=0;i<seq_man->get_number_of_sequences();i++){
		std::free(qfmat[i]);
		std::free(qf_submat[i]);
		std::free(qfmat_state[i]);
	}
	std::free(qf_submat);
	std::free(qfmat);
	std::free(qfmat_state);
	std::free(genomes_block_count);
	std::free(order_offsets);
	std::free(order_offsets_after_concat);
	std::free(pairs_diff);
	std::free(pairs_diff_integer);
	std::free(indel_distance);
	std::free(indel_kept);
	std::free(indel_type);
	std::free(cons_order_A_B_T1);
	std::free(cons_order_A_B_T2);
	std::free(cons_order_B_C_T1);
	std::free(cons_order_B_C_T2);
	std::free(genomes_affected);
	std::free(seq_for_reverse);
	std::free(mc);
	std::free(f0);
	std::free(f1);
	
	delete sm_A;
	delete sm_B;
	delete sm_C;
	//delete sm_D;
	//delete sm_E;

	delete operations_queue;

	delete mp;
	delete words_dictionary;
}

