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
			//ptr->b.id = ptr->b.order; // Order will be modified when applying events, but the id will not
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
	printf("DIFFTYPES: %"PRId64", %"PRId64"\n", diff_type_1, diff_type_2);
	/*
	for(i=0;i<n_sequences;i++){
		if(cons_order_T1[i] != NULL) {printf("in T1@%"PRIu64": ",i); printBlockWriteMode(cons_order_T1[i]);}
		if(cons_order_T2[i] != NULL) {printf("in T2@%"PRIu64": ",i); printBlockWriteMode(cons_order_T2[i]);}
	}
	*/
	
	if(diff_type_1 == 1){
		//Switch them so that const_order_T1 has the one with higher diff
		Block ** aux_list = cons_order_T1;
		cons_order_T1 = cons_order_T2;
		cons_order_T2 = aux_list;
	}
	
	for(i=0;i<n_sequences;i++){
		if(cons_order_T1[i] != NULL) {printf("in T1@%"PRIu64": ",i); printBlockWriteMode(cons_order_T1[i]);}
		if(cons_order_T2[i] != NULL) {printf("in T2@%"PRIu64": ",i); printBlockWriteMode(cons_order_T2[i]);}
	}

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
					printf("enter at@@@@@@@@@: "); printBlock(the_pointer_to_retrieve_synteny);
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
	
	//printSyntenyBlock(A->sb);
	//getchar();

	//Remove intermediate synteny block list
	(*A)->next = (*C)->next;
	if((*C)->next != NULL) (*A)->next->prev = (*A);
	(*B) = (*A)->next;
	if((*B) != NULL) (*C) = (*B)->next; else (*B) = NULL;
	//B and C cant be accessed now. Its not dangling because of the mempool.

	printf("Situation after concat\n");
	if(A != NULL){ printSyntenyBlock((*A)->sb); printf("=was A=======000000\n");}
	if((*B) != NULL){ printSyntenyBlock((*B)->sb); printf("=was B=======000000\n");}
	if((*C) != NULL){ printSyntenyBlock((*C)->sb); printf("=was C=======000000\n");}
	getchar();

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
				
				//Change list of fragments
				fl_ptr = sb_ptr->b->f_list;
				while(fl_ptr != NULL){
					if(fl_ptr->f->strand == 'r') fl_ptr->f->strand = 'f';
					fl_ptr = fl_ptr->next;
				}

			}else{
				throw "Could not find sequence and/or block to reversion";
			}
		}
	}
}

void reverse_duplication(Synteny_list * B, Synteny_list * C, Block * dup, hash_table * ht, events_queue * operations_queue, uint64_t last_s_id){
	if(dup == NULL) return;
	//Modify DNA sequence by removing block and shifting char bytes
	char * dna_ptr = dup->genome->seq;
	memmove(&dna_ptr[dup->start], &dna_ptr[dup->end], dup->genome->len - dup->end);
	//Change max len
	dup->genome->len -= (dup->genome->len - dup->end);

	//Remove references from synteny list
	Synteny_block * sb_ptr = B->sb;
	Synteny_block * last = NULL;
	while(sb_ptr != NULL){

		if(isBlockEqualToWithOrder(sb_ptr->b, dup)){
			sb_ptr = sb_ptr->next;
			if(last != NULL){
				last->next = sb_ptr;
			}
			
		}
		last = sb_ptr;
		if(sb_ptr != NULL) sb_ptr = sb_ptr->next;
	}

	//Remove block from ht
	ht->remove_block(dup);

	//Add operation to queue
	//Coordinates and order
	rearrangement _r = { -((int64_t)(dup->end - dup->start)), -1, B->id, 0xFFFFFFFFFFFFFFFF, B->id, '0', dup->genome->id};

	operations_queue->insert_event(_r);

	//in case there were blocks after (now without order)
	while(sb_ptr != NULL){
		if(isBlockEqualTo(sb_ptr->b, dup)){
			//Remove what was added
			sb_ptr->b->order = (uint64_t)((int64_t)sb_ptr->b->order + _r.mod_order);
			sb_ptr->b->start = (uint64_t)((int64_t)sb_ptr->b->start + _r.mod_coordinates);
			sb_ptr->b->end = (uint64_t)((int64_t)sb_ptr->b->end + _r.mod_coordinates);
		}
		sb_ptr = sb_ptr->next;
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

void reverse_tranposition(Synteny_list * A, Synteny_list * B, Synteny_list * C, Synteny_list * K1, Synteny_list * K2, Block ** cons_order_A_B_T1, Block ** cons_order_A_B_T2, bool * genomes_affected, uint64_t n_sequences, events_queue * operations_queue){
	uint64_t i;

	// Remember: T1 holds those blocks with furthest distance (i.e. contained in K1 K2)

	for(i=0;i<n_sequences;i++){
		if(cons_order_A_B_T1[i] != NULL && genomes_affected[cons_order_A_B_T1[i]->genome->id] == true){
			//Block has to be moved to its counterpart
			// Block cons_order_A_B_T1[i] gets order from A post

			//Find the block from A that has the same genome
			Synteny_block * sb_ptr = A->sb;
			while(sb_ptr != NULL){
				if(sb_ptr->b->genome->id == cons_order_A_B_T1[i]->genome->id) break;
				sb_ptr = sb_ptr->next;
			}
			
		}
		if(cons_order_A_B_T2[i] != NULL && genomes_affected[cons_order_A_B_T2[i]->genome->id] == true){
			//Block has to be moved to its counterpart
			// Block cons_order_A_B_T2[i] gets order from K1 post
		}
	}
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
	if(indel_distance == NULL) terror("Could not allocate indels distance vector");

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
	memory_pool * mp = new memory_pool(1, 2*POOL_SIZE);
	/*
	memory_pool * mp = new memory_pool(1,
	 seq_man->get_number_of_sequences()*sizeofSlist()*2 +
	  2*seq_man->get_number_of_sequences()*sizeof(unsigned char) +
	  seq_man->get_number_of_sequences()*sizeof(Slist *));
	*/

	
	// For handling rearragement operations
	rearrangement * current_rea;
	events_queue * operations_queue = new events_queue(n_sequences);

	//Lists of synteny blocks to address evolutionary events
	Synteny_list * A, * B = NULL, * C = NULL; //, * D = NULL, * E = NULL;
	Synteny_block * sb_ptr;

	//To have some statistics
	uint64_t current_step = 0, current_concats = 0, t_concats = 0;
	uint64_t t_inversions = 0, t_duplications = 0, t_transpositions = 0;


	while(!stop_criteria){
		
		//Display current iteration
		printf("\nAfter %"PRIu64" step(s):\n\tTotal concats: %"PRIu64", this round: %"PRIu64"\n", current_step++, t_concats, current_concats);
		printf("\tTotal inversions: %"PRIu64".\n\tTotal duplications: %"PRIu64"\n", t_inversions, t_duplications);
		printf("\tTotal transpositions: %"PRIu64"\n", t_transpositions);
		current_concats = 0;
		getchar();

		//In case nothing gets done, stop iterating
		stop_criteria = true;

		//Get head of synteny list
		A = sbl; 

		//Set offset orders to zero
		memset(order_offsets, 0, n_sequences*sizeof(int64_t));

		//Copy pointers of first consecutive blocks
		if(A != NULL) B = A->next; else return; // A must exist (minimum)
		if(B != NULL) C = B->next;
		//if(C != NULL) D = C->next;
		//if(D != NULL) E = D->next;


		//Generate their strand matrices
		sm_A->reset();
		sm_B->reset();
		sm_C->reset();
		sm_A->add_fragment_strands(A);
		sm_B->add_fragment_strands(B);
		sm_C->add_fragment_strands(C);
		//sm_D->add_fragment_strands(D);
		//sm_E->add_fragment_strands(E);

		while(A != NULL){ // AT least one to detect duplications



			//Traverse rearrangement queue and apply modifications
			if(had_modifying_event){
				if(B != NULL){
					sb_ptr = B->sb;
					while(sb_ptr != NULL){
						current_rea = operations_queue->get_aggregated_event(sb_ptr->b, B->id);
						if(current_rea != NULL){
							sb_ptr->b->order = (uint64_t)((int64_t) sb_ptr->b->order + current_rea->mod_order);
							sb_ptr->b->start = (uint64_t)((int64_t) sb_ptr->b->start + current_rea->mod_coordinates);
							sb_ptr->b->end = (uint64_t)((int64_t) sb_ptr->b->end + current_rea->mod_coordinates);
						}
						sb_ptr = sb_ptr->next;
					}
				}
			}
			if(C != NULL){
				
				sb_ptr = C->sb;
				while(sb_ptr != NULL){
					current_rea = operations_queue->get_aggregated_event(sb_ptr->b, C->id);
					if(current_rea != NULL){
						sb_ptr->b->order = (uint64_t)((int64_t) sb_ptr->b->order + current_rea->mod_order);
						sb_ptr->b->start = (uint64_t)((int64_t) sb_ptr->b->start + current_rea->mod_coordinates);
						sb_ptr->b->end = (uint64_t)((int64_t) sb_ptr->b->end + current_rea->mod_coordinates);
					}
					sb_ptr = sb_ptr->next;
				}
				
				
				
				/*
				while((current_rea = operations_queue->get_next_element(C->id-2)) != NULL){//Minus 2 because we check on C
					//apply current_rea
					printf("Applying: \n"); printRearrangement(current_rea); printf("CUZ I AM: %"PRIu64"\n", C->id);
					printf("TO\n");
					
					if(A != NULL){ printSyntenyBlock(A->sb); printf("=was A=======000000\n");}
					if(B != NULL){ printSyntenyBlock(B->sb); printf("=was B=======000000\n");}
					if(C != NULL){ printSyntenyBlock(C->sb); printf("=was C=======000000\n");}


					if(had_modifying_event) apply_queue_operation(current_rea, B); //Because B enters new
					apply_queue_operation(current_rea, C);


					printf("Yields\n");
					if(A != NULL){ printSyntenyBlock(A->sb); printf("=was A=======000000\n");}
					if(B != NULL){ printSyntenyBlock(B->sb); printf("=was B=======000000\n");}
					if(C != NULL){ printSyntenyBlock(C->sb); printf("=was C=======000000\n");}
					getchar();

				}
				*/
			}
			//rearrangement r = {1,2,3,4};
			//operations_queue->insert_event(r);
			//getchar();
			


			//Recompute order of the last one added to the list (because of concatenation)
			//recompute_orders_from_offset(order_offsets, 1, C);
			//if(!had_modifying_event) recompute_orders_from_offset(order_offsets, 1, B);
			
			//printDebugBlockOrderByGenome(E, 0);
			
			had_modifying_event = false;
			
			printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
			printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
			if(A != NULL){ printSyntenyBlock(A->sb); printf("=was A=======000000\n");}
			if(B != NULL){ printSyntenyBlock(B->sb); printf("=was B=======000000\n");}
			if(C != NULL){ printSyntenyBlock(C->sb); printf("=was C=======000000\n");}
			
			//if(D != NULL){ printSyntenyBlock(D->sb); printf("=was D=======000000\n");}
			//if(E != NULL){ printSyntenyBlock(E->sb); printf("=was E=======000000\n");}
			getchar();

			// Transpositions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			
			
			if(B != NULL && C != NULL && synteny_level_across_lists(3, A, B, C)){
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

							printf("USING: "); printBlock(retrieve_synteny);

							//Retrieve the syteny from the block
							Synteny_list * sl_prev = NULL, * sl_after = NULL;
							Block * aux;
							aux = retrieve_synteny->prev;// blocks_ht->get_previous_block(retrieve_synteny);
							if(aux != NULL) sl_prev = aux->present_in_synteny;
							aux = retrieve_synteny->next;//blocks_ht->get_next_block(retrieve_synteny);
							if(aux != NULL) sl_after = aux->present_in_synteny;

							printf("##############################\n");
							if(sl_prev != NULL) printSyntenyBlock(sl_prev->sb);
							printf("##############################\n");
							if(sl_after != NULL) printSyntenyBlock(sl_after->sb);
							printf("##############################\n");

							if(sl_prev != A && sl_after != C && synteny_level_across_lists(5, A, B, C, sl_prev, sl_after)){
								printSyntenyBlock(sl_prev->sb);
								printSyntenyBlock(sl_after->sb);
								printf("Detected transposition at B\n");

								//To reverse the transposition we have to align the B synteny block
								//To find out which block moved first
						
								memset(genomes_affected, false, n_sequences*sizeof(bool));
								read_words_from_synteny_block_and_align(seq_man, B, kmer_size, words_dictionary, qfmat, qfmat_state);
								mp->reset_to(0,0);
								//Note: The "genomes_affected" should hold which one are the blocks that moved (i.e. genome ids)
								UPGMA_joining_clustering(qfmat, qf_submat, qfmat_state, seq_man->get_number_of_sequences(), mp, genomes_affected);

								//Now we know which blocks moved
								//cons_A_B_T1 has the "further" blocks
								//whereas cons_A_B_T2 has the closest
								

								t_transpositions++;
								//getchar();
								stop_criteria = false;
							}else{
								printf("Different synteny in ALL for transposition\n");
							}
						}else{
							printf("Retrieved synteny block is null for transposition\n");
						}
					}else{
						printf("Wrong consecutive order cant be discriminated for transposition\n");
					}
				}else{
					printf("Sorry, genomes involved differ in transposition\n");
				}
			}else{
				printf("A,B,C have different synteny in transposition\n");
			}
			
			//getchar();


			// Duplications %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			// Left and right synteny must have same synteny level
			if(B != NULL && C != NULL && synteny_level_across_lists(2, A, C) > 0){

				memset(genomes_block_count, 0, n_sequences*sizeof(uint64_t));
				//If there is not the same number of genomes involved
				if(!genomes_involved_in_synteny(genomes_block_count, n_sequences, 1, B)){
					//There are duplications in B
					//Find those that have more synteny level
					for(i=0;i<n_sequences;i++){
						if(genomes_block_count[i] > 1){
							// Genome i has duplications
							printf("Duplications in %"PRIu64"\n", i);

							// Now we would know which blocks are duplications
							// list_of_dups = who_is_dup(sbl ...)

							//And reverse it

							// REMOVE THIS PART
							Block * current_dup = NULL;
							
							if(i == 0){
								current_dup = B->sb->next->b;
								printf("USING: "); printBlock(current_dup);
							}
							
							if(i == 2){
								current_dup = B->sb->next->next->b;
								printf("USING: "); printBlock(current_dup);
							}
							
							if(current_dup != NULL) reverse_duplication(B, C, current_dup, blocks_ht, operations_queue, *last_s_id);
							// UNTIL HERE

							t_duplications++;
							//getchar();
							stop_criteria = false;
						}
					}
					//getchar();
				}else{
					printf("Genomes involved not qualifying for duplication\n");
				}
			}else{
				printf("Wrong synteny for duplication\n");
			}
			
			// Inversions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
			if(B != NULL && C != NULL && synteny_level_across_lists(3, A, B, C) > 0){
				if(sm_A->get_strands_type() != MIXED && 
				sm_A->get_strands_type() == sm_C->get_strands_type() &&
				sm_B->get_strands_type() == MIXED){

					memset(genomes_block_count, 0, n_sequences*sizeof(uint64_t));
					if(genomes_involved_in_synteny(genomes_block_count, n_sequences, 3, A, B, C)){
						if(consecutive_block_order(pairs_diff, 3, A, B, C)){
							printf("Attention: this looks like a reversion\n");

							//Clear out array of genomes affected
							memset(genomes_affected, false, n_sequences*sizeof(bool));

							read_words_from_synteny_block_and_align(seq_man, B, kmer_size, words_dictionary, qfmat, qfmat_state);
							mp->reset_to(0,0);
							UPGMA_joining_clustering(qfmat, qf_submat, qfmat_state, seq_man->get_number_of_sequences(), mp, genomes_affected);
							//getchar();
							
							// IMPORTANT: UPGMA should modify "genomes_affected" to tell which genomes (blocks) have the reversion in B

							//genomes_affected[0] = true;
							reverse_reversion(B, seq_man, genomes_affected);
							//Recalculate strand matrix (in case there is a concatenation)
							sm_B->reset();
							sm_B->add_fragment_strands(B);

							printf("AFTER\n:");printSyntenyBlock(B->sb);

							t_inversions++;
							stop_criteria = false;
						}else{
							printf("Not consecutive order for inversion\n");
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
			if(B != NULL && C != NULL && synteny_level_across_lists(3, A, B, C) > 0){
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

							//First check if there are any indels
							memset(indel_distance, 0, n_sequences*sizeof(uint64_t));
							distance_between_blocks(indel_distance, A, B);
							for(i=0;i<n_sequences;i++){
								printf("D[%"PRIu64"] -:: %"PRIu64"\n", i, indel_distance[i]);
							}
							printf("Got some concat here\n");
							if(A != NULL){ printSyntenyBlock(A->sb); printf("=was A=======000000\n");}
							if(B != NULL){ printSyntenyBlock(B->sb); printf("=was B=======000000\n");}
							if(C != NULL){ printSyntenyBlock(C->sb); printf("=was C=======000000\n");}
							getchar();
							concat_synteny_blocks(&A, &B, &C);
							printf("Just in case after\n");
							if(A != NULL){ printSyntenyBlock(A->sb); printf("=was A=======000000\n");}
							if(B != NULL){ printSyntenyBlock(B->sb); printf("=was B=======000000\n");}
							if(C != NULL){ printSyntenyBlock(C->sb); printf("=was C=======000000\n");}
							getchar();
							t_concats++;
							current_concats++;
							
							//Add offset to orders
							Synteny_block * sb_ptr = A->sb;
							while(sb_ptr != NULL){
								//If genome was involved we have to add offset to the concat
								if(genomes_block_count[sb_ptr->b->genome->id] != 0){
																	
									rearrangement _r = {0, -2, sb_ptr->b->order, 0xFFFFFFFFFFFFFFFF, A->id, '0', sb_ptr->b->genome->id};
									operations_queue->insert_event(_r);

									//order_offsets[i] += 2; //Two because two blocks are shrinked into one
								}
								sb_ptr = sb_ptr->next;
							}
							//Make the machine dont stop
							had_modifying_event = true;
							stop_criteria = false;
						}else{
							printf("Non consecutive order in blocks for concat...\n");
						}
					}else{
						printf("Genomes involved different number concat...\n"); //getchar();
					}
				}else{
					printf("Frags differ in strand for concat...\n"); //getchar();
				}	
			}else{
				printf("Different synteny levels for concat...\n"); //getchar();
			}

			//Advance pointers
			/*
			printf("$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$$\n");
			printf("BEFORE CALLING THE ADVANCE$\n");
			if(A != NULL){ printSyntenyBlock(A->sb); printf("=was A=======000000\n");}
			if(B != NULL){ printSyntenyBlock(B->sb); printf("=was B=======000000\n");}
			if(C != NULL){ printSyntenyBlock(C->sb); printf("=was C=======000000\n");}
			getchar();
			*/
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
			
	
			
		}
	}

	printf("\nAfter %"PRIu64" step(s):\n\tTotal concats: %"PRIu64", this round: %"PRIu64"\n", current_step++, t_concats, current_concats);
	printf("\tTotal inversions: %"PRIu64".\n\tTotal duplications: %"PRIu64"\n", t_inversions, t_duplications);
	printf("\tTotal transpositions: %"PRIu64"\n", t_transpositions);


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
	std::free(pairs_diff_integer);
	std::free(indel_distance);
	std::free(cons_order_A_B_T1);
	std::free(cons_order_A_B_T2);
	std::free(cons_order_B_C_T1);
	std::free(cons_order_B_C_T2);
	std::free(genomes_affected);
	
	delete sm_A;
	delete sm_B;
	delete sm_C;
	//delete sm_D;
	//delete sm_E;

	delete operations_queue;

	delete mp;
	delete words_dictionary;
}

