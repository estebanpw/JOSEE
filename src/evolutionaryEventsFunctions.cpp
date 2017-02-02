#include <stdio.h>
#include <cstdlib>
#include <sys/time.h>
#include <inttypes.h>
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
	free(seq_orders);
	
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
	Synteny_block * curr_sb = NULL;


	for(i=0;i<ht->get_size();i++){
		ptr = ht->get_key_at(i);

		while(ptr != NULL){
			//Each block is here
			//For each block, add the blocks linked by the fragments
			memset(had_genome_bitmask, 0, n_seqs); //Reset genome counters
			synteny_level = 0; //Restart synteny level
			flptr = ptr->b.f_list;
			while(flptr != NULL){
				//Each fragment in the current block
				//Only if we did not already have the genome in the frags
				if(had_genome_bitmask[flptr->f->seqX] == (unsigned char)0) aux_block = ht->get_block_from_frag(flptr->f, 0);

				//Insert frag_x
				if(aux_block != NULL && aux_block->present_in_synteny == 0){
					aux_block->present_in_synteny = 1;
					Synteny_block * aux_sb = (Synteny_block *) mp->request_bytes(pre_comp_sb);
					aux_sb->b = aux_block; //insert at the head
					aux_sb->next = curr_sb;
					curr_sb = aux_sb;
					aux_block = NULL;
					had_genome_bitmask[flptr->f->seqX] = 1;
					synteny_level++;
					//printf("\t"); printBlock(aux_sb->b);
				}

				if(had_genome_bitmask[flptr->f->seqY] == (unsigned char)0) aux_block = ht->get_block_from_frag(flptr->f, 1);

				//Insert frag_y
				if(aux_block != NULL && aux_block->present_in_synteny == 0){
					aux_block->present_in_synteny = 1;
					Synteny_block * aux_sb = (Synteny_block *) mp->request_bytes(pre_comp_sb);
					aux_sb->b = aux_block; //insert at the head
					aux_sb->next = curr_sb;
					curr_sb = aux_sb;
					aux_block = NULL;
					had_genome_bitmask[flptr->f->seqY] = 1;
					synteny_level++;
					//printf("\t"); printBlock(aux_sb->b);
				}
				
				flptr = flptr->next;
			}

			// End synteny block
			if(synteny_level > 1){
				curr_sbl->sb = curr_sb;
				curr_sbl->synteny_level = synteny_level;
				curr_sbl->next = (Synteny_list *) mp->request_bytes(pre_comp_sbl);
				curr_sbl = curr_sbl->next;
				curr_sbl->next = NULL;
			}
			if(synteny_level == 1){
				//Since there is a minimum synteny, restore its level so that it can be used
				curr_sb->b->present_in_synteny = 0;
				mp->reset_n_bytes(pre_comp_sb);
			}
			curr_sb = NULL;
			//Go to next block
			ptr = ptr->next;
			
			
			//printf("broke stnyteny ---------------------------\n");
			//getchar();
		}
	}
	free(had_genome_bitmask);
	return sbl;
}



void has_reversion_in_truple(Synteny_block * a, Synteny_block * b, Synteny_block * c){
	/*

	Block * A, * B, * C;

	A = &a->b;
	B = &b->b;
	C = &c->b;

	//At this point we have the reference to each block A,B,C
	//We need to compute whether any of the fragments that generate B is reversed
	//Also that both A and C have no reversion themselves
	//Also that A,B,C have the same syteny

	if(A->synteny_level != B->synteny_level || A->synteny_level != C->synteny_level) return; //No equal synteny

	Frags_list * ptr;
	char b_strand;

	//Checking for no reversion in A
	ptr = a->b.f_list;
	while(ptr != NULL){
		if(ptr->f->strand == 'r') return;
		ptr = ptr->next;
	}
	//Checking for no reversion in C
	ptr = c->b.f_list;
	while(ptr != NULL){
		if(ptr->f->strand == 'r') return;
		ptr = ptr->next;
	}


	//Current truple
	printBlock(A);
	printBlock(B);
	

	//Check that there is a reversion in one or two fragments in B
	ptr = b->b.f_list;
	b_strand = ptr->f->strand; //First strand
	
	while(ptr != NULL){

		//Print fragments of b
		printf("\t"); printFragment(ptr->f);

		
		if(b_strand != ptr->f->strand){ //If there is one not equal to the first strand
			//Found reversion
			fprintf(stdout, "Found reversion between:\n");
			printBlock(A);
			printBlock(B);
			printBlock(C);
			getchar();
		}
		
		ptr = ptr->next;
	}
	printBlock(C); //So that it prints in order!
	getchar();
	*/
}


void detect_evolutionary_event(Synteny_list * sbl, sequence_manager * seq_man, uint32_t kmer_size){
	
	dictionary_hash * words_dictionary = new dictionary_hash(seq_man->get_maximum_length()/TABLE_RATE, seq_man->get_maximum_length(), kmer_size);
	Quickfrag ** qfmat = (Quickfrag **) std::malloc(seq_man->get_number_of_sequences()*seq_man->get_number_of_sequences()*sizeof(Quickfrag *));
	unsigned char ** qfmat_state = (unsigned char **) std::malloc(seq_man->get_number_of_sequences()*seq_man->get_number_of_sequences()*sizeof(unsigned char *));
	if(qfmat == NULL || qfmat_state == NULL) terror("Could not allocate pairwise alignment matrix (1)");
	uint64_t i;
	for(i=0;i<seq_man->get_number_of_sequences();i++){
		qfmat[i] = (Quickfrag *) std::malloc(seq_man->get_number_of_sequences()*sizeofQuickfrag());
		qfmat_state[i] = (unsigned char *) std::malloc(seq_man->get_number_of_sequences()*sizeof(unsigned char));
		if(qfmat[i] == NULL || qfmat_state[i] == NULL) terror("Could not allocate pairwsie alignment matrix (2)");
	}


	//Lists of synteny blocks to address evolutionary events
	Synteny_list * A, * B = NULL, * C = NULL, * D = NULL, * E = NULL;
	A = sbl;

	//Copy pointers of first consecutive blocks
	if(A != NULL) B = A->next; else return;
	if(B != NULL) C = B->next; else return; //Three at least
	if(C != NULL) D = C->next;
	if(D != NULL) E = D->next;

	while(A != NULL && B != NULL && C != NULL){ // AT least three


		

		//Evolutionary events
		
		//Work only with those that share the same synteny level 
		//Level 3 synteny
		if(A->synteny_level > 2 && A->synteny_level == B->synteny_level && B->synteny_level == C->synteny_level){

			//5-level
			if(C->synteny_level == D->synteny_level && D->synteny_level == E->synteny_level){

				if(C != NULL) printSyntenyBlock(C->sb);//sm_C->print_strand_matrix();printf("\nNEXT\n");
				//Lets say we pass list C as argument since its the one in the middle
				read_words_from_synteny_block_and_align(seq_man, C, kmer_size, words_dictionary, qfmat, qfmat_state);

				/*
				if(A != NULL) printSyntenyBlock(A->sb);//sm_A->print_strand_matrix();printf("\nNEXT\n");
				if(B != NULL) printSyntenyBlock(B->sb);//sm_B->print_strand_matrix();printf("\nNEXT\n");
				if(C != NULL) printSyntenyBlock(C->sb);//sm_C->print_strand_matrix();printf("\nNEXT\n");
				if(D != NULL) printSyntenyBlock(D->sb);//sm_D->print_strand_matrix();printf("\nNEXT\n");
				if(E != NULL) printSyntenyBlock(E->sb);//sm_E->print_strand_matrix();printf("\nNEXT\n");
				*/

			}else{

				if(B != NULL) printSyntenyBlock(B->sb);//sm_B->print_strand_matrix();printf("\nNEXT\n");
				read_words_from_synteny_block_and_align(seq_man, B, kmer_size, words_dictionary, qfmat, qfmat_state);

				/*
				if(A != NULL) printSyntenyBlock(A->sb);//sm_A->print_strand_matrix();printf("\nNEXT\n");
				if(B != NULL) printSyntenyBlock(B->sb);//sm_B->print_strand_matrix();printf("\nNEXT\n");
				if(C != NULL) printSyntenyBlock(C->sb);//sm_C->print_strand_matrix();printf("\nNEXT\n");
				*/
				
			}

			//getchar();

		}

		//advance pointers
		A = B;
		B = C;	
		C = D;
		D = E;

		//next iteration
		if(E != NULL) E = E->next;

	}

	for(i=0;i<seq_man->get_number_of_sequences();i++){
		free(qfmat[i]);
		free(qfmat_state[i]);
	}
	free(qfmat);
	free(qfmat_state);

	delete words_dictionary;
}

