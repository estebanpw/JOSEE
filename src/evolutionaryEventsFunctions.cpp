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



void map_frags_to_genomes(unsigned char ** map_table, struct FragFile * frags, uint64_t n_frags){

	uint64_t i, j, from, to, seq;
	//For all frags
	for(i=0;i<n_frags;i++){

		seq = frags[i].seqX;
		from = frags[i].xStart;
		to = frags[i].xEnd;

		//Map coordinates in frag for seqX, which is always forward
		for(j=from;j<to;j++){
			if(j == from) map_table[seq][j] = OPENFRAG;
			if(j == to) map_table[seq][j] = CLOSEFRAG;
			
			if(map_table[seq][j] == NOFRAG){
				map_table[seq][j] = COVERFRAG;
			}
		}

		//Map coordinates in frag for seqY, which might be reversed
		//Remember RAMGECKO coordinates are global respective to forward and with Ystart > Yend when reversed
		//So reversed should only be switched

		if(frags[i].strand == 'f'){
			seq = frags[i].seqY;
			from = frags[i].yStart;
			to = frags[i].yEnd;
		}else{
			seq = frags[i].seqY;
			from = frags[i].yEnd;
			to = frags[i].yStart;	
		}

		for(j=from;j<to;j++){
			if(map_table[seq][j] == NOFRAG || map_table[seq][j] == COVERFRAG){
				if(j == from) map_table[seq][j] = OPENFRAG;
				else if(j == to) map_table[seq][j] = CLOSEFRAG;
				else map_table[seq][j] = COVERFRAG;
			}
		}	
	}

	//At this point all coordinates have been mapped for the current fragments
}



inline void copyFragWithNewCoordinates(struct FragFile * destination, struct FragFile * source, uint64_t xStart, uint64_t yStart, uint64_t xEnd, uint64_t yEnd, uint64_t len){
    destination->xStart = xStart;
    destination->xEnd = xEnd;

	destination->yStart = yStart;
    destination->yEnd = yEnd;	
	/*
    if(source->strand == 'f'){
		destination->yStart = yStart;
    	destination->yEnd = yEnd;		
	}else{
		destination->yStart = yEnd;
    	destination->yEnd = yStart;
	}
	*/

    destination->length = len;
    destination->seqX = source->seqX;
    destination->seqY = source->seqY;
    destination->strand = source->strand;
}

struct FragFile * trim_fragments_and_map(unsigned char ** map_table, struct FragFile * frags, uint64_t * n_frags, uint64_t min_len, Sequence * sequences){

	struct FragFile new_frag;
	struct FragFile * list_new_frags;
	uint64_t list_reallocs = 1;
	uint64_t new_frags_number = 0;
	uint64_t size_fragment = sizeofFragment(); //To not compute it every time

	
	//Allocate memory
	list_new_frags = (struct FragFile *) std::malloc(INIT_TRIM_FRAGS*sizeofFragment());
	if(list_new_frags == NULL) terror("Could not allocate memory for list of new fragments in trimming");

	//Start trimming process
	uint64_t i, jX, jY, fromX, fromY, toX, toY, seqX, seqY, frag_len;
	char strand;
	uint64_t cur_new_len;

	for(i=0; i<*n_frags; i++){
		
		//Copy frag values
		fromX = frags[i].xStart; 
		toX = frags[i].xEnd; 
		
		
		strand = frags[i].strand;
		/*
		if(strand == 'f'){
			fromY = frags[i].yStart;
			toY = frags[i].yEnd;	
		}else{
			fromY = frags[i].yEnd;	
			toY = frags[i].yStart;
		}
		*/
		fromY = frags[i].yStart;
		toY = frags[i].yEnd;	
		
		

		

		seqX = frags[i].seqX; seqY = frags[i].seqY;

		frag_len = frags[i].length;
		cur_new_len = 1;
		jX = fromX+1;
		if(strand == 'f') jY = fromY+1; else jY = fromY-1;

		/*
		if(seqX == 0 && frags[i].xEnd == 5165){
			printFragment(&frags[i]);
			getchar();
		}
		*/

		while(jX < toX && jY != toY){
			//Check how long until there is a break (by starting or ending of frag)
			while(cur_new_len < frag_len && jX < sequences[seqX].len && jY < sequences[seqY].len && jY >= 0){

				if(map_table[seqX][jX] != COVERFRAG) break;
				if(map_table[seqY][jY] != COVERFRAG) break;

				cur_new_len++; jX++;
				if(strand == 'f'){ jY++; }else{ if(jY > 0) jY--; else break;} //To scape the buffer overflow of uints64
			}

			//At this point, jX and jY hold the ending coordinates of the new fragment
			//And fromX and fromY hold the starting coordinates
			if(cur_new_len >= min_len){ //Filtering

				//The fragment must be snipped out and saved
				copyFragWithNewCoordinates(&new_frag, &frags[i], fromX, fromY, jX, jY, cur_new_len);
				memcpy(&list_new_frags[new_frags_number], &new_frag, size_fragment);
				new_frags_number++;
				//Check if we need to realloc the list of new frags
				if(new_frags_number == list_reallocs*INIT_TRIM_FRAGS){
					list_reallocs++;
					list_new_frags = (struct FragFile *) std::realloc(list_new_frags, list_reallocs*INIT_TRIM_FRAGS*size_fragment);
					if(list_new_frags == NULL) terror("Could not realloc fragments on the trimming process");
				}

				


				//And set the mapping grid to the new values

				/*
				if(seqY == 3 && jY == 209 && fromY == 2){
					printf("Got it with fromY: %"PRIu64" jY: %"PRIu64" toY: %"PRIu64"\n", fromY, jY, toY);
					printFragment(&frags[i]);
				}
				*/

				//Close where you finished
				map_table[seqX][jX] = CLOSEFRAG;
				map_table[seqY][jY] = CLOSEFRAG;

				//Set open for this one
				
				map_table[seqX][fromX] = OPENFRAG;
				map_table[seqY][fromY] = OPENFRAG;
				

				//Open next if it was cut in between
				/*
				if(jX+1 < sequences[seqX].len && jX+1 < toX) map_table[seqX][jX+1] = OPENFRAG;
				if(strand == 'f'){
					if(jY+1 < sequences[seqY].len && jY+1 < toY) map_table[seqY][jY+1] = OPENFRAG;	
				}else{
					if(jY > 0 && jY-1 < fromY) map_table[seqY][jY-1] = OPENFRAG;
				}
				*/
				//Set the fromX and fromY to 1 (start frag) again in case this is not the first time we split

				/* DEBUG PURPOSES */
				/*
				if(frags[i].seqY == 3 && frags[i].yEnd == 5375){
					printf("Created: "); printFragment(&new_frag);
					print_maptable_portion(map_table, frags[i].yStart - 20, frags[i].yEnd + 20, 50, frags[i].seqY);
					getchar();
				}
				*/

			}
			//If you are here, either the fragment was too short, or was written correctly or we are at the end of the frag
			//Just keep going
			//Copy frag values
			
			/*
			fromX = jX+1; //One to move from an ending 3 to an opening 1
			if(strand == 'f') fromY = jY+1; else fromY = jY-1; //Same
			*/
			//NEW::::::::
			fromX = jX;
			fromY = jY;
			
			//And one more to skip the opening 1
			
			jX = jX+1;
			if(strand == 'f') jY = jY+1; else fromY = jY-1; 
			

			cur_new_len = 1;

			//End of outside while
		}


	}

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
	uint64_t i;
	Bucket * ptr;
	Block * aux;
	Frags_list * flptr;
	uint64_t pre_comp_sb = sizeofSyntenyBlock();
	uint64_t pre_comp_sbl = sizeofSyntenyList();

	Synteny_list * sbl = (Synteny_list *) mp->request_bytes(pre_comp_sbl);
	Synteny_list * curr_sbl = sbl;
	Synteny_block * curr_sb = NULL;
	/* ATTENTION SOMETHING IS NOT RIGHT HERE REDO!! */
	for(i=0;i<ht->get_size();i++){
		ptr = ht->get_key_at(i);
		curr_sbl->next = NULL;
		while(ptr != NULL){
			//Each block is here
			//For each block, add the blocks linked by the fragments
			curr_sbl->sb = NULL;
			flptr = ptr->b.f_list;
			while(flptr != NULL){
				//Each fragment in the current block
				aux = ht->get_block_from_frag(flptr->f);
				if(aux != NULL){
					Synteny_block * aux_sb = (Synteny_block *) mp->request_bytes(pre_comp_sb);
					aux_sb->b = aux; //insert at the head
					aux_sb->next = curr_sb;
					curr_sb = aux_sb;
					printf("Added: "); printBlock(curr_sb->b); getchar();
				}
				flptr = flptr->next;
			}

			// End block
			curr_sbl->sb = curr_sb;
			curr_sbl->next = (Synteny_list *) mp->request_bytes(pre_comp_sbl);
			printf("\t break synteny!!\n");
			curr_sbl = curr_sbl->next;
			ptr = ptr->next;
		}
	}
	return sbl;
}

void has_reversion_in_truple(Bucket * a, Bucket * b, Bucket * c){
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

void detect_evolutionary_event(hash_table * ht, uint64_t max_len_seq){
	/*
	Bucket * a = NULL, * b = NULL, * c = NULL, * d = NULL, * e = NULL;

	uint64_t i = 0;

	//Get first 5 consecutive synteny blocks
	while(a == NULL) {
		a = ht->get_key_at(i);
		if(a == NULL) i++;
	}
	b = a->next;
	while(b == NULL) {
		i++;
		b = ht->get_key_at(i);
	}
	c = b->next;
	while(c == NULL) {
		i++;
		c = ht->get_key_at(i);
	}
	d = c->next;
	while(d == NULL) {
		i++;
		d = ht->get_key_at(i);
	}
	e = d->next;
	while(e == NULL) {
		i++;
		e = ht->get_key_at(i);
	}

	Frags_list * fl;
	printf("Showing consecutive synteny blocks\n");
	printf("\t"); fl = a->b.f_list; printBlock(&a->b); while(fl != NULL){ printFragment(fl->f); fl = fl->next; }
	printf("\t"); fl = b->b.f_list; printBlock(&b->b); while(fl != NULL){ printFragment(fl->f); fl = fl->next; }
	printf("\t"); fl = c->b.f_list; printBlock(&c->b); while(fl != NULL){ printFragment(fl->f); fl = fl->next; }
	printf("\t"); fl = d->b.f_list; printBlock(&d->b); while(fl != NULL){ printFragment(fl->f); fl = fl->next; }
	printf("\t"); fl = e->b.f_list; printBlock(&e->b); while(fl != NULL){ printFragment(fl->f); fl = fl->next; }
	*/

	/*
	while(i<max_len_seq){

		//TODO
		if(a != NULL && b != NULL && c != NULL && d != NULL && e != NULL){
			//Got 5 consecutive synteny blocks
		}
			
	}
	*/
}

