#include <stdio.h>
#include <cstdlib>
#include "structs.h"
#include "commonFunctions.h"
#include "comparisonFunctions.h"

memory_pool::memory_pool(uint64_t max_pools)
{
	this->current_pool = 0;
	this->max_pools = max_pools;
	this->mem_pool = (char **) std::calloc(max_pools, sizeof(char *));
	this->base = (uint64_t *) std::malloc(max_pools * sizeof(uint64_t));
	this->base[0] = 0;
	if (this->mem_pool == NULL) terror("Could not allocate memory pools");
	this->mem_pool[0] = (char *) std::calloc(POOL_SIZE, sizeof(char));
	if (this->mem_pool[0] == NULL) terror("Could not allocate initial memory pool");
}

void * memory_pool::request_bytes(uint64_t n_bytes)
{
	void * ptr;
	if (this->base[this->current_pool] + n_bytes >= POOL_SIZE) {
		this->current_pool++;
		if(this->current_pool == this->max_pools) terror("Reached maximum number of pools. Exiting.");
		this->mem_pool[this->current_pool] = (char *) std::calloc(POOL_SIZE, sizeof(char));
		if (this->mem_pool[this->current_pool] == NULL) terror("Could not allocate memory pool");
		this->base[this->current_pool] = 0;
	}
	
	ptr = &this->mem_pool[this->current_pool][0] + this->base[this->current_pool];
	this->base[this->current_pool] = this->base[this->current_pool] + n_bytes;
	
	return ptr;
}

void memory_pool::reset_n_bytes(uint64_t bytes){
	//Makes no checks, assuming you allocated some a priori
	if(bytes >= this->base[current_pool]){
		this->base[current_pool] = this->base[current_pool] - bytes;
	}
}

memory_pool::~memory_pool()
{
	for (uint64_t i = 0; i <= this->current_pool; i++) {
		free(this->mem_pool[i]);
	}
	free(this->mem_pool);
	free(this->base);
}


hash_table::hash_table(memory_pool * main_mem_pool, uint64_t init_size, Sequence * sequences, uint64_t highest_key){
	this->mp = main_mem_pool;
	this->ht_size = init_size;
	this->ht = (Bucket **) this->mp->request_bytes(init_size*sizeof(Bucket *));
	this->sequences = sequences;

	uint64_t i;
	for(i=0;i<init_size;i++) this->ht[i] = NULL;

    computed_sizeof_block = sizeofBucket(); //Avoid overcomputing
    computed_sizeof_frags_list = sizeofFrags_list(); //Avoid overcomputing

	//Just in case the init size is larger than the highest genome
	if(init_size < highest_key){
		this->key_factor = (double)(init_size)/(highest_key);
	}else{
		this->key_factor = 1.0;
	}
	

}

uint64_t hash_table::compute_hash(uint64_t key){
	return (uint64_t) (key_factor * key); //Partitionates the table from 0 to max genome length
}

void hash_table::insert_block(struct FragFile * f){
	
	this->insert_x_side(f);
	this->insert_y_side(f);

}

void hash_table::insert_x_side(struct FragFile * f){
	uint64_t hash_x = compute_hash(f->xStart);
	

	//Condition to insert in the frags list
	int insert_on_list = 0;
	int exit = 0;

	//Get memory for buckets
	Bucket * bkt_x = (Bucket *) this->mp->request_bytes(this->computed_sizeof_block);

	//Fill data of block
	bkt_x->next = NULL;
	bkt_x->b.start = f->xStart;
	bkt_x->b.end = f->xEnd;
	bkt_x->b.order = 0; //Order will be assigned later
	bkt_x->b.genome = &this->sequences[f->seqX];

	//Insertions
	Bucket * ptr = ht[hash_x];
	Bucket * theoretical_position = NULL;

	
	while(ptr != NULL){
		if(isBlockEqualTo(&bkt_x->b, &ptr->b)){
			
			this->mp->reset_n_bytes(this->computed_sizeof_block); //First reset the bytes taken for the block
			if(idNotInList(ptr->b.f_list, f)){
				//The block exists but not linked to this fragment, so add it to the list

				insert_on_list = 1;	
			}else{
				//If the block already exists for this genome and for this fragment then it is a repetition
				//(Only varies its y-coordinates)
				//What do here?
				
			}
			//Exit since the block exists
			exit = 1;
		}

		//If the block is not equal but a block for this genome is already contained in the bucket
		//Then we should insert ordered
		//if(bkt_x->b.genome == ptr->b.genome){
		if(bkt_x->b.start > ptr->b.start){ //Only if its bigger so that it keeps the reference to the previous
			theoretical_position = ptr;
		}
		//}
		if(exit == 1) break;
		ptr = ptr->next;	
	}

	//Actual insertion: If null pointer then the block was not contained in the set
	//If not null pointer, it was already contained and thus we have to reset the bytes requested
	if(ptr == NULL){
		
		if(ht[hash_x] == NULL) this->n_entries++; 

		Frags_list * frag_pointer = (Frags_list *) this->mp->request_bytes(this->computed_sizeof_frags_list);
		
		//Insert between theoretical position and its next
		
		if(theoretical_position == NULL){
			bkt_x->next = ht[hash_x];  //Insert at the head
			ht[hash_x] = bkt_x;
		}else{
			bkt_x->next = theoretical_position->next;
			theoretical_position->next = bkt_x;
		}
		


		//Insert frag into list
		bkt_x->b.f_list = frag_pointer;
		bkt_x->b.f_list->next = NULL; 
		bkt_x->b.f_list->f = f;
		
		bkt_x->b.present_in_synteny = 0;

		this->n_buckets++;

	}

	if(ptr != NULL && insert_on_list == 1){
		
		Frags_list * frag_pointer = (Frags_list *) this->mp->request_bytes(this->computed_sizeof_frags_list);
		frag_pointer->next = ptr->b.f_list;
		frag_pointer->f = f;
		ptr->b.f_list = frag_pointer;


	}
	

}

void hash_table::insert_y_side(struct FragFile * f){
	uint64_t hash_y = compute_hash(f->yStart);
		
	//Condition to insert in the frags list
	int insert_on_list = 0;
	int exit = 0;

	//Get memory for buckets
	Bucket * bkt_y = (Bucket *) this->mp->request_bytes(this->computed_sizeof_block);

	//Fill data of block
	bkt_y->next = NULL;
	bkt_y->b.start = f->yStart;
	bkt_y->b.end = f->yEnd;
	bkt_y->b.order = 0; //Order will be assigned later
	bkt_y->b.genome = &this->sequences[f->seqY];

	//Insertions
	Bucket * ptr = ht[hash_y];
	Bucket * theoretical_position = NULL;


	while(ptr != NULL){
		
		if(isBlockEqualTo(&bkt_y->b, &ptr->b)){

			this->mp->reset_n_bytes(this->computed_sizeof_block); //First reset the bytes taken for the block
			if(idNotInList(ptr->b.f_list, f)){
				//The block exists but not linked to this fragment, so add it to the list
				insert_on_list = 1;	
			}else{
				//If the block already exists for this genome and for this fragment then it is a repetition
				//(Only varies its y-coordinates)
				//What do here?
				
			}
			//Exit since the block exists
			exit = 1;
		}

		//If the block is not equal but a block for this genome is already contained in the bucket
		//Then we should insert ordered
		//if(bkt_y->b.genome == ptr->b.genome){
		if(bkt_y->b.start > ptr->b.start){
			theoretical_position = ptr;
		}
		//}
		if(exit == 1) break;

		ptr = ptr->next;	
	}

	//Actual insertion: If null pointer then the block was not contained in the set
	//If not null pointer, it was already contained and thus we have to reset the bytes requested
	if(ptr == NULL){
		
		if(ht[hash_y] == NULL) this->n_entries++; 

		Frags_list * frag_pointer = (Frags_list *) this->mp->request_bytes(this->computed_sizeof_frags_list);
		
		//Insert between theoretical position and its next
		if(theoretical_position == NULL){
			bkt_y->next = ht[hash_y];  //Insert at the head
			ht[hash_y] = bkt_y;
		}else{
			bkt_y->next = theoretical_position->next;
			theoretical_position->next = bkt_y;
		}

		

		//Insert frag into list
		bkt_y->b.f_list = frag_pointer;
		bkt_y->b.f_list->next = NULL; 
		bkt_y->b.f_list->f = f;

		bkt_y->b.present_in_synteny = 0;

		this->n_buckets++;

	}

	if(ptr != NULL && insert_on_list == 1){
		
		Frags_list * frag_pointer = (Frags_list *) this->mp->request_bytes(this->computed_sizeof_frags_list);
		frag_pointer->next = ptr->b.f_list;
		frag_pointer->f = f;
		ptr->b.f_list = frag_pointer;	

	}
	
}

void hash_table::print_hash_table(int print){
	uint64_t i, bck_counter, total_buckets = 0, block_len_verifier;
	Bucket * ptr;
	Frags_list * fl;
	int had_reversed = 0;
	for(i=0;i<this->ht_size;i++){
		bck_counter = 0;
		ptr = this->ht[i];
		had_reversed = 0;
		while(ptr != NULL){ 
			if(print == 2){
				printBlock(&ptr->b);
				had_reversed = 0;
				block_len_verifier = ptr->b.end - ptr->b.start;
				fl = ptr->b.f_list;
				while(fl != NULL){
					fprintf(stdout, "\t"); printFragment(fl->f);
					if(block_len_verifier != fl->f->length) terror("Found different length of fragment in block");
					if(fl->f->strand == 'r') had_reversed = 1;
					fl = fl->next;
				}
				//getchar();
			}
			bck_counter++; ptr = ptr->next; 
		}
		if(print >= 1){
			fprintf(stdout, "Entry %"PRIu64" contains %"PRIu64" buckets\n", i, bck_counter);
			if(had_reversed == 1) getchar();
		}
		total_buckets += bck_counter;
	}
	fprintf(stdout, "%"PRIu64" buckets.\n", total_buckets);
}

Bucket * hash_table::get_value_at(uint64_t pos){
	Bucket * ptr = this->get_key_at(compute_hash(pos));
	while(ptr != NULL && ptr->b.start != pos) ptr = ptr->next;

	return ptr;
}

Block * hash_table::get_block_from_frag(struct FragFile * f, int x_or_y){
	Bucket * ptr = this->get_key_at(compute_hash(f->xStart));
	
	while(ptr != NULL){
		if(x_or_y == 0 && ptr->b.start == f->xStart && ptr->b.end == f->xEnd
		&& ptr->b.genome->id == f->seqX) { return &ptr->b;} 
		if(x_or_y == 1 && ptr->b.start == f->yStart && ptr->b.end == f->yEnd
		&& ptr->b.genome->id == f->seqY) { return &ptr->b;} 

		ptr = ptr->next;
	}
	return NULL;
}

void hash_table::write_blocks_and_breakpoints_to_file(FILE * out_blocks, FILE * out_breakpoints, uint64_t n_sequences){
	uint64_t i, block_counts = 0;
	Bucket * ptr;
	uint64_t * bps_from = (uint64_t *) std::calloc(n_sequences, sizeof(uint64_t));
	uint64_t * bps_to = (uint64_t *) std::calloc(n_sequences, sizeof(uint64_t));
	if(bps_from == NULL || bps_to == NULL) terror("Could not allocate breakpoint coordinates");

	fprintf(out_blocks, "id\tseq\torder\tstart\tend\tlength\n");
	fprintf(out_breakpoints, "id\tseq\tstart\tend\tlength\n");
	for(i=0;i<this->ht_size;i++){
		ptr = this->ht[i];
		while(ptr != NULL){ 
			
			bps_to[ptr->b.genome->id] = ptr->b.start;
			

			fprintf(out_breakpoints, "%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\n", 
			block_counts, ptr->b.genome->id, bps_from[ptr->b.genome->id], bps_to[ptr->b.genome->id], bps_to[ptr->b.genome->id] - bps_from[ptr->b.genome->id]);
			
			fprintf(out_blocks, "%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\n", 
			block_counts, ptr->b.genome->id, ptr->b.order, ptr->b.start, ptr->b.end, ptr->b.end-ptr->b.start);
			
			bps_from[ptr->b.genome->id] = ptr->b.end;

			block_counts++;
			ptr = ptr->next; 
		}
	}
	std::free(bps_from);
	std::free(bps_to);
}


strand_matrix::strand_matrix(uint64_t sequences){
	uint64_t i;
	this->n_seqs = sequences;
	this->squared_sequences = sequences * sequences;
	this->sm = (unsigned char **) std::calloc(sequences, sizeof(unsigned char *));
	if(this->sm == NULL) terror("Could not allocate strand matrix");
	for(i=0;i<sequences;i++){
		this->sm[i] = (unsigned char *) std::calloc(sequences, sizeof(unsigned char));
		if(this->sm[i] == NULL) terror("Could not allocate strand matrix subdimensions");
	}
}

void strand_matrix::add_fragment_strands(Synteny_list * sbl){
	Synteny_block * sb_ptr;
	Frags_list * fl;
	if(sbl != NULL){
		sb_ptr = sbl->sb;
		//printf("A block...\n");
		while(sb_ptr != NULL){

			fl = sb_ptr->b->f_list;
			while(fl != NULL){
				//printf("A frag...\n");
				//TODO
				//Check if the fragment has to be added or not to the strand matrix
				//I think all fragments should be added and that the current issues are bugs
				if(1==1){ //REPLACE: has_to_be_added(sb_ptr->b->)
					(fl->f->strand == 'f') ? this->sm[fl->f->seqX][fl->f->seqY] = FORWARD : this->sm[fl->f->seqX][fl->f->seqY] = REVERSE;
				}

				fl = fl->next;
			}
			

			sb_ptr = sb_ptr->next;
		}
	}
	//printf("==========================================\n");
}

void strand_matrix::print_strand_matrix(){
	uint64_t i,j;
	for(i=0;i<n_seqs;i++){
		for(j=0;j<n_seqs;j++){
			fprintf(stdout, "%u\t", sm[i][j]);
		}
		fprintf(stdout, "\n");
	}
}

int strand_matrix::is_block_reversed(uint64_t block_number){
	uint64_t i, n_rev = 0, n_for = 0;
	for(i=0; i<block_number; i++){
		if(this->sm[block_number][i] == FORWARD) n_for++;
		if(this->sm[block_number][i] == REVERSE) n_rev++;
	}
	for(i=block_number+1;i<this->n_seqs;i++){
		if(this->sm[i][block_number] == FORWARD) n_for++;
		if(this->sm[i][block_number] == REVERSE) n_rev++;
	}

	//Now we have to check in case that some are forward because the others are reversed
	

	return (n_rev > n_for) ? 1 : -1;
}

strand_matrix::~strand_matrix(){
	uint64_t i;
	for (i = 0; i < this->n_seqs; i++) {
		std::free(this->sm[i]);
	}
	std::free(this->sm);
}
