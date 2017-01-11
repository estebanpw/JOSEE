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
			break;
		}

		//If the block is not equal but a block for this genome is already contained in the bucket
		//Then we should insert ordered
		//if(bkt_x->b.genome == ptr->b.genome){
			if(bkt_x->b.start > ptr->b.start){ //Only if its bigger so that it keeps the reference to the previous
				theoretical_position = ptr;
			}
		//}

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
		ht[hash_x]->b.f_list = frag_pointer;
		ht[hash_x]->b.f_list->next = NULL; 
		ht[hash_x]->b.f_list->f = f;
		
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
		if(ht[hash_y]->b.start == 52 && ht[hash_y]->b.end == 209 && f->yStart == 85 && f->yEnd == 186){
			printf("\nresult: %d\n", isBlockEqualTo(&bkt_y->b, &ptr->b));
			printf("Trying to insert : "); printBlock(&bkt_y->b);
			printf("Compared to      : "); printBlock(&ptr->b);

		}
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
			break;
		}

		//If the block is not equal but a block for this genome is already contained in the bucket
		//Then we should insert ordered
		//if(bkt_y->b.genome == ptr->b.genome){
			if(bkt_y->b.start > ptr->b.start){
				theoretical_position = ptr;
			}
		//}

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
		
		if(ht[hash_y]->b.start == 52 && ht[hash_y]->b.end == 209 && f->yStart == 85 && f->yEnd == 186){
			printf("im inserting it!!"); //BUg here
		}

		//Insert frag into list
		ht[hash_y]->b.f_list = frag_pointer;
		ht[hash_y]->b.f_list->next = NULL; 
		ht[hash_y]->b.f_list->f = f;

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
	uint64_t i, bck_counter, total_buckets = 0;
	Bucket * ptr;
	Frags_list * fl;
	for(i=0;i<this->ht_size;i++){
		bck_counter = 0;
		ptr = this->ht[i];
		while(ptr != NULL){ 
			if(print == 2){
				printBlock(&ptr->b);
				fl = ptr->b.f_list;
				while(fl != NULL){
					fprintf(stdout, "\t"); printFragment(fl->f);
					if(fl->f->strand == 'r') getchar();
					fl = fl->next;
				}
				getchar();
			}
			bck_counter++; ptr = ptr->next; 
		}
		if(print >= 1){
			fprintf(stdout, "Entry %"PRIu64" contains %"PRIu64" buckets\n", i, bck_counter);
			getchar();
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




