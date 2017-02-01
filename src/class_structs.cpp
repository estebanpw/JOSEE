#include <stdio.h>
#include <cstdlib>
#include <string.h>
#include <ctype.h>
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


hash_table::hash_table(memory_pool * main_mem_pool, uint64_t init_size, sequence_manager * sequences, uint64_t highest_key){
	this->mp = main_mem_pool;
	this->ht_size = init_size;
	this->ht = (Bucket **) this->mp->request_bytes(init_size*sizeof(Bucket *));
	this->sequences = sequences;
	this->n_buckets = init_size;

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
	bkt_x->b.genome = this->sequences->get_sequence_by_label(f->seqX);

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
				//Insert in the fragment list; the block will get inserted on the other coordinate anyway
				//insert_on_list = 1;
				
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
	bkt_y->b.genome = this->sequences->get_sequence_by_label(f->seqY);

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
				//insert_on_list = 1;
				
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
				block_len_verifier = ptr->b.end - ptr->b.start + 1;
				fl = ptr->b.f_list;
				while(fl != NULL){
					fprintf(stdout, "\t"); printFragment(fl->f);
					if(block_len_verifier != fl->f->length) terror("Found different length of fragment in block");
					if(fl->f->strand == 'r') had_reversed = 1;
					fl = fl->next;
				}
				getchar();
			}
			bck_counter++; ptr = ptr->next; 
		}
		if(print >= 1){
			fprintf(stdout, "Entry %"PRIu64" contains %"PRIu64" buckets\n", i, bck_counter);
			//if(had_reversed == 1) getchar();
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

void hash_table::write_blocks_and_breakpoints_to_file(FILE * out_blocks, FILE * out_breakpoints){
	uint64_t i, block_counts = 0;
	Bucket * ptr;
	uint64_t * bps_from = (uint64_t *) std::calloc(this->sequences->get_number_of_sequences(), sizeof(uint64_t));
	uint64_t * bps_to = (uint64_t *) std::calloc(this->sequences->get_number_of_sequences(), sizeof(uint64_t));
	if(bps_from == NULL || bps_to == NULL) terror("Could not allocate breakpoint coordinates");


	unsigned char print_only_noncoding = 1;

	fprintf(out_blocks, "id\tseq\torder\tstart\tend\tlength\n");
	fprintf(out_breakpoints, "id\tseq\tstart\tend\tlength\n");


	for(i=0;i<this->ht_size;i++){
		ptr = this->ht[i];
		while(ptr != NULL){ 
			
			bps_to[ptr->b.genome->id] = ptr->b.start;
			
			//If there actually is a breakpoint
			if(bps_to[ptr->b.genome->id] > bps_from[ptr->b.genome->id]+1){
				
				switch(print_only_noncoding){
					
					case 1:
					{
						//Only if there is no matching gene
						/*
						Annotation * az = binary_search_annotations(bps_from[ptr->b.genome->id]+1, bps_to[ptr->b.genome->id]-1, this->sequences->get_annotation_list(ptr->b.genome->id), this->sequences->get_annotations_number_in_list(ptr->b.genome->id) - 1);
						if(az == NULL) {
							printf("%"PRIu64", %"PRIu64" in %"PRIu64"--> \n", bps_from[ptr->b.genome->id]+1, bps_to[ptr->b.genome->id]-1, ptr->b.genome->id);
							getchar();
						}
						*/
						
						
						

						if(NULL == binary_search_annotations(bps_from[ptr->b.genome->id]+1, bps_to[ptr->b.genome->id]-1, this->sequences->get_annotation_list(ptr->b.genome->id), this->sequences->get_annotations_number_in_list(ptr->b.genome->id) - 1)){
							fprintf(out_breakpoints, "%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\n", 
							block_counts, ptr->b.genome->id, bps_from[ptr->b.genome->id]+1, bps_to[ptr->b.genome->id]-1, bps_to[ptr->b.genome->id] - bps_from[ptr->b.genome->id] + 1);
						}
					}
					break;
					default: fprintf(out_breakpoints, "%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\n", 
							block_counts, ptr->b.genome->id, bps_from[ptr->b.genome->id]+1, bps_to[ptr->b.genome->id]-1, bps_to[ptr->b.genome->id] - bps_from[ptr->b.genome->id] + 1);
				}				
			}
			
			switch(print_only_noncoding){
				case 1:
				{
					if(NULL == binary_search_annotations(ptr->b.start, ptr->b.end, this->sequences->get_annotation_list(ptr->b.genome->id), this->sequences->get_annotations_number_in_list(ptr->b.genome->id)-1)){
						fprintf(out_blocks, "%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\n", 
						block_counts, ptr->b.genome->id, ptr->b.order, ptr->b.start, ptr->b.end, ptr->b.end-ptr->b.start + 1);
					}
				}
				break;
				default: fprintf(out_blocks, "%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\t%"PRIu64"\n", 
			block_counts, ptr->b.genome->id, ptr->b.order, ptr->b.start, ptr->b.end, ptr->b.end-ptr->b.start + 1);

			}
			
			
			
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

sequence_manager::sequence_manager(){
	this->annotation_lists = NULL;
	this->n_annotations = NULL;
	this->n_sequences = 0;
}

uint64_t sequence_manager::load_sequences_descriptors(FILE * lengths_file){


    //Calculate number of sequences according to size of lengths file
    fseeko(lengths_file, 0L, SEEK_END);
    this->n_sequences = ftello(lengths_file)/sizeof(uint64_t);
    fseeko(lengths_file, 0L, SEEK_SET);

    
    //Allocate heap for sequences struct to hold lengths and ids
	this->sequences = (Sequence *) std::malloc(n_sequences*sizeofSequence());


    if(this->sequences == NULL) terror("Could not allocate memory for sequence descriptors");

    //Load sequence data into sequences descriptors
    uint64_t i=0, acum = 0;
    while(i<this->n_sequences){
        this->sequences[i].id = i;
        this->sequences[i].acum = acum;
        if(1 != fread(&this->sequences[i].len, sizeof(uint64_t), 1, lengths_file)) terror("Wrong number of sequences or sequence file corrupted");
        acum += this->sequences[i].len;
        //fprintf(stdout, "[INFO] Sequence %"PRIu64" has length %"PRIu64"\n", i, st[i].len);
        i++;
    }

	return this->n_sequences;
}

uint64_t sequence_manager::get_maximum_length(){
    uint64_t i;
    uint64_t m_len = 0;
    for(i=0;i<this->n_sequences;i++){
        if(m_len < this->get_sequence_by_label(i)->len){
            m_len = this->get_sequence_by_label(i)->len;
        }
    }
    return m_len;
}


void sequence_manager::print_sequences_data(){
    uint64_t i;
    fprintf(stdout, "[INFO] Sequences data:\n");
    for(i=0;i<this->n_sequences;i++){
        fprintf(stdout, "\t(%"PRIu64")\tL:%"PRIu64"\tC:%"PRIu32"\tF:%"PRIu64"\n", i, this->sequences[i].len, this->sequences[i].coverage, this->sequences[i].n_frags);
    }
}

Sequence * sequence_manager::get_sequence_by_label(uint64_t label){
	if(label < this->n_sequences) return &this->sequences[label]; else return NULL;
}

void sequence_manager::read_dna_sequences(char * paths_to_files){
    
    uint64_t i;
    
    FILE * lf = fopen64(paths_to_files, "rt");
    if(lf == NULL) terror("Could not open list of genomic files");


    char ** all_sequences_names = (char **) std::malloc (this->n_sequences*sizeof(char *));
    for(i=0;i<this->n_sequences;i++){
        all_sequences_names[i] = (char *) std::malloc(READLINE*sizeof(char));
        if(all_sequences_names[i] == NULL) terror("Could not allocate paths to files");
    }


    i = 0;
    while(i < this->n_sequences && !feof(lf)){
        if(fgets(all_sequences_names[i], READLINE, lf) > 0){
            if(all_sequences_names[i][0] != '\0' && all_sequences_names[i][0] != '\n'){
                all_sequences_names[i][strlen(all_sequences_names[i])-1] = '\0';
                i++;
            }
        }
    }
    if(i != this->n_sequences) { printf("%"PRIu64"\n", i);terror("Something went wrong. Incorrect number of files"); }

    fclose(lf);
    
    
    
    //Char to hold all sequences
    char ** all_sequences = (char **) std::calloc(this->n_sequences, sizeof(char *));
    if(all_sequences == NULL) terror("Could not allocate sequences pointer");

    //Vector to tell for sequence reallocs
    uint64_t * n_reallocs = (uint64_t *) std::calloc(this->n_sequences, sizeof(uint64_t));
    if(n_reallocs == NULL) terror("Could not allocate realloc count vector");

    //Read using buffered fgetc
    uint64_t idx = 0, r = 0, curr_pos;
    char * temp_seq_buffer = NULL;
    if ((temp_seq_buffer = (char *) std::calloc(READBUF, sizeof(char))) == NULL) {
        terror("Could not allocate memory for read buffer");
    }
    //To force reading from the buffer
    idx = READBUF + 1;
    char c;
    
    
    FILE * current; 

    //Read sequences and load into array
    for(i=0;i<this->n_sequences;i++){
        current = fopen64(all_sequences_names[i], "rt");
        all_sequences[i] = (char *) std::calloc(SEQ_REALLOC, sizeof(char));
        if(all_sequences[i] == NULL) terror("Could not allocate genome sequence");
        if(current == NULL) terror("Could not open fasta file");

        curr_pos = 0;
        idx = READBUF + 1;
        r = 0;

        c = buffered_fgetc(temp_seq_buffer, &idx, &r, current);
        while((!feof(current) || (feof(current) && idx < r))){

            if(c == '>'){
                while(c != '\n') c = buffered_fgetc(temp_seq_buffer, &idx, &r, current); //Skip id

                while(c != '>' && (!feof(current) || (feof(current) && idx < r))){ //Until next id
                    c = buffered_fgetc(temp_seq_buffer, &idx, &r, current);
                    c = toupper(c);
                    if(c >= 'A' && c <= 'Z'){
                        all_sequences[i][curr_pos++] = c;
                        if(curr_pos >= SEQ_REALLOC*n_reallocs[i]){
                            n_reallocs[i]++;
                            all_sequences[i] = (char *) std::realloc(all_sequences[i], n_reallocs[i]*SEQ_REALLOC);
                            if(all_sequences[i] == NULL) terror("Could not realloc sequence");
                        }
                    }
                }
                curr_pos++; //one for the *

            }else{
                c = buffered_fgetc(temp_seq_buffer, &idx, &r, current);
            }
            
        }
        //Realloc final size
        all_sequences[i] = (char *) std::realloc(all_sequences[i], curr_pos);
        if(all_sequences[i] == NULL) terror("Could not realloc sequence");
        this->sequences[i].seq = all_sequences[i]; //Assign the current sequence to its correspondent


        fclose(current);
    }

    std::free(temp_seq_buffer);
    std::free(n_reallocs);

    for(i=0;i<this->n_sequences;i++){
        std::free(all_sequences_names[i]);
    }
    std::free(all_sequences_names);
    std::free(all_sequences); //But not the individual pointers, which are pointed by SEQUENCEs
}

void sequence_manager::read_annotations(){
	//Only if a path was given
	if(this->path_annotations == NULL){
		terror("Requested to load annotations but no annotation file was given");
	}

	//Open file
	FILE * gbconcat = fopen64(this->path_annotations, "rt");
	if(gbconcat == NULL) terror("Could not open annotations file");

	//Allocate in memory pool space for annotations
	this->annotation_lists = (Annotation **) std::malloc(this->n_sequences*sizeof(Annotation *));
	this->n_annotations = (uint64_t *) std::calloc(this->n_sequences, sizeof(uint64_t));
	if(this->n_annotations == NULL) terror("Could not allocate annotations number");
	uint64_t * n_reallocs = (uint64_t *) std::malloc(this->n_sequences * sizeof(uint64_t));
	if(n_reallocs == NULL) terror("Could not allocate realloc annotations number");

	//Buffer line
	char line[READLINE], nullString[READLINE];
	if(0 == fgets(line, READLINE, gbconcat)) terror("Could not read annotations file first line");

	//Current label for genome
	int64_t current_label = -1;

	//Temporary annotations
	Annotation an1, an2;

	//Read gene annotations file
	while(!feof(gbconcat)){
		while(!strncmp(line, "VERSION", 7) == 0){ //Until finding a version
			if(0 == fgets(line, READLINE, gbconcat)) break;
		}
		current_label++; //Annotation file corresponding to label

		this->annotation_lists[current_label] = (Annotation *) std::malloc(INIT_ANNOTS*sizeofAnnotation());
		if(this->annotation_lists[current_label] == NULL) terror("Could not allocate annotation sublists");
		n_reallocs[current_label] = 1;

		if(0 == fgets(line, READLINE, gbconcat)) break;

		//Now until finding another "VERSION"
		while(!strncmp(line, "VERSION", 7) == 0){
			if(0 == fgets(line, READLINE, gbconcat)) break;

			if(strncmp(line, "     gene", 9) == 0){
				//Fill gene positions
				if(strstr(line, "join") != NULL){
					//gene            join(839406..839615,1..1215)
					sscanf(line, "%[^(](%"PRIu64"%[.>]%"PRIu64",%"PRIu64"%[.>]%"PRIu64")", nullString, &an1.start, nullString, &an1.end, &an2.start, nullString, &an2.end);
					an1.strand = 'f';
					an2.strand = 'f';
					an1.product = NULL;
					an2.product = NULL;

					//See if it still has space to add the two annotations
					if(this->n_annotations[current_label]+2 >= n_reallocs[current_label]*INIT_ANNOTS){
						n_reallocs[current_label]++;
						this->annotation_lists[current_label] = (Annotation *) std::realloc(this->annotation_lists[current_label], n_reallocs[current_label]*INIT_ANNOTS);
						if(this->annotation_lists == NULL) terror("Could not realloc list of annotations");
					}
					memcpy(&this->annotation_lists[current_label][this->n_annotations[current_label]], &an1, sizeofAnnotation()); 
					this->n_annotations[current_label]++;
					
					memcpy(&this->annotation_lists[current_label][this->n_annotations[current_label]], &an2, sizeofAnnotation()); 
					this->n_annotations[current_label]++;
					

				}else if(strstr(line, "complement") != NULL){
					//Its complemented
					//     gene            complement(16694..16957)
				
					sscanf(line, "%[^(](%"PRIu64"%[.>]%"PRIu64")", nullString, &an1.start, nullString, &an1.end);
					an1.strand = 'r';
					an1.product = NULL;
					//See if it still has space to add the annotation
					if(this->n_annotations[current_label] == n_reallocs[current_label]*INIT_ANNOTS){
						n_reallocs[current_label]++;
						this->annotation_lists[current_label] = (Annotation *) std::realloc(this->annotation_lists[current_label], n_reallocs[current_label]*INIT_ANNOTS);
						if(this->annotation_lists == NULL) terror("Could not realloc list of annotations");
					}
					memcpy(&this->annotation_lists[current_label][this->n_annotations[current_label]], &an1, sizeofAnnotation()); 
					this->n_annotations[current_label]++;


				}else{
					//Straight
					//     gene            1..1392
					sscanf(line, "%s%"PRIu64"%[.>]%"PRIu64"", nullString, &an1.start, nullString, &an1.end);
					an1.strand = 'f';
					
					//See if it still has space to add the annotation
					if(this->n_annotations[current_label] == n_reallocs[current_label]*INIT_ANNOTS){
						n_reallocs[current_label]++;
						this->annotation_lists[current_label] = (Annotation *) std::realloc(this->annotation_lists[current_label], n_reallocs[current_label]*INIT_ANNOTS);
						if(this->annotation_lists == NULL) terror("Could not realloc list of annotations");
					}
					memcpy(&this->annotation_lists[current_label][this->n_annotations[current_label]], &an1, sizeofAnnotation()); 
					this->n_annotations[current_label]++;

				}

			}
		}
		quick_sort_annotations(this->annotation_lists[current_label], 0, this->n_annotations[current_label]-1);
	}

	std::free(n_reallocs);
	fclose(gbconcat);
}

void sequence_manager::print_annotations(){
	uint64_t i, j;
	for(i=0;i<this->n_sequences;i++){
		for(j=0;j<this->n_annotations[i];j++){
			printAnnotation(&this->annotation_lists[i][j]);
		}
	}
}

sequence_manager::~sequence_manager(){
	uint64_t i;
	for(i=0;i<this->n_sequences;i++){
		if(this->sequences[i].seq != NULL) std::free(this->sequences[i].seq);
	}
	if(this->n_annotations != NULL) std::free(this->n_annotations);
	if(this->annotation_lists != NULL){
		for(i=0;i<this->n_sequences;i++){
			std::free(this->annotation_lists[i]);
		}
		std::free(this->annotation_lists);
	}
}

dictionary_hash::dictionary_hash(uint64_t init_size, uint64_t highest_key, uint32_t kmer_size){
	this->words = (Wordbucket *) std::calloc(init_size, sizeofWordbucket());
	this->ht_size = init_size;
	this->kmer_size = kmer_size;

	//Just in case the init size is larger than the highest genome
	if(init_size < highest_key){
		this->key_factor = (double)(init_size)/(highest_key);
	}else{
		this->key_factor = 1.0;
	}
}


uint64_t dictionary_hash::compute_hash(char * kmer){
	return hashOfWord(kmer, this->kmer_size);
}

dictionary_hash::~dictionary_hash(){
	std::free(this->words);
}