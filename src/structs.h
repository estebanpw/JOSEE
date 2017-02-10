#ifndef STRUCTS_H
#define STRUCTS_H

#include <inttypes.h>


#pragma pack(push, 1)

#define SEQ_REALLOC 5000000
#define INIT_SEQS 20
#define INIT_ANNOTS 5000
#define READLINE 2000
#define READBUF 50000000 //50MB
#define INIT_TRIM_FRAGS 10000
#define TABLE_RATE 100 //hash table lengh divisor
#define INIT_CANDIDATES_ALIGN 100 //For wordbucket list

#define MAX_MEM_POOLS 256
#define POOL_SIZE 1024*1024*128 //128 MB

#define NOFRAG 0
#define OPENFRAG 1
#define COVERFRAG 2
#define CLOSEFRAG 3

#define FORWARD 1
#define REVERSE 2
#define MIXED 3

#define POINT 4

class memory_pool;
class hash_table;
class sequence_manager;

//Struct for FragHits, af2png and leeFrag programs
struct FragFile {
    //Diagonal where the frag is located
    //This value is calculated as:
    //posX - posY
    int64_t diag;
    //Start position in sequence X
    uint64_t xStart;
    //Start position in Sequence Y
    uint64_t yStart;
    //End position in Sequence X
    uint64_t xEnd;
    //End position in Sequence Y
    uint64_t yEnd;
    //Fragment Length
    //For ungaped aligment is:
    //xEnd-xStart+1
    uint64_t length;
    //Number of identities in the
    //fragment
    uint64_t ident;
    //Score of the fragment. This
    //depends on the score matrix
    //used
    uint64_t score;
    //Percentage of similarity. This
    //is calculated as score/scoreMax
    //Where score max is the maximum
    //score possible
    float similarity;
    //sequence number in the 'X' file
    uint64_t seqX;
    //sequence number in the 'Y' file
    uint64_t seqY;
    //synteny block id
    int64_t block;
    //'f' for the forward strain and 'r' for the reverse
    char strand;
    //E-value of fragment
    long double evalue;
};

//Sequence descriptor
typedef struct sequence{
    uint64_t id;    //Label of the sequence
    uint64_t len;   //Length in nucleotides of the sequence
    uint64_t acum;  //Accumulated length from the sequences found before in the file (if any)
    uint32_t coverage; //The percentage of bases covered by fragments of a minimum length
    uint64_t n_frags; //Number of fragments that the sequence had
    char * seq; //DNA sequence
} Sequence;

typedef struct quickfrag{
    uint64_t x_start;
    uint64_t y_start;
    uint64_t t_len;
    int64_t diag;
    long double sim;
    Sequence * x;
    Sequence * y;
} Quickfrag;

typedef struct frags_list{
    struct FragFile * f;
    struct frags_list * next;
} Frags_list;

typedef struct linked_list_pos{
    uint64_t pos;
    struct linked_list_pos * next;
} llpos;

//A block that belongs to a genome and that has some synteny level (conserved block)
typedef struct block{
    uint64_t start;     //Starting coordinate
    uint64_t end;       //Ending coordinate
    uint64_t order;     //Order of block according to the genome
    Frags_list * f_list;    //List of fragments that compose it
    Sequence * genome;      //A pointer to the genome to which it belongs
    unsigned char present_in_synteny;   //To tell whether it has already been used in a synteny block
    unsigned char strand_in_synteny;    //The strand that it has at the synteny block
    //unsigned char repetition;   //To tell if the block is a repetition or not
} Block;

typedef struct word{
    uint64_t hash;
    uint64_t pos;
    char strand;
    Block * b;
} Word;

//A synteny block is a collection of blocks
typedef struct synteny_block{
    Block * b;
    struct synteny_block * next;
} Synteny_block;

typedef struct synteny_list{
    Synteny_block * sb;
    uint64_t synteny_level;
    struct synteny_list * next;
} Synteny_list;


//Class for allocating memory only once and requesting particular amounts of bytes
class memory_pool{

private:
    char ** mem_pool;
    uint64_t * base;
    uint64_t current_pool;
    uint64_t max_pools;
    uint64_t pool_size;

public:
    memory_pool(uint64_t max_pools, uint64_t pool_size);
    void * request_bytes(uint64_t bytes);
    void reset_n_bytes(uint64_t bytes);
    void reset_to(uint64_t pool, uint64_t position){ this->current_pool = 0; this->base[current_pool] = 0;}
    void full_reset();
    ~memory_pool();
};


typedef struct bucket {
	Block b;
	struct bucket * next;
} Bucket;

typedef struct wordbucket{
    Word w;
    struct wordbucket * next;
} Wordbucket;



typedef struct annotation{
    uint64_t start;
    uint64_t end;
    char strand;
    char * product; //Requires hard copy
} Annotation;


class sequence_manager
{
private:
    Sequence * sequences;   //A pointer to the sequences
    uint64_t n_sequences;   //Number of sequences
    char * path_annotations;
    Annotation ** annotation_lists;
    uint64_t * n_annotations;

public:
    sequence_manager();
    void set_path_annotations(char * p){ path_annotations = p; }
    char * get_path_annotations(){ return this->path_annotations; }
    uint64_t load_sequences_descriptors(FILE * lengths_file);
    Sequence * get_sequence_by_label(uint64_t label);
    uint64_t get_maximum_length();
    uint64_t get_number_of_sequences() { return n_sequences; }
    void print_sequences_data();
    void print_annotations();
    Annotation * get_annotation_list(uint64_t label){ return this->annotation_lists[label]; }
    uint64_t get_annotations_number_in_list(uint64_t label){ return this->n_annotations[label]; }
    void read_dna_sequences(char * paths_to_files);
    void read_annotations();
    void print_sequence_region(uint64_t label, uint64_t from, uint64_t to);
    ~sequence_manager();
};

//Hash-table class for Blocks
class hash_table
{

private:
	Bucket ** ht; // the table itself
	memory_pool * mp;
    uint64_t ht_size; //Size for the hash table //Init size
    uint64_t n_buckets; //To compute the load factor
    uint64_t n_entries; //Used up entries
    double key_factor; //To partitionate the space by the largest genome
    uint64_t computed_sizeof_block; //Avoid overcomputing
    uint64_t computed_sizeof_frags_list; //Avoid overcomputing
    sequence_manager * sequences;

public:
    hash_table(memory_pool * main_mem_pool, uint64_t init_size, sequence_manager * sequences, uint64_t highest_key);
	void insert_block(struct FragFile * f);
    Bucket * get_key_at(uint64_t pos){ if(pos < ht_size && pos >= 0) return ht[pos]; else return NULL; } //Returns a reference to the key by absolute position
    Bucket * get_value_at(uint64_t pos); //Returns a reference to the key computed from the hash of x_pos
    Block * get_block_from_frag(struct FragFile * f, int x_or_y);
    //double get_load_factor(){ return (double)ht_size/n_buckets;}
    uint64_t get_size(){ return ht_size; }
    void print_hash_table(int print);
    Bucket * get_iterator(){ return ht[0];}
    void write_blocks_and_breakpoints_to_file(FILE * out_blocks, FILE * out_breakpoints);

private:
	uint64_t compute_hash(uint64_t key);
    void insert_x_side(struct FragFile * f);
    void insert_y_side(struct FragFile * f);
};

class dictionary_hash{
private:
    Wordbucket ** words;
    Wordbucket ** list; //To retrieve all possible sequences to align with
    uint64_t list_allocs; //how many times the list was reallocated
    uint64_t n_list_pointers; //The number of pointers in the current list
    uint64_t ht_size;
    uint32_t kmer_size;
    uint64_t computed_sizeofwordbucket;
    double key_factor;
    memory_pool * mp;
public:
    dictionary_hash(uint64_t init_size, uint64_t highest_key, uint32_t kmer_size);
    Wordbucket ** put_and_hit(char * kmer, char strand, uint64_t position, Block * b);
    uint64_t get_candidates_number(){ return this->n_list_pointers;}
    void clear();
    ~dictionary_hash();
private:
    uint64_t compute_hash(char * kmer);
};

typedef struct slist{
    Sequence * s;
    struct slist * next;
} Slist;



// There will be one strand matrix per synteny block
class strand_matrix
{

public:
    unsigned char ** sm;
    uint64_t n_seqs;
    uint64_t squared_sequences;
    uint64_t acu_frags_forward;
    uint64_t acu_frags_reverse;

    strand_matrix(uint64_t sequences);
    //Unsure about this one below
    int is_block_reversed(uint64_t block_number);
    //
    int get_strands_type();
    void reset();
    void add_fragment_strands(Synteny_list * sbl);
    void print_strand_matrix();
    ~strand_matrix();
};

#endif
