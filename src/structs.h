#ifndef STRUCTS_H
#define STRUCTS_H

#include <inttypes.h>


#pragma pack(push, 1)

#define SEQ_REALLOC 5000000
#define INIT_SEQS 20
#define READLINE 2000
#define READBUF 50000000 //50MB
#define INIT_TRIM_FRAGS 10000

#define MAX_MEM_POOLS 256
#define POOL_SIZE 1024*1024*256 //256 MB

#define NOFRAG 0
#define OPENFRAG 1
#define COVERFRAG 2
#define CLOSEFRAG 3

#define FORWARD 1
#define REVERSE 2

class memory_pool;
class hash_table;

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
    char * seq; //DNA sequence
} Sequence;

typedef struct frags_list{
    struct FragFile * f;
    struct frags_list * next;
} Frags_list;

//A block that belongs to a genome and that has some synteny level (conserved block)
typedef struct block{
    uint64_t start;     //Starting coordinate
    uint64_t end;       //Ending coordinate
    uint64_t order;     //Order of block according to the genome
    Frags_list * f_list;    //List of fragments that compose it
    Sequence * genome;      //A pointer to the genome to which it belongs
    unsigned char present_in_synteny;   //To tell whether it has already been used in a synteny block
    unsigned char strand_in_synteny;    //The strand that it has at the synteny block
} Block;

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

public:
    memory_pool(uint64_t max_pools);
    void * request_bytes(uint64_t bytes);
    void reset_n_bytes(uint64_t bytes);
    ~memory_pool();
};


typedef struct bucket {
	Block b;
	struct bucket * next;
} Bucket;

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
    Sequence * sequences;

public:
    hash_table(memory_pool * main_mem_pool, uint64_t init_size, Sequence * sequences, uint64_t highest_key);
	void insert_block(struct FragFile * f);
    Bucket * get_key_at(uint64_t pos){ if(pos < ht_size && pos >= 0) return ht[pos]; else return NULL; } //Returns a reference to the key by absolute position
    Bucket * get_value_at(uint64_t pos); //Returns a reference to the key computed from the hash of x_pos
    Block * get_block_from_frag(struct FragFile * f, int x_or_y);
    //double get_load_factor(){ return (double)ht_size/n_buckets;}
    uint64_t get_size(){ return ht_size; }
    void print_hash_table(int print);
    Bucket * get_iterator(){ return ht[0];}
    void write_blocks_and_breakpoints_to_file(FILE * out_blocks, FILE * out_breakpoints, uint64_t n_sequences);

private:
	uint64_t compute_hash(uint64_t key);
    void insert_x_side(struct FragFile * f);
    void insert_y_side(struct FragFile * f);
};

// There will be one strand matrix per synteny block
class strand_matrix
{

public:
    unsigned char ** sm;
    uint64_t n_seqs;
    uint64_t squared_sequences;

    strand_matrix(uint64_t sequences);
    int is_block_reversed(uint64_t block_number);
    void reset() { for(uint64_t i=0;i<n_seqs;i++){ for(uint64_t j=0;j<n_seqs;j++){ this->sm[i][j] = 0; } } }
    void add_fragment_strands(Synteny_list * sbl);
    void print_strand_matrix();
    ~strand_matrix();
};

#endif
