#ifndef STRUCTS_H
#define STRUCTS_H

#include <inttypes.h>

#pragma pack(push, 1)

#define READLINE 2000
#define READBUF 50000000 //50MB
#define INIT_TRIM_FRAGS 10000

#define MAX_MEM_POOLS 256
#define POOL_SIZE 1024*1024*256 //256 MB

#define NOFRAG 0
#define OPENFRAG 1
#define COVERFRAG 2
#define CLOSEFRAG 3

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

typedef struct sequence{
    uint64_t id;    //Label of the sequence
    uint64_t len;   //Length in nucleotides of the sequence
    uint64_t acum;  //Accumulated length from the sequences found before in the file (if any)
} Sequence;

//A block that belongs to a genome and that has some synteny level (conserved block)
typedef struct block{
    uint64_t start;     //Starting coordinate
    uint64_t end;       //Ending coordinate
    uint64_t order;     //Order of block according to the genome
    uint64_t synteny_level; //Number of other genomes where this block is present
    Sequence * genome;    //A pointer to the genome to which it belongs
} Block;




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

typedef struct frags_list{
    struct FragFile * f;
    struct frags_list * next;
} Frags_list;


typedef struct bucket {
	Block b;
    Frags_list * f_list;
	struct bucket * next;
} Bucket;

//Hash-table class for Blocks
class hash_table
{

private:
	Bucket * ht; // the table itself
	memory_pool * mp;
    uint64_t ht_size; //Size for the hash table
    double key_factor; //To partitionate the space by the largest genome
    uint64_t computed_sizeof_block; //Avoid overcomputing
    uint64_t computed_sizeof_frags_list; //Avoid overcomputing
    Sequence * sequences;

public:
    hash_table(memory_pool * main_mem_pool, uint64_t init_size, Sequence * sequences, uint64_t highest_key);
	//Bucket * getBucketAt(uint64_t index);
	//Chunk * getChunkByKey(const Vec3GLui& key);
	void insert_block(struct FragFile * f);
    Bucket * keys_iterator(){ return &ht[0]; }
	~hash_table();

private:
	uint64_t compute_hash(uint64_t key);
    void insert_x_side(struct FragFile * f);
    void insert_y_side(struct FragFile * f);
};

#endif
