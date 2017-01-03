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

//A block that belongs to a genome and that has some synteny level (conserved block)
typedef struct block{
    uint64_t start;
    uint64_t end;
    uint64_t order;
    uint32_t synteny_level;
} Block;

//A table of length equal to a genome with pointers to blocks per position
typedef struct blocktable{
    Block * blocks;
} Block_table;

typedef struct sequence{
    uint64_t id;    //Label of the sequence
    uint64_t len;   //Length in nucleotides of the sequence
    uint64_t acum;  //Accumulated length from the sequences found before in the file (if any)
    Block_table * block_table; //A table for quick access (to be substituted by hash table) to blocks
} Sequence;

class memory_pool{

private:
    char ** mem_pool;
    uint64_t * base;
    uint64_t current_pool;

public:
    memory_pool();
    void * request_bytes(uint64_t bytes);
    ~memory_pool();
};

#endif
