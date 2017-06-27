#ifndef STRUCTS_H
#define STRUCTS_H

#include <inttypes.h>
#include <iostream>
#include <list>
#include <vector>

#pragma pack(push, 1)

#define SEQ_REALLOC 5000000
//#define MAX_READ_SIZE 5000 //Maximum length of read to have a portion of the table allocated
#define MAX_WINDOW_SIZE 70 //Maximum window length to explore NW table
#define SYN_TABLE_REALLOC 10000
#define INIT_SEQS 20
#define INIT_ANNOTS 5000
#define INIT_BLOCKS_TO_ADD 50
#define INIT_REARRANGEMENT_CAPACITY 100
#define WRITE_BUFFER_CAPACITY 100000 //100 KB
#define READLINE 2000
#define READBUF 50000000 //50MB
#define INIT_TRIM_FRAGS 10000
#define TABLE_RATE 100 //hash table lengh divisor
#define INIT_CANDIDATES_ALIGN 100 //For wordbucket list
#define PRINT_RATE 70
#define MAX_LINE 2048

#define MAX_MEM_POOLS 256
#define POOL_SIZE 1024*1024*128 //128 MB

#define NOFRAG 0
#define OPENFRAG 1
#define COVERFRAG 2
#define CLOSEFRAG 3

#define FORWARD 1
#define REVERSE 2
#define MIXED 3

#define INSERTION 1
#define DELETION 2
#define NOTHING 3

#define POINT 4

#define UINT64_T_MAX 0xFFFFFFFFFFFFFFFF
#define SEQUENCE_INDELS_LEN 1000*1000*2 // 1 M

#define IGAP -24
#define EGAP -8

#define EVENTS_WRITE_TIME 5

#define PRINT_BLOCKS 1
#define PRINT_BLOCKS_AND_FRAGS 2

#define MIN(X, Y) (((X) < (Y)) ? (X) : (Y))
#define MAX(X, Y) (((X) > (Y)) ? (X) : (Y))

// To globally measure RAM usage
extern uint64_t total_bytes_in_use;

//Class prototypes
class memory_pool;
class hash_table;
class sequence_manager;

//Struct prototypes
struct synteny_list;

//Enum of events
enum Event { inversion, duplication, transposition, insertion, deletion, indel, none, concatenation}; // Indel is generic for insertion or deletion, and has to be decided

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

//Struct for NW two rows
struct cell{
	int64_t score;
	uint64_t xe;
	uint64_t ye;
	uint64_t xs;
	uint64_t ys;
	uint64_t igaps;
	uint64_t egaps;
	uint64_t ident;
};

//Struct for full NW
struct cell_f{
    int64_t score;
    uint64_t xfrom;
    uint64_t yfrom;
};
struct positioned_cell{
    int64_t score;
    uint64_t xpos;
    uint64_t ypos;
};
struct best_cell{
    struct positioned_cell c;
    uint64_t j_prime;
};
typedef struct{
    uint64_t identities;
    uint64_t length;
    uint64_t igaps;
    uint64_t egaps;
} BasicAlignment;

//Sequence descriptor
typedef struct sequence{
    uint64_t id;    //Label of the sequence
    uint64_t len;   //Length in nucleotides of the sequence
    uint64_t acum;  //Accumulated length from the sequences found before in the file (if any)
    uint32_t coverage; //The percentage of bases covered by fragments of a minimum length
    uint64_t n_frags; //Number of fragments that the sequence had
    char * seq; //DNA sequence
} Sequence;

// Struct used to compute the alignment between two hits
typedef struct quickfrag{
    uint64_t x_start;
    uint64_t y_start;
    uint64_t t_len;
    int64_t diag;
    long double sim;
    Sequence * x;
    Sequence * y;
} Quickfrag;

// A collection of fragments 
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
    struct synteny_list * present_in_synteny;   //To tell whether it has already been used in a synteny block
    unsigned char strand_in_synteny;    //The strand that it has at the synteny block
    struct block * prev;
    struct block * next;
    uint64_t id;
} Block;

//A struct to hold arrays of blocks
typedef struct holder{
    Block ** cons_AB;
    Block ** cons_BC;
    uint64_t n_sequences;
} Holder;

//Word struct that identifies a kmer in a sequence
typedef struct word{
    uint64_t hash;
    uint64_t pos;
    char strand;
    Block * b;
} Word;

typedef struct point{
    uint64_t x;
    uint64_t y;
} Point;

//A synteny block is a collection of blocks
typedef struct synteny_block{
    Block * b;
    struct synteny_block * next;
} Synteny_block;

//A synteny list is a collection of synteny blocks
typedef struct synteny_list{
    Synteny_block * sb;
    uint64_t synteny_level;
    struct synteny_list * next;
    struct synteny_list * prev;
    uint64_t id;
} Synteny_list;

struct triplet{
    Synteny_list * A;
    Synteny_list * B;
    Synteny_list * C;
    struct triplet * next;
    Event etype;
};
//Class for allocating memory only once and requesting particular amounts of bytes
class memory_pool{

private:
    std::vector<char *> * mem_pool;
    std::vector<uint64_t> * base;
    uint64_t current_pool;
    uint64_t max_pools;
    uint64_t pool_size;

public:
    memory_pool(uint64_t pool_size);
    void * request_bytes(uint64_t bytes);
    void reset_n_bytes(uint64_t bytes);
    void reset_to(uint64_t pool, uint64_t position){ this->current_pool = 0; this->base->at(this->current_pool) = 0;}
    void full_reset();
    uint64_t get_pools_used(){ return current_pool; }
    ~memory_pool();
};

//Bucket structs to hold blocks in hash table (1) and words in dictionary (2)
typedef struct bucket {
	Block b;
	struct bucket * next;
} Bucket;

typedef struct wordbucket{
    Word w;
    struct wordbucket * next;
} Wordbucket;

// For pthreads alignment of NW
struct alignment_arguments{
    char * seq_A;
    Sequence * s_x;
    Sequence * s_y;
    uint64_t start_A;
    uint64_t end_A;
    char * seq_B;
    uint64_t start_B;
    uint64_t end_B;
    int64_t iGap;
    int64_t eGap;
    struct cell * mc;
    struct cell * f0;
    struct cell * f1;
    struct cell * alignment;
    struct cell * alignment_reverse;
    char * seq_for_reverse;
};

typedef struct annotation{
    uint64_t start;
    uint64_t end;
    char strand;
    char * product; //Requires hard copy
} Annotation;

// Sequence class management
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
    void print_sequence_region(FILE * fout, uint64_t label, uint64_t from, uint64_t to);
    ~sequence_manager();
};

//Hash-table class for Blocks management
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
    Block * get_previous_block(Block * b);
    Block * get_next_block(Block * b);
    void release_temp_last_blocks();
    void remove_block(Block * b);
    void write_blocks_and_breakpoints_to_file(FILE * out_blocks, FILE * out_breakpoints);
    Block ** last_blocks; // Used to generate the double link between blocks

private:
	uint64_t compute_hash(uint64_t key);
    void insert_x_side(struct FragFile * f);
    void insert_y_side(struct FragFile * f);
};

// Dictionary class to handle words and hits between sequences
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

class markdown_event_hash{
private:
    triplet ** array;
    uint64_t size;
    memory_pool * mp;
public:
    markdown_event_hash(uint64_t size);
    void put(triplet * k);
    triplet * find_triplet(triplet * k);
    void remove(triplet * k);
    ~markdown_event_hash();

private:
    uint64_t compute_hash(triplet * k);
};

typedef struct slist{
    Sequence * s;
    struct slist * next;
} Slist;



// There will be one strand matrix per synteny block to compute if a synteny list might have a reversion
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
    int do_forwards_require_less_changes(uint64_t genome);
    int get_strands_type();
    unsigned char get_strands(uint64_t l1, uint64_t l2);
    void set_strands(uint64_t l1, uint64_t l2, unsigned char v);
    uint64_t get_frags_forward(){ return this->acu_frags_forward; }
    uint64_t get_frags_reverse(){ return this->acu_frags_reverse; }
    void reset();
    void add_fragment_strands(Synteny_list * sbl);
    void print_strand_matrix();
    ~strand_matrix();
};

// Structs to pass events and handle them in one function
struct e_inversion{
    Block * inv; //A block with inversion
};

struct e_duplication{
    Block * orig; //The original block from which the duplication "duplicated"
    Block * dup; //A block that has been detected as duplication
};

struct e_transposition{
    Block * before_trans; 
    Block * transposed;    
};

struct e_insertion{
    Block * insertion;
};

struct e_deletion{
    Block * deletion;
    uint64_t removed;
};

struct e_concatenation{
    Block * involved;
    uint64_t coords1;
    uint64_t coords2;
};

// Struct to modify blocks given previous rearrangments
struct rearrangement{
    int64_t mod_coordinates; //coordinate offset to add
    int64_t mod_order; //order offset to add 
    uint64_t b1_id; //Range start of blocks affected
    uint64_t b2_id; //Range end of blocks affected
    uint64_t ends_at; //Ending id
    uint64_t type;
    uint64_t affects_who; //-1 for all, either specify genome label
};

// Queue-like class to handle all operations that have to be done 
class events_queue{
private:
    std::list<rearrangement> * rea_queue;
    std::list<rearrangement>::iterator rea_itera;
    rearrangement * aggregated_r;
    uint64_t n_sequences;
public:
    events_queue(uint64_t n_sequences);
    void insert_event(rearrangement r);
    uint64_t get_queue_size(){ return this->rea_queue->size(); }
    rearrangement * get_aggregated_event(Block * b, uint64_t s_id);
    void begin_iterator(){ this->rea_itera = this->rea_queue->begin(); }
    void print_queue();
    ~events_queue();
};

// Class to register events and write them to file
class ee_log{
private:
    FILE * logfile;
    char * write_buffer;
    uint64_t write_index;
    uint64_t event_count;
    char tmp[256];
    char path[MAX_LINE];

public:
    ee_log(FILE * logfile, char * path);
    void register_event(Event e, void * event_data); 
    void force_write();   
    ~ee_log();

private:
    void write(const char * data);

};

struct Event_handling{
    Event type_of_event;
    bool * genomes_affected;
};
struct Indel_handling{
    Synteny_block * concat;
    uint64_t * bp_lengths;
    uint64_t n_sequences;
};

// Struct for reducing parameters in multiple alignment function 
struct arguments_multiple_alignment{
    sequence_manager * seq_man;
    char ** seq_for_reverse;
    Synteny_list * sbl;
    Quickfrag ** qfmat;
    unsigned char ** qfmat_state;
    int iGap;
    int eGap;
    struct cell ** mc;
    struct cell ** f0;
    struct cell ** f1;
    pthread_t * threads;
    double ** submat;
    memory_pool * mp;
    // For full NW
    char * aux_dummy_sequence; // The one that is not used in the backtrackings
    char ** recon_X;
    uint64_t * sequence_ids;  // Tells which one of the seqs in seq_X is after aligning
    char ** recon_Y;
    char ** recon_Z; 
    char ** seq_X;
    char ** seq_Y;
    char ** seq_Z; 
    int64_t * cell_path_y;
    struct positioned_cell * mc_f;
    struct cell_f ** table_f;
    char * writing_buffer_alignment;
    long double window;
    uint64_t n_sequences;

};


#endif
