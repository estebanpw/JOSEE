#include <cstdio>
#include <inttypes.h>
#include <cstring>
#include <cstdlib>
#include <float.h>
#include <math.h>
#include <cstdint>
#include "structs.h"
/*
    Returns the score for a nucl
*/

int valOfNucl(char c);

/*
    Compares two nucleotides and produces a score
*/
int64_t compare_letters(char a, char b);

/*
    Gets maximum length of a group of sequences (uses strlen)
*/
uint64_t get_max_length_of_sequences(char ** a, uint64_t n_x);


/*
    Gets index of maximum length of a group of sequences (uses strlen)
*/
uint64_t get_id_of_longest_sequence(char ** a, uint64_t n_x);

/*
    Extracts a row/column of superposed nucleotides
*/
void pile_chars_up(char ** list_seq_X, char ** list_seq_Y, char * list_piled_X, char * list_piled_Y, uint64_t nx, uint64_t ny, uint64_t posx, uint64_t posy);

/*
    Compares rows/columns of superposed nucleotides to yield a general score
*/
int64_t compare_piled_up_chars(char * list_piled_X, char * list_piled_Y, uint64_t nx, uint64_t ny);

/*
    Builds an anchored alignment for a given group of sequences 
*/
void build_unanchored_alignment(int64_t * cell_path_y, char ** seq_piles_up_X, char ** seq_piles_up_Y, uint64_t n_x, uint64_t n_y);

/*
    Computes the multiple alignment of n to m sequences 
*/
void build_multiple_alignment(char ** reconstruct_X, char ** reconstruct_Y, char ** my_x, char ** my_y, uint64_t nx, uint64_t ny, struct cell_f ** table, struct positioned_cell * mc, char * writing_buffer_alignment, uint64_t xlen, uint64_t ylen, int64_t * cell_path_y, long double * window, int64_t iGap, int64_t eGap, BasicAlignment * ba, char * aux);

/*
    Calculates a path of positions to travel on a NW matrix
*/
void calculate_y_cell_path(Point p0, Point p1, Point p2, Point p3, int64_t * y_points);

/*
    Gets the multiple alignment of a NW table with n and m sequences
*/
void backtrackingNW(char * X, uint64_t Xstart, uint64_t Xend, char * Y, uint64_t Ystart, uint64_t Yend, struct cell_f ** table, char * rec_X, char * rec_Y, struct best_cell * bc, uint64_t * ret_head_x, uint64_t * ret_head_y, int64_t * cell_path_y, uint64_t window_size, BasicAlignment * ba);
/*
    Calculates NW table with two rows and stores a cellpath of scores, identities, gaps and starting and ending positions
*/

struct cell NWscore2rows(char * X, uint64_t Xstart, uint64_t Xend, char * Y, uint64_t Ystart, uint64_t Yend, int64_t iGap, int64_t eGap, struct cell * mc, struct cell * f0, struct cell * f1);