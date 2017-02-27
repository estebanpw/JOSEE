
#include <inttypes.h>
#include <cstring>
#include <iostream>
#include <list>
#include <vector>
#include <random>
#include "commonFunctions.h"

#pragma pack(push, 1)

static char nts[] = {'A', 'C', 'G', 'T'};

struct a_sequence{
    char * s;
};

class dna_mutation;
class dna_duplication;
class dna_insertion;
class dna_deletion;
class dna_inversion;

void dna_generator_gc(uint64_t l, char * s, std::uniform_int_distribution<uint64_t> * u, std::default_random_engine * r);



class dna_mutation{
private:
    a_sequence * sequence;
    uint64_t * s_len;
    std::default_random_engine generator;
    std::uniform_real_distribution<long double> d_r_unif;
    std::uniform_int_distribution<uint64_t> d_u_unif;
    long double p;

public:
    dna_mutation(long double p, a_sequence * sequence, uint64_t * s_len, uint64_t seed);
    void set_p(long double p){ this->p = p; }
    long double get_p(){ return this->p; }
    void step();

};

class dna_duplication{
private:
    a_sequence * sequence;
    uint64_t * s_len;
    uint64_t d_len;
    std::default_random_engine generator;
    std::uniform_real_distribution<long double> d_r_unif;
    std::uniform_int_distribution<uint64_t> d_u_unif;
    long double p;

public:
    dna_duplication(long double p, a_sequence * sequence, uint64_t * s_len, uint64_t d_len, uint64_t seed);
    void set_p(long double p){ this->p = p; }
    long double get_p(){ return this->p; }
    void step();

};

class dna_insertion{
private:
    a_sequence * sequence;
    uint64_t * s_len;
    uint64_t i_len;
    std::default_random_engine generator;
    std::uniform_real_distribution<long double> d_r_unif;
    std::uniform_int_distribution<uint64_t> d_u_unif;
    long double p;

public:
    dna_insertion(long double p, a_sequence * sequence, uint64_t * s_len, uint64_t i_len, uint64_t seed);
    void set_p(long double p){ this->p = p; }
    long double get_p(){ return this->p; }
    void step();

};

class dna_deletion{
private:
    a_sequence * sequence;
    uint64_t * s_len;
    uint64_t d_len;
    std::default_random_engine generator;
    std::uniform_real_distribution<long double> d_r_unif;
    std::uniform_int_distribution<uint64_t> d_u_unif;
    long double p;

public:
    dna_deletion(long double p, a_sequence * sequence, uint64_t * s_len, uint64_t d_len, uint64_t seed);
    void set_p(long double p){ this->p = p; }
    long double get_p(){ return this->p; }
    void step();

};

class dna_inversion{
private:
    a_sequence * sequence;
    uint64_t * s_len;
    uint64_t i_len;
    std::default_random_engine generator;
    std::uniform_real_distribution<long double> d_r_unif;
    std::uniform_int_distribution<uint64_t> d_u_unif;
    long double p;

public:
    dna_inversion(long double p, a_sequence * sequence, uint64_t * s_len, uint64_t i_len, uint64_t seed);
    void set_p(long double p){ this->p = p; }
    long double get_p(){ return this->p; }
    void step();

};

class dna_transposition{
private:
    a_sequence * sequence;
    uint64_t * s_len;
    uint64_t t_len;
    std::default_random_engine generator;
    std::uniform_real_distribution<long double> d_r_unif;
    std::uniform_int_distribution<uint64_t> d_u_unif;
    long double p;

public:
    dna_transposition(long double p, a_sequence * sequence, uint64_t * s_len, uint64_t t_len, uint64_t seed);
    void set_p(long double p){ this->p = p; }
    long double get_p(){ return this->p; }
    void step();

};