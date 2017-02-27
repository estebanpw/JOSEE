
#include <inttypes.h>
#include <cstring>
#include <iostream>
#include <list>
#include <vector>
#include <random>

#pragma pack(push, 1)

static char nts[] = {'A', 'C', 'G', 'T'};

struct a_sequence{
    char * s;
};

void dna_generator_gc(uint64_t l, char * s, std::uniform_int_distribution<uint64_t> * u, std::default_random_engine * r);

class mutation{
private:
    a_sequence * sequence;
    uint64_t * s_len;
    std::default_random_engine generator;
    std::uniform_real_distribution<long double> d_r_unif;
    std::uniform_int_distribution<uint64_t> d_u_unif;
    long double p;

public:
    mutation(long double p, a_sequence * sequence, uint64_t * s_len, uint64_t seed);
    void step();

};

class duplication{
private:
    a_sequence * sequence;
    uint64_t * s_len;
    uint64_t d_len;
    std::default_random_engine generator;
    std::uniform_real_distribution<long double> d_r_unif;
    std::uniform_int_distribution<uint64_t> d_u_unif;
    long double p;

public:
    duplication(long double p, a_sequence * sequence, uint64_t * s_len, uint64_t d_len, uint64_t seed);
    void step();

};

class insertion{
private:
    a_sequence * sequence;
    uint64_t * s_len;
    uint64_t i_len;
    std::default_random_engine generator;
    std::uniform_real_distribution<long double> d_r_unif;
    std::uniform_int_distribution<uint64_t> d_u_unif;
    long double p;

public:
    insertion(long double p, a_sequence * sequence, uint64_t * s_len, uint64_t i_len, uint64_t seed);
    void step();

};

class deletion{
private:
    a_sequence * sequence;
    uint64_t * s_len;
    uint64_t d_len;
    std::default_random_engine generator;
    std::uniform_real_distribution<long double> d_r_unif;
    std::uniform_int_distribution<uint64_t> d_u_unif;
    long double p;

public:
    deletion(long double p, a_sequence * sequence, uint64_t * s_len, uint64_t d_len, uint64_t seed);
    void step();

};

