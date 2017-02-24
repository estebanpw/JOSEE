
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

void dna_generator_gc(uint64_t l, char * s, std::uniform_int_distribution<uint64_t> * u, std::default_random_engine * r){
    for(uint64_t i=0;i<l;i++){
        s[i] = nts[(*u)(*r)];
    }
}

class mutation{
private:
    a_sequence * sequence;
    uint64_t * s_len;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> d_r_unif;
    std::uniform_int_distribution<uint64_t> d_u_unif;
    double p;

public:
    mutation(double p, a_sequence * sequence, uint64_t * s_len);
    void step();

};

class duplication{
private:
    a_sequence * sequence;
    uint64_t * s_len;
    uint64_t d_len;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> d_r_unif;
    std::uniform_int_distribution<uint64_t> d_u_unif;
    double p;

public:
    duplication(double p, a_sequence * sequence, uint64_t * s_len, uint64_t d_len);
    void step();

};

class insertion{
private:
    a_sequence * sequence;
    uint64_t * s_len;
    uint64_t i_len;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> d_r_unif;
    std::uniform_int_distribution<uint64_t> d_u_unif;
    double p;

public:
    insertion(double p, a_sequence * sequence, uint64_t * s_len, uint64_t i_len);
    void step();

};

class deletion{
private:
    a_sequence * sequence;
    uint64_t * s_len;
    uint64_t d_len;
    std::default_random_engine generator;
    std::uniform_real_distribution<double> d_r_unif;
    std::uniform_int_distribution<uint64_t> d_u_unif;
    double p;

public:
    deletion(double p, a_sequence * sequence, uint64_t * s_len, uint64_t d_len);
    void step();

};

