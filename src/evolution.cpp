#include "evolution.h"

void dna_generator_gc(uint64_t l, char * s, std::uniform_int_distribution<uint64_t> * u, std::default_random_engine * r){
    for(uint64_t i=0;i<l;i++){
        s[i] = nts[(*u)(*r)];
    }
}

dna_mutation::dna_mutation(long double p, a_sequence * sequence, uint64_t * s_len, uint64_t seed){
    this->sequence = sequence;
    this->p = p;
    this->s_len = s_len;
    this->d_r_unif = std::uniform_real_distribution<long double>(0.0, 1.0);
    this->d_u_unif = std::uniform_int_distribution<uint64_t>(0, 3);
    this->generator = std::default_random_engine(seed);
    
}

void dna_mutation::step(){
    uint64_t i;
    char c;
    

    for(i=0;i<*this->s_len;i++){
        if(this->d_r_unif(this->generator) <= this->p){
            // Mutate individual nucleotide
            c = this->sequence->s[i];
            while(c == this->sequence->s[i]){
                this->sequence->s[i] = nts[this->d_u_unif(this->generator)];
            }
        }
    }
}


dna_duplication::dna_duplication(long double p, a_sequence * sequence, uint64_t * s_len, uint64_t d_len, uint64_t seed){
    this->sequence = sequence;
    this->p = p;
    this->s_len = s_len;
    this->d_len = d_len;
    this->d_r_unif = std::uniform_real_distribution<long double>(0.0, 1.0);
    this->d_u_unif = std::uniform_int_distribution<uint64_t>(0, 3);
    this->generator = std::default_random_engine(seed);
}

void dna_duplication::step(){
    uint64_t i;    
    for(i=0;i<(*this->s_len - this->d_len);i++){
        if(this->d_r_unif(this->generator) <= this->p){

            // Generate insertion
            // First add space
            this->sequence->s = (char *) realloc(this->sequence->s, (*this->s_len+this->d_len)*sizeof(char));
            // Displace region
            memmove(&this->sequence->s[i+this->d_len], &this->sequence->s[i], *this->s_len-i);
            // And generate the new region
            memcpy(&this->sequence->s[i], &this->sequence->s[i+this->d_len], this->d_len);
            // And increase the size of the sequence for all 
            *this->s_len += this->d_len;
        }
    }
}

dna_insertion::dna_insertion(long double p, a_sequence * sequence, uint64_t * s_len, uint64_t i_len, uint64_t seed){
    this->sequence = sequence;
    this->p = p;
    this->s_len = s_len;
    this->i_len = i_len;
    this->d_r_unif = std::uniform_real_distribution<long double>(0.0, 1.0);
    this->d_u_unif = std::uniform_int_distribution<uint64_t>(0, 3);
    this->generator = std::default_random_engine(seed);
}

void dna_insertion::step(){
    uint64_t i;
    
    for(i=0;i<*this->s_len;i++){
        if(this->d_r_unif(this->generator) <= this->p){
            // Generate insertion
            // First add space
            this->sequence->s = (char *) realloc(this->sequence->s, (*this->s_len+this->i_len)*sizeof(char));
            // Displace region
            memmove(&this->sequence->s[i+this->i_len], &this->sequence->s[i], *this->s_len-i);
            // And generate the new region
            dna_generator_gc(this->i_len, &this->sequence->s[i], &this->d_u_unif, &this->generator);
            // And increase the size of the sequence for all 
            *this->s_len += this->i_len;
        }
    }
}

dna_deletion::dna_deletion(long double p, a_sequence * sequence, uint64_t * s_len, uint64_t d_len, uint64_t seed){
    this->sequence = sequence;
    this->p = p;
    this->s_len = s_len;
    this->d_len = d_len;
    this->d_r_unif = std::uniform_real_distribution<long double>(0.0, 1.0);
    this->d_u_unif = std::uniform_int_distribution<uint64_t>(0, 3);
    this->generator = std::default_random_engine(seed);
}

void dna_deletion::step(){
    uint64_t i;
    
    for(i=this->d_len;i<*this->s_len;i++){
        if(this->d_r_unif(this->generator) <= this->p){
            // Generate deletion
            
            
            // Displace region
            memmove(&this->sequence->s[i-this->d_len], &this->sequence->s[i], *this->s_len-i);
            // Remove space
            this->sequence->s = (char *) realloc(this->sequence->s, (*this->s_len-this->d_len)*sizeof(char));
            // And decrease the size of the sequence for all 
            *this->s_len -= this->d_len;
        }
    }
}