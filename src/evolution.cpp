#include "evolution.h"

mutation::mutation(double p, a_sequence * sequence, uint64_t * s_len){
    this->sequence = sequence;
    this->p = p;
    this->s_len = s_len;
    this->d_r_unif = std::uniform_real_distribution<double>(0.0, 1.0);
    this->d_u_unif = std::uniform_int_distribution<uint64_t>(0, 3);
    
}

void mutation::step(){
    uint64_t i;
    char c;
    
    for(i=0;i<*this->s_len;i++){
        if(this->d_r_unif(this->generator) > this->p){
            // Mutate individual nucleotide
            c = this->sequence->s[i];
            while(c == this->sequence->s[i]){
                this->sequence->s[i] = nts[this->d_u_unif(this->generator)];
            }
        }
    }
}


duplication::duplication(double p, a_sequence * sequence, uint64_t * s_len, uint64_t d_len){
    this->sequence = sequence;
    this->p = p;
    this->s_len = s_len;
    this->d_len = d_len;
    this->d_r_unif = std::uniform_real_distribution<double>(0.0, 1.0);
    this->d_u_unif = std::uniform_int_distribution<uint64_t>(0, 3);
    
}

void duplication::step(){
    uint64_t i;
    
    for(i=0;i<*this->s_len;i++){
        if(this->d_r_unif(this->generator) > this->p){
            // Generate insertion
            // First add space
            this->sequence->s = (char *) realloc(this->sequence->s, (*this->s_len+2*this->d_len)*sizeof(char));
            // Displace region
            memmove(&this->sequence->s[i], &this->sequence->s[i+this->d_len], this->d_len);
            // And generate the new region
            dna_generator_gc(this->d_len, &this->sequence->s[i], &this->d_u_unif, &this->generator);
            // And increase the size of the sequence for all 
            *this->s_len += this->d_len;
        }
    }
}

insertion::insertion(double p, a_sequence * sequence, uint64_t * s_len, uint64_t i_len){
    this->sequence = sequence;
    this->p = p;
    this->s_len = s_len;
    this->i_len = i_len;
    this->d_r_unif = std::uniform_real_distribution<double>(0.0, 1.0);
    this->d_u_unif = std::uniform_int_distribution<uint64_t>(0, 3);
    
}

void insertion::step(){
    uint64_t i;
    
    for(i=0;i<*this->s_len;i++){
        if(this->d_r_unif(this->generator) > this->p){
            // Generate insertion
            // First add space
            this->sequence->s = (char *) realloc(this->sequence->s, (*this->s_len+this->i_len)*sizeof(char));
            // Displace region
            memmove(&this->sequence->s[i], &this->sequence->s[i+this->i_len], this->i_len);
            // And generate the new region
            dna_generator_gc(this->i_len, &this->sequence->s[i], &this->d_u_unif, &this->generator);
            // And increase the size of the sequence for all 
            *this->s_len += this->i_len;
        }
    }
}

deletion::deletion(double p, a_sequence * sequence, uint64_t * s_len, uint64_t d_len){
    this->sequence = sequence;
    this->p = p;
    this->s_len = s_len;
    this->d_len = d_len;
    this->d_r_unif = std::uniform_real_distribution<double>(0.0, 1.0);
    this->d_u_unif = std::uniform_int_distribution<uint64_t>(0, 3);
    
}

void deletion::step(){
    uint64_t i;
    
    for(i=this->d_len;i<*this->s_len;i++){
        if(this->d_r_unif(this->generator) > this->p){
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