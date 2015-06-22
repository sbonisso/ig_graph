/*
 * Encoding.h
 *
 *  Created on: 
 *      Author: sbonisso
 */

#ifndef ENCODING_H_
#define ENCODING_H_

#include <string>
#include <unordered_map>
#include <iostream>
#include <fstream>
#include <sstream>
#include <map>
#include <unordered_map>

class Encoding {

public:
    Encoding() {}
    virtual ~Encoding() {}
    
    static unsigned int encode_kmer(std::string kmer);
    static std::string decode_kmer(unsigned int kmer_val, int k);
    static  unsigned int update_kmer(unsigned int kmer_val, char bp, int k);
};


#endif
