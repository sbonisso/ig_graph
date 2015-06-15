#include "Encoding.hpp"
/**
 * given a k-mer string, return the encoding as an unsigned int
 * can have max k-mer length of 16
 */
unsigned int Encoding::encode_kmer(std::string kmer) {
    int len = (int)kmer.size()-1;
    unsigned val = 0;
    for(int i = len; i >= 0; i--) {
	if(kmer[i] == 'A') { val += 0; }
	else if(kmer[i] == 'C') { val += 1; }
	else if(kmer[i] == 'G') { val += 2; }
	else if(kmer[i] == 'T') { val += 3; }
	else {}
	if(i > 0) { val = val << 2; }
    }
    return val;
}
/**
 * given a value for a k-mer, and a length, return the encoded k-mer as a string
 */
std::string Encoding::decode_kmer(unsigned int kmer_val, int k) {
    std::string kmer = "";
    for(int i = 0; i < k; i++) {
	unsigned int prev_val = kmer_val;
	kmer_val = kmer_val >> 2;
	kmer_val = kmer_val << 2;
	unsigned int diff = (prev_val - kmer_val);
	if(diff == 0) { kmer += "A"; }
	else if(diff == 1) { kmer += "C"; }
	else if(diff == 2) { kmer += "G"; }
	else if(diff == 3) { kmer += "T"; }
	else {}
	kmer_val = kmer_val >> 2;
    }
    return kmer;
}
/**
 * given an encoded k-mer, next bp, and length, return new k-mer value of next
 * shifted, k-mer. For use in creating k-mers from longer reads efficiently
 */
unsigned int Encoding::update_kmer(unsigned int kmer_val, char bp, int k) {
    // remove last char
    kmer_val = kmer_val >> 2;  
    // add new char to most sig end
    unsigned int tmp_val = 0;
    if(bp == 'A') { tmp_val += 0; }
    else if(bp == 'C') { tmp_val += 1; }
    else if(bp == 'G') { tmp_val += 2; }
    else if(bp == 'T') { tmp_val += 3; }
    else {}
    tmp_val = tmp_val << (2*(k-1));
    return (kmer_val | tmp_val);
}
