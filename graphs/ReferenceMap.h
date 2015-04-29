#ifndef REFERENCEMAP_H_
#define REFERENCEMAP_H_

#include <iostream>
#include <string>
#include <vector>
#include <algorithm>
#include <ctime>
#include <time.h>
#include <unordered_map>

using namespace std;

class ReferenceMap {

public: 
    ReferenceMap();
    
    bool addReference(string ref_id, string seq);
    bool addKmerToReference(string ref_id, string kmer, int ind);
    
    int getIndex(string ref_id, string kmer);
    string getSeq(string ref_id);
    vector<string> getIDs();
    
    int size();
    
private:
    unordered_map<string, unordered_map<string,int> > h_;
    unordered_map<string, string> ref_h_;
};

#endif
