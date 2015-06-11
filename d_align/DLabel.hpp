#ifndef D_LABEL_HPP
#define D_LABEL_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <limits>
#include <algorithm>    // std::reverse
#include <cmath>
#include <vector>

using namespace std;

class DLabel {
public:
    DLabel();
    DLabel(std::string s, int sc, int iden);
    virtual ~DLabel();
    
    bool operator > (const DLabel& dl) const;
    bool operator < (const DLabel& dl) const;
    
    friend ostream& operator<<(ostream& out, const DLabel& dl);
    
    std::string label;
    int score;
    int ident;
};


#endif
