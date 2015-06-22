#ifndef D_CLASSIFY_HPP
#define D_CLASSIFY_HPP

#include <iostream>
#include <fstream>
#include <string>
#include <limits>
#include <algorithm>    // std::reverse
#include <cmath>
#include <vector>
#include <functional>

#include "prettyprint.hpp"

#include "file_io/FastaParser.hpp"
#include "file_io/FastaRefID.hpp"

#include "d_align/DLabel.hpp"

#include "tests/TestDClassify.hpp"

#include <seqan/align.h>

typedef seqan::String<char> TSequence; // sequence type
typedef seqan::StringSet<TSequence> TStringSet; // container for strings
typedef seqan::StringSet<TSequence, seqan::Dependent<> > TDepStringSet; // dependent string set
typedef seqan::Graph<seqan::Alignment<TDepStringSet> > TAlignGraph; // alignment graph
typedef seqan::Align<TSequence,seqan::ArrayGaps> TAlign;
typedef seqan::Row<TAlign>::Type TRow;

#ifdef DEBUG
#define DCLASS_DEBUG_PRINT(EXPR) DEBUG_PRINT("[D_CLASSIFY]", EXPR)
#else
#define DCLASS_DEBUG_PRINT(x) do {} while (0)
#endif

class DClassify {
public:
    DClassify();
    DClassify(std::string dref_fasta);
    virtual ~DClassify();

    std::vector<DLabel> classify_d(std::string seq,
				   std::pair<int,int> v_part,
				   std::pair<int,int> j_part);
    
    void set_num_d_report(int num_d);
    
    friend class TestDClassify; // for unit testing
    
protected:
    
    void read_in_refs(std::string ref_fasta, 
		      std::map<std::string,std::string> &ref_h);

    std::vector<DLabel> score_d(std::string seq,
				std::pair<int,int> v_part,
				std::pair<int,int> j_part);
    
    
    std::map<std::string,std::string> drefs_h_;
    
    int num_report_;
    int match_s_; 
    int mismatch_s_;
    int indel_s_;
    //int get_local_score(std::string ref_seq, std::string read_substr);
    pair<int,int> get_local_score(std::string ref_seq, std::string read_substr);
};

#endif
