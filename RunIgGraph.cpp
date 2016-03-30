#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>

#include <math.h>
#include <thread>
#include <mutex>

#include <tclap/CmdLine.h> 	  // for parsing command line args
#include "ProcessMemUsage.h" 	  // for checking process memory usage
#include "Utilities.hpp"          // for debug printing macro

#include "file_io/FastaParser.hpp"  // for parsing FASTA files
#include "file_io/FastaRefID.hpp"

#include "graphs/CanonicalAntibodyGraph.hpp"
#include "graphs/MultiKCanonAbGraph.hpp"
#include "graphs/ColorProfileMatrix.hpp"
#include "graphs/CreateProfile.hpp"

#include "d_align/DClassify.hpp"

using namespace std;
using namespace TCLAP;

#ifdef DEBUG
#define MAIN_DEBUG_PRINT(EXPR) DEBUG_PRINT("[MAIN]", EXPR)
#else
#define MAIN_DEBUG_PRINT(x) do {} while (0)
#endif

std::mutex mu_cerr;
std::mutex mu_out;
/**
 * worker thread function
 */
void process_fasta(ostream &out, CreateProfile cp, string &readFasta, 
		   int start_index, int end_index, int thread_index) {
    FastaRefID<FastaParser> fp(readFasta);
    if(!fp.openFile()) {
	cerr<<readFasta<<" not openend properly!\n";
	exit(1);
    }
    int index = -1;
    int tenth_val = end_index < 10 ? 1 : (end_index-start_index)/10;
    int percent_done = 0;
    while(fp.hasNextSequence()) {
    	index++;
    	pair<string,string> entry = fp.getNextEntry();
    	string readID = fp.getCurrIDLine().substr(1);  // remove the '>'
    	string seq = entry.second;
    	if(index < start_index || index > end_index) { continue; }
	// output each 10%	
	if( (index-start_index) % tenth_val == 0) {
	    mu_cerr.lock();
	    cerr<<"[T"<<thread_index<<"] "<<percent_done<<"%\t"<<index<<endl; 
	    mu_cerr.unlock();
	    percent_done += 10;
	}
    	MAIN_DEBUG_PRINT("READID:\t"<<readID);	
    	MAIN_DEBUG_PRINT("SEQ:\t"<<seq);
    	// one-stop shop for computing everything for sequence
    	cp.compute(seq);
    	// output to stream - mutex lock it
    	mu_out.lock();
    	out<<readID<<"\t"<<cp;
    	mu_out.unlock();
    }
    fp.closeFile();
}
/**
 *
 */
int main(int argc, char **argv) {
    string readFasta, vRefFasta, dRefFasta, jRefFasta, outFile, paramDir;
    int k = 21, maxKeep = 3, scoringScheme = 0;
    int v_k = -1; int d_k = -1; int j_k = -1;
    int num_thread = 1;
    bool no_cdr3 = false;
    bool out_scores = false;
    bool fill_in_d = false;
    string homeDir(getenv("HOME"));
    string default_paramDir("data/model.bin");
    try {
	CmdLine cmd("Command description message", ' ', "0.0");

	ValueArg<std::string> 
	    readNameArg("r","read_file","Name to print1",true,"homer","string");
	cmd.add( readNameArg );

	ValueArg<std::string> 
	    vRefNameArg("v","V_ref_file","Name to print2",true,"blah","string");
	cmd.add( vRefNameArg );

	ValueArg<std::string> 
	    dRefNameArg("d","D_ref_file","Name to print2",true,"blah","string");
	cmd.add( dRefNameArg );

	ValueArg<std::string> 
	    jRefNameArg("j","J_ref_file","Name to print2",true,"blah","string");
	cmd.add( jRefNameArg );

	ValueArg<int> kSizeArg("k","k_param","size of k",false,-1,"int");
	cmd.add( kSizeArg );
	
	ValueArg<int> vkSizeArg("V","v_k_param","size of k for V",false,-1,"int");
	cmd.add( vkSizeArg );
	ValueArg<int> dkSizeArg("D","d_k_param","size of k for D",false,-1,"int");
	cmd.add( dkSizeArg );
	ValueArg<int> jkSizeArg("J","j_k_param","size of k for J",false,-1,"int");
	cmd.add( jkSizeArg );

	ValueArg<int> 
	    maxReportArg("m", "max_report", 
			 "maximum gene segements to report, [1,2,..], default = 2", 
			 false, 2, "int");
	cmd.add(maxReportArg);

	ValueArg<int> 
	    scoringArg("s", "scoring", "[0 = standard, 1 = prob]", false, 1, "int");
	cmd.add(scoringArg);

	SwitchArg cdr3Arg("c", "no_cdr3", "Omit computing CDR3 sequence", false);
	cmd.add(cdr3Arg);

	SwitchArg 
	    outScoresArg("S", "output_scores", 
			 "flag to output top scorig of each V, D, and J", false);
	cmd.add(outScoresArg);

	ValueArg<int> nThreadArg("t","num_thread","number of threads",false,1,"int");
	cmd.add( nThreadArg );

	ValueArg<std::string> 
	    paramDirArg("p", "param_dir", "path to parameter directory", 
			false, default_paramDir, "string");
	cmd.add(paramDirArg);

	ValueArg<std::string> 
	    outFileNameArg("o","output_file","Name of output file",
			   false,"","string");
	cmd.add( outFileNameArg );

	SwitchArg 
	    fillInDArg("f", "fill_in_d", 
		       "flag to perform additional alignment for filling in missing D gene-segments", false);
	cmd.add(fillInDArg);
	
	// Parse the argv array.
	cmd.parse( argc, argv );
	// Get the value parsed by each arg.
	readFasta = readNameArg.getValue();
	vRefFasta = vRefNameArg.getValue();
	dRefFasta = dRefNameArg.getValue();
	jRefFasta = jRefNameArg.getValue();
	outFile = outFileNameArg.getValue();
	// parameter arg values
	k = kSizeArg.getValue();
	v_k = vkSizeArg.getValue();
	d_k = dkSizeArg.getValue();
	j_k = jkSizeArg.getValue();
	maxKeep = maxReportArg.getValue();
	scoringScheme = scoringArg.getValue();
	no_cdr3 = cdr3Arg.getValue();
	out_scores = outScoresArg.getValue();
	// experimental
	fill_in_d = fillInDArg.getValue();
	// parameter dir
	paramDir = paramDirArg.getValue();
	// num threads
	num_thread = nThreadArg.getValue();

	cerr<<"READS =\t"<<readFasta<<endl;
	cerr<<"REFS =\t"<<vRefFasta<<endl;
	cerr<<"OUTFILE =\t"<<outFile<<endl;
	cerr<<"K =\t"<<k<<endl;
	cerr<<"K_v = "<<v_k<<"\tK_d = "<<d_k<<"\tJ_k = "<<j_k<<endl;
	cerr<<"MAX REPORT =\t"<<maxKeep<<endl;
	cerr<<"NO CDR3 =\t"<<no_cdr3<<endl;
	cerr<<"NUM THREAD =\t"<<num_thread<<endl;
    }catch(ArgException &e)  // catch any exceptions
    { cerr << "error: " << e.error() << " for arg " << e.argId() << endl; }

    if(k < 0 && v_k < 0) { 
	cerr<<"Must specify a value for k or v_k/d_k/j_k"<<endl; 
	exit(1);
    }
    bool standardScoring = scoringScheme == 0;
    cerr<<"SCORING = "<<(standardScoring ? "STANDARD" : "PROBABILISTIC")<<endl;
    if(!standardScoring) {
	MAIN_DEBUG_PRINT(paramDir);
	(MutationNBModel::getInstance()).setParamDir(paramDir);
#ifdef DEBUG
	int lmerLen = (MutationNBModel::getInstance()).getLMerLen();
#endif
	MAIN_DEBUG_PRINT("L-MER LEN:\t"<<lmerLen);
    }
    //
    ofstream out_s;
    if(outFile != "") { 
	out_s.open(outFile, std::ios::out);
	if(!out_s.is_open()) {
	    cerr<<outFile<<" not opened... exiting\n";
	    exit(1);
	}
    }
    ostream &out = (outFile == "") ? std::cout : out_s;
    
    //CanonicalAntibodyGraph cab(k);
    MultiKCanonAbGraph cab( (k < 0 ? v_k : k), 
			    (k < 0 ? d_k : k),
			    (k < 0 ? j_k : k) );
    
    cab.addVReferences(vRefFasta);
    cab.addDReferences(dRefFasta);
    cab.addJReferences(jRefFasta);   
    //
    CreateProfile cp(&cab, !no_cdr3, out_scores);
    cp.set_scoring(scoringScheme);
    DClassify dc(dRefFasta);
    if(fill_in_d) {
	cp.fill_in_d();
	cp.set_d_classify(&dc);
    }
    MAIN_DEBUG_PRINT("SIZE:\t"<<cp.getProfileSize());
    //
    std::ifstream inFile(readFasta); 
    int num_lines = std::count(std::istreambuf_iterator<char>(inFile), 
			       std::istreambuf_iterator<char>(), '\n');
    std::cerr<<"NUM FASTA LINES:\t"<<num_lines<<std::endl;
    int chunk_size = ceil( (double)((double)num_lines/2.0)/((double)num_thread) );
    
    std::thread *tt = new std::thread[num_thread];
    // create the threads for work ...
    for(int i = 0; i < num_thread; i++) {
    	int strt_i = i*chunk_size;
    	int end_i = (i+1)*chunk_size - 1;
    	tt[i] = std::thread(process_fasta, 
			    std::ref(out), 
			    cp,
			    std::ref(readFasta), 
			    strt_i, end_i, i);
    }
    // ... now join the threads
    for(int i = 0; i < num_thread; i++) { 
    	tt[i].join();
    }
    delete [] tt;   
    
    if(out_s.is_open()) { out_s.close(); }
    
    // checks process memory usage
    double vm, rss;
    process_mem_usage(vm, rss);
    cerr<<"VM: " << vm <<"\tRSS: " << rss<<endl;
    
    return 0;
}

