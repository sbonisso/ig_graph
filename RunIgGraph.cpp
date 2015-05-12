
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <utility>
#include "Utils.h"

#include <tclap/CmdLine.h> 	  // for parsing command line args
#include "ProcessMemUsage.h" 	  // for checking process memory usage

#include "file_io/FastaParser.h"  // for parsing FASTA files
#include "file_io/FastaMult.h"
#include "file_io/FastaRefID.h"
#include "file_io/FastaWriter.h"  // for writing FASTA files

#include "graphs/CanonicalAntibodyGraph.h"
#include "graphs/MultiKCanonAbGraph.h"
#include "graphs/ColorProfileMatrix.h"
#include "graphs/CreateProfile.h"

using namespace std;
using namespace TCLAP;

#ifdef DEBUG
#define MAIN_DEBUG_PRINT(EXPR) DEBUG_PRINT("[MAIN]", EXPR)
#else
#define MAIN_DEBUG_PRINT(x) do {} while (0)
#endif

int main(int argc, char **argv) {
    string readFasta, vRefFasta, dRefFasta, jRefFasta, outFile, paramDir;
    int k = 21, maxKeep = 3, scoringScheme = 0;
    int v_k = -1; int d_k = -1; int j_k = -1;
    string homeDir(getenv("HOME"));
    string default_paramDir(homeDir+"/.nb_params/4mer_amp/");
    try {
	CmdLine cmd("Command description message", ' ', "0.9");

	ValueArg<std::string> readNameArg("r","read_file","Name to print1",true,"homer","string");
	cmd.add( readNameArg );

	ValueArg<std::string> vRefNameArg("v","V_ref_file","Name to print2",true,"blah","string");
	cmd.add( vRefNameArg );

	ValueArg<std::string> dRefNameArg("d","D_ref_file","Name to print2",true,"blah","string");
	cmd.add( dRefNameArg );

	ValueArg<std::string> jRefNameArg("j","J_ref_file","Name to print2",true,"blah","string");
	cmd.add( jRefNameArg );

	ValueArg<int> kSizeArg("k","k_param","size of k",false,-1,"int");
	cmd.add( kSizeArg );
	
	ValueArg<int> vkSizeArg("V","v_k_param","size of k for V",false,-1,"int");
	cmd.add( vkSizeArg );
	ValueArg<int> dkSizeArg("D","d_k_param","size of k for D",false,-1,"int");
	cmd.add( dkSizeArg );
	ValueArg<int> jkSizeArg("J","j_k_param","size of k for J",false,-1,"int");
	cmd.add( jkSizeArg );

	ValueArg<int> maxReportArg("m", "max_report", "maximum gene segements to report, [1,2,..,10], default = 2", false, 2, "int");
	cmd.add(maxReportArg);

	ValueArg<int> scoringArg("s", "scoring", "[0 = standard, 1 = prob]", false, 1, "int");
	cmd.add(scoringArg);

	ValueArg<std::string> paramDirArg("p", "param_dir", "path to parameter directory", false, default_paramDir, "string");
	cmd.add(paramDirArg);

	ValueArg<std::string> outFileNameArg("o","output_file","Name of output file",false,"","string");
	cmd.add( outFileNameArg );
	
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
	// parameter dir
	paramDir = paramDirArg.getValue();

	cerr<<"READS =\t"<<readFasta<<endl;
	cerr<<"REFS =\t"<<vRefFasta<<endl;
	cerr<<"OUTFILE =\t"<<outFile<<endl;
	cerr<<"K =\t"<<k<<endl;
	cerr<<"K_v = "<<v_k<<"\tK_d = "<<d_k<<"\tJ_k = "<<j_k<<endl;
	cerr<<"MAX REPORT =\t"<<maxKeep<<endl;
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
	(MutationNBProbabilities::getInstance()).setParamDir(paramDir);
#ifdef DEBUG
	int lmerLen = (MutationNBProbabilities::getInstance()).getLMerLen();
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
    CreateProfile cp(&cab);
    MAIN_DEBUG_PRINT("SIZE:\t"<<cp.getProfileSize());
    //
    FastaRefID<FastaParser> fp(readFasta);
    fp.openFile();
    int index = 0;
    while(fp.hasNextSequence()) {
	pair<string,string> entry = fp.getNextEntry();
	string readID = fp.getCurrIDLine().substr(1);  // remove the '>'
	string seq = entry.second;
	if(index % 100 == 0) { cerr<<index<<endl; }
	MAIN_DEBUG_PRINT("READID:\t"<<readID);	
	MAIN_DEBUG_PRINT("SEQ:\t"<<seq);
	
	// one-stop shop for computing everything for sequence
	cp.compute(seq);
	
	//ColorProfileMatrix cp_mat = cp.getColorProfile(seq);
        //cp.computeCDR3(seq, cp_mat);
	
	// vector<double> vscores = cp.getPredictedVScores();
	// vector<double> dscores = cp.getPredictedDScores();
	// vector<double> jscores = cp.getPredictedJScores();

	// cout<<"D:\t"<<cp.getDScores()<<endl;
	// cout<<cp_mat.getDColorProfile()<<endl;
	// MAIN_DEBUG_PRINT("V:\t"<<cp.getVScores());
	// MAIN_DEBUG_PRINT("D:\t"<<cp.getDScores());
	// MAIN_DEBUG_PRINT("J:\t"<<cp.getJScores());
	
	// vector<string> vpred = cp.getPredictedV(2);
	// vector<string> dpred = cp.getPredictedD(2);
	// vector<string> jpred = cp.getPredictedJ(2);	
	
	// out<<readID<<"\t"<<vpred<<"\t"<<dpred<<"\t"<<jpred<<"\t"
	//    <<vscores<<"\t"<<dscores<<"\t"<<jscores<<endl;
	out<<readID<<"\t"<<cp;
	
	// MAIN_DEBUG_PRINT("V:\t"<<vpred);
	// MAIN_DEBUG_PRINT("D:\t"<<dpred);
	// MAIN_DEBUG_PRINT("J:\t"<<jpred);
	
	index++;
    }
    fp.closeFile();
    
    if(out_s.is_open()) { out_s.close(); }
    
    // checks process memory usage
    double vm, rss;
    process_mem_usage(vm, rss);
    cerr<<"VM: " << vm <<"\tRSS: " << rss<<endl;
    
    return 0;
}

