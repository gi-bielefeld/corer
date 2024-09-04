#ifndef IO_HPP
#define IO_HPP

#include <bifrost/CompactedDBG.hpp>
#include <bifrost/ColoredCDBG.hpp>

#include <htslib/sam.h>

#include "CoreInfo.h"

#define OPTIONS "i:c:o:f:q:d:t:sah"
//Minimum required parameters are the graph files (sequences and colors), and the output prefix
#define MIN_PARAM_NB 3
#define BIFROST_VERBOSE_MODE false
#define DEFAULT_NB_THREADS 1
#define OUTPUT_CORE_SNIPPETS_DEFAULT false

//This function prints usage infos
inline void dspHlp(){
	cerr << "Corer [-hsa] [-q QUORUM] [-d DELTA] [-i Graph_File] [-c Graph_Color_File] [-o Output_File_Prefix] [-f Core_K-mer_F" <<\
	"ile] [-t Nb_Threads]" << 
	endl << endl;
	cerr << "Extracting a pangenome's core." << endl << endl;
	cerr << "Required parameters:" << endl;
	cerr << "   -i   --igraph  Input graph file in gfa(.gz) or bfg format" << endl;
	cerr << "   -c   --cgraph  Input graph color file in color.bfg format" << endl;
	cerr << "   -o   --ograph  Output graph file prefix" << endl << endl;
	cerr << "Optional parameters with required argument:" << endl;
	cerr << "   -q   --quorum     Absolute quorum defining the core (default is 90\% of pangenome size)" << endl;
	cerr << "   -d   --delta      Maximum distance between core k-mers (default is 50)" << endl;
	cerr << "   -f   --coreKmers  FASTA file containing sequences which k-mers should be used as core k-mers (neglects quorum)" << 
	endl;
	cerr << "   -t   --threads    Number of threads (default is 1)" << endl << endl;
	cerr << "Optional parameters without argument:" << endl;
	cerr << "   -a   --approx    Search for inexact core k-mers (only active if core k-mer file (-f) given)" << endl;
	cerr << "   -s   --snippets  Output unitig core snippets to stdout" << endl;
	cerr << "   -h   --help      Display this help message" << endl;
}

//This function cleans up all temporary data structures created for reading in a FASTA file
inline void cleanUpFASTAparser(sam_hdr_t *&samhdr, samFile *&file, bam1_t *&bamdata){
	if(samhdr) sam_hdr_destroy(samhdr);

	if(file) sam_close(file);

	if(bamdata) bam_destroy1(bamdata);
}

//This function parses the program parameters. Returns false if given arguments are not valid
const bool prsArgs(int& nArgs, char** argList, string& inGfl, string& inCfl, string& outPref, string& iKfl, uint32_t& qrm, uint32_t&
	 dlt, size_t& nThrds, bool& oSnps, bool& apprxSrch);

//This function iterates over the given graph and outputs all core and bridging parts as snippets
void outputSnippets(const ColoredCDBG<CoreInfo>& cdbg);

//This function constructs a graph only consisting of a detected core and writes it to the specified output file
void genCoreGraph(ColoredCDBG<CoreInfo>& cdbg, const string& oName, const size_t& thrds);

//This function reads in a FASTA file and returns whether successful. Sequences of read FASTA records are appended to vector seqs
bool readFasta(const char *filename, vector<string>& seqs);

#endif
