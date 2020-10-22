#ifndef IO_HPP
#define IO_HPP

#include <bifrost/CompactedDBG.hpp>
#include <bifrost/ColoredCDBG.hpp>

#include "CoreInfo.h"

#define OPTIONS "i:o:q:d:t:sh"
#define MIN_PARAM_NB 2
#define GFA_FILE_ENDING ".gfa"
#define COLOR_FILE_ENDING ".bfg_colors"
#define DEFAULT_NB_THREADS 1
#define OUTPUT_CORE_SNIPPETS_DEFAULT false

//This function prints usage infos
inline void dspHlp(){
	cerr << "Corer [-hs] [-q QUORUM] [-d DELTA] [-i Graph_File_Prefix] [-o Output_File_Prefix] [-t Nb_Threads]" << endl << endl;
	cerr << "Extracting a pangenome's core." << endl << endl;
	cerr << "Required parameters:" << endl;
	cerr << "   -i   --igraph  Input graph file prefix" << endl;
	cerr << "   -o   --ograph  Output graph file prefix" << endl << endl;
	cerr << "Optional parameters with required argument:" << endl;
	cerr << "   -q   --quorum   Absolute quorum defining the core (default is 90\% of pangenome size)" << endl;
	cerr << "   -d   --delta    Maximum distance between core k-mers (default is 50)" << endl;
	cerr << "   -t   --threads  Number of threads (default is 1)" << endl << endl;
	cerr << "Optional parameters without argument:" << endl;
	cerr << "   -s   --snippets  Output unitig core snippets to stdout" << endl;
	cerr << "   -h   --help      Display this help message" << endl;
}

//This function parses the program parameters. Returns false if given arguments are not valid
const bool prsArgs(int& nArgs, char** argList, string& inPref, string& outPref, uint32_t& qrm, uint32_t& dlt, size_t& nThrds, bool& oSnps);

//This function iterates over the given graph and outputs all core and bridging parts as snippets
void outputSnippets(const ColoredCDBG<CoreInfo>& cdbg);

#endif
