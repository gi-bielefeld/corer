#ifndef IO_HPP
#define IO_HPP

#include <bifrost/CompactedDBG.hpp>
#include <bifrost/ColoredCDBG.hpp>

#define OPTIONS "g:q:d:t:h"
#define MIN_PARAM_NB 2
#define GFA_FILE_ENDING ".gfa"
#define COLOR_FILE_ENDING ".bfg_colors"
#define DEFAULT_NB_THREADS 1

//This function parses the program parameters. Returns false if given arguments are not valid
const bool prsArgs(int& nArgs, char** argList, string& filePref, uint32_t& qrm, uint32_t& dlt, size_t& nThrds);

//This function prints usage infos
inline void dspHlp(){
	cerr << "Corer [-h] [-q QUORUM] [-d DELTA] [-g Graph_File_Prefix] [-t Nb_Threads]" << endl << endl;
	cerr << "Extract a pangenome's core." << endl << endl;
	cerr << "Required parameters:" << endl;
	cerr << "   -g   --graph   Graph file prefix" << endl << endl;
	cerr << "Optional parameters with required argument:" << endl;
	cerr << "   -q   --quorum   Absolute quorum defining the core (default is 90\% of pangenome size)" << endl << endl;
	cerr << "   -d   --delta    Maximum distance between core k-mers (default is 50)" << endl << endl;
	cerr << "   -t   --threads  Number of threads (default is 1)" << endl << endl;
	cerr << "Optional parameters without argument:" << endl;
	cerr << "   -h   --help   Display help message" << endl;
}

#endif
