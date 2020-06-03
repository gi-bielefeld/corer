#ifndef IO_HPP
#define IO_HPP

#define OPTIONS "g:q:d:h"
#define MIN_PARAM_NB 2

//This function parses the program parameters. Returns false if given arguments are not valid
const bool prsArgs(int& nb_args, char** argList, string& filePref, uint32_t& qrm, uint32_t& dlt);

//This function prints usage infos
inline void dspHlp(){
	cerr << "Corer [-h] [-q QUORUM] [-d DELTA] [-g Graph_File_Prefix]" << endl << endl;
	cerr << "Extract a pangenome's core." << endl << endl;
	cerr << "Required parameters:" << endl;
	cerr << "   -g   --graph   Graph file prefix" << endl << endl;
	cerr << "Optional parameters with required argument:" << endl;
	cerr << "   -q   --quorum   Absolute quorum defining the core (default is 90\% of pangenome size)" << endl << endl;
	cerr << "   -d   --delta    Maximum distance between core k-mers (default is 50)" << endl << endl;
	cerr << "Optional parameters without argument:" << endl;
	cerr << "   -h   --help   Display help message" << endl;
}

#endif
