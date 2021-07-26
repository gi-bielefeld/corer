#include <bifrost/CompactedDBG.hpp>
#include <bifrost/ColoredCDBG.hpp>

#include "Core.cpp"
#include "Traversal.h"
#include "IO.cpp"
#include "Bridging.cpp"

int main(int argc, char **argv){
	bool oSnps = OUTPUT_CORE_SNIPPETS_DEFAULT;
	uint32_t qrm = 0;
	uint32_t dlt = DEFAULT_DELTA;
	size_t thrds = DEFAULT_NB_THREADS;
	string iFilePref;
	string oFilePref;
	ColoredCDBG<CoreInfo> cdbg = ColoredCDBG<CoreInfo>();

	//Parse arguments
	if(!prsArgs(argc, argv, iFilePref, oFilePref, qrm, dlt, thrds, oSnps)){
		//Display help message
		dspHlp();
		return EXIT_FAILURE;
	}

	//Load graph
	if(!cdbg.read(iFilePref + GFA_FILE_ENDING, iFilePref + COLOR_FILE_ENDING, thrds, false)){
		cerr << "ERROR: Graph could not be loaded" << endl;
		return EXIT_FAILURE;
	}

	//Set quorum if not already done
	if(qrm == 0){
		qrm = max((uint32_t) MIN_QUORUM, (uint32_t) (cdbg.getNbColors() * DEFAULT_CORE_RATIO));
		cerr << "NOTE: No quorum value given; quorum is set to " << qrm << endl;
	}

	//Walk through the graph and mark all core parts within each unitig
	markCore(cdbg, qrm, dlt);
	//Walk through the graph and mark all bridging k-mers within each unitig
	detectBrdg(cdbg, dlt);

	//Testing
	// return 0;

	//Check if unitig snippet output is requested
	if(oSnps){
		//Output core as unitig snippets
		outputSnippets(cdbg);
	}

	//Construct and write core graph
	genCoreGraph(cdbg, oFilePref, thrds);

	return EXIT_SUCCESS;
}
