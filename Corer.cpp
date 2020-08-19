#include <bifrost/CompactedDBG.hpp>
#include <bifrost/ColoredCDBG.hpp>

#include "Core.cpp"
#include "Traversal.h"
#include "IO.cpp"

int main(int argc, char **argv){
	uint32_t qrm = 0;
	uint32_t dlt = DEFAULT_DELTA;
	string gFilePref;
	ColoredCDBG<CoreInfo> cdbg = ColoredCDBG<CoreInfo>();

	//Parse arguments
	if(!prsArgs(argc, argv, gFilePref, qrm, dlt)){
		//Display help message
		dspHlp();
		return EXIT_FAILURE;
	}

	//Testing
	// cout << "Loaded parameters" << endl;
	// cout << "Parameters are gFilePref: " << gFilePref << " qrm: " << qrm << " dlt: " << dlt << endl;

	//Load graph
	if(!cdbg.read(gFilePref + GFA_FILE_ENDING, gFilePref + COLOR_FILE_ENDING, true)){
		cerr << "ERROR: Graph could not be loaded" << endl;
		return EXIT_FAILURE;
	}

	//Set quorum if not already done
	if(qrm == 0){
		qrm = min((uint32_t) MIN_QUORUM, (uint32_t) (cdbg.getNbColors() * DEFAULT_CORE_RATIO));
		cerr << "NOTE: No quorum value given; quorum is set to " << qrm << endl;
	}

	//Testing
	// cout << "Colors in graph are " << cdbg.getNbColors() << endl;
	// cout << "First unitig is " << cdbg.begin()->mappedSequenceToString() << endl;
	// cout << "Quorum is " << qrm << endl;

	//Walk through the graph and mark all core parts within each unitig
	markCore(cdbg, qrm, dlt);
	//Walk through the graph and mark all bridging k-mers within each unitig
	detectBrdg(cdbg, dlt);//TODO: Implement this function!

	//Output results...

	return EXIT_SUCCESS;
}