#include <bifrost/CompactedDBG.hpp>
#include <bifrost/ColoredCDBG.hpp>

#include "Core.h"
#include "Traversal.h"
#include "IO.cpp"
#include "CoreInfo.h"

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
	cout << "Loaded parameters" << endl;
	cout << "Parameters are gFilePref: " << gFilePref << " qrm: " << qrm << " dlt: " << dlt << endl;

	//Load graph
	if(!cdbg.read(gFilePref + GFA_FILE_ENDING, gFilePref + COLOR_FILE_ENDING, true)){
		cerr << "ERROR: Graph could not be loaded" << endl;
		return EXIT_FAILURE;
	}

	return EXIT_SUCCESS;
}