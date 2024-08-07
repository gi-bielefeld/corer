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
	string iGFile, iCFile, iKFile;
	string oFilePref;
	vector<string> seqList;
	ColoredCDBG<CoreInfo> cdbg = ColoredCDBG<CoreInfo>();
	TravTrackQueue queue;

	//Parse arguments
	if(!prsArgs(argc, argv, iGFile, iCFile, oFilePref, iKFile, qrm, dlt, thrds, oSnps)){
		//Display help message
		dspHlp();
		return EXIT_FAILURE;
	}

	//Load graph
	if(!cdbg.read(iGFile, iCFile, thrds, BIFROST_VERBOSE_MODE)){
		cerr << "ERROR: Graph could not be loaded" << endl;
		return EXIT_FAILURE;
	}

	//Check if a core k-mer file is given
	if(!iKFile.empty()){
		if(!readFasta(iKFile.c_str(), seqList)) return EXIT_FAILURE;

		//Mark all occurring k-mers as core
		markKmers(cdbg, seqList, dlt);
		//Initialize queue for graph traversal
		queue = initializeQueue(cdbg);
	} else{
		//Set quorum if not already done
		if(qrm == 0){
			qrm = max((uint32_t) MIN_QUORUM, (uint32_t) (cdbg.getNbColors() * DEFAULT_CORE_RATIO));
			cerr << "NOTE: No quorum value given; quorum is set to " << qrm << endl;
		}

		//Detect all core k-mers
		queue = detectCore(cdbg, qrm, dlt);
	}

	//Annotate unitigs with distances to next core k-mers
	annotateDists(cdbg, queue, dlt);

	//Walk through the graph and mark all bridging k-mers within each unitig
	markBrdg(cdbg, dlt);

	//Check if unitig snippet output is requested
	if(oSnps){
		//Output core as unitig snippets
		outputSnippets(cdbg);
	}

	//Construct and write core graph
	genCoreGraph(cdbg, oFilePref, thrds);

	return EXIT_SUCCESS;
}
