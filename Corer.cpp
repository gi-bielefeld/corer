#include <bifrost/CompactedDBG.hpp>
#include <bifrost/ColoredCDBG.hpp>

#include "Core.cpp"
#include "Traversal.h"
#include "IO.cpp"
#include "Bridging.cpp"

int main(int argc, char **argv){
	uint32_t qrm = 0;
	uint32_t dlt = DEFAULT_DELTA;
	size_t thrds = DEFAULT_NB_THREADS;
	string gFilePref;
	ColoredCDBG<CoreInfo> cdbg = ColoredCDBG<CoreInfo>();

	//Parse arguments
	if(!prsArgs(argc, argv, gFilePref, qrm, dlt, thrds)){
		//Display help message
		dspHlp();
		return EXIT_FAILURE;
	}

	//Testing
	// cout << "Loaded parameters" << endl;
	// cout << "Parameters are gFilePref: " << gFilePref << " qrm: " << qrm << " dlt: " << dlt << endl;

	//Load graph
	if(!cdbg.read(gFilePref + GFA_FILE_ENDING, gFilePref + COLOR_FILE_ENDING, thrds, false)){
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

	//Testing
	// cout << "Our graph is:" << endl;
	// for(ColoredCDBG<CoreInfo>::iterator i = cdbg.begin(); i != cdbg.end(); ++i){
	// 	cout << i->referenceUnitigToString() << " Colors are:" << endl;

	// 	for(UnitigColors::const_iterator j = i->getData()->getUnitigColors(*i)->begin(*i); j != i->getData()->getUnitigColors(*i)->end(); ++j) cout << "Pos:" << j.getKmerPosition() << " ID:" << j.getColorID() << endl;

	// 	cout << "coreList:" << endl;
	// 	for(list<pair<uint32_t, uint32_t>>::const_iterator k = i->getData()->getData(*i)->coreList.begin(); k != i->getData()->getData(*i)->coreList.end(); ++k) cout << "[" << k->first << "," << k->second << "]" << endl;
	// }

	//Walk through the graph and mark all bridging k-mers within each unitig
	detectBrdg(cdbg, dlt);

	//Output results...

	//For a first sanity check it should be fine to output the core simply as unitig snippets
	outputSnippets(cdbg);

	return EXIT_SUCCESS;
}