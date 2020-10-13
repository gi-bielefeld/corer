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
	if(!cdbg.read(gFilePref + GFA_FILE_ENDING, gFilePref + COLOR_FILE_ENDING, thrds, true)){
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
	detectBrdg(cdbg, dlt);

	//Output results...

	//For a first sanity check it should be fine to output the core simply as unitig snippets//
	size_t start, end;

	//Iterate over all unitigs
	for(UnitigColorMap<CoreInfo>::iterator i = cdbg.begin(); i != cdbg.end(); ++i){
		//Check if unitig has no core k-mers
		if(i->getData()->getData(*i)->coreList.empty()){
			//Check if unitig's sequence is marked as bridging
			if(i->getData()->getData(*i)->preBrdg || i->getData()->getData(*i)->sufBrdg){
				//Output the complete sequence
				cout << i->mappedSequenceToString() << endl;
			}
		} else{
			list<pair<uint32_t, uint32_t>>::const_iterator intvl = i->getData()->getData(*i)->coreList.begin();

			//Check if unitig's beginning is marked as bridging
			if(i->getData()->getData(*i)->preBrdg){
				//The first substring we have to output starts at the sequence's beginning
				start = 0;
			} else{
				//The first substring we have to output starts at the first core interval's beginning
				start = intvl->first;
			}

			//We assume that our substring will end at the first core interval's end
			end = intvl->second;
			//Move to next interval
			++intvl;

			//Keep outputting substings as long as intervals are left
			while(intvl != i->getData()->getData(*i)->coreList.end()){
				//Output last substring
				cout << i->mappedSequenceToString().substr(start, end - start + 1) << endl;
				//The next substring starts at the current interval
				start = intvl->first;
				//We assume it ends with the current interval
				end = intvl->second;
				//Move to next interval
				++intvl;
			}

			//Check if unitig's suffix is marked as bridging
			if(i->getData()->getData(*i)->sufBrdg){
				//Output last substring reaching to sequence's end
				cout << i->mappedSequenceToString().substr(start) << endl;
			}
		}
	}

	return EXIT_SUCCESS;
}