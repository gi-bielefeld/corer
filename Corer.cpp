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

	//Testing
	// cout << "Loaded graph:" << endl;
	// for(ColoredCDBG<CoreInfo>::iterator u = cdbg.begin(); u != cdbg.end(); ++u) cout << u->mappedSequenceToString() << endl;
	// exit(EXIT_SUCCESS);

	//Set quorum if not already done
	if(qrm == 0){
		qrm = max((uint32_t) MIN_QUORUM, (uint32_t) (cdbg.getNbColors() * DEFAULT_CORE_RATIO));
		cerr << "NOTE: No quorum value given; quorum is set to " << qrm << endl;
	}

	//Testing
	// cout << "Detecting core k-mers...";

	//Walk through the graph and mark all core parts within each unitig
	markCore(cdbg, qrm, dlt);

	//Testing
	// cout << "done" << endl;
	// cout << "Detecting bridging k-mers..." << endl;
	// cout "Our graph is:" << endl;
	// for(ColoredCDBG<CoreInfo>::iterator i = cdbg.begin(); i != cdbg.end(); ++i){
	// 	if(!i->mappedSequenceToString().compare("ATGCTGTTTAA") || !i->mappedSequenceToString().compare("TTAAACAGCAT")) cerr << "main: Found unitig " << i->mappedSequenceToString() << endl << "main: Suffix bridging flag is " << (i->getData()->getData(*i)->sufBrdg ? "" : "not ") << "set" << endl;
	//	cout << i->referenceUnitigToString() << endl;//" Colors are:" << endl;
	//
	// 	for(UnitigColors::const_iterator j = i->getData()->getUnitigColors(*i)->begin(*i); j != i->getData()->getUnitigColors(*i)->end(); ++j) cout << "Pos:" << j.getKmerPosition() << " ID:" << j.getColorID() << endl;
	//
	// 	// cout << "coreList:" << endl;
	// 	// for(list<pair<uint32_t, uint32_t>>::const_iterator k = i->getData()->getData(*i)->coreList.begin(); k != i->getData()->getData(*i)->coreList.end(); ++k) cout << "[" << k->first << "," << k->second << "]" << endl;
	// }
	// return EXIT_SUCCESS;
	//First, only output core k-mers (important for validation script)
	// cout << "Core k-mers:" << endl;
	// outputSnippets(cdbg);

	//Walk through the graph and mark all bridging k-mers within each unitig
	detectBrdg(cdbg, dlt);

	//Testing
	// cerr << "main: sucCoreDist of unitig TGTTAAACAGC: " << cdbg.find(Kmer("TGTTAAACAGC")).getData()->getData(cdbg.find(Kmer("TGTTAAACAGC")))->sucCoreDist << endl;
	// cerr << "main: predCoreDist of unitig TGTTAAACAGC: " << cdbg.find(Kmer("TGTTAAACAGC")).getData()->getData(cdbg.find(Kmer("TGTTAAACAGC")))->predCoreDist << endl;
	// cout << "done" << endl << "Construct output graph...";
	// for(ColoredCDBG<CoreInfo>::iterator i = cdbg.begin(); i != cdbg.end(); ++i){
	// 	cout << "Unitig: " << i->mappedSequenceToString() << " CoreInfo: preBrdg: " << (i->getData()->getData(*i)->preBrdg ? "true" : "false") << " sufBrdg: " << (i->getData()->getData(*i)->sufBrdg ? "true" : "false") << " sucCoreDist: " << i->getData()->getData(*i)->sucCoreDist << " predCoreDist: " << i->getData()->getData(*i)->predCoreDist << " coreList:" << endl;
	// 	for(list<pair<uint32_t, uint32_t>>::const_iterator interval = i->getData()->getData(*i)->coreList.begin(); interval != i->getData()->getData(*i)->coreList.end(); ++interval){
	// 		cout << "[" << interval->first << "," << interval->second << "]" << endl;
	// 	}
	// }
	// return EXIT_SUCCESS;

	//Check if unitig snippet output is requested
	if(oSnps){
		//Testing
		// //For a first sanity check it should be fine to output the core simply as unitig snippets
		// cout << "Whole core:" << endl;

		//Output core as unitig snippets
		outputSnippets(cdbg);
	}

	//Construct and write core graph
	genCoreGraph(cdbg, oFilePref, thrds);

	//Testing
	// cout << "done" << endl;

	return EXIT_SUCCESS;
}
