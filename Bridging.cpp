#include "Bridging.h"
#include "Traversal.cpp"

//This function checks if the distance from the end of a unitig to either its closest core k-mer or its beginning is already delta or more
const bool lCrTooFar(const size_t& ulen, const list<pair<uint32_t, uint32_t>>& coreList, const uint32_t& dlt){
	if(!coreList.empty()) return ulen - coreList.back().second > dlt;

	return ulen >= dlt;
}

//This function checks if the distance from the beginning of a unitig to either its closest core k-mer or its end is already delta or more
const bool rCrTooFar(const size_t& ulen, const list<pair<uint32_t, uint32_t>>& coreList, const uint32_t& dlt){
	if(!coreList.empty()) return coreList.front().first >= dlt;

	return ulen >= dlt;
}

//This function takes a list of paths obtained by a BFS and a flag indicating whether paths are successive, and marks all involved non-core k-mers as bridging
void markBrdg(const list<Path>& pths, const bool& sucPths){
	//Iterate over paths
	for(list<Path>::const_iterator p = pths.begin(); p != pths.end(); ++p){
		//Testing
		// cout << p->first << ":" << endl;

		//The first unitig of a successor's path potenially has a bridging suffix and vice versa
		list<UnitigColorMap<CoreInfo>>::const_iterator u = p->second.begin();

		if(sucPths){
			u->getData()->getData(*u)->sufBrdg = true;
		} else{
			u->getData()->getData(*u)->preBrdg = true;
		}

		//Testing
		// cout << u->referenceUnitigToString() << endl;
		// if(!u->mappedSequenceToString().compare("CACAATAAAAAA") || !u->mappedSequenceToString().compare("TTTTTTATTGTG")){
		// 	cout << "markBrdg: Found unitig " << u->mappedSequenceToString() << endl << "markBrdg: Suffix bridging flag of unitig ATGCTGTTTAA is " << (u->getGraph()->find(Kmer("ATGCTGTTTAA"))->getData()-getData(u->getGraph()->find(Kmer("ATGCTGTTTAA")))->sufBrdg ? "" : "not ") << "set" << endl;
		// }
		// if(!u->mappedSequenceToString().compare("ATGCTGTTTAA") || !u->mappedSequenceToString().compare("TTAAACAGCAT")){
		// 	cerr << "markBrdg: Found the unitig of interest. It belongs to the path:" << endl;
		// 	for(list<UnitigColorMap<CoreInfo>>::const_iterator v = p->second.begin(); v != p->second.end(); ++v) cerr << v->mappedSequenceToString() << endl;
		// 	cerr << "markBrdg: sucPths is " << (sucPths ? "" : "not ") << "set" << endl;
		// 	// exit(EXIT_SUCCESS);
		// }

		++u;

		//Iterate over unitigs in path
		for(; u != p->second.end(); ++u){
			//Check which bridging flag we have to set
			if(sucPths ^ u->strand){
				//Mark k-mers at unitig's end (reference strand orientation) as bridging
				u->getData()->getData(*u)->sufBrdg = true;
			} else{
				//Mark k-mers at unitig's beginning (reference strand orientation) as bridging
				u->getData()->getData(*u)->preBrdg = true;
			}

			//Testing
			// if(!u->mappedSequenceToString().compare("ATGCTGTTTAA") || !u->mappedSequenceToString().compare("TTAAACAGCAT")){
			// 	cerr << "markBrdg: Found the unitig of interest. It belongs to the path:" << endl;
			// 	for(list<UnitigColorMap<CoreInfo>>::const_iterator v = p->second.begin(); v != p->second.end(); ++v) cerr << v->mappedSequenceToString() << endl;
			// 	cerr << "markBrdg: sucPths is " << (sucPths ? "" : "not ") << "set" << endl;
			// 	// exit(EXIT_SUCCESS);
			// }
		}
	}
}

//This function detects and marks all bridging k-mers between core parts on different unitigs in the graph
void detectBrdg(ColoredCDBG<CoreInfo>& cdbg, const uint32_t& dlt){
	//Length of the path from a core k-mer to the current unitigs beginning
	uint32_t exstPthLen;
	//Pointer to current unitig's CoreInfo object
	CoreInfo* cInfo;
	//List of paths leading to reachable core k-mers on successive unitigs
	list<Path> sucPaths;
	//List of paths leading to reachable core k-mers on predecessive unitigs
	list<Path> predPaths;

	//Testing
	// UnitigColorMap<CoreInfo> lastUnitig;

	//Iterate over unitigs
	for(ColoredCDBG<CoreInfo>::iterator i = cdbg.begin(); i != cdbg.end(); ++i){
		//Get CoreInfo object
		cInfo = i->getData()->getData(*i);

		//Testing
		// if(cdbg.find(Kmer("ATGCTGTTTAA")).getData()->getData(cdbg.find(Kmer("ATGCTGTTTAA")))->sufBrdg){
		// 	cerr << "detectBrdg: Suffix bridging flag has been set" << endl << "Last unitig was " << lastUnitig.mappedSequenceToString() << endl;
		// 	exit(EXIT_SUCCESS);
		// }
		// lastUnitig = *i;
		// if(!i->mappedSequenceToString().compare("ATGCTGTTTAA") || !i->mappedSequenceToString().compare("TTAAACAGCAT")){
		// 	cerr << "detectBrdg: Found unitig " << i->mappedSequenceToString() << endl;
		// 	cerr << "detectBrdg: Suffix bridging flag is " << (cInfo->sufBrdg ? "" : "not ") << "set" << endl;

		// 	// exit(EXIT_SUCCESS);
		// }

		//Check if last k-mer on unitig is neither marked as bridging nor core and ensure that the distance we have to bridge to the left side (i.e. the distance to the closest core k-mer on this unitig or the unitig's beginning) is not already too large
		if((cInfo->coreList.empty() || cInfo->coreList.back().second < i->len - 1) && !cInfo->sufBrdg && !lCrTooFar(i->len, cInfo->coreList, dlt)){
			//Testing
			// if(!i->mappedSequenceToString().compare("ATGCTGTTTAA") || !i->mappedSequenceToString().compare("TTAAACAGCAT")){
			// 	cerr << "detectBrdg: Do BFS on successors" << endl;
			// }

			//Do BFS on successive unitigs and check if we need to try a BFS on predecessors as well (which is the case only if there is a core k-mer on the current unitig or the BFS on successive unitigs was successful)
			if(!doSucBFS(*i, (dlt + 1) / 2, sucPaths) && cInfo->coreList.empty()) continue;
		}

		//Testing
		// uint32_t counter = 0;
		// if(!i->mappedSequenceToString().compare("ATGCTGTTTAA") || !i->mappedSequenceToString().compare("TTAAACAGCAT")){
		// 	cerr << "detectBrdg: " << sucPaths.size() << " paths on successors have been found" << endl;
			
		// 	for(list<Path>::const_iterator p = sucPaths.begin(); p != sucPaths.end(); ++p){
		// 		cerr << "Path " << ++counter << " length (in k-mers): " << p->first << endl;

		// 		for(list<UnitigColorMap<CoreInfo>>::const_iterator u = p->second.begin(); u != p->second.end(); ++u) cerr << u->mappedSequenceToString() << " coreList (min start, max end): " << (u->getData()->getData(*u)->coreList.empty() ? -1 : u->getData()->getData(*u)->coreList.front().first) << " " << (u->getData()->getData(*u)->coreList.empty() ? -1 : u->getData()->getData(*u)->coreList.back().second) << endl;
		// 	}
		// 	exit(EXIT_SUCCESS);
		// }
			

		//Check if first k-mer on unitig is neither marked as bridging nor core and ensure that the distance we have to bridge to the right side (i.e. the distance to the closest core kmer on this unitig or the unitig's end) is not already too large
		if((cInfo->coreList.empty() || cInfo->coreList.front().first > 0) && !cInfo->preBrdg && !rCrTooFar(i->len, cInfo->coreList, dlt)){
			//Check if there are no core k-mers on this unitig
			if(cInfo->coreList.empty()){
				//The length of the existing path the minimum length path found during BFS on successive unitigs in addition to all k-mers on the current unitig (except for the one we start at)
				exstPthLen = findMinPthLen(sucPaths) + i->len - 1;
			} else{
				//A path has to reach the leftmost core k-mer on this unitig only
				exstPthLen = cInfo->coreList.front().first;
			}

			//Do BFS on predecessive unitigs and mark all bridging k-mers if necessary
			if(doPredBFS(*i, min((dlt + 1) / 2, (uint32_t) (dlt - exstPthLen)), predPaths) || !cInfo->coreList.empty()){
				//Testing
				// cout << i->referenceUnitigToString() << endl;
				// cout << sucPaths.size() << endl;
				// if(!i->mappedSequenceToString().compare("ATGCTGTTTAA") || !i->mappedSequenceToString().compare("TTAAACAGCAT")){
				// 	// cerr << "detectBrdg: Found unitig " << i->mappedSequenceToString() << endl;
				// 	// cerr << "detectBrdg: " << sucPaths.size() << " paths on successors have been found" << endl;
				// 	// uint32_t counter = 0;
			
				// 	// for(list<Path>::const_iterator p = sucPaths.begin(); p != sucPaths.end(); ++p){
				// 	// 	cerr << "Path " << ++counter << " length (in k-mers): " << p->first << endl;

				// 	// 	for(list<UnitigColorMap<CoreInfo>>::const_iterator u = p->second.begin(); u != p->second.end(); ++u) cerr << u->mappedSequenceToString() << " coreList (min start, max end): " << (u->getData()->getData(*u)->coreList.empty() ? -1 : u->getData()->getData(*u)->coreList.front().first) << " " << (u->getData()->getData(*u)->coreList.empty() ? -1 : u->getData()->getData(*u)->coreList.back().second) << endl;
				// 	// }

				// 	cerr << "detectBrdg: " << predPaths.size() << " paths on predecessors have been found" << endl;
				// 	counter = 0;
			
				// 	for(list<Path>::const_iterator p = predPaths.begin(); p != predPaths.end(); ++p){
				// 		cerr << "Path " << ++counter << " length (in k-mers): " << p->first << endl;

				// 		for(list<UnitigColorMap<CoreInfo>>::const_iterator u = p->second.begin(); u != p->second.end(); ++u) cerr << u->mappedSequenceToString() << " coreList (min start, max end): " << (u->getData()->getData(*u)->coreList.empty() ? -1 : u->getData()->getData(*u)->coreList.front().first) << " " << (u->getData()->getData(*u)->coreList.empty() ? -1 : u->getData()->getData(*u)->coreList.back().second) << endl;
				// 	}

				// 	exit(EXIT_SUCCESS);
				// }

				//Mark k-mers of successive result paths as bridging
				markBrdg(sucPaths, true);
				//Mark k-mers of predecessive result paths as bridging
				markBrdg(predPaths, false);
			}
		}

		//Clear path lists
		sucPaths.clear();
		predPaths.clear();
	}
}