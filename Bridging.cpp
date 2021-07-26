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
		//The first unitig of a successor's path potenially has a bridging suffix and vice versa
		list<UnitigColorMap<CoreInfo>>::const_iterator u = p->second.begin();

		//Testing
		// if(!u->referenceUnitigToString().compare("CACCCGTCTGCGGCCAGATTATTCGCCTGTCGCATCGCCCG")){
		// 	cout << "markBrdg: Found a path in which k-mer CACCCGTCTGCGGCCAGATTA appears:" << endl;
		// 	for(list<UnitigColorMap<CoreInfo>>::const_iterator v = p->second.begin(); v != p->second.end(); ++v) cout << v->mappedSequenceToString() << endl;
		// }

		if(sucPths){
			u->getData()->getData(*u)->sufBrdg = true;
		} else{
			u->getData()->getData(*u)->preBrdg = true;
		}

		++u;

		//Iterate over unitigs in path
		for(; u != p->second.end(); ++u){
			//Testing
			// if(!u->referenceUnitigToString().compare("CACCCGTCTGCGGCCAGATTATTCGCCTGTCGCATCGCCCG")){
			// 	cout << "markBrdg: Found a path in which k-mer CACCCGTCTGCGGCCAGATTA appears:" << endl;
			// 	for(list<UnitigColorMap<CoreInfo>>::const_iterator v = p->second.begin(); v != p->second.end(); ++v) cout << v->mappedSequenceToString() << endl;
			// }

			//Check which bridging flag we have to set
			if(sucPths ^ u->strand){
				//Mark k-mers at unitig's end (reference strand orientation) as bridging
				u->getData()->getData(*u)->sufBrdg = true;
			} else{
				//Mark k-mers at unitig's beginning (reference strand orientation) as bridging
				u->getData()->getData(*u)->preBrdg = true;
			}
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

	//Iterate over unitigs
	for(ColoredCDBG<CoreInfo>::iterator i = cdbg.begin(); i != cdbg.end(); ++i){
		//Get CoreInfo object
		cInfo = i->getData()->getData(*i);

		//Check if last k-mer on unitig is neither marked as bridging nor core and ensure that the distance we have to bridge to the left side (i.e. the distance to the closest core k-mer on this unitig or the unitig's beginning) is not already too large
		if((cInfo->coreList.empty() || cInfo->coreList.back().second < i->len - 1) && !cInfo->sufBrdg && !lCrTooFar(i->len, cInfo->coreList, dlt)){
			//Do BFS on successive unitigs and check if we need to try a BFS on predecessors as well (which is the case only if there is a core k-mer on the current unitig or the BFS on successive unitigs was successful)
			if(!doSucBFS(*i, (dlt + 1) / 2, sucPaths) && cInfo->coreList.empty()) continue;
		}

		//Check if first k-mer on unitig is neither marked as bridging nor core and ensure that the distance we have to bridge to the right side (i.e. the distance to the closest core kmer on this unitig or the unitig's end) is not already too large
		if((cInfo->coreList.empty() || cInfo->coreList.front().first > 0) && !cInfo->preBrdg && !rCrTooFar(i->len, cInfo->coreList, dlt)){
			//Check if there are no core k-mers on this unitig
			if(cInfo->coreList.empty()){
				//The length of the existing path is the minimum length path found during BFS on successive unitigs in addition to all k-mers on the current unitig (except for
				//the one we start at)
				exstPthLen = findMinPthLen(sucPaths) + i->len - 1;
			} else{
				//A path has to reach the leftmost core k-mer on this unitig only
				exstPthLen = cInfo->coreList.front().first;
			}

			//We have to catch this to prevent an overflow
			if(dlt < exstPthLen) exstPthLen = dlt;

			//Testing
			// if(!i->referenceUnitigToString().compare("CACCCGTCTGCGGCCAGATTATTCGCCTGTCGCATCGCCCG")){
			// 	cout << "detectBrdg: We are dealing with unitig " << i->referenceUnitigToString() << endl;
			// 	cout << "detectBrdg: exstPthLen: " << exstPthLen << " (uint32_t) (dlt - exstPthLen):" << (uint32_t) (dlt - exstPthLen) << endl;
			// }

			//Do BFS on predecessive unitigs and mark all bridging k-mers if necessary
			if(doPredBFS(*i, min((dlt + 1) / 2, (uint32_t) (dlt - exstPthLen)), predPaths) || !cInfo->coreList.empty()){
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